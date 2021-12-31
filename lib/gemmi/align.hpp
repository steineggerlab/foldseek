// Copyright 2020 Global Phasing Ltd.
//
// Sequence alignment, label_seq_id assignment, structure superposition.

#ifndef GEMMI_ALIGN_HPP_
#define GEMMI_ALIGN_HPP_

#include "model.hpp"
#include "seqalign.hpp"  // for align_sequences
#include "qcp.hpp"       // for superpose_positions
#include "polyheur.hpp"  // for are_connected3

namespace gemmi {

// Sequence alignment and label_seq_id assignment

// helper function for sequence alignment
inline std::vector<bool> prepare_free_gapo(const ConstResidueSpan& polymer,
                                           PolymerType polymer_type) {
  std::vector<bool> gaps;
  gaps.reserve(polymer.size());
  gaps.push_back(true); // free gap opening at the beginning of sequence
  if (!is_polypeptide(polymer_type) && !is_polynucleotide(polymer_type))
    return gaps;
  auto first_conformer = polymer.first_conformer();
  auto res = first_conformer.begin();
  for (auto next_res = res; ++next_res != first_conformer.end(); res = next_res)
    gaps.push_back(!are_connected3(*res, *next_res, polymer_type));
  return gaps;
}

// pre: !!polymer
inline AlignmentResult align_sequence_to_polymer(
                                     const std::vector<std::string>& full_seq,
                                     const ConstResidueSpan& polymer,
                                     PolymerType polymer_type,
                                     const AlignmentScoring& scoring) {
  std::map<std::string, std::uint8_t> encoding;
  for (const std::string& res_name : scoring.matrix_encoding)
    encoding.emplace(res_name, (std::uint8_t)encoding.size());
  for (const Residue& res : polymer)
    encoding.emplace(res.name, (std::uint8_t)encoding.size());
  for (const std::string& mon_list : full_seq)
    encoding.emplace(Entity::first_mon(mon_list), (std::uint8_t)encoding.size());
  if (encoding.size() > 255)
    return AlignmentResult();

  std::vector<std::uint8_t> encoded_full_seq(full_seq.size());
  for (size_t i = 0; i != full_seq.size(); ++i)
    encoded_full_seq[i] = encoding.at(Entity::first_mon(full_seq[i]));

  std::vector<std::uint8_t> encoded_model_seq;
  encoded_model_seq.reserve(polymer.size());
  for (const Residue& res : polymer.first_conformer())
    encoded_model_seq.push_back(encoding.at(res.name));

  return align_sequences(encoded_full_seq, encoded_model_seq,
                         prepare_free_gapo(polymer, polymer_type),
                         (std::uint8_t)encoding.size(), scoring);
}

// check for exact match between model sequence and full sequence (SEQRES)
inline bool seqid_matches_seqres(const ConstResidueSpan& polymer,
                                 const Entity& ent) {
  if (ent.full_sequence.size() != polymer.size())
    return false;
  int idx = 0;
  for (const Residue& res : polymer) {
    if (ent.full_sequence[idx] != res.name ||
        ++idx != *res.seqid.num || res.seqid.has_icode())
      return false;
  }
  return true;
}

// Uses sequence alignment (model to SEQRES) to assign label_seq.
// force: assign label_seq even if full sequence is not known (assumes no gaps)
inline void assign_label_seq_to_polymer(ResidueSpan& polymer,
                                        const Entity* ent, bool force) {
  AlignmentResult result;

  // sequence not known
  if (!ent || ent->full_sequence.empty()) {
    if (!force)
      return;
    PolymerType ptype = (ent && ent->polymer_type != PolymerType::Unknown
                         ? ent->polymer_type
                         : check_polymer_type(polymer));
    const Residue* prev = nullptr;
    for (const Residue& res : polymer.first_conformer()) {
      if (prev && !are_connected3(*prev, res, ptype))
        result.push_cigar(1, 1);  // assume a single insertion
      result.push_cigar(0, 1);
      prev = &res;
    }

  // exact match - common case that doesn't require alignment
  } else if (seqid_matches_seqres(polymer, *ent)) {
    result.push_cigar(0, (int)ent->full_sequence.size());

  // sequence alignment
  } else {
    AlignmentScoring scoring;
    result = align_sequence_to_polymer(ent->full_sequence, polymer,
                                       ent->polymer_type, scoring);
  }

  auto res_group = polymer.first_conformer().begin();
  int id = 1;
  for (AlignmentResult::Item item : result.cigar) {
    switch (item.op()) {
      case 'I':
        id += item.len();
        break;
      case 'D':  // leaving label_seq as it is
        for (uint32_t i = 0; i < item.len(); ++i)
          res_group++;
        break;
      case 'M':  // not checking for mismatches
        for (uint32_t i = 0; i < item.len(); ++i, ++id)
          for (Residue* res = &*res_group++; res != &*res_group; ++res)
            res->label_seq = id;
        break;
    }
  }
}

inline void clear_label_seq_id(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (Residue& res : chain.residues)
        res.label_seq = Residue::OptionalNum();
}

inline void assign_label_seq_id(Structure& st, bool force) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      if (ResidueSpan polymer = chain.get_polymer())
        if (!polymer.front().label_seq || !polymer.back().label_seq) {
          const Entity* ent = st.get_entity_of(polymer);
          assign_label_seq_to_polymer(polymer, ent, force);
        }
}


// superposition

enum class SupSelect {
  CaP,  // only Ca (aminoacids) or P (nucleotides) atoms
  MainChain,  // only main chain atoms
  All
};

inline void prepare_positions_for_superposition(std::vector<Position>& pos1,
                                                std::vector<Position>& pos2,
                                                ConstResidueSpan fixed,
                                                ConstResidueSpan movable,
                                                PolymerType ptype,
                                                SupSelect sel,
                                                char altloc='\0') {
  AlignmentScoring scoring;
  AlignmentResult result = align_sequence_to_polymer(fixed.extract_sequence(),
                                                     movable, ptype, scoring);
  auto it1 = fixed.first_conformer().begin();
  auto it2 = movable.first_conformer().begin();
  std::vector<AtomNameElement> used_atoms;
  if (sel == SupSelect::CaP) {
    if (is_polynucleotide(ptype))
      used_atoms.push_back({"P", El::P});
    else
      used_atoms.push_back({"CA", El::C});
  } else if (sel == SupSelect::MainChain) {
    used_atoms = get_mainchain_atoms(ptype);
  }
  for (AlignmentResult::Item item : result.cigar) {
    char op = item.op();
    for (uint32_t i = 0; i < item.len(); ++i) {
      if (op == 'M' && it1->name == it2->name) {
        if (!used_atoms.empty()) {
          for (const AtomNameElement& ane : used_atoms) {
            const Atom* a1 = it1->find_atom(ane.atom_name, altloc, ane.el);
            const Atom* a2 = it2->find_atom(ane.atom_name, altloc, ane.el);
            if (a1 && a2) {
              pos1.push_back(a1->pos);
              pos2.push_back(a2->pos);
            }
          }
        } else {
          for (const Atom& a1 : it1->atoms)
            if (a1.altloc_matches(altloc))
              if (const Atom* a2 = it2->find_atom(a1.name, altloc, a1.element)) {
                pos1.push_back(a1.pos);
                pos2.push_back(a2->pos);
              }
        }
      }
      if (op == 'M' || op == 'I')
        ++it1;
      if (op == 'M' || op == 'D')
        ++it2;
    }
  }
}

inline SupResult calculate_current_rmsd(ConstResidueSpan fixed,
                                        ConstResidueSpan movable,
                                        PolymerType ptype,
                                        SupSelect sel,
                                        char altloc='\0') {
  std::vector<Position> pos1, pos2;
  prepare_positions_for_superposition(pos1, pos2, fixed, movable, ptype, sel, altloc);
  SupResult r;
  r.count = pos1.size();
  double sd = 0;
  for (size_t i = 0; i != pos1.size(); ++i)
    sd += pos1[i].dist_sq(pos2[i]);
  r.rmsd = std::sqrt(sd / r.count);
  return r;
}

inline SupResult calculate_superposition(ConstResidueSpan fixed,
                                         ConstResidueSpan movable,
                                         PolymerType ptype,
                                         SupSelect sel,
                                         int trim_cycles=0,
                                         double trim_cutoff=2.0,
                                         char altloc='\0') {
  std::vector<Position> pos1, pos2;
  prepare_positions_for_superposition(pos1, pos2, fixed, movable, ptype, sel, altloc);
  const double* weights = nullptr;
  size_t len = pos1.size();
  SupResult sr = superpose_positions(pos1.data(), pos2.data(), len, weights);

  for (int n = 0; n < trim_cycles; ++n) {
    double max_dist_sq = sq(trim_cutoff * sr.rmsd);
    size_t p = 0;
    for (size_t i = 0; i != len; ++i) {
      Vec3 m2 = sr.transform.apply(pos2[i]);
      if (m2.dist_sq(pos1[i]) <= max_dist_sq) {
        if (i != p) {
          pos1[p] = pos1[i];
          pos2[p] = pos2[i];
        }
        ++p;
      }
    }
    if (p == len)
      break;
    len = p;
    if (len < 3)
      fail("in calculate_superposition(): only ", std::to_string(len),
           " atoms after trimming");
    sr = superpose_positions(pos1.data(), pos2.data(), len, weights);
  }

  return sr;
}

} // namespace gemmi
#endif
