// Copyright 2018 Global Phasing Ltd.
//
// Topo(logy) - restraints (from a monomer library) applied to a model.

#ifndef GEMMI_TOPO_HPP_
#define GEMMI_TOPO_HPP_

#include "chemcomp.hpp"  // for ChemComp
#include "monlib.hpp"    // for MonLib
#include "model.hpp"     // for Residue, Atom
#include "calculate.hpp" // for calculate_angle, calculate_dihedral
#include "polyheur.hpp"  // for are_connected

namespace gemmi {

struct Topo {
  // We have internal pointers in this class (pointers setup in
  // apply_restraints() that point to ResInfo::chemcomp.rt),
  // disable copying this class.
  Topo() = default;
  Topo(Topo const&) = delete;
  Topo& operator=(Topo const&) = delete;

  struct Bond {
    const Restraints::Bond* restr;
    std::array<Atom*, 2> atoms;
    double calculate() const { return atoms[0]->pos.dist(atoms[1]->pos); }
    double calculate_z() const {
      return std::abs(calculate() - restr->value) / restr->esd;
    }
  };
  struct Angle {
    const Restraints::Angle* restr;
    std::array<Atom*, 3> atoms;
    double calculate() const {
      return calculate_angle(atoms[0]->pos, atoms[1]->pos, atoms[2]->pos);
    }
    double calculate_z() const { return angle_z(calculate(), *restr); }
  };
  struct Torsion {
    const Restraints::Torsion* restr;
    std::array<Atom*, 4> atoms;
    double calculate() const {
      return calculate_dihedral(atoms[0]->pos, atoms[1]->pos,
                                atoms[2]->pos, atoms[3]->pos);
    }
    double calculate_z() const {
      return angle_z(calculate(), *restr, 360. / std::max(1, restr->period));
    }
  };
  struct Chirality {
    const Restraints::Chirality* restr;
    std::array<Atom*, 4> atoms;
    double calculate() const {
      return calculate_chiral_volume(atoms[0]->pos, atoms[1]->pos,
                                     atoms[2]->pos, atoms[3]->pos);
    }
    bool check() const { return !restr->is_wrong(calculate()); }
  };
  struct Plane {
    const Restraints::Plane* restr;
    std::vector<Atom*> atoms;
    bool has(const Atom* atom) const {
      return in_vector(const_cast<Atom*>(atom), atoms);
    }
  };

  enum class Provenance { None, PrevLink, Monomer, NextLink, ExtraLink };
  enum class RKind { Bond, Angle, Torsion, Chirality, Plane };
  struct Force {
    Provenance provenance;
    RKind rkind;
    size_t index; // index in the respective vector (bonds, ...) in Topo
  };

  struct ResInfo {
    Residue* res;
    struct Prev {
      std::string link;
      int idx; // relative to current index, 0 means n/a
      const ResInfo* get(const ResInfo* t) const { return t + idx; }
      ResInfo* get(ResInfo* t) const { return t + idx; }
    };
    std::vector<Prev> prev;
    std::vector<std::string> mods;
    ChemComp chemcomp;
    std::vector<Force> forces;

    ResInfo(Residue* r) : res(r) {}
    void add_mod(const std::string& m) {
      if (!m.empty())
        mods.push_back(m);
    }
  };

  struct ChainInfo {
    std::string name;
    std::string entity_id;
    bool polymer;
    PolymerType polymer_type;
    std::vector<ResInfo> res_infos;

    void initialize(ResidueSpan& subchain, const Entity* ent);
    void setup_polymer_links();
    void add_refmac_builtin_modifications();
    struct RGroup {
      std::vector<ResInfo>::iterator begin, end;
    };
    RGroup group_from(std::vector<ResInfo>::iterator b) const {
      auto e = b + 1;
      while (e != res_infos.end() && e->res->group_key() == b->res->group_key())
        ++e;
      return RGroup{b, e};
    }
  };

  struct ExtraLink {
    Residue* res1;
    Residue* res2;
    char alt1 = '\0';
    char alt2 = '\0';
    std::string link_id;
    std::vector<Force> forces;
  };

  template<typename T>
  static int has_atom(const Atom* a, const T& t) {
    for (int i = 0; (size_t) i != t.atoms.size(); ++i)
      if (t.atoms[i] == a)
        return i;
    return -1;
  }

  std::vector<ChainInfo> chain_infos;
  std::vector<ExtraLink> extras;

  // Restraints applied to Model
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::vector<Plane> planes;

  ResInfo* find_resinfo(const Residue* res) {
    for (ChainInfo& ci : chain_infos)
      for (ResInfo& ri : ci.res_infos)
        if (ri.res == res)
          return &ri;
    return nullptr;
  }

  const Restraints::Bond* take_bond(const Atom* a, const Atom* b) const {
    for (const Bond& bond : bonds)
      if ((bond.atoms[0] == a && bond.atoms[1] == b) ||
          (bond.atoms[0] == b && bond.atoms[1] == a))
        return bond.restr;
    return nullptr;
  }

  const Restraints::Angle* take_angle(const Atom* a,
                                      const Atom* b,
                                      const Atom* c) const {
    for (const Angle& ang : angles)
      if (ang.atoms[1] == b && ((ang.atoms[0] == a && ang.atoms[2] == c) ||
                                (ang.atoms[0] == c && ang.atoms[2] == a)))
        return ang.restr;
    return nullptr;
  }

  const Chirality* get_chirality(const Atom* ctr) const {
    for (const Chirality& chir : chirs)
      if (chir.atoms[0] == ctr)
        return &chir;
    return nullptr;
  }

  std::vector<Force> apply_restraints(const Restraints& rt,
                                      Residue& res, Residue* res2,
                                      char altloc='*') {
    std::string altlocs;
    if (altloc == '*') {
      // find all distinct altlocs
      add_distinct_altlocs(res, altlocs);
      if (res2)
        add_distinct_altlocs(*res2, altlocs);
    }
    if (altlocs.empty())
      altlocs += altloc;

    std::vector<Force> forces;
    Provenance pro = Provenance::None;
    for (const Restraints::Bond& bond : rt.bonds)
      for (char alt : altlocs)
        if (Atom* at1 = bond.id1.get_from(res, res2, alt))
          if (Atom* at2 = bond.id2.get_from(res, res2, alt)) {
            forces.push_back({pro, RKind::Bond, bonds.size()});
            bonds.push_back({&bond, {{at1, at2}}});
            if (!at1->altloc && !at2->altloc)
              break;
          }
    for (const Restraints::Angle& angle : rt.angles)
      for (char alt : altlocs)
        if (Atom* at1 = angle.id1.get_from(res, res2, alt))
          if (Atom* at2 = angle.id2.get_from(res, res2, alt))
            if (Atom* at3 = angle.id3.get_from(res, res2, alt)) {
              forces.push_back({pro, RKind::Angle, angles.size()});
              angles.push_back({&angle, {{at1, at2, at3}}});
              if (!at1->altloc && !at2->altloc && !at3->altloc)
                break;
            }
    for (const Restraints::Torsion& tor : rt.torsions)
      for (char alt : altlocs)
        if (Atom* at1 = tor.id1.get_from(res, res2, alt))
          if (Atom* at2 = tor.id2.get_from(res, res2, alt))
            if (Atom* at3 = tor.id3.get_from(res, res2, alt))
              if (Atom* at4 = tor.id4.get_from(res, res2, alt)) {
                forces.push_back({pro, RKind::Torsion, torsions.size()});
                torsions.push_back({&tor, {{at1, at2, at3, at4}}});
                if (!at1->altloc && !at2->altloc &&
                    !at3->altloc && !at4->altloc)
                  break;
          }
    for (const Restraints::Chirality& chir : rt.chirs)
      for (char alt : altlocs)
        if (Atom* at1 = chir.id_ctr.get_from(res, res2, alt))
          if (Atom* at2 = chir.id1.get_from(res, res2, alt))
            if (Atom* at3 = chir.id2.get_from(res, res2, alt))
              if (Atom* at4 = chir.id3.get_from(res, res2, alt)) {
                forces.push_back({pro, RKind::Chirality, chirs.size()});
                chirs.push_back({&chir, {{at1, at2, at3, at4}}});
                if (!at1->altloc && !at2->altloc &&
                    !at3->altloc && !at4->altloc)
                  break;
              }
    for (const Restraints::Plane& plane : rt.planes)
      for (char alt : altlocs) {
        std::vector<Atom*> atoms;
        for (const Restraints::AtomId& id : plane.ids)
          if (Atom* atom = id.get_from(res, res2, alt))
            atoms.push_back(atom);
        if (atoms.size() >= 4) {
          forces.push_back({pro, RKind::Plane, planes.size()});
          planes.push_back({&plane, atoms});
        }
        if (std::all_of(atoms.begin(), atoms.end(),
                        [](Atom* a) { return !a->altloc; }))
          break;
      }
    return forces;
  }

  void apply_internal_restraints_to_residue(ResInfo& ri) {
    auto forces = apply_restraints(ri.chemcomp.rt, *ri.res, nullptr);
    for (const auto& f : forces)
      ri.forces.push_back({Provenance::Monomer, f.rkind, f.index});
  }

  void apply_restraints_to_residue(ResInfo& ri, const MonLib& monlib) {
    for (ResInfo::Prev& prev : ri.prev) {
      if (const ChemLink* link = monlib.find_link(prev.link)) {
        ResInfo* prev_ri = prev.get(&ri);
        auto forces = apply_restraints(link->rt, *prev_ri->res, ri.res);
        for (const auto& f : forces) {
          ri.forces.push_back({Provenance::PrevLink, f.rkind, f.index});
          prev_ri->forces.push_back({Provenance::NextLink, f.rkind, f.index});
        }
      }
    }
    apply_internal_restraints_to_residue(ri);
  }

  void apply_restraints_to_extra_link(ExtraLink& link, const MonLib& monlib) {
    const ChemLink* cl = monlib.find_link(link.link_id);
    if (!cl) {
      printf("Warning: ignoring link '%s' as it is not in the monomer library",
             link.link_id.c_str());
      return;
    }
    if (link.alt1 && link.alt2 && link.alt1 != link.alt2)
      printf("Warning: LINK between different conformers %c and %c.",
             link.alt1, link.alt2);
    char alt = link.alt1 ? link.alt1 : link.alt2;
    ResInfo* ri1 = find_resinfo(link.res1);
    ResInfo* ri2 = find_resinfo(link.res2);
    auto forces = apply_restraints(cl->rt, *link.res1, link.res2, alt);
    for (Force& f : forces) {
      f.provenance = Provenance::ExtraLink;
      link.forces.push_back(f);
      if (ri1)
        ri1->forces.push_back(f);
      if (ri2)
        ri2->forces.push_back(f);
    }
  }

  // Model is non-const b/c we store non-const pointers to residues in Topo.
  // Because of the pointers, don't add or remove residues after this step.
  // Monlib may get modified by addition of extra links from the model.
  void initialize_refmac_topology(const Structure& st, Model& model0,
                                  MonLib& monlib);

  // This step stores pointers to gemmi::Atom's from model0,
  // so after this step don't add or remove atoms.
  // monlib is needed only for links.
  void finalize_refmac_topology(const MonLib& monlib) {
    for (ChainInfo& chain_info : chain_infos)
      for (ResInfo& ri : chain_info.res_infos)
        apply_restraints_to_residue(ri, monlib);
    for (ExtraLink& link : extras)
      apply_restraints_to_extra_link(link, monlib);
  }
};

inline void Topo::ChainInfo::initialize(ResidueSpan& subchain, const Entity* ent) {
  res_infos.reserve(subchain.size());
  name = subchain.at(0).subchain;
  if (ent) {
    entity_id = ent->name;
    polymer = ent->entity_type == EntityType::Polymer;
    polymer_type = ent->polymer_type;
  } else {
    polymer = false;
    polymer_type = PolymerType::Unknown;
  }
  for (Residue& res : subchain)
    res_infos.emplace_back(&res);
}

inline void Topo::ChainInfo::setup_polymer_links() {
  if (!polymer || res_infos.empty())
    return;
  //ResidueGroup residue_groups()
  RGroup prev_group = group_from(res_infos.begin());
  while (prev_group.end != res_infos.end()) {
    RGroup group = group_from(prev_group.end);
    for (auto ri = group.begin; ri != group.end; ++ri)
      for (auto prev_ri = prev_group.begin; prev_ri != prev_group.end; ++prev_ri) {
        ResInfo::Prev p{{}, 0};
        if (are_connected(*prev_ri->res, *ri->res, polymer_type)) {
          p.idx = int(prev_ri - ri);
          if (is_polypeptide(polymer_type)) {
            if (ri->chemcomp.group == "P-peptide")
              p.link = "P";  // PCIS, PTRANS
            else if (ri->chemcomp.group == "M-peptide")
              p.link = "NM"; // NMCIS, NMTRANS
            p.link += prev_ri->res->is_cis ? "CIS" : "TRANS";
          } else if (is_polynucleotide(polymer_type)) {
            p.link = "p";
          } else {
            p.link = "?";
          }
        } else {
          p.link = "gap";
        }
        ri->prev.push_back(p);
      }
    prev_group = group;
  }
}

inline void Topo::ChainInfo::add_refmac_builtin_modifications() {
  if (polymer && !res_infos.empty()) {
    // we try to get exactly the same numbers that makecif produces
    for (Topo::ResInfo& ri : res_infos)
      if (polymer_type == PolymerType::PeptideL)
        ri.mods.emplace_back("AA-STAND");
    Topo::ResInfo& front = res_infos.front();
    Topo::ResInfo& back = res_infos.back();
    if (is_polypeptide(polymer_type)) {
      if (front.chemcomp.group == "P-peptide")
        front.mods.emplace_back("NH2");
      else
        front.mods.emplace_back("NH3");
      back.mods.emplace_back(back.res->find_atom("OXT", '*') ? "COO" : "TERMINUS");
    } else if (is_polynucleotide(polymer_type)) {
      front.mods.emplace_back("5*END");
      back.mods.emplace_back("TERMINUS");
    }
  }
}


// Model is non-const b/c we store non-const pointers to residues in Topo.
inline void Topo::initialize_refmac_topology(const Structure& st, Model& model0,
                                             MonLib& monlib) {
  // initialize chains and residues
  for (Chain& chain : model0.chains)
    for (ResidueSpan& sub : chain.subchains()) {
      const Entity* ent = st.get_entity_of(sub);
      chain_infos.emplace_back();
      chain_infos.back().initialize(sub, ent);
    }
  for (ChainInfo& ci : chain_infos) {
    // copy monomer description
    for (ResInfo& ri : ci.res_infos) {
      auto it = monlib.monomers.find(ri.res->name);
      if (it != monlib.monomers.end())
        ri.chemcomp = it->second;
    }

    ci.setup_polymer_links();

    ci.add_refmac_builtin_modifications();

    // add modifications from standard links
    for (ResInfo& ri : ci.res_infos)
      for (ResInfo::Prev& prev : ri.prev)
        if (const ChemLink* link = monlib.find_link(prev.link)) {
          prev.get(&ri)->add_mod(link->side1.mod);
          ri.add_mod(link->side2.mod);
        }
  }
  // add extra links
  for (const Connection& conn : st.connections) {
    // ignoring hydrogen bonds and metal coordination
    if (conn.type == Connection::Hydrog || conn.type == Connection::MetalC)
      continue;
    ExtraLink extra;
    extra.res1 = model0.find_cra(conn.partner1, true).residue;
    extra.res2 = model0.find_cra(conn.partner2, true).residue;
    if (!extra.res1 || !extra.res2)
      continue;
    extra.alt1 = conn.partner1.altloc;
    extra.alt2 = conn.partner2.altloc;
    const ChemLink* match =
        monlib.match_link(extra.res1->name, conn.partner1.atom_name,
                          extra.res2->name, conn.partner2.atom_name);
    if (!match) {
      match = monlib.match_link(extra.res2->name, conn.partner2.atom_name,
                                extra.res1->name, conn.partner1.atom_name);
      if (match) {
        std::swap(extra.res1, extra.res2);
        std::swap(extra.alt1, extra.alt2);
      }
    }
    if (match) {
      extra.link_id = match->id;
      // add modifications from the link
      find_resinfo(extra.res1)->add_mod(match->side1.mod);
      find_resinfo(extra.res2)->add_mod(match->side2.mod);
    } else {
      ChemLink cl;
      cl.side1.comp = extra.res1->name;
      cl.side2.comp = extra.res2->name;
      cl.id = cl.side1.comp + "-" + cl.side2.comp;
      Restraints::Bond bond;
      bond.id1 = Restraints::AtomId{1, conn.partner1.atom_name};
      bond.id2 = Restraints::AtomId{2, conn.partner2.atom_name};
      bond.type = BondType::Unspec;
      bond.aromatic = false;
      bond.value = conn.reported_distance;
      bond.esd = 0.02;
      cl.rt.bonds.push_back(bond);
      monlib.ensure_unique_link_name(cl.id);
      monlib.links.emplace(cl.id, cl);
      extra.link_id = cl.id;
    }
    extras.push_back(extra);
  }

  for (ChainInfo& chain_info : chain_infos)
    for (ResInfo& ri : chain_info.res_infos) {
      // apply modifications
      for (const std::string& modif : ri.mods) {
        if (const ChemMod* chem_mod = monlib.find_mod(modif))
          try {
            chem_mod->apply_to(ri.chemcomp);
          } catch(std::runtime_error& e) {
            printf("Failed to apply modification %s to %s: %s\n",
                   chem_mod->id.c_str(), ri.res->name.c_str(), e.what());
          }
        else
          printf("Modification not found: %s\n", modif.c_str());
      }
    }
}

} // namespace gemmi
#endif
