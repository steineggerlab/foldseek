// Copyright 2017-2021 Global Phasing Ltd.
//
// Modify various properties of the model.

// For modifications that depend on entities or connectivity see polyheur.hpp.

#ifndef GEMMI_MODIFY_HPP_
#define GEMMI_MODIFY_HPP_

#include "model.hpp"
#include "util.hpp"      // for vector_remove_if
#include <set>

namespace gemmi {

// Remove alternative conformations.
template<class T> void remove_alternative_conformations(T& obj) {
  for (auto& child : obj.children())
    remove_alternative_conformations(child);
}
template<> inline void remove_alternative_conformations(Chain& chain) {
  std::set<SeqId> seqids;
  for (size_t i = 0; i < chain.residues.size(); ) {
    if (seqids.insert(chain.residues[i].seqid).second)
      ++i;
    else
      chain.residues.erase(chain.residues.begin() + i);
  }
  for (Residue& residue : chain.residues) {
    std::set<std::string> names;
    for (size_t i = 0; i < residue.atoms.size(); ) {
      Atom& atom = residue.atoms[i];
      atom.altloc = '\0';
      if (names.insert(atom.name).second)
        ++i;
      else
        residue.atoms.erase(residue.atoms.begin() + i);
    }
  }
}

// Remove hydrogens.
template<class T> void remove_hydrogens(T& obj) {
  for (auto& child : obj.children())
    remove_hydrogens(child);
}
template<> inline void remove_hydrogens(Residue& res) {
  vector_remove_if(res.atoms, [](const Atom& a) {
    return a.element == El::H || a.element == El::D;
  });
}

// Remove anisotropic ADP
template<class T> void remove_anisou(T& obj) {
  for (auto& child : obj.children())
    remove_anisou(child);
}
template<> inline void remove_anisou(Atom& atom) {
  atom.aniso = {0, 0, 0, 0, 0, 0};
}

// Set absent ANISOU to value from B_iso
template<class T> void ensure_anisou(T& obj) {
  for (auto& child : obj.children())
    ensure_anisou(child);
}
template<> inline void ensure_anisou(Atom& atom) {
  if (!atom.aniso.nonzero()) {
    float u = float(1. / gemmi::u_to_b() * atom.b_iso);
    atom.aniso = {u, u, u, 0.f, 0.f, 0.f};
  }
}

// apply Transform to both atom's position and ADP
template<class T> void transform_pos_and_adp(T& obj, const Transform& tr) {
  for (auto& child : obj.children())
    transform_pos_and_adp(child, tr);
}
template<> inline void transform_pos_and_adp(Atom& atom, const Transform& tr) {
  atom.pos = Position(tr.apply(atom.pos));
  if (atom.aniso.nonzero())
    atom.aniso = atom.aniso.transformed_by<float>(tr.mat);
}

} // namespace gemmi
#endif
