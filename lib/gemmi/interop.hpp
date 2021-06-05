// Copyright 2020 Global Phasing Ltd.
//
// Interoperability between Model (MX) and SmallStructure (SX).

#ifndef GEMMI_INTEROP_HPP_
#define GEMMI_INTEROP_HPP_

#include "model.hpp"
#include "small.hpp"
#include "calculate.hpp"  // for count_atom_sites

namespace gemmi {

inline SmallStructure::Site atom_to_site(const Atom& atom, const UnitCell& cell) {
  SmallStructure::Site site;
  site.label = atom.name;
  site.type_symbol = atom.element.name();
  site.fract = cell.fractionalize(atom.pos);
  site.occ = atom.occ;
  // occupancy may need to be adjusted if the atom is on special position
  if (atom.occ <= 0.5) {
    int n_mates = cell.is_special_position(atom.pos);
    if (n_mates > 0 && atom.occ * (n_mates + 1) <= 1.0)
      site.occ = atom.occ * (n_mates + 1);
  }
  site.u_iso = atom.b_iso / u_to_b();
  if (atom.aniso.nonzero()) {
    if (cell.alpha == 90. || cell.beta == 90. || cell.gamma == 90.) {
      site.aniso = atom.aniso.scaled(1.0);
    } else {
      SMat33<double> t = atom.aniso.transformed_by<>(cell.frac.mat);
      Vec3 v = {1.0 / cell.ar, 1.0 / cell.br, 1.0 / cell.cr};
      site.aniso = {t.u11 * v.x * v.x,
                    t.u22 * v.y * v.y,
                    t.u33 * v.z * v.z,
                    t.u12 * v.x * v.y,
                    t.u13 * v.x * v.z,
                    t.u23 * v.y * v.z};
    }
  }
  site.element = atom.element;
  site.charge = atom.charge;
  return site;
}

inline SmallStructure mx_to_sx_structure(const Structure& st, int n=0) {
  const Model& model = st.models.at(n);
  SmallStructure small;
  small.name = st.name;
  small.cell = st.cell;
  small.spacegroup_hm = st.spacegroup_hm;
  small.sites.reserve(count_atom_sites(model));
  for (const Chain& chain : model.chains)
    for (const Residue& residue : chain.residues)
      for (const Atom& atom : residue.atoms)
        small.sites.push_back(atom_to_site(atom, st.cell));
  return small;
}

} // namespace gemmi
#endif
