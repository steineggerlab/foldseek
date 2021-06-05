// Copyright 2018 Global Phasing Ltd.
//
// Read small molecule CIF file into SmallStructure (from small.hpp).

#ifndef GEMMI_SMCIF_HPP_
#define GEMMI_SMCIF_HPP_

#include "small.hpp"     // SmallStructure
#include "cifdoc.hpp"
#include "numb.hpp"      // for as_number
#include "symmetry.hpp"  // SpaceGroup

namespace gemmi {

inline
SmallStructure make_small_structure_from_block(const cif::Block& block_) {
  using cif::as_number;
  using cif::as_string;
  cif::Block& block = const_cast<cif::Block&>(block_);
  SmallStructure st;
  st.name = block.name;

  // unit cell and symmetry
  cif::Table cell = block.find("_cell_",
                               {"length_a", "length_b", "length_c",
                                "angle_alpha", "angle_beta", "angle_gamma"});
  if (cell.ok()) {
    auto c = cell.one();
    if (!cif::is_null(c[0]) && !cif::is_null(c[1]) && !cif::is_null(c[2]))
      st.cell.set(as_number(c[0]), as_number(c[1]), as_number(c[2]),
                  as_number(c[3]), as_number(c[4]), as_number(c[5]));
  }
  for (const char* tag : {"_space_group_name_H-M_alt",
                          "_symmetry_space_group_name_H-M"})
    if (const std::string* val = block.find_value(tag)) {
      st.spacegroup_hm = as_string(*val);
      break;
    }

  enum { kLabel, kSymbol, kX, kY, kZ, kUiso, kBiso, kOcc, kDisorderGroup };
  cif::Table atom_table = block.find("_atom_site_",
                                     {"label",
                                      "?type_symbol",
                                      "?fract_x",
                                      "?fract_y",
                                      "?fract_z",
                                      "?U_iso_or_equiv",
                                      "?B_iso_or_equiv",
                                      "?occupancy",
                                      "?disorder_group"});
  for (auto row : atom_table) {
    SmallStructure::Site site;
    site.label = as_string(row[kLabel]);
    if (row.has(kSymbol))
      site.type_symbol = as_string(row[kSymbol]);
    else
      site.type_symbol = site.label;
    if (row.has(kX))
      site.fract.x = as_number(row[kX]);
    if (row.has(kY))
      site.fract.y = as_number(row[kY]);
    if (row.has(kZ))
      site.fract.z = as_number(row[kZ]);
    if (row.has(kUiso))
      site.u_iso = as_number(row[kUiso], 0.0);
    else if (row.has(kBiso))
      site.u_iso = as_number(row[kBiso], 0.0) / u_to_b();
    if (row.has(kOcc))
      site.occ = as_number(row[kOcc], 1.0);
    if (row.has2(kDisorderGroup))
      // ignore non-integer _atom_site_disorder_group
      site.disorder_group = string_to_int(row[kDisorderGroup], false);
    split_element_and_charge(site.type_symbol, &site);
    st.sites.push_back(site);
  }

  std::vector<SmallStructure::Site>::iterator aniso_site = st.sites.begin();
  for (auto row : block.find("_atom_site_aniso_",
                             {"label", "U_11", "U_22", "U_33",
                              "U_12", "U_13", "U_23"})) {
    std::string aniso_label = row.str(0);
    if (aniso_site == st.sites.end() || aniso_site->label != aniso_label) {
      aniso_site = std::find_if(st.sites.begin(), st.sites.end(),
              [&](SmallStructure::Site& x) { return x.label == aniso_label; });
      if (aniso_site == st.sites.end())
        continue;
    }
    aniso_site->aniso.u11 = as_number(row[1], 0.0);
    aniso_site->aniso.u22 = as_number(row[2], 0.0);
    aniso_site->aniso.u33 = as_number(row[3], 0.0);
    aniso_site->aniso.u12 = as_number(row[4], 0.0);
    aniso_site->aniso.u13 = as_number(row[5], 0.0);
    aniso_site->aniso.u23 = as_number(row[6], 0.0);
    ++aniso_site;
  }

  for (auto row : block.find("_atom_type_",
                             {"symbol",
                              "scat_dispersion_real",
                              "scat_dispersion_imag"})) {
    SmallStructure::AtomType atom_type;
    atom_type.symbol = row.str(0);
    atom_type.dispersion_real = cif::as_number(row[1]);
    atom_type.dispersion_imag = cif::as_number(row[2]);
    split_element_and_charge(atom_type.symbol, &atom_type);
    st.atom_types.push_back(atom_type);
  }
  if (cif::Column w_col = block.find_values("_diffrn_radiation_wavelength"))
    st.wavelength = cif::as_number(w_col.at(0));

  const SpaceGroup* sg = find_spacegroup_by_name(st.spacegroup_hm,
                                                 st.cell.alpha, st.cell.gamma);
  st.cell.set_cell_images_from_spacegroup(sg);
  return st;
}

} // namespace gemmi
#endif
