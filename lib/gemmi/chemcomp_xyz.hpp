// Copyright 2018 Global Phasing Ltd.
//
// Reading coordinates from chemical component or Refmac monomer library files.

#ifndef GEMMI_CHEMCOMP_XYZ_HPP_
#define GEMMI_CHEMCOMP_XYZ_HPP_

#include <array>
#include "cifdoc.hpp" // for Block, etc
#include "fail.hpp"   // for fail
#include "numb.hpp"   // for as_number
#include "model.hpp"  // for Atom, Residue, etc

namespace gemmi {

// Reading chemical component as a coordinate file.
enum class ChemCompModel {
  Xyz,     // _chem_comp_atom.x, etc
  Example, // _chem_comp_atom.model_Cartn_x
  Ideal    // _chem_comp_atom.pdbx_model_Cartn_x_ideal
};

inline Residue make_residue_from_chemcomp_block(const cif::Block& block,
                                                ChemCompModel kind) {
  std::array<std::string, 3> xyz_tags;
  switch (kind) {
    case ChemCompModel::Xyz:
      xyz_tags = {{"x", "y", "z"}};
      break;
    case ChemCompModel::Example:
      xyz_tags = {{"model_Cartn_x", "model_Cartn_y", "model_Cartn_z"}};
      break;
    case ChemCompModel::Ideal:
      xyz_tags = {{"pdbx_model_Cartn_x_ideal",
                   "pdbx_model_Cartn_y_ideal",
                   "pdbx_model_Cartn_z_ideal"}};
      break;
  }
  Residue res;
  cif::Column col =
    const_cast<cif::Block&>(block).find_values("_chem_comp_atom.comp_id");
  if (col && col.length() > 0)
    res.name = col[0];
  else
    res.name = block.name.substr(starts_with(block.name, "comp_") ? 5 : 0);
  cif::Table table = const_cast<cif::Block&>(block).find("_chem_comp_atom.",
          {"atom_id", "type_symbol", "?charge",
           xyz_tags[0], xyz_tags[1], xyz_tags[2]});
  res.atoms.resize(table.length());
  int n = 0;
  for (auto row : table) {
    Atom& atom = res.atoms[n++];
    atom.name = row.str(0);
    atom.element = Element(row.str(1));
    if (row.has2(2))
      // Charge is defined as integer, but some cif files in the wild have
      // trailing '.000', so we read it as floating-point number.
      atom.charge = (signed char) std::round(cif::as_number(row[2]));
    atom.pos = Position(cif::as_number(row[3]),
                        cif::as_number(row[4]),
                        cif::as_number(row[5]));
  }
  return res;
}

inline Model make_model_from_chemcomp_block(const cif::Block& block,
                                            ChemCompModel kind) {
  std::string name;
  switch (kind) {
    case ChemCompModel::Xyz: name = "xyz"; break;
    case ChemCompModel::Example: name = "example_xyz"; break;
    case ChemCompModel::Ideal: name = "ideal_xyz"; break;
  }
  Model model(name);
  model.chains.emplace_back("");
  model.chains[0].residues.push_back(
      make_residue_from_chemcomp_block(block, kind));
  return model;
}

// For CCD input - returns a structure with two single-residue models:
// example (model_Cartn_x) and ideal (pdbx_model_Cartn_x_ideal).
// For Refmac dictionary (monomer library) files returns structure with
// a single model.
inline Structure make_structure_from_chemcomp_block(const cif::Block& block) {
  Structure st;
  st.input_format = CoorFormat::ChemComp;
  if (const std::string* name = block.find_value("_chem_comp.id"))
    st.name = *name;
  if (block.has_any_value("_chem_comp_atom.x"))
    st.models.push_back(
        make_model_from_chemcomp_block(block, ChemCompModel::Xyz));
  if (block.has_any_value("_chem_comp_atom.model_Cartn_x"))
    st.models.push_back(
        make_model_from_chemcomp_block(block, ChemCompModel::Example));
  if (block.has_any_value("_chem_comp_atom.pdbx_model_Cartn_x_ideal"))
    st.models.push_back(
        make_model_from_chemcomp_block(block, ChemCompModel::Ideal));
  return st;
}

// a helper function for use with make_structure_from_chemcomp_block():
//   int n = check_chemcomp_block_number(doc);
//   if (n != -1)
//     Structure st = make_structure_from_chemcomp_block(doc.blocks[n]);
inline int check_chemcomp_block_number(const cif::Document& doc) {
  // monomer library file without global_
  if (doc.blocks.size() == 2 && doc.blocks[0].name == "comp_list")
    return 1;
  // monomer library file with global_
  if (doc.blocks.size() == 3 && doc.blocks[0].name.empty() &&
      doc.blocks[1].name == "comp_list")
    return 2;
  // CCD file
  if (doc.blocks.size() == 1 &&
      !doc.blocks[0].has_tag("_atom_site.id") &&
      doc.blocks[0].has_tag("_chem_comp_atom.atom_id"))
    return 0;
  return -1;
}

inline Structure make_structure_from_chemcomp_doc(const cif::Document& doc) {
  int n = check_chemcomp_block_number(doc);
  if (n == -1)
    fail("Not a chem_comp format.");
  return make_structure_from_chemcomp_block(doc.blocks[n]);
}

} // namespace gemmi
#endif
