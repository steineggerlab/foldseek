// Copyright 2019 Global Phasing Ltd.
//
// Function used in both mmcif.hpp and refln.hpp (for coordinate and
// reflection mmCIF files).

#ifndef GEMMI_MMCIF_IMPL_HPP_
#define GEMMI_MMCIF_IMPL_HPP_

#include "cifdoc.hpp"    // for cif::Block
#include "numb.hpp"      // for cif::as_number
#include "unitcell.hpp"  // for UnitCell
#include "symmetry.hpp"  // for SpaceGroup

namespace gemmi {
namespace impl {

inline void set_cell_from_mmcif(cif::Block& block, UnitCell& cell,
                                bool mmcif=true) {
  cif::Table tab = block.find((mmcif ? "_cell." : "_cell_"),
                              {"length_a", "length_b", "length_c",
                               "angle_alpha", "angle_beta", "angle_gamma"});
  if (tab.ok()) {
    auto c = tab.one();
    using namespace gemmi::cif;
    if (!is_null(c[0]) && !is_null(c[1]) && !is_null(c[2]))
      cell.set(as_number(c[0]), as_number(c[1]), as_number(c[2]),
               as_number(c[3]), as_number(c[4]), as_number(c[5]));
  }
}

inline const std::string* find_spacegroup_hm_value(const cif::Block& block) {
  const char* hm_tag = "_symmetry.space_group_name_H-M";
  return block.find_value(hm_tag);
}

} // namespace impl
} // namespace gemmi

#endif
