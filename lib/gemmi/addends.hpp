// Copyright 2020 Global Phasing Ltd.
//
// Addends to scattering form factors used in DensityCalculator
// and in StructureFactorCalculator.

#ifndef GEMMI_ADDENDS_HPP_
#define GEMMI_ADDENDS_HPP_

#include <array>
#include "elem.hpp"  // for El, Element

namespace gemmi {

struct Addends {
  std::array<float, (int)El::END> values = {};

  void set(Element el, float val) { values[el.ordinal()] = val; }
  float get(Element el) const { return values[el.ordinal()]; }
  size_t size() const { return values.size(); }
  void clear() {
    for (size_t i = 0; i != size(); ++i)
      values[i] = 0.;
  }
  void subtract_z(bool except_hydrogen=false) {
    for (int z = 2; z < (int)El::D; ++z)
      values[z] -= z;
    if (!except_hydrogen) {
      values[(int)El::H] -= 1;
      values[(int)El::D] -= 1;
    }
  }
};

} // namespace gemmi
#endif
