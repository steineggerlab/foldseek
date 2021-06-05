// Copyright 2020 Global Phasing Ltd.
//
// Reciprocal space utilities.

#ifndef GEMMI_RECIPROC_HPP_
#define GEMMI_RECIPROC_HPP_

#include <array>
#include "fail.hpp"      // for fail
#include "symmetry.hpp"  // for SpaceGroup
#include "unitcell.hpp"  // for UnitCell

namespace gemmi {

// dmin should include a tiny margin for numerical errors
template<typename Func>
void for_all_reflections(Func func,
                         const UnitCell& cell, const SpaceGroup* spacegroup,
                         double dmin, double dmax=0., bool unique=true) {
  Miller lim = cell.get_hkl_limits(dmin);
  double inv_dmin2 = 1. / sq(dmin);
  double inv_dmax2 = dmax > 0 ? 1. / sq(dmax) : 0;
  ReciprocalAsu asu(spacegroup);
  GroupOps gops = spacegroup->operations();
  Miller hkl;
  for (hkl[0] = -lim[0]; hkl[0] <= lim[0]; ++hkl[0])
    for (hkl[1] = -lim[1]; hkl[1] <= lim[1]; ++hkl[1])
      for (hkl[2] = -lim[2]; hkl[2] <= lim[2]; ++hkl[2])
        if (!unique || asu.is_in(hkl)) {
          double inv_d2 = cell.calculate_1_d2(hkl);
          if (inv_d2 <= inv_dmin2 && inv_d2 > inv_dmax2 &&
              !gops.is_systematically_absent(hkl))
            func(hkl);
        }
}

// dmin should include a tiny margin for numerical errors
inline int count_reflections(const UnitCell& cell, const SpaceGroup* spacegroup,
                             double dmin, double dmax=0., bool unique=true) {
  int counter = 0;
  for_all_reflections([&counter](const Miller&) { ++counter; },
                      cell, spacegroup, dmin, dmax, unique);
  return counter;
}

inline std::vector<Miller>
make_miller_vector(const UnitCell& cell, const SpaceGroup* spacegroup,
                   double dmin, double dmax=0., bool unique=true) {
  std::vector<Miller> hkls;
  for_all_reflections([&hkls](const Miller& hkl) { hkls.push_back(hkl); },
                      cell, spacegroup, dmin, dmax, unique);
  return hkls;
}


} // namespace gemmi
#endif
