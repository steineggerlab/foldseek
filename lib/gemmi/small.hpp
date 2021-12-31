// Copyright 2018 Global Phasing Ltd.
//
// Representation of small molecule or inorganic crystal.
// Flat list of atom sites. Minimal functionality.

#ifndef GEMMI_SMALL_HPP_
#define GEMMI_SMALL_HPP_

#include <cctype>        // for isalpha
#include <algorithm>     // for any_of
#include <bitset>
#include <string>
#include <vector>
#include "elem.hpp"      // Element
#include "math.hpp"      // SMat33
#include "symmetry.hpp"  // find_spacegroup_by_name
#include "unitcell.hpp"  // UnitCell, Fractional
#include "util.hpp"      // vector_remove_if

namespace gemmi {

struct SmallStructure {
  struct Site {
    std::string label;
    std::string type_symbol;
    Fractional fract;
    double occ = 1.0;
    double u_iso = 0.;
    SMat33<double> aniso = {0, 0, 0, 0, 0, 0};
    int disorder_group = 0;
    Element element = El::X;
    signed char charge = 0;  // [-8, +8]

    Position orth(const gemmi::UnitCell& cell_) const {
      return cell_.orthogonalize(fract);
    }

    std::string element_and_charge_symbol() const {
      std::string s = element.name();
      if (charge != 0) {
        s += std::to_string(std::abs(charge));
        s += charge > 0 ? '+' : '-';
      }
      return s;
    }
  };

  struct AtomType {
    std::string symbol;
    Element element = El::X;
    signed char charge = 0;  // [-8, +8]
    double dispersion_real;
    double dispersion_imag;
  };

  std::string name;
  UnitCell cell;
  std::string spacegroup_hm;
  std::vector<Site> sites;
  std::vector<AtomType> atom_types;
  double wavelength = 0.; // the first wavelength if multiple

  std::vector<Site> get_all_unit_cell_sites() const;

  const SpaceGroup* find_spacegroup() const {
    return find_spacegroup_by_name(spacegroup_hm, cell.alpha, cell.gamma);
  }

  const AtomType* get_atom_type(const std::string& symbol) const {
    for (const AtomType& at : atom_types)
      if (at.symbol == symbol)
        return &at;
    return nullptr;
  }

  // similar to Model::present_elements() from model.hpp
  std::bitset<(size_t)El::END> present_elements() const {
    std::bitset<(size_t)El::END> table;
    for (const Site& atom : sites)
      table.set((size_t)atom.element.elem);
    return table;
  }

  void remove_hydrogens() {
    vector_remove_if(sites, [](const Site& a) { return a.element.is_hydrogen(); });
  }

  // pre: atoms on special positions have "chemical" occupancy (i.e. not divided
  // by n for n-fold symmetry)
  void change_occupancies_to_crystallographic(double max_dist=0.4) {
    for (Site& site : sites) {
      int n_mates = cell.is_special_position(site.fract, max_dist);
      if (n_mates != 0)
        site.occ /= (n_mates + 1);
    }
  }

  void setup_cell_images() {
    cell.set_cell_images_from_spacegroup(find_spacegroup());
  }
};

template<typename T>
inline void split_element_and_charge(const std::string& label, T* dest) {
  int len = label.size() > 1 && std::isalpha(label[1]) ? 2 : 1;
  dest->element = len == 1 ? impl::find_single_letter_element(label[0] & ~0x20)
                           : find_element(label.c_str());
  if (dest->element != El::X && (label.back() == '+' || label.back() == '-')) {
    int sign = label.back() == '+' ? 1 : -1;
    if (label.size() - len == 1)
      dest->charge = sign;
    else if (label.size() - len == 2 && label[len] >= '0' && label[len] <= '9')
      dest->charge = sign * (label[len] - '0');
  }
}

inline std::vector<SmallStructure::Site>
SmallStructure::get_all_unit_cell_sites() const {
  const double SPECIAL_POS_TOL = 0.4;
  std::vector<Site> all;
  for (const Site& site : sites) {
    size_t start = all.size();
    all.push_back(site);
    for (const FTransform& image : cell.images) {
      Fractional fpos = image.apply(site.fract);
      if (std::any_of(all.begin() + start, all.end(), [&](const Site& other) {
            return cell.distance_sq(fpos, other.fract) < sq(SPECIAL_POS_TOL);
          }))
        continue;
      all.push_back(site);
      all.back().fract = fpos;
    }
  }
  return all;
}

} // namespace gemmi
#endif
