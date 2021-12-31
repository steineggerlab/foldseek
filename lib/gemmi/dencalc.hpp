// Copyright 2019 Global Phasing Ltd.
//
// Tools to prepare a grid with values of electron density of a model.

#ifndef GEMMI_DENCALC_HPP_
#define GEMMI_DENCALC_HPP_

#include <cassert>
#include "addends.hpp"  // for Addends
#include "formfact.hpp" // for ExpSum
#include "grid.hpp"     // for Grid
#include "model.hpp"    // for Structure, ...

namespace gemmi {

template <typename PrecalExpSum>
double determine_cutoff_radius(double x1, const PrecalExpSum& precal, double cutoff_level) {
  double y1, dy;
  std::tie(y1, dy) = precal.calculate_with_derivative(x1);
  // Generally, density is supposed to decrease with radius.
  // But if we have addends (in particular -Z for Mott-Bothe),
  // it can first rise, then decrease. We want to be after the maximum.
  while (dy > 0) { // unlikely
    x1 += 1.0;
    std::tie(y1, dy) = precal.calculate_with_derivative(x1);
  }
  double x2 = x1;
  double y2 = y1;
  if (y1 < cutoff_level) {
    while (y1 < cutoff_level) {
      x2 = x1;
      y2 = y1;
      x1 -= 0.5f;
      std::tie(y1, dy) = precal.calculate_with_derivative(x1);
      // with addends it's possible to land on the left side of the maximum
      if (dy > 0) { // unlikely
        while (dy > 0 && x1 + 0.1 < x2) {
          x1 += 0.1;
          std::tie(y1, dy) = precal.calculate_with_derivative(x1);
        }
        if (y1 < cutoff_level)
          return x1;
        break;
      }
      if (x1 < 0) { // unlikely
        x1 = 0;
        y1 = precal.calculate(x1 * x1);
        break;
      }
    }
  } else {
    while (y2 > cutoff_level) {
      x1 = x2;
      y1 = y2;
      x2 += 0.5f;
      y2 = precal.calculate(x2 * x2);
    }
  }

  return x1 + (x1 - x2) / (y1 - y2) * (cutoff_level - y1);
}

// approximated radius of electron density (IT92) above cutoff=1e-5 for C
inline double it92_radius_approx(double b) {
  return (8.5 + 0.075 * b) / (2.4 + 0.0045 * b);
}

inline double get_minimum_b_iso(const Model& model) {
  double b_min = 1000.;
  for (const Chain& chain : model.chains)
    for (const Residue& residue : chain.residues)
      for (const Atom& atom : residue.atoms)
        if (atom.b_iso < b_min)
          b_min = atom.b_iso;
  return b_min;
}

// Usual usage:
// - set d_min and optionally also other parameters,
// - set addends to f' values for your wavelength (see fprime.hpp)
// - use set_grid_cell_and_spacegroup() to set grid's unit cell and space group
// - check that Table has SF coefficients for all elements that are to be used
// - call put_model_density_on_grid()
// - do FFT using transform_map_to_f_phi()
// - if blur is used, multiply the SF by reciprocal_space_multiplier()
template <typename Table, typename Real>
struct DensityCalculator {
  Grid<Real> grid;
  double d_min = 0.;
  double rate = 1.5;
  double blur = 0.;
  float cutoff = 1e-5f;
  Addends addends;

  using coef_type = typename Table::Coef::coef_type;

  double requested_grid_spacing() const { return d_min / (2 * rate); }

  void set_refmac_compatible_blur(const Model& model) {
    double b_min = get_minimum_b_iso(model);
    blur = std::max(u_to_b() / 1.1 * sq(requested_grid_spacing()) - b_min, 0.);
  }

  // pre: check if Table::has(atom.element)
  void add_atom_density_to_grid(const Atom& atom) {
    Element el = atom.element;
    do_add_atom_density_to_grid(atom, Table::get(el), addends.get(el));
  }

  // Parameter c is a constant factor and has the same meaning as either addend
  // or c in scattering factor coefficients (a1, b1, ..., c).
  void add_c_contribution_to_grid(const Atom& atom, float c) {
    do_add_atom_density_to_grid(atom, GaussianCoef<0, 1, coef_type>{0}, c);
  }

  template<int N>
  double estimate_radius(const ExpSum<N, coef_type>& precal, double b) const {
    if (N == 1)
      return std::sqrt(std::log(cutoff / std::abs(precal.a[0])) / precal.b[0]);
    double x1 = it92_radius_approx(b);
    return determine_cutoff_radius(x1, precal, cutoff);
  }

  template<typename Coef>
  void do_add_atom_density_to_grid(const Atom& atom, const Coef& coef, float addend) {
    Fractional fpos = grid.unit_cell.fractionalize(atom.pos);
    if (!atom.aniso.nonzero()) {
      // isotropic
      double b = atom.b_iso + blur;
      auto precal = coef.precalculate_density_iso(b, addend);
      double radius = estimate_radius(precal, b);
      grid.template use_points_around<true>(fpos, radius, [&](Real& point, double r2) {
          point += Real(atom.occ * precal.calculate((Real)r2));
      }, /*fail_on_too_large_radius=*/false);
    } else {
      // anisotropic
      SMat33<double> aniso_b = atom.aniso.scaled(u_to_b()).added_kI(blur);
      // rough estimate, so we don't calculate eigenvalues
      double b_max = std::max(std::max(aniso_b.u11, aniso_b.u22), aniso_b.u33);
      auto precal_iso = coef.precalculate_density_iso(b_max, addend);
      double radius = estimate_radius(precal_iso, b_max);
      auto precal = coef.precalculate_density_aniso_b(aniso_b, addend);
      int du = (int) std::ceil(radius / grid.spacing[0]);
      int dv = (int) std::ceil(radius / grid.spacing[1]);
      int dw = (int) std::ceil(radius / grid.spacing[2]);
      grid.template use_points_in_box<true>(fpos, du, dv, dw,
                             [&](Real& point, const Position& delta) {
        if (delta.length_sq() < radius * radius)
          point += Real(atom.occ * precal.calculate(delta));
      }, false);
    }
  }

  void initialize_grid() {
    grid.data.clear();
    grid.set_size_from_spacing(requested_grid_spacing(), true);
  }

  void add_model_density_to_grid(const Model& model) {
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& atom : res.atoms)
          add_atom_density_to_grid(atom);
  }

  void put_model_density_on_grid(const Model& model) {
    initialize_grid();
    add_model_density_to_grid(model);
    grid.symmetrize_sum();
  }

  void set_grid_cell_and_spacegroup(const Structure& st) {
    grid.unit_cell = st.cell;
    grid.spacegroup = st.find_spacegroup();
  }

  // The argument is 1/d^2 - as outputted by unit_cell.calculate_1_d2(hkl).
  double reciprocal_space_multiplier(double inv_d2) const {
    return std::exp(blur * 0.25 * inv_d2);
  }

  double mott_bethe_factor(const Miller& hkl) const {
    double inv_d2 = grid.unit_cell.calculate_1_d2(hkl);
    double factor = -1. / (2 * pi() * pi() * bohrradius()) / inv_d2;
    return blur == 0 ? factor : factor * reciprocal_space_multiplier(inv_d2);
  }
};

} // namespace gemmi
#endif
