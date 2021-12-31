// Copyright 2019 Global Phasing Ltd.
//
// Direct calculation of structure factors.
//
// It does not use optimizations described in the literature,
// cf. Bourhis et al (2014) https://doi.org/10.1107/S2053273314022207,
// because direct calculations are not used in MX if performance is important.
// For FFT-based calculations see dencalc.hpp + fourier.hpp.

#ifndef GEMMI_SFCALC_HPP_
#define GEMMI_SFCALC_HPP_

#include <complex>
#include "addends.hpp" // for Addends
#include "model.hpp"   // for Structure, ...
#include "small.hpp"   // for SmallStructure

namespace gemmi {

// calculate part of the structure factor: exp(2 pi i r * s)
inline std::complex<double> calculate_sf_part(const Fractional& fpos,
                                              const Miller& hkl) {
  double arg = 2 * pi() * (hkl[0]*fpos.x + hkl[1]*fpos.y + hkl[2]*fpos.z);
  return std::complex<double>{std::cos(arg), std::sin(arg)};
}

template <typename Table>
class StructureFactorCalculator {
public:
  StructureFactorCalculator(const UnitCell& cell) : cell_(cell) {}

  void set_stol2_and_scattering_factors(const Miller& hkl) {
    stol2_ = cell_.calculate_stol_sq(hkl);
    scattering_factors_.clear();
    scattering_factors_.resize(addends.size(), 0.);
  }

  double get_scattering_factor(Element element) {
    double& sfactor = scattering_factors_[element.ordinal()];
    if (sfactor == 0.) {
      if (!Table::has(element.elem))
        fail("Missing scattering factor for ", element.name());
      sfactor = Table::get(element.elem).calculate_sf(stol2_) + addends.get(element);
    }
    return sfactor;
  }

  // Calculation of Debye-Waller factor with isotropic ADPs
  double dwf_iso(const SmallStructure::Site& site) const {
    return std::exp(-u_to_b() * stol2_ * site.u_iso);
  }
  double dwf_iso(const Atom& atom) const {
    return std::exp(-stol2_ * atom.b_iso);
  }

  // Calculation of Debye-Waller factor exp(-2 pi^2 s.U.s)
  // cf. B. Rupp's book, p. 641 or RWGK & Adams 2002, J. Appl. Cryst. 35, 477
  // Small molecule and macromolecular anisotropic U's are defined differently,
  // so we have two functions.
  double dwf_aniso(const SmallStructure::Site& site, const Vec3& hkl) const {
    Vec3 arh(cell_.ar * hkl.x, cell_.br * hkl.y, cell_.cr * hkl.z);
    return std::exp(-2 * pi() * pi() * site.aniso.r_u_r(arh));
  }
  double dwf_aniso(const Atom& atom, const Vec3& hkl) const {
    return std::exp(-2 * pi() * pi() *
                    atom.aniso.transformed_by<>(cell_.frac.mat).r_u_r(hkl));
  }

  template<typename Site>
  std::complex<double> calculate_sf_from_atom_sf(const Fractional& fract,
                                                 const Site& site,
                                                 const Miller& hkl,
                                                 double sf) {
    double oc_sf = site.occ * sf;
    std::complex<double> sum = calculate_sf_part(fract, hkl);
    if (!site.aniso.nonzero()) {
      for (const FTransform& image : cell_.images)
        sum += calculate_sf_part(image.apply(fract), hkl);
      return oc_sf * dwf_iso(site) * sum;
    }
    Vec3 vhkl(hkl[0], hkl[1], hkl[2]);
    sum *= dwf_aniso(site, vhkl);
    for (const FTransform& image : cell_.images)
      sum += calculate_sf_part(image.apply(fract), hkl) *
             dwf_aniso(site, image.mat.left_multiply(vhkl));
    return oc_sf * sum;
  }

  template<typename Site>
  std::complex<double> calculate_sf_from_atom(const Fractional& fract,
                                              const Site& site,
                                              const Miller& hkl) {
    return calculate_sf_from_atom_sf(fract, site, hkl, get_scattering_factor(site.element));
  }

  std::complex<double> calculate_sf_from_model(const Model& model, const Miller& hkl) {
    std::complex<double> sf = 0.;
    set_stol2_and_scattering_factors(hkl);
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& site : res.atoms)
          sf += calculate_sf_from_atom(cell_.fractionalize(site.pos), site, hkl);
    return sf;
  }

  // Z part of Mott-Bethe formula (when need to use different model)
  std::complex<double> calculate_mb_z(const Model& model, const Miller& hkl, bool only_h) {
    std::complex<double> sf = 0.;
    stol2_ = cell_.calculate_stol_sq(hkl);
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& site : res.atoms)
          if (!only_h || site.element.is_hydrogen())
            sf += calculate_sf_from_atom_sf(cell_.fractionalize(site.pos), site, hkl, -1.*site.element.atomic_number());
    return sf;
  }

  double mott_bethe_factor() const {
    return -1. / (8 * pi() * pi() * bohrradius()) / stol2_;
  }

  // The occupancy is assumed to take into account symmetry,
  // i.e. to be fractional if the atom is on special position.
  std::complex<double> calculate_sf_from_small_structure(const SmallStructure& small,
                                                         const Miller& hkl) {
    std::complex<double> sf = 0.;
    set_stol2_and_scattering_factors(hkl);
    for (const SmallStructure::Site& site : small.sites)
      sf += calculate_sf_from_atom(site.fract, site, hkl);
    return sf;
  }

private:
  const UnitCell& cell_;
  double stol2_;
  std::vector<double> scattering_factors_;
public:
  Addends addends;  // usually f' for X-rays
};

} // namespace gemmi
#endif
