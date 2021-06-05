// Copyright 2019 Global Phasing Ltd.

// Calculation of atomic form factors approximated by a sum of Gaussians.
// Tables with numeric coefficient are in it92.hpp and c4322.hpp.

#ifndef GEMMI_FORMFACT_HPP_
#define GEMMI_FORMFACT_HPP_

#include <cmath>     // for exp, sqrt
#include <utility>   // for pair
#include "math.hpp"  // for pi()

namespace gemmi {

// precalculated density of an isotropic atom
template<int N, typename Real>
struct ExpSum {
  Real a[N], b[N];

  Real calculate(Real r2) const {
    Real density = 0;
    for (int i = 0; i < N; ++i)
      density += a[i] * std::exp(b[i] * r2);
    return density;
  }

  std::pair<Real,Real> calculate_with_derivative(Real r) const {
    Real density = 0;
    Real derivative = 0;
    for (int i = 0; i < N; ++i) {
      Real y = a[i] * std::exp(b[i] * r * r);
      density += y;
      derivative += 2 * b[i] * r * y;
    }
    return std::make_pair(density, derivative);
  }
};

// precalculated density of an anisotropic atom
template<int N, typename Real>
struct ExpAnisoSum {
  Real a[N];
  SMat33<Real> b[N];

  Real calculate(const Vec3& r) const {
    Real density = 0;
    for (int i = 0; i < N; ++i)
      density += a[i] * std::exp(b[i].r_u_r(r));
    return density;
  }
};

// Gaussian coefficients with functions to calculate sf and density.
template<int N, int WithC, typename Real>
struct GaussianCoef {
  using coef_type = Real;
  std::array<Real, 2*N+WithC> coefs;
  Real a(int n) const { return coefs[n]; }
  Real b(int n) const { return coefs[N+n]; }
  Real c() const { return WithC ? coefs[2*N] : 0; }

  void set_coefs(const std::array<Real, 2*N+WithC>& c) { coefs = c; }

  // argument: (sin(theta)/lambda)^2
  Real calculate_sf(Real stol2) const {
    Real sf = c();
    for (int i = 0; i < N; ++i)
      sf += a(i) * std::exp(-b(i)*stol2);
    return sf;
  }

  static constexpr Real pow15(Real x) { return x * std::sqrt(x); }

  Real calculate_density_iso(Real r2, Real B) const {
    constexpr Real _4pi = 4 * pi();
    Real r2pi = r2 * pi();
    Real density = c() * pow15(_4pi / B) * std::exp(-(_4pi / B) * r2pi);
    for (int i = 0; i < N; ++i) {
      Real t = _4pi / (b(i)+B);
      density += a(i) * pow15(t) * std::exp(-t*r2pi);
    }
    return density;
  }

  // note: addend is considered only if WithC (addend is usually dispersion f')
  ExpSum<N+WithC,Real> precalculate_density_iso(Real B, Real addend=0) const {
    ExpSum<N+WithC,Real> prec;
    constexpr Real _4pi = 4 * pi();
    for (int i = 0; i < N; ++i) {
      Real t = _4pi / (b(i)+B);
      prec.a[i] = a(i) * pow15(t);
      prec.b[i] = -t * pi();
    }
    if (WithC) {
      Real t = _4pi / B;
      prec.a[N] = (c() + addend) * pow15(t);
      prec.b[N] = -t * pi();
    }
    return prec;
  }

  Real calculate_density_aniso(const Vec3& r, const SMat33<float>& U) const {
    constexpr Real pi2 = sq(pi());
    const SMat33<Real> B = U.scaled(8 * pi2);
    Real density = c() * pow15(4 * pi()) / std::sqrt(B.determinant()) *
                   std::exp(-4 * pi2 * B.inverse().r_u_r(r));
    for (int i = 0; i < N; ++i) {
      SMat33<Real> Bb = B.added_kI(b(i));
      density += a(i) * pow15(4 * pi()) / std::sqrt(Bb.determinant()) *
                 std::exp(-4 * pi2 * Bb.inverse().r_u_r(r));
    }
    return density;
  }

  ExpAnisoSum<N+WithC,Real> precalculate_density_aniso_b(const SMat33<Real>& B,
                                                         Real addend=0) const {
    constexpr Real m4pi2 = -4 * sq(pi());
    ExpAnisoSum<N+WithC,Real> prec;
    for (int i = 0; i < N; ++i) {
      SMat33<Real> Bb = B.added_kI(b(i));
      prec.a[i] = a(i) * pow15(4 * pi()) / std::sqrt(Bb.determinant());
      prec.b[i] = Bb.inverse().scaled(m4pi2);
    }
    if (WithC) {
      prec.a[N] = (c() + addend) * pow15(4 * pi()) / std::sqrt(B.determinant());
      prec.b[N] = B.inverse().scaled(m4pi2);
    }
    return prec;
  }

  ExpAnisoSum<N+WithC,Real> precalculate_density_aniso_u(const SMat33<float>& U,
                                                         Real addend=0) const {
    constexpr Real UtoB = 8 * sq(pi());
    return precalculate_density_aniso_b(U.scaled(UtoB), addend);
  }
};

} // namespace gemmi
#endif
