// Ratio of modified Bessel functions I1(x)/I0(x).
//
// When gemmi requires C++17 we might use std::cyl_bessel_if.
// Crystallographic codes (including Refmac and cctbx) often use polynomial
// approximation of I0 and I1 from p. 378 of Abramowitz and Stegun.
// Gemmi uses approximation based on polynomial coefficients from bessel_i0
// and bessel_i1(float) from Boost.Math:
// https://www.boost.org/doc/libs/1_76_0/libs/math/doc/html/math_toolkit/bessel/mbessel.html
// This approximation was derived in 2017 by John Maddock,
// building on the work of Pavel Holoborodko:
// https://www.advanpix.com/2015/11/11/rational-approximations-for-the-modified-bessel-function-of-the-first-kind-i0-computations-double-precision/
// The efficiency is similar to that of scitbx.math.bessel_i1_over_i0.

#ifndef GEMMI_BESSEL_HPP_
#define GEMMI_BESSEL_HPP_

#include <cmath>

namespace gemmi {

template<int N>
inline float evaluate_polynomial(const float(&poly)[N], float x) {
  static_assert(N > 1, "");
  float result = poly[N-1];
  for (int i = N-2; i >= 0; --i)
    result = result * x + poly[i];
  return result;
}

inline float bessel_i1_over_i0(float x) {
  static const float P1[] = {
    8.333333221e-02f,
    6.944453712e-03f,
    3.472097211e-04f,
    1.158047174e-05f,
    2.739745142e-07f,
    5.135884609e-09f,
    5.262251502e-11f,
    1.331933703e-12f
  };
  static const float Q1[] = {
    1.00000003928615375e+00f,
    2.49999576572179639e-01f,
    2.77785268558399407e-02f,
    1.73560257755821695e-03f,
    6.96166518788906424e-05f,
    1.89645733877137904e-06f,
    4.29455004657565361e-08f,
    3.90565476357034480e-10f,
    1.48095934745267240e-11f
  };
  static const float P2[] = {
    3.98942115977513013e-01f,
    -1.49581264836620262e-01f,
    -4.76475741878486795e-02f,
    -2.65157315524784407e-02f,
    -1.47148600683672014e-01f
  };
  static const float Q2[] = {
    3.98942651588301770e-01f,
    4.98327234176892844e-02f,
    2.91866904423115499e-02f,
    1.35614940793742178e-02f,
    1.31409251787866793e-01f
  };
  static const float Q3[] = {
    3.98942391532752700e-01f,
    4.98455950638200020e-02f,
    2.94835666900682535e-02f
  };

  if (x < 0)
    return -bessel_i1_over_i0(-x);

  if (x < 7.75) {
     float a = x * x / 4;
     float bessel0 = a * evaluate_polynomial(Q1, a) + 1;
     float R[3] = { 1, 0.5f, evaluate_polynomial(P1, a) };
     float bessel1 = x * evaluate_polynomial(R, a) / 2;
     return bessel1 / bessel0;
  }

  float p = evaluate_polynomial(P2, 1 / x);
  float q = x < 50 ? evaluate_polynomial(Q2, 1 / x)
                   : evaluate_polynomial(Q3, 1 / x);
  return p / q;
}

} // namespace gemmi
#endif
