// Copyright 2021 Global Phasing Ltd.
//
// Unit cell reductions: Buerger, Niggli, Selling-Delaunay.

#ifndef GEMMI_CELLRED_HPP_
#define GEMMI_CELLRED_HPP_

#include <cmath>
#include <array>
#include "math.hpp"  // for deg

namespace gemmi {

struct SellingVector;

// GruberVector contains G6 vector (G for Gruber) and cell reduction algorithms.
// Originally, in B. Gruber, Acta Cryst. A29, 433 (1973), the vector was called
// "characteristic" of a lattice/cell.
// Functions that take epsilon as a parameter use it for comparisons,
// as proposed in Grosse-Kunstleve et al, Acta Cryst. (2004) A60, 1.
struct GruberVector {
  //    a.a  b.b c.c 2b.c 2a.c 2a.b
  double A, B, C, xi, eta, zeta;  // the 1973 paper uses names A B C ξ η ζ

  // m - orthogonalization matrix of a primitive cell
  explicit GruberVector(const Mat33& m)
    : A(m.column_dot(0,0)),
      B(m.column_dot(1,1)),
      C(m.column_dot(2,2)),
      xi(2 * m.column_dot(1,2)),
      eta(2 * m.column_dot(0,2)),
      zeta(2 * m.column_dot(0,1)) {}

  explicit GruberVector(const std::array<double,6>& g6)
    : A(g6[0]), B(g6[1]), C(g6[2]), xi(g6[3]), eta(g6[4]), zeta(g6[5]) {}

  std::array<double,6> parameters() const { return {A, B, C, xi, eta, zeta}; }
  std::array<double,6> cell_parameters() const {
    // inverse of UnitCell::g6()
    double a = std::sqrt(A);
    double b = std::sqrt(B);
    double c = std::sqrt(C);
    return {a, b, c,
            deg(std::acos(xi/(2*b*c))),
            deg(std::acos(eta/(2*a*c))),
            deg(std::acos(zeta/(2*a*b)))};
  }

  SellingVector selling() const;

  bool is_normalized() const {
    // eq(3) from Gruber 1973
    return A <= B && B <= C &&
           (A != B || std::abs(xi) <= std::abs(eta)) &&
           (B != C || std::abs(eta) <= std::abs(zeta)) &&
           (xi > 0) == (eta > 0) && (xi > 0) == (zeta > 0);
  }

  bool is_buerger(double epsilon=1e-9) const {
    return is_normalized() &&
           // eq (4) from Gruber 1973
           std::abs(xi) <= B + epsilon &&
           std::abs(eta) <= A + epsilon &&
           std::abs(zeta) <= A + epsilon;
  }

  // Algorithm N from Gruber (1973).
  // Returns branch taken in N3.
  void normalize(double eps=1e-9) {
    if (A - B > eps || (A - B >= -eps && std::abs(xi) > std::abs(eta) + eps)) { // N1
      std::swap(A, B);
      std::swap(xi, eta);
    }
    if (B - C > eps || (B - C >= -eps && std::abs(eta) > std::abs(zeta) + eps)) { // N2
      std::swap(B, C);
      std::swap(eta, zeta);
      // To make it faster, instead of "go to the point N1" we repeat N1 once
      // (which is equivalent - three swaps are sufficient to reorder ABC).
      if (A - B > eps || (A - B >= -eps && std::abs(xi) > std::abs(eta) + eps)) { // N1
        std::swap(A, B);
        std::swap(xi, eta);
      }
    }
    // N3
    // xi * eta * zeta > 0 <=> positive count is 1 or 3 and no zeros
    int pos_count = (xi > eps) + (eta > eps) + (zeta > eps);
    int nonneg_count = (xi >= -eps) + (eta >= -eps) + (zeta >= -eps);
    double sgn = (pos_count == nonneg_count && pos_count % 2 == 1) ? 1 : -1;
    xi = std::copysign(xi, sgn);
    eta = std::copysign(eta, sgn);
    zeta = std::copysign(zeta, sgn);
  }

  // Algorithm B from Gruber (1973).
  // Returns true if no change was needed.
  bool buerger_step() {
    if (std::abs(xi) > B) { // B2
      double j = std::floor(0.5*xi/B + 0.5);
      C += j * (j*B - xi);
      xi -= 2 * j * B;
      eta -= j * zeta;
    } else if (std::abs(eta) > A) { // B3
      double j = std::floor(0.5*eta/A + 0.5);
      C += j * (j*A - eta);
      xi -= j * zeta;
      eta -= 2 * j * A;
    } else if (std::abs(zeta) > A) { // B4
      double j = std::floor(0.5*zeta/A + 0.5);
      B += j * (j*A - zeta);
      xi -= j * eta;
      zeta -= 2 * j * A;
    } else if (xi + eta + zeta + A + B < 0) { // B5
      double j = std::floor(0.5 * (xi + eta) / (A + B + zeta) + 0.5);
      C += j * (j * (A + B + zeta) - (xi + eta));
      xi -= j * (2*B + zeta);
      eta -= j * (2*A + zeta);
    } else {
      return true;
    }
    return false;
  }

  // Returns number of iterations.
  int buerger_reduce() {
    int n = 0;
    double prev_sum = -1;
    int stall_count = 0;
    for (;;) {
      normalize();
      // In rare cases numerical errors push the algorithm into infinite loop,
      // as described in Grosse-Kunstleve et al, Acta Cryst. (2004) A60, 1.
      // Ad-hoc solution: stop if a+b+c is stalled for 5 iterations.
      if (++n > 8) {  // don't waste time during the first few iterations
        double sum = std::sqrt(A) + std::sqrt(B) + std::sqrt(C);
        if (std::abs(sum - prev_sum) < sum * 1e-6) {
          if (++stall_count == 5)
            break;
        } else {
          stall_count = 0;
        }
        prev_sum = sum;
      }
      if (buerger_step())
        break;
    }
    return n;
  }

  // To be called after normalize() or is_normalized().
  // Returns true if it already was Niggli cell.
  // Algorithm from Krivy & Gruber, Acta Cryst. (1976) A32, 297.
  bool niggli_step(double epsilon=1e-9) {
    if (std::abs(xi) > B + epsilon ||  // step 5. from Krivy-Gruber (1976)
        (xi >= B - epsilon && 2 * eta < zeta - epsilon) ||
        (xi <= -(B - epsilon) && zeta < -epsilon)) {
      double sign_xi = xi >= 0 ? 1 : -1;
      C += B - xi * sign_xi;
      eta -= zeta * sign_xi;
      xi -= 2 * B * sign_xi;
    } else if (std::abs(eta) > A + epsilon ||  // step 6.
               (eta >= A - epsilon && 2 * xi < zeta - epsilon) ||
               (eta <= -(A - epsilon) && zeta < -epsilon)) {
      double sign_eta = eta >= 0 ? 1 : -1;
      C += A - eta * sign_eta;
      xi -= zeta * sign_eta;
      eta -= 2 * A * sign_eta;
    } else if (std::abs(zeta) > A + epsilon ||  // step 7.
               (zeta >= A - epsilon && 2 * xi < eta - epsilon) ||
               (zeta <= -(A - epsilon) && eta < -epsilon)) {
      double sign_zeta = zeta >= 0 ? 1 : -1;
      B += A - zeta * sign_zeta;
      xi -= eta * sign_zeta;
      zeta -= 2 * A * sign_zeta;
    } else if (xi + eta + zeta + A + B < -epsilon || // step 8.
               (xi + eta + zeta + A + B <= epsilon && 2 * (A + eta) + zeta > epsilon)) {
      C += A + B + xi + eta + zeta;
      xi += 2 * B + zeta;
      eta += 2 * A + zeta;
    } else {
      return true;
    }
    return false;
  }

  // Returns number of iterations.
  int niggli_reduce(double epsilon=1e-9, int iteration_limit=100) {
    int n = 0;
    for (;;) {
      normalize(epsilon);
      if (++n == iteration_limit || niggli_step(epsilon))
        break;
    }
    return n;
  }

  bool is_niggli(double epsilon=1e-9) const {
    return is_normalized() && GruberVector(*this).niggli_step(epsilon);
  }
};


// Selling-Delaunay reduction. Based on:
// - chapter "Delaunay reduction and standardization" in
//   International Tables for Crystallography vol. A (2016), sec. 3.1.2.3.
//   https://onlinelibrary.wiley.com/iucr/itc/Ac/ch3o1v0001/
// - Patterson & Love (1957), Acta Cryst. 10, 111,
//   "Remarks on the Delaunay reduction", doi:10.1107/s0365110x57000328
// - Andrews et al (2019), Acta Cryst. A75, 115,
//   "Selling reduction versus Niggli reduction for crystallographic lattices".
struct SellingVector {
  // b.c a.c a.b a.d b.d c.d
  std::array<double,6> s;

  explicit SellingVector(const std::array<double,6>& s_) : s(s_) {}

  explicit SellingVector(const Mat33& orth) {
    Vec3 b[4];
    for (int i = 0; i < 3; ++i)
      b[i] = orth.column_copy(i);
    b[3]= -b[0] - b[1] - b[2];
    s[0] = b[1].dot(b[2]);
    s[1] = b[0].dot(b[2]);
    s[2] = b[0].dot(b[1]);
    s[3] = b[0].dot(b[3]);
    s[4] = b[1].dot(b[3]);
    s[5] = b[2].dot(b[3]);
  }

  // The reduction minimizes the sum b_i^2 which is equal to -2 sum s_i.
  double sum_b_squared() const {
    return -2 * (s[0] + s[1] + s[2] + s[3] + s[4] + s[5]);
  }

  bool is_reduced(double eps=1e-9) const {
    return std::all_of(s.begin(), s.end(), [eps](double x) { return x <= eps; });
  }

  bool reduce_step(double eps=1e-9) {
    //printf(" s = %g %g %g %g %g %g  sum=%g\n",
    //       s[0], s[1], s[2], s[3], s[4], s[5], sum_b_squared());
    const int table[6][5] = {
      // When negating s[n] we need to apply operations from table[n]:
      // 2 x add, subtract, 2 x swap&add
      {2, 4, 3, 1, 5},  // 0
      {2, 3, 4, 0, 5},  // 1
      {1, 3, 5, 0, 4},  // 2
      {1, 2, 0, 4, 5},  // 3
      {0, 2, 1, 3, 5},  // 4
      {0, 1, 2, 3, 4},  // 5
    };

    double max_s = eps;
    int max_s_pos = -1;
    for (int i = 0; i < 6; ++i)
      if (s[i] > max_s) {
        max_s = s[i];
        max_s_pos = i;
      }
    if (max_s_pos < 0)
      return false;
    const int (&indices)[5] = table[max_s_pos];
    s[max_s_pos] = -max_s;
    s[indices[0]] += max_s;
    s[indices[1]] += max_s;
    s[indices[2]] -= max_s;
    std::swap(s[indices[3]], s[indices[4]]);
    s[indices[3]] += max_s;
    s[indices[4]] += max_s;
    //printf("  s[%d]=%g  sum: %g\n", max_s_pos, max_s, sum_b_squared());
    return true;
  }

  // Returns number of iterations.
  int reduce(double eps=1e-9, int iteration_limit=100) {
    int n = 0;
    while (++n != iteration_limit)
      if (!reduce_step(eps))
        break;
    return n;
  }

  std::array<double,6> g6_parameters() const {
    return {-s[1]-s[2]-s[3], -s[0]-s[2]-s[4], -s[0]-s[1]-s[5], 2*s[0], 2*s[1], 2*s[2]};
  }

  GruberVector gruber() const { return GruberVector(g6_parameters()); }

  // Swap values to make a <= b <= c <= d
  void sort(double eps=1e-9) {
    double abcd_sq_neg[4] = {
      // -a^2, -b^2, -c^2, -d^2 (negated - to be sorted in descending order)
      s[1]+s[2]+s[3], s[0]+s[2]+s[4], s[0]+s[1]+s[5], s[3]+s[4]+s[5]
    };
    // First, make sure that d >= a,b,c (therefore -d^2 <= -a^2,...).
    int min_idx = 3;
    for (int i = 0; i < 3; ++i)
      if (abcd_sq_neg[i] < abcd_sq_neg[min_idx] - eps)
        min_idx = i;
    switch (min_idx) {
      case 0:  // a <-> d
        std::swap(s[1], s[5]);
        std::swap(s[2], s[4]);
        break;
      case 1:  // b <-> d
        std::swap(s[0], s[5]);
        std::swap(s[2], s[3]);
        break;
      case 2:  // c <-> d
        std::swap(s[0], s[4]);
        std::swap(s[1], s[3]);
        break;
    }
    // we could stop here and not care about the order of a,b,c.
    std::swap(abcd_sq_neg[min_idx], abcd_sq_neg[3]);
    if (abcd_sq_neg[0] < abcd_sq_neg[1] - eps) {  // a <-> b
      std::swap(s[0], s[1]);
      std::swap(s[3], s[4]);
      std::swap(abcd_sq_neg[0], abcd_sq_neg[1]);
    }
    if (abcd_sq_neg[1] < abcd_sq_neg[2] - eps) {  // b <-> c
      std::swap(s[1], s[2]);
      std::swap(s[4], s[5]);
      std::swap(abcd_sq_neg[1], abcd_sq_neg[2]);
    }
    if (abcd_sq_neg[0] < abcd_sq_neg[1] - eps) {  // a <-> b
      std::swap(s[0], s[1]);
      std::swap(s[3], s[4]);
      //std::swap(abcd_sq_neg[0], abcd_sq_neg[1]);
    }
  }

  std::array<double,6> cell_parameters() const {
    return gruber().cell_parameters();
  }
};

SellingVector GruberVector::selling() const {
  double s0 = 0.5 * xi;
  double s1 = 0.5 * eta;
  double s2 = 0.5 * zeta;
  return SellingVector({s0, s1, s2, -A - s1 - s2, -B - s0 - s2, -C - s0 - s1});
}

} // namespace gemmi
#endif
