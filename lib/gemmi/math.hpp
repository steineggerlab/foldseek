// Copyright 2018 Global Phasing Ltd.
//
// Math utilities. 3D linear algebra.

#ifndef GEMMI_MATH_HPP_
#define GEMMI_MATH_HPP_

#include <cmath>      // for fabs, cos, sqrt, round
#include <cstdio>     // for snprintf
#include <algorithm>  // for min
#include <array>
#include <stdexcept>  // for out_of_range
#include <string>
#include <vector>

namespace gemmi {

constexpr double pi() { return 3.1415926535897932384626433832795029; }

// The value used in converting between energy[eV] and wavelength[Angstrom].
// $ units -d15 'h * c / eV / angstrom'
constexpr double hc() { return 12398.4197386209; }

// The Bohr radius (a0) in Angstroms.
constexpr double bohrradius() { return 0.529177210903; }

// Used in conversion of ADPs (atomic displacement parameters).
constexpr double u_to_b() { return 8 * pi() * pi(); }

constexpr double deg(double angle) { return 180.0 / pi() * angle; }
constexpr double rad(double angle) { return pi() / 180.0 * angle; }

constexpr float sq(float x) { return x * x; }
constexpr double sq(double x) { return x * x; }

inline int iround(double d) { return static_cast<int>(std::round(d)); }

inline double angle_abs_diff(double a, double b, double full=360.0) {
  double d = std::abs(a - b);
  if (d > full)
    d -= std::floor(d / full) * full;
  return std::min(d, full - d);
}

struct Vec3 {
  double x, y, z;

  Vec3() : x(0), y(0), z(0) {}
  Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
  explicit Vec3(std::array<int, 3> h) : x(h[0]), y(h[1]), z(h[2]) {}

  double& at(int i) {
    switch (i) {
      case 0: return x;
      case 1: return y;
      case 2: return z;
      default: throw std::out_of_range("Vec3 index must be 0, 1 or 2.");
    }
  }
  double at(int i) const { return const_cast<Vec3*>(this)->at(i); }

  Vec3 operator-() const { return {-x, -y, -z}; }
  Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
  Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
  Vec3 operator*(double d) const { return {x*d, y*d, z*d}; }
  Vec3 operator/(double d) const { return *this * (1.0/d); }
  Vec3& operator-=(const Vec3& o) { *this = *this - o; return *this; }
  Vec3& operator+=(const Vec3& o) { *this = *this + o; return *this; }
  Vec3& operator*=(double d) { *this = *this * d; return *this; }
  Vec3& operator/=(double d) { return operator*=(1.0/d); }

  Vec3 negated() const { return {-x, -y, -z}; }
  double dot(const Vec3& o) const { return x*o.x + y*o.y + z*o.z; }
  Vec3 cross(const Vec3& o) const {
    return {y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x};
  }
  double length_sq() const { return x * x + y * y + z * z; }
  double length() const { return std::sqrt(length_sq()); }
  Vec3 changed_magnitude(double m) const { return operator*(m / length()); }
  Vec3 normalized() const { return changed_magnitude(1.0); }
  double dist_sq(const Vec3& o) const { return (*this - o).length_sq(); }
  double dist(const Vec3& o) const { return std::sqrt(dist_sq(o)); }
  double angle(const Vec3& o) const {
    return std::acos(dot(o) / std::sqrt(length_sq() * o.length_sq()));
  }
  bool approx(const Vec3& o, double epsilon) const {
    return std::fabs(x - o.x) <= epsilon &&
           std::fabs(y - o.y) <= epsilon &&
           std::fabs(z - o.z) <= epsilon;
  }
  std::string str() const {
    using namespace std;
    char buf[64] = {0};
    snprintf(buf, 63, "[%g %g %g]", x, y, z);
    return buf;
  }
};

inline Vec3 operator*(double d, const Vec3& v) { return v * d; }

struct Mat33 {
  double a[3][3] = { {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} };

  // make it accessible with ".a"
  typedef double row_t[3];
  const row_t& operator[](int i) const { return a[i]; }
  row_t& operator[](int i) { return a[i]; }

  Mat33() = default;
  explicit Mat33(double d) : a{{d, d, d}, {d, d, d}, {d, d, d}} {}
  Mat33(double a1, double a2, double a3, double b1, double b2, double b3,
        double c1, double c2, double c3)
  : a{{a1, a2, a3}, {b1, b2, b3}, {c1, c2, c3}} {}

  Vec3 row_copy(int i) const {
    if (i < 0 || i > 2)
      throw std::out_of_range("Mat33 row index must be 0, 1 or 2.");
    return Vec3(a[i][0], a[i][1], a[i][2]);
  }

  Vec3 column_copy(int i) const {
    if (i < 0 || i > 2)
      throw std::out_of_range("Mat33 column index must be 0, 1 or 2.");
    return Vec3(a[0][i], a[1][i], a[2][i]);
  }

  Vec3 multiply(const Vec3& p) const {
    return {a[0][0] * p.x + a[0][1] * p.y + a[0][2] * p.z,
            a[1][0] * p.x + a[1][1] * p.y + a[1][2] * p.z,
            a[2][0] * p.x + a[2][1] * p.y + a[2][2] * p.z};
  }
  Vec3 left_multiply(const Vec3& p) const {
    return {a[0][0] * p.x + a[1][0] * p.y + a[2][0] * p.z,
            a[0][1] * p.x + a[1][1] * p.y + a[2][1] * p.z,
            a[0][2] * p.x + a[1][2] * p.y + a[2][2] * p.z};
  }
  // p has elements from the main diagonal of a 3x3 diagonal matrix
  Mat33 multiply_by_diagonal(const Vec3& p) const {
    return Mat33(a[0][0] * p.x, a[0][1] * p.y, a[0][2] * p.z,
                 a[1][0] * p.x, a[1][1] * p.y, a[1][2] * p.z,
                 a[2][0] * p.x, a[2][1] * p.y, a[2][2] * p.z);
  }
  Mat33 multiply(const Mat33& b) const {
    Mat33 r;
    for (int i = 0; i != 3; ++i)
      for (int j = 0; j != 3; ++j)
        r[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
    return r;
  }
  Mat33 transpose() const {
    return Mat33(a[0][0], a[1][0], a[2][0],
                 a[0][1], a[1][1], a[2][1],
                 a[0][2], a[1][2], a[2][2]);
  }
  double trace() const { return a[0][0] + a[1][1] + a[2][2]; }

  bool approx(const Mat33& other, double epsilon) const {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        if (std::fabs(a[i][j] - other.a[i][j]) > epsilon)
          return false;
    return true;
  }
  double determinant() const {
    return a[0][0] * (a[1][1]*a[2][2] - a[2][1]*a[1][2]) +
           a[0][1] * (a[1][2]*a[2][0] - a[2][2]*a[1][0]) +
           a[0][2] * (a[1][0]*a[2][1] - a[2][0]*a[1][1]);
  }
  Mat33 inverse() const {
    Mat33 inv;
    double inv_det = 1.0 / determinant();
    inv[0][0] = inv_det * (a[1][1] * a[2][2] - a[2][1] * a[1][2]);
    inv[0][1] = inv_det * (a[0][2] * a[2][1] - a[0][1] * a[2][2]);
    inv[0][2] = inv_det * (a[0][1] * a[1][2] - a[0][2] * a[1][1]);
    inv[1][0] = inv_det * (a[1][2] * a[2][0] - a[1][0] * a[2][2]);
    inv[1][1] = inv_det * (a[0][0] * a[2][2] - a[0][2] * a[2][0]);
    inv[1][2] = inv_det * (a[1][0] * a[0][2] - a[0][0] * a[1][2]);
    inv[2][0] = inv_det * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
    inv[2][1] = inv_det * (a[2][0] * a[0][1] - a[0][0] * a[2][1]);
    inv[2][2] = inv_det * (a[0][0] * a[1][1] - a[1][0] * a[0][1]);
    return inv;
  }
  bool is_identity() const {
    return a[0][0] == 1 && a[0][1] == 0 && a[0][2] == 0 &&
           a[1][0] == 0 && a[1][1] == 1 && a[1][2] == 0 &&
           a[2][0] == 0 && a[2][1] == 0 && a[2][2] == 1;
  }

  double column_dot(int i, int j) const {
    return a[0][i] * a[0][j] + a[1][i] * a[1][j] + a[2][i] * a[2][j];
  }
};

// Symmetric matrix 3x3. Used primarily for an ADP tensor.
template<typename T> struct SMat33 {
  T u11, u22, u33, u12, u13, u23;

  // The PDB ANISOU record has the above order, but in a different context
  // (such as metric tensor) the order of Voigt notation may be preferred.
  std::array<T, 6> elements_pdb() const   { return {{u11, u22, u33, u12, u13, u23}}; }
  std::array<T, 6> elements_voigt() const { return {{u11, u22, u33, u23, u13, u12}}; }

  Mat33 as_mat33() const {
    return Mat33(u11, u12, u13, u12, u22, u23, u13, u23, u33);
  }

  T trace() const { return u11 + u22 + u33; }
  bool nonzero() const { return trace() != 0; }

  bool all_zero() const {
    return u11 == 0 && u22 == 0 && u33 == 0 && u12 == 0 && u13 == 0 && u23 == 0;
  }

  void scale(T s) const {
    u11 *= s; u22 *= s; u33 *= s; u12 *= s; u13 *= s; u23 *= s;
  };

  template<typename Real>
  SMat33<Real> scaled(Real s) const {
    return SMat33<Real>{u11*s, u22*s, u33*s, u12*s, u13*s, u23*s};
  }

  // returns U + kI
  SMat33<T> added_kI(T k) const {
    return {u11+k, u22+k, u33+k, u12, u13, u23};
  }

  // returns squared norm r^T U r where U is this matrix and vector r is arg
  double r_u_r(const Vec3& r) const {
    return r.x * r.x * u11 + r.y * r.y * u22 + r.z * r.z * u33 +
      2 * (r.x * r.y * u12 + r.x * r.z * u13 + r.y * r.z * u23);
  }
  double r_u_r(const std::array<int,3>& h) const {
    // it's faster to first convert ints to doubles (Vec3)
    return r_u_r(Vec3(h));
  }

  Vec3 multiply(const Vec3& p) const {
    return {u11 * p.x + u12 * p.y + u13 * p.z,
            u12 * p.x + u22 * p.y + u23 * p.z,
            u13 * p.x + u23 * p.y + u33 * p.z};
  }

  SMat33 operator-(const SMat33& o) const {
    return {u11-o.u11, u22-o.u22, u33-o.u33, u12-o.u12, u13-o.u13, u23-o.u23};
  }
  SMat33 operator+(const SMat33& o) const {
    return {u11+o.u11, u22+o.u22, u33+o.u33, u12+o.u12, u13+o.u13, u23+o.u23};
  }

  // return M U M^T
  template<typename Real=double>
  SMat33<Real> transformed_by(const Mat33& m) const {
    // slightly faster than m.multiply(as_mat33()).multiply(m.transpose());
    auto elem = [&](int i, int j) {
      return static_cast<Real>(
          m[i][0] * (m[j][0] * u11 + m[j][1] * u12 + m[j][2] * u13) +
          m[i][1] * (m[j][0] * u12 + m[j][1] * u22 + m[j][2] * u23) +
          m[i][2] * (m[j][0] * u13 + m[j][1] * u23 + m[j][2] * u33));
    };
    return SMat33<Real>{elem(0, 0), elem(1, 1), elem(2, 2),
                        elem(0, 1), elem(0, 2), elem(1, 2)};
  }

  T determinant() const {
    return u11 * (u22*u33 - u23*u23) +
           u12 * (u23*u13 - u33*u12) +
           u13 * (u12*u23 - u13*u22);
  }

  SMat33 inverse() const {
    SMat33 inv;
    T inv_det = 1.0f / determinant();
    inv.u11 = inv_det * (u22 * u33 - u23 * u23);
    inv.u22 = inv_det * (u11 * u33 - u13 * u13);
    inv.u33 = inv_det * (u11 * u22 - u12 * u12);
    inv.u12 = inv_det * (u13 * u23 - u12 * u33);
    inv.u13 = inv_det * (u12 * u23 - u13 * u22);
    inv.u23 = inv_det * (u12 * u13 - u11 * u23);
    return inv;
  }

  // Based on https://en.wikipedia.org/wiki/Eigenvalue_algorithm
  std::array<double, 3> calculate_eigenvalues() const {
    double p1 = u12*u12 + u13*u13 + u23*u23;
    if (p1 == 0)
      return {{u11, u22, u33}};
    double q = (1./3.) * trace();
    SMat33<double> b{u11 - q, u22 - q, u33 - q, u12, u13, u23};
    double p2 = sq(b.u11) + sq(b.u22) + sq(b.u33) + 2 * p1;
    double p = std::sqrt((1./6.) * p2);
    double r = b.determinant() / ((1./3.) * p2 * p);
    double phi = 0;
    if (r <= -1)
      phi = (1./3.) * pi();
    else if (r < 1)
      phi = (1./3.) * std::acos(r);
    double eig1 = q + 2 * p * std::cos(phi);
    double eig3 = q + 2 * p * std::cos(phi + 2./3.*pi());
    return {{eig1, 3 * q - eig1 - eig3, eig3}};
  }

  // Assumes one of the eigenvalue calculate above.
  // May not work if eigenvalues are not distinct.
  Vec3 calculate_eigenvector(double eigenvalue) const {
    Vec3 r0(u11 - eigenvalue, u12, u13);
    Vec3 r1(u12, u22 - eigenvalue, u23);
    Vec3 r2(u13, u23, u33 - eigenvalue);
    Vec3 cr[3] = {r0.cross(r1), r0.cross(r2), r1.cross(r2)};
    int idx = 0;
    double lensq = 0;
    for (int i = 0; i < 3; ++i) {
      double tmp = cr[i].length_sq();
      if (tmp > lensq) {
        idx = i;
        lensq = tmp;
      }
    }
    if (lensq == 0)
      return Vec3(0, 0, 1); // an arbitrary choice for the special case
    return cr[idx] / std::sqrt(lensq);
  }
};

struct Transform {
  Mat33 mat;
  Vec3 vec;

  Transform inverse() const {
    Mat33 minv = mat.inverse();
    return {minv, minv.multiply(vec).negated()};
  }

  Vec3 apply(const Vec3& x) const { return mat.multiply(x) + vec; }

  Transform combine(const Transform& b) const {
    return {mat.multiply(b.mat), vec + mat.multiply(b.vec)};
  }

  bool is_identity() const {
    return mat.is_identity() && vec.x == 0. && vec.y == 0. && vec.z == 0.;
  }
  void set_identity() { mat = Mat33(); vec = Vec3(); }

  bool approx(const Transform& o, double epsilon) const {
    return mat.approx(o.mat, epsilon) && vec.approx(o.vec, epsilon);
  }
};

template<typename Pos>
struct Box {
  Pos minimum = Pos(INFINITY, INFINITY, INFINITY);
  Pos maximum = Pos(-INFINITY, -INFINITY, -INFINITY);
  void extend(const Pos& p) {
    if (p.x < minimum.x) minimum.x = p.x;
    if (p.y < minimum.y) minimum.y = p.y;
    if (p.z < minimum.z) minimum.z = p.z;
    if (p.x > maximum.x) maximum.x = p.x;
    if (p.y > maximum.y) maximum.y = p.y;
    if (p.z > maximum.z) maximum.z = p.z;
  }
  Pos get_size() const { return maximum - minimum; }
  void add_margins(const Pos& p) { minimum -= p; maximum += p; }
  void add_margin(double m) { add_margins(Pos(m, m, m)); }
};

// popular single-pass algorithm for calculating variance and mean
struct Variance {
  int n = 0;
  double sum_sq = 0.;
  double mean_x = 0.;

  Variance() = default;
  template <typename T> Variance(T begin, T end) : Variance() {
    for (auto i = begin; i != end; ++i)
      add_point(*i);
  }
  void add_point(double x) {
    ++n;
    double dx = x - mean_x;
    mean_x += dx / n;
    sum_sq += dx * (x - mean_x);
  }
  double for_sample() const { return sum_sq / (n - 1); }
  double for_population() const { return sum_sq / n; }
};

struct Covariance : Variance {
  double mean_y = 0.;
  void add_point(double x, double y) {
    ++n;
    double dx = x - mean_x;
    mean_x += dx / n;
    mean_y += (y - mean_y) / n;
    sum_sq += dx * (y - mean_y);
  }
};

struct Correlation {
  int n = 0;
  double sum_xx = 0.;
  double sum_yy = 0.;
  double sum_xy = 0.;
  double mean_x = 0.;
  double mean_y = 0.;
  void add_point(double x, double y) {
    ++n;
    double weight = (n - 1.0) / n;
    double dx = x - mean_x;
    double dy = y - mean_y;
    sum_xx += weight * dx * dx;
    sum_yy += weight * dy * dy;
    sum_xy += weight * dx * dy;
    mean_x += dx / n;
    mean_y += dy / n;
  }
  double coefficient() const { return sum_xy / std::sqrt(sum_xx * sum_yy); }
  double x_variance() const { return sum_xx / n; }
  double y_variance() const { return sum_yy / n; }
  double covariance() const { return sum_xy / n; }
  double mean_ratio() const { return mean_y / mean_x; }
  // the regression line
  double slope() const { return sum_xy / sum_xx; }
  double intercept() const { return mean_y - slope() * mean_x; }
};


struct DataStats {
  double dmin = NAN;
  double dmax = NAN;
  double dmean = NAN;
  double rms = NAN;
  size_t nan_count = 0;
};

template<typename T>
DataStats calculate_data_statistics(const std::vector<T>& data) {
  DataStats stats;
  double sum = 0;
  double sq_sum = 0;
  stats.dmin = INFINITY;
  stats.dmax = -INFINITY;
  for (double d : data) {
    if (std::isnan(d)) {
      stats.nan_count++;
      continue;
    }
    sum += d;
    sq_sum += d * d;
    if (d < stats.dmin)
      stats.dmin = d;
    if (d > stats.dmax)
      stats.dmax = d;
  }
  if (stats.nan_count != data.size()) {
    stats.dmean = sum / (data.size() - stats.nan_count);
    stats.rms = std::sqrt(sq_sum / (data.size() - stats.nan_count) - stats.dmean * stats.dmean);
  } else {
    stats.dmin = NAN;
    stats.dmax = NAN;
  }
  return stats;
}

// internally used functions
namespace impl {
template<typename T> bool is_same(T a, T b) { return a == b; }
template<> inline bool is_same(float a, float b) {
  return std::isnan(b) ? std::isnan(a) : a == b;
}
template<> inline bool is_same(double a, double b) {
  return std::isnan(b) ? std::isnan(a) : a == b;
}
} // namespace impl

} // namespace gemmi
#endif
