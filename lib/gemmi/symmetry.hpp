// Copyright 2017-2019 Global Phasing Ltd.
//
// Crystallographic Symmetry. Space Groups. Coordinate Triplets.
//
// If this is all that you need from Gemmi you can just copy this file,
// fail.hpp and LICENSE.txt to your project.

#ifndef GEMMI_SYMMETRY_HPP_
#define GEMMI_SYMMETRY_HPP_

#include <cstdint>
#include <cstdlib>    // for strtol
#include <cstring>    // for memchr, strchr
#include <array>
#include <algorithm>  // for count, sort, remove
#include <functional> // for hash
#include <stdexcept>  // for runtime_error, invalid_argument
#include <string>
#include <tuple>      // for tie
#include <vector>

#include "fail.hpp"   // for fail, unreachable

// we use brace elision with std:array's
#ifdef __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
#endif

namespace gemmi {

// UTILS

namespace impl {

// copied a helper function from atox.hpp to keep it a two-header lib
inline const char* skip_blank(const char* p) {
  if (p)
    while (*p == ' ' || *p == '\t' || *p == '_') // '_' can be used as space
      ++p;
  return p;
}

} // namespace impl


// OP

// Op is a symmetry operation, or a change-of-basic transformation,
// or a different operation of similar kind.
// Both "rotation" matrix and translation vector are fractional, with DEN
// used as the denominator.
struct Op {
  static constexpr int DEN = 24;  // 24 to handle 1/8 in change-of-basis
  typedef std::array<std::array<int, 3>, 3> Rot;
  typedef std::array<int, 3> Tran;

  Rot rot;
  Tran tran;

  std::string triplet() const;

  Op inverse() const;

  // If the translation points outside of the unit cell, wrap it.
  Op& wrap() {
    for (int i = 0; i != 3; ++i) {
      if (tran[i] >= DEN) // elements need to be in [0,DEN)
        tran[i] %= DEN;
      else if (tran[i] < 0)
        tran[i] = ((tran[i] + 1) % DEN) + DEN - 1;
    }
    return *this;
  }

  Op& translate(const Tran& a) {
    for (int i = 0; i != 3; ++i)
      tran[i] += a[i];
    return *this;
  }

  Op translated(const Tran& a) const { return Op(*this).translate(a); }

  Op add_centering(const Tran& a) const { return translated(a).wrap(); }

  Rot negated_rot() const {
    return { -rot[0][0], -rot[0][1], -rot[0][2],
             -rot[1][0], -rot[1][1], -rot[1][2],
             -rot[2][0], -rot[2][1], -rot[2][2] };
  }

  Op negated() const { return { negated_rot(), { -tran[0], -tran[1], -tran[2] } }; }

  Rot transposed_rot() const {
    return { rot[0][0], rot[1][0], rot[2][0],
             rot[0][1], rot[1][1], rot[2][1],
             rot[0][2], rot[1][2], rot[2][2] };
  }

  // DEN^3 for rotation, -DEN^3 for rotoinversion
  int det_rot() const {
    return rot[0][0] * (rot[1][1] * rot[2][2] - rot[1][2] * rot[2][1])
         - rot[0][1] * (rot[1][0] * rot[2][2] - rot[1][2] * rot[2][0])
         + rot[0][2] * (rot[1][0] * rot[2][1] - rot[1][1] * rot[2][0]);
  }

  Op combine(const Op& b) const {
    Op r;
    for (int i = 0; i != 3; ++i) {
      r.tran[i] = tran[i] * Op::DEN;
      for (int j = 0; j != 3; ++j) {
        r.rot[i][j] = (rot[i][0] * b.rot[0][j] +
                       rot[i][1] * b.rot[1][j] +
                       rot[i][2] * b.rot[2][j]) / Op::DEN;
        r.tran[i] += rot[i][j] * b.tran[j];
      }
      r.tran[i] /= Op::DEN;
    }
    return r;
  }

  std::array<double, 3> apply_to_xyz(const std::array<double, 3>& xyz) const {
    std::array<double, 3> out;
    for (int i = 0; i != 3; ++i)
      out[i] = (rot[i][0] * xyz[0] + rot[i][1] * xyz[1] + rot[i][2] * xyz[2] +
                tran[i]) / Op::DEN;
    return out;
  }

  // Miller is defined in the same way in namespace gemmi in unitcell.hpp
  using Miller = std::array<int, 3>;

  Miller apply_to_hkl_without_division(const Miller& hkl) const {
    Miller r;
    for (int i = 0; i != 3; ++i)
      r[i] = (rot[0][i] * hkl[0] + rot[1][i] * hkl[1] + rot[2][i] * hkl[2]);
    return r;
  }
  static Miller divide_hkl_by_DEN(const Miller& hkl) {
    return {{ hkl[0] / DEN, hkl[1] / DEN, hkl[2] / DEN }};
  }
  Miller apply_to_hkl(const Miller& hkl) const {
    return divide_hkl_by_DEN(apply_to_hkl_without_division(hkl));
  }

  double phase_shift(const Miller& hkl) const {
    constexpr double mult = -2 * 3.1415926535897932384626433832795 / Op::DEN;
    return mult * (hkl[0] * tran[0] + hkl[1] * tran[1] + hkl[2] * tran[2]);
  }

  std::array<std::array<int, 4>, 4> int_seitz() const {
    std::array<std::array<int, 4>, 4> t;
    for (int i = 0; i < 3; ++i)
      t[i] = { rot[i][0], rot[i][1], rot[i][2], tran[i] };
    t[3] = { 0, 0, 0, 1 };
    return t;
  }

  std::array<std::array<double, 4>, 4> float_seitz() const {
    std::array<std::array<double, 4>, 4> t;
    double m = 1.0 / Op::DEN;
    for (int i = 0; i < 3; ++i)
      t[i] = { m * rot[i][0], m * rot[i][1], m * rot[i][2], m * tran[i] };
    t[3] = { 0., 0., 0., 1. };
    return t;
  }

  static constexpr Op identity() {
    return {{DEN,0,0, 0,DEN,0, 0,0,DEN}, {0,0,0}};
  }
  bool operator<(const Op& rhs) const {
    return std::tie(rot, tran) < std::tie(rhs.rot, rhs.tran);
  }
};

inline bool operator==(const Op& a, const Op& b) {
  return a.rot == b.rot && a.tran == b.tran;
}
inline bool operator!=(const Op& a, const Op& b) { return !(a == b); }

inline Op operator*(const Op& a, const Op& b) { return a.combine(b).wrap(); }
inline Op& operator*=(Op& a, const Op& b) { a = a * b; return a; }

inline Op Op::inverse() const {
  int detr = det_rot();
  if (detr == 0)
    fail("cannot invert matrix: " + Op{rot, {0,0,0}}.triplet());
  int d2 = Op::DEN * Op::DEN;
  Op inv;
  inv.rot[0][0] = d2 * (rot[1][1] * rot[2][2] - rot[2][1] * rot[1][2]) / detr;
  inv.rot[0][1] = d2 * (rot[0][2] * rot[2][1] - rot[0][1] * rot[2][2]) / detr;
  inv.rot[0][2] = d2 * (rot[0][1] * rot[1][2] - rot[0][2] * rot[1][1]) / detr;
  inv.rot[1][0] = d2 * (rot[1][2] * rot[2][0] - rot[1][0] * rot[2][2]) / detr;
  inv.rot[1][1] = d2 * (rot[0][0] * rot[2][2] - rot[0][2] * rot[2][0]) / detr;
  inv.rot[1][2] = d2 * (rot[1][0] * rot[0][2] - rot[0][0] * rot[1][2]) / detr;
  inv.rot[2][0] = d2 * (rot[1][0] * rot[2][1] - rot[2][0] * rot[1][1]) / detr;
  inv.rot[2][1] = d2 * (rot[2][0] * rot[0][1] - rot[0][0] * rot[2][1]) / detr;
  inv.rot[2][2] = d2 * (rot[0][0] * rot[1][1] - rot[1][0] * rot[0][1]) / detr;
  for (int i = 0; i != 3; ++i)
    inv.tran[i] = (-tran[0] * inv.rot[i][0]
                   -tran[1] * inv.rot[i][1]
                   -tran[2] * inv.rot[i][2]) / Op::DEN;
  return inv;
}


// TRIPLET -> OP

inline int interpret_miller_character(char c, const std::string& s) {
  static const signed char values[] =
    //a  b  c  d  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z
    { 1, 2, 3, 0, 0, 0, 0, 1, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3 };
  size_t idx = size_t((c | 0x20) - 'a');  // "|0x20" = to lower
  if (idx >= sizeof(values) || values[idx] == 0)
    fail("unexpected character '", c, "' in: ", s);
  return values[idx] - 1;
}

inline std::array<int, 4> parse_triplet_part(const std::string& s) {
  std::array<int, 4> r = { 0, 0, 0, 0 };
  int num = Op::DEN;
  const char* c = s.c_str();
  while (*(c = impl::skip_blank(c))) {
    if (*c == '+' || *c == '-') {
      num = (*c == '+' ? Op::DEN : -Op::DEN);
      c = impl::skip_blank(++c);
    }
    if (num == 0)
      fail("wrong or unsupported triplet format: " + s);
    int r_idx;
    int den = 1;
    if (*c >= '0' && *c <= '9') {
      // syntax examples in this branch: "1", "-1/2", "+2*x", "1/2 * b"
      char* endptr;
      num *= std::strtol(c, &endptr, 10);
      if (*endptr == '/')
        den = std::strtol(endptr + 1, &endptr, 10);
      if (*endptr == '*') {
        c = impl::skip_blank(endptr + 1);
        r_idx = interpret_miller_character(*c, s);
        ++c;
      } else {
        c = endptr;
        r_idx = 3;
      }
    } else {
      // syntax examples in this branch: "x", "+a", "-k/3"
      r_idx = interpret_miller_character(*c, s);
      c = impl::skip_blank(++c);
      if (*c == '/') {
        char* endptr;
        den = std::strtol(c + 1, &endptr, 10);
        c = endptr;
      }
    }
    if (den != 1) {
      if (den <= 0 || Op::DEN % den != 0)
        fail("Wrong denominator " + std::to_string(den) + " in: " + s);
      num /= den;
    }
    r[r_idx] += num;
    num = 0;
  }
  if (num != 0)
    fail("trailing sign in: " + s);
  return r;
}

inline Op parse_triplet(const std::string& s) {
  if (std::count(s.begin(), s.end(), ',') != 2)
    fail("expected exactly two commas in triplet");
  size_t comma1 = s.find(',');
  size_t comma2 = s.find(',', comma1 + 1);
  auto a = parse_triplet_part(s.substr(0, comma1));
  auto b = parse_triplet_part(s.substr(comma1 + 1, comma2 - (comma1 + 1)));
  auto c = parse_triplet_part(s.substr(comma2 + 1));
  Op::Rot rot = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]};
  Op::Tran tran = {a[3], b[3], c[3]};
  return { rot, tran };
}


// OP -> TRIPLET

namespace impl {

// much faster than s += std::to_string(n) for n in 0 ... 99
inline void append_small_number(std::string& s, int n) {
  if (n < 0 || n >= 100) {
    s += std::to_string(n);
  } else if (n < 10) {
    s += char('0' + n);
  } else { // 10 ... 99
    int tens = n / 10;
    s += char('0' + tens);
    s += char('0' + n - 10 * tens);
  }
}

inline void append_sign_of(std::string& s, int n) {
  if (n < 0)
    s += '-';
  else if (!s.empty())
    s += '+';
}

// append w/DEN fraction reduced to the lowest terms
inline std::pair<int,int> get_op_fraction(int w) {
  // Op::DEN == 24 == 2 * 2 * 2 * 3
  int denom = 1;
  for (int i = 0; i != 3; ++i)
    if (w % 2 == 0)  // 2, 2, 2
      w /= 2;
    else
      denom *= 2;
  if (w % 3 == 0)    // 3
    w /= 3;
  else
    denom *= 3;
  return {w, denom};
}

} // namespace impl

inline std::string make_triplet_part(const std::array<int, 3>& xyz, int w,
                                     char style='x') {
  std::string s;
  for (int i = 0; i != 3; ++i)
    if (xyz[i] != 0) {
      impl::append_sign_of(s, xyz[i]);
      int a = std::abs(xyz[i]);
      if (a != Op::DEN) {
        std::pair<int,int> frac = impl::get_op_fraction(a);
        if (frac.first == 1) {  // e.g. "x/3"
          s += char(style + i);
          s += '/';
          impl::append_small_number(s, frac.second);
        } else {  // e.g. "2/3*x"
          impl::append_small_number(s, frac.first);
          if (frac.second != 1) {
            s += '/';
            impl::append_small_number(s, frac.second);
          }
          s += '*';
          s += char(style + i);
        }
      } else {
        s += char(style + i);
      }
    }
  if (w != 0) {
    impl::append_sign_of(s, w);
    std::pair<int,int> frac = impl::get_op_fraction(std::abs(w));
    impl::append_small_number(s, frac.first);
    if (frac.second != 1) {
      s += '/';
      impl::append_small_number(s, frac.second);
    }
  }
  return s;
}

inline std::string Op::triplet() const {
  return make_triplet_part(rot[0], tran[0]) +
   "," + make_triplet_part(rot[1], tran[1]) +
   "," + make_triplet_part(rot[2], tran[2]);
}


// GROUPS OF OPERATIONS

// corresponds to Table A1.4.2.2 in ITfC vol.B (edition 2010)
inline std::vector<Op::Tran> centring_vectors(char centring_type) {
  constexpr int h = Op::DEN / 2;
  constexpr int t = Op::DEN / 3;
  constexpr int d = 2 * t;
  // note: find_centering() depends on the order of operations in vector
  switch (centring_type & ~0x20) {
    case 'P': return {{0, 0, 0}};
    case 'A': return {{0, 0, 0}, {0, h, h}};
    case 'B': return {{0, 0, 0}, {h, 0, h}};
    case 'C': return {{0, 0, 0}, {h, h, 0}};
    case 'I': return {{0, 0, 0}, {h, h, h}};
    case 'R': return {{0, 0, 0}, {d, t, t}, {t, d, d}};
    // hall_symbols.html has no H, ITfC 2010 has no S and T
    case 'H': return {{0, 0, 0}, {d, t, 0}, {t, d, 0}};
    case 'S': return {{0, 0, 0}, {t, t, d}, {d, t, d}};
    case 'T': return {{0, 0, 0}, {t, d, t}, {d, t, d}};
    case 'F': return {{0, 0, 0}, {0, h, h}, {h, 0, h}, {h, h, 0}};
    default: fail("not a centring type: ", centring_type);
  }
}


struct GroupOps {
  std::vector<Op> sym_ops;
  std::vector<Op::Tran> cen_ops;

  int order() const { return static_cast<int>(sym_ops.size()*cen_ops.size()); }

  void add_missing_elements();

  char find_centering() const {
    if (cen_ops.size() == 1 && cen_ops[0] == Op::Tran{0, 0, 0})
      return 'P';
    std::vector<Op::Tran> trans = cen_ops;
    std::sort(trans.begin(), trans.end());
    for (char c : {'A', 'B', 'C', 'I', 'F', 'R', 'H', 'S', 'T'}) {
      std::vector<Op::Tran> c_vectors = centring_vectors(c);
      if (c == 'R' || c == 'H') // these two are returned not sorted
        std::swap(c_vectors[1], c_vectors[2]);
      if (trans == c_vectors)
        return c;
    }
    return 0;
  }

  Op* find_by_rotation(const Op::Rot& r) {
    for (Op& op : sym_ops)
      if (op.rot == r)
        return &op;
    return nullptr;
  }

  const Op* find_by_rotation(const Op::Rot& r) const {
    return const_cast<GroupOps*>(this)->find_by_rotation(r);
  }

  bool is_centric() const {
    return find_by_rotation({-Op::DEN,0,0, 0,-Op::DEN,0, 0,0,-Op::DEN}) != nullptr;
  }

  bool is_reflection_centric(const Op::Miller& hkl) const {
    Op::Miller mhkl = {{-Op::DEN * hkl[0], -Op::DEN * hkl[1], -Op::DEN * hkl[2]}};
    for (const Op& op : sym_ops)
      if (op.apply_to_hkl_without_division(hkl) == mhkl)
        return true;
    return false;
  }

  int epsilon_factor_without_centering(const Op::Miller& hkl) const {
    Op::Miller denh = {{Op::DEN * hkl[0], Op::DEN * hkl[1], Op::DEN * hkl[2]}};
    int epsilon = 0;
    for (const Op& op : sym_ops)
      if (op.apply_to_hkl_without_division(hkl) == denh)
        ++epsilon;
    return epsilon;
  }
  int epsilon_factor(const Op::Miller& hkl) const {
    return epsilon_factor_without_centering(hkl) * (int) cen_ops.size();
  }

  static bool has_phase_shift(const Op::Tran& c, const Op::Miller& hkl) {
    return (hkl[0] * c[0] + hkl[1] * c[1] + hkl[2] * c[2]) % Op::DEN != 0;
  }

  bool is_systematically_absent(const Op::Miller& hkl) const {
    for (auto i = cen_ops.begin() + 1; i != cen_ops.end(); ++i)
      if (has_phase_shift(*i, hkl))
        return true;
    Op::Miller denh = {{Op::DEN * hkl[0], Op::DEN * hkl[1], Op::DEN * hkl[2]}};
    for (auto op = sym_ops.begin() + 1; op != sym_ops.end(); ++op)
      if (op->apply_to_hkl_without_division(hkl) == denh) {
        for (const Op::Tran& c : cen_ops)
          if (has_phase_shift({{op->tran[0] + c[0],
                                op->tran[1] + c[1],
                                op->tran[2] + c[2]}}, hkl))
            return true;
      }
    return false;
  }

  void change_basis_impl(const Op& cob, const Op& inv) {
    if (sym_ops.empty() || cen_ops.empty())
      return;

    // Apply change-of-basis to sym_ops.
    // Ignore the first item in sym_ops -- it's identity.
    for (auto op = sym_ops.begin() + 1; op != sym_ops.end(); ++op)
      *op = cob.combine(*op).combine(inv).wrap();

    // The number of centering vectors may be different.
    // As an ad-hoc method (not proved to be robust) add lattice points
    // from a super-cell.
    int idet = inv.det_rot() / (Op::DEN * Op::DEN * Op::DEN);
    if (idet > 1) {
      std::vector<Op::Tran> new_cen_ops;
      new_cen_ops.reserve(cen_ops.size() * idet * idet * idet);
      for (int i = 0; i < idet; ++i)
        for (int j = 0; j < idet; ++j)
          for (int k = 0; k < idet; ++k)
            for (Op::Tran& cen : cen_ops)
              new_cen_ops.push_back({i * Op::DEN + cen[0],
                                     j * Op::DEN + cen[1],
                                     k * Op::DEN + cen[2]});
      cen_ops.swap(new_cen_ops);
    }

    // Apply change-of-basis to centering vectors
    Op cvec = Op::identity();
    for (auto tr = cen_ops.begin() + 1; tr != cen_ops.end(); ++tr) {
      cvec.tran = *tr;
      *tr = cob.combine(cvec).combine(inv).wrap().tran;
    }

    // Remove redundant centering vectors.
    for (int i = static_cast<int>(cen_ops.size()) - 1; i > 0; --i)
      for (int j = i - 1; j >= 0; --j)
        if (cen_ops[i] == cen_ops[j]) {
          cen_ops.erase(cen_ops.begin() + i);
          break;
        }
  }

  void change_basis_forward(const Op& cob) { change_basis_impl(cob, cob.inverse()); }
  void change_basis_backward(const Op& inv) { change_basis_impl(inv.inverse(), inv); }

  std::vector<Op> all_ops_sorted() const {
    std::vector<Op> ops;
    ops.reserve(sym_ops.size() * cen_ops.size());
    for (const Op& so : sym_ops)
      for (const Op::Tran& co : cen_ops)
        ops.push_back(so.add_centering(co));
    std::sort(ops.begin(), ops.end());
    return ops;
  }

  Op get_op(int n) const {
    int n_cen = n / (int) sym_ops.size();
    int n_sym = n % (int) sym_ops.size();
    return sym_ops.at(n_sym).add_centering(cen_ops.at(n_cen));
  }

  bool is_same_as(const GroupOps& other) const {
    if (cen_ops.size() != other.cen_ops.size() ||
        sym_ops.size() != other.sym_ops.size())
      return false;
    return all_ops_sorted() == other.all_ops_sorted();
  }

  bool has_same_centring(const GroupOps& other) const {
    if (cen_ops.size() != other.cen_ops.size())
      return false;
    if (std::is_sorted(cen_ops.begin(), cen_ops.end()) &&
        std::is_sorted(other.cen_ops.begin(), other.cen_ops.end()))
      return cen_ops == other.cen_ops;
    std::vector<Op::Tran> v1 = cen_ops;
    std::vector<Op::Tran> v2 = other.cen_ops;
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    return v1 == v2;
  }

  bool has_same_rotations(const GroupOps& other) const {
    if (sym_ops.size() != other.sym_ops.size())
      return false;
    auto sorted_rotations = [](const GroupOps& g) {
      std::vector<Op::Rot> r(g.sym_ops.size());
      for (size_t i = 0; i != r.size(); ++i)
        r[i] = g.sym_ops[i].rot;
      std::sort(r.begin(), r.end());
      return r;
    };
    return sorted_rotations(*this) == sorted_rotations(other);
  }

  // minimal multiplicity for real-space grid in each direction
  // examples: 1,2,1 for P21, 1,1,6 for P61
  std::array<int, 3> find_grid_factors() const {
    const int T = Op::DEN;
    int r[3] = {T, T, T};
    for (Op op : *this)
      for (int i = 0; i != 3; ++i)
        if (op.tran[i] != 0 && op.tran[i] < r[i])
          r[i] = op.tran[i];
    return {T / r[0], T / r[1], T / r[2]};
  }

  bool are_directions_symmetry_related(int u, int v) const {
    for (const Op& op : sym_ops)
      if (op.rot[u][v] != 0)
        return true;
    return false;
  }

  struct Iter {
    const GroupOps& gops;
    int n_sym, n_cen;
    void operator++() {
      if (++n_sym == (int) gops.sym_ops.size()) {
        ++n_cen;
        n_sym = 0;
      }
    }
    Op operator*() const {
      return gops.sym_ops.at(n_sym).translated(gops.cen_ops.at(n_cen)).wrap();
    }
    bool operator==(const Iter& other) const {
      return n_sym == other.n_sym && n_cen == other.n_cen;
    }
    bool operator!=(const Iter& other) const { return !(*this == other); }
  };

  Iter begin() const { return {*this, 0, 0}; };
  Iter end() const { return {*this, 0, (int) cen_ops.size()}; };
};

inline void GroupOps::add_missing_elements() {
  // We always keep identity as sym_ops[0].
  if (sym_ops.empty() || sym_ops[0] != Op::identity())
    fail("oops");
  if (sym_ops.size() == 1)
    return;
  auto check_size = [&]() {
    if (sym_ops.size() > 1023)
      fail("1000+ elements in the group should not happen");
  };
  // Below we assume that all centring vectors are already known (in cen_ops)
  // so when checking for a new element we compare only the 3x3 matrix.
  // Dimino's algorithm. https://physics.stackexchange.com/a/351400/95713
  std::vector<Op> gen(sym_ops.begin() + 1, sym_ops.end());
  sym_ops.resize(2);
  const Op::Rot idrot = Op::identity().rot;
  for (Op g = sym_ops[1] * sym_ops[1]; g.rot != idrot; g *= sym_ops[1]) {
    sym_ops.push_back(g);
    check_size();
  }
  for (size_t i = 1; i < gen.size(); ++i) {
    std::vector<Op> coset_repr(1, Op::identity());
    size_t init_size = sym_ops.size();
    for (;;) {
      size_t len = coset_repr.size();
      for (size_t j = 0; j != len; ++j) {
        for (size_t n = 0; n != i + 1; ++n) {
          Op sg = gen[n] * coset_repr[j];
          if (find_by_rotation(sg.rot) == nullptr) {
            sym_ops.push_back(sg);
            for (size_t k = 1; k != init_size; ++k)
              sym_ops.push_back(sg * sym_ops[k]);
            coset_repr.push_back(sg);
          }
        }
      }
      if (len == coset_repr.size())
        break;
      check_size();
    }
  }
}

// Create GroupOps from Ops by separating centering vectors
inline GroupOps split_centering_vectors(const std::vector<Op>& ops) {
  const Op identity = Op::identity();
  GroupOps go;
  go.sym_ops.push_back(identity);
  for (const Op& op : ops)
    if (Op* old_op = go.find_by_rotation(op.rot)) {
      if (op.rot == identity.rot)  // pure shift
        go.cen_ops.push_back(op.tran);
      if (op.tran == identity.tran)  // or rather |op.tran| < |old_op->tran| ?
        old_op->tran = op.tran;
    } else {
      go.sym_ops.push_back(op);
    }
  return go;
}

// INTERPRETING HALL SYMBOLS
// based on both ITfC vol.B ch.1.4 (2010)
// and http://cci.lbl.gov/sginfo/hall_symbols.html

// matrices for Nz from Table 3 and 4 from hall_symbols.html
inline Op::Rot hall_rotation_z(int N) {
  constexpr int d = Op::DEN;
  switch (N) {
    case 1: return {d,0,0,  0,d,0,  0,0,d};
    case 2: return {-d,0,0, 0,-d,0, 0,0,d};
    case 3: return {0,-d,0, d,-d,0, 0,0,d};
    case 4: return {0,-d,0, d,0,0,  0,0,d};
    case 6: return {d,-d,0, d,0,0,  0,0,d};
    case '\'': return {0,-d,0, -d,0,0, 0,0,-d};
    case '"':  return {0,d,0,   d,0,0, 0,0,-d};
    case '*':  return {0,0,d,   d,0,0, 0,d,0};
    default: fail("incorrect axis definition");
  }
}
inline Op::Tran hall_translation_from_symbol(char symbol) {
  constexpr int h = Op::DEN / 2;
  constexpr int q = Op::DEN / 4;
  switch (symbol) {
    case 'a': return {h, 0, 0};
    case 'b': return {0, h, 0};
    case 'c': return {0, 0, h};
    case 'n': return {h, h, h};
    case 'u': return {q, 0, 0};
    case 'v': return {0, q, 0};
    case 'w': return {0, 0, q};
    case 'd': return {q, q, q};
    default: fail(std::string("unknown symbol: ") + symbol);
  }
}

inline Op hall_matrix_symbol(const char* start, const char* end,
                             int pos, int& prev) {
  Op op = Op::identity();
  bool neg = (*start == '-');
  const char* p = (neg ? start + 1 : start);
  if (*p < '1' || *p == '5' || *p > '6')
    fail("wrong n-fold order notation: " + std::string(start, end));
  int N = *p++ - '0';
  int fractional_tran = 0;
  char principal_axis = '\0';
  char diagonal_axis = '\0';
  for (; p < end; ++p) {
    if (*p >= '1' && *p <= '5') {
      if (fractional_tran != '\0')
        fail("two numeric subscripts");
      fractional_tran = *p - '0';
    } else if (*p == '\'' || *p == '"' || *p == '*') {
      if (N != (*p == '*' ? 3 : 2))
        fail("wrong symbol: " + std::string(start, end));
      diagonal_axis = *p;
    } else if (*p == 'x' || *p == 'y' || *p == 'z') {
      principal_axis = *p;
    } else {
      op.translate(hall_translation_from_symbol(*p));
    }
  }
  // fill in implicit values
  if (!principal_axis && !diagonal_axis) {
    if (pos == 1) {
      principal_axis = 'z';
    } else if (pos == 2 && N == 2) {
      if (prev == 2 || prev == 4)
        principal_axis = 'x';
      else if (prev == 3 || prev == 6)
        diagonal_axis = '\'';
    } else if (pos == 3 && N == 3) {
      diagonal_axis = '*';
    } else if (N != 1) {
      fail("missing axis");
    }
  }
  // get the operation
  op.rot = hall_rotation_z(diagonal_axis ? diagonal_axis : N);
  if (neg)
    op.rot = op.negated_rot();
  auto alter_order = [](const Op::Rot& r, int i, int j, int k) {
    return Op::Rot{ r[i][i], r[i][j], r[i][k],
                    r[j][i], r[j][j], r[j][k],
                    r[k][i], r[k][j], r[k][k] };
  };
  if (principal_axis == 'x')
    op.rot = alter_order(op.rot, 2, 0, 1);
  else if (principal_axis == 'y')
    op.rot = alter_order(op.rot, 1, 2, 0);
  if (fractional_tran)
    op.tran[principal_axis - 'x'] += Op::DEN / N * fractional_tran;
  prev = N;
  return op;
}

// Parses either short (0 0 1) or long notation (x,y,z+1/12)
// but without multpliers (such as 1/2x) to keep things simple for now.
inline Op parse_hall_change_of_basis(const char* start, const char* end) {
  if (std::memchr(start, ',', end - start) != nullptr) // long symbol
    return parse_triplet(std::string(start, end));
  // short symbol (0 0 1)
  Op cob = Op::identity();
  char* endptr;
  for (int i = 0; i != 3; ++i) {
    cob.tran[i] = std::strtol(start, &endptr, 10) % 12 * (Op::DEN / 12);
    start = endptr;
  }
  if (endptr != end)
    fail("unexpected change-of-basis format: " + std::string(start, end));
  return cob;
}

inline GroupOps generators_from_hall(const char* hall) {
  auto find_blank = [](const char* p) {
    while (*p != '\0' && *p != ' ' && *p != '\t' && *p != '_') // '_' == ' '
      ++p;
    return p;
  };
  if (hall == nullptr)
    fail("null");
  hall = impl::skip_blank(hall);
  GroupOps ops;
  ops.sym_ops.emplace_back(Op::identity());
  bool centrosym = (hall[0] == '-');
  const char* lat = impl::skip_blank(centrosym ? hall + 1 : hall);
  if (!lat)
    fail("not a hall symbol: " + std::string(hall));
  ops.cen_ops = centring_vectors(*lat);
  int counter = 0;
  int prev = 0;
  const char* part = impl::skip_blank(lat + 1);
  while (*part != '\0' && *part != '(') {
    const char* space = find_blank(part);
    ++counter;
    if (part[0] != '1' || (part[1] != ' ' && part[1] != '\0')) {
      Op op = hall_matrix_symbol(part, space, counter, prev);
      ops.sym_ops.emplace_back(op);
    }
    part = impl::skip_blank(space);
  }
  if (centrosym)
    ops.sym_ops.emplace_back(Op::identity().negated());
  if (*part == '(') {
    const char* rb = std::strchr(part, ')');
    if (!rb)
      fail("missing ')': " + std::string(hall));
    if (ops.sym_ops.empty())
      fail("misplaced translation: " + std::string(hall));
    ops.change_basis_forward(parse_hall_change_of_basis(part + 1, rb));

    if (*impl::skip_blank(find_blank(rb + 1)) != '\0')
      fail("unexpected characters after ')': " + std::string(hall));
  }
  return ops;
}

inline GroupOps symops_from_hall(const char* hall) {
  GroupOps ops = generators_from_hall(hall);
  ops.add_missing_elements();
  return ops;
}

// CRYSTAL SYSTEMS, POINT GROUPS AND LAUE CLASSES

enum class CrystalSystem : unsigned char {
  Triclinic=0, Monoclinic, Orthorhombic, Tetragonal, Trigonal, Hexagonal, Cubic
};

inline const char* crystal_system_str(CrystalSystem system) {
  static const char* names[7] = {
    "triclinic", "monoclinic", "orthorhombic", "tetragonal",
    "trigonal", "hexagonal", "cubic"
  };
  return names[static_cast<int>(system)];
}

enum class PointGroup : unsigned char {
  C1=0, Ci, C2, Cs, C2h, D2, C2v, D2h, C4, S4, C4h, D4, C4v, D2d, D4h, C3,
  C3i, D3, C3v, D3d, C6, C3h, C6h, D6, C6v, D3h, D6h, T, Th, O, Td, Oh
};

inline const char* point_group_hm(PointGroup pg) {
  static const char hm_pointgroup_names[32][6] = {
    "1", "-1", "2", "m", "2/m", "222", "mm2", "mmm",
    "4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm", "3",
    "-3", "32", "3m", "-3m", "6", "-6", "6/m", "622",
    "6mm", "-62m", "6/mmm", "23", "m-3", "432", "-43m", "m-3m",
  };
  return hm_pointgroup_names[static_cast<int>(pg)];
}

// http://reference.iucr.org/dictionary/Laue_class
enum class Laue : unsigned char {
  L1=0, L2m, Lmmm, L4m, L4mmm, L3, L3m, L6m, L6mmm, Lm3, Lm3m
};

inline Laue pointgroup_to_laue(PointGroup pg) {
  static const Laue laue[32] = {
    Laue::L1, Laue::L1,
    Laue::L2m, Laue::L2m, Laue::L2m,
    Laue::Lmmm, Laue::Lmmm, Laue::Lmmm,
    Laue::L4m, Laue::L4m, Laue::L4m,
    Laue::L4mmm, Laue::L4mmm, Laue::L4mmm, Laue::L4mmm,
    Laue::L3, Laue::L3,
    Laue::L3m, Laue::L3m, Laue::L3m,
    Laue::L6m, Laue::L6m, Laue::L6m,
    Laue::L6mmm, Laue::L6mmm, Laue::L6mmm, Laue::L6mmm,
    Laue::Lm3, Laue::Lm3,
    Laue::Lm3m, Laue::Lm3m, Laue::Lm3m,
  };
  return laue[static_cast<int>(pg)];
}

// return centrosymmetric pointgroup from the Laue class
inline PointGroup laue_to_pointgroup(Laue laue) {
  static const PointGroup pg[11] = {
    PointGroup::Ci, PointGroup::C2h, PointGroup::D2h, PointGroup::C4h,
    PointGroup::D4h, PointGroup::C3i, PointGroup::D3d, PointGroup::C6h,
    PointGroup::D6h, PointGroup::Th, PointGroup::Oh
  };
  return pg[static_cast<int>(laue)];
}

inline const char* laue_class_str(Laue laue) {
  return point_group_hm(laue_to_pointgroup(laue));
}

inline CrystalSystem crystal_system(Laue laue) {
  static const CrystalSystem crystal_systems[11] = {
    CrystalSystem::Triclinic,
    CrystalSystem::Monoclinic,
    CrystalSystem::Orthorhombic,
    CrystalSystem::Tetragonal, CrystalSystem::Tetragonal,
    CrystalSystem::Trigonal,   CrystalSystem::Trigonal,
    CrystalSystem::Hexagonal,  CrystalSystem::Hexagonal,
    CrystalSystem::Cubic,      CrystalSystem::Cubic
  };
  return crystal_systems[static_cast<int>(laue)];
}

inline CrystalSystem crystal_system(PointGroup pg) {
  return crystal_system(pointgroup_to_laue(pg));
}

inline unsigned char point_group_index_and_category(int space_group_number) {
  // 0x20=Sohncke, 0x40=enantiomorphic, 0x80=symmorphic
  enum : unsigned char { S=0x20, E=(0x20|0x40), Y=0x80, Z=(0x20|0x80) };
  static unsigned char indices[230] = {
     0|Z,  1|Y,  2|Z,  2|S,  2|Z,  3|Y,  3,    3|Y,  3,    4|Y,  // 1-10
     4,    4|Y,  4,    4,    4,    5|Z,  5|S,  5|S,  5|S,  5|S,  // 11-20
     5|Z,  5|Z,  5|Z,  5|S,  6|Y,  6,    6,    6,    6,    6,    // 21-30
     6,    6,    6,    6,    6|Y,  6,    6,    6|Y,  6,    6,    // 31-40
     6,    6|Y,  6,    6|Y,  6,    6,    7|Y,  7,    7,    7,    // 41-50
     7,    7,    7,    7,    7,    7,    7,    7,    7,    7,    // 51-60
     7,    7,    7,    7,    7|Y,  7,    7,    7,    7|Y,  7,    // 61-70
     7|Y,  7,    7,    7,    8|Z,  8|E,  8|S,  8|E,  8|Z,  8|S,  // 71-80
     9|Y,  9|Y, 10|Y, 10,   10,   10,   10|Y, 10,   11|Z, 11|S,  // 81-90
    11|E, 11|E, 11|S, 11|S, 11|E, 11|E, 11|Z, 11|S, 12|Y, 12,    // 91-100
    12,   12,   12,   12,   12,   12,   12|Y, 12,   12,   12,    // 101-110
    13|Y, 13,   13,   13,   13|Y, 13,   13,   13,   13|Y, 13,    // 111-120
    13|Y, 13,   14|Y, 14,   14,   14,   14,   14,   14,   14,    // 121-130
    14,   14,   14,   14,   14,   14,   14,   14,   14|Y, 14,    // 131-140
    14,   14,   15|Z, 15|E, 15|E, 15|Z, 16|Y, 16|Y, 17|Z, 17|Z,  // 141-150
    17|E, 17|E, 17|E, 17|E, 17|Z, 18|Y, 18|Y, 18,   18,   18|Y,  // 151-160
    18,   19|Y, 19,   19|Y, 19,   19|Y, 19,   20|Z, 20|E, 20|E,  // 161-170
    20|E, 20|E, 20|S, 21|Y, 22|Y, 22,   23|Z, 23|E, 23|E, 23|E,  // 171-180
    23|E, 23|S, 24|Y, 24,   24,   24,   25|Y, 25,   25|Y, 25,    // 181-190
    26|Y, 26,   26,   26,   27|Z, 27|Z, 27|Z, 27|S, 27|S, 28|Y,  // 191-200
    28,   28|Y, 28,   28|Y, 28,   28,   29|Z, 29|S, 29|Z, 29|S,  // 201-210
    29|Z, 29|E, 29|E, 29|S, 30|Y, 30|Y, 30|Y, 30,   30,   30,    // 211-220
    31|Y, 31,   31,   31,   31|Y, 31,   31,   31,   31|Y, 31     // 221-230
  };
  return indices[space_group_number-1];
}

inline PointGroup point_group(int space_group_number) {
  auto n = point_group_index_and_category(space_group_number);
  return static_cast<PointGroup>(n & 0x1f);
}

// true for 65 Sohncke (non-enantiogenic) space groups
inline bool is_sohncke(int space_group_number) {
  return (point_group_index_and_category(space_group_number) & 0x20) != 0;
}

// true for 22 space groups (11 enantiomorphic pairs)
inline bool is_enantiomorphic(int space_group_number) {
  return (point_group_index_and_category(space_group_number) & 0x40) != 0;
}

// true for 73 space groups
inline bool is_symmorphic(int space_group_number) {
  return (point_group_index_and_category(space_group_number) & 0x80) != 0;
}

// Generated by tools/gen_sg_table.py.
inline const char* get_basisop(int basisop_idx) {
  static const char* basisops[49] = {
    "x,y,z",
    "z,x,y",
    "y,z,x",
    "z,y,-x",
    "x,y,-x+z",
    "-x,z,y",
    "-x+z,x,y",
    "y,-x,z",
    "y,-x+z,x",
    "x-z,y,z",
    "z,x-z,y",
    "y,z,x-z",
    "z,y,-x+z",
    "x+z,y,-x",
    "x+1/4,y+1/4,z",
    "-x+z,z,y",
    "-x,x+z,y",
    "y,-x+z,z",
    "y,-x,x+z",
    "x+1/4,y-1/4,z",
    "x-1/4,y-1/4,z-1/4",
    "x-1/4,y-1/4,z",
    "z,x-1/4,y-1/4",
    "y-1/4,z,x-1/4",
    "x-1/2,y-1/4,z+1/4",
    "z+1/4,x-1/2,y-1/4",
    "y-1/4,z+1/4,x-1/2",
    "x+1/8,y+1/8,z+1/8",
    "x+1/4,y-1/4,z+1/4",
    "x-1/4,y+1/4,z",
    "x+1/4,y+1/4,z+1/4",
    "x,y+1/4,z+1/8",
    "x-1/4,y+1/4,z+1/4",
    "x-1/4,y+1/4,z-1/4",
    "x-1/2,y+1/4,z+1/8",
    "x-1/2,y+1/4,z-3/8",
    "-y+z,x+z,-x+y+z",
    "x-1/8,y-1/8,z-1/8",
    "x+1/4,y+1/4,-x+z-1/4",
    "x+1/4,y,z",
    "x,y,z+1/4",
    "-x,-y/2+z/2,y/2+z/2",
    "-x/2+z/2,-y,x/2+z/2",
    "x/2+y/2,x/2-y/2,-z",
    "y/2+z/2,x/2+z/2,x/2+y/2",
    "-x/2+y/2+z/2,x/2-y/2+z/2,x/2+y/2-z/2",
    "-x/2+z,x/2,y",
    "x-z/2,y,z/2",
    "x/2+y/2,-x/2+y/2,z",
  };
  return basisops[basisop_idx];
}

// Returns a change-of-basis operator for centred -> primitive transformation.
// The same operator as inverse of z2p_op in sgtbx.
inline Op::Rot centred_to_primitive(char centring_type) {
  constexpr int D = Op::DEN;
  constexpr int H = Op::DEN / 2;
  constexpr int T = Op::DEN / 3;
  switch (centring_type) {
    case 'P': return {D,0,0, 0,D,0, 0,0,D};
    case 'A': return {-D,0,0, 0,-H,H, 0,H,H};
    case 'B': return {-H,0,H, 0,-D,0, H,0,H};
    case 'C': return {H,H,0, H,-H,0, 0,0,-D};
    case 'I': return {-H,H,H, H,-H,H, H,H,-H};
    case 'R': return {2*T,-T,-T, T,T,-2*T, T,T,T};
    case 'H': return {2*T,-T,0, T,T,0, 0,0,D};  // not used normally
    case 'F': return {0,H,H, H,0,H, H,H,0};
    default: fail("not a centring type: ", centring_type);
  }
}


// LIST OF CRYSTALLOGRAPHIC SPACE GROUPS

struct SpaceGroup { // typically 44 bytes
  int number;
  int ccp4;
  char hm[11];  // Hermann-Mauguin (international) notation
  char ext;
  char qualifier[5];
  char hall[15];
  int basisop_idx;

  std::string xhm() const {
    std::string ret = hm;
    if (ext) {
      ret += ':';
      ret += ext;
    }
    return ret;
  }

  char centring_type() const { return ext == 'R' ? 'P' : hm[0]; }

  // (old) CCP4 spacegroup names start with H for hexagonal setting
  char ccp4_lattice_type() const { return ext == 'H' ? 'H' : hm[0]; }

  // P 1 2 1 -> P2, but P 1 1 2 -> P112. R 3:H -> H3.
  std::string short_name() const {
    std::string s(hm);
    size_t len = s.size();
    if (len > 6 && s[2] == '1' && s[len - 2] == ' ' && s[len - 1] == '1')
      s = s[0] + s.substr(4, len - 4 - 2);
    if (ext == 'H')
      s[0] = 'H';
    s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
    return s;
  }

  bool is_sohncke() const { return gemmi::is_sohncke(number); }
  bool is_enantiomorphic() const { return gemmi::is_enantiomorphic(number); }
  bool is_symmorphic() const { return gemmi::is_symmorphic(number); }
  PointGroup point_group() const { return gemmi::point_group(number); }
  const char* point_group_hm() const {
    return gemmi::point_group_hm(point_group());
  }
  Laue laue_class() const { return pointgroup_to_laue(point_group()); }
  const char* laue_str() const { return laue_class_str(laue_class()); }
  CrystalSystem crystal_system() const {
    return gemmi::crystal_system(point_group());
  }
  const char* crystal_system_str() const {
    return gemmi::crystal_system_str(crystal_system());
  }

  const char* basisop_str() const { return get_basisop(basisop_idx); }
  Op basisop() const { return parse_triplet(basisop_str()); }
  bool is_reference_setting() const { return basisop_idx == 0; }

  Op centred_to_primitive() const {
    return {gemmi::centred_to_primitive(centring_type()), {0,0,0}};
  }

  GroupOps operations() const { return symops_from_hall(hall); }
};

struct SpaceGroupAltName {
  char hm[11];
  char ext;
  int pos;
};

// the template here is only to substitute C++17 inline variables
// https://stackoverflow.com/questions/38043442/how-do-inline-variables-work
namespace impl {

template<class Dummy>
struct Tables_
{
  static const SpaceGroup main[557];
  static const SpaceGroupAltName alt_names[28];
  static const unsigned char ccp4_hkl_asu[230];
};

template<class Dummy>
const SpaceGroup Tables_<Dummy>::main[557] = {
  // This table was generated by tools/gen_sg_table.py.
  // First 530 entries in the same order as in SgInfo, sgtbx and ITB.
  // Note: spacegroup 68 has three duplicates with different H-M names.
  {  1,    1, "P 1"       ,   0,     "", "P 1"           , 0 }, //   0
  {  2,    2, "P -1"      ,   0,     "", "-P 1"          , 0 }, //   1
  {  3,    3, "P 1 2 1"   ,   0,    "b", "P 2y"          , 0 }, //   2
  {  3, 1003, "P 1 1 2"   ,   0,    "c", "P 2"           , 1 }, //   3
  {  3,    0, "P 2 1 1"   ,   0,    "a", "P 2x"          , 2 }, //   4
  {  4,    4, "P 1 21 1"  ,   0,    "b", "P 2yb"         , 0 }, //   5
  {  4, 1004, "P 1 1 21"  ,   0,    "c", "P 2c"          , 1 }, //   6
  {  4,    0, "P 21 1 1"  ,   0,    "a", "P 2xa"         , 2 }, //   7
  {  5,    5, "C 1 2 1"   ,   0,   "b1", "C 2y"          , 0 }, //   8
  {  5, 2005, "A 1 2 1"   ,   0,   "b2", "A 2y"          , 3 }, //   9
  {  5, 4005, "I 1 2 1"   ,   0,   "b3", "I 2y"          , 4 }, //  10
  {  5,    0, "A 1 1 2"   ,   0,   "c1", "A 2"           , 1 }, //  11
  {  5, 1005, "B 1 1 2"   ,   0,   "c2", "B 2"           , 5 }, //  12
  {  5,    0, "I 1 1 2"   ,   0,   "c3", "I 2"           , 6 }, //  13
  {  5,    0, "B 2 1 1"   ,   0,   "a1", "B 2x"          , 2 }, //  14
  {  5,    0, "C 2 1 1"   ,   0,   "a2", "C 2x"          , 7 }, //  15
  {  5,    0, "I 2 1 1"   ,   0,   "a3", "I 2x"          , 8 }, //  16
  {  6,    6, "P 1 m 1"   ,   0,    "b", "P -2y"         , 0 }, //  17
  {  6, 1006, "P 1 1 m"   ,   0,    "c", "P -2"          , 1 }, //  18
  {  6,    0, "P m 1 1"   ,   0,    "a", "P -2x"         , 2 }, //  19
  {  7,    7, "P 1 c 1"   ,   0,   "b1", "P -2yc"        , 0 }, //  20
  {  7,    0, "P 1 n 1"   ,   0,   "b2", "P -2yac"       , 9 }, //  21
  {  7,    0, "P 1 a 1"   ,   0,   "b3", "P -2ya"        , 3 }, //  22
  {  7,    0, "P 1 1 a"   ,   0,   "c1", "P -2a"         , 1 }, //  23
  {  7,    0, "P 1 1 n"   ,   0,   "c2", "P -2ab"        , 10}, //  24
  {  7, 1007, "P 1 1 b"   ,   0,   "c3", "P -2b"         , 5 }, //  25
  {  7,    0, "P b 1 1"   ,   0,   "a1", "P -2xb"        , 2 }, //  26
  {  7,    0, "P n 1 1"   ,   0,   "a2", "P -2xbc"       , 11}, //  27
  {  7,    0, "P c 1 1"   ,   0,   "a3", "P -2xc"        , 7 }, //  28
  {  8,    8, "C 1 m 1"   ,   0,   "b1", "C -2y"         , 0 }, //  29
  {  8,    0, "A 1 m 1"   ,   0,   "b2", "A -2y"         , 3 }, //  30
  {  8,    0, "I 1 m 1"   ,   0,   "b3", "I -2y"         , 4 }, //  31
  {  8,    0, "A 1 1 m"   ,   0,   "c1", "A -2"          , 1 }, //  32
  {  8, 1008, "B 1 1 m"   ,   0,   "c2", "B -2"          , 5 }, //  33
  {  8,    0, "I 1 1 m"   ,   0,   "c3", "I -2"          , 6 }, //  34
  {  8,    0, "B m 1 1"   ,   0,   "a1", "B -2x"         , 2 }, //  35
  {  8,    0, "C m 1 1"   ,   0,   "a2", "C -2x"         , 7 }, //  36
  {  8,    0, "I m 1 1"   ,   0,   "a3", "I -2x"         , 8 }, //  37
  {  9,    9, "C 1 c 1"   ,   0,   "b1", "C -2yc"        , 0 }, //  38
  {  9,    0, "A 1 n 1"   ,   0,   "b2", "A -2yab"       , 12}, //  39
  {  9,    0, "I 1 a 1"   ,   0,   "b3", "I -2ya"        , 13}, //  40
  {  9,    0, "A 1 a 1"   ,   0,  "-b1", "A -2ya"        , 3 }, //  41
  {  9,    0, "C 1 n 1"   ,   0,  "-b2", "C -2yac"       , 14}, //  42
  {  9,    0, "I 1 c 1"   ,   0,  "-b3", "I -2yc"        , 4 }, //  43
  {  9,    0, "A 1 1 a"   ,   0,   "c1", "A -2a"         , 1 }, //  44
  {  9,    0, "B 1 1 n"   ,   0,   "c2", "B -2ab"        , 15}, //  45
  {  9,    0, "I 1 1 b"   ,   0,   "c3", "I -2b"         , 16}, //  46
  {  9, 1009, "B 1 1 b"   ,   0,  "-c1", "B -2b"         , 5 }, //  47
  {  9,    0, "A 1 1 n"   ,   0,  "-c2", "A -2ab"        , 10}, //  48
  {  9,    0, "I 1 1 a"   ,   0,  "-c3", "I -2a"         , 6 }, //  49
  {  9,    0, "B b 1 1"   ,   0,   "a1", "B -2xb"        , 2 }, //  50
  {  9,    0, "C n 1 1"   ,   0,   "a2", "C -2xac"       , 17}, //  51
  {  9,    0, "I c 1 1"   ,   0,   "a3", "I -2xc"        , 18}, //  52
  {  9,    0, "C c 1 1"   ,   0,  "-a1", "C -2xc"        , 7 }, //  53
  {  9,    0, "B n 1 1"   ,   0,  "-a2", "B -2xab"       , 11}, //  54
  {  9,    0, "I b 1 1"   ,   0,  "-a3", "I -2xb"        , 8 }, //  55
  { 10,   10, "P 1 2/m 1" ,   0,    "b", "-P 2y"         , 0 }, //  56
  { 10, 1010, "P 1 1 2/m" ,   0,    "c", "-P 2"          , 1 }, //  57
  { 10,    0, "P 2/m 1 1" ,   0,    "a", "-P 2x"         , 2 }, //  58
  { 11,   11, "P 1 21/m 1",   0,    "b", "-P 2yb"        , 0 }, //  59
  { 11, 1011, "P 1 1 21/m",   0,    "c", "-P 2c"         , 1 }, //  60
  { 11,    0, "P 21/m 1 1",   0,    "a", "-P 2xa"        , 2 }, //  61
  { 12,   12, "C 1 2/m 1" ,   0,   "b1", "-C 2y"         , 0 }, //  62
  { 12,    0, "A 1 2/m 1" ,   0,   "b2", "-A 2y"         , 3 }, //  63
  { 12,    0, "I 1 2/m 1" ,   0,   "b3", "-I 2y"         , 4 }, //  64
  { 12,    0, "A 1 1 2/m" ,   0,   "c1", "-A 2"          , 1 }, //  65
  { 12, 1012, "B 1 1 2/m" ,   0,   "c2", "-B 2"          , 5 }, //  66
  { 12,    0, "I 1 1 2/m" ,   0,   "c3", "-I 2"          , 6 }, //  67
  { 12,    0, "B 2/m 1 1" ,   0,   "a1", "-B 2x"         , 2 }, //  68
  { 12,    0, "C 2/m 1 1" ,   0,   "a2", "-C 2x"         , 7 }, //  69
  { 12,    0, "I 2/m 1 1" ,   0,   "a3", "-I 2x"         , 8 }, //  70
  { 13,   13, "P 1 2/c 1" ,   0,   "b1", "-P 2yc"        , 0 }, //  71
  { 13,    0, "P 1 2/n 1" ,   0,   "b2", "-P 2yac"       , 9 }, //  72
  { 13,    0, "P 1 2/a 1" ,   0,   "b3", "-P 2ya"        , 3 }, //  73
  { 13,    0, "P 1 1 2/a" ,   0,   "c1", "-P 2a"         , 1 }, //  74
  { 13,    0, "P 1 1 2/n" ,   0,   "c2", "-P 2ab"        , 10}, //  75
  { 13, 1013, "P 1 1 2/b" ,   0,   "c3", "-P 2b"         , 5 }, //  76
  { 13,    0, "P 2/b 1 1" ,   0,   "a1", "-P 2xb"        , 2 }, //  77
  { 13,    0, "P 2/n 1 1" ,   0,   "a2", "-P 2xbc"       , 11}, //  78
  { 13,    0, "P 2/c 1 1" ,   0,   "a3", "-P 2xc"        , 7 }, //  79
  { 14,   14, "P 1 21/c 1",   0,   "b1", "-P 2ybc"       , 0 }, //  80
  { 14, 2014, "P 1 21/n 1",   0,   "b2", "-P 2yn"        , 9 }, //  81
  { 14, 3014, "P 1 21/a 1",   0,   "b3", "-P 2yab"       , 3 }, //  82
  { 14,    0, "P 1 1 21/a",   0,   "c1", "-P 2ac"        , 1 }, //  83
  { 14,    0, "P 1 1 21/n",   0,   "c2", "-P 2n"         , 10}, //  84
  { 14, 1014, "P 1 1 21/b",   0,   "c3", "-P 2bc"        , 5 }, //  85
  { 14,    0, "P 21/b 1 1",   0,   "a1", "-P 2xab"       , 2 }, //  86
  { 14,    0, "P 21/n 1 1",   0,   "a2", "-P 2xn"        , 11}, //  87
  { 14,    0, "P 21/c 1 1",   0,   "a3", "-P 2xac"       , 7 }, //  88
  { 15,   15, "C 1 2/c 1" ,   0,   "b1", "-C 2yc"        , 0 }, //  89
  { 15,    0, "A 1 2/n 1" ,   0,   "b2", "-A 2yab"       , 12}, //  90
  { 15,    0, "I 1 2/a 1" ,   0,   "b3", "-I 2ya"        , 13}, //  91
  { 15,    0, "A 1 2/a 1" ,   0,  "-b1", "-A 2ya"        , 3 }, //  92
  { 15,    0, "C 1 2/n 1" ,   0,  "-b2", "-C 2yac"       , 19}, //  93
  { 15,    0, "I 1 2/c 1" ,   0,  "-b3", "-I 2yc"        , 4 }, //  94
  { 15,    0, "A 1 1 2/a" ,   0,   "c1", "-A 2a"         , 1 }, //  95
  { 15,    0, "B 1 1 2/n" ,   0,   "c2", "-B 2ab"        , 15}, //  96
  { 15,    0, "I 1 1 2/b" ,   0,   "c3", "-I 2b"         , 16}, //  97
  { 15, 1015, "B 1 1 2/b" ,   0,  "-c1", "-B 2b"         , 5 }, //  98
  { 15,    0, "A 1 1 2/n" ,   0,  "-c2", "-A 2ab"        , 10}, //  99
  { 15,    0, "I 1 1 2/a" ,   0,  "-c3", "-I 2a"         , 6 }, // 100
  { 15,    0, "B 2/b 1 1" ,   0,   "a1", "-B 2xb"        , 2 }, // 101
  { 15,    0, "C 2/n 1 1" ,   0,   "a2", "-C 2xac"       , 17}, // 102
  { 15,    0, "I 2/c 1 1" ,   0,   "a3", "-I 2xc"        , 18}, // 103
  { 15,    0, "C 2/c 1 1" ,   0,  "-a1", "-C 2xc"        , 7 }, // 104
  { 15,    0, "B 2/n 1 1" ,   0,  "-a2", "-B 2xab"       , 11}, // 105
  { 15,    0, "I 2/b 1 1" ,   0,  "-a3", "-I 2xb"        , 8 }, // 106
  { 16,   16, "P 2 2 2"   ,   0,     "", "P 2 2"         , 0 }, // 107
  { 17,   17, "P 2 2 21"  ,   0,     "", "P 2c 2"        , 0 }, // 108
  { 17, 1017, "P 21 2 2"  ,   0,  "cab", "P 2a 2a"       , 1 }, // 109
  { 17, 2017, "P 2 21 2"  ,   0,  "bca", "P 2 2b"        , 2 }, // 110
  { 18,   18, "P 21 21 2" ,   0,     "", "P 2 2ab"       , 0 }, // 111
  { 18, 3018, "P 2 21 21" ,   0,  "cab", "P 2bc 2"       , 1 }, // 112
  { 18, 2018, "P 21 2 21" ,   0,  "bca", "P 2ac 2ac"     , 2 }, // 113
  { 19,   19, "P 21 21 21",   0,     "", "P 2ac 2ab"     , 0 }, // 114
  { 20,   20, "C 2 2 21"  ,   0,     "", "C 2c 2"        , 0 }, // 115
  { 20,    0, "A 21 2 2"  ,   0,  "cab", "A 2a 2a"       , 1 }, // 116
  { 20,    0, "B 2 21 2"  ,   0,  "bca", "B 2 2b"        , 2 }, // 117
  { 21,   21, "C 2 2 2"   ,   0,     "", "C 2 2"         , 0 }, // 118
  { 21,    0, "A 2 2 2"   ,   0,  "cab", "A 2 2"         , 1 }, // 119
  { 21,    0, "B 2 2 2"   ,   0,  "bca", "B 2 2"         , 2 }, // 120
  { 22,   22, "F 2 2 2"   ,   0,     "", "F 2 2"         , 0 }, // 121
  { 23,   23, "I 2 2 2"   ,   0,     "", "I 2 2"         , 0 }, // 122
  { 24,   24, "I 21 21 21",   0,     "", "I 2b 2c"       , 0 }, // 123
  { 25,   25, "P m m 2"   ,   0,     "", "P 2 -2"        , 0 }, // 124
  { 25,    0, "P 2 m m"   ,   0,  "cab", "P -2 2"        , 1 }, // 125
  { 25,    0, "P m 2 m"   ,   0,  "bca", "P -2 -2"       , 2 }, // 126
  { 26,   26, "P m c 21"  ,   0,     "", "P 2c -2"       , 0 }, // 127
  { 26,    0, "P c m 21"  ,   0, "ba-c", "P 2c -2c"      , 7 }, // 128
  { 26,    0, "P 21 m a"  ,   0,  "cab", "P -2a 2a"      , 1 }, // 129
  { 26,    0, "P 21 a m"  ,   0, "-cba", "P -2 2a"       , 3 }, // 130
  { 26,    0, "P b 21 m"  ,   0,  "bca", "P -2 -2b"      , 2 }, // 131
  { 26,    0, "P m 21 b"  ,   0, "a-cb", "P -2b -2"      , 5 }, // 132
  { 27,   27, "P c c 2"   ,   0,     "", "P 2 -2c"       , 0 }, // 133
  { 27,    0, "P 2 a a"   ,   0,  "cab", "P -2a 2"       , 1 }, // 134
  { 27,    0, "P b 2 b"   ,   0,  "bca", "P -2b -2b"     , 2 }, // 135
  { 28,   28, "P m a 2"   ,   0,     "", "P 2 -2a"       , 0 }, // 136
  { 28,    0, "P b m 2"   ,   0, "ba-c", "P 2 -2b"       , 7 }, // 137
  { 28,    0, "P 2 m b"   ,   0,  "cab", "P -2b 2"       , 1 }, // 138
  { 28,    0, "P 2 c m"   ,   0, "-cba", "P -2c 2"       , 3 }, // 139
  { 28,    0, "P c 2 m"   ,   0,  "bca", "P -2c -2c"     , 2 }, // 140
  { 28,    0, "P m 2 a"   ,   0, "a-cb", "P -2a -2a"     , 5 }, // 141
  { 29,   29, "P c a 21"  ,   0,     "", "P 2c -2ac"     , 0 }, // 142
  { 29,    0, "P b c 21"  ,   0, "ba-c", "P 2c -2b"      , 7 }, // 143
  { 29,    0, "P 21 a b"  ,   0,  "cab", "P -2b 2a"      , 1 }, // 144
  { 29,    0, "P 21 c a"  ,   0, "-cba", "P -2ac 2a"     , 3 }, // 145
  { 29,    0, "P c 21 b"  ,   0,  "bca", "P -2bc -2c"    , 2 }, // 146
  { 29,    0, "P b 21 a"  ,   0, "a-cb", "P -2a -2ab"    , 5 }, // 147
  { 30,   30, "P n c 2"   ,   0,     "", "P 2 -2bc"      , 0 }, // 148
  { 30,    0, "P c n 2"   ,   0, "ba-c", "P 2 -2ac"      , 7 }, // 149
  { 30,    0, "P 2 n a"   ,   0,  "cab", "P -2ac 2"      , 1 }, // 150
  { 30,    0, "P 2 a n"   ,   0, "-cba", "P -2ab 2"      , 3 }, // 151
  { 30,    0, "P b 2 n"   ,   0,  "bca", "P -2ab -2ab"   , 2 }, // 152
  { 30,    0, "P n 2 b"   ,   0, "a-cb", "P -2bc -2bc"   , 5 }, // 153
  { 31,   31, "P m n 21"  ,   0,     "", "P 2ac -2"      , 0 }, // 154
  { 31,    0, "P n m 21"  ,   0, "ba-c", "P 2bc -2bc"    , 7 }, // 155
  { 31,    0, "P 21 m n"  ,   0,  "cab", "P -2ab 2ab"    , 1 }, // 156
  { 31,    0, "P 21 n m"  ,   0, "-cba", "P -2 2ac"      , 3 }, // 157
  { 31,    0, "P n 21 m"  ,   0,  "bca", "P -2 -2bc"     , 2 }, // 158
  { 31,    0, "P m 21 n"  ,   0, "a-cb", "P -2ab -2"     , 5 }, // 159
  { 32,   32, "P b a 2"   ,   0,     "", "P 2 -2ab"      , 0 }, // 160
  { 32,    0, "P 2 c b"   ,   0,  "cab", "P -2bc 2"      , 1 }, // 161
  { 32,    0, "P c 2 a"   ,   0,  "bca", "P -2ac -2ac"   , 2 }, // 162
  { 33,   33, "P n a 21"  ,   0,     "", "P 2c -2n"      , 0 }, // 163
  { 33,    0, "P b n 21"  ,   0, "ba-c", "P 2c -2ab"     , 7 }, // 164
  { 33,    0, "P 21 n b"  ,   0,  "cab", "P -2bc 2a"     , 1 }, // 165
  { 33,    0, "P 21 c n"  ,   0, "-cba", "P -2n 2a"      , 3 }, // 166
  { 33,    0, "P c 21 n"  ,   0,  "bca", "P -2n -2ac"    , 2 }, // 167
  { 33,    0, "P n 21 a"  ,   0, "a-cb", "P -2ac -2n"    , 5 }, // 168
  { 34,   34, "P n n 2"   ,   0,     "", "P 2 -2n"       , 0 }, // 169
  { 34,    0, "P 2 n n"   ,   0,  "cab", "P -2n 2"       , 1 }, // 170
  { 34,    0, "P n 2 n"   ,   0,  "bca", "P -2n -2n"     , 2 }, // 171
  { 35,   35, "C m m 2"   ,   0,     "", "C 2 -2"        , 0 }, // 172
  { 35,    0, "A 2 m m"   ,   0,  "cab", "A -2 2"        , 1 }, // 173
  { 35,    0, "B m 2 m"   ,   0,  "bca", "B -2 -2"       , 2 }, // 174
  { 36,   36, "C m c 21"  ,   0,     "", "C 2c -2"       , 0 }, // 175
  { 36,    0, "C c m 21"  ,   0, "ba-c", "C 2c -2c"      , 7 }, // 176
  { 36,    0, "A 21 m a"  ,   0,  "cab", "A -2a 2a"      , 1 }, // 177
  { 36,    0, "A 21 a m"  ,   0, "-cba", "A -2 2a"       , 3 }, // 178
  { 36,    0, "B b 21 m"  ,   0,  "bca", "B -2 -2b"      , 2 }, // 179
  { 36,    0, "B m 21 b"  ,   0, "a-cb", "B -2b -2"      , 5 }, // 180
  { 37,   37, "C c c 2"   ,   0,     "", "C 2 -2c"       , 0 }, // 181
  { 37,    0, "A 2 a a"   ,   0,  "cab", "A -2a 2"       , 1 }, // 182
  { 37,    0, "B b 2 b"   ,   0,  "bca", "B -2b -2b"     , 2 }, // 183
  { 38,   38, "A m m 2"   ,   0,     "", "A 2 -2"        , 0 }, // 184
  { 38,    0, "B m m 2"   ,   0, "ba-c", "B 2 -2"        , 7 }, // 185
  { 38,    0, "B 2 m m"   ,   0,  "cab", "B -2 2"        , 1 }, // 186
  { 38,    0, "C 2 m m"   ,   0, "-cba", "C -2 2"        , 3 }, // 187
  { 38,    0, "C m 2 m"   ,   0,  "bca", "C -2 -2"       , 2 }, // 188
  { 38,    0, "A m 2 m"   ,   0, "a-cb", "A -2 -2"       , 5 }, // 189
  { 39,   39, "A b m 2"   ,   0,     "", "A 2 -2b"       , 0 }, // 190
  { 39,    0, "B m a 2"   ,   0, "ba-c", "B 2 -2a"       , 7 }, // 191
  { 39,    0, "B 2 c m"   ,   0,  "cab", "B -2a 2"       , 1 }, // 192
  { 39,    0, "C 2 m b"   ,   0, "-cba", "C -2a 2"       , 3 }, // 193
  { 39,    0, "C m 2 a"   ,   0,  "bca", "C -2a -2a"     , 2 }, // 194
  { 39,    0, "A c 2 m"   ,   0, "a-cb", "A -2b -2b"     , 5 }, // 195
  { 40,   40, "A m a 2"   ,   0,     "", "A 2 -2a"       , 0 }, // 196
  { 40,    0, "B b m 2"   ,   0, "ba-c", "B 2 -2b"       , 7 }, // 197
  { 40,    0, "B 2 m b"   ,   0,  "cab", "B -2b 2"       , 1 }, // 198
  { 40,    0, "C 2 c m"   ,   0, "-cba", "C -2c 2"       , 3 }, // 199
  { 40,    0, "C c 2 m"   ,   0,  "bca", "C -2c -2c"     , 2 }, // 200
  { 40,    0, "A m 2 a"   ,   0, "a-cb", "A -2a -2a"     , 5 }, // 201
  { 41,   41, "A b a 2"   ,   0,     "", "A 2 -2ab"      , 0 }, // 202
  { 41,    0, "B b a 2"   ,   0, "ba-c", "B 2 -2ab"      , 7 }, // 203
  { 41,    0, "B 2 c b"   ,   0,  "cab", "B -2ab 2"      , 1 }, // 204
  { 41,    0, "C 2 c b"   ,   0, "-cba", "C -2ac 2"      , 3 }, // 205
  { 41,    0, "C c 2 a"   ,   0,  "bca", "C -2ac -2ac"   , 2 }, // 206
  { 41,    0, "A c 2 a"   ,   0, "a-cb", "A -2ab -2ab"   , 5 }, // 207
  { 42,   42, "F m m 2"   ,   0,     "", "F 2 -2"        , 0 }, // 208
  { 42,    0, "F 2 m m"   ,   0,  "cab", "F -2 2"        , 1 }, // 209
  { 42,    0, "F m 2 m"   ,   0,  "bca", "F -2 -2"       , 2 }, // 210
  { 43,   43, "F d d 2"   ,   0,     "", "F 2 -2d"       , 0 }, // 211
  { 43,    0, "F 2 d d"   ,   0,  "cab", "F -2d 2"       , 1 }, // 212
  { 43,    0, "F d 2 d"   ,   0,  "bca", "F -2d -2d"     , 2 }, // 213
  { 44,   44, "I m m 2"   ,   0,     "", "I 2 -2"        , 0 }, // 214
  { 44,    0, "I 2 m m"   ,   0,  "cab", "I -2 2"        , 1 }, // 215
  { 44,    0, "I m 2 m"   ,   0,  "bca", "I -2 -2"       , 2 }, // 216
  { 45,   45, "I b a 2"   ,   0,     "", "I 2 -2c"       , 0 }, // 217
  { 45,    0, "I 2 c b"   ,   0,  "cab", "I -2a 2"       , 1 }, // 218
  { 45,    0, "I c 2 a"   ,   0,  "bca", "I -2b -2b"     , 2 }, // 219
  { 46,   46, "I m a 2"   ,   0,     "", "I 2 -2a"       , 0 }, // 220
  { 46,    0, "I b m 2"   ,   0, "ba-c", "I 2 -2b"       , 7 }, // 221
  { 46,    0, "I 2 m b"   ,   0,  "cab", "I -2b 2"       , 1 }, // 222
  { 46,    0, "I 2 c m"   ,   0, "-cba", "I -2c 2"       , 3 }, // 223
  { 46,    0, "I c 2 m"   ,   0,  "bca", "I -2c -2c"     , 2 }, // 224
  { 46,    0, "I m 2 a"   ,   0, "a-cb", "I -2a -2a"     , 5 }, // 225
  { 47,   47, "P m m m"   ,   0,     "", "-P 2 2"        , 0 }, // 226
  { 48,   48, "P n n n"   , '1',     "", "P 2 2 -1n"     , 20}, // 227
  { 48,    0, "P n n n"   , '2',     "", "-P 2ab 2bc"    , 0 }, // 228
  { 49,   49, "P c c m"   ,   0,     "", "-P 2 2c"       , 0 }, // 229
  { 49,    0, "P m a a"   ,   0,  "cab", "-P 2a 2"       , 1 }, // 230
  { 49,    0, "P b m b"   ,   0,  "bca", "-P 2b 2b"      , 2 }, // 231
  { 50,   50, "P b a n"   , '1',     "", "P 2 2 -1ab"    , 21}, // 232
  { 50,    0, "P b a n"   , '2',     "", "-P 2ab 2b"     , 0 }, // 233
  { 50,    0, "P n c b"   , '1',  "cab", "P 2 2 -1bc"    , 22}, // 234
  { 50,    0, "P n c b"   , '2',  "cab", "-P 2b 2bc"     , 1 }, // 235
  { 50,    0, "P c n a"   , '1',  "bca", "P 2 2 -1ac"    , 23}, // 236
  { 50,    0, "P c n a"   , '2',  "bca", "-P 2a 2c"      , 2 }, // 237
  { 51,   51, "P m m a"   ,   0,     "", "-P 2a 2a"      , 0 }, // 238
  { 51,    0, "P m m b"   ,   0, "ba-c", "-P 2b 2"       , 7 }, // 239
  { 51,    0, "P b m m"   ,   0,  "cab", "-P 2 2b"       , 1 }, // 240
  { 51,    0, "P c m m"   ,   0, "-cba", "-P 2c 2c"      , 3 }, // 241
  { 51,    0, "P m c m"   ,   0,  "bca", "-P 2c 2"       , 2 }, // 242
  { 51,    0, "P m a m"   ,   0, "a-cb", "-P 2 2a"       , 5 }, // 243
  { 52,   52, "P n n a"   ,   0,     "", "-P 2a 2bc"     , 0 }, // 244
  { 52,    0, "P n n b"   ,   0, "ba-c", "-P 2b 2n"      , 7 }, // 245
  { 52,    0, "P b n n"   ,   0,  "cab", "-P 2n 2b"      , 1 }, // 246
  { 52,    0, "P c n n"   ,   0, "-cba", "-P 2ab 2c"     , 3 }, // 247
  { 52,    0, "P n c n"   ,   0,  "bca", "-P 2ab 2n"     , 2 }, // 248
  { 52,    0, "P n a n"   ,   0, "a-cb", "-P 2n 2bc"     , 5 }, // 249
  { 53,   53, "P m n a"   ,   0,     "", "-P 2ac 2"      , 0 }, // 250
  { 53,    0, "P n m b"   ,   0, "ba-c", "-P 2bc 2bc"    , 7 }, // 251
  { 53,    0, "P b m n"   ,   0,  "cab", "-P 2ab 2ab"    , 1 }, // 252
  { 53,    0, "P c n m"   ,   0, "-cba", "-P 2 2ac"      , 3 }, // 253
  { 53,    0, "P n c m"   ,   0,  "bca", "-P 2 2bc"      , 2 }, // 254
  { 53,    0, "P m a n"   ,   0, "a-cb", "-P 2ab 2"      , 5 }, // 255
  { 54,   54, "P c c a"   ,   0,     "", "-P 2a 2ac"     , 0 }, // 256
  { 54,    0, "P c c b"   ,   0, "ba-c", "-P 2b 2c"      , 7 }, // 257
  { 54,    0, "P b a a"   ,   0,  "cab", "-P 2a 2b"      , 1 }, // 258
  { 54,    0, "P c a a"   ,   0, "-cba", "-P 2ac 2c"     , 3 }, // 259
  { 54,    0, "P b c b"   ,   0,  "bca", "-P 2bc 2b"     , 2 }, // 260
  { 54,    0, "P b a b"   ,   0, "a-cb", "-P 2b 2ab"     , 5 }, // 261
  { 55,   55, "P b a m"   ,   0,     "", "-P 2 2ab"      , 0 }, // 262
  { 55,    0, "P m c b"   ,   0,  "cab", "-P 2bc 2"      , 1 }, // 263
  { 55,    0, "P c m a"   ,   0,  "bca", "-P 2ac 2ac"    , 2 }, // 264
  { 56,   56, "P c c n"   ,   0,     "", "-P 2ab 2ac"    , 0 }, // 265
  { 56,    0, "P n a a"   ,   0,  "cab", "-P 2ac 2bc"    , 1 }, // 266
  { 56,    0, "P b n b"   ,   0,  "bca", "-P 2bc 2ab"    , 2 }, // 267
  { 57,   57, "P b c m"   ,   0,     "", "-P 2c 2b"      , 0 }, // 268
  { 57,    0, "P c a m"   ,   0, "ba-c", "-P 2c 2ac"     , 7 }, // 269
  { 57,    0, "P m c a"   ,   0,  "cab", "-P 2ac 2a"     , 1 }, // 270
  { 57,    0, "P m a b"   ,   0, "-cba", "-P 2b 2a"      , 3 }, // 271
  { 57,    0, "P b m a"   ,   0,  "bca", "-P 2a 2ab"     , 2 }, // 272
  { 57,    0, "P c m b"   ,   0, "a-cb", "-P 2bc 2c"     , 5 }, // 273
  { 58,   58, "P n n m"   ,   0,     "", "-P 2 2n"       , 0 }, // 274
  { 58,    0, "P m n n"   ,   0,  "cab", "-P 2n 2"       , 1 }, // 275
  { 58,    0, "P n m n"   ,   0,  "bca", "-P 2n 2n"      , 2 }, // 276
  { 59,   59, "P m m n"   , '1',     "", "P 2 2ab -1ab"  , 21}, // 277
  { 59, 1059, "P m m n"   , '2',     "", "-P 2ab 2a"     , 0 }, // 278
  { 59,    0, "P n m m"   , '1',  "cab", "P 2bc 2 -1bc"  , 22}, // 279
  { 59,    0, "P n m m"   , '2',  "cab", "-P 2c 2bc"     , 1 }, // 280
  { 59,    0, "P m n m"   , '1',  "bca", "P 2ac 2ac -1ac", 23}, // 281
  { 59,    0, "P m n m"   , '2',  "bca", "-P 2c 2a"      , 2 }, // 282
  { 60,   60, "P b c n"   ,   0,     "", "-P 2n 2ab"     , 0 }, // 283
  { 60,    0, "P c a n"   ,   0, "ba-c", "-P 2n 2c"      , 7 }, // 284
  { 60,    0, "P n c a"   ,   0,  "cab", "-P 2a 2n"      , 1 }, // 285
  { 60,    0, "P n a b"   ,   0, "-cba", "-P 2bc 2n"     , 3 }, // 286
  { 60,    0, "P b n a"   ,   0,  "bca", "-P 2ac 2b"     , 2 }, // 287
  { 60,    0, "P c n b"   ,   0, "a-cb", "-P 2b 2ac"     , 5 }, // 288
  { 61,   61, "P b c a"   ,   0,     "", "-P 2ac 2ab"    , 0 }, // 289
  { 61,    0, "P c a b"   ,   0, "ba-c", "-P 2bc 2ac"    , 3 }, // 290
  { 62,   62, "P n m a"   ,   0,     "", "-P 2ac 2n"     , 0 }, // 291
  { 62,    0, "P m n b"   ,   0, "ba-c", "-P 2bc 2a"     , 7 }, // 292
  { 62,    0, "P b n m"   ,   0,  "cab", "-P 2c 2ab"     , 1 }, // 293
  { 62,    0, "P c m n"   ,   0, "-cba", "-P 2n 2ac"     , 3 }, // 294
  { 62,    0, "P m c n"   ,   0,  "bca", "-P 2n 2a"      , 2 }, // 295
  { 62,    0, "P n a m"   ,   0, "a-cb", "-P 2c 2n"      , 5 }, // 296
  { 63,   63, "C m c m"   ,   0,     "", "-C 2c 2"       , 0 }, // 297
  { 63,    0, "C c m m"   ,   0, "ba-c", "-C 2c 2c"      , 7 }, // 298
  { 63,    0, "A m m a"   ,   0,  "cab", "-A 2a 2a"      , 1 }, // 299
  { 63,    0, "A m a m"   ,   0, "-cba", "-A 2 2a"       , 3 }, // 300
  { 63,    0, "B b m m"   ,   0,  "bca", "-B 2 2b"       , 2 }, // 301
  { 63,    0, "B m m b"   ,   0, "a-cb", "-B 2b 2"       , 5 }, // 302
  { 64,   64, "C m c a"   ,   0,     "", "-C 2ac 2"      , 0 }, // 303
  { 64,    0, "C c m b"   ,   0, "ba-c", "-C 2ac 2ac"    , 7 }, // 304
  { 64,    0, "A b m a"   ,   0,  "cab", "-A 2ab 2ab"    , 1 }, // 305
  { 64,    0, "A c a m"   ,   0, "-cba", "-A 2 2ab"      , 3 }, // 306
  { 64,    0, "B b c m"   ,   0,  "bca", "-B 2 2ab"      , 2 }, // 307
  { 64,    0, "B m a b"   ,   0, "a-cb", "-B 2ab 2"      , 5 }, // 308
  { 65,   65, "C m m m"   ,   0,     "", "-C 2 2"        , 0 }, // 309
  { 65,    0, "A m m m"   ,   0,  "cab", "-A 2 2"        , 1 }, // 310
  { 65,    0, "B m m m"   ,   0,  "bca", "-B 2 2"        , 2 }, // 311
  { 66,   66, "C c c m"   ,   0,     "", "-C 2 2c"       , 0 }, // 312
  { 66,    0, "A m a a"   ,   0,  "cab", "-A 2a 2"       , 1 }, // 313
  { 66,    0, "B b m b"   ,   0,  "bca", "-B 2b 2b"      , 2 }, // 314
  { 67,   67, "C m m a"   ,   0,     "", "-C 2a 2"       , 0 }, // 315
  { 67,    0, "C m m b"   ,   0, "ba-c", "-C 2a 2a"      , 14}, // 316
  { 67,    0, "A b m m"   ,   0,  "cab", "-A 2b 2b"      , 1 }, // 317
  { 67,    0, "A c m m"   ,   0, "-cba", "-A 2 2b"       , 3 }, // 318
  { 67,    0, "B m c m"   ,   0,  "bca", "-B 2 2a"       , 2 }, // 319
  { 67,    0, "B m a m"   ,   0, "a-cb", "-B 2a 2"       , 5 }, // 320
  { 68,   68, "C c c a"   , '1',     "", "C 2 2 -1ac"    , 24}, // 321
  { 68,    0, "C c c a"   , '2',     "", "-C 2a 2ac"     , 0 }, // 322
  { 68,    0, "C c c b"   , '1', "ba-c", "C 2 2 -1ac"    , 24}, // 323
  { 68,    0, "C c c b"   , '2', "ba-c", "-C 2a 2c"      , 21}, // 324
  { 68,    0, "A b a a"   , '1',  "cab", "A 2 2 -1ab"    , 25}, // 325
  { 68,    0, "A b a a"   , '2',  "cab", "-A 2a 2b"      , 1 }, // 326
  { 68,    0, "A c a a"   , '1', "-cba", "A 2 2 -1ab"    , 25}, // 327
  { 68,    0, "A c a a"   , '2', "-cba", "-A 2ab 2b"     , 3 }, // 328
  { 68,    0, "B b c b"   , '1',  "bca", "B 2 2 -1ab"    , 26}, // 329
  { 68,    0, "B b c b"   , '2',  "bca", "-B 2ab 2b"     , 2 }, // 330
  { 68,    0, "B b a b"   , '1', "a-cb", "B 2 2 -1ab"    , 26}, // 331
  { 68,    0, "B b a b"   , '2', "a-cb", "-B 2b 2ab"     , 5 }, // 332
  { 69,   69, "F m m m"   ,   0,     "", "-F 2 2"        , 0 }, // 333
  { 70,   70, "F d d d"   , '1',     "", "F 2 2 -1d"     , 27}, // 334
  { 70,    0, "F d d d"   , '2',     "", "-F 2uv 2vw"    , 0 }, // 335
  { 71,   71, "I m m m"   ,   0,     "", "-I 2 2"        , 0 }, // 336
  { 72,   72, "I b a m"   ,   0,     "", "-I 2 2c"       , 0 }, // 337
  { 72,    0, "I m c b"   ,   0,  "cab", "-I 2a 2"       , 1 }, // 338
  { 72,    0, "I c m a"   ,   0,  "bca", "-I 2b 2b"      , 2 }, // 339
  { 73,   73, "I b c a"   ,   0,     "", "-I 2b 2c"      , 0 }, // 340
  { 73,    0, "I c a b"   ,   0, "ba-c", "-I 2a 2b"      , 28}, // 341
  { 74,   74, "I m m a"   ,   0,     "", "-I 2b 2"       , 0 }, // 342
  { 74,    0, "I m m b"   ,   0, "ba-c", "-I 2a 2a"      , 28}, // 343
  { 74,    0, "I b m m"   ,   0,  "cab", "-I 2c 2c"      , 1 }, // 344
  { 74,    0, "I c m m"   ,   0, "-cba", "-I 2 2b"       , 3 }, // 345
  { 74,    0, "I m c m"   ,   0,  "bca", "-I 2 2a"       , 2 }, // 346
  { 74,    0, "I m a m"   ,   0, "a-cb", "-I 2c 2"       , 5 }, // 347
  { 75,   75, "P 4"       ,   0,     "", "P 4"           , 0 }, // 348
  { 76,   76, "P 41"      ,   0,     "", "P 4w"          , 0 }, // 349
  { 77,   77, "P 42"      ,   0,     "", "P 4c"          , 0 }, // 350
  { 78,   78, "P 43"      ,   0,     "", "P 4cw"         , 0 }, // 351
  { 79,   79, "I 4"       ,   0,     "", "I 4"           , 0 }, // 352
  { 80,   80, "I 41"      ,   0,     "", "I 4bw"         , 0 }, // 353
  { 81,   81, "P -4"      ,   0,     "", "P -4"          , 0 }, // 354
  { 82,   82, "I -4"      ,   0,     "", "I -4"          , 0 }, // 355
  { 83,   83, "P 4/m"     ,   0,     "", "-P 4"          , 0 }, // 356
  { 84,   84, "P 42/m"    ,   0,     "", "-P 4c"         , 0 }, // 357
  { 85,   85, "P 4/n"     , '1',     "", "P 4ab -1ab"    , 29}, // 358
  { 85,    0, "P 4/n"     , '2',     "", "-P 4a"         , 0 }, // 359
  { 86,   86, "P 42/n"    , '1',     "", "P 4n -1n"      , 30}, // 360
  { 86,    0, "P 42/n"    , '2',     "", "-P 4bc"        , 0 }, // 361
  { 87,   87, "I 4/m"     ,   0,     "", "-I 4"          , 0 }, // 362
  { 88,   88, "I 41/a"    , '1',     "", "I 4bw -1bw"    , 31}, // 363
  { 88,    0, "I 41/a"    , '2',     "", "-I 4ad"        , 0 }, // 364
  { 89,   89, "P 4 2 2"   ,   0,     "", "P 4 2"         , 0 }, // 365
  { 90,   90, "P 4 21 2"  ,   0,     "", "P 4ab 2ab"     , 0 }, // 366
  { 91,   91, "P 41 2 2"  ,   0,     "", "P 4w 2c"       , 0 }, // 367
  { 92,   92, "P 41 21 2" ,   0,     "", "P 4abw 2nw"    , 0 }, // 368
  { 93,   93, "P 42 2 2"  ,   0,     "", "P 4c 2"        , 0 }, // 369
  { 94,   94, "P 42 21 2" ,   0,     "", "P 4n 2n"       , 0 }, // 370
  { 95,   95, "P 43 2 2"  ,   0,     "", "P 4cw 2c"      , 0 }, // 371
  { 96,   96, "P 43 21 2" ,   0,     "", "P 4nw 2abw"    , 0 }, // 372
  { 97,   97, "I 4 2 2"   ,   0,     "", "I 4 2"         , 0 }, // 373
  { 98,   98, "I 41 2 2"  ,   0,     "", "I 4bw 2bw"     , 0 }, // 374
  { 99,   99, "P 4 m m"   ,   0,     "", "P 4 -2"        , 0 }, // 375
  {100,  100, "P 4 b m"   ,   0,     "", "P 4 -2ab"      , 0 }, // 376
  {101,  101, "P 42 c m"  ,   0,     "", "P 4c -2c"      , 0 }, // 377
  {102,  102, "P 42 n m"  ,   0,     "", "P 4n -2n"      , 0 }, // 378
  {103,  103, "P 4 c c"   ,   0,     "", "P 4 -2c"       , 0 }, // 379
  {104,  104, "P 4 n c"   ,   0,     "", "P 4 -2n"       , 0 }, // 380
  {105,  105, "P 42 m c"  ,   0,     "", "P 4c -2"       , 0 }, // 381
  {106,  106, "P 42 b c"  ,   0,     "", "P 4c -2ab"     , 0 }, // 382
  {107,  107, "I 4 m m"   ,   0,     "", "I 4 -2"        , 0 }, // 383
  {108,  108, "I 4 c m"   ,   0,     "", "I 4 -2c"       , 0 }, // 384
  {109,  109, "I 41 m d"  ,   0,     "", "I 4bw -2"      , 0 }, // 385
  {110,  110, "I 41 c d"  ,   0,     "", "I 4bw -2c"     , 0 }, // 386
  {111,  111, "P -4 2 m"  ,   0,     "", "P -4 2"        , 0 }, // 387
  {112,  112, "P -4 2 c"  ,   0,     "", "P -4 2c"       , 0 }, // 388
  {113,  113, "P -4 21 m" ,   0,     "", "P -4 2ab"      , 0 }, // 389
  {114,  114, "P -4 21 c" ,   0,     "", "P -4 2n"       , 0 }, // 390
  {115,  115, "P -4 m 2"  ,   0,     "", "P -4 -2"       , 0 }, // 391
  {116,  116, "P -4 c 2"  ,   0,     "", "P -4 -2c"      , 0 }, // 392
  {117,  117, "P -4 b 2"  ,   0,     "", "P -4 -2ab"     , 0 }, // 393
  {118,  118, "P -4 n 2"  ,   0,     "", "P -4 -2n"      , 0 }, // 394
  {119,  119, "I -4 m 2"  ,   0,     "", "I -4 -2"       , 0 }, // 395
  {120,  120, "I -4 c 2"  ,   0,     "", "I -4 -2c"      , 0 }, // 396
  {121,  121, "I -4 2 m"  ,   0,     "", "I -4 2"        , 0 }, // 397
  {122,  122, "I -4 2 d"  ,   0,     "", "I -4 2bw"      , 0 }, // 398
  {123,  123, "P 4/m m m" ,   0,     "", "-P 4 2"        , 0 }, // 399
  {124,  124, "P 4/m c c" ,   0,     "", "-P 4 2c"       , 0 }, // 400
  {125,  125, "P 4/n b m" , '1',     "", "P 4 2 -1ab"    , 21}, // 401
  {125,    0, "P 4/n b m" , '2',     "", "-P 4a 2b"      , 0 }, // 402
  {126,  126, "P 4/n n c" , '1',     "", "P 4 2 -1n"     , 20}, // 403
  {126,    0, "P 4/n n c" , '2',     "", "-P 4a 2bc"     , 0 }, // 404
  {127,  127, "P 4/m b m" ,   0,     "", "-P 4 2ab"      , 0 }, // 405
  {128,  128, "P 4/m n c" ,   0,     "", "-P 4 2n"       , 0 }, // 406
  {129,  129, "P 4/n m m" , '1',     "", "P 4ab 2ab -1ab", 29}, // 407
  {129,    0, "P 4/n m m" , '2',     "", "-P 4a 2a"      , 0 }, // 408
  {130,  130, "P 4/n c c" , '1',     "", "P 4ab 2n -1ab" , 29}, // 409
  {130,    0, "P 4/n c c" , '2',     "", "-P 4a 2ac"     , 0 }, // 410
  {131,  131, "P 42/m m c",   0,     "", "-P 4c 2"       , 0 }, // 411
  {132,  132, "P 42/m c m",   0,     "", "-P 4c 2c"      , 0 }, // 412
  {133,  133, "P 42/n b c", '1',     "", "P 4n 2c -1n"   , 32}, // 413
  {133,    0, "P 42/n b c", '2',     "", "-P 4ac 2b"     , 0 }, // 414
  {134,  134, "P 42/n n m", '1',     "", "P 4n 2 -1n"    , 33}, // 415
  {134,    0, "P 42/n n m", '2',     "", "-P 4ac 2bc"    , 0 }, // 416
  {135,  135, "P 42/m b c",   0,     "", "-P 4c 2ab"     , 0 }, // 417
  {136,  136, "P 42/m n m",   0,     "", "-P 4n 2n"      , 0 }, // 418
  {137,  137, "P 42/n m c", '1',     "", "P 4n 2n -1n"   , 32}, // 419
  {137,    0, "P 42/n m c", '2',     "", "-P 4ac 2a"     , 0 }, // 420
  {138,  138, "P 42/n c m", '1',     "", "P 4n 2ab -1n"  , 33}, // 421
  {138,    0, "P 42/n c m", '2',     "", "-P 4ac 2ac"    , 0 }, // 422
  {139,  139, "I 4/m m m" ,   0,     "", "-I 4 2"        , 0 }, // 423
  {140,  140, "I 4/m c m" ,   0,     "", "-I 4 2c"       , 0 }, // 424
  {141,  141, "I 41/a m d", '1',     "", "I 4bw 2bw -1bw", 34}, // 425
  {141,    0, "I 41/a m d", '2',     "", "-I 4bd 2"      , 0 }, // 426
  {142,  142, "I 41/a c d", '1',     "", "I 4bw 2aw -1bw", 35}, // 427
  {142,    0, "I 41/a c d", '2',     "", "-I 4bd 2c"     , 0 }, // 428
  {143,  143, "P 3"       ,   0,     "", "P 3"           , 0 }, // 429
  {144,  144, "P 31"      ,   0,     "", "P 31"          , 0 }, // 430
  {145,  145, "P 32"      ,   0,     "", "P 32"          , 0 }, // 431
  {146,  146, "R 3"       , 'H',     "", "R 3"           , 0 }, // 432
  {146, 1146, "R 3"       , 'R',     "", "P 3*"          , 36}, // 433
  {147,  147, "P -3"      ,   0,     "", "-P 3"          , 0 }, // 434
  {148,  148, "R -3"      , 'H',     "", "-R 3"          , 0 }, // 435
  {148, 1148, "R -3"      , 'R',     "", "-P 3*"         , 36}, // 436
  {149,  149, "P 3 1 2"   ,   0,     "", "P 3 2"         , 0 }, // 437
  {150,  150, "P 3 2 1"   ,   0,     "", "P 3 2\""       , 0 }, // 438
  {151,  151, "P 31 1 2"  ,   0,     "", "P 31 2 (0 0 4)", 0 }, // 439
  {152,  152, "P 31 2 1"  ,   0,     "", "P 31 2\""      , 0 }, // 440
  {153,  153, "P 32 1 2"  ,   0,     "", "P 32 2 (0 0 2)", 0 }, // 441
  {154,  154, "P 32 2 1"  ,   0,     "", "P 32 2\""      , 0 }, // 442
  {155,  155, "R 3 2"     , 'H',     "", "R 3 2\""       , 0 }, // 443
  {155, 1155, "R 3 2"     , 'R',     "", "P 3* 2"        , 36}, // 444
  {156,  156, "P 3 m 1"   ,   0,     "", "P 3 -2\""      , 0 }, // 445
  {157,  157, "P 3 1 m"   ,   0,     "", "P 3 -2"        , 0 }, // 446
  {158,  158, "P 3 c 1"   ,   0,     "", "P 3 -2\"c"     , 0 }, // 447
  {159,  159, "P 3 1 c"   ,   0,     "", "P 3 -2c"       , 0 }, // 448
  {160,  160, "R 3 m"     , 'H',     "", "R 3 -2\""      , 0 }, // 449
  {160, 1160, "R 3 m"     , 'R',     "", "P 3* -2"       , 36}, // 450
  {161,  161, "R 3 c"     , 'H',     "", "R 3 -2\"c"     , 0 }, // 451
  {161, 1161, "R 3 c"     , 'R',     "", "P 3* -2n"      , 36}, // 452
  {162,  162, "P -3 1 m"  ,   0,     "", "-P 3 2"        , 0 }, // 453
  {163,  163, "P -3 1 c"  ,   0,     "", "-P 3 2c"       , 0 }, // 454
  {164,  164, "P -3 m 1"  ,   0,     "", "-P 3 2\""      , 0 }, // 455
  {165,  165, "P -3 c 1"  ,   0,     "", "-P 3 2\"c"     , 0 }, // 456
  {166,  166, "R -3 m"    , 'H',     "", "-R 3 2\""      , 0 }, // 457
  {166, 1166, "R -3 m"    , 'R',     "", "-P 3* 2"       , 36}, // 458
  {167,  167, "R -3 c"    , 'H',     "", "-R 3 2\"c"     , 0 }, // 459
  {167, 1167, "R -3 c"    , 'R',     "", "-P 3* 2n"      , 36}, // 460
  {168,  168, "P 6"       ,   0,     "", "P 6"           , 0 }, // 461
  {169,  169, "P 61"      ,   0,     "", "P 61"          , 0 }, // 462
  {170,  170, "P 65"      ,   0,     "", "P 65"          , 0 }, // 463
  {171,  171, "P 62"      ,   0,     "", "P 62"          , 0 }, // 464
  {172,  172, "P 64"      ,   0,     "", "P 64"          , 0 }, // 465
  {173,  173, "P 63"      ,   0,     "", "P 6c"          , 0 }, // 466
  {174,  174, "P -6"      ,   0,     "", "P -6"          , 0 }, // 467
  {175,  175, "P 6/m"     ,   0,     "", "-P 6"          , 0 }, // 468
  {176,  176, "P 63/m"    ,   0,     "", "-P 6c"         , 0 }, // 469
  {177,  177, "P 6 2 2"   ,   0,     "", "P 6 2"         , 0 }, // 470
  {178,  178, "P 61 2 2"  ,   0,     "", "P 61 2 (0 0 5)", 0 }, // 471
  {179,  179, "P 65 2 2"  ,   0,     "", "P 65 2 (0 0 1)", 0 }, // 472
  {180,  180, "P 62 2 2"  ,   0,     "", "P 62 2 (0 0 4)", 0 }, // 473
  {181,  181, "P 64 2 2"  ,   0,     "", "P 64 2 (0 0 2)", 0 }, // 474
  {182,  182, "P 63 2 2"  ,   0,     "", "P 6c 2c"       , 0 }, // 475
  {183,  183, "P 6 m m"   ,   0,     "", "P 6 -2"        , 0 }, // 476
  {184,  184, "P 6 c c"   ,   0,     "", "P 6 -2c"       , 0 }, // 477
  {185,  185, "P 63 c m"  ,   0,     "", "P 6c -2"       , 0 }, // 478
  {186,  186, "P 63 m c"  ,   0,     "", "P 6c -2c"      , 0 }, // 479
  {187,  187, "P -6 m 2"  ,   0,     "", "P -6 2"        , 0 }, // 480
  {188,  188, "P -6 c 2"  ,   0,     "", "P -6c 2"       , 0 }, // 481
  {189,  189, "P -6 2 m"  ,   0,     "", "P -6 -2"       , 0 }, // 482
  {190,  190, "P -6 2 c"  ,   0,     "", "P -6c -2c"     , 0 }, // 483
  {191,  191, "P 6/m m m" ,   0,     "", "-P 6 2"        , 0 }, // 484
  {192,  192, "P 6/m c c" ,   0,     "", "-P 6 2c"       , 0 }, // 485
  {193,  193, "P 63/m c m",   0,     "", "-P 6c 2"       , 0 }, // 486
  {194,  194, "P 63/m m c",   0,     "", "-P 6c 2c"      , 0 }, // 487
  {195,  195, "P 2 3"     ,   0,     "", "P 2 2 3"       , 0 }, // 488
  {196,  196, "F 2 3"     ,   0,     "", "F 2 2 3"       , 0 }, // 489
  {197,  197, "I 2 3"     ,   0,     "", "I 2 2 3"       , 0 }, // 490
  {198,  198, "P 21 3"    ,   0,     "", "P 2ac 2ab 3"   , 0 }, // 491
  {199,  199, "I 21 3"    ,   0,     "", "I 2b 2c 3"     , 0 }, // 492
  {200,  200, "P m -3"    ,   0,     "", "-P 2 2 3"      , 0 }, // 493
  {201,  201, "P n -3"    , '1',     "", "P 2 2 3 -1n"   , 20}, // 494
  {201,    0, "P n -3"    , '2',     "", "-P 2ab 2bc 3"  , 0 }, // 495
  {202,  202, "F m -3"    ,   0,     "", "-F 2 2 3"      , 0 }, // 496
  {203,  203, "F d -3"    , '1',     "", "F 2 2 3 -1d"   , 27}, // 497
  {203,    0, "F d -3"    , '2',     "", "-F 2uv 2vw 3"  , 0 }, // 498
  {204,  204, "I m -3"    ,   0,     "", "-I 2 2 3"      , 0 }, // 499
  {205,  205, "P a -3"    ,   0,     "", "-P 2ac 2ab 3"  , 0 }, // 500
  {206,  206, "I a -3"    ,   0,     "", "-I 2b 2c 3"    , 0 }, // 501
  {207,  207, "P 4 3 2"   ,   0,     "", "P 4 2 3"       , 0 }, // 502
  {208,  208, "P 42 3 2"  ,   0,     "", "P 4n 2 3"      , 0 }, // 503
  {209,  209, "F 4 3 2"   ,   0,     "", "F 4 2 3"       , 0 }, // 504
  {210,  210, "F 41 3 2"  ,   0,     "", "F 4d 2 3"      , 0 }, // 505
  {211,  211, "I 4 3 2"   ,   0,     "", "I 4 2 3"       , 0 }, // 506
  {212,  212, "P 43 3 2"  ,   0,     "", "P 4acd 2ab 3"  , 0 }, // 507
  {213,  213, "P 41 3 2"  ,   0,     "", "P 4bd 2ab 3"   , 0 }, // 508
  {214,  214, "I 41 3 2"  ,   0,     "", "I 4bd 2c 3"    , 0 }, // 509
  {215,  215, "P -4 3 m"  ,   0,     "", "P -4 2 3"      , 0 }, // 510
  {216,  216, "F -4 3 m"  ,   0,     "", "F -4 2 3"      , 0 }, // 511
  {217,  217, "I -4 3 m"  ,   0,     "", "I -4 2 3"      , 0 }, // 512
  {218,  218, "P -4 3 n"  ,   0,     "", "P -4n 2 3"     , 0 }, // 513
  {219,  219, "F -4 3 c"  ,   0,     "", "F -4a 2 3"     , 0 }, // 514
  {220,  220, "I -4 3 d"  ,   0,     "", "I -4bd 2c 3"   , 0 }, // 515
  {221,  221, "P m -3 m"  ,   0,     "", "-P 4 2 3"      , 0 }, // 516
  {222,  222, "P n -3 n"  , '1',     "", "P 4 2 3 -1n"   , 20}, // 517
  {222,    0, "P n -3 n"  , '2',     "", "-P 4a 2bc 3"   , 0 }, // 518
  {223,  223, "P m -3 n"  ,   0,     "", "-P 4n 2 3"     , 0 }, // 519
  {224,  224, "P n -3 m"  , '1',     "", "P 4n 2 3 -1n"  , 30}, // 520
  {224,    0, "P n -3 m"  , '2',     "", "-P 4bc 2bc 3"  , 0 }, // 521
  {225,  225, "F m -3 m"  ,   0,     "", "-F 4 2 3"      , 0 }, // 522
  {226,  226, "F m -3 c"  ,   0,     "", "-F 4a 2 3"     , 0 }, // 523
  {227,  227, "F d -3 m"  , '1',     "", "F 4d 2 3 -1d"  , 27}, // 524
  {227,    0, "F d -3 m"  , '2',     "", "-F 4vw 2vw 3"  , 0 }, // 525
  {228,  228, "F d -3 c"  , '1',     "", "F 4d 2 3 -1ad" , 37}, // 526
  {228,    0, "F d -3 c"  , '2',     "", "-F 4ud 2vw 3"  , 0 }, // 527
  {229,  229, "I m -3 m"  ,   0,     "", "-I 4 2 3"      , 0 }, // 528
  {230,  230, "I a -3 d"  ,   0,     "", "-I 4bd 2c 3"   , 0 }, // 529
  // And extra entries from syminfo.lib
  {  5, 5005, "I 1 21 1"  ,   0,     "", "I 2yb"         , 38}, // 530
  {  5, 3005, "C 1 21 1"  ,   0,     "", "C 2yb"         , 14}, // 531
  { 18, 1018, "P 21212(a)",   0,     "", "P 2ab 2a"      , 14}, // 532
  { 20, 1020, "C 2 2 21a)",   0,     "", "C 2ac 2"       , 39}, // 533
  { 21, 1021, "C 2 2 2a"  ,   0,     "", "C 2ab 2b"      , 14}, // 534
  { 22, 1022, "F 2 2 2a"  ,   0,     "", "F 2 2c"        , 40}, // 535
  { 23, 1023, "I 2 2 2a"  ,   0,     "", "I 2ab 2bc"     , 33}, // 536
  { 94, 1094, "P 42 21 2a",   0,     "", "P 4bc 2a"      , 20}, // 537
  {197, 1197, "I 2 3a"    ,   0,     "", "I 2ab 2bc 3"   , 30}, // 538
  // And extra entries from Open Babel, double checked with Crystallographic
  // Space Group Diagrams and Tables at http://img.chem.ucl.ac.uk/sgp/
  // triclinic - enlarged unit cells
  {  1,    0, "A 1"       ,   0,     "", "A 1"           , 41}, // 539
  {  1,    0, "B 1"       ,   0,     "", "B 1"           , 42}, // 540
  {  1,    0, "C 1"       ,   0,     "", "C 1"           , 43}, // 541
  {  1,    0, "F 1"       ,   0,     "", "F 1"           , 44}, // 542
  {  1,    0, "I 1"       ,   0,     "", "I 1"           , 45}, // 543
  {  2,    0, "A -1"      ,   0,     "", "-A 1"          , 41}, // 544
  {  2,    0, "B -1"      ,   0,     "", "-B 1"          , 42}, // 545
  {  2,    0, "C -1"      ,   0,     "", "-C 1"          , 43}, // 546
  {  2,    0, "F -1"      ,   0,     "", "-F 1"          , 44}, // 547
  {  2,    0, "I -1"      ,   0,     "", "-I 1"          , 45}, // 548
  // monoclinic
  {  4,    0, "C 1 1 21"  ,   0,     "", "C 2c"          , 46}, // 549
  { 12,    0, "F 1 2/m 1" ,   0,     "", "-F 2y"         , 47}, // 550
  { 64,    0, "A b a m"   ,   0,     "", "-A 2 2ab"      , 3 }, // 551
  // tetragonal - enlarged C- and F-centred unit cells
  {117,    0, "C -4 2 b"  ,   0,     "", "C -4 2ya"      , 48}, // 552
  { 97,    0, "F 4 2 2" ,     0,     "", "F 4 2"         , 48}, // 553
  {139,    0, "F 4/m m m" ,   0,     "", "-F 4 2"        , 48}, // 554
  { 89,    0, "C 4 2 2" ,     0,     "", "C 4 2"         , 48}, // 555
  { 90,    0, "C 4 2 21" ,    0,     "", "C 4a 2"        , 48}, // 556
};

template<class Dummy>
const SpaceGroupAltName Tables_<Dummy>::alt_names[28] = {
  // In 1990's ITfC vol.A changed some of the standard names, introducing
  // symbols 'e' and 'g'. sgtbx interprets these new symbols with
  // option ad_hoc_1992. spglib uses only the new symbols.
  {"A e m 2",   0, 190}, // A b m 2
  {"B m e 2",   0, 191}, // B m a 2
  {"B 2 e m",   0, 192}, // B 2 c m
  {"C 2 m e",   0, 193}, // C 2 m b
  {"C m 2 e",   0, 194}, // C m 2 a
  {"A e 2 m",   0, 195}, // A c 2 m
  {"A e a 2",   0, 202}, // A b a 2
  {"B b e 2",   0, 203}, // B b a 2
  {"B 2 e b",   0, 204}, // B 2 c b
  {"C 2 c e",   0, 205}, // C 2 c b
  {"C c 2 e",   0, 206}, // C c 2 a
  {"A e 2 a",   0, 207}, // A c 2 a
  {"C m c e",   0, 303}, // C m c a
  {"C c m e",   0, 304}, // C c m b
  {"A e m a",   0, 305}, // A b m a
  {"A e a m",   0, 306}, // A c a m
  {"B b e m",   0, 307}, // B b c m
  {"B m e b",   0, 308}, // B m a b
  {"C m m e",   0, 315}, // C m m a
  {"A e m m",   0, 317}, // A b m m
  {"B m e m",   0, 319}, // B m c m
  {"C c c e", '1', 321}, // C c c a
  {"C c c e", '2', 322}, // C c c a
  {"A e a a", '1', 325}, // A b a a
  {"A e a a", '2', 326}, // A b a a
  {"B b e b", '1', 329}, // B b c b
  {"B b e b", '2', 330}, // B b c b
  // help with  parsing of unusual setting names that are present in the PDB
  {"P 21 21 2a", 0, 532}, // P 21212(a)
};


// This table was generated by tools/gen_reciprocal_asu.py.
template<class Dummy>
const unsigned char Tables_<Dummy>::ccp4_hkl_asu[230] = {
  0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 7, 6, 7, 6, 7, 7, 7,
  6, 7, 6, 7, 7, 6, 6, 7, 7, 7, 7, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9
};

} // namespace impl

using spacegroup_tables = impl::Tables_<void>;

inline const SpaceGroup* find_spacegroup_by_number(int ccp4) noexcept {
  if (ccp4 == 0)
    return &spacegroup_tables::main[0];
  for (const SpaceGroup& sg : spacegroup_tables::main)
    if (sg.ccp4 == ccp4)
      return &sg;
  return nullptr;
}

inline const SpaceGroup& get_spacegroup_by_number(int ccp4) {
  const SpaceGroup* sg = find_spacegroup_by_number(ccp4);
  if (sg == nullptr)
    throw std::invalid_argument("Invalid space-group number: "
                                + std::to_string(ccp4));
  return *sg;
}

inline const SpaceGroup& get_spacegroup_reference_setting(int number) {
  for (const SpaceGroup& sg : spacegroup_tables::main)
    if (sg.number == number && sg.is_reference_setting())
      return sg;
  throw std::invalid_argument("Invalid space-group number: "
                              + std::to_string(number));
}

// the angles alpha and gamma are optional. If provided they are only used
// to distinguish hexagonal and rhombohedral settings (e.g. for "R 3").
inline const SpaceGroup* find_spacegroup_by_name(std::string name,
                                  double alpha=0., double gamma=0.) noexcept {
  const char* p = impl::skip_blank(name.c_str());
  if (*p >= '0' && *p <= '9') { // handle numbers
    char *endptr;
    long n = std::strtol(p, &endptr, 10);
    return *endptr == '\0' ? find_spacegroup_by_number(n) : nullptr;
  }
  char first = *p & ~0x20; // to uppercase
  if (first == '\0')
    return nullptr;
  if (first == 'H')
    first = 'R';
  p = impl::skip_blank(p+1);
  size_t start = p - name.c_str();
  // change letters to lower case, except the letter after :
  for (size_t i = start; i < name.size(); ++i) {
    if (name[i] >= 'A' && name[i] <= 'Z')
      name[i] |= 0x20;  // to lowercase
    else if (name[i] == ':')
      while (++i < name.size())
        if (name[i] >= 'a' && name[i] <= 'z')
          name[i] &= ~0x20;  // to uppercase
  }
  // The string that const char* p points to was just modified.
  // This confuses some compilers (GCC 4.8), so let's re-assign p.
  p = name.c_str() + start;

  for (const SpaceGroup& sg : spacegroup_tables::main)
    if (sg.hm[0] == first) {
      if (sg.hm[2] == *p) {
        const char* a = impl::skip_blank(p + 1);
        const char* b = impl::skip_blank(sg.hm + 3);
        while (*a == *b && *b != '\0') {
          a = impl::skip_blank(a+1);
          b = impl::skip_blank(b+1);
        }
        if (*b == '\0' &&
            (*a == '\0' || (*a == ':' && *impl::skip_blank(a+1) == sg.ext))) {
          // Change hexagonal settings to rhombohedral if the unit cell angles
          // are more consistent with the latter.
          // We have possible ambiguity in the hexagonal crystal family.
          // For instance, "R 3" may mean "R 3:H" (hexagonal setting) or
          // "R 3:R" (rhombohedral setting). The :H symbols come first
          // in the table and are used by default. The ratio gamma:alpha
          // is 120:90 in the hexagonal system and 1:1 in rhombohedral.
          // We assume that the 'R' entry follows directly the 'H' entry.
          if (*a == '\0' && sg.ext == 'H' && gamma < 1.125 * alpha)
            return &sg + 1;
          return &sg;
        }
      } else if (sg.hm[2] == '1' && sg.hm[3] == ' ') {
        // check monoclinic short names, matching P2 to "P 1 2 1";
        // as an exception "B 2" == "B 1 1 2" (like in the PDB)
        const char* b = sg.hm + 4;
        if (*b != '1' || (first == 'B' && *++b == ' ' && *++b != '1')) {
          char end = (b == sg.hm + 4 ? ' ' : '\0');
          const char* a = impl::skip_blank(p);
          while (*a == *b && *b != end) {
            ++a;
            ++b;
          }
          if (*impl::skip_blank(a) == '\0' && *b == end)
            return &sg;
        }
      }
    }
  for (const SpaceGroupAltName& sg : spacegroup_tables::alt_names)
    if (sg.hm[0] == first && sg.hm[2] == *p) {
      const char* a = impl::skip_blank(p + 1);
      const char* b = impl::skip_blank(sg.hm + 3);
      while (*a == *b && *b != '\0') {
        a = impl::skip_blank(a+1);
        b = impl::skip_blank(b+1);
      }
      if (*b == '\0' &&
          (*a == '\0' || (*a == ':' && *impl::skip_blank(a+1) == sg.ext)))
        return &spacegroup_tables::main[sg.pos];
    }
  return nullptr;
}

inline const SpaceGroup& get_spacegroup_by_name(const std::string& name) {
  const SpaceGroup* sg = find_spacegroup_by_name(name);
  if (sg == nullptr)
    throw std::invalid_argument("Unknown space-group name: " + name);
  return *sg;
}

inline const SpaceGroup& get_spacegroup_p1() {
  return spacegroup_tables::main[0];
}

inline const SpaceGroup* find_spacegroup_by_ops(const GroupOps& gops) {
  char c = gops.find_centering();
  for (const SpaceGroup& sg : spacegroup_tables::main)
    if ((c == sg.hall[0] || c == sg.hall[1]) &&
        gops.is_same_as(sg.operations()))
      return &sg;
  return nullptr;
}

// Reciprocal space asu (asymmetric unit).
// The same 12 choices of ASU as in CCP4 symlib and cctbx.
struct ReciprocalAsu {
  int idx;
  Op::Rot rot;
  bool is_ref;

  ReciprocalAsu(const SpaceGroup* sg) {
    if (sg == nullptr)
      fail("Missing space group");
    idx = spacegroup_tables::ccp4_hkl_asu[sg->number - 1];
    is_ref = sg->is_reference_setting();
    if (!is_ref)
      rot = sg->basisop().rot;
  }

  bool is_in(const Op::Miller& hkl) const {
    if (is_ref)
      return is_in_reference_setting(hkl[0], hkl[1], hkl[2]);
    Op::Miller r;
    for (int i = 0; i != 3; ++i)
      r[i] = rot[0][i] * hkl[0] + rot[1][i] * hkl[1] + rot[2][i] * hkl[2];
    return is_in_reference_setting(r[0], r[1], r[2]);
  }

  bool is_in_reference_setting(int h, int k, int l) const {
    switch (idx) {
      case 0: return l>0 || (l==0 && (h>0 || (h==0 && k>=0)));
      case 1: return k>=0 && (l>0 || (l==0 && h>=0));
      case 2: return h>=0 && k>=0 && l>=0;
      case 3: return l>=0 && ((h>=0 && k>0) || (h==0 && k==0));
      case 4: return h>=k && k>=0 && l>=0;
      case 5: return (h>=0 && k>0) || (h==0 && k==0 && l>=0);
      case 6: return h>=k && k>=0 && (k>0 || l>=0);
      case 7: return h>=k && k>=0 && (h>k || l>=0);
      case 8: return h>=0 && ((l>=h && k>h) || (l==h && k==h));
      case 9: return k>=l && l>=h && h>=0;
    }
    unreachable();
  }

  const char* condition_str() const {
    switch (idx) {
      case 0: return "l>0 or (l=0 and (h>0 or (h=0 and k>=0)))";
      case 1: return "k>=0 and (l>0 or (l=0 and h>=0))";
      case 2: return "h>=0 and k>=0 and l>=0";
      case 3: return "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))";
      case 4: return "h>=k and k>=0 and l>=0";
      case 5: return "(h>=0 and k>0) or (h=0 and k=0 and l>=0)";
      case 6: return "h>=k and k>=0 and (k>0 or l>=0)";
      case 7: return "h>=k and k>=0 and (h>k or l>=0)";
      case 8: return "h>=0 and ((l>=h and k>h) or (l=h and k=h))";
      case 9: return "k>=l and l>=h and h>=0";
    }
    unreachable();
  }

  // Returns hkl in asu and MTZ ISYM - 2*n-1 for reflections in the positive
  // asu (I+ of a Friedel pair), 2*n for reflections in the negative asu (I-).
  std::pair<Op::Miller,int> to_asu(const Op::Miller& hkl, const GroupOps& gops) const {
    int isym = 0;
    for (const Op& op : gops.sym_ops) {
      ++isym;
      Op::Miller new_hkl = op.apply_to_hkl_without_division(hkl);
      if (is_in(new_hkl))
        return {Op::divide_hkl_by_DEN(new_hkl), isym};
      ++isym;
      Op::Miller negated_new_hkl{{-new_hkl[0], -new_hkl[1], -new_hkl[2]}};
      if (is_in(negated_new_hkl))
        return {Op::divide_hkl_by_DEN(negated_new_hkl), isym};
    }
    fail("Oops, maybe inconsistent GroupOps?");
  }
};

} // namespace gemmi

namespace std {
template<> struct hash<gemmi::Op> {
  size_t operator()(const gemmi::Op& op) const {
    size_t h = 0;
    for (int i = 0; i != 3; ++i)
      for (int j = 0; j != 3; ++j)
        h = (h << 2) ^ (op.rot[i][j] + 1);
    for (int i = 0; i != 3; ++i)
      h = (h << 5) ^ op.tran[i];
    return h;
  }
};
} // namespace std

#ifdef __clang__
# pragma clang diagnostic pop  // ignored -Wmissing-braces
#endif

#endif
