// Copyright 2020 Global Phasing Ltd.
//
// AsuData for storing reflection data.

#ifndef GEMMI_ASUDATA_HPP_
#define GEMMI_ASUDATA_HPP_

#include <complex>       // for arg, abs
#include <tuple>         // for tie
#include <algorithm>     // for sort, is_sorted
#include "unitcell.hpp"
#include "symmetry.hpp"

namespace gemmi {

struct ComplexCorrelation {
  int n = 0;
  double sum_xx = 0.;
  double sum_yy = 0.;
  std::complex<double> sum_xy = 0.;
  std::complex<double> mean_x = 0.;
  std::complex<double> mean_y = 0.;
  void add_point(std::complex<double> x, std::complex<double> y) {
    ++n;
    double inv_n = 1.0 / n;
    double weight = (n - 1.0) * inv_n;
    std::complex<double> dx = x - mean_x;
    std::complex<double> dy = y - mean_y;
    sum_xx += weight * std::norm(dx);
    sum_yy += weight * std::norm(dy);
    sum_xy += weight * (dx * std::conj(dy));
    mean_x += dx * inv_n;
    mean_y += dy * inv_n;
  }
  void add_point(std::complex<float> x, std::complex<float> y) {
    add_point(std::complex<double>(x), std::complex<double>(y));
  }
  std::complex<double> coefficient() const { return sum_xy / std::sqrt(sum_xx * sum_yy); }
  double mean_ratio() const { return std::abs(mean_y) / std::abs(mean_x); }
};


// pre: both are sorted
template<typename Func, typename T>
void for_matching_reflections(const std::vector<T>& a,
                              const std::vector<T>& b,
                              const Func& func) {
  auto r1 = a.begin();
  auto r2 = b.begin();
  while (r1 != a.end() && r2 != b.end()) {
    if (r1->hkl == r2->hkl) {
      func(*r1, *r2);
      ++r1;
      ++r2;
    } else if (std::tie(r1->hkl[0], r1->hkl[1], r1->hkl[2]) <
               std::tie(r2->hkl[0], r2->hkl[1], r2->hkl[2])) {
      ++r1;
    } else {
      ++r2;
    }
  }
}

// pre: both are sorted
template<typename T>
Correlation calculate_hkl_value_correlation(const std::vector<T>& a,
                                            const std::vector<T>& b) {
  Correlation cor;
  for_matching_reflections(a, b, [&cor](const T& x, const T& y) {
      cor.add_point(x.value, y.value);
  });
  return cor;
}

// pre: both are sorted
template<typename T>
ComplexCorrelation calculate_hkl_complex_correlation(const std::vector<T>& a,
                                                     const std::vector<T>& b) {
  ComplexCorrelation cor;
  for_matching_reflections(a, b, [&cor](const T& x, const T& y) {
      cor.add_point(x.value, y.value);
  });
  return cor;
}

// pre: both are sorted
template<typename T>
int count_equal_values(const std::vector<T>& a, const std::vector<T>& b) {
  int count = 0;
  for_matching_reflections(a, b, [&count](const T& x, const T& y) {
      if (x.value == y.value)
        ++count;
  });
  return count;
}

template<typename T>
struct HklValue {
  Miller hkl;
  T value;

  bool operator<(const Miller& m) const {
    return std::tie(hkl[0], hkl[1], hkl[2]) < std::tie(m[0], m[1], m[2]);
  }
  bool operator<(const HklValue& o) const { return operator<(o.hkl); }
};

template<typename T>
struct ValueSigma {
  using value_type = T;
  T value, sigma;

  bool operator==(const ValueSigma& o) const {
    return value == o.value && sigma == o.sigma;
  }
};

namespace impl {
template<typename T>
void move_to_asu(const GroupOps&, const Miller& hkl, int, HklValue<T>& hkl_value) {
  hkl_value.hkl = hkl;
}

template<typename R>
void move_to_asu(const GroupOps& gops, const Miller& hkl, int isym,
                 HklValue<std::complex<R>>& hkl_value) {
  hkl_value.hkl = hkl;
  // cf. Mtz::ensure_asu()
  const Op& op = gops.sym_ops[(isym - 1) / 2];
  double shift = op.phase_shift(hkl);
  if (shift != 0) {
    if (isym % 2 == 0)
      shift = -shift;
    double phase = std::arg(hkl_value.value) + shift;
    hkl_value.value = std::polar(std::abs(hkl_value.value), (R)phase);
  }
}
} // namespace impl

template<typename T>
struct AsuData {
  std::vector<HklValue<T>> v;
  UnitCell unit_cell_;
  const SpaceGroup* spacegroup_ = nullptr;
  // function defining FPhiProxy interface
  size_t stride() const { return 1; }
  size_t size() const { return v.size(); }
  Miller get_hkl(size_t n) const { return v[n].hkl; }
  double get_f(size_t n) const { return std::abs(v[n].value); }
  double get_phi(size_t n) const { return std::arg(v[n].value); }
  const UnitCell& unit_cell() const { return unit_cell_; }
  const SpaceGroup* spacegroup() const { return spacegroup_; }
  void ensure_sorted() {
    if (!std::is_sorted(v.begin(), v.end()))
      std::sort(v.begin(), v.end());
  }

  void ensure_asu() {
    if (!spacegroup_)
      fail("AsuData::ensure_asu(): space group not set");
    GroupOps gops = spacegroup_->operations();
    ReciprocalAsu asu(spacegroup_);
    for (HklValue<T>& hkl_value : v) {
      const Miller& hkl = hkl_value.hkl;
      if (asu.is_in(hkl))
        continue;
      auto result = asu.to_asu(hkl, gops);
      impl::move_to_asu(gops, result.first, result.second, hkl_value);
    }
  }

  // load values from one column
  template<typename DataProxy>
  void load_values(const DataProxy& proxy, const std::string& label,
                   bool as_is=false) {
    std::size_t col = proxy.column_index(label);
    unit_cell_ = proxy.unit_cell();
    spacegroup_ = proxy.spacegroup();
    for (size_t i = 0; i < proxy.size(); i += proxy.stride()) {
      auto num = proxy.get_num(i + col);
      if (!std::isnan(num))
        v.push_back({proxy.get_hkl(i), (T)num});
    }
    if (!as_is) {
      ensure_asu();
      ensure_sorted();
    }
  }

  // load values from two or more columns
  template<int N, typename DataProxy>
  void load_values(const DataProxy& proxy, const std::array<std::string,N>& labels,
                   bool as_is=false) {
    std::array<std::size_t, N> cols;
    for (int i = 0; i < N; ++i)
      cols[i] = proxy.column_index(labels[i]);
    unit_cell_ = proxy.unit_cell();
    spacegroup_ = proxy.spacegroup();
    using Val = typename T::value_type;
    for (size_t i = 0; i < proxy.size(); i += proxy.stride()) {
      std::array<Val, N> nums;
      for (int j = 0; j < N; ++j)
        nums[j] = (Val) proxy.get_num(i + cols[j]);
      if (!std::any_of(nums.begin(), nums.end(), [](Val f) { return std::isnan(f); })) {
        v.emplace_back();
        v.back().hkl = proxy.get_hkl(i);
        set_value_from_array(v.back().value, nums);
      }
    }
    if (!as_is) {
      ensure_asu();
      ensure_sorted();
    }
  }

private:
  // for T being a number, std::array and std::complex, respectively:
  static void set_value_from_array(T& val, const std::array<T,1>& nums) { val = nums[0]; }
  static void set_value_from_array(T& val, const T& nums) { val = nums; }
  template<typename R>
  static void set_value_from_array(std::complex<R>& val, const std::array<R,2>& nums) {
    R theta = (R)rad(nums[1]);
    val = {nums[0] * std::cos(theta), nums[0] * std::sin(theta)};
  }
  template<typename R>
  static void set_value_from_array(ValueSigma<R>& val, const std::array<R,2>& nums) {
    val.value = nums[0];
    val.sigma = nums[1];
  }
};

template<typename T, int N, typename Data>
AsuData<T> make_asu_data(const Data& data, const std::array<std::string,N>& labels,
                         bool as_is=false) {
  AsuData<T> asu_data;
  asu_data.template load_values<N>(data_proxy(data), labels, as_is);
  return asu_data;
}
template<typename T, typename Data>
AsuData<T> make_asu_data(const Data& data, const std::string& label, bool as_is) {
  AsuData<T> asu_data;
  asu_data.load_values(data_proxy(data), label, as_is);
  return asu_data;
}

} // namespace gemmi
#endif
