// Copyright 2020 Global Phasing Ltd.
//
// ReciprocalGrid -- grid for reciprocal space data.

#ifndef GEMMI_RECGRID_HPP_
#define GEMMI_RECGRID_HPP_

#include "asudata.hpp"
#include "grid.hpp"

namespace gemmi {

template<typename T>
struct ReciprocalGrid : GridBase<T> {
  bool half_l = false; // hkl grid that stores only l>=0
  bool has_index(int u, int v, int w) const {
    bool half_u = (half_l && this->axis_order == AxisOrder::ZYX);
    bool half_w = (half_l && this->axis_order != AxisOrder::ZYX);
    return std::abs(half_u ? u : 2 * u) < this->nu &&
           std::abs(2 * v) < this->nv &&
           std::abs(half_w ? w : 2 * w) < this->nw;
  }
  void check_index(int u, int v, int w) const {
    if (!has_index(u, v, w))
      throw std::out_of_range("ReciprocalGrid: index out of grid.");
  }
  // Similar to Grid::index_n(), but works only for -nu <= u < nu, etc.
  size_t index_n(int u, int v, int w) const {
    return this->index_q(u >= 0 ? u : u + this->nu,
                         v >= 0 ? v : v + this->nv,
                         w >= 0 ? w : w + this->nw);
  }
  size_t index_checked(int u, int v, int w) const {
    check_index(u, v, w);
    return index_n(u, v, w);
  }
  T get_value(int u, int v, int w) const {
    return this->data[index_checked(u, v, w)];
  }
  T get_value_or_zero(int u, int v, int w) const {
    return has_index(u, v, w) ? this->data[index_n(u, v, w)] : T{};
  }
  void set_value(int u, int v, int w, T x) {
    this->data[index_checked(u, v, w)] = x;
  }
  Miller to_hkl(const typename GridBase<T>::Point& point) const {
    Miller hkl{{point.u, point.v, point.w}};
    if (2 * point.u >= this->nu &&
        !(half_l && this->axis_order == AxisOrder::ZYX))
      hkl[0] -= this->nu;
    if (2 * point.v >= this->nv)
      hkl[1] -= this->nv;
    if (2 * point.w >= this->nw &&
        !(half_l && this->axis_order != AxisOrder::ZYX))
      hkl[2] -= this->nw;
    if (this->axis_order == AxisOrder::ZYX)
      std::swap(hkl[0], hkl[2]);
    return hkl;
  }

  double calculate_1_d2(const typename GridBase<T>::Point& point) const {
    return this->unit_cell.calculate_1_d2(to_hkl(point));
  }
  double calculate_d(const typename GridBase<T>::Point& point) const {
    return this->unit_cell.calculate_d(to_hkl(point));
  }

  // the result is always sorted by h,k,l
  template <typename R=T>
  AsuData<R> prepare_asu_data(double dmin=0, double unblur=0,
                              bool with_000=false, bool with_sys_abs=false,
                              bool mott_bethe=false) {
    AsuData<R> asu_data;
    if (this->axis_order == AxisOrder::ZYX)
      fail("get_asu_values(): ZYX order is not supported yet");
    int max_h = (this->nu - 1) / 2;
    int max_k = (this->nv - 1) / 2;
    int max_l = half_l ? this->nw - 1 : (this->nw - 1) / 2;
    double max_1_d2 = 0.;
    if (dmin != 0.) {
      max_1_d2 = 1. / (dmin * dmin);
      max_h = std::min(max_h, int(1. / (dmin * this->unit_cell.ar)));
      max_k = std::min(max_k, int(1. / (dmin * this->unit_cell.br)));
      max_l = std::min(max_l, int(1. / (dmin * this->unit_cell.cr)));
    }
    gemmi::ReciprocalAsu asu(this->spacegroup);
    std::unique_ptr<GroupOps> gops;
    if (!with_sys_abs && this->spacegroup)
      gops.reset(new GroupOps(this->spacegroup->operations()));
    Miller hkl;
    for (hkl[0] = -max_h; hkl[0] <= max_h; ++hkl[0]) {
      int hi = hkl[0] >= 0 ? hkl[0] : hkl[0] + this->nu;
      for (hkl[1] = -max_k; hkl[1] <= max_k; ++hkl[1]) {
        int ki = hkl[1] >= 0 ? hkl[1] : hkl[1] + this->nv;
        for (hkl[2] = (half_l ? 0 : -max_l); hkl[2] <= max_l; ++hkl[2])
          if (asu.is_in(hkl) &&
              (max_1_d2 == 0. || this->unit_cell.calculate_1_d2(hkl) < max_1_d2) &&
              (with_sys_abs || !gops->is_systematically_absent(hkl)) &&
              (with_000 || !(hkl[0] == 0 && hkl[1] == 0 && hkl[2] == 0))) {
            int li = hkl[2] >= 0 ? hkl[2] : hkl[2] + this->nw;
            asu_data.v.push_back({hkl, this->get_value_q(hi, ki, li)});
          }
      }
    }
    if (unblur != 0. || mott_bethe)
      for (HklValue<R>& hv : asu_data.v) {
        double inv_d2 = this->unit_cell.calculate_1_d2(hv.hkl);
        double mult = 1;
        if (unblur != 0)
          // cf. reciprocal_space_multiplier()
          mult = std::exp(unblur * 0.25 * inv_d2);
        if (mott_bethe)
          // cf. mott_bethe_factor
          mult *= -1. / (2 * pi() * pi() * bohrradius()) / inv_d2;
        hv.value *= static_cast<decltype(std::abs(hv.value))>(mult);
      }
    asu_data.unit_cell_ = this->unit_cell;
    asu_data.spacegroup_ = this->spacegroup;
    return asu_data;
  }
};

template<typename T> using FPhiGrid = ReciprocalGrid<std::complex<T>>;

} // namespace gemmi
#endif
