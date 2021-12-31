// Copyright 2020 Global Phasing Ltd.
//
// Least-squares fitting - Levenberg-Marquardt method.
//
// Based on the code from fityk (but here it's under MPL 2.0).

#ifndef GEMMI_LEVMAR_HPP_
#define GEMMI_LEVMAR_HPP_

#include <cassert>
#include <cmath>      // for fabs
#include <algorithm>  // for min
#include <vector>
#include "fail.hpp"   // for fail

//#define GEMMI_DEBUG_LEVMAR

namespace gemmi {

/// This function solves a set of linear algebraic equations using
/// Gauss-Jordan elimination with partial pivoting.
///
/// A * x = b
///
/// a is n x n matrix (in vector)
/// b is vector of length n,
/// This function returns vector x[] in b[], and 1-matrix in a[].
inline void jordan_solve(std::vector<double>& a, std::vector<double>& b) {
  assert(a.size() == b.size() * b.size());
  int n = (int) b.size();
  for (int i = 0; i < n; i++) {
    // looking for a pivot element
    int maxnr = -1;
    double amax = 0;
    for (int j = i; j < n; j++) {
      double aji = std::fabs(a[n * j + i]);
      if (aji > amax) {
        maxnr = j;
        amax = aji;
      }
    }
    // handle singular matrix
    if (maxnr == -1) {
      // i-th column has only zeros.
      // If it's the same about i-th row, and b[i]==0, let x[i]==0.
      for (int j = i; j < n; j++)
        if (a[n * i + j] != 0. || b[i] != 0.)
          fail("Trying to reverse singular matrix. Column ", i, " is zeroed.");
      continue; // x[i]=b[i], b[i]==0
    }
    // interchanging rows
    if (maxnr != i) {
      for (int j = i; j < n; j++)
        std::swap(a[n * maxnr + j], a[n * i + j]);
      std::swap(b[i], b[maxnr]);
    }
    // divide by a_ii -- to get a_ii=1
    double c = 1.0 / a[i * n + i];
    for (int j = i; j < n; j++)
      a[i * n + j] *= c;
    b[i] *= c;
    // subtract -- to zero all remaining elements of this row
    for (int k = 0; k < n; k++)
      if (k != i) {
        double d = a[k * n + i];
        for (int j = i; j < n; j++)
          a[k * n + j] -= a[i * n + j] * d;
        b[k] -= b[i] * d;
      }
  }
}


#ifdef GEMMI_DEBUG_LEVMAR
inline void debug_print(const std::string& name, std::vector<double> &a) {
  fprintf(stderr, " %s:", name.c_str());
  for (double& x : a)
    fprintf(stderr, " %g", x);
  fprintf(stderr, "\n");
}
#else
inline void debug_print(const std::string&, std::vector<double>&) {}
#endif

struct LevMar {
  // termination criteria
  int eval_limit = 100;
  double lambda_limit = 1e+15;
  double stop_rel_change = 1e-5;

  // adjustable parameters (normally the default values work fine)
  double lambda_up_factor = 10;
  double lambda_down_factor = 0.1;
  double lambda_start = 0.001;

  // values set in fit() that can be inspected later
  double initial_wssr;
  int eval_count;  // number of function evaluations

  // arrays used during refinement
  std::vector<double> alpha; // matrix
  std::vector<double> beta;  // vector
  std::vector<double> temp_alpha, temp_beta; // working arrays


  template<typename Target>
  double fit(Target& target) {
    eval_count = 0;
    std::vector<double> initial_a = target.get_parameters();
    debug_print("ini", initial_a);
    initial_wssr = this->compute_wssr(target.compute_values(), target.points);
    std::vector<double> best_a = initial_a;
    size_t na = initial_a.size();

    double lambda = lambda_start;
    alpha.resize(na * na);
    beta.resize(na);

    double wssr = initial_wssr;
    this->compute_derivatives(target);

    int small_change_counter = 0;
    for (int iter = 0; ; iter++) {
      if (eval_limit > 0 && eval_count >= eval_limit)
        break;

      // prepare next parameters -> temp_beta
      temp_alpha = alpha;
      for (size_t j = 0; j < na; j++)
        temp_alpha[na * j + j] *= (1.0 + lambda);
      temp_beta = beta;

      // Matrix solution (Ax=b)  temp_alpha * da == temp_beta
      jordan_solve(temp_alpha, temp_beta);

      for (size_t i = 0; i < na; i++)
        // put new a[] into temp_beta[]
        temp_beta[i] += best_a[i];

      target.set_parameters(temp_beta);
      double new_wssr = this->compute_wssr(target.compute_values(), target.points);
#ifdef GEMMI_DEBUG_LEVMAR
      fprintf(stderr, " #%d WSSR=%.8g %+g%% (%+.4g%%) lambda=%g\n",
              iter, new_wssr, 100. * (new_wssr / initial_wssr - 1.),
              100. * (new_wssr / wssr - 1.), lambda);
      if (new_wssr < wssr)
        debug_print("", temp_beta);
#endif
      if (new_wssr < wssr) {
        double rel_change = (wssr - new_wssr) / wssr;
        wssr = new_wssr;
        best_a = temp_beta;

        if (wssr == 0)
          break;
        // termination criterium: negligible change of wssr
        if (rel_change < stop_rel_change) {
          if (++small_change_counter >= 2)
            break;
        } else {
          small_change_counter = 0;
        }
        this->compute_derivatives(target);
        lambda *= lambda_down_factor;
      } else { // worse fitting
        if (lambda > lambda_limit) // termination criterium: large lambda
          break;
        lambda *= lambda_up_factor;
      }
    }

    target.set_parameters(wssr < initial_wssr ? best_a : initial_a);
    return wssr;
  }

private:
  template<typename Target>
  void compute_derivatives(const Target& target) {
    assert(alpha.size() == beta.size() * beta.size());
    int na = (int)beta.size();
    assert(na != 0);
#ifdef GEMMI_DEBUG_LEVMAR
    check_derivatives(const_cast<Target&>(target));
#endif
    std::fill(alpha.begin(), alpha.end(), 0.0);
    std::fill(beta.begin(), beta.end(), 0.0);
    // Iterating over points is tiled to limit memory usage. It's also a little
    // faster than a single loop over all points for large number of points.
    const size_t kMaxTileSize = 1024;
    std::vector<double> dy_da;
    size_t n = target.points.size();
    for (size_t tstart = 0; tstart < n; tstart += kMaxTileSize) {
      size_t tsize = std::min(n - tstart, kMaxTileSize);
      std::vector<double> yy(tsize, 0.);
      dy_da.resize(tsize * na);
      std::fill(dy_da.begin(), dy_da.end(), 0.);
      target.compute_values_and_derivatives(tstart, tsize, yy, dy_da);
      for (size_t i = 0; i != tsize; ++i) {
        double weight = target.points[tstart + i].get_weight();
        double dy_sig = weight * (target.points[tstart + i].get_y() - yy[i]);
        double* t = &dy_da[i * na];
        for (int j = 0; j != na; ++j) {
          if (t[j] != 0) {
            t[j] *= weight;
            for (int k = j; k != -1; --k)
              alpha[na * j + k] += t[j] * t[k];
            beta[j] += dy_sig * t[j];
          }
        }
      }
    }

    // Only half of the alpha matrix was filled above. Fill the rest.
    for (int j = 1; j < na; j++)
      for (int k = 0; k < j; k++)
        alpha[na * k + j] = alpha[na * j + k];
  }

#ifdef GEMMI_DEBUG_LEVMAR
  template<typename Target>
  void check_derivatives(Target& target) {
    assert(!beta.empty());
    assert(alpha.size() == beta.size() * beta.size());
    std::vector<double> yy(target.points.size(), 0.);
    std::vector<double> dy_da(target.points.size() * beta.size());
    target.compute_values_and_derivatives(0, target.points.size(), yy, dy_da);
    std::vector<double> yy2 = target.compute_values();
    assert(yy.size() == yy2.size());
    for (size_t i = 0; i != yy.size(); ++i) {
      double m = std::max(std::fabs(yy[i]), std::fabs(yy2[i]));
      if (m > 1e-5 && std::fabs(yy[i] - yy2[i]) > 1e-6 * m)
        fprintf(stderr, "!! value %zu: %g != %g\n", i, yy[i], yy2[i]);
    }
    const double numerical_h = 1e-3;
    const double small_number = 1e-6; // prevents h==0
    std::vector<double> aa = target.get_parameters();
    assert(aa.size() == beta.size());
    for (size_t k = 0; k < aa.size(); k++) {
      double acopy = aa[k];
      double h = std::max(std::fabs(acopy), small_number) * numerical_h;
      aa[k] -= h;
      target.set_parameters(aa);
      std::vector<double> y_left = target.compute_values();
      aa[k] = acopy + h;
      target.set_parameters(aa);
      std::vector<double> y_right = target.compute_values();
      aa[k] = acopy;
      for (size_t i = 0; i != target.points.size(); ++i) {
        double symbolic = dy_da[i * aa.size() + k];
        double numeric = (y_right[i] - y_left[i]) / (2 * h);
        double m = std::max(std::fabs(symbolic), std::fabs(numeric));
        if (m > 1e-3 && std::fabs(symbolic - numeric) > 0.02 * m)
          fprintf(stderr, "!! dy[%zu]/da[%zu]: %g != %g\n", i, k, symbolic, numeric);
        if (i == 30)
          break;
      }
    }
    target.set_parameters(aa);
  }
#endif

  template<typename Point>
  double compute_wssr(const std::vector<double>& yy, const std::vector<Point>& data) {
    // long double here notably increases the accuracy of calculations
    long double wssr = 0;
    for (int j = 0; j < (int)yy.size(); j++) {
      double dy = data[j].get_weight() * (data[j].get_y() - yy[j]);
      wssr += dy * dy;
    }
    ++eval_count;
    return (double) wssr;
  }
};

} // namespace gemmi
#endif
