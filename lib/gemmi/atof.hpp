// Copyright 2020 Global Phasing Ltd.
//
// Functions that convert string to floating-point number ignoring locale.
// Simple wrappers around fastfloat::from_chars().

#ifndef GEMMI_ATOF_HPP_
#define GEMMI_ATOF_HPP_

#include "atox.hpp"   // for is_space
#include "third_party/fast_float/fast_float.h"

namespace gemmi {

using fast_float::from_chars_result;

inline from_chars_result fast_from_chars(const char* start, const char* end, double& d) {
  while (start < end && is_space(*start))
    ++start;
  if (start < end && *start == '+')
    ++start;
  return fast_float::from_chars(start, end, d);
}

inline from_chars_result fast_from_chars(const char* start, double& d) {
  while (is_space(*start))
    ++start;
  if (*start == '+')
    ++start;
  return fast_float::from_chars(start, start + std::strlen(start), d);
}

inline double fast_atof(const char* p, const char** endptr=nullptr) {
  double d = 0;
  auto result = fast_from_chars(p, d);
  if (endptr)
    *endptr = result.ptr;
  return d;
}

} // namespace gemmi
#endif
