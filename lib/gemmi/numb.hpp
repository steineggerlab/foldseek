// Copyright 2017 Global Phasing Ltd.
//
// Utilities for parsing CIF numbers (the CIF spec calls it 'numb').
//
// Numb - the numeric type in CIF - is a number with optional
// standard uncertainty (s.u.) in brackets: 1.23(8).
// Mmcif file do not use s.u. though - they define own numeric categories.

#ifndef GEMMI_NUMB_HPP_
#define GEMMI_NUMB_HPP_

#include <cmath>   // for NAN
#include <string>
#include "third_party/fast_float/fast_float.h"

namespace gemmi {
namespace cif {

inline double as_number(const std::string& s, double nan=NAN) {
  const char* start = s.data();
  const char* end = s.data() + s.size();
  if (*start == '+')
    ++start;
  // NaN, Inf and -Inf are not allowed in CIF
  char f = start[int(*start == '-')] | 0x20;
  if (f == 'i' || f == 'n')
    return nan;

  double d;
  auto result = fast_float::from_chars(start, end, d);
  if (result.ec != std::errc())
    return nan;
  if (*result.ptr == '(') {
    const char* p = result.ptr + 1;
    while (*p >= '0' && *p <= '9')
      ++p;
    if (*p == ')')
      result.ptr = p + 1;
  }
  return result.ptr == end ? d : nan;
}

inline bool is_numb(const std::string& s) {
  return !std::isnan(as_number(s));
}


// for use in templates (see also as_any() functions in cifdoc.hpp)
inline float as_any(const std::string& s, float null) {
  return (float) as_number(s, null);
}
inline double as_any(const std::string& s, double null) {
  return as_number(s, null);
}

} // namespace cif
} // namespace gemmi
#endif
