// Copyright 2018 Global Phasing Ltd.
//
// Locale-independent functions that convert string to integer,
// equivalents of standard isspace and isdigit, and a few helper functions.
//
// This file is named similarly to the standard functions atoi() and atof().
// But the functions here are not meant to be equivalent to the standard
// library functions. They are locale-independent (a good thing when reading
// numbers from files). They don't set errno, don't signal overflow and
// underflow. Due to the limited scope these functions tend to be faster
// than the standard-library ones.

#ifndef GEMMI_ATOX_HPP_
#define GEMMI_ATOX_HPP_

#include <stdexcept>  // for invalid_argument
#include <string>

namespace gemmi {

// equivalent of std::isspace for C locale (no handling of EOF)
inline bool is_space(char c) {
  static const std::uint8_t table[256] = { // 1 for 9-13 and 32
    0,0,0,0,0,0,0,0, 0,1,1,1,1,1,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0
  };
  return table[(std::uint8_t)c] != 0;
}

// equivalent of std::isblank for C locale (no handling of EOF)
inline bool is_blank(char c) {
  return c == ' ' || c == '\t';
}

// equivalent of std::isdigit for C locale (no handling of EOF)
inline bool is_digit(char c) {
  return c >= '0' && c <= '9';
}

inline const char* skip_blank(const char* p) {
  if (p)
    while (is_blank(*p))
      ++p;
  return p;
}

inline const char* skip_word(const char* p) {
  if (p)
    while (*p != '\0' && !is_space(*p))
      ++p;
  return p;
}

inline std::string read_word(const char* line) {
  line = skip_blank(line);
  return std::string(line, skip_word(line));
}

inline std::string read_word(const char* line, const char** endptr) {
  line = skip_blank(line);
  *endptr = skip_word(line);
  return std::string(line, *endptr);
}

// no checking for overflow
inline int string_to_int(const char* p, bool checked, size_t length=0) {
  int mult = -1;
  int n = 0;
  size_t i = 0;
  while ((length == 0 || i < length) && is_space(p[i]))
    ++i;
  if (p[i] == '-') {
    mult = 1;
    ++i;
  } else if (p[i] == '+') {
    ++i;
  }
  bool has_digits = false;
  // use negative numbers because INT_MIN < -INT_MAX
  for (; (length == 0 || i < length) && is_digit(p[i]); ++i) {
    n = n * 10 - (p[i] - '0');
    has_digits = true;
  }
  if (checked) {
    while ((length == 0 || i < length) && is_space(p[i]))
      ++i;
    if (!has_digits || p[i] != '\0')
      throw std::invalid_argument("not an integer: " +
                                  std::string(p, length ? length : i+1));
  }
  return mult * n;
}

inline int string_to_int(const std::string& str, bool checked) {
  return string_to_int(str.c_str(), checked);
}

inline int simple_atoi(const char* p, const char** endptr=nullptr) {
  int mult = -1;
  int n = 0;
  while (is_space(*p))
    ++p;
  if (*p == '-') {
    mult = 1;
    ++p;
  } else if (*p == '+') {
    ++p;
  }
  for (; is_digit(*p); ++p)
    n = n * 10 - (*p - '0'); // use negative numbers because INT_MIN < -INT_MAX
  if (endptr)
    *endptr = p;
  return mult * n;
}

inline int no_sign_atoi(const char* p, const char** endptr=nullptr) {
  int n = 0;
  while (is_space(*p))
    ++p;
  for (; is_digit(*p); ++p)
    n = n * 10 + (*p - '0');
  if (endptr)
    *endptr = p;
  return n;
}

} // namespace gemmi
#endif
