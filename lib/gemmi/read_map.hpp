// Copyright 2021 Global Phasing Ltd.
//
// Functions for reading possibly gzipped CCP4 map files.
// Trivial wrappers that can make compilation faster.
//
// The definitions are guarded by #ifdef *_IMPLEMENTATION.
// If you use this header in your program, you need in exactly one file:
// #define GEMMI_READ_MAP_IMPLEMENTATION
// #include <gemmi/read_map.hpp>

#ifndef GEMMI_GZREAD_MAP_HPP_
#define GEMMI_GZREAD_MAP_HPP_

#include "ccp4.hpp"  // for Ccp4

namespace gemmi {

Ccp4<float> read_ccp4_map(const std::string& path, bool setup);
Ccp4<int8_t> read_ccp4_mask(const std::string& path, bool setup);

} // namespace gemmi

#endif


#ifdef GEMMI_READ_MAP_IMPLEMENTATION

#include "gz.hpp"  // for MaybeGzipped

namespace gemmi {

Ccp4<float> read_ccp4_map(const std::string& path, bool setup) {
  Ccp4<float> ccp4;
  ccp4.read_ccp4(MaybeGzipped(path));
  if (setup)
    ccp4.setup(GridSetup::Full, NAN);
  return ccp4;
}

Ccp4<int8_t> read_ccp4_mask(const std::string& path, bool setup) {
  Ccp4<int8_t> ccp4;
  ccp4.read_ccp4(MaybeGzipped(path));
  if (setup)
    ccp4.setup(GridSetup::Full, -1);
  return ccp4;
}

} // namespace gemmi

#endif
