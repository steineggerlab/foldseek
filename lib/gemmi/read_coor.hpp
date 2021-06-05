// Copyright 2021 Global Phasing Ltd.
//
// Functions for reading possibly gzipped coordinate files.
// Trivial wrappers that can make compilation faster.
//
// The definitions are guarded by #ifdef *_IMPLEMENTATION.
// If you use this header in your program, you need in exactly one file:
// #define GEMMI_READ_COOR_IMPLEMENTATION
// #include <gemmi/read_coor.cpp>
//
// See src/input.cpp for an example.

#ifndef GEMMI_READ_COOR_HPP_
#define GEMMI_READ_COOR_HPP_

#include "model.hpp"  // for Structure

namespace gemmi {

Structure read_structure_gz(const std::string& path, CoorFormat format=CoorFormat::Unknown);

Structure read_pdb_gz(const std::string& path, PdbReadOptions options=PdbReadOptions());

CoorFormat coor_format_from_ext_gz(const std::string& path);

} // namespace gemmi

#endif


#ifdef GEMMI_READ_COOR_IMPLEMENTATION

#include "mmread.hpp" // for read_structure
#include "pdb.hpp"    // for read_pdb
#include "gz.hpp"     // for MaybeGzipped

namespace gemmi {

Structure read_pdb_gz(const std::string& path, PdbReadOptions options) {
  return read_pdb(MaybeGzipped(path), options);
}

Structure read_structure_gz(const std::string& path, CoorFormat format) {
  return read_structure(MaybeGzipped(path), format);
}

CoorFormat coor_format_from_ext_gz(const std::string& path) {
  return coor_format_from_ext(MaybeGzipped(path).basepath());
}

} // namespace gemmi
#endif
