// Copyright 2021 Global Phasing Ltd.
//
// Functions for reading possibly gzipped CIF files.
// Trivial wrappers that can make compilation faster.
//
// The definitions are guarded by #ifdef *_IMPLEMENTATION.
// If you use this header in your program, you need in exactly one file:
// #define GEMMI_READ_CIF_IMPLEMENTATION
// #include <gemmi/read_cif.hpp>

#ifndef GEMMI_READ_CIF_HPP_
#define GEMMI_READ_CIF_HPP_

#include "cifdoc.hpp" // for Document
#include "input.hpp"  // for CharArray

namespace gemmi {

cif::Document read_cif_gz(const std::string& path);
cif::Document read_mmjson_gz(const std::string& path);
CharArray read_into_buffer_gz(const std::string& path);
cif::Document read_cif_from_buffer(const CharArray& buffer, const char* name);

inline cif::Document read_cif_or_mmjson_gz(const std::string& path) {
  if (giends_with(path, "json") || giends_with(path, "js"))
    return read_mmjson_gz(path);
  return read_cif_gz(path);
}

} // namespace gemmi

#endif


#ifdef GEMMI_READ_CIF_IMPLEMENTATION

#include "cif.hpp"    // for cif::read
#include "json.hpp"   // for cif::read_mmjson
#include "gz.hpp"     // for MaybeGzipped

namespace gemmi {

cif::Document read_cif_gz(const std::string& path) {
  return cif::read(MaybeGzipped(path));
}

cif::Document read_mmjson_gz(const std::string& path) {
  return cif::read_mmjson(MaybeGzipped(path));
}

CharArray read_into_buffer_gz(const std::string& path) {
  return cif::read_into_buffer(MaybeGzipped(path));
}

cif::Document read_cif_from_buffer(const CharArray& buffer, const char* name) {
  return cif::read_memory(buffer.data(), buffer.size(), name);
}

} // namespace gemmi

#endif
