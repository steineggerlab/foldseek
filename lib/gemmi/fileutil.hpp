// Copyright 2018 Global Phasing Ltd.
//
// File-related utilities.

#ifndef GEMMI_FILEUTIL_HPP_
#define GEMMI_FILEUTIL_HPP_

#include <cctype>    // for isdigit, isalnum
#include <cstdio>    // for FILE, fopen, fclose
#include <cstdlib>   // getenv
#include <cstring>   // strlen
#include <initializer_list>
#include <memory>    // for unique_ptr
#include <string>
#include "fail.hpp"  // for fail, sys_fail
#include "util.hpp"  // for to_lower

#if defined(_WIN32) && !defined(GEMMI_USE_FOPEN)
#include "utf.hpp"
#endif

namespace gemmi {

// strip directory and suffixes from filename
inline std::string path_basename(const std::string& path,
                                 std::initializer_list<const char*> exts) {
  size_t pos = path.find_last_of("\\/");
  std::string basename = pos == std::string::npos ? path : path.substr(pos + 1);
  for (const char* ext : exts) {
    size_t len = std::strlen(ext);
    if (basename.size() > len &&
        basename.compare(basename.length() - len, len, ext, len) == 0)
      basename.resize(basename.length() - len);
  }
  return basename;
}

// file operations
typedef std::unique_ptr<std::FILE, decltype(&std::fclose)> fileptr_t;

inline fileptr_t file_open(const char* path, const char* mode) {
  std::FILE* file;
#if defined(_WIN32) && !defined(GEMMI_USE_FOPEN)
  std::wstring wpath = UTF8_to_wchar(path);
  std::wstring wmode = UTF8_to_wchar(mode);
  if ((file = ::_wfopen(wpath.c_str(), wmode.c_str())) == nullptr)
#else
  if ((file = std::fopen(path, mode)) == nullptr)
#endif
    sys_fail(std::string("Failed to open ") + path +
             (*mode == 'w' ? " for writing" : ""));
  return fileptr_t(file, &std::fclose);
}

// helper function for treating "-" as stdin or stdout
inline fileptr_t file_open_or(const char* path, const char* mode,
                              std::FILE* dash_stream) {
  if (path[0] == '-' && path[1] == '\0')
    return fileptr_t(dash_stream, [](std::FILE*) { return 0; });
  return file_open(path, mode);
}

inline std::size_t file_size(std::FILE* f, const std::string& path) {
  if (std::fseek(f, 0, SEEK_END) != 0)
    sys_fail(path + ": fseek failed");
  long length = std::ftell(f);
  if (length < 0)
    sys_fail(path + ": ftell failed");
  if (std::fseek(f, 0, SEEK_SET) != 0)
    sys_fail(path + ": fseek failed");
  return length;
}

inline bool is_pdb_code(const std::string& str) {
  return str.length() == 4 && std::isdigit(str[0]) && std::isalnum(str[1]) &&
                              std::isalnum(str[2]) && std::isalnum(str[3]);
}

// Call it after checking the code with gemmi::is_pdb_code(code).
// The convention for $PDB_DIR is the same as in BioJava, see the docs.
// type is the requested file type: 'M' for mmCIF or 'P' for PDB.
inline std::string expand_pdb_code_to_path(const std::string& code, char type) {
  std::string path;
  if (const char* pdb_dir = std::getenv("PDB_DIR")) {
    std::string lc = to_lower(code);
    path = pdb_dir;
    path += "/structures/divided/";
    path += (type == 'M' ? "mmCIF/" : "pdb/");
    path += lc.substr(1, 2) + "/";
    path += (type == 'M' ?  lc + ".cif.gz" : "pdb" + lc + ".ent.gz");
  }
  return path;
}

// type must be 'M' for mmCIF or 'P' for PDB
inline std::string expand_if_pdb_code(const std::string& input, char type='M') {
  std::string path;
  if (is_pdb_code(input)) {
    path = gemmi::expand_pdb_code_to_path(input, type);
    if (path.empty())
      fail(input + " is a PDB code, but $PDB_DIR is not set.");
  } else {
    path = input;
  }
  return path;
}

// helper function for working with binary files
inline bool is_little_endian() {
  std::uint32_t x = 1;
  return *reinterpret_cast<char *>(&x) == 1;
}

inline void swap_two_bytes(void* start) {
  char* bytes = static_cast<char*>(start);
  std::swap(bytes[0], bytes[1]);
}

inline void swap_four_bytes(void* start) {
  char* bytes = static_cast<char*>(start);
  std::swap(bytes[0], bytes[3]);
  std::swap(bytes[1], bytes[2]);
}

} // namespace gemmi
#endif
