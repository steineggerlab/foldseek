// Copyright 2018 Global Phasing Ltd.
//
// Classes for iterating files in a directory tree, top-down,
// in an alphabetical order.  It wraps the tinydir library (as we cannot
// depend on C++17 <filesystem> yet).

// DirWalk<> iterates through all files and directories.
// CifWalk yields only cif files (either files that end with .cif or .cif.gz,
// or files that look like SF mmCIF files from wwPDB, e.g. r3aaasf.ent.gz).
// It's good for traversing a local copy of the wwPDB archive.
// PdbWalk: .pdb or .ent (optionally with .gz) except r????sf.ent
// CoorFileWalk: .cif, .pdb or .ent (optionally with .gz)
//               except r????sf.ent and *-sf.cif
//
// Usage:
//   for (const std::string& file : gemmi::DirWalk<>(top_dir))
//     do_something(file);
// or
//   for (const std::string& file : gemmi::CifWalk(top_dir))
//     do_something(file);
// You should also catch std::runtime_error.

#ifndef GEMMI_DIRWALK_HPP_
#define GEMMI_DIRWALK_HPP_

#include <stdexcept>  // for runtime_error
#include <string>
#include <vector>
#include <cassert>
#if defined(_MSC_VER) && !defined(NOMINMAX)
# define NOMINMAX
#endif
#include "third_party/tinydir.h"

#include "util.hpp"  // for giends_with
#include "fail.hpp"  // for sys_fail
#if defined(_WIN32) && defined(_UNICODE)
 #include "utf.hpp"
#endif

namespace gemmi {

inline std::string as_utf8(const _tinydir_char_t* path) {
#if defined(_WIN32) && defined(_UNICODE)
  return wchar_to_UTF8(path);
#else
  return path;
#endif
}

// linear-time glob matching: https://research.swtch.com/glob
inline bool glob_match(const std::string& pattern, const std::string& str) {
  size_t pat_next = 0;
  size_t str_next = std::string::npos;
  size_t pat_pos = 0;
  size_t str_pos = 0;
  while (pat_pos < pattern.size() || str_pos < str.size()) {
    if (pat_pos < pattern.size()) {
      char c = pattern[pat_pos];
      if (c == '*') {
        pat_next = pat_pos;
        str_next = str_pos + 1;
        pat_pos++;
        continue;
      } else if (str_pos < str.size() && (c == '?' || c == str[str_pos])) {
        pat_pos++;
        str_pos++;
        continue;
      }
    }
    if (str_next > str.size())
      return false;
    pat_pos = pat_next;
    str_pos = str_next;
  }
  return true;
}


namespace impl {
// the SF mmCIF files from PDB have names such as
// divided/structure_factors/aa/r3aaasf.ent.gz
inline bool is_rxsf_ent_filename(const std::string& filename) {
  return filename[0] == 'r' && giends_with(filename, "sf.ent")
         && filename.find('.') >= 4;
}

struct IsMmCifFile { // actually we don't know what kind of cif file it is
  static bool check(const std::string& filename) {
    return giends_with(filename, ".cif") || giends_with(filename, ".mmcif");
  }
};

struct IsCifFile {
  static bool check(const std::string& filename) {
    return giends_with(filename, ".cif") || is_rxsf_ent_filename(filename);
  }
};

struct IsPdbFile {
  static bool check(const std::string& filename) {
    return giends_with(filename, ".pdb") ||
           (giends_with(filename, ".ent") && !is_rxsf_ent_filename(filename));
  }
};

struct IsCoordinateFile {
  static bool check(const std::string& filename) {
    // the SF mmCIF files from RCSB website have names such as 3AAA-sf.cif
    return IsPdbFile::check(filename) ||
           (IsMmCifFile::check(filename) && !giends_with(filename, "-sf.cif"));
  }
};

struct IsAnyFile {
  static bool check(const std::string&) { return true; }
};

struct IsMatchingFile {
  bool check(const std::string& filename) const {
    return glob_match(pattern, filename);
  }
  std::string pattern;
};

} // namespace impl


template<bool FileOnly=true, typename Filter=impl::IsAnyFile>
class DirWalk {
public:
  explicit DirWalk(const char* path) {
#if defined(_WIN32) && defined(_UNICODE)
    std::wstring str = UTF8_to_wchar(path);
    const _tinydir_char_t* xpath = str.c_str();
#else
    const char* xpath = path;
#endif
    if (tinydir_file_open(&top_, xpath) == -1)
      sys_fail("Cannot open " + std::string(path));
  }
  explicit DirWalk(const std::string& path) : DirWalk(path.c_str()) {}
  ~DirWalk() {
    for (auto& d : dirs_)
      tinydir_close(&d.second);
  }
  void push_dir(size_t cur_pos, const _tinydir_char_t* path) {
    dirs_.emplace_back();
    dirs_.back().first = cur_pos;
    if (tinydir_open_sorted(&dirs_.back().second, path) == -1)
      sys_fail("Cannot open directory " + as_utf8(path));
  }
  size_t pop_dir() {
    assert(!dirs_.empty());
    size_t old_pos = dirs_.back().first;
    tinydir_close(&dirs_.back().second);
    dirs_.pop_back();
    return old_pos;
  }

  struct Iter {
    DirWalk& walk;
    size_t cur;

    const tinydir_dir& get_dir() const { return walk.dirs_.back().second; }

    const tinydir_file& get() const {
      if (walk.dirs_.empty())
        return walk.top_;
      assert(cur < get_dir().n_files);
      return get_dir()._files[cur];
    }

    std::string operator*() const { return as_utf8(get().path); }

    // checks for "." and ".."
    bool is_special(const _tinydir_char_t* name) const {
      return name[0] == '.' && (name[1] == '\0' ||
                                (name[1] == '.' && name[2] == '\0'));
    }

    size_t depth() const { return walk.dirs_.size(); }

    void next() { // depth first
      const tinydir_file& tf = get();
      if (tf.is_dir) {
        walk.push_dir(cur, tf.path);
        cur = 0;
      } else {
        cur++;
      }
      while (!walk.dirs_.empty()) {
        if (cur == get_dir().n_files)
          cur = walk.pop_dir() + 1;
        else if (is_special(get_dir()._files[cur].name))
          cur++;
        else
          break;
      }
    }

    void operator++() {
      for (;;) {
        next();
        const tinydir_file& f = get();
        if ((!FileOnly && f.is_dir)
            || (!f.is_dir && walk.filter.check(as_utf8(f.name)))
            || walk.is_single_file()
            || (depth() == 0 && cur == 1))
          break;
      }
    }

    // == and != is used only to compare with end()
    bool operator==(const Iter& o) const { return depth()==0 && cur == o.cur; }
    bool operator!=(const Iter& o) const { return !operator==(o); }
  };

  Iter begin() {
    Iter it{*this, 0};
    if (FileOnly && !is_single_file()) // i.e. the top item is a directory
      ++it;
    return it;
  }

  Iter end() { return Iter{*this, 1}; }
  bool is_single_file() { return !top_.is_dir; }

private:
  friend struct Iter;
  tinydir_file top_;
  std::vector<std::pair<size_t, tinydir_dir>> dirs_;
protected:
  Filter filter;
};

using CifWalk = DirWalk<true, impl::IsCifFile>;
using MmCifWalk = DirWalk<true, impl::IsMmCifFile>;
using PdbWalk = DirWalk<true, impl::IsPdbFile>;
using CoorFileWalk = DirWalk<true, impl::IsCoordinateFile>;

struct GlobWalk : public DirWalk<true, impl::IsMatchingFile> {
  GlobWalk(const std::string& path, const std::string& glob) : DirWalk(path) {
    filter.pattern = glob;
  }
};

} // namespace gemmi
#endif
