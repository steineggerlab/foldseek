// Copyright 2019 Global Phasing Ltd.
//
// Ofstream and Ifstream: wrappers around std::ofstream and std::ifstream.
//
// They offer two extra features:
//  - on MSVC supports Unicode filenames (the filename is passed in UTF-8),
//  - optionally, filename "-" can be interpretet as stdout or stderr.

#ifndef GEMMI_OFSTREAM_HPP_
#define GEMMI_OFSTREAM_HPP_

#if defined(_MSC_VER)
# include "utf.hpp"
#elif defined(_WIN32) && defined(__has_include)
# if __has_include(<filesystem>)
#  include <filesystem>
#  include "utf.hpp"
# endif
#endif
#include <fstream>
#include <memory>
#include "fail.hpp"

namespace gemmi {

template<typename T>
inline void open_stream_from_utf8_path(T& ptr, const std::string& filename) {
#if defined(_MSC_VER)
    std::wstring wfilename = UTF8_to_wchar(filename.c_str());
    ptr->open(wfilename.c_str());
#elif defined(_WIN32) && defined(__cpp_lib_filesystem)
    std::wstring wfilename = UTF8_to_wchar(filename.c_str());
    ptr->open(std::filesystem::path(wfilename));
#else
    ptr->open(filename);
#endif
}

// note: move of std::ofstream doesn't work in GCC 4.8.

struct Ofstream {
  Ofstream(const std::string& filename, std::ostream* dash=nullptr) {
    if (filename.size() == 1 && filename[0] == '-' && dash) {
      ptr_ = dash;
      return;
    }
    keeper_.reset(new std::ofstream);
    open_stream_from_utf8_path(keeper_, filename);
    if (!*keeper_)
      sys_fail("Failed to open " + filename + " for writing");
    ptr_ = keeper_.get();
  }

  std::ostream* operator->() { return ptr_; }
  std::ostream& ref() { return *ptr_; }

private:
  std::unique_ptr<std::ofstream> keeper_;
  std::ostream* ptr_;
};

struct Ifstream {
  Ifstream(const std::string& filename, std::istream* dash=nullptr) {
    if (filename.size() == 1 && filename[0] == '-' && dash) {
      ptr_ = dash;
      return;
    }
    keeper_.reset(new std::ifstream);
    open_stream_from_utf8_path(keeper_, filename);
    if (!*keeper_)
      sys_fail("Failed to open " + filename);
    ptr_ = keeper_.get();
  }

  std::istream* operator->() { return ptr_; }
  std::istream& ref() { return *ptr_; }

private:
  std::unique_ptr<std::ifstream> keeper_;
  std::istream* ptr_;
};



} // namespace gemmi
#endif
