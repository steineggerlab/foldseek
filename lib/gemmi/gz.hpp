// Copyright 2017 Global Phasing Ltd.
//
// Functions for transparent reading of gzipped files. Uses zlib.

#ifndef GEMMI_GZ_HPP_
#define GEMMI_GZ_HPP_
#include <cassert>
#include <cstdio>       // fseek, ftell, fread
#include <climits>      // INT_MAX
#include <memory>
#include <string>
#include <zlib.h>
#include "fail.hpp"     // fail, sys_fail
#include "fileutil.hpp" // file_open
#include "input.hpp"    // BasicInput
#include "util.hpp"     // iends_with

namespace gemmi {

// Throws if the size is not found or if it is suspicious.
// Anything outside of the arbitrary limits from 1 to 10x of the compressed
// size looks suspicious to us.
// **This function should not be relied upon.**
// In particular, if the return values is >= 4GiB - it's only a guess.
inline size_t estimate_uncompressed_size(const std::string& path) {
  fileptr_t f = file_open(path.c_str(), "rb");
  if (std::fseek(f.get(), -4, SEEK_END) != 0)
    sys_fail("fseek() failed (empty file?): " + path);
  long pos = std::ftell(f.get());
  if (pos <= 0)
    sys_fail("ftell() failed on " + path);
  size_t gzipped_size = pos + 4;
  unsigned char buf[4];
  if (std::fread(buf, 1, 4, f.get()) != 4)
    sys_fail("Failed to read last 4 bytes of: " + path);
  unsigned orig_size = (buf[3] << 24) | (buf[2] << 16) | (buf[1] << 8) | buf[0];
  if (orig_size + 100 < gzipped_size || orig_size > 100 * gzipped_size) {
    // The size is stored as 32-bit number. If the original size exceeds 4GiB,
    // the stored number is modulo 4 GiB. So we just guess...
    if (gzipped_size > 1073741824)
      return 4294967295U + (sizeof(size_t) > 4 ? orig_size : 0);
    fail("Cannot determine uncompressed size of " + path +
         "\nWould it be " + std::to_string(gzipped_size) + " -> " +
         std::to_string(orig_size) + " bytes?");
  }
  return orig_size;
}

inline size_t big_gzread(gzFile file, void* buf, size_t len) {
  // In zlib >= 1.2.9 we could use gzfread()
  // return gzfread(buf, len, 1, f) == 1;
  size_t read_bytes = 0;
  while (len > INT_MAX) {
    int ret = gzread(file, buf, INT_MAX);
    read_bytes += ret;
    if (ret != INT_MAX)
      return read_bytes;
    len -= INT_MAX;
    buf = (char*) buf + INT_MAX;
  }
  read_bytes += gzread(file, buf, (unsigned) len);
  return read_bytes;
}

class MaybeGzipped : public BasicInput {
public:
  struct GzStream {
    gzFile f;
    char* gets(char* line, int size) { return gzgets(f, line, size); }
    int getc() { return gzgetc(f); }
    bool read(void* buf, size_t len) { return big_gzread(f, buf, len) == len; }
  };

  explicit MaybeGzipped(const std::string& path)
    : BasicInput(path), file_(nullptr) {}
  ~MaybeGzipped() {
    if (file_)
#if ZLIB_VERNUM >= 0x1235
      gzclose_r(file_);
#else
      gzclose(file_);
#endif
  }

  size_t gzread_checked(void* buf, size_t len) {
    size_t read_bytes = big_gzread(file_, buf, len);
    if (read_bytes != len && !gzeof(file_)) {
      int errnum;
      std::string err_str = gzerror(file_, &errnum);
      if (errnum == Z_ERRNO)
        sys_fail("failed to read " + path());
      if (errnum)
        fail("Error reading " + path() + ": " + err_str);
    }
    if (read_bytes > len)  // should never happen
      fail("Error reading " + path());
    return read_bytes;
  }

  bool is_compressed() const { return iends_with(path(), ".gz"); }
  std::string basepath() const {
    return is_compressed() ? path().substr(0, path().size() - 3) : path();
  }

  CharArray uncompress_into_buffer(size_t limit=0) {
    if (!is_compressed())
      return BasicInput::uncompress_into_buffer();
    size_t size = (limit == 0 ? estimate_uncompressed_size(path()) : limit);
    open();
    if (size > 3221225471)
      // if this exception is changed adjust src/cif2mtz.cpp
      fail("For now gz files above 3 GiB uncompressed are not supported.\n"
           "To read " + path() + " first uncompress it.");
    CharArray mem(size);
    size_t read_bytes = gzread_checked(mem.data(), size);
    // if the file is shorter than the size from header, adjust size
    if (read_bytes < size) {
      mem.set_size(read_bytes);  // should we call resize() here
    } else if (limit == 0) { // read_bytes == size
    // if the file is longer than the size from header, read in the rest
      int next_char;
      while (!gzeof(file_) && (next_char = gzgetc(file_)) != -1) {
        if (mem.size() > 3221225471)
          fail("For now gz files above 3 GiB uncompressed are not supported.\n"
               "To read " + path() + " first uncompress it.");
        gzungetc(next_char, file_);
        size_t old_size = mem.size();
        mem.resize(2 * old_size);
        size_t n = gzread_checked(mem.data() + old_size, old_size);
        mem.set_size(old_size + n);
      }
    }
    return mem;
  }

  GzStream get_uncompressing_stream() {
    assert(is_compressed());
    open();
#if ZLIB_VERNUM >= 0x1235
    gzbuffer(file_, 64*1024);
#endif
    return GzStream{file_};
  }

private:
  gzFile file_;

  void open() {
    file_ = gzopen(path().c_str(), "rb");
    if (!file_)
      sys_fail("Failed to gzopen " + path());
  }
};

} // namespace gemmi

#endif
