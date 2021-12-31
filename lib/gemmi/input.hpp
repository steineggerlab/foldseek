// Copyright 2018 Global Phasing Ltd.
//
// Input abstraction.
// Used to decouple file reading and uncompression.

#ifndef GEMMI_INPUT_HPP_
#define GEMMI_INPUT_HPP_

#include <cassert>
#include <cstddef> // for ptrdiff_t
#include <cstdio>  // for FILE, fseek, fread
#include <cstdlib> // for malloc, realloc
#include <cstring> // for memchr
#include <memory>  // for unique_ptr
#include <string>
#include "fail.hpp"  // for unreachable

namespace gemmi {

// provides the same methods as GzStream
struct FileStream {
  std::FILE* f;
  // used in pdb.hpp
  char* gets(char* line, int size) { return std::fgets(line, size, f); }
  int getc() { return std::fgetc(f); }
  // used in ccp4.hpp
  bool read(void* buf, size_t len) { return std::fread(buf, len, 1, f) == 1; }

  // used in mtz.hpp
  std::string read_rest() {
    std::string ret;
    int c = std::fgetc(f);
    if (c != EOF) {
      ret += (char)c;
      char buf[512];
      for (;;) {
        size_t n = std::fread(buf, 1, sizeof(buf), f);
        ret.append(buf, n);
        if (n != sizeof(buf))
          break;
      }
    }
    return ret;
  }

  bool seek(std::ptrdiff_t offset) {
#if defined(_MSC_VER)
    return _fseeki64(f, offset, SEEK_SET) == 0;
#elif defined(__MINGW32__)
    return fseeko(f, (_off_t)offset, SEEK_SET) == 0;
#else
    return std::fseek(f, (long)offset, SEEK_SET) == 0;
#endif
  }
};

struct MemoryStream {
  MemoryStream(const char* start_, size_t size)
    : start(start_), end(start_ + size), cur(start_) {}

  char* gets(char* line, int size) {
    --size; // fgets reads in at most one less than size characters
    if (cur >= end)
      return nullptr;
    if (size > end - cur)
      size = int(end - cur);
    const char* nl = (const char*) std::memchr(cur, '\n', size);
    size_t len = nl ? nl - cur + 1 : size;
    std::memcpy(line, cur, len);
    line[len] = '\0';
    cur += len;
    return line;
  }
  int getc() { return cur < end ? *cur++ : EOF; }

  bool read(void* buf, size_t len) {
    if (cur + len > end)
      return false;
    std::memcpy(buf, cur, len);
    cur += len;
    return true;
  }

  std::string read_rest() {
    const char* last = cur;
    cur = end;
    return std::string(last, end);
  }

  int seek(std::ptrdiff_t offset) {
    cur = start + offset;
    return cur < end;
  }

private:
  const char* const start;
  const char* const end;
  const char* cur;
};

class CharArray {
  std::unique_ptr<char, decltype(&std::free)> ptr_;
  size_t size_;
public:
  CharArray() : ptr_(nullptr, &std::free), size_(0) {}
  explicit CharArray(size_t n) : ptr_((char*)std::malloc(n), &std::free), size_(n) {};
  explicit operator bool() const { return (bool)ptr_; }
  char* data() { return ptr_.get(); }
  const char* data() const { return ptr_.get(); }
  size_t size() const { return size_; }
  void set_size(size_t n) { size_ = n; }

  MemoryStream stream() const { return MemoryStream(data(), size()); }

  void resize(size_t n) {
    char* new_ptr = (char*) std::realloc(ptr_.get(), n);
    if (!new_ptr && n != 0)
      fail("Out of memory.");
    (void) ptr_.release();
    ptr_.reset(new_ptr);
    size_ = n;
  }

  // Remove first n bytes making space for more text at the returned position.
  char* roll(size_t n) {
    assert(n <= size());
    std::memmove(data(), data() + n, n);
    return data() + n;
  }
};

class BasicInput {
public:
  explicit BasicInput(const std::string& path) : path_(path) {}

  const std::string& path() const { return path_; };
  const std::string& basepath() const { return path_; };

  // Does the path stands for stdin?
  // Each reading function needs to call it (some functions use stdin
  // and some std::cin, so we don't try to unify it here).
  bool is_stdin() const { return path() == "-"; };

  // providing the same interface as MaybeGzipped
  bool is_compressed() const { return false; }
  FileStream get_uncompressing_stream() const { assert(0); unreachable(); }
  // for reading (uncompressing into memory) the whole file at once
  CharArray uncompress_into_buffer(size_t=0) { return {}; }

private:
  std::string path_;
};

template<typename Input>
inline size_t copy_line_from_stream(char* line, int size, Input&& in) {
  if (!in.gets(line, size))
    return 0;
  size_t len = std::strlen(line);
  // If a line is longer than size we discard the rest of it.
  if (len > 0 && line[len-1] != '\n')
    for (int c = in.getc(); c > 0 /* not 0 nor EOF */ && c != '\n'; c = in.getc())
      continue;
  return len;
}

} // namespace gemmi
#endif
