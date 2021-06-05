// Conversion between UTF-8 and wchar. Used only for file names on Windows.

#ifndef GEMMI_UTF_HPP_
#define GEMMI_UTF_HPP_

#include <stdexcept>  // for runtime_error
#include <string>

namespace gemmi {

// from Mark Ransom's answer
// https://stackoverflow.com/questions/148403/utf8-to-from-wide-char-conversion-in-stl/148766#148766
inline std::wstring UTF8_to_wchar(const char* in) {
  std::wstring out;
  unsigned int codepoint = 0;
  while (*in != 0) {
    unsigned char ch = static_cast<unsigned char>(*in);
    if (ch <= 0x7f)
      codepoint = ch;
    else if (ch <= 0xbf)
      codepoint = (codepoint << 6) | (ch & 0x3f);
    else if (ch <= 0xdf)
      codepoint = ch & 0x1f;
    else if (ch <= 0xef)
      codepoint = ch & 0x0f;
    else
      codepoint = ch & 0x07;
    ++in;
    if ((*in & 0xc0) != 0x80 && codepoint <= 0x10ffff) {
      if (sizeof(wchar_t) > 2) {
        out.append(1, static_cast<wchar_t>(codepoint));
      } else if (codepoint > 0xffff) {
        out.append(1, static_cast<wchar_t>(0xd800 + (codepoint >> 10)));
        out.append(1, static_cast<wchar_t>(0xdc00 + (codepoint & 0x03ff)));
      } else if (codepoint < 0xd800 || codepoint >= 0xe000) {
        out.append(1, static_cast<wchar_t>(codepoint));
      }
    }
  }
  return out;
}

inline std::string wchar_to_UTF8(const wchar_t* in) {
  std::string out;
  unsigned int codepoint = 0;
  while (*in != 0) {
    if (*in >= 0xd800 && *in <= 0xdbff) {
      codepoint = ((*in - 0xd800) << 10) + 0x10000;
    } else {
      if (*in >= 0xdc00 && *in <= 0xdfff)
        codepoint |= *in - 0xdc00;
      else
        codepoint = *in;
      if (codepoint <= 0x7f) {
        out += static_cast<char>(codepoint);
      } else if (codepoint <= 0x7ff) {
        out += static_cast<char>(0xc0 | ((codepoint >> 6) & 0x1f));
        out += static_cast<char>(0x80 | (codepoint & 0x3f));
      } else if (codepoint <= 0xffff) {
        out += static_cast<char>(0xe0 | ((codepoint >> 12) & 0x0f));
        out += static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f));
        out += static_cast<char>(0x80 | (codepoint & 0x3f));
      } else {
        out += static_cast<char>(0xf0 | ((codepoint >> 18) & 0x07));
        out += static_cast<char>(0x80 | ((codepoint >> 12) & 0x3f));
        out += static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f));
        out += static_cast<char>(0x80 | (codepoint & 0x3f));
      }
      codepoint = 0;
    }
    ++in;
  }
  return out;
}

} // namespace gemmi
#endif
