// Copyright 2020 Global Phasing Ltd.
//
// Read sequence from PIR or FASTA format.

#ifndef GEMMI_PIRFASTA_HPP_
#define GEMMI_PIRFASTA_HPP_

#include <cctype>   // for isspace
#include <istream>
#include "fail.hpp"
#include "util.hpp"

namespace gemmi {

inline std::string read_pir_or_fasta(std::istream& is) {
  std::string sequence;
  std::string line;
  std::getline(is, line);
  if (line[0] != '>')
    fail("PIR/FASTA files start with '>'");
  const bool pir_format = (line.size() > 3 && line[3] == ';');
  bool ensure_pir_format = false;
  bool desc_line = pir_format;
  size_t discard_n = 0;
  while (is.peek() != '>' && std::getline(is, line)) {
    for (char c : line) {
      if (('A' <= c && c <= 'Z') || ('a' <= c && c <= 'z')) {
        sequence += c;
      } else if (!std::isspace(c) && c != '-') {
        if (desc_line) {
          ensure_pir_format = true;
          sequence.clear();
          break;
        }
        if (pir_format && c == '*')
          return sequence.substr(discard_n);
        fail("unexpected character in sequence file: ", c);
      }
    }
    if (desc_line)
      discard_n = sequence.size();
    desc_line = false;
  }
  if (ensure_pir_format)
    fail("PIR sequence must end with '*'");
  return sequence;
}

} // namespace gemmi
#endif
