// Copyright 2017 Global Phasing Ltd.

// Writing cif::Document or its parts to std::ostream.

#ifndef GEMMI_TO_CIF_HPP_
#define GEMMI_TO_CIF_HPP_

#include <ostream>
#include "cifdoc.hpp"

namespace gemmi {
namespace cif {

enum class Style {
  Simple,
  NoBlankLines,
  PreferPairs,  // write single-row loops as pairs
  Pdbx,         // PreferPairs + put '#' (empty comments) between categories
  Indent35,     // start values in pairs from 35th column
};

// CIF files are read in binary mode. It makes difference only for text fields.
// If the text field with \r\n would be written as is in text mode on Windows
// \r would get duplicated. As a workaround, here we convert \r\n to \n.
// Hopefully \r that gets removed here is never meaningful.
inline void write_text_field(std::ostream& os, const std::string& value) {
  for (size_t pos = 0, end = 0; end != std::string::npos; pos = end + 1) {
    end = value.find("\r\n", pos);
    size_t len = (end == std::string::npos ? value.size() : end) - pos;
    os.write(value.c_str() + pos, len);
  }
}

inline void write_out_pair(std::ostream& os, const std::string& name,
                           const std::string& value, Style style) {
  os << name;
  if (is_text_field(value)) {
    os.put('\n');
    write_text_field(os, value);
  } else {
    if (name.size() + value.size() > 120)
      os.put('\n');
    else if (style == Style::Indent35 && name.size() < 34)
      os.write("                                  ", 34 - name.size());
    else
      os.put(' ');
    os << value;
  }
  os.put('\n');
}

inline void write_out_loop_values(std::ostream& os, const Loop& loop) {
  size_t ncol = loop.tags.size();
  size_t col = 0;
  for (const std::string& val : loop.values) {
    bool text_field = is_text_field(val);
    os.put(col++ == 0 || text_field ? '\n' : ' ');
    if (text_field)
      write_text_field(os, val);
    else
      os << val;
    if (col == ncol)
      col = 0;
  }
}

inline void write_out_loop(std::ostream& os, const Loop& loop, Style style) {
  if (loop.values.empty())
    return;
  if ((style == Style::PreferPairs || style == Style::Pdbx) &&
      loop.length() == 1) {
    for (size_t i = 0; i != loop.tags.size(); ++i)
      write_out_pair(os, loop.tags[i], loop.values[i], style);
    return;
  }
  os << "loop_";
  for (const std::string& tag : loop.tags)
    os << '\n' << tag;
  write_out_loop_values(os, loop);
  os.put('\n');
}

inline void write_out_item(std::ostream& os, const Item& item, Style style) {
  switch (item.type) {
    case ItemType::Pair:
      write_out_pair(os, item.pair[0], item.pair[1], style);
      break;
    case ItemType::Loop:
      write_out_loop(os, item.loop, style);
      break;
    case ItemType::Frame:
      os << "save_" << item.frame.name << '\n';
      for (const Item& inner_item : item.frame.items)
        write_out_item(os, inner_item, style);
      os << "save_\n";
      break;
    case ItemType::Comment:
      os << item.pair[1] << '\n';
      break;
    case ItemType::Erased:
      break;
  }
}

inline bool should_be_separated_(const Item& a, const Item& b) {
  if (a.type == ItemType::Comment || b.type == ItemType::Comment)
    return false;
  if (a.type != ItemType::Pair || b.type != ItemType::Pair)
    return true;
  // check if we have mmcif-like tags from different categories
  auto adot = a.pair[0].find('.');
  if (adot == std::string::npos)
    return false;
  auto bdot = b.pair[0].find('.');
  return adot != bdot || a.pair[0].compare(0, adot, b.pair[0], 0, adot) != 0;
}

inline void write_cif_block_to_stream(std::ostream& os, const Block& block,
                                      Style style=Style::Simple) {
    os << "data_" << block.name << '\n';
    if (style == Style::Pdbx)
      os << "#\n";
    const Item* prev = nullptr;
    for (const Item& item : block.items)
      if (item.type != ItemType::Erased) {
        if (prev && style != Style::NoBlankLines &&
            should_be_separated_(*prev, item))
          os << (style == Style::Pdbx ? "#\n" : "\n");
        write_out_item(os, item, style);
        prev = &item;
      }
    if (style == Style::Pdbx)
      os << "#\n";
}

inline void write_cif_to_stream(std::ostream& os, const Document& doc,
                                Style style=Style::Simple) {
  bool first = true;
  for (const Block& block : doc.blocks) {
    if (!first)
      os.put('\n'); // extra blank line for readability
    write_cif_block_to_stream(os, block, style);
    first = false;
  }
}

} // namespace cif
} // namespace gemmi

#endif
