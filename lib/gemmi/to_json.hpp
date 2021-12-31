// Copyright 2017 Global Phasing Ltd.

// Writing cif::Document or its parts as JSON (mmJSON, CIF-JSON, etc).

#ifndef GEMMI_TO_JSON_HPP_
#define GEMMI_TO_JSON_HPP_
#include <cctype>    // for isdigit
#include <ostream>   // for ostream
#include <set>       // for set
#include <string>    // for string
#include <vector>    // for vector
#include "cifdoc.hpp"
#include "numb.hpp"  // for is_numb
#include "util.hpp"  // for starts_with

namespace gemmi {
namespace cif {

class JsonWriter {
public:
  bool comcifs = false;  // conform to the COMCIFS CIF-JSON draft
  bool group_ddl2_categories = false;  // for mmJSON
  bool with_data_keyword = false;  // for mmJSON
  bool bare_tags = false;  // "tag" instead of "_tag"
  bool values_as_arrays = false;  // "_tag": ["value"]
  bool lowercase_names = true; // write case-insensitive names as lower case
  int quote_numbers = 1;  // 0=never (no s.u.), 1=mix, 2=always
  std::string cif_dot = "null";  // how to convert '.' from CIF
  explicit JsonWriter(std::ostream& os) : os_(os), linesep_("\n ") {}
  void write_json(const Document& d);
  void set_comcifs() {
    comcifs = true;
    values_as_arrays = true;
    quote_numbers = 2;
    cif_dot = "false";
  }
  void set_mmjson() {
    group_ddl2_categories = true;
    with_data_keyword = true;
    bare_tags = true;
    values_as_arrays = true;
    lowercase_names = false;
    quote_numbers = 0;
  }

private:
  std::ostream& os_;
  std::string linesep_;

  void change_indent(int n) { linesep_.resize(linesep_.size() + n, ' '); }

  // returns category with trailing dot
  std::string get_tag_category(const std::string& tag) const {
    if (!group_ddl2_categories)
      return std::string{};
    size_t pos = tag.find('.');
    if (pos == std::string::npos)
      return std::string{};
    return tag.substr(0, pos + 1);
  }

  std::string get_loop_category(const Loop& loop) const {
    if (loop.tags.empty())
      return std::string{};
    std::string cat = get_tag_category(loop.tags[0]);
    for (size_t i = 1; i < loop.tags.size(); ++i)
      if (!starts_with(loop.tags[i], cat))
        return std::string{};
    return cat;
  }

  // based on tao/json/internal/escape.hpp
  static void escape(std::ostream& os, const std::string& s, size_t pos,
                     bool to_lower) {
    static const char* h = "0123456789abcdef";
    const char* p = s.data() + pos;
    const char* l = p;
    const char* const e = s.data() + s.size();
    while (p != e) {
      const unsigned char c = *p;
      if (c == '\\') {
        os.write(l, p - l);
        l = ++p;
        os << "\\\\";
      } else if (c == '"') {
        os.write(l, p - l);
        l = ++p;
        os << "\\\"";
      } else if (c < 32) {
        os.write(l, p - l);
        l = ++p;
        switch ( c ) {
          case '\b': os << "\\b"; break;
          case '\f': os << "\\f"; break;
          case '\n': os << "\\n"; break;
          case '\r': os << "\\r"; break;
          case '\t': os << "\\t"; break;
          default: os << "\\u00" << h[(c & 0xf0) >> 4] << h[c & 0x0f];
        }
      } else if (to_lower && c >= 'A' && c <= 'Z') {
        os.write(l, p - l);
        l = ++p;
        os.put(c + 32);
      } else if (c == 127) {
        os.write(l, p - l);
        l = ++p;
        os << "\\u007f";
      } else {
        ++p;
      }
    }
    os.write(l, p - l);
  }

  void write_string(const std::string& s, size_t pos=0, bool to_lower=false) {
    os_.put('"');
    escape(os_, s, pos, to_lower);
    os_.put('"');
  }

  void write_as_number(const std::string& value) {
    // if we are here, value is not empty
    if (value[0] == '.') // in JSON numbers cannot start with dot
      os_.put('0');
    // in JSON the number cannot start with +
    size_t pos = 0;
    if (value[pos] == '+') {
      pos = 1;
    } else if (value[pos] == '-') { // make handling -001 easier
      os_.put('-');
      pos = 1;
    }
    // in JSON left-padding with 0s is not allowed
    while (value[pos] == '0' && std::isdigit(value[pos+1]))
      ++pos;
    // in JSON dot must be followed by digit
    size_t dotpos = value.find('.');
    if (dotpos != std::string::npos && !std::isdigit(value[dotpos+1])) {
      os_ << value.substr(pos, dotpos+1-pos) << '0';
      pos = dotpos + 1;
    }
    if (value.back() != ')')
      os_ << value.c_str() + pos;
    else
      os_ << value.substr(pos, value.find('(', pos) - pos);
  }

  void write_value(const std::string& value) {
    if (value == "?")
      os_ << "null";
    else if (value == ".")
      os_ << cif_dot;
    else if (quote_numbers < 2 && is_numb(value) &&
             // exception: 012 (but not 0.12) is assumed to be a string
             (value[0] != '0' || value[1] == '.' || value[1] == '\0') &&
             (quote_numbers == 0 || value.back() != ')'))
      write_as_number(value);
    else
      write_string(as_string(value));
  }

  void open_cat(const std::string& cat, size_t* tag_pos) {
    if (!cat.empty()) {
      change_indent(+1);
      write_string(cat.substr(0, cat.size() - 1), bare_tags ? 1 : 0, lowercase_names);
      os_ << ": {" << linesep_;
      *tag_pos += cat.size() - 1;
    }
  }

  void close_cat(std::string& cat, size_t* tag_pos) {
    if (!cat.empty()) {
      change_indent(-1);
      os_ << linesep_ << '}';
      *tag_pos -= cat.size() - 1;
      cat.clear();
    }
  }

  void write_loop(const Loop& loop) {
    size_t ncol = loop.tags.size();
    const auto& vals = loop.values;
    std::string cat = get_loop_category(loop);
    size_t tag_pos = bare_tags ? 1 : 0;
    open_cat(cat, &tag_pos);
    for (size_t i = 0; i < ncol; i++) {
      if (i != 0)
        os_ << "," << linesep_;
      write_string(loop.tags[i], tag_pos, lowercase_names);
      os_ << ": [";
      for (size_t j = i; j < vals.size(); j += ncol) {
        if (j != i)
          os_.put(',');
        write_value(vals[j]);
      }
      os_.put(']');
    }
    close_cat(cat, &tag_pos);
  }


  // works for both block and frame
  void write_map(const std::string& name, const std::vector<Item>& items) {
    write_string(name, 0, lowercase_names);
    os_ << ": ";
    change_indent(+1);
    char first = '{';
    bool has_frames = false;
    std::string cat;
    size_t tag_pos = bare_tags ? 1 : 0;
    // When grouping into categories, only consecutive tags are grouped.
    std::set<std::string> seen_cats;
    for (const Item& item : items) {
      switch (item.type) {
        case ItemType::Pair:
          if (!cat.empty() && !starts_with(item.pair[0], cat))
            close_cat(cat, &tag_pos);
          os_ << first << linesep_;
          if (group_ddl2_categories && cat.empty()) {
            cat = get_tag_category(item.pair[0]);
            if (seen_cats.insert(cat).second)
              open_cat(cat, &tag_pos);
          }
          write_string(item.pair[0], tag_pos, lowercase_names);
          os_ << ": ";
          if (values_as_arrays)
            os_.put('[');
          write_value(item.pair[1]);
          if (values_as_arrays)
            os_.put(']');
          first = ',';
          break;
        case ItemType::Loop:
          close_cat(cat, &tag_pos);
          os_ << first << linesep_;
          write_loop(item.loop);
          first = ',';
          break;
        case ItemType::Frame:
          has_frames = true;
          break;
        case ItemType::Comment:
          break;
        case ItemType::Erased:
          break;
      }
    }
    if (has_frames) {  // usually, we don't have any frames
      os_ << first << linesep_ << "\"Frames\": ";
      change_indent(+1);
      first = '{';
      for (const Item& item : items)
        if (item.type == ItemType::Frame) {
          os_ << first << linesep_;
          write_map(item.frame.name, item.frame.items);
          first = ',';
        }
      change_indent(-1);
      os_ << linesep_ << '}';
    }
    close_cat(cat, &tag_pos);
    change_indent(-1);
    os_ << linesep_ << '}';
  }
};

inline void JsonWriter::write_json(const Document& d) {
  os_.put('{');
  if (comcifs) {
    os_ << R"(
 "CIF-JSON": {
  "Metadata": {
   "cif-version": "2.0",
   "schema-name": "CIF-JSON",
   "schema-version": "1.0.0",
   "schema-uri": "http://www.iucr.org/resources/cif/cif-json.json"
  },)";
    change_indent(+1);
  }
  for (const Block& block : d.blocks) {
    if (&block != &d.blocks[0])
      os_.put(',');
    // start mmJSON with {"data_ so it can be easily recognized
    if (&block != &d.blocks[0] || comcifs || !with_data_keyword)
      os_ << linesep_;
    write_map((with_data_keyword ? "data_" : "") + block.name, block.items);
  }
  if (comcifs)
    os_ << "\n }";
  os_ << "\n}\n";
}

inline void write_mmjson_to_stream(std::ostream& os, const Document& doc) {
  cif::JsonWriter writer(os);
  writer.set_mmjson();
  writer.write_json(doc);
}

} // namespace cif
} // namespace gemmi
#endif
