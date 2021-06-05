// Copyright 2017 Global Phasing Ltd.
//
// struct Document that represents the CIF file (but can be also
// read from JSON file, such as CIF-JSON or mmJSON).

#ifndef GEMMI_CIFDOC_HPP_
#define GEMMI_CIFDOC_HPP_
#include "iterator.hpp"  // for StrideIter, IndirectIter
#include "atox.hpp"  // for string_to_int
#include "fail.hpp"  // for fail
#include "util.hpp"  // for starts_with, to_lower
#include "tostr.hpp"  // for tostr
#include <cassert>
#include <cstring>   // for memchr
#include <algorithm> // for move, find_if, all_of, min, rotate
#include <array>
#include <initializer_list>
#include <iosfwd>    // for size_t, ptrdiff_t
#include <new>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#if defined(_MSC_VER)
#pragma warning(push)
// warning C4244: an integer type is converted to a smaller integer type
#pragma warning(disable: 4244)
// warning C4267: conversion from 'size_t' to 'type', possible loss of data
#pragma warning(disable: 4267)
#endif

namespace gemmi {
namespace cif {
using std::size_t;
using gemmi::fail;

enum class ItemType : unsigned char {
  Pair,
  Loop,
  Frame,
  Comment,
  Erased,
};

inline uint8_t char_table(char c) {
  static const uint8_t table[256] = {
   // 0  1  2  3  4  5  6  7  8  9  A  B  C  D  E  F
      0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 2, 0, 0, // 0
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 1
      2, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, // 2
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, // 3
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, // 4
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, // 5
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, // 6
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, // 7
   // 128-255
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  };
  return table[static_cast<unsigned char>(c)];
}

inline void assert_tag(const std::string& tag) {
  if (tag[0] != '_')
    fail("Tag should start with '_', got: " + tag);
}

inline void ensure_mmcif_category(std::string& cat) {
  if (cat[0] != '_')
    fail("Category should start with '_', got: " + cat);
  if (*(cat.end() - 1) != '.')
    cat += '.';
}

inline bool is_null(const std::string& value) {
  return value.size() == 1 && (value[0] == '?' || value[0] == '.');
}

inline std::string as_string(const std::string& value) {
  if (value.empty() || is_null(value))
    return "";
  if (value[0] == '"' || value[0] == '\'')
    return std::string(value.begin() + 1, value.end() - 1);
  if (value[0] == ';' && value.size() > 2 && *(value.end() - 2) == '\n') {
    bool crlf = *(value.end() - 3) == '\r';
    return std::string(value.begin() + 1, value.end() - (crlf ? 3 : 2));
  }
  return value;
}

inline std::string as_string(const std::string* value) {
  return value ? as_string(*value) : std::string();
}

inline char as_char(const std::string& value, char null) {
  if (is_null(value))
    return null;
  if (value.size() < 2)
    return value[0];
  const std::string s = as_string(value);
  if (s.size() < 2)
    return s[0];
  fail("Not a single character: " + value);
}

inline int as_int(const std::string& str) {
  return string_to_int(str, true);
}

inline int as_int(const std::string& str, int null) {
  return is_null(str) ? null : as_int(str);
}

// for use in templates (see also as_any() functions in numb.hpp)
inline int as_any(const std::string& s, int null) { return as_int(s, null); }
inline char as_any(const std::string& s, char null) { return as_char(s, null); }


using Pair = std::array<std::string, 2>;

// used only as arguments when creating Item
struct LoopArg {};
struct FrameArg { std::string str; };
struct CommentArg { std::string str; };

struct Loop {
  std::vector<std::string> tags;
  std::vector<std::string> values;

  // search and access
  int find_tag(std::string tag) const {
    tag = gemmi::to_lower(tag);
    auto f = std::find_if(tags.begin(), tags.end(),
               [&tag](const std::string& t) { return gemmi::iequal(t, tag); });
    return f == tags.end() ? -1 : f - tags.begin();
  }
  bool has_tag(const std::string& tag) const { return find_tag(tag) != -1; }
  size_t width() const { return tags.size(); }
  size_t length() const { return values.size() / tags.size(); }
  const std::string& val(size_t row, size_t col) const {
    return values[row * tags.size() + col];
  }

  void clear() { tags.clear(); values.clear(); }

  template <typename T> void add_row(T new_values, int pos=-1) {
    if (new_values.size() != tags.size())
      fail("add_row(): wrong row length.");
    auto it = values.end();
    if (pos >= 0 && pos * width() < values.size())
      it = values.begin() + pos * tags.size();
    values.insert(it, new_values.begin(), new_values.end());
  }
  void add_row(std::initializer_list<std::string> new_values, int pos=-1) {
    add_row<std::initializer_list<std::string>>(new_values, pos);
  }
  // comments are added relying on how cif writing works
  void add_comment_and_row(std::initializer_list<std::string> ss) {
    if (ss.size() != tags.size() + 1)
      fail("add_comment_and_row(): wrong row length.");
    std::vector<std::string> vec(ss.begin() + 1, ss.end());
    vec[0] = "#" + *ss.begin() + "\n" + vec[0];
    return add_row(vec);
  }
  void pop_row() {
    if (values.size() < tags.size())
      fail("pop_row() called on empty Loop");
    values.resize(values.size() - tags.size());
  }

  void set_all_values(std::vector<std::vector<std::string>> columns);
};


struct Item;
struct Block;

// Accessor to a specific loop column, or to a single value from a Pair.
class Column {
public:
  Column() : item_(nullptr) {}
  Column(Item* item, size_t col) : item_(item), col_(col) {}
  using iterator = StrideIter<std::string>;
  iterator begin();
  iterator end();
  using const_iterator = StrideIter<const std::string>;
  const_iterator begin() const { return const_cast<Column*>(this)->begin(); }
  const_iterator end() const { return const_cast<Column*>(this)->end(); }

  Loop* get_loop() const;
  std::string* get_tag();
  const std::string* get_tag() const {
    return const_cast<Column*>(this)->get_tag();
  }
  int length() const {
    if (const Loop* loop = get_loop())
      return loop->length();
    return item_ ? 1 : 0;
  }
  explicit operator bool() const { return item_ != nullptr ; }
  std::string& operator[](int n);
  std::string& at(int n) {
    if (n < 0)
      n += length();
    if (n < 0 || n >= length())
      throw std::out_of_range("Cannot access element " + std::to_string(n) +
          " in Column with length " + std::to_string(length()));
    return operator[](n);
  }
  const std::string& at(int n) const {
    return const_cast<Column*>(this)->at(n);
  }

  std::string str(int n) const { return as_string(at(n)); }
  const Item* item() const { return item_; }
  Item* item() { return item_; }
  size_t col() const { return col_; }

private:
  Item* item_;
  size_t col_;  // for loop this is a column index in item_->loop
};

// Some values can be given either in loop or as tag-value pairs.
// The latter case is equivalent to a loop with a single row.
// We optimized for loops, and in case of tag-values we copy the values
// into the `values` vector.
struct Table {
  Item* loop_item;
  Block& bloc;
  std::vector<int> positions;
  size_t prefix_length;

  struct Row {
    Table& tab;
    int row_index;

    std::string& value_at_unsafe(int pos);
    std::string& value_at(int pos) {
      if (pos == -1)
        throw std::out_of_range("Cannot access missing optional tag.");
      return value_at_unsafe(pos);
    }
    const std::string& value_at(int pos) const {
      return const_cast<Row*>(this)->value_at(pos);
    }

    std::string& at(int n) {
      return value_at(tab.positions.at(n < 0 ? n + size() : n));
    }
    const std::string& at(int n) const { return const_cast<Row*>(this)->at(n); }

    std::string& operator[](int n);
    const std::string& operator[](int n) const {
      return const_cast<Row*>(this)->operator[](n);
    }

    std::string* ptr_at(int n) {
      int pos = tab.positions.at(n < 0 ? n + size() : n);
      return pos >= 0 ? &value_at(pos) : nullptr;
    }
    const std::string* ptr_at(int n) const {
      return const_cast<Row*>(this)->ptr_at(n);
    }

    bool has(int n) const { return tab.positions.at(n) >= 0; }
    bool has2(int n) const { return has(n) && !cif::is_null(operator[](n)); }

    const std::string& one_of(int n1, int n2) const {
      static const std::string nul(1, '.');
      if (has2(n1))
       return operator[](n1);
      if (has(n2))
       return operator[](n2);
      return nul;
    }

    size_t size() const { return tab.width(); }

    std::string str(int n) const { return as_string(at(n)); }

    using iterator = IndirectIter<Row, std::string>;
    using const_iterator = IndirectIter<const Row, const std::string>;
    iterator begin() { return iterator({this, tab.positions.begin()}); }
    iterator end() { return iterator({this, tab.positions.end()}); }
    const_iterator begin() const {
      return const_iterator({this, tab.positions.begin()});
    }
    const_iterator end() const {
      return const_iterator({this, tab.positions.end()});
    }
  };

  Loop* get_loop();
  bool ok() const { return !positions.empty(); }
  size_t width() const { return positions.size(); }
  size_t length() const;
  size_t size() const { return length(); }
  bool has_column(int n) const { return ok() && positions.at(n) >= 0; }
  int first_of(int n1, int n2) const { return positions.at(n1) >= 0 ? n1 : n2; }
  Row tags() { return Row{*this, -1}; }
  Row operator[](int n) { return Row{*this, n}; }

  Row at(int n) {
    if (n < 0)
      n += length();
    if (n < 0 || static_cast<size_t>(n) >= length())
      throw std::out_of_range("No row with index " + std::to_string(n));
    return (*this)[n];
  }

  Row one() {
    if (length() != 1)
      fail("Expected one value, found " + std::to_string(length()));
    return (*this)[0];
  }

  std::string get_prefix() const {
    for (int pos : positions)
      if (pos >= 0)
        return const_cast<Table*>(this)->tags()
               .value_at(pos).substr(0, prefix_length);
    fail("The table has no columns.");
  }

  Row find_row(const std::string& s);

  template <typename T> void append_row(T new_values);
  void append_row(std::initializer_list<std::string> new_values) {
    append_row<std::initializer_list<std::string>>(new_values);
  }
  void remove_row(int row_index) { remove_rows(row_index, row_index+1); }
  void remove_rows(int start, int end);
  Column column_at_pos(int pos);
  Column column(int n) {
    int pos = positions.at(n);
    if (pos == -1)
      fail("Cannot access absent column");
    return column_at_pos(pos);
  }

  // prefix is optional
  int find_column_position(const std::string& tag) const {
    Row tag_row = const_cast<Table*>(this)->tags();
    for (int pos : positions) {
      const std::string& v = tag_row.value_at_unsafe(pos);
      if (v == tag || v.compare(prefix_length, std::string::npos, tag) == 0)
        return pos;
    }
    fail("Column name not found: " + tag);
  }

  Column find_column(const std::string& tag) {
    return column_at_pos(find_column_position(tag));
  }

  void erase();

  void convert_pair_to_loop();

  // It is not a proper input iterator, but just enough for using range-for.
  struct iterator {
    Table& parent;
    int index;
    void operator++() { index++; }
    bool operator==(const iterator& o) const { return index == o.index; }
    bool operator!=(const iterator& o) const { return index != o.index; }
    Row operator*() { return parent[index]; }
    const std::string& get(int n) const { return parent[index].at(n); }
  };
  iterator begin() { return iterator{*this, 0}; }
  iterator end() { return iterator{*this, (int)length()}; }
};

struct Block {
  std::string name;
  std::vector<Item> items;

  explicit Block(const std::string& name_) : name(name_) {}
  Block() {}

  void swap(Block& o) { name.swap(o.name); items.swap(o.items); }
  // access functions
  const Item* find_pair_item(const std::string& tag) const;
  const Pair* find_pair(const std::string& tag) const;
  const std::string* find_value(const std::string& tag) const {
    const Pair* pair = find_pair(tag);
    return pair ? &(*pair)[1] : nullptr;
  }
  Column find_loop(const std::string& tag);
  const Item* find_loop_item(const std::string& tag) const;
  Column find_values(const std::string& tag);
  bool has_tag(const std::string& tag) const {
    return const_cast<Block*>(this)->find_values(tag).item() != nullptr;
  }
  bool has_any_value(const std::string& tag) const {
    Column c = const_cast<Block*>(this)->find_values(tag);
    return c.item() != nullptr && !std::all_of(c.begin(), c.end(), is_null);
  }
  Table find(const std::string& prefix,
             const std::vector<std::string>& tags);
  Table find(const std::vector<std::string>& tags) { return find({}, tags); }
  Table find_any(const std::string& prefix,
                 const std::vector<std::string>& tags);
  Table find_or_add(const std::string& prefix, std::vector<std::string> tags) {
    Table t = find(prefix, tags);
    if (!t.ok()) {
      for (int i = 0; i != (int) tags.size(); ++i)
        t.positions.push_back(i);
      t.loop_item = &setup_loop_item(find_any(prefix, tags), prefix,
                                     std::move(tags));
    }
    return t;
  }
  Block* find_frame(std::string name);
  Table item_as_table(Item& item);

  size_t get_index(const std::string& tag) const;

  // modifying functions
  void set_pair(const std::string& tag, const std::string& value);

  Loop& init_loop(const std::string& prefix, std::vector<std::string> tags) {
    return setup_loop(find_any(prefix, tags), prefix, std::move(tags));
  }

  void move_item(int old_pos, int new_pos);

  // mmCIF specific functions
  std::vector<std::string> get_mmcif_category_names() const;
  Table find_mmcif_category(std::string cat);

  Loop& init_mmcif_loop(std::string cat, std::vector<std::string> tags) {
    ensure_mmcif_category(cat);
    return setup_loop(find_mmcif_category(cat), cat, std::move(tags));
  }

private:
  Item& setup_loop_item(Table&& tab, const std::string& prefix,
                        std::vector<std::string>&& tags);
  Loop& setup_loop(Table&& tab, const std::string& prefix,
                   std::vector<std::string>&& tags);
};

struct Item {
  ItemType type;
  int line_number = -1;
  union {
    Pair pair;
    Loop loop;
    Block frame;
  };

  explicit Item(LoopArg)
    : type{ItemType::Loop}, loop{} {}
  explicit Item(std::string&& t)
    : type{ItemType::Pair}, pair{{std::move(t), std::string()}} {}
  Item(const std::string& t, const std::string& v)
    : type{ItemType::Pair}, pair{{t, v}} {}
  explicit Item(FrameArg&& frame_arg)
    : type{ItemType::Frame}, frame(frame_arg.str) {}
  explicit Item(CommentArg&& comment)
    : type{ItemType::Comment}, pair{{std::string(), std::move(comment.str)}} {}

  Item(Item&& o) noexcept
      : type(o.type), line_number(o.line_number) {
    move_value(std::move(o));
  }
  Item(const Item& o)
      : type(o.type), line_number(o.line_number) {
    copy_value(o);
  }

  Item& operator=(Item o) { set_value(std::move(o)); return *this; }

  ~Item() { destruct(); }

  void erase() {
    destruct();
    type = ItemType::Erased;
  }

  bool has_prefix(const std::string& prefix) const {
    return (type == ItemType::Pair && gemmi::starts_with(pair[0], prefix)) ||
           (type == ItemType::Loop && !loop.tags.empty() &&
            gemmi::starts_with(loop.tags[0], prefix));
  }

  void set_value(Item&& o) {
    if (type == o.type) {
      switch (type) {
        case ItemType::Pair: pair = std::move(o.pair); break;
        case ItemType::Loop: loop = std::move(o.loop); break;
        case ItemType::Frame: frame = std::move(o.frame); break;
        case ItemType::Comment: pair = std::move(o.pair); break;
        case ItemType::Erased: break;
      }
    } else {
      destruct();
      type = o.type;
      move_value(std::move(o));
    }
  }

private:
  void destruct() {
    switch (type) {
      case ItemType::Pair: pair.~Pair(); break;
      case ItemType::Loop: loop.~Loop(); break;
      case ItemType::Frame: frame.~Block(); break;
      case ItemType::Comment: pair.~Pair(); break;
      case ItemType::Erased: break;
    }
  }

  void copy_value(const Item& o) {
    if (o.type == ItemType::Pair || o.type == ItemType::Comment)
      new (&pair) Pair(o.pair);
    else if (o.type == ItemType::Loop)
      new (&loop) Loop(o.loop);
    else if (o.type == ItemType::Frame)
      new (&frame) Block(o.frame);
  }

  void move_value(Item&& o) {
    if (o.type == ItemType::Pair || o.type == ItemType::Comment)
      new (&pair) Pair(std::move(o.pair));
    else if (o.type == ItemType::Loop)
      new (&loop) Loop(std::move(o.loop));
    else if (o.type == ItemType::Frame)
      new (&frame) Block(std::move(o.frame));
  }
};

inline void Loop::set_all_values(std::vector<std::vector<std::string>> columns){
  size_t w = columns.size();
  if (w != width())
    fail(tostr("set_all_values(): expected ", width(), " columns, got ", w));
  if (w == 0)
    return;
  size_t h = columns[0].size();
  for (auto& col : columns)
    if (col.size() != h)
      fail("set_all_values(): all columns must have the same length");
  values.resize(w * h);
  for (size_t i = 0; i != h; ++i)
    for (size_t j = 0; j != w; ++j)
      values[w * i + j] = std::move(columns[j][i]);
}

inline std::string* Column::get_tag() {
  if (!item_)
    return nullptr;
  if (Loop* loop = get_loop())
    return &loop->tags.at(col_);
  return &item_->pair[0];
}

inline Loop* Column::get_loop() const {
  return item_ && item_->type == ItemType::Loop ? &item_->loop : nullptr;
}
inline Column::iterator Column::begin() {
  if (Loop* loop = get_loop())
    return iterator({loop->values.data(), col_, (unsigned) loop->width()});
  if (item_ && item_->type == ItemType::Pair)
    return iterator({&item_->pair[1], 0, 1});
  return iterator();
}

inline Column::iterator Column::end() {
  if (Loop* loop = get_loop())
    return iterator({loop->values.data() + loop->values.size(),
                    col_, (unsigned) loop->width()});
  if (item_ && item_->type == ItemType::Pair)
    return iterator({&item_->pair[1] + 1, 0, 1});
  return iterator();
}

inline std::string& Column::operator[](int n) {
  if (Loop* loop = get_loop())
    return loop->values[n * loop->width() + col_];
  return item_->pair[1];
}

inline std::string& Table::Row::operator[](int n) {
  int pos = tab.positions[n];
  if (Loop* loop = tab.get_loop()) {
    if (row_index == -1) // tags
      return loop->tags[pos];
    return loop->values[loop->width() * row_index + pos];
  }
  return tab.bloc.items[pos].pair[row_index == -1 ? 0 : 1];
}

inline std::string& Table::Row::value_at_unsafe(int pos) {
  Loop* loop = tab.get_loop();
  if (row_index == -1) { // tags
    if (loop)
      return loop->tags.at(pos);
    return tab.bloc.items[pos].pair[0];
  }
  if (loop)
    return loop->values.at(loop->width() * row_index + pos);
  return tab.bloc.items[pos].pair[1];
}

inline Loop* Table::get_loop() {
  return loop_item ? &loop_item->loop : nullptr;
}

inline size_t Table::length() const {
  return loop_item ? loop_item->loop.length() : (positions.empty() ? 0 : 1);
}

inline Table::Row Table::find_row(const std::string& s) {
  int pos = positions.at(0);
  if (const Loop* loop = get_loop()) {
    for (size_t i = 0; i < loop->values.size(); i += loop->width())
      if (as_string(loop->values[i + pos]) == s)
        return Row{*this, static_cast<int>(i / loop->width())};
  } else if (as_string(bloc.items[pos].pair[1]) == s) {
    return Row{*this, 0};
  }
  fail("Not found in " + *column_at_pos(pos).get_tag() + ": " + s);
}

template <typename T> void Table::append_row(T new_values) {
  if (!ok())
    fail("append_row(): table not found");
  if (new_values.size() != width())
    fail("append_row(): wrong row length");
  if (!loop_item)
    convert_pair_to_loop();
  Loop& loop = loop_item->loop;
  size_t cur_size = loop.values.size();
  loop.values.resize(cur_size + loop.width(), ".");
  int n = 0;
  for (const auto& value : new_values)
    loop.values[cur_size + positions[n++]] = value;
}

inline void Table::remove_rows(int start, int end) {
  if (!ok())
    // this function is used mostly through remove_row()
    fail("remove_row(): table not found");
  if (!loop_item)
    convert_pair_to_loop();
  Loop& loop = loop_item->loop;
  size_t start_pos = start * loop.width();
  size_t end_pos = end * loop.width();
  if (start_pos >= end_pos || end_pos > loop.values.size())
    throw std::out_of_range("remove_row(): invalid index");
  loop.values.erase(loop.values.begin() + start_pos,
                    loop.values.begin() + end_pos);
}

inline Column Table::column_at_pos(int pos) {
  if (loop_item)
    return Column(loop_item, pos);
  return Column(&bloc.items[pos], 0);
}

inline void Table::erase() {
  if (loop_item)
    loop_item->erase();
  else
    for (int pos : positions)
      bloc.items[pos].erase();
}

inline void Table::convert_pair_to_loop() {
  assert(loop_item == nullptr);
  Item new_item(LoopArg{});
  new_item.loop.tags.resize(positions.size());
  new_item.loop.values.resize(positions.size());
  for (size_t i = 0; i != positions.size(); ++i) {
    Item& item = bloc.items[positions[i]];
    new_item.loop.tags[i].swap(item.pair[0]);
    new_item.loop.values[i].swap(item.pair[1]);
    item.erase();
  }
  loop_item = &bloc.items.at(positions[0]);
  loop_item->set_value(std::move(new_item));
}

inline const Item* Block::find_pair_item(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Pair && i.pair[0] == tag)
      return &i;
  return nullptr;
}

inline const Pair* Block::find_pair(const std::string& tag) const {
  const Item* item = find_pair_item(tag);
  return item ? &item->pair : nullptr;
}

inline void Block::set_pair(const std::string& tag, const std::string& value) {
  assert_tag(tag);
  for (Item& i : items) {
    if (i.type == ItemType::Pair && i.pair[0] == tag) {
      i.pair[1] = value;
      return;
    }
    if (i.type == ItemType::Loop && i.loop.find_tag(tag) != -1) {
      i.set_value(Item(tag, value));
      return;
    }
  }
  items.emplace_back(tag, value);
}

inline Column Block::find_loop(const std::string& tag) {
  Column c = find_values(tag);
  return c.item() && c.item()->type == ItemType::Loop ? c : Column();
}

inline const Item* Block::find_loop_item(const std::string& tag) const {
  for (const Item& i : items)
    if (i.type == ItemType::Loop && i.loop.find_tag(tag) != -1)
      return &i;
  return nullptr;
}

inline Column Block::find_values(const std::string& tag) {
  for (Item& i : items)
    if (i.type == ItemType::Loop) {
      int pos = i.loop.find_tag(tag);
      if (pos != -1)
        return Column{&i, static_cast<size_t>(pos)};
    } else if (i.type == ItemType::Pair) {
      if (i.pair[0] == tag)
        return Column{&i, 0};
    }
  return Column{nullptr, 0};
}

inline Block* Block::find_frame(std::string frame_name) {
  frame_name = gemmi::to_lower(frame_name);
  for (Item& i : items)
    if (i.type == ItemType::Frame && gemmi::iequal(i.frame.name, frame_name))
      return &i.frame;
  return nullptr;
}

inline Table Block::item_as_table(Item& item) {
  if (item.type != ItemType::Loop)
    fail("item_as_table: item is not Loop");
  std::vector<int> indices(item.loop.tags.size());
  for (size_t j = 0; j != indices.size(); ++j)
    indices[j] = (int) j;
  return Table{&item, *this, indices, 0};
}

inline size_t Block::get_index(const std::string& tag) const {
  for (size_t i = 0; i != items.size(); ++i) {
    const Item& item = items[i];
    if ((item.type == ItemType::Pair && item.pair[0] == tag) ||
        (item.type == ItemType::Loop && item.loop.find_tag(tag) != -1))
      return i;
  }
  fail(tag + " not found in block");
}

inline void Block::move_item(int old_pos, int new_pos) {
  if (old_pos < 0)
    old_pos += items.size();
  if ((size_t) old_pos >= items.size())
    fail("move_item: old_pos out of range");
  if (new_pos < 0)
    new_pos += items.size();
  if ((size_t) new_pos >= items.size())
    fail("move_item: new_pos out of range");
  auto src = items.begin() + old_pos;
  auto dst = items.begin() + new_pos;
  if (src < dst)
    std::rotate(src, src+1, dst+1);
  else
    std::rotate(dst, src, src+1);
}

inline std::vector<std::string> Block::get_mmcif_category_names() const {
  std::vector<std::string> cats;
  for (const Item& item : items) {
    const std::string* tag = nullptr;
    if (item.type == ItemType::Pair)
      tag = &item.pair[0];
    else if (item.type == ItemType::Loop && !item.loop.tags.empty())
      tag = &item.loop.tags[0];
    if (tag)
      for (auto j = cats.rbegin(); j != cats.rend(); ++j)
        if (gemmi::starts_with(*tag, *j)) {
          tag = nullptr;
          break;
        }
    if (tag) {
      size_t dot = tag->find('.');
      if (dot != std::string::npos)
        cats.emplace_back(*tag, 0, dot + 1);
    }
  }
  return cats;
}

inline Item& Block::setup_loop_item(Table&& tab, const std::string& prefix,
                                    std::vector<std::string>&& tags) {
  Item *item;
  if (tab.loop_item) {
    item = tab.loop_item;
    item->loop.clear();
  } else if (tab.ok()) {
    item = &tab.bloc.items.at(tab.positions[0]);
    tab.erase();
    item->set_value(Item(LoopArg{}));
  } else {
    items.emplace_back(LoopArg{});
    item = &items.back();
  }
  for (std::string& tag : tags) {
    tag.insert(0, prefix);
    assert_tag(tag);
  }
  item->loop.tags = std::move(tags);
  return *item;
}

inline Loop& Block::setup_loop(Table&& tab, const std::string& prefix,
                               std::vector<std::string>&& tags) {
  return setup_loop_item(std::move(tab), prefix, std::move(tags)).loop;
}

inline Table Block::find(const std::string& prefix,
                         const std::vector<std::string>& tags) {
  Item* loop_item = nullptr;
  if (!tags.empty()) {
    if (tags[0][0] == '?')
      fail("The first tag in find() cannot be ?optional.");
    loop_item = find_loop(prefix + tags[0]).item();
  }

  std::vector<int> indices;
  indices.reserve(tags.size());
  if (loop_item) {
    for (const std::string& tag : tags) {
      std::string full_tag = prefix + (tag[0] != '?' ? tag : tag.substr(1));
      int idx = loop_item->loop.find_tag(full_tag);
      if (idx == -1 && tag[0] != '?') {
        loop_item = nullptr;
        indices.clear();
        break;
      }
      indices.push_back(idx);
    }
  } else {
    for (const std::string& tag : tags) {
      std::string full_tag = prefix + (tag[0] != '?' ? tag : tag.substr(1));
      if (const Item* p = find_pair_item(full_tag)) {
        indices.push_back(p - items.data());
      } else if (tag[0] == '?') {
        indices.push_back(-1);
      } else {
        indices.clear();
        break;
      }
    }
  }
  return Table{loop_item, *this, indices, prefix.length()};
}

inline Table Block::find_any(const std::string& prefix,
                             const std::vector<std::string>& tags) {
  std::vector<int> indices;
  for (auto tag = tags.begin(); tag != tags.end(); ++tag) {
    Column column = find_values(prefix + *tag);
    if (Item* item = column.item()) {
      if (item->type == ItemType::Loop) {
        indices.push_back(column.col());
        while (++tag != tags.end()) {
          int idx = item->loop.find_tag(prefix + *tag);
          if (idx != -1)
            indices.push_back(idx);
        }
        return Table{item, *this, indices, prefix.length()};
      } else {
        indices.push_back(item - items.data());
        while (++tag != tags.end())
          if (const Item* p = find_pair_item(prefix + *tag))
            indices.push_back(p - items.data());
        return Table{nullptr, *this, indices, prefix.length()};
      }
    }
  }
  return Table{nullptr, *this, indices, prefix.length()};
}

inline Table Block::find_mmcif_category(std::string cat) {
  ensure_mmcif_category(cat);
  std::vector<int> indices;
  for (Item& i : items)
    if (i.has_prefix(cat)) {
      if (i.type == ItemType::Loop) {
        indices.resize(i.loop.tags.size());
        for (size_t j = 0; j != indices.size(); ++j) {
          indices[j] = j;
          const std::string& tag = i.loop.tags[j];
          if (!starts_with(tag, cat))
            fail("Tag " + tag + " in loop with " + cat);
        }
        return Table{&i, *this, indices, cat.length()};
      } else {
        indices.push_back(&i - items.data());
      }
    }
  return Table{nullptr, *this, indices, cat.length()};
}


struct Document {
  std::string source;
  std::vector<Block> blocks;

  // implementation detail: items of the currently parsed block or frame
  std::vector<Item>* items_ = nullptr;

  Block& add_new_block(const std::string& name, int pos=-1) {
    if (find_block(name))
      fail("Block with such name already exists: " + name);
    if (pos > 0 && static_cast<size_t>(pos) > blocks.size())
      throw std::out_of_range("add_new_block(): invalid position");
    return *blocks.emplace(pos < 0 ? blocks.end() : blocks.begin() + pos, name);
  }

  void clear() noexcept {
    source.clear();
    blocks.clear();
    items_ = nullptr;
  }

  // returns blocks[0] if the document has exactly one block (like mmCIF)
  Block& sole_block() {
    if (blocks.size() > 1)
      fail("single data block expected, got " + std::to_string(blocks.size()));
    return blocks.at(0);
  }
  const Block& sole_block() const {
    return const_cast<Document*>(this)->sole_block();
  }

  Block* find_block(const std::string& name) {
    for (Block& b : blocks)
      if (b.name == name)
        return &b;
    return nullptr;
  }
  const Block* find_block(const std::string& name) const {
    return const_cast<Document*>(this)->find_block(name);
  }
};


[[noreturn]]
inline void cif_fail(const std::string& source, const Block& b,
                     const Item& item, const std::string& s) {
  fail(tostr(source, ':', item.line_number, " in data_", b.name, ": ", s));
}

inline void check_for_missing_values_in_block(const Block& block,
                                              const std::string& source) {
  for (const Item& item : block.items) {
    if (item.type == ItemType::Pair) {
      if (item.pair[1].empty())
        cif_fail(source, block, item, item.pair[0] + " has no value");
    } else if (item.type == ItemType::Frame) {
      check_for_missing_values_in_block(item.frame, source);
    }
  }
}

// Throw an error if any item (pair) value is missing
inline void check_for_missing_values(const Document& d) {
  for (const Block& block : d.blocks)
    check_for_missing_values_in_block(block, d.source);
}

// Throw an error if any block name, frame name or tag is duplicated.
inline void check_for_duplicates(const Document& d) {
  // check for duplicate block names (except empty "" which is global_)
  std::unordered_set<std::string> names;
  for (const Block& block : d.blocks) {
    bool ok = names.insert(gemmi::to_lower(block.name)).second;
    if (!ok && !block.name.empty())
      fail(d.source + ": duplicate block name: ", block.name);
  }
  // check for dups inside each block
  std::unordered_set<std::string> frame_names;
  for (const Block& block : d.blocks) {
    names.clear();
    frame_names.clear();
    for (const Item& item : block.items) {
      if (item.type == ItemType::Pair) {
        bool ok = names.insert(gemmi::to_lower(item.pair[0])).second;
        if (!ok)
          cif_fail(d.source, block, item, "duplicate tag " + item.pair[0]);
      } else if (item.type == ItemType::Loop) {
        for (const std::string& t : item.loop.tags) {
          bool ok = names.insert(gemmi::to_lower(t)).second;
          if (!ok)
            cif_fail(d.source, block, item, "duplicate tag " + t);
        }
      } else if (item.type == ItemType::Frame) {
        bool ok = frame_names.insert(gemmi::to_lower(item.frame.name)).second;
        if (!ok)
          cif_fail(d.source, block, item, "duplicate save_" + item.frame.name);
      }
    }
  }
}

inline bool is_text_field(const std::string& val) {
  size_t len = val.size();
  return len > 2 && val[0] == ';' && (val[len-2] == '\n' || val[len-2] == '\r');
}

inline std::string quote(std::string v) {
  if (std::all_of(v.begin(), v.end(), [](char c) { return char_table(c) == 1; })
      && !v.empty() && !is_null(v))
    return v;
  char q = ';';
  if (std::memchr(v.c_str(), '\n', v.size()) == nullptr) {
    if (std::memchr(v.c_str(), '\'', v.size()) == nullptr)
      q = '\'';
    else if (std::memchr(v.c_str(), '"', v.size()) == nullptr)
      q = '"';
  }
  v.insert(v.begin(), q);
  if (q == ';')
    v += '\n';
  v += q;
  return v;
}

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

} // namespace cif
} // namespace gemmi
#endif
