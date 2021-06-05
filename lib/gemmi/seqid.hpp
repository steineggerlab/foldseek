// Copyright 2017 Global Phasing Ltd.
//
// SeqId -- residue number and insertion code together.

#ifndef GEMMI_SEQID_HPP_
#define GEMMI_SEQID_HPP_

#include <cstdlib>    // for strtol
#include <stdexcept>  // for invalid_argument
#include <string>

namespace gemmi {

// Optional int value. N is a special value that means not-set.
template<int N> struct OptionalInt {
  enum { None=N };
  int value = None;

  OptionalInt() = default;
  OptionalInt(int n) : value(n) {}
  bool has_value() const { return value != None; }
  std::string str(char null='?') const {
    return has_value() ? std::to_string(value) : std::string(1, null);
  }
  OptionalInt& operator=(int n) { value = n; return *this; }
  bool operator==(const OptionalInt& o) const { return value == o.value; }
  bool operator!=(const OptionalInt& o) const { return value != o.value; }
  bool operator<(const OptionalInt& o) const {
    return has_value() && o.has_value() && value < o.value;
  }
  bool operator==(int n) const { return value == n; }
  bool operator!=(int n) const { return value != n; }
  OptionalInt operator+(OptionalInt o) const {
    return OptionalInt(has_value() && o.has_value() ? value + o.value : N);
  }
  OptionalInt operator-(OptionalInt o) const {
    return OptionalInt(has_value() && o.has_value() ? value - o.value : N);
  }
  OptionalInt& operator+=(int n) { if (has_value()) value += n; return *this; }
  OptionalInt& operator-=(int n) { if (has_value()) value -= n; return *this; }
  explicit operator int() const { return value; }
  explicit operator bool() const { return has_value(); }
  // these are defined for partial compatibility with C++17 std::optional
  using value_type = int;
  int& operator*() { return value; }
  const int& operator*() const { return value; }
  int& emplace(int n) { value = n; return value; }
};

struct SeqId {
  using OptionalNum = OptionalInt<-999>;

  OptionalNum num; // sequence number
  char icode = ' ';  // insertion code

  SeqId() = default;
  SeqId(int num_, char icode_) { num = num_; icode = icode_; }
  SeqId(OptionalNum num_, char icode_) { num = num_; icode = icode_; }
  explicit SeqId(const std::string& str) {
    char* endptr;
    num = std::strtol(str.c_str(), &endptr, 10);
    if (endptr == str.c_str() || (*endptr != '\0' && endptr[1] != '\0'))
      throw std::invalid_argument("Not a seqid: " + str);
    icode = (*endptr | 0x20);
  }

  bool operator==(const SeqId& o) const {
    return num == o.num && (icode | 0x20) == (o.icode | 0x20);
  }
  bool operator!=(const SeqId& o) const { return !operator==(o); }
  bool operator<(const SeqId& o) const {
    return (*num * 256 + icode) < (*o.num * 256 + o.icode);
  }

  char has_icode() const { return icode != ' '; }

  std::string str() const {
    std::string r = num.str();
    if (icode != ' ')
      r += icode;
    return r;
  }
};

// Sequence ID (sequence number + insertion code) + residue name + segment ID
struct ResidueId {
  SeqId seqid;
  std::string segment; // segid - up to 4 characters in the PDB file
  std::string name;

  // used for first_conformation iterators, etc.
  SeqId group_key() const { return seqid; }

  bool matches(const ResidueId& o) const {
    return seqid == o.seqid && segment == o.segment && name == o.name;
  }
  bool matches_noseg(const ResidueId& o) const {
    return seqid == o.seqid && name == o.name;
  }
  bool operator==(const ResidueId& o) const { return matches(o); }
  std::string str() const { return seqid.str() + "(" + name + ")"; }
};

inline std::string atom_str(const std::string& chain_name,
                            const ResidueId& res_id,
                            const std::string& atom_name,
                            char altloc) {
  std::string r = chain_name;
  r += '/';
  r += res_id.name;
  r += ' ';
  r += res_id.seqid.str();
  r += '/';
  r += atom_name;
  if (altloc) {
    r += '.';
    r += altloc;
  }
  return r;
}

struct AtomAddress {
  std::string chain_name;
  ResidueId res_id;
  std::string atom_name;
  char altloc = '\0';

  AtomAddress() = default;
  AtomAddress(const std::string& ch, const ResidueId& resid,
              const std::string& atom, char alt='\0')
    : chain_name(ch), res_id(resid), atom_name(atom), altloc(alt) {}
  AtomAddress(const std::string& ch, const SeqId& seqid, const std::string& res,
              const std::string& atom, char alt='\0')
    : chain_name(ch), res_id({seqid, "", res}), atom_name(atom), altloc(alt) {}
  bool operator==(const AtomAddress& o) const {
    return chain_name == o.chain_name && res_id.matches(o.res_id) &&
           atom_name == o.atom_name && altloc == o.altloc;
  }

  std::string str() const {
    return atom_str(chain_name, res_id, atom_name, altloc);
  }
};

} // namespace gemmi

namespace std {
template <> struct hash<gemmi::ResidueId> {
  size_t operator()(const gemmi::ResidueId& r) const {
    size_t seqid_hash = (*r.seqid.num << 7) + (r.seqid.icode | 0x20);
    return seqid_hash ^ hash<string>()(r.segment) ^ hash<string>()(r.name);
  }
};
} // namespace std

#endif
