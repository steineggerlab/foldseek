// Copyright 2018 Global Phasing Ltd.
//
// Selections.

#ifndef GEMMI_SELECT_HPP_
#define GEMMI_SELECT_HPP_

#include <string>
#include <cstdlib>   // for strtol
#include <cctype>    // for isalpha
#include <climits>   // for INT_MIN, INT_MAX
#include "fail.hpp"  // for fail
#include "util.hpp"  // for is_in_list
#include "model.hpp" // for Model, Chain, etc
#include "iterator.hpp" // for FilterProxy
#include "sprintf.hpp" // for to_str
#include "atof.hpp"  // for fast_from_chars

namespace gemmi {

// from http://www.ccp4.ac.uk/html/pdbcur.html
// Specification of the selection sets:
// either
//     /mdl/chn/s1.i1-s2.i2/at[el]:aloc
// or
//     /mdl/chn/*(res).ic/at[el]:aloc
//

struct Selection {
  struct List {
    bool all = true;
    bool inverted = false;
    std::string list;  // comma-separated

    std::string str() const {
      if (all)
        return "*";
      return inverted ? "!" + list : list;
    }

    bool has(const std::string& name) const {
      if (all)
        return true;
      bool found = is_in_list(name, list);
      return inverted ? !found : found;
    }
  };

  struct FlagList {
    std::string pattern;
    bool has(char flag) const {
      if (pattern.empty())
        return true;
      bool invert = (pattern[0] == '!');
      bool found = (pattern.find(flag, invert ? 1 : 0) != std::string::npos);
      return invert ? !found : found;
    }
  };

  struct SequenceId {
    int seqnum;
    char icode;

    bool empty() const {
      return seqnum == INT_MIN || seqnum == INT_MAX;
    }

    std::string str() const {
      std::string s;
      if (!empty()) {
        s = std::to_string(seqnum);
        if (icode != '*') {
          s += '.';
          if (icode != ' ')
            s += icode;
        }
      }
      return s;
    }

    int compare(const SeqId& seqid) const {
      if (seqnum != *seqid.num)
        return seqnum < *seqid.num ? -1 : 1;
      if (icode != '*' && icode != seqid.icode)
        return icode < seqid.icode ? -1 : 1;
      return 0;
    }
  };

  struct AtomInequality {
    char property;
    int relation;
    double value;

    bool matches(const Atom& a) const {
      double atom_value = 0.;
      if (property == 'q')
        atom_value = a.occ;
      else if (property == 'b')
        atom_value = a.b_iso;
      if (relation < 0)
        return atom_value < value;
      if (relation > 0)
        return atom_value > value;
      return atom_value == value;
    }

    std::string str() const {
      std::string r = ";";
      r += property;
      r += relation == 0 ? '=' : relation < 0 ? '<' : '>';
      r += to_str(value);
      return r;
    }
  };

  int mdl = 0;            // 0 = all
  List chain_ids;
  SequenceId from_seqid = {INT_MIN, '*'};
  SequenceId to_seqid = {INT_MAX, '*'};
  List residue_names;
  List atom_names;
  List elements;
  List altlocs;
  FlagList residue_flags;
  FlagList atom_flags;
  std::vector<AtomInequality> atom_inequalities;

  Selection() = default;
  Selection(const std::string& cid);

  std::string str() const {
    std::string cid = "/";
    if (mdl != 0)
      cid += std::to_string(mdl);
    cid += '/';
    cid += chain_ids.str();
    cid += '/';
    cid += from_seqid.str();
    if (!residue_names.all) {
      cid += '(';
      cid += residue_names.str();
      cid += ')';
    }
    if (!from_seqid.empty() || !to_seqid.empty()) {
      cid += '-';
      cid += to_seqid.str();
    }
    cid += '/';
    if (!atom_names.all)
      cid += atom_names.str();
    if (!elements.all) {
      cid += '[';
      cid += elements.str();
      cid += ']';
    }
    if (!altlocs.all) {
      cid += ':';
      cid += altlocs.str();
    }
    for (const AtomInequality& ai : atom_inequalities)
      cid += ai.str();
    return cid;
  }

  bool matches(const Model& model) const {
    return mdl == 0 || std::to_string(mdl) == model.name;
  }
  bool matches(const Chain& chain) const {
    return chain_ids.has(chain.name);
  }
  bool matches(const Residue& res) const {
    return residue_names.has(res.name) &&
           from_seqid.compare(res.seqid) <= 0 &&
           to_seqid.compare(res.seqid) >= 0 &&
           residue_flags.has(res.flag);
  }
  bool matches(const Atom& a) const {
    return atom_names.has(a.name) &&
           (elements.all || elements.has(a.element.uname())) &&
           (altlocs.all || altlocs.has(std::string(a.altloc ? 1 : 0, a.altloc))) &&
           atom_flags.has(a.flag) &&
           std::all_of(atom_inequalities.begin(), atom_inequalities.end(),
                       [&](const AtomInequality& i) { return i.matches(a); });
  }
  bool matches(const CRA& cra) const {
    return (cra.chain == nullptr || matches(*cra.chain)) &&
           (cra.residue == nullptr || matches(*cra.residue)) &&
           (cra.atom == nullptr || matches(*cra.atom));
  }

  FilterProxy<Selection, Model> models(Structure& st) const {
    return {*this, st.models};
  }
  FilterProxy<Selection, Chain> chains(Model& model) const {
    return {*this, model.chains};
  }
  FilterProxy<Selection, Residue> residues(Chain& chain) const {
    return {*this, chain.residues};
  }
  FilterProxy<Selection, Atom> atoms(Residue& residue) const {
    return {*this, residue.atoms};
  }

  CRA first_in_model(Model& model) const {
    if (matches(model))
      for (Chain& chain : model.chains) {
        if (matches(chain))
          for (Residue& res : chain.residues) {
            if (matches(res))
              for (Atom& atom : res.atoms) {
                if (matches(atom))
                  return {&chain, &res, &atom};
              }
          }
      }
    return {nullptr, nullptr, nullptr};
  }

  std::pair<Model*, CRA> first(Structure& st) const {
    for (Model& model : st.models) {
      CRA cra = first_in_model(model);
      if (cra.chain)
        return {&model, cra};
    }
    return {nullptr, {nullptr, nullptr, nullptr}};
  }

  template<typename T>
  void add_matching_children(const T& orig, T& target) const {
    for (const auto& orig_child : orig.children())
      if (matches(orig_child)) {
        target.children().push_back(orig_child.empty_copy());
        add_matching_children(orig_child, target.children().back());
      }
  }
  void add_matching_children(const Atom&, Atom&) const {}

  Selection& set_residue_flags(const std::string& pattern) {
    residue_flags.pattern = pattern;
    return *this;
  }
  Selection& set_atom_flags(const std::string& pattern) {
    atom_flags.pattern = pattern;
    return *this;
  }

  template<typename T>
  T copy_selection(const T& orig) const {
    T copied = orig.empty_copy();
    add_matching_children(orig, copied);
    return copied;
  }

  template<typename T>
  void remove_selected(T& t) const {
    for (auto& child : t.children())
      if (matches(child))
        remove_selected(child);
    vector_remove_if(t.children(),
                     [&](typename T::child_type& c) { return c.children().empty(); });
  }
  void remove_selected(Residue& res) const {
    if (atom_names.all && elements.all && altlocs.all &&
        atom_flags.pattern.empty() && atom_inequalities.empty())
      res.atoms.clear();
    else
      vector_remove_if(res.atoms, [&](Atom& c) { return matches(c); });
  }

  template<typename T>
  void remove_not_selected(T& t) const {
    vector_remove_if(t.children(), [&](typename T::child_type& c) { return !matches(c); });
    for (auto& child : t.children())
      remove_not_selected(child);
  }
  void remove_not_selected(Atom&) const {}
};

namespace impl {

inline int determine_omitted_cid_fields(const std::string& cid) {
  if (cid[0] == '/')
    return 0; // model
  if (std::isdigit(cid[0]) || cid[0] == '.' || cid[0] == '(' || cid[0] == '-')
    return 2; // residue
  size_t sep = cid.find_first_of("/([:;");
  if (sep == std::string::npos || cid[sep] == '/')
    return 1; // chain
  if (cid[sep] == '(')
    return 2; // residue
  return 3;  // atom
}

inline Selection::List make_cid_list(const std::string& cid, size_t pos, size_t end) {
  Selection::List list;
  list.all = (cid[pos] == '*');
  if (cid[pos] == '!') {
    list.inverted = true;
    ++pos;
  }
  list.list = cid.substr(pos, end - pos);
  return list;
}

inline Selection::SequenceId parse_cid_seqid(const std::string& cid, size_t& pos,
                                             int default_seqnum) {
  size_t initial_pos = pos;
  int seqnum = default_seqnum;
  char icode = ' ';
  if (cid[pos] == '*') {
    ++pos;
    icode = '*';
  } else if (std::isdigit(cid[pos])) {
    char* endptr;
    seqnum = std::strtol(&cid[pos], &endptr, 10);
    pos = endptr - &cid[0];
  }
  if (cid[pos] == '.')
    ++pos;
  if (initial_pos != pos && (std::isalpha(cid[pos]) || cid[pos] == '*'))
    icode = cid[pos++];
  return {seqnum, icode};
}

inline Selection::AtomInequality parse_atom_inequality(const std::string& cid,
                                                       size_t pos, size_t end) {
  Selection::AtomInequality r;
  while (cid[pos] == ' ')
    ++pos;
  if (cid[pos] != 'q' && cid[pos] != 'b')
    fail("Invalid selection syntax (at ", cid[pos], "): ", cid);
  r.property = cid[pos];
  ++pos;
  while (cid[pos] == ' ')
    ++pos;
  if (cid[pos] == '<')
    r.relation = -1;
  else if (cid[pos] == '>')
    r.relation = 1;
  else if (cid[pos] == '=')
    r.relation = 0;
  else
    fail("Invalid selection syntax (at ", cid[pos], "): ", cid);
  ++pos;
  auto result = fast_from_chars(cid.c_str() + pos, r.value);
  if (result.ec != std::errc())
    fail("Invalid selection syntax (number expected at '", cid.substr(pos), "'): ", cid);
  pos = size_t(result.ptr - cid.c_str());
  while (cid[pos] == ' ')
    ++pos;
  if (pos != std::min(end, cid.size()))
    fail("Invalid selection syntax (at ", cid[pos], "): ", cid);
  return r;
}

inline void parse_cid(const std::string& cid, Selection& sel) {
  if (cid.empty() || (cid.size() == 1 && cid[0] == '*'))
    return;
  int omit = determine_omitted_cid_fields(cid);
  size_t sep = 0;
  // model
  if (omit == 0) {
    sep = cid.find('/', 1);
    if (sep != 1 && cid[1] != '*') {
      char* endptr;
      sel.mdl = std::strtol(&cid[1], &endptr, 10);
      size_t end_pos = endptr - &cid[0];
      if (end_pos != sep && end_pos != cid.length())
        fail("Expected model number first: " + cid);
    }
  }

  // chain
  if (omit <= 1 && sep != std::string::npos) {
    size_t pos = (sep == 0 ? 0 : sep + 1);
    sep = cid.find('/', pos);
    sel.chain_ids = make_cid_list(cid, pos, sep);
  }

  // residue; MMDB CID syntax: s1.i1-s2.i2 or *(res).ic
  // In gemmi both 14.a and 14a are accepted.
  // *(ALA). and *(ALA) and (ALA). can be used instead of (ALA) for
  // compatibility with MMDB.
  if (omit <= 2 && sep != std::string::npos) {
    size_t pos = (sep == 0 ? 0 : sep + 1);
    if (cid[pos] != '(')
      sel.from_seqid = parse_cid_seqid(cid, pos, INT_MIN);
    if (cid[pos] == '(') {
      ++pos;
      size_t right_br = cid.find(')', pos);
      sel.residue_names = make_cid_list(cid, pos, right_br);
      pos = right_br + 1;
    }
    // allow "(RES)." and "(RES).*" and "(RES)*"
    if (cid[pos] == '.')
      ++pos;
    if (cid[pos] == '*')
      ++pos;
    if (cid[pos] == '-') {
      ++pos;
      sel.to_seqid = parse_cid_seqid(cid, pos, INT_MAX);
    }
    sep = pos;
  }

  // atom;  at[el]:aloc
  if (sep < cid.size()) {
    if (sep != 0 && cid[sep] != '/')
      fail("Invalid selection syntax: " + cid);
    size_t pos = (sep == 0 ? 0 : sep + 1);
    size_t end = cid.find_first_of("[:;", pos);
    if (end != pos) {
      sel.atom_names = make_cid_list(cid, pos, end);
      // Chain name can be empty, but not atom name,
      // so we interpret empty atom name as *.
      if (!sel.atom_names.inverted && sel.atom_names.list.empty())
        sel.atom_names.all = true;
    }
    if (end != std::string::npos) {
      if (cid[end] == '[') {
        pos = end + 1;
        end = cid.find(']', pos);
        if (end == std::string::npos)
          fail("Invalid selection syntax (no matching ']'): " + cid);
        sel.elements = make_cid_list(cid, pos, end);
        sel.elements.list = to_upper(sel.elements.list);
        ++end;
      }
      if (cid[end] == ':') {
        pos = end + 1;
        end = cid.find(';', pos);
        sel.altlocs = make_cid_list(cid, pos, end);
      }
      while (end != std::string::npos && cid[end] == ';') {
        pos = end + 1;
        end = cid.find(';', pos);
        sel.atom_inequalities.push_back(parse_atom_inequality(cid, pos, end));
      }
      if (end < cid.length())
        fail("Invalid selection syntax (atom properties): " + cid);
    }
  }
}

} // namespace impl


inline Selection::Selection(const std::string& cid) {
  impl::parse_cid(cid, *this);
}

} // namespace gemmi
#endif
