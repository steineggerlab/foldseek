// Copyright 2018 Global Phasing Ltd.
//
// ChemComp - chemical component that represents a monomer from Refmac
// monomer library, or from PDB CCD.

#ifndef GEMMI_CHEMCOMP_HPP_
#define GEMMI_CHEMCOMP_HPP_

#include <string>
#include <vector>
#include "cifdoc.hpp"
#include "elem.hpp"  // for Element
#include "fail.hpp"  // for fail, unreachable
#include "numb.hpp"  // for as_number
#include "util.hpp"  // for istarts_with, join_str
#include "model.hpp" // for Residue, Atom

namespace gemmi {

enum class BondType {
  Unspec, Single, Double, Triple, Aromatic, Deloc, Metal
};
enum class ChiralityType { Positive, Negative, Both };

struct Restraints {
  struct AtomId {
    int comp;
    std::string atom;

    bool operator==(const AtomId& o) const {
      return comp == o.comp && atom == o.atom;
    }
    bool operator!=(const AtomId& o) const { return !operator==(o); }

    bool operator==(const std::string& name) const { return atom == name; }
    bool operator!=(const std::string& name) const { return atom != name; }

    bool operator<(const AtomId& o) const {
      return comp == o.comp ? atom < o.atom : comp < o.comp;
    }

    const Atom* get_from(const Residue& res1, const Residue* res2,
                         char altloc) const {
      const Residue* residue;
      if (comp == 1 || res2 == nullptr)
        residue = &res1;
      else if (comp == 2)
        residue = res2;
      else
        throw std::out_of_range("Unexpected component ID");
      if (const Atom* ret = residue->find_atom(atom, altloc))
        // skip riding hydrogens, they won't be restrained (to be revised)
        if (ret->calc_flag != CalcFlag::Calculated || !ret->is_hydrogen())
          return ret;
      return nullptr;
    }
    Atom* get_from(Residue& res1, Residue* res2, char altloc) const {
      const Residue& cres1 = res1;
      return const_cast<Atom*>(get_from(cres1, res2, altloc));
    }
  };

  static std::string lexicographic_str(const std::string& name1,
                                       const std::string& name2) {
    return name1 < name2 ? name1 + "-" + name2 : name2 + "-" + name1;
  }

  enum class DistanceOf { ElectronCloud, Nucleus };

  struct Bond {
    AtomId id1, id2;
    BondType type;
    bool aromatic;
    double value;
    double esd;
    double value_nucleus;
    double esd_nucleus;
    std::string str() const { return id1.atom + "-" + id2.atom; }
    std::string lexicographic_str() const {
      return Restraints::lexicographic_str(id1.atom, id2.atom);
    }
    double distance(DistanceOf of) const {
      return of == DistanceOf::ElectronCloud ? value : value_nucleus;
    }
  };

  struct Angle {
    AtomId id1, id2, id3;
    double value;
    double esd;
    double radians() const { return rad(value); }
    std::string str() const {
      return id1.atom + "-" + id2.atom + "-" + id3.atom;
    }
  };

  struct Torsion {
    std::string label;
    AtomId id1, id2, id3, id4;
    double value;
    double esd;
    int period;
    std::string str() const {
      return id1.atom + "-" + id2.atom + "-" + id3.atom + "-" + id4.atom;
    }
  };

  struct Chirality {
    AtomId id_ctr, id1, id2, id3;
    ChiralityType sign;

    bool is_wrong(double volume) const {
      return (sign == ChiralityType::Positive && volume < 0) ||
             (sign == ChiralityType::Negative && volume > 0);
    }
    std::string str() const {
      return id_ctr.atom + "," + id1.atom + "," + id2.atom + "," + id3.atom;
    }
  };

  struct Plane {
    std::string label;
    std::vector<AtomId> ids;
    double esd;
    std::string str() const {
      return join_str(ids, ',', [](const AtomId& a) { return a.atom; });
    }
  };

  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::vector<Plane> planes;

  bool empty() const {
    return bonds.empty() && angles.empty() && torsions.empty() &&
           chirs.empty() && planes.empty();
  }

  template<typename T>
  std::vector<Bond>::iterator find_bond(const T& a1, const T& a2) {
    return std::find_if(bonds.begin(), bonds.end(), [&](const Bond& b) {
        return (b.id1 == a1 && b.id2 == a2) || (b.id1 == a2 && b.id2 == a1);
    });
  }
  template<typename T>
  std::vector<Bond>::const_iterator find_bond(const T& a1, const T& a2) const {
    return const_cast<Restraints*>(this)->find_bond(a1, a2);
  }
  const Bond& get_bond(const AtomId& a1, const AtomId& a2) const {
    auto it = find_bond(a1, a2);
    if (it == bonds.end())
      fail("Bond restraint not found: ", a1.atom, '-', a2.atom);
    return *it;
  }

  template<typename T>
  bool are_bonded(const T& a1, const T& a2) const {
    return find_bond(a1, a2) != bonds.end();
  }

  template<typename T>
  const AtomId* first_bonded_atom(const T& a) const {
    for (const Bond& bond : bonds) {
      if (bond.id1 == a)
        return &bond.id2;
      if (bond.id2.atom == a)
        return &bond.id1;
    }
    return nullptr;
  }

  template<typename A, typename T>
  void for_each_bonded_atom(const A& a, const T& func) const {
    for (const Bond& bond : bonds) {
      const AtomId* other = nullptr;
      if (bond.id1 == a)
        other = &bond.id2;
      else if (bond.id2 == a)
        other = &bond.id1;
      if (other != nullptr)
        if (!func(*other))
          break;
    }
  }

  // BFS
  std::vector<AtomId> find_shortest_path(const AtomId& a, const AtomId& b,
                                         std::vector<AtomId> visited) const {
    int start = (int) visited.size();
    int end = -1;
    visited.push_back(b);
    std::vector<int> parent(visited.size(), -1);
    for (int n = start; end == -1 && n != (int) visited.size(); ++n) {
      for_each_bonded_atom(AtomId(visited[n]), [&](const AtomId& id) {
          if (id == a)
            end = (int) visited.size();
          if (!in_vector(id, visited)) {
            visited.push_back(id);
            parent.push_back(n);
          }
          return true;
      });
    }
    std::vector<AtomId> path;
    for (int n = end; n != -1; n = parent[n])
      path.push_back(visited[n]);
    return path;
  }

  template<typename T>
  std::vector<Angle>::iterator find_angle(const T& a, const T& b, const T& c) {
    return std::find_if(angles.begin(), angles.end(), [&](const Angle& ang) {
        return ang.id2 == b && ((ang.id1 == a && ang.id3 == c) ||
                                (ang.id1 == c && ang.id3 == a));
    });
  }
  template<typename T> std::vector<Angle>::const_iterator
  find_angle(const T& a, const T& b, const T& c) const {
    return const_cast<Restraints*>(this)->find_angle(a, b, c);
  }
  const Angle& get_angle(const AtomId& a, const AtomId& b, const AtomId& c)
                                                                        const {
    auto it = const_cast<Restraints*>(this)->find_angle(a, b, c);
    if (it == angles.end())
      fail("Angle restraint not found: ", a.atom, '-', b.atom, '-', c.atom);
    return *it;
  }

  template<typename T>
  std::vector<Torsion>::iterator find_torsion(const T& a, const T& b,
                                              const T& c, const T& d) {
    return std::find_if(torsions.begin(), torsions.end(),
                        [&](const Torsion& t) {
        return (t.id1 == a && t.id2 == b && t.id3 == c && t.id4 == d) ||
               (t.id1 == d && t.id2 == c && t.id3 == b && t.id4 == a);
    });
  }
  template<typename T> std::vector<Torsion>::const_iterator
  find_torsion(const T& a, const T& b, const T& c, const T& d) const {
    return const_cast<Restraints*>(this)->find_torsion(a, b, c, d);
  }

  template<typename T>
  std::vector<Chirality>::iterator find_chir(const T& ctr, const T& a,
                                             const T& b, const T& c) {
    return std::find_if(chirs.begin(), chirs.end(), [&](const Chirality& t) {
        return t.id_ctr == ctr && ((t.id1 == a && t.id2 == b && t.id3 == c) ||
                                   (t.id1 == b && t.id2 == c && t.id3 == a) ||
                                   (t.id1 == c && t.id2 == a && t.id3 == b));
    });
  }
  template<typename T> std::vector<Chirality>::const_iterator
  find_chir(const T& ctr, const T& a, const T& b, const T& c) const {
    return const_cast<Restraints*>(this)->find_chir(ctr, a, b, c);
  }

  double chiral_abs_volume(const Restraints::Chirality& ch) const {
    double mult = get_bond(ch.id_ctr, ch.id1).value *
                  get_bond(ch.id_ctr, ch.id2).value *
                  get_bond(ch.id_ctr, ch.id3).value;
    double x = 1;
    double y = 2;
    for (double a : {get_angle(ch.id1, ch.id_ctr, ch.id2).value,
                     get_angle(ch.id2, ch.id_ctr, ch.id3).value,
                     get_angle(ch.id3, ch.id_ctr, ch.id1).value}) {
      double cosine = a == 90. ? 0. : std::cos(rad(a));
      x -= cosine * cosine;
      y *= cosine;
    }
    return mult * std::sqrt(x + y);
  }

  std::vector<Plane>::iterator get_plane(const std::string& label) {
    return std::find_if(planes.begin(), planes.end(),
                        [&label](const Plane& p) { return p.label == label; });
  }

  Plane& get_or_add_plane(const std::string& label) {
    std::vector<Plane>::iterator it = get_plane(label);
    if (it != planes.end())
      return *it;
    planes.push_back(Plane{label, {}, 0.0});
    return planes.back();
  }
};

template<typename Restr>
double angle_z(double value_rad, const Restr& restr, double full=360.) {
  return angle_abs_diff(deg(value_rad), restr.value, full) / restr.esd;
}

struct ChemComp {
  struct Atom {
    std::string id;
    Element el;
    // _chem_comp_atom.partial_charge can be non-integer,
    // _chem_comp_atom.charge is always integer (but sometimes has format
    //  '0.000' which is not correct but we ignore it).
    float charge;
    std::string chem_type;

    gemmi::Atom to_full_atom() const {
      gemmi::Atom atom;
      atom.name = id;
      atom.calc_flag = CalcFlag::Calculated;
      atom.occ = 0.0f;
      atom.b_iso = 0.0f;
      atom.element = el;
      atom.charge = static_cast<signed char>(std::round(charge));
      return atom;
    }
    bool is_hydrogen() const { return gemmi::is_hydrogen(el); }
  };

  std::string name;
  std::string group;
  std::vector<Atom> atoms;
  Restraints rt;

  std::vector<Atom>::iterator find_atom(const std::string& atom_id) {
    return std::find_if(atoms.begin(), atoms.end(),
                        [&](const Atom& a) { return a.id == atom_id; });
  }
  std::vector<Atom>::const_iterator find_atom(const std::string& atom_id) const{
    return const_cast<ChemComp*>(this)->find_atom(atom_id);
  }

  int get_atom_index(const std::string& atom_id) const {
    auto it = find_atom(atom_id);
    if (it == atoms.end())
      fail("Chemical componenent ", name, " has no atom ", atom_id);
    return int(it - atoms.begin());
  }

  const Atom& get_atom(const std::string& atom_id) const {
    return atoms[get_atom_index(atom_id)];
  }

  void remove_nonmatching_restraints() {
    vector_remove_if(rt.bonds, [&](const Restraints::Bond& x) {
      return find_atom(x.id1.atom) == atoms.end() ||
             find_atom(x.id2.atom) == atoms.end();
    });
    vector_remove_if(rt.angles, [&](const Restraints::Angle& x) {
      return find_atom(x.id1.atom) == atoms.end() ||
             find_atom(x.id2.atom) == atoms.end() ||
             find_atom(x.id3.atom) == atoms.end();
    });
    vector_remove_if(rt.torsions, [&](const Restraints::Torsion& x) {
      return find_atom(x.id1.atom) == atoms.end() ||
             find_atom(x.id2.atom) == atoms.end() ||
             find_atom(x.id3.atom) == atoms.end() ||
             find_atom(x.id4.atom) == atoms.end();
    });
    vector_remove_if(rt.chirs, [&](const Restraints::Chirality& x) {
      return find_atom(x.id_ctr.atom) == atoms.end() ||
             find_atom(x.id1.atom) == atoms.end() ||
             find_atom(x.id2.atom) == atoms.end() ||
             find_atom(x.id3.atom) == atoms.end();
    });
    for (Restraints::Plane& plane : rt.planes)
      vector_remove_if(plane.ids, [&](const Restraints::AtomId& x) {
        return find_atom(x.atom) == atoms.end();
      });
  }

  ChemComp& remove_hydrogens() {
    vector_remove_if(atoms, [](const ChemComp::Atom& a) {
      return a.is_hydrogen();
    });
    remove_nonmatching_restraints();
    return *this;
  }
};

inline BondType bond_type_from_string(const std::string& s) {
  if (istarts_with(s, "sing"))
    return BondType::Single;
  if (istarts_with(s, "doub"))
    return BondType::Double;
  if (istarts_with(s, "trip"))
    return BondType::Triple;
  if (istarts_with(s, "arom"))
    return BondType::Aromatic;
  if (istarts_with(s, "metal"))
    return BondType::Metal;
  if (istarts_with(s, "delo") || s == "1.5")
    return BondType::Deloc;
  if (cif::is_null(s))
    return BondType::Unspec;
  // program PDB2TNT produces a restraint file with bond type 'coval'
  if (s == "coval")
    return BondType::Unspec;
  throw std::out_of_range("Unexpected bond type: " + s);
}

inline const char* bond_type_to_string(BondType btype) {
  switch (btype) {
    case BondType::Unspec: return ".";
    case BondType::Single: return "single";
    case BondType::Double: return "double";
    case BondType::Triple: return "triple";
    case BondType::Aromatic: return "aromatic";
    case BondType::Deloc: return "deloc";
    case BondType::Metal: return "metal";
  }
  unreachable();
}

inline float order_of_bond_type(BondType btype) {
  switch (btype) {
    case BondType::Single: return 1.0f;
    case BondType::Double: return 2.0f;
    case BondType::Triple: return 3.0f;
    case BondType::Aromatic: return 1.5f;
    case BondType::Deloc: return 1.5f;
    case BondType::Metal: return 1.0f;
    case BondType::Unspec: return 0.0f;
  }
  unreachable();
}

// it doesn't handle crossN types from the monomer library
inline ChiralityType chirality_from_string(const std::string& s) {
  switch (s[0] | 0x20) {
    case 'p': return ChiralityType::Positive;
    case 'n': return ChiralityType::Negative;
    case 'b': return ChiralityType::Both;
    default: throw std::out_of_range("Unexpected chirality: " + s);
  }
}

inline const char* chirality_to_string(ChiralityType chir_type) {
  switch (chir_type) {
    case ChiralityType::Positive: return "positive";
    case ChiralityType::Negative: return "negative";
    case ChiralityType::Both: return "both";
  }
  unreachable();
}

inline ChemComp make_chemcomp_from_block(const cif::Block& block_) {
  ChemComp cc;
  cc.name = block_.name.substr(starts_with(block_.name, "comp_") ? 5 : 0);
  cif::Block& block = const_cast<cif::Block&>(block_);
  cif::Column group_col = block.find_values("_chem_comp.group");
  if (group_col)
    cc.group = group_col.str(0);
  for (auto row : block.find("_chem_comp_atom.",
                             {"atom_id", "type_symbol", "?type_energy",
                             "?charge", "?partial_charge"}))
    cc.atoms.push_back({row.str(0), Element(row.str(1)),
                        (float) cif::as_number(row.one_of(3, 4), 0.0),
                        row.has(2) ? row.str(2) : ""});
  for (auto row : block.find("_chem_comp_bond.",
                             {"atom_id_1", "atom_id_2",              // 0, 1
                              "?type", "?value_order",               // 2, 3
                              "?aromatic", "?pdbx_aromatic_flag",    // 4, 5
                              "?value_dist", "?value_dist_esd",      // 6, 7
                              "?value_dist_nucleus", "?value_dist_nucleus_esd"})) { // 8, 9
    bool aromatic_flag = (row.one_of(4, 5)[0] | 0x20) == 'y';
    double dist = row.has(6) ? cif::as_number(row[6]) : NAN;
    double esd = row.has(7) ? cif::as_number(row[7]) : NAN;
    double dist_nucl = row.has(8) ? cif::as_number(row[8]) : NAN;
    double esd_nucl = row.has(9) ? cif::as_number(row[9]) : NAN;
    cc.rt.bonds.push_back({{1, row.str(0)}, {1, row.str(1)},
                          bond_type_from_string(row.one_of(2, 3)),
                          aromatic_flag, dist, esd, dist_nucl, esd_nucl});
  }
  for (auto row : block.find("_chem_comp_angle.",
                             {"atom_id_1", "atom_id_2", "atom_id_3",
                              "value_angle", "value_angle_esd"}))
    cc.rt.angles.push_back({{1, row.str(0)}, {1, row.str(1)}, {1, row.str(2)},
                            cif::as_number(row[3]), cif::as_number(row[4])});
  for (auto row : block.find("_chem_comp_tor.",
                             {"id",
                              "atom_id_1", "atom_id_2",
                              "atom_id_3", "atom_id_4",
                              "value_angle", "value_angle_esd",
                              "period"}))
    cc.rt.torsions.push_back({row.str(0),
                              {1, row.str(1)}, {1, row.str(2)},
                              {1, row.str(3)}, {1, row.str(4)},
                              cif::as_number(row[5]), cif::as_number(row[6]),
                              cif::as_int(row[7])});
  for (auto row : block.find("_chem_comp_chir.",
                             {"atom_id_centre",
                              "atom_id_1", "atom_id_2", "atom_id_3",
                              "volume_sign"}))
    if (row[4][0] != 'c') // ignore crossN
      cc.rt.chirs.push_back({{1, row.str(0)},
                             {1, row.str(1)}, {1, row.str(2)}, {1, row.str(3)},
                             chirality_from_string(row[4])});
  for (auto row : block.find("_chem_comp_plane_atom.",
                             {"plane_id", "atom_id" , "dist_esd"})) {
    Restraints::Plane& plane = cc.rt.get_or_add_plane(row.str(0));
    if (plane.esd == 0.0)
      plane.esd = cif::as_number(row[2]);
    plane.ids.push_back({1, row.str(1)});
  }
  return cc;
}

inline ChemComp make_chemcomp_from_cif(const std::string& name,
                                       const cif::Document& doc) {
  const cif::Block* block = doc.find_block("comp_" + name);
  if (!block)
    block = doc.find_block(name);
  if (!block)
    throw std::runtime_error("data_comp_" + name + " not in the cif file");
  ChemComp cc = make_chemcomp_from_block(*block);
  if (cc.group.empty())
    if (const cif::Block* list = doc.find_block("comp_list")) {
      cif::Table tab =
        const_cast<cif::Block*>(list)->find("_chem_comp.", {"id", "group"});
      if (tab.ok())
        cc.group = tab.find_row(name).str(1);
    }
  return cc;
}

} // namespace gemmi
#endif
