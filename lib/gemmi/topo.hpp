// Copyright 2018 Global Phasing Ltd.
//
// Topo(logy) - restraints (from a monomer library) applied to a model.

#ifndef GEMMI_TOPO_HPP_
#define GEMMI_TOPO_HPP_

#include <map>           // for multimap
#include <ostream>       // for ostream
#include "chemcomp.hpp"  // for ChemComp
#include "monlib.hpp"    // for MonLib
#include "model.hpp"     // for Residue, Atom
#include "calculate.hpp" // for calculate_angle, calculate_dihedral
#include "polyheur.hpp"  // for are_connected

namespace gemmi {

struct Topo {
  // We have internal pointers in this class (pointers setup in
  // apply_restraints() that point to ResInfo::chemcomp.rt),
  // disable copying this class.
  Topo() = default;
  Topo(Topo const&) = delete;
  Topo& operator=(Topo const&) = delete;

  struct Bond {
    const Restraints::Bond* restr;
    std::array<Atom*, 2> atoms;
    double calculate() const { return atoms[0]->pos.dist(atoms[1]->pos); }
    double calculate_z() const {
      return std::abs(calculate() - restr->value) / restr->esd;
    }
  };
  struct Angle {
    const Restraints::Angle* restr;
    std::array<Atom*, 3> atoms;
    double calculate() const {
      return calculate_angle(atoms[0]->pos, atoms[1]->pos, atoms[2]->pos);
    }
    double calculate_z() const { return angle_z(calculate(), *restr); }
  };
  struct Torsion {
    const Restraints::Torsion* restr;
    std::array<Atom*, 4> atoms;
    double calculate() const {
      return calculate_dihedral(atoms[0]->pos, atoms[1]->pos,
                                atoms[2]->pos, atoms[3]->pos);
    }
    double calculate_z() const {
      return angle_z(calculate(), *restr, 360. / std::max(1, restr->period));
    }
  };
  struct Chirality {
    const Restraints::Chirality* restr;
    std::array<Atom*, 4> atoms;
    double calculate() const {
      return calculate_chiral_volume(atoms[0]->pos, atoms[1]->pos,
                                     atoms[2]->pos, atoms[3]->pos);
    }
    double calculate_z(double ideal_abs_vol, double esd) const {
      double calc = calculate();
      if (restr->sign == ChiralityType::Negative ||
          (restr->sign == ChiralityType::Both && calc < 0))
        ideal_abs_vol *= -1;
      return std::abs(calc - ideal_abs_vol) / esd;
    }
    bool check() const { return !restr->is_wrong(calculate()); }
  };
  struct Plane {
    const Restraints::Plane* restr;
    std::vector<Atom*> atoms;
    bool has(const Atom* atom) const {
      return in_vector(const_cast<Atom*>(atom), atoms);
    }
  };

  enum class RKind { Bond, Angle, Torsion, Chirality, Plane };
  struct Rule {
    RKind rkind;
    size_t index; // index in the respective vector (bonds, ...) in Topo
  };

  struct Link {
    std::string link_id;
    Residue* res1 = nullptr;
    Residue* res2 = nullptr;
    std::vector<Rule> link_rules;
    // altloc is used only for ChainInfo::extras, not for ResInfo::prev
    char alt1 = '\0';
    char alt2 = '\0';
  };

  struct ResInfo {
    Residue* res;
    // in case of microheterogeneity we may have 2+ previous residues
    std::vector<Link> prev;
    std::vector<std::string> mods;
    ChemComp chemcomp;
    std::vector<Rule> monomer_rules;

    ResInfo(Residue* r) : res(r) {}
    void add_mod(const std::string& m) {
      if (!m.empty())
        mods.push_back(m);
    }
  };

  // corresponds to a sub-chain
  struct ChainInfo {
    const Chain& chain_ref;
    std::string subchain_name;
    std::string entity_id;
    bool polymer;
    PolymerType polymer_type;
    std::vector<ResInfo> res_infos;

    ChainInfo(ResidueSpan& subchain, const Chain& chain, const Entity* ent);
    void setup_polymer_links();
    void add_refmac_builtin_modifications();
    struct RGroup {
      std::vector<ResInfo>::iterator begin, end;
    };
    RGroup group_from(std::vector<ResInfo>::iterator b) const {
      auto e = b + 1;
      while (e != res_infos.end() && e->res->group_key() == b->res->group_key())
        ++e;
      return RGroup{b, e};
    }
  };

  template<typename T>
  static int has_atom(const Atom* a, const T& t) {
    for (int i = 0; (size_t) i != t.atoms.size(); ++i)
      if (t.atoms[i] == a)
        return i;
    return -1;
  }

  std::ostream* warnings = nullptr;
  std::vector<ChainInfo> chain_infos;
  std::vector<Link> extras;

  // Restraints applied to Model
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::vector<Plane> planes;

  std::multimap<const Atom*, Bond*> bond_index;
  std::multimap<const Atom*, Angle*> angle_index;

  ResInfo* find_resinfo(const Residue* res) {
    for (ChainInfo& ci : chain_infos)
      for (ResInfo& ri : ci.res_infos)
        if (ri.res == res)
          return &ri;
    return nullptr;
  }

  Bond* first_bond_in_link(Link& link) {
    for (const Rule& rule : link.link_rules)
      if (rule.rkind == RKind::Bond)
        return &bonds[rule.index];
    return nullptr;
  }

  const Restraints::Bond* take_bond(const Atom* a, const Atom* b) const {
    auto range = bond_index.equal_range(a);
    for (auto i = range.first; i != range.second; ++i) {
      const Bond* bond = i->second;
      if ((bond->atoms[0] == b && bond->atoms[1] == a) ||
          (bond->atoms[1] == b && bond->atoms[0] == a))
        return bond->restr;
    }
    return nullptr;
  }

  const Restraints::Angle* take_angle(const Atom* a,
                                      const Atom* b,
                                      const Atom* c) const {
    auto range = angle_index.equal_range(b);
    for (auto i = range.first; i != range.second; ++i) {
      const Angle* ang = i->second;
      if ((ang->atoms[0] == a && ang->atoms[2] == c) ||
          (ang->atoms[0] == c && ang->atoms[2] == a))
        return ang->restr;
    }
    return nullptr;
  }

  const Chirality* get_chirality(const Atom* ctr) const {
    for (const Chirality& chir : chirs)
      if (chir.atoms[0] == ctr)
        return &chir;
    return nullptr;
  }

  double ideal_chiral_abs_volume(const Chirality &ch) const {
    const Restraints::Bond* bond_c1 = take_bond(ch.atoms[0], ch.atoms[1]);
    const Restraints::Bond* bond_c2 = take_bond(ch.atoms[0], ch.atoms[2]);
    const Restraints::Bond* bond_c3 = take_bond(ch.atoms[0], ch.atoms[3]);
    const Restraints::Angle* angle_1c2 = take_angle(ch.atoms[1], ch.atoms[0], ch.atoms[2]);
    const Restraints::Angle* angle_2c3 = take_angle(ch.atoms[2], ch.atoms[0], ch.atoms[3]);
    const Restraints::Angle* angle_3c1 = take_angle(ch.atoms[3], ch.atoms[0], ch.atoms[1]);
    if (bond_c1 && bond_c2 && bond_c3 && angle_1c2 && angle_2c3 && angle_3c1)
      return chiral_abs_volume(bond_c1->value, bond_c2->value, bond_c3->value,
                               angle_1c2->value, angle_2c3->value, angle_3c1->value);
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::vector<Rule> apply_restraints(const Restraints& rt,
                                     Residue& res, Residue* res2,
                                     char altloc='*') {
    std::string altlocs;
    if (altloc == '*') {
      // find all distinct altlocs
      add_distinct_altlocs(res, altlocs);
      if (res2)
        add_distinct_altlocs(*res2, altlocs);
    }
    if (altlocs.empty())
      altlocs += altloc;

    std::vector<Rule> rules;
    for (const Restraints::Bond& bond : rt.bonds)
      for (char alt : altlocs)
        if (Atom* at1 = bond.id1.get_from(res, res2, alt))
          if (Atom* at2 = bond.id2.get_from(res, res2, alt)) {
            rules.push_back({RKind::Bond, bonds.size()});
            bonds.push_back({&bond, {{at1, at2}}});
            if (!at1->altloc && !at2->altloc)
              break;
          }
    for (const Restraints::Angle& angle : rt.angles)
      for (char alt : altlocs)
        if (Atom* at1 = angle.id1.get_from(res, res2, alt))
          if (Atom* at2 = angle.id2.get_from(res, res2, alt))
            if (Atom* at3 = angle.id3.get_from(res, res2, alt)) {
              rules.push_back({RKind::Angle, angles.size()});
              angles.push_back({&angle, {{at1, at2, at3}}});
              if (!at1->altloc && !at2->altloc && !at3->altloc)
                break;
            }
    for (const Restraints::Torsion& tor : rt.torsions)
      for (char alt : altlocs)
        if (Atom* at1 = tor.id1.get_from(res, res2, alt))
          if (Atom* at2 = tor.id2.get_from(res, res2, alt))
            if (Atom* at3 = tor.id3.get_from(res, res2, alt))
              if (Atom* at4 = tor.id4.get_from(res, res2, alt)) {
                rules.push_back({RKind::Torsion, torsions.size()});
                torsions.push_back({&tor, {{at1, at2, at3, at4}}});
                if (!at1->altloc && !at2->altloc &&
                    !at3->altloc && !at4->altloc)
                  break;
          }
    for (const Restraints::Chirality& chir : rt.chirs)
      for (char alt : altlocs)
        if (Atom* at1 = chir.id_ctr.get_from(res, res2, alt))
          if (Atom* at2 = chir.id1.get_from(res, res2, alt))
            if (Atom* at3 = chir.id2.get_from(res, res2, alt))
              if (Atom* at4 = chir.id3.get_from(res, res2, alt)) {
                rules.push_back({RKind::Chirality, chirs.size()});
                chirs.push_back({&chir, {{at1, at2, at3, at4}}});
                if (!at1->altloc && !at2->altloc &&
                    !at3->altloc && !at4->altloc)
                  break;
              }
    for (const Restraints::Plane& plane : rt.planes)
      for (char alt : altlocs) {
        std::vector<Atom*> atoms;
        for (const Restraints::AtomId& id : plane.ids)
          if (Atom* atom = id.get_from(res, res2, alt))
            atoms.push_back(atom);
        if (atoms.size() >= 4) {
          rules.push_back({RKind::Plane, planes.size()});
          planes.push_back({&plane, atoms});
        }
        if (std::all_of(atoms.begin(), atoms.end(),
                        [](Atom* a) { return !a->altloc; }))
          break;
      }
    return rules;
  }

  void apply_restraints_to_residue(ResInfo& ri, const MonLib& monlib) {
    // link restraints
    for (Link& prev : ri.prev)
      if (const ChemLink* link = monlib.find_link(prev.link_id)) {
        auto rules = apply_restraints(link->rt, *prev.res1, ri.res);
        vector_move_extend(prev.link_rules, std::move(rules));
      }
    // monomer restraints
    auto rules = apply_restraints(ri.chemcomp.rt, *ri.res, nullptr);
    vector_move_extend(ri.monomer_rules, std::move(rules));
  }

  void apply_restraints_to_extra_link(Link& link, const MonLib& monlib) {
    const ChemLink* cl = monlib.find_link(link.link_id);
    if (!cl) {
      err("ignoring link '" + link.link_id + "' as it is not in the monomer library");
      return;
    }
    if (link.alt1 && link.alt2 && link.alt1 != link.alt2)
      err(tostr("LINK between different conformers ", link.alt1, " and ", link.alt2, '.'));
    char alt = link.alt1 ? link.alt1 : link.alt2;
    auto rules = apply_restraints(cl->rt, *link.res1, link.res2, alt);
    vector_move_extend(link.link_rules, std::move(rules));
  }

  // Model is non-const b/c we store non-const pointers to residues in Topo.
  // Because of the pointers, don't add or remove residues after this step.
  // Monlib may get modified by addition of extra links from the model.
  void initialize_refmac_topology(const Structure& st, Model& model0,
                                  MonLib& monlib, bool ignore_unknown_links=false);

  // This step stores pointers to gemmi::Atom's from model0,
  // so after this step don't add or remove atoms.
  // monlib is needed only for links.
  void finalize_refmac_topology(const MonLib& monlib) {
    for (ChainInfo& chain_info : chain_infos)
      for (ResInfo& ri : chain_info.res_infos)
        apply_restraints_to_residue(ri, monlib);
    for (Link& link : extras)
      apply_restraints_to_extra_link(link, monlib);

    // create indices
    for (Bond& bond : bonds) {
      bond_index.emplace(bond.atoms[0], &bond);
      if (bond.atoms[1] != bond.atoms[0])
        bond_index.emplace(bond.atoms[1], &bond);
    }
    for (Angle& ang : angles)
      angle_index.emplace(ang.atoms[1], &ang);
  }

  GEMMI_NOINLINE void err(const std::string& msg) const {
    if (warnings == nullptr)
      fail(msg);
    *warnings << "Warning: " << msg << std::endl;
  }
};

inline Topo::ChainInfo::ChainInfo(ResidueSpan& subchain,
                                  const Chain& chain, const Entity* ent)
  : chain_ref(chain) {
  subchain_name = subchain.at(0).subchain;
  res_infos.reserve(subchain.size());
  if (ent) {
    entity_id = ent->name;
    polymer = ent->entity_type == EntityType::Polymer;
    polymer_type = ent->polymer_type;
  } else {
    polymer = false;
    polymer_type = PolymerType::Unknown;
  }
  for (Residue& res : subchain)
    res_infos.emplace_back(&res);
}

inline void Topo::ChainInfo::setup_polymer_links() {
  if (!polymer || res_infos.empty())
    return;
  //ResidueGroup residue_groups()
  RGroup prev_group = group_from(res_infos.begin());
  while (prev_group.end != res_infos.end()) {
    RGroup group = group_from(prev_group.end);
    for (auto ri = group.begin; ri != group.end; ++ri)
      for (auto prev_ri = prev_group.begin; prev_ri != prev_group.end; ++prev_ri) {
        Link p;
        p.res1 = prev_ri->res;
        p.res2 = ri->res;
        assert(prev_ri - ri == p.res1 - p.res2);
        if (are_connected(*prev_ri->res, *ri->res, polymer_type)) {
          if (is_polypeptide(polymer_type)) {
            if (ri->chemcomp.group == "P-peptide")
              p.link_id = "P";  // PCIS, PTRANS
            else if (ri->chemcomp.group == "M-peptide")
              p.link_id = "NM"; // NMCIS, NMTRANS
            p.link_id += prev_ri->res->is_cis ? "CIS" : "TRANS";
          } else if (is_polynucleotide(polymer_type)) {
            p.link_id = "p";
          } else {
            p.link_id = "?";
          }
        } else {
          p.link_id = "gap";
        }
        ri->prev.push_back(p);
      }
    prev_group = group;
  }
}

inline void Topo::ChainInfo::add_refmac_builtin_modifications() {
  if (polymer && !res_infos.empty()) {
    // we try to get exactly the same numbers that makecif produces
    for (Topo::ResInfo& ri : res_infos)
      if (polymer_type == PolymerType::PeptideL)
        ri.mods.emplace_back("AA-STAND");
    Topo::ResInfo& front = res_infos.front();
    Topo::ResInfo& back = res_infos.back();
    if (is_polypeptide(polymer_type)) {
      if (front.chemcomp.group == "P-peptide")
        front.mods.emplace_back("NH2");
      else
        front.mods.emplace_back("NH3");
      back.mods.emplace_back(back.res->find_atom("OXT", '*') ? "COO" : "TERMINUS");
    } else if (is_polynucleotide(polymer_type)) {
      front.mods.emplace_back(front.res->find_atom("P", '*') ? "p5*END" : "5*END");
      back.mods.emplace_back("TERMINUS");
    }
  }
}

inline Restraints::Bond bond_restraint_from_connection(const Connection& conn) {
  Restraints::Bond bond;
  bond.id1 = Restraints::AtomId{1, conn.partner1.atom_name};
  bond.id2 = Restraints::AtomId{2, conn.partner2.atom_name};
  bond.type = BondType::Unspec;
  bond.aromatic = false;
  bond.value = conn.reported_distance;
  bond.esd = 0.02;
  return bond;
}

// Model is non-const b/c we store non-const pointers to residues in Topo.
inline void Topo::initialize_refmac_topology(const Structure& st, Model& model0,
                                             MonLib& monlib, bool ignore_unknown_links) {
  // initialize chains and residues
  for (Chain& chain : model0.chains)
    for (ResidueSpan& sub : chain.subchains()) {
      const Entity* ent = st.get_entity_of(sub);
      chain_infos.emplace_back(sub, chain, ent);
    }
  for (ChainInfo& ci : chain_infos) {
    // copy monomer description
    for (ResInfo& ri : ci.res_infos) {
      auto it = monlib.monomers.find(ri.res->name);
      if (it != monlib.monomers.end())
        ri.chemcomp = it->second;
      else
        err("unknown chemical component " + ri.res->name
            + " in chain " +  ci.chain_ref.name);
    }

    ci.setup_polymer_links();

    ci.add_refmac_builtin_modifications();

    // add modifications from standard links
    for (ResInfo& ri : ci.res_infos)
      for (Link& prev : ri.prev)
        if (const ChemLink* link = monlib.find_link(prev.link_id)) {
          ResInfo* ri_prev = &ri + (prev.res1 - prev.res2);
          ri_prev->add_mod(link->side1.mod);
          ri.add_mod(link->side2.mod);
        }
  }
  // add extra links
  for (const Connection& conn : st.connections) {
    // ignoring hydrogen bonds and metal coordination
    if (conn.type == Connection::Hydrog || conn.type == Connection::MetalC)
      continue;
    Link extra;
    extra.res1 = model0.find_cra(conn.partner1, true).residue;
    extra.res2 = model0.find_cra(conn.partner2, true).residue;
    if (!extra.res1 || !extra.res2)
      continue;
    extra.alt1 = conn.partner1.altloc;
    extra.alt2 = conn.partner2.altloc;

    // first try to find ChemLink by name (and check if it matches)
    const ChemLink* match = monlib.find_link(conn.link_id);
    if (match && (
          match->rt.bonds.empty() ||
          match->rt.bonds[0].id1.atom != conn.partner1.atom_name ||
          match->rt.bonds[0].id2.atom != conn.partner2.atom_name ||
          !monlib.link_side_matches_residue(match->side1, extra.res1->name) ||
          !monlib.link_side_matches_residue(match->side2, extra.res2->name)))
      match = nullptr;
    // if ChemLink was not found, use the first matching link (if any)
    if (!match)
      match = monlib.match_link(extra.res1->name, conn.partner1.atom_name,
                                extra.res2->name, conn.partner2.atom_name);
    if (!match) {
      match = monlib.match_link(extra.res2->name, conn.partner2.atom_name,
                                extra.res1->name, conn.partner1.atom_name);
      if (match) {
        std::swap(extra.res1, extra.res2);
        std::swap(extra.alt1, extra.alt2);
      }
    }

    if (!match && ignore_unknown_links)
      continue;

    if (match) {
      extra.link_id = match->id;
      // add modifications from the link
      find_resinfo(extra.res1)->add_mod(match->side1.mod);
      find_resinfo(extra.res2)->add_mod(match->side2.mod);
    } else {
      // create a new ChemLink and add it to the monomer library
      ChemLink cl;
      cl.side1.comp = extra.res1->name;
      cl.side2.comp = extra.res2->name;
      cl.id = cl.side1.comp + "-" + cl.side2.comp;
      cl.rt.bonds.push_back(bond_restraint_from_connection(conn));

      monlib.ensure_unique_link_name(cl.id);
      monlib.links.emplace(cl.id, cl);
      extra.link_id = cl.id;
    }
    extras.push_back(extra);
  }

  for (ChainInfo& chain_info : chain_infos)
    for (ResInfo& ri : chain_info.res_infos) {
      // apply modifications
      for (const std::string& modif : ri.mods) {
        if (const ChemMod* chem_mod = monlib.find_mod(modif))
          try {
            chem_mod->apply_to(ri.chemcomp);
          } catch(std::runtime_error& e) {
            err("failed to apply modification " + chem_mod->id
                + " to " + ri.res->name + ": " + e.what());
          }
        else
          err("modification not found: " + modif);
      }
    }
}

} // namespace gemmi
#endif
