// Copyright 2019 Global Phasing Ltd.
//
// Searching for links based on the _chem_link table from monomer dictionary.

#ifndef GEMMI_LINKHUNT_HPP_
#define GEMMI_LINKHUNT_HPP_

#include <map>
#include <unordered_map>
#include "elem.hpp"
#include "model.hpp"
#include "monlib.hpp"
#include "neighbor.hpp"
#include "contact.hpp"

namespace gemmi {

struct LinkHunt {
  struct Match {
    const ChemLink* chem_link = nullptr;
    int chem_link_count = 0;
    int score = -1000;
    CRA cra1;
    CRA cra2;
    bool same_image;
    float bond_length = 0.f;
    Connection* conn = nullptr;
  };

  double global_max_dist = 2.34; // ZN-CYS
  std::multimap<std::string, const ChemLink*> links;
  std::unordered_map<std::string, ChemLink::Group> res_group;

  void index_chem_links(const MonLib& monlib) {
    for (const auto& iter : monlib.links) {
      const ChemLink& link = iter.second;
      if (link.rt.bonds.empty())
        continue;
      if (link.rt.bonds.size() > 1)
        fprintf(stderr, "Note: considering only the first bond in %s\n",
                link.id.c_str());
      if (link.side1.comp.empty() && link.side2.comp.empty())
        if (link.side1.group == ChemLink::Group::Null ||
            link.side2.group == ChemLink::Group::Null ||
            link.id == "SS")
          continue;
      const Restraints::Bond& bond = link.rt.bonds[0];
      if (bond.value > global_max_dist)
        global_max_dist = bond.value;
      links.emplace(bond.lexicographic_str(), &link);
    }
    for (const auto& ri : monlib.residue_infos)
      res_group.emplace(ri.first, ChemLink::group_from_residue_info(ri.second));
  }

  bool match_link_side(const ChemLink::Side& side,
                       const std::string& resname) const {
    if (!side.comp.empty())
      return side.comp == resname;
    if (side.group == ChemLink::Group::Null)
      return false;
    auto iter = res_group.find(resname);
    return iter != res_group.end() && side.matches_group(iter->second);
  }

  std::vector<Match> find_possible_links(Structure& st,
                                         double bond_margin,
                                         double radius_margin,
                                         ContactSearch::Ignore ignore) {
    std::vector<Match> results;
    Model& model = st.first_model();
    double search_radius = std::max(global_max_dist * bond_margin,
                                    /*max r1+r2 ~=*/3.0 * radius_margin);
    NeighborSearch ns(model, st.cell, std::max(5.0, search_radius));
    ns.populate();

    ContactSearch contacts((float) search_radius);
    contacts.ignore = ignore;
    contacts.for_each_contact(ns, [&](const CRA& cra1, const CRA& cra2,
                                      int image_idx, float dist_sq) {
        Match match;

        // search for a match in chem_links
        if (bond_margin > 0) {
          auto range = links.equal_range(Restraints::lexicographic_str(
                                            cra1.atom->name, cra2.atom->name));
          for (auto iter = range.first; iter != range.second; ++iter) {
            const ChemLink& link = *iter->second;
            const Restraints::Bond& bond = link.rt.bonds[0];
            if (dist_sq > sq(bond.value * bond_margin))
              continue;
            bool order1;
            if (bond.id1.atom == cra1.atom->name &&
                match_link_side(link.side1, cra1.residue->name) &&
                match_link_side(link.side2, cra2.residue->name))
              order1 = true;
            else if (bond.id2.atom == cra1.atom->name &&
                match_link_side(link.side2, cra1.residue->name) &&
                match_link_side(link.side1, cra2.residue->name))
              order1 = false;
            else
              continue;
            int link_score = link.calculate_score(
                    order1 ? *cra1.residue : *cra2.residue,
                    order1 ? cra2.residue : cra1.residue,
                    cra1.atom->altloc ? cra1.atom->altloc : cra2.atom->altloc);
            match.chem_link_count++;
            if (link_score > match.score) {
              match.chem_link = &link;
              match.score = link_score;
              if (order1) {
                match.cra1 = cra1;
                match.cra2 = cra2;
              } else {
                match.cra1 = cra2;
                match.cra2 = cra1;
              }
            }
          }
        }

        // potential other links according to covalent radii
        if (!match.chem_link) {
          float r1 = cra1.atom->element.covalent_r();
          float r2 = cra2.atom->element.covalent_r();
          if (dist_sq > sq((r1 + r2) * radius_margin))
            return;
          match.cra1 = cra1;
          match.cra2 = cra2;
        }

        // finalize
        match.same_image = !image_idx;
        match.bond_length = std::sqrt(dist_sq);
        results.push_back(match);
    });

    // add references to st.connections
    for (Match& match : results)
      match.conn = st.find_connection_by_cra(match.cra1, match.cra2);

    return results;
  }
};

} // namespace gemmi
#endif
