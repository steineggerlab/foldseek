// Copyright 2018 Global Phasing Ltd.
//
// Cell-linked lists method for atom searching (a.k.a. grid search, binning,
// bucketing, cell technique for neighbor search, etc).

#ifndef GEMMI_NEIGHBOR_HPP_
#define GEMMI_NEIGHBOR_HPP_

#include <vector>
#include <cmath>  // for INFINITY, sqrt

#include "fail.hpp"      // for fail
#include "grid.hpp"
#include "model.hpp"
#include "small.hpp"

namespace gemmi {

struct NeighborSearch {

  struct Mark {
    float x, y, z;
    char altloc;
    Element element;
    int image_idx;
    int chain_idx;
    int residue_idx;
    int atom_idx;

    Mark(const Position& p, char alt, El el, int im, int ch, int res, int atom)
    : x(float(p.x)), y(float(p.y)), z(float(p.z)), altloc(alt), element(el),
      image_idx(im), chain_idx(ch), residue_idx(res), atom_idx(atom) {}

    Position pos() const { return {x, y, z}; }

    CRA to_cra(Model& mdl) const {
      Chain& c = mdl.chains.at(chain_idx);
      Residue& r = c.residues.at(residue_idx);
      Atom& a = r.atoms.at(atom_idx);
      return {&c, &r, &a};
    }
    const_CRA to_cra(const Model& mdl) const {
      const Chain& c = mdl.chains.at(chain_idx);
      const Residue& r = c.residues.at(residue_idx);
      const Atom& a = r.atoms.at(atom_idx);
      return {&c, &r, &a};
    }
    SmallStructure::Site& to_site(SmallStructure& small) const {
      return small.sites.at(atom_idx);
    }
    const SmallStructure::Site& to_site(const SmallStructure& small) const {
      return small.sites.at(atom_idx);
    }
    float dist_sq(const Position& p) const {
      return sq((float)p.x - x) + sq((float)p.y - y) + sq((float)p.z - z);
    }
  };

  Grid<std::vector<Mark>> grid;
  double radius_specified = 0.;
  Model* model = nullptr;
  SmallStructure* small_structure = nullptr;
  bool include_h = true;

  NeighborSearch() = default;
  // Model is not const so it can be modified in for_each_contact()
  NeighborSearch(Model& model_, const UnitCell& cell, double max_radius) {
    initialize(model_, cell, max_radius);
  }
  NeighborSearch(SmallStructure& small, double max_radius) {
    small_structure = &small;
    radius_specified = max_radius;
    grid.unit_cell = small.cell;
    set_grid_size();
  }
  void initialize(Model& model, const UnitCell& cell, double max_radius);
  NeighborSearch& populate(bool include_h_=true);
  void add_chain(const Chain& chain, bool include_h_=true);
  void add_chain_n(const Chain& chain, int n_ch);
  void add_atom(const Atom& atom, int n_ch, int n_res, int n_atom);
  void add_site(const SmallStructure::Site& site, int n);

  // assumes data in [0, 1), but uses index_n to handle numeric deviations
  std::vector<Mark>& get_subcell(const Fractional& fr) {
    return grid.data[grid.index_n(int(fr.x * grid.nu),
                                  int(fr.y * grid.nv),
                                  int(fr.z * grid.nw))];
  }

  template<typename Func>
  void for_each_cell(const Position& pos, const Func& func);
  template<typename Func>
  void for_each(const Position& pos, char alt, float radius, const Func& func);

  // with radius==0 it uses radius_specified
  std::vector<Mark*> find_atoms(const Position& pos, char alt, float radius) {
    if (radius == 0.f)
      radius = (float) radius_specified;
    std::vector<Mark*> out;
    for_each(pos, alt, radius, [&out](Mark& a, float) { out.push_back(&a); });
    return out;
  }

  // with max_dist==0 it uses radius_specified
  std::vector<Mark*> find_neighbors_(const Position& pos, char altloc,
                                     float min_dist, float max_dist) {
    std::vector<Mark*> out;
    if (max_dist == 0.f)
      max_dist = (float) radius_specified;
    for_each(pos, altloc, max_dist, [&](Mark& a, float dist_sq) {
        if (dist_sq >= sq(min_dist))
          out.push_back(&a);
    });
    return out;
  }
  std::vector<Mark*> find_neighbors(const Atom& atom,
                                    float min_dist, float max_dist) {
    return find_neighbors_(atom.pos, atom.altloc, min_dist, max_dist);
  }
  std::vector<Mark*> find_site_neighbors(const SmallStructure::Site& site,
                                         float min_dist, float max_dist) {
    Position pos = grid.unit_cell.orthogonalize(site.fract);
    return find_neighbors_(pos, '\0', min_dist, max_dist);
  }

  Mark* find_nearest_atom(const Position& pos) {
    Mark* mark = nullptr;
    float nearest_dist_sq = float(radius_specified * radius_specified);
    for_each_cell(pos, [&](std::vector<Mark>& marks, const Fractional& fr) {
        Position p = grid.unit_cell.orthogonalize(fr);
        for (Mark& m : marks) {
          float dist_sq = m.dist_sq(p);
          if (dist_sq < nearest_dist_sq) {
            mark = &m;
            nearest_dist_sq = dist_sq;
          }
        }
    });
    return mark;
  }

  float dist_sq(const Position& pos1, const Position& pos2) const {
    return (float) grid.unit_cell.distance_sq(pos1, pos2);
  }
  float dist(const Position& pos1, const Position& pos2) const {
    return std::sqrt(dist_sq(pos1, pos2));
  }

  FTransform get_image_transformation(int image_idx) const {
    // 0 is for identity, other indices are shifted by one.
    if (image_idx == 0)
      return Transform{};
    if ((size_t)image_idx <= grid.unit_cell.images.size())
      return grid.unit_cell.images[image_idx-1];
    fail("No such image index: " + std::to_string(image_idx));
  }

private:
  void set_grid_size() {
    grid.set_size_from_spacing(radius_specified, false);
    if (grid.nu < 3 || grid.nv < 3 || grid.nw < 3)
      grid.set_size_without_checking(std::max(grid.nu, 3),
                                     std::max(grid.nv, 3),
                                     std::max(grid.nw, 3));
  }
};


inline void NeighborSearch::initialize(Model& model_, const UnitCell& cell,
                                       double max_radius) {
  model = &model_;
  radius_specified = max_radius;
  if (cell.is_crystal()) {
    grid.unit_cell = cell;
  } else {
    Box<Position> box;
    for (const Chain& chain : model->chains)
      for (const Residue& res : chain.residues)
        for (const Atom& atom : res.atoms)
          box.extend(atom.pos);
    // We need to take into account strict NCS from MTRIXn.
    // To avoid additional function parameter that would pass Structure::ncs,
    // here we reconstruct ncs transforms from cell.images.
    // images store fractional transforms, but for non-crystal it should be
    // the same as Cartesian transform.
    std::vector<Transform> ncs;
    for (size_t n = cell.cs_count; n < cell.images.size(); n += cell.cs_count + 1)
      ncs.push_back(cell.images[n]);
    // The box needs include all NCS images as well.
    if (!ncs.empty()) {
      for (CRA cra : model->all())
        for (const Transform& tr : ncs)
          box.extend(Position(tr.apply(cra.atom->pos)));
    }
    box.add_margin(1.5 * max_radius);  // much more than needed
    Position size = box.get_size();
    grid.unit_cell.set(size.x, size.y, size.z, 90, 90, 90);
    for (const Transform& tr : ncs) {
      UnitCell& c = grid.unit_cell;
      // cf. add_ncs_images_to_cs_images()
      c.images.push_back(c.frac.combine(tr.combine(c.orth)));
    }
  }
  set_grid_size();
}

inline NeighborSearch& NeighborSearch::populate(bool include_h_) {
  include_h = include_h_;
  if (model) {
    for (int n_ch = 0; n_ch != (int) model->chains.size(); ++n_ch)
      add_chain_n(model->chains[n_ch], n_ch);
  } else if (small_structure) {
    for (int n = 0; n != (int) small_structure->sites.size(); ++n) {
      SmallStructure::Site& site = small_structure->sites[n];
      if (include_h || !site.element.is_hydrogen())
        add_site(site, n);
    }
  } else {
    fail("NeighborSearch not initialized");
  }
  return *this;
}

inline void NeighborSearch::add_chain(const Chain& chain, bool include_h_) {
  if (!model)
    fail("NeighborSearch.add_chain(): model not initialized yet");
  // to be safe avoid (&chain - model.chains[0]) which could be UB
  for (int n_ch = 0; n_ch != (int) model->chains.size(); ++n_ch)
    if (&model->chains[n_ch] == &chain) {
      include_h = include_h_;
      add_chain_n(chain, n_ch);
      return;
    }
  fail("NeighborSearch.add_chain(): chain not in this model");
}

inline void NeighborSearch::add_chain_n(const Chain& chain, int n_ch) {
  for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
    const Residue& res = chain.residues[n_res];
    for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
      const Atom& atom = res.atoms[n_atom];
      if (include_h || !atom.is_hydrogen())
        add_atom(atom, n_ch, n_res, n_atom);
    }
  }
}

inline void NeighborSearch::add_atom(const Atom& atom,
                                     int n_ch, int n_res, int n_atom) {
  const UnitCell& gcell = grid.unit_cell;
  Fractional frac0 = gcell.fractionalize(atom.pos);
  {
    Fractional frac = frac0.wrap_to_unit();
    Position pos = gcell.orthogonalize(frac);
    get_subcell(frac).emplace_back(pos, atom.altloc, atom.element.elem,
                                   0, n_ch, n_res, n_atom);
  }
  for (int n_im = 0; n_im != (int) gcell.images.size(); ++n_im) {
    Fractional frac = gcell.images[n_im].apply(frac0).wrap_to_unit();
    Position pos = gcell.orthogonalize(frac);
    get_subcell(frac).emplace_back(pos, atom.altloc, atom.element.elem,
                                   n_im + 1, n_ch, n_res, n_atom);
  }
}

// We exclude special position images of atoms here, but not in add_atom.
// This choice is somewhat arbitrary, but it also reflects the fact that
// in MX files occupances of atoms on special positions are (almost always)
// fractional and all images are to be taken into account.
inline void NeighborSearch::add_site(const SmallStructure::Site& site, int n) {
  const double SPECIAL_POS_TOL = 0.4;
  const UnitCell& gcell = grid.unit_cell;
  std::vector<Fractional> others;
  others.reserve(gcell.images.size());
  Fractional frac0 = site.fract.wrap_to_unit();
  {
    Position pos = gcell.orthogonalize(frac0);
    get_subcell(frac0).emplace_back(pos, '\0', site.element.elem, 0, -1, -1, n);
  }
  for (int n_im = 0; n_im != (int) gcell.images.size(); ++n_im) {
    Fractional frac = gcell.images[n_im].apply(site.fract).wrap_to_unit();
    if (gcell.distance_sq(frac, frac0) < sq(SPECIAL_POS_TOL) ||
        std::any_of(others.begin(), others.end(), [&](const Fractional& f) {
          return gcell.distance_sq(frac, f) < sq(SPECIAL_POS_TOL);
        }))
      continue;
    Position pos = gcell.orthogonalize(frac);
    get_subcell(frac).emplace_back(pos, '\0', site.element.elem,
                                   n_im + 1, -1, -1, n);
    others.push_back(frac);
  }
}

template<typename Func>
void NeighborSearch::for_each_cell(const Position& pos, const Func& func) {
  Fractional fr = grid.unit_cell.fractionalize(pos).wrap_to_unit();
  const int u0 = int(fr.x * grid.nu);
  const int v0 = int(fr.y * grid.nv);
  const int w0 = int(fr.z * grid.nw);
  const int uend = u0 + std::min(3, grid.nu) - 1;
  const int vend = v0 + std::min(3, grid.nv) - 1;
  const int wend = w0 + std::min(3, grid.nw) - 1;
  for (int w = w0 - 1; w < wend; ++w) {
    int dw = w >= grid.nw ? -1 : w < 0 ? 1 : 0;
    for (int v = v0 - 1; v < vend; ++v) {
      int dv = v >= grid.nv ? -1 : v < 0 ? 1 : 0;
      for (int u = u0 - 1; u < uend; ++u) {
        int du = u >= grid.nu ? -1 : u < 0 ? 1 : 0;
        size_t idx = grid.index_q(u + du * grid.nu,
                                  v + dv * grid.nv,
                                  w + dw * grid.nw);
        func(grid.data[idx], Fractional(fr.x + du, fr.y + dv, fr.z + dw));
      }
    }
  }
}

template<typename Func>
void NeighborSearch::for_each(const Position& pos, char alt, float radius,
                              const Func& func) {
  if (radius <= 0.f)
    return;
  for_each_cell(pos, [&](std::vector<Mark>& marks, const Fractional& fr) {
      Position p = grid.unit_cell.orthogonalize(fr);
      for (Mark& m : marks) {
        float dist_sq = m.dist_sq(p);
        if (m.dist_sq(p) < sq(radius) && is_same_conformer(alt, m.altloc))
          func(m, dist_sq);
      }
  });
}


inline void remove_cras(Model& model, std::vector<CRA>& vec) {
  // sort in reverse order, so items can be erased without invalidating pointers
  std::sort(vec.begin(), vec.end(), [](const CRA& a, const CRA& b) {
      return std::tie(a.chain, a.residue, a.atom) > std::tie(b.chain, b.residue, b.atom);
  });
  const Atom* prev_a = nullptr;
  for (CRA& cra : vec) {
    if (cra.atom == prev_a)
      continue;
    prev_a = cra.atom;
    auto atom_idx = cra.atom - cra.residue->atoms.data();
    cra.residue->atoms.erase(cra.residue->atoms.begin() + atom_idx);
    if (cra.residue->atoms.empty()) {
      auto res_idx = cra.residue - cra.chain->residues.data();
      cra.chain->residues.erase(cra.chain->residues.begin() + res_idx);
      if (cra.chain->residues.empty()) {
        auto chain_idx = cra.chain - model.chains.data();
        model.chains.erase(model.chains.begin() + chain_idx);
      }
    }
  }
}

// To be used after expand_ncs() and make_assembly().
// Searches and merges overlapping equivalent atoms from different chains.
inline void merge_atoms_in_expanded_model(Model& model, const UnitCell& cell,
                                          double max_dist=0.2) {
  using Mark = NeighborSearch::Mark;
  NeighborSearch ns(model, cell, 4.0);
  ns.populate(true);
  std::vector<CRA> to_be_deleted;
  for (int n_ch = 0; n_ch != (int) model.chains.size(); ++n_ch) {
    Chain& chain = model.chains[n_ch];
    for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
      Residue& res = chain.residues[n_res];
      for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
        Atom& atom = res.atoms[n_atom];
        std::vector<CRA> equiv;
        ns.for_each_cell(atom.pos, [&](std::vector<Mark>& marks, const Fractional& fr) {
            for (Mark& m : marks) {
              // We look for the same atoms, but copied to a different chain.
              // First quick check that filters out most of non-matching pairs.
              if (m.altloc != atom.altloc || m.element != atom.element ||
                  m.chain_idx == n_ch || m.atom_idx != n_atom)
                continue;
              // Now check if everything else matches.
              CRA cra = m.to_cra(model);
              if (cra.atom &&
                  cra.atom->serial == atom.serial &&
                  cra.atom->name == atom.name &&
                  cra.atom->b_iso == atom.b_iso &&
                  cra.residue->matches_noseg(res) &&
                  m.dist_sq(ns.grid.unit_cell.orthogonalize(fr)) < sq(max_dist))
              equiv.push_back(cra);
            }
        });
        if (!equiv.empty()) {
          // calculate new occupancy and position for the first atom
          int n = 1;
          if (cell.is_crystal())
            n += cell.is_special_position(atom.pos, max_dist);
          double occ_sum = n * atom.occ;
          if (occ_sum > 0) {
            atom.pos *= occ_sum;
            for (CRA& cra : equiv) {
              atom.pos += cra.atom->occ * cra.atom->pos;
              occ_sum += cra.atom->occ;
            }
            atom.pos /= occ_sum;
            atom.occ = float(occ_sum / n);
          }
          // discard overlapping equivalent atoms
          for (CRA& cra : equiv) {
            // change atoms to avoid processing them again
            cra.atom->serial = -1;
            cra.atom->name.clear();
            // deleting now would invalidate indices in NeighborSearch
            to_be_deleted.push_back(cra);
          }
        }
      }
    }
  }
  remove_cras(model, to_be_deleted);
}


// simple map alignment - this function may change in the future
template<typename T>
void interpolate_grid_of_aligned_model(Grid<T>& dest, const Grid<T>& src,
                                       const Transform& tr, NeighborSearch& ns,
                                       double radius=0.) {
  if (radius <= 0.)
    radius = ns.radius_specified;
  else if (radius > ns.radius_specified)
    fail("set_grid_values_interpolated_from(): radius exceeds NeighborSearch radius");
  // mask model or its part that was used for NeighborSearch
  std::vector<bool> mask(dest.data.size(), false);
  for (const std::vector<NeighborSearch::Mark>& marks : ns.grid.data)
    for (const NeighborSearch::Mark& mark : marks)
      if (mark.image_idx == 0)  // leave out symmetry mates
        dest.template use_points_around<true>(
            ns.grid.unit_cell.fractionalize(mark.pos()),
            radius,
            [&](T& point, double) { mask[&point - dest.data.data()] = true; });
  size_t idx = 0;
  for (int w = 0; w != dest.nw; ++w)
    for (int v = 0; v != dest.nv; ++v)
      for (int u = 0; u != dest.nu; ++u, ++idx) {
        if (!mask[idx])
          continue;
        Position pos2 = dest.get_position(u, v, w);
        // in contacts use only nodes that are nearer to the original molecule
        NeighborSearch::Mark* mark = ns.find_nearest_atom(pos2);
        if (mark && mark->image_idx == 0) {
          Position delta = mark->to_cra(*ns.model).atom->pos - mark->pos();
          Position pos1 = Position(tr.apply(pos2 + delta));
          dest.data[idx] = src.interpolate_value(pos1);
        }
      }
}

} // namespace gemmi
#endif
