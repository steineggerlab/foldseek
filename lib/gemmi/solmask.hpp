// Copyright 2021 Global Phasing Ltd.
//
// Flat bulk solvent mask. With helper tools that modify data on grid.

#ifndef GEMMI_SOLMASK_HPP_
#define GEMMI_SOLMASK_HPP_

#include "grid.hpp"      // for Grid
#include "floodfill.hpp" // for FloodFill
#include "model.hpp"     // for Model, Atom, ...

namespace gemmi {

enum class AtomicRadiiSet { VanDerWaals, Cctbx, Refmac, Constant };

// data from cctbx/eltbx/van_der_waals_radii.py used to generate identical mask
inline float cctbx_vdw_radius(El el) {
  static constexpr float radii[] = {
    /*X*/  1.00f,
    /*H*/  1.20f, /*He*/ 1.40f,
    /*Li*/ 1.82f, /*Be*/ 0.63f, /*B*/  1.75f, /*C*/ 1.775f, /*N*/ 1.50f,
    /*O*/  1.45f, /*F*/  1.47f, /*Ne*/ 1.54f,
    /*Na*/ 2.27f, /*Mg*/ 1.73f, /*Al*/ 1.50f, /*Si*/ 2.10f, /*P*/ 1.90f,
    /*S*/  1.80f, /*Cl*/ 1.75f, /*Ar*/ 1.88f,
    /*K*/  2.75f, /*Ca*/ 1.95f, /*Sc*/ 1.32f, /*Ti*/ 1.95f, /*V*/  1.06f,
    /*Cr*/ 1.13f, /*Mn*/ 1.19f, /*Fe*/ 1.26f, /*Co*/ 1.13f, /*Ni*/ 1.63f,
    /*Cu*/ 1.40f, /*Zn*/ 1.39f, /*Ga*/ 1.87f, /*Ge*/ 1.48f, /*As*/ 0.83f,
    /*Se*/ 1.90f, /*Br*/ 1.85f, /*Kr*/ 2.02f,
    /*Rb*/ 2.65f, /*Sr*/ 2.02f, /*Y*/  1.61f, /*Zr*/ 1.42f, /*Nb*/ 1.33f,
    /*Mo*/ 1.75f, /*Tc*/ 2.00f, /*Ru*/ 1.20f, /*Rh*/ 1.22f, /*Pd*/ 1.63f,
    /*Ag*/ 1.72f, /*Cd*/ 1.58f, /*In*/ 1.93f, /*Sn*/ 2.17f, /*Sb*/ 1.12f,
    /*Te*/ 1.26f, /*I*/  1.98f, /*Xe*/ 2.16f,
    /*Cs*/ 3.01f, /*Ba*/ 2.41f, /*La*/ 1.83f, /*Ce*/ 1.86f, /*Pr*/ 1.62f,
    /*Nd*/ 1.79f, /*Pm*/ 1.76f, /*Sm*/ 1.74f, /*Eu*/ 1.96f, /*Gd*/ 1.69f,
    /*Tb*/ 1.66f, /*Dy*/ 1.63f, /*Ho*/ 1.61f, /*Er*/ 1.59f, /*Tm*/ 1.57f,
    /*Yb*/ 1.54f, /*Lu*/ 1.53f, /*Hf*/ 1.40f, /*Ta*/ 1.22f, /*W*/  1.26f,
    /*Re*/ 1.30f, /*Os*/ 1.58f, /*Ir*/ 1.22f, /*Pt*/ 1.72f, /*Au*/ 1.66f,
    /*Hg*/ 1.55f, /*Tl*/ 1.96f, /*Pb*/ 2.02f, /*Bi*/ 1.73f, /*Po*/ 1.21f,
    /*At*/ 1.12f, /*Rn*/ 2.30f,
    /*Fr*/ 3.24f, /*Ra*/ 2.57f, /*Ac*/ 2.12f, /*Th*/ 1.84f, /*Pa*/ 1.60f,
    /*U*/  1.75f, /*Np*/ 1.71f, /*Pu*/ 1.67f, /*Am*/ 1.66f, /*Cm*/ 1.65f,
    /*Bk*/ 1.64f, /*Cf*/ 1.63f, /*Es*/ 1.62f, /*Fm*/ 1.61f, /*Md*/ 1.60f,
    /*No*/ 1.59f, /*Lr*/ 1.58f, /*Rf*/ 1.00f, /*Db*/ 1.00f, /*Sg*/ 1.00f,
    /*Bh*/ 1.00f, /*Hs*/ 1.00f, /*Mt*/ 1.00f, /*Ds*/ 1.00f, /*Rg*/ 1.00f,
    /*Cn*/ 1.00f, /*Nh*/ 1.00f, /*Fl*/ 1.00f, /*Mc*/ 1.00f, /*Lv*/ 1.00f,
    /*Ts*/ 1.00f, /*Og*/ 1.00f,
    /*D*/  1.20f, /*END*/0.f
  };
  static_assert(radii[static_cast<int>(El::D)] == 1.2f, "Hmm");
  static_assert(sizeof(radii) / sizeof(radii[0]) ==
                static_cast<int>(El::END) + 1, "Hmm");
  return radii[static_cast<int>(el)];
}

// Data from Refmac's ener_lib.cif: ionic radius - 0.2A or vdW radius + 0.2A.
// For full compatibility use r_probe=1.0A and r_shrink=0.8A.
inline float refmac_radius_for_bulk_solvent(El el) {
#if 0
  static constexpr float radii[] = {
    /*X*/  1.00f,
    /*H*/  1.40f, /*He*/ 1.60f,
    /*Li*/ 0.53f, /*Be*/ 0.21f, /*B*/  0.05f, /*C*/  1.90f, /*N*/  1.12f,
    /*O*/  1.08f, /*F*/  0.99f, /*Ne*/ 0.92f,
    /*Na*/ 0.93f, /*Mg*/ 0.51f, /*Al*/ 0.33f, /*Si*/ 0.20f, /*P*/  0.39f,
    /*S*/  0.20f, /*Cl*/ 1.47f, /*Ar*/ 1.34f,
    /*K*/  1.31f, /*Ca*/ 0.94f, /*Sc*/ 0.69f, /*Ti*/ 0.36f, /*V*/  0.48f,
    /*Cr*/ 0.33f, /*Mn*/ 0.26f, /*Fe*/ 0.48f, /*Co*/ 0.34f, /*Ni*/ 0.43f,
    /*Cu*/ 0.51f, /*Zn*/ 0.54f, /*Ga*/ 0.41f, /*Ge*/ 0.20f, /*As*/ 0.28f,
    /*Se*/ 0.22f, /*Br*/ 0.53f, /*Kr*/ 1.49f,
    /*Rb*/ 1.28f, /*Sr*/ 1.12f, /*Y*/  0.84f, /*Zr*/ 0.53f, /*Nb*/ 0.42f,
    /*Mo*/ 0.35f, /*Tc*/ 0.31f, /*Ru*/ 0.32f, /*Rh*/ 0.49f, /*Pd*/ 0.58f,
    /*Ag*/ 0.61f, /*Cd*/ 0.72f, /*In*/ 0.56f, /*Sn*/ 0.49f, /*Sb*/ 0.70f,
    /*Te*/ 0.37f, /*I*/  0.36f, /*Xe*/ 1.70f,
    /*Cs*/ 1.61f, /*Ba*/ 1.29f, /*La*/ 0.97f, /*Ce*/ 0.81f, /*Pr*/ 0.79f,
    /*Nd*/ 0.92f, /*Pm*/ 0.91f, /*Sm*/ 0.90f, /*Eu*/ 0.89f, /*Gd*/ 0.88f,
    /*Tb*/ 0.70f, /*Dy*/ 0.85f, /*Ho*/ 0.84f, /*Er*/ 0.83f, /*Tm*/ 0.82f,
    /*Yb*/ 0.81f, /*Lu*/ 0.80f, /*Hf*/ 0.52f, /*Ta*/ 0.58f, /*W*/  0.36f,
    /*Re*/ 0.32f, /*Os*/ 0.33f, /*Ir*/ 0.51f, /*Pt*/ 0.51f, /*Au*/ 0.51f,
    /*Hg*/ 0.90f, /*Tl*/ 0.69f, /*Pb*/ 0.59f, /*Bi*/ 0.70f, /*Po*/ 0.61f,
    /*At*/ 0.56f, /*Rn*/ 1.80f,
    /*Fr*/ 1.74f, /*Ra*/ 1.42f, /*Ac*/ 1.06f, /*Th*/ 0.88f, /*Pa*/ 0.72f,
    /*U*/  0.46f, /*Np*/ 0.65f, /*Pu*/ 0.65f, /*Am*/ 0.79f, /*Cm*/ 0.79f,
    /*Bk*/ 0.77f, /*Cf*/ 0.76f, /*Es*/ 1.00f, /*Fm*/ 1.00f, /*Md*/ 1.00f,
    /*No*/ 1.00f, /*Lr*/ 1.00f, /*Rf*/ 1.00f, /*Db*/ 1.00f, /*Sg*/ 1.00f,
    /*Bh*/ 1.00f, /*Hs*/ 1.00f, /*Mt*/ 1.00f, /*Ds*/ 1.00f, /*Rg*/ 1.00f,
    /*Cn*/ 1.00f, /*Nh*/ 1.00f, /*Fl*/ 1.00f, /*Mc*/ 1.00f, /*Lv*/ 1.00f,
    /*Ts*/ 1.00f, /*Og*/ 1.00f,
    /*D*/  1.40f, /*END*/0.f
  };
  static_assert(radii[static_cast<int>(El::D)] == 1.40f, "Hmm");
  return radii[static_cast<int>(el)];
#else
  // temporary solution used in Refmac
  switch (el) {
    case El::H: return 1.4f;
    case El::O: return 1.08f;
    case El::C: return 2.0f;
    case El::N: return 1.12f;
    default: return 1.6f;
  };
#endif
}

// mask utilities
template<typename T>
void mask_points_in_constant_radius(Grid<T>& mask, const Model& model,
                                    double radius, T value) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        mask.set_points_around(atom.pos, radius, value);
}

template<typename T>
void mask_points_in_varied_radius(Grid<T>& mask, const Model& model,
                                  AtomicRadiiSet atomic_radii_set,
                                  double r_probe, T value) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms) {
        El elem = atom.element.elem;
        double r = 0;
        switch (atomic_radii_set) {
          case AtomicRadiiSet::VanDerWaals: r = vdw_radius(elem); break;
          case AtomicRadiiSet::Cctbx: r = cctbx_vdw_radius(elem); break;
          case AtomicRadiiSet::Refmac: r = refmac_radius_for_bulk_solvent(elem); break;
          case AtomicRadiiSet::Constant: assert(0); break;
        }
        mask.set_points_around(atom.pos, r + r_probe, value);
      }
}

// All points != value in a distance < r from value are set to margin_value
template<typename T>
void set_margin_around(Grid<T>& mask, double r, T value, T margin_value) {
  int du = (int) std::floor(r / mask.spacing[0]);
  int dv = (int) std::floor(r / mask.spacing[1]);
  int dw = (int) std::floor(r / mask.spacing[2]);
  double max_spacing2 = sq(std::max(std::max(mask.unit_cell.a / mask.nu,
                                             mask.unit_cell.b / mask.nv),
                                             mask.unit_cell.c / mask.nw)) + 1e-6;
  if (2 * du >= mask.nu || 2 * dv >= mask.nv || 2 * dw >= mask.nw)
    fail("grid operation failed: radius bigger than half the unit cell?");
  std::vector<std::array<int,3>> stencil1;
  std::vector<std::array<int,3>> stencil2;
  for (int w = -dw; w <= dw; ++w)
    for (int v = -dv; v <= dv; ++v)
      for (int u = -du; u <= du; ++u) {
        Fractional fdelta = mask.get_fractional(u, v, w);
        double r2 = mask.unit_cell.orthogonalize_difference(fdelta).length_sq();
        if (r2 <= r * r && r2 != 0.) {
          std::array<int,3> wvu{{w <= 0 ? w : w - mask.nw,
                                 v <= 0 ? v : v - mask.nv,
                                 u <= 0 ? u : u - mask.nu}};
          if (r2 < max_spacing2)
            stencil1.push_back(wvu);
          else
            stencil2.push_back(wvu);
        }
      }
  if (stencil2.empty()) {
    for (const typename Grid<T>::Point& p : mask)
      if (*p.value != value) {
        for (const auto& wvu : stencil1) {
          size_t idx = mask.index_near_zero(p.u + wvu[2], p.v + wvu[1], p.w + wvu[0]);
          if (mask.data[idx] == value) {
            *p.value = margin_value;
            break;
          }
        }
      }
  } else {
    for (const typename Grid<T>::Point& p : mask) {
      if (*p.value == value) {
        bool found = false;
        for (const auto& wvu : stencil1) {
          size_t idx = mask.index_near_zero(p.u + wvu[2], p.v + wvu[1], p.w + wvu[0]);
          if (mask.data[idx] != value) {
            mask.data[idx] = margin_value;
            found = true;
          }
        }
        if (found)
          for (const auto& wvu : stencil2) {
            size_t idx = mask.index_near_zero(p.u + wvu[2], p.v + wvu[1], p.w + wvu[0]);
            if (mask.data[idx] != value)
              mask.data[idx] = margin_value;
          }
      }
    }
  }
  //printf("stencil sizes: %zu\n", stencil1.size(), stencil2.size());
  //printf("margin: %zu\n", std::count(mask.data.begin(), mask.data.end(), margin_value));
}

struct SolventMasker {
  // parameters for used only in put_solvent_mask_on_grid()
  AtomicRadiiSet atomic_radii_set = AtomicRadiiSet::VanDerWaals;
  double rprobe;
  double rshrink;
  double island_min_volume;
  double constant_r;

  SolventMasker(AtomicRadiiSet choice, double constant_r_=0.) {
    set_radii(choice, constant_r_);
  }

  void set_radii(AtomicRadiiSet choice, double constant_r_=0.) {
    atomic_radii_set = choice;
    constant_r = constant_r_;
    switch (choice) {
      case AtomicRadiiSet::VanDerWaals:
        rprobe = 1.0;
        rshrink = 1.1;
        island_min_volume = 0.;
        break;
      case AtomicRadiiSet::Cctbx:
        rprobe = 1.11;
        rshrink = 0.9;
        island_min_volume = 0.;
        break;
      case AtomicRadiiSet::Refmac:
        rprobe = 1.0;
        rshrink = 0.8;
        island_min_volume = 50;  // the exact value used in Refmac is yet to be found
        break;
      case AtomicRadiiSet::Constant:
        rprobe = 0;
        rshrink = 0;
        island_min_volume = 0.;
        break;
    }
  }

  template<typename T> void clear(Grid<T>& grid) const { grid.fill((T)1); }

  template<typename T> void mask_points(Grid<T>& grid, const Model& model) const {
    if (atomic_radii_set == AtomicRadiiSet::Constant)
      mask_points_in_constant_radius(grid, model, constant_r + rprobe, (T)0);
    else
      mask_points_in_varied_radius(grid, model, atomic_radii_set, rprobe, (T)0);
  }

  template<typename T> void symmetrize(Grid<T>& grid) const {
    grid.symmetrize([&](T a, T b) { return a == (T)0 || b == (T)0 ? (T)0 : (T)1; });
  }

  template<typename T> void shrink(Grid<T>& grid) const {
    set_margin_around(grid, rshrink, (T)1, (T)-1);
    grid.change_values((T)-1, (T)1);
  }

  template<typename T> void invert(Grid<T>& grid) const {
    for (auto& v : grid.data)
      v = (T)1 - v;
  }


  // Removes small islands of Land=1 in the sea of 0. Uses flood fill.
  // cf. find_blobs_by_flood_fill()
  template<typename T> int remove_islands(Grid<T>& grid) const {
    if (island_min_volume <= 0)
      return 0;
    size_t limit = static_cast<size_t>(island_min_volume * grid.point_count()
                                       / grid.unit_cell.volume);
    int counter = 0;
    FloodFill<T,1> flood_fill{grid};
    flood_fill.for_each_islands([&](typename FloodFill<T,1>::Result& r) {
        //printf("island %d: %zu in %zu (limit: %zu)\n",
        //       counter, r.point_count(), lines.size(), limit);
        if (r.point_count() <= limit) {
          ++counter;
          flood_fill.set_volume_values(r, (T)0);
        }
    });
    return counter;
  }

  template<typename T> void put_mask_on_grid(Grid<T>& grid, const Model& model) const {
    clear(grid);
    assert(!grid.data.empty());
    mask_points(grid, model);
    symmetrize(grid);
    shrink(grid);
    remove_islands(grid);
  }

  void set_to_zero(Grid<float>& grid, const Model& model) const {
    mask_points(grid, model);
    grid.symmetrize([&](float a, float b) { return b == 0.f ? 0.f : a; });
  }

#if 0
  template<typename T> void put_mask_on_grid(Grid<T>& grid, const Model& model) {
    // use twice finer grid for solvent mask
    Grid<std::int8_t> mask;
    mask.copy_metadata_from(grid);
    mask.set_size(2*grid.nu, 2*grid.nv, 2*grid.nw);
    mask.data.resize(8 * grid.data.size(), 1);
    put_mask_on_grid(mask, model);
    for (int w = 0, idx = 0; w < grid.nw; ++w)
      for (int v = 0; v < grid.nv; ++v)
        for (int u = 0; u < grid.nu; ++u, ++idx) {
          grid.data[idx] = 0;
          for (int wa = 0; wa < 2; ++wa)
            for (int va = 0; va < 2; ++va)
              for (int ua = 0; ua < 2; ++ua)
                grid.data[idx] += mask.get_value_q(2*u + ua, 2*v + va, 2*w + wa);
          grid.data[idx] *= 1. / 8;
        }
#endif
};

} // namespace gemmi
#endif
