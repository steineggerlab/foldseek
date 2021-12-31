// Copyright 2019 Global Phasing Ltd.
//
// Reads reflection data from the mmCIF format.

#ifndef GEMMI_REFLN_HPP_
#define GEMMI_REFLN_HPP_

#include <array>
#include "cifdoc.hpp"
#include "fail.hpp"       // for fail
#include "mmcif_impl.hpp" // for set_cell_from_mmcif, read_spacegroup_from_block
#include "numb.hpp"       // for as_number
#include "symmetry.hpp"   // for SpaceGroup
#include "unitcell.hpp"   // for UnitCell

namespace gemmi {

struct ReflnBlock {
  cif::Block block;
  std::string entry_id;
  UnitCell cell;
  const SpaceGroup* spacegroup = nullptr;
  double wavelength;
  cif::Loop* refln_loop = nullptr;
  cif::Loop* diffrn_refln_loop = nullptr;
  cif::Loop* default_loop = nullptr;

  ReflnBlock() = default;
  ReflnBlock(ReflnBlock&& rblock_) = default;
  ReflnBlock(cif::Block&& block_) : block(block_) {
    entry_id = cif::as_string(block.find_value("_entry.id"));
    impl::set_cell_from_mmcif(block, cell);
    if (const std::string* hm = impl::find_spacegroup_hm_value(block))
      spacegroup = find_spacegroup_by_name(cif::as_string(*hm),
                                           cell.alpha, cell.gamma);
    cell.set_cell_images_from_spacegroup(spacegroup);
    const char* wave_tag = "_diffrn_radiation_wavelength.wavelength";
    cif::Column wave_col = block.find_values(wave_tag);
    wavelength = wave_col.length() == 1 ? cif::as_number(wave_col[0]) : 0.;
    refln_loop = block.find_loop("_refln.index_h").get_loop();
    diffrn_refln_loop = block.find_loop("_diffrn_refln.index_h").get_loop();
    default_loop = refln_loop ? refln_loop : diffrn_refln_loop;
  }
  ReflnBlock& operator=(ReflnBlock&&) = default;
  ReflnBlock& operator=(const ReflnBlock& o) {
    block = o.block;
    entry_id = o.entry_id;
    cell = o.cell;
    spacegroup = o.spacegroup;
    wavelength = o.wavelength;
    if (o.refln_loop)
      refln_loop = block.find_loop("_refln.index_h").get_loop();
    if (o.diffrn_refln_loop)
      diffrn_refln_loop = block.find_loop("_diffrn_refln.index_h").get_loop();
    default_loop = refln_loop ? refln_loop : diffrn_refln_loop;
    return *this;
  }

  bool ok() const { return default_loop != nullptr; }
  void check_ok() const { if (!ok()) fail("Invalid ReflnBlock"); }

  // position after "_refln." or "_diffrn_refln."
  size_t tag_offset() const { return refln_loop ? 7 : 14; }

  void use_unmerged(bool unmerged) {
    default_loop = unmerged ? diffrn_refln_loop : refln_loop;
  }
  bool is_unmerged() const { return ok() && default_loop == diffrn_refln_loop; }

  std::vector<std::string> column_labels() const {
    check_ok();
    std::vector<std::string> labels(default_loop->tags.size());
    for (size_t i = 0; i != labels.size(); ++i)
      labels[i].assign(default_loop->tags[i], tag_offset(), std::string::npos);
    return labels;
  }

  int find_column_index(const std::string& tag) const {
    if (!ok())
      return -1;
    size_t name_pos = tag_offset();
    for (int i = 0; i != (int) default_loop->tags.size(); ++i)
      if (default_loop->tags[i].compare(name_pos, std::string::npos, tag) == 0)
        return i;
    return -1;
  }

  size_t get_column_index(const std::string& tag) const {
    int idx = find_column_index(tag);
    if (idx == -1)
      fail("Column not found: " + tag);
    return idx;
  }

  template<typename T>
  std::vector<T> make_vector(const std::string& tag, T null) const {
    size_t n = get_column_index(tag);
    std::vector<T> v(default_loop->length());
    for (size_t j = 0; j != v.size(); n += default_loop->width(), ++j)
      v[j] = cif::as_any(default_loop->values[n], null);
    return v;
  }

  std::array<size_t,3> get_hkl_column_indices() const {
    return {{get_column_index("index_h"),
             get_column_index("index_k"),
             get_column_index("index_l")}};
  }

  std::vector<Miller> make_miller_vector() const {
    auto hkl_idx = get_hkl_column_indices();
    std::vector<Miller> v(default_loop->length());
    for (size_t j = 0, n = 0; j != v.size(); j++, n += default_loop->width())
      for (int i = 0; i != 3; ++i)
        v[j][i] = cif::as_int(default_loop->values[n + hkl_idx[i]]);
    return v;
  }

  std::vector<double> make_1_d2_vector() const {
    if (!cell.is_crystal() || cell.a <= 0)
      fail("Unit cell is not known");
    auto hkl_idx = get_hkl_column_indices();
    std::vector<double> r(default_loop->length());
    for (size_t j = 0, n = 0; j != r.size(); j++, n += default_loop->width()) {
      Miller hkl;
      for (int i = 0; i != 3; ++i)
        hkl[i] = cif::as_int(default_loop->values[n + hkl_idx[i]]);
      r[j] = cell.calculate_1_d2(hkl);
    }
    return r;
  }

  std::vector<double> make_d_vector() const {
    std::vector<double> vec = make_1_d2_vector();
    for (double& d : vec)
      d = 1.0 / std::sqrt(d);
    return vec;
  }
};

// moves blocks from the argument to the return value
inline
std::vector<ReflnBlock> as_refln_blocks(std::vector<cif::Block>&& blocks) {
  std::vector<ReflnBlock> r;
  r.reserve(blocks.size());
  for (cif::Block& block : blocks)
    r.emplace_back(std::move(block));
  blocks.clear();
  // Some blocks miss space group tag, try to fill it in.
  const SpaceGroup* first_sg = nullptr;
  for (ReflnBlock& rblock : r)
    if (!first_sg)
      first_sg = rblock.spacegroup;
    else if (!rblock.spacegroup)
      rblock.spacegroup = first_sg;
  return r;
}

// Get the first (merged) block with required labels.
// Optionally, block name can be specified.
inline ReflnBlock get_refln_block(std::vector<cif::Block>&& blocks,
                                  const std::vector<std::string>& labels,
                                  const char* block_name=nullptr) {
  const SpaceGroup* first_sg = nullptr;
  for (cif::Block& block : blocks) {
    if (!first_sg)
      if (const std::string* hm = impl::find_spacegroup_hm_value(block)) {
        first_sg = find_spacegroup_by_name(cif::as_string(*hm));
        if (first_sg && first_sg->ext == 'H') {
          UnitCell cell;
          impl::set_cell_from_mmcif(block, cell);
          first_sg = find_spacegroup_by_name(cif::as_string(*hm),
                                             cell.alpha, cell.gamma);
        }
      }
    if (block_name && block.name != block_name)
      continue;
    if (cif::Loop* loop = block.find_loop("_refln.index_h").get_loop())
      if (std::all_of(labels.begin(), labels.end(),
            [&](const std::string& s) { return loop->has_tag("_refln."+s); })) {
        ReflnBlock rblock(std::move(block));
        if (!rblock.spacegroup && first_sg)
          rblock.spacegroup = first_sg;
        return rblock;
      }
  }
  fail("Required block or tags not found in the SF-mmCIF file.");
}

inline ReflnBlock hkl_cif_as_refln_block(cif::Block& block) {
  ReflnBlock rb;
  rb.block.swap(block);
  rb.entry_id = rb.block.name;
  impl::set_cell_from_mmcif(rb.block, rb.cell, /*mmcif=*/false);
  const char* hm_tag = "_symmetry_space_group_name_H-M";
  if (const std::string* hm = block.find_value(hm_tag))
    rb.spacegroup = find_spacegroup_by_name(cif::as_string(*hm),
                                            rb.cell.alpha, rb.cell.gamma);
  rb.cell.set_cell_images_from_spacegroup(rb.spacegroup);
  rb.refln_loop = rb.block.find_loop("_refln_index_h").get_loop();
  rb.default_loop = rb.refln_loop;
  return rb;
}

// Abstraction of data source, cf. MtzDataProxy.
struct ReflnDataProxy {
  const ReflnBlock& rb_;
  std::array<size_t,3> hkl_cols_;
  explicit ReflnDataProxy(const ReflnBlock& rb)
    : rb_(rb), hkl_cols_(rb_.get_hkl_column_indices()) {}
  size_t stride() const { return loop().tags.size(); }
  size_t size() const { return loop().values.size(); }
  using num_type = double;
  double get_num(size_t n) const { return cif::as_number(loop().values[n]); }
  const UnitCell& unit_cell() const { return rb_.cell; }
  const SpaceGroup* spacegroup() const { return rb_.spacegroup; }
  Miller get_hkl(size_t offset) const {
    return {{get_int(offset + hkl_cols_[0]),
             get_int(offset + hkl_cols_[1]),
             get_int(offset + hkl_cols_[2])}};
  }
  size_t column_index(const std::string& label) const { return rb_.get_column_index(label); }
private:
  const cif::Loop& loop() const { rb_.check_ok(); return *rb_.default_loop; }
  int get_int(size_t n) const { return cif::as_int(loop().values[n]); }
};

inline ReflnDataProxy data_proxy(const ReflnBlock& rb) { return ReflnDataProxy(rb); }

} // namespace gemmi
#endif
