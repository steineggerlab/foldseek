// Copyright 2020 Global Phasing Ltd.
//
// A class for converting MTZ (merged or unmerged) to SF-mmCIF

// TODO:
//  - cell parameters may be different in CELL and DCELL records, check for it
//  - check that the FP column is not from Refmac
//  - should we allow for repeated column names in MTZ?

#ifndef GEMMI_MTZ2CIF_HPP_
#define GEMMI_MTZ2CIF_HPP_

#include <ostream>
#include <unordered_map>
#include <climits>       // for INT_MIN
#include <set>
#include <algorithm>     // for all_of
#include "asudata.hpp"   // for calculate_hkl_value_correlation
#include "eig3.hpp"      // for eigen_decomposition
#include "mtz.hpp"       // for Mtz
#include "xds_ascii.hpp" // for XdsAscii
#include "atox.hpp"      // for read_word
#include "merge.hpp"     // for Intensities, read_unmerged_intensities_from_mtz
#include "sprintf.hpp"   // for gf_snprintf, to_str
#include "version.hpp"   // for GEMMI_VERSION

namespace gemmi {

class MtzToCif {
public:
  enum Var { Dot=-1, Qmark=-2, Counter=-3, DatasetId=-4, Image=-5 };

  // options that can be set directly
  std::vector<std::string> spec_lines; // conversion specification (cf. default_spec)
  const char* block_name = nullptr;  // NAME in data_NAME
  std::string entry_id = "xxxx";     // _entry.id
  bool with_comments = true;         // write comments
  bool with_history = true;          // write MTZ history in comments
  bool skip_empty = false;           // skip reflections with no values
  bool enable_UB = false;            // write _diffrn_orient_matrix.UB
  bool write_staraniso_tensor = true; // write _reflns.pdbx_aniso_B_tensor_*
  bool write_special_marker_for_pdb = false;
  int less_anomalous = 0;            // skip (+)/(-) columns even if in spec
  std::string skip_empty_cols;       // columns used to determine "emptiness"
  double wavelength = NAN;           // user-specified wavelength
  int trim = 0;                      // output only reflections -N<=h,k,l<=N
  int free_flag_value = -1;          // -1 = auto: 0 or (if we have >50% of 0's) 1
  std::string staraniso_version;     // for _software.version in "special_marker"

  static const char** default_spec(bool for_merged) {
    static const char* merged[] = {
      "H                          H index_h",
      "K                          H index_k",
      "L                          H index_l",
      "? IMEAN|I|IOBS|I-obs       J intensity_meas",
      "& SIG{prev}                Q intensity_sigma",
      "? I(+)|IOBS(+)|I-obs(+)    K pdbx_I_plus",
      "& SIG{prev}                M pdbx_I_plus_sigma",
      "? I(-)|IOBS(-)|I-obs(-)    K pdbx_I_minus",
      "& SIG{prev}                M pdbx_I_minus_sigma",
      // TODO: FP from Refmac should show warning or error
      "? F|FP|FOBS|F-obs          F F_meas_au",
      "& SIG{prev}                Q F_meas_sigma_au",
      "? F(+)|FOBS(+)|F-obs(+)    G pdbx_F_plus",
      "& SIG{prev}                L pdbx_F_plus_sigma",
      "? F(-)|FOBS(-)|F-obs(-)    G pdbx_F_minus",
      "& SIG{prev}                L pdbx_F_minus_sigma",
      "? FREE|RFREE|FREER|FreeR_flag|R-free-flags I status S",
      "? FWT|2FOFCWT              F pdbx_FWT",
      "& PHWT|PH2FOFCWT           P pdbx_PHWT .3f",
      "? DELFWT|FOFCWT            F pdbx_DELFWT",
      "& DELPHWT|PHDELWT|PHFOFCWT P pdbx_DELPHWT .3f",
      nullptr
    };
    static const char* unmerged[] = {
      "$dataset    diffrn_id",  // diffrn_id - sweep id deduced from BATCH
      "$counter    id",         // reflection counter (1, 2, ...)
      "H         H index_h",
      "K         H index_k",
      "L         H index_l",
      "? I       J intensity_net",
      "& SIGI    Q intensity_sigma .5g",
      // new! https://github.com/wwpdb-dictionaries/mmcif_pdbx/pull/33
      "?ROT      R pdbx_scan_angle",
      "$image      pdbx_image_id",
      nullptr
    };
    return for_merged ? merged : unmerged;
  }

  void write_cif(const Mtz& mtz, const Mtz* mtz2,
                 SMat33<double>* staraniso_b, std::ostream& os);
  void write_cif_from_xds(const XdsAscii& xds, std::ostream& os);

private:
  // describes which MTZ column is to be translated to what mmCIF column
  struct Trans {
    int col_idx;
    bool is_status = false;
    std::string tag;  // excluding category
    std::string format = "%g";
    int min_width = 0;
  };

  // data corresponding to one sweep (dataset) in unmerged MTZ file
  struct SweepData {
    int id = 0;
    int batch_count = 0;
    int offset = 0;
    int crystal_id = 1;
    const Mtz::Batch* first_batch = nullptr;
    const Mtz::Dataset* dataset = nullptr;
  };

  std::vector<SweepData> sweeps;
  std::unordered_map<int, int> sweep_indices;

  std::vector<Trans> recipe;

  const Trans* get_status_translation() const {
    for (const Trans& t: recipe)
      if (t.is_status)
        return &t;
    return nullptr;
  }

  bool gather_sweep_data(const Mtz& mtz) {
    bool ok = true;
    int prev_number = INT_MIN;
    int prev_dataset = INT_MIN;
    SweepData sweep;
    for (const Mtz::Batch& batch : mtz.batches) {
      int dataset_id = batch.dataset_id();
      if (dataset_id != prev_dataset || batch.number != prev_number + 1) {
        prev_dataset = dataset_id;
        if (sweep.id != 0)
          sweeps.push_back(sweep);
        sweep.id++;
        sweep.batch_count = 0;
        sweep.first_batch = &batch;
        sweep.offset = batch.number - (batch.number % 100);
        try {
          sweep.dataset = &mtz.dataset(dataset_id);
        } catch (std::runtime_error&) {
          sweep.dataset = nullptr;
          ok = false;
          mtz.warn("Reference to absent dataset: " + std::to_string(dataset_id));
        }
      }
      sweep_indices.emplace(batch.number, sweep.id - 1);
      sweep.batch_count++;
      prev_number = batch.number;
    }
    if (sweep.id != 0)
      sweeps.push_back(sweep);
    return ok;
  }

  // Get the first (non-zero) DWAVEL corresponding to an intensity, amplitude
  // or sigma column from the template.
  static double get_wavelength(const Mtz& mtz, const std::vector<Trans>& spec) {
    for (const Trans& tr : spec) {
      if (tr.col_idx < 0)
        continue;
      const Mtz::Column& col = mtz.columns.at(tr.col_idx);
      if (col.type == 'F' || col.type == 'J' ||
          col.type == 'K' || col.type == 'G' ||
          col.type == 'Q' || col.type == 'M' || col.type == 'L') { // sigma
        double wavelength = mtz.dataset(col.dataset_id).wavelength;
        if (wavelength != 0.)
          return wavelength;
      }
    }
    return 0.;
  }

  static int find_column_index(const std::string& column, const Mtz& mtz) {
    int idx = -1;
    for (const std::string& label : split_str(column, '|')) {
      for (size_t i = 0; i != mtz.columns.size(); ++i) {
        if (mtz.columns[i].label == label) {
          if (idx == -1)
            idx = (int) i;
          else
            mtz.warn("Column label duplicated: " + label);
        }
      }
      if (idx != -1)
        break;
    }
    return idx;
  }

  static int check_format(const std::string& fmt) {
    // expected format: [#_+-]?\d*(\.\d+)?[fFgGeEc]
    int min_width = 0;
    if (fmt.find('%') != std::string::npos)
      fail("Specify format without %. Got: " + fmt);
    const char* p = fmt.c_str();
    if (*p == '_' || *p == '+' || *p == '-' || *p == '#')
     ++p;
    if (is_digit(*p)) {
      min_width = *p++ - '0';
      if (is_digit(*p)) // two digits of width number max
        min_width = min_width * 10 + (*p++ - '0');
    }
    if (*p == '.' && is_digit(*(p+1))) {
      p += 2;
      if (is_digit(*p)) // two digits of precision numbers max
        ++p;
    }
    if (!std::isalpha(*p) || *(p+1) != '\0')
      fail("wrong format : " + fmt + "\nCorrect examples: g, .4f, 12.5e");
    char c = alpha_up(*p);
    if (c != 'F' && c != 'G' && c != 'E')
      fail("expected floating-point format, got: " + fmt);
    if (min_width > 32)
      fail("the width exceeds 32: " + fmt);
    return min_width;
  }

  // state for parse_spec_line
  struct SpecParserState {
    size_t verified_spec_size = 0;
    bool discard_next_line = false;
  };

  void prepare_recipe(const Mtz& mtz) {
    recipe.clear();
    SpecParserState state;
    if (!spec_lines.empty()) {
      for (const std::string& line : spec_lines)
        parse_spec_line(line.c_str(), mtz, state);
    } else {
      const char** lines = default_spec(/*for_merged=*/mtz.batches.empty());
      for (; *lines != nullptr; ++lines)
        parse_spec_line(*lines, mtz, state);
    }
    if (recipe.empty())
      fail("empty translation recipe");
    for (size_t i = 0; i != recipe.size(); ++i)
      for (size_t j = i + 1; j != recipe.size(); ++j)
        if (recipe[i].tag == recipe[j].tag)
          fail("duplicated output tag: " + recipe[i].tag);
    // H, K, L must be the first columns in MTZ and are required in _refln
    for (int i = 2; i != -1; --i)
      if (!in_vector_f([&](const Trans& t) { return t.col_idx == i; }, recipe)) {
        Trans tr;
        tr.col_idx = i;
        tr.tag = "index_";
        tr.tag += "hkl"[i]; // h, k or l
        recipe.insert(recipe.begin(), tr);
      }
  }

  // adds results to recipe
  void parse_spec_line(const char* line, const Mtz& mtz, SpecParserState& state) {
    Trans tr;
    const char* p = line;
    if (*p == '&') {
      if (state.discard_next_line)
        return;
    } else {
      state.verified_spec_size = recipe.size();
      state.discard_next_line = false;
    }
    bool optional = (*p == '?' || *p == '&');
    if (optional)
      ++p;
    std::string column = read_word(p, &p);

    char ctype = '\0';
    if (column[0] != '$') {
      std::string type = read_word(p, &p);
      if (type.size() != 1)
        fail("Spec error: MTZ type '" + type + "' is not one character,"
             "\nin line: " + line);
      ctype = type[0];
      if (less_anomalous > 0) {
        bool is_i_ano = (ctype == 'K' || ctype == 'M');
        if (less_anomalous == 1
            ? is_i_ano && mtz.count_type('J') != 0 && mtz.count_type('G') > 1 &&
                          staraniso_version.empty()
            : is_i_ano || ctype == 'G' || ctype == 'D' || ctype == 'L') {
          recipe.resize(state.verified_spec_size);
          state.discard_next_line = true;
          return;
        }
      }
    }

    tr.tag = read_word(p, &p);
    if (tr.tag[0] == '_' || tr.tag.find('.') != std::string::npos)
      fail("Spec error: expected tag part after _refln., got: " +
           tr.tag + "\nin line: " + line);

    if (column[0] == '$') {
      if (column.size() == 2 && column[1] == '.')
        tr.col_idx = Var::Dot;
      else if (column.size() == 2 && column[1] == '?')
        tr.col_idx = Var::Qmark;
      else if (mtz.is_merged())
        fail("Variables other than $. and $? can be used only for unmerged files");
      else if (column.compare(1, 7, "counter") == 0)
        tr.col_idx = Var::Counter;
      else if (column.compare(1, 7, "dataset") == 0)
        tr.col_idx = Var::DatasetId;
      else if (column.compare(1, 5, "image") == 0)
        tr.col_idx = Var::Image;
      else
        fail("Unknown variable in the spec: " + column);
    } else {
      if (column.find('{') != std::string::npos &&
          !recipe.empty() && recipe.back().col_idx >= 0)
        replace_all(column, "{prev}", mtz.columns[recipe.back().col_idx].label);
      tr.col_idx = find_column_index(column, mtz);
      if (tr.col_idx == -1) {
        if (!optional)
          fail("Column not found: " + column);
        recipe.resize(state.verified_spec_size);
        state.discard_next_line = true;
        return;
      }
      const Mtz::Column& col = mtz.columns[tr.col_idx];
      if (ctype != '*' && col.type != ctype)
        fail("Column ", col.label, " has type ", col.type, " not ", ctype);
    }

    std::string fmt = read_word(p, &p);
    if (!fmt.empty()) {
      if (fmt.size() == 1 && fmt[0] == 'S') {
        tr.is_status = true;
      } else {
        tr.min_width = check_format(fmt);
        tr.format = "%" + fmt;
        if (tr.format[1] == '_')
          tr.format[1] = ' ';
      }
    }

    recipe.push_back(tr);
  }

  void write_special_marker_if_requested(std::ostream& os, bool merged) const {
    if (!write_special_marker_for_pdb)
      return;
    os << "### IF YOU MODIFY THIS FILE, REMOVE THIS SIGNATURE: ###\n";
    if (!merged || staraniso_version.empty()) {
      os << "_software.pdbx_ordinal 1\n"
            "_software.classification 'data extraction'\n"
            "_software.name gemmi\n"
            "_software.version " GEMMI_VERSION "\n";
    } else {
      os << "loop_\n"
            "_software.pdbx_ordinal\n"
            "_software.classification\n"
            "_software.name\n"
            "_software.version\n"
            "1 'data extraction' gemmi " GEMMI_VERSION "\n";
      // STARANISO here tells that intensities were scaled anisotropically.
      os << "2 'data scaling' STARANISO '" << staraniso_version << "'\n";
    }
    os << "_pdbx_audit_conform.dict_name mmcif_pdbx.dic\n"
          "_pdbx_audit_conform.dict_version 5.339\n"
          "_pdbx_audit_conform.dict_location "
          "https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic\n"
          "### END OF SIGNATURE ###\n\n";
  }

  void write_cell_and_symmetry(const UnitCell& cell, double* rmsds,
                               const SpaceGroup* sg,
                               char* buf, std::ostream& os) const;

  void write_main_loop(const Mtz& mtz, char* buf, std::ostream& os);
};

inline bool validate_merged_mtz_deposition_columns(const Mtz& mtz, std::ostream& out) {
  bool ok = true;
  if (!mtz.rfree_column()) {
    out << "ERROR. Merged file is missing free-set flag.\n";
    ok = false;
  }
  if (!mtz.imean_column() && !mtz.iplus_column()) {
    out << "ERROR. Merged file is missing intensities.\n";
    ok = false;
  }
  if (!mtz.column_with_one_of_labels({"F", "FP", "FOBS", "F-obs",
                                      "F(+)", "FOBS(+)", "F-obs(+)"})) {
    out << "Merged file is missing amplitudes\n"
           "(which is fine if intensities were used for refinement)\n";
  }
  if (!ok) {
    out << "Columns in the merged file:";
    for (const Mtz::Column& col : mtz.columns)
      out << ' ' << col.label;
    out << '\n';
  }
  return ok;
}

// note: both mi and ui get modified
inline bool validate_merged_intensities(Intensities& mi, Intensities& ui,
                                        bool relaxed_check, std::ostream& out) {
  // XDS files have 4 significant digits. Using accuracy 5x the precision.
  const double max_diff = 0.005;
  out << "Checking if both files match...\n";
  bool ok = true;
  if (ui.spacegroup == mi.spacegroup) {
    out << "The same space group: " << mi.spacegroup_str() << '\n';
  } else {
    GroupOps gops1 = ui.spacegroup->operations();
    GroupOps gops2 = mi.spacegroup->operations();
    if (!gops1.has_same_centring(gops2) || !gops1.has_same_rotations(gops2))
      ok = false;
    out << (ok ? "WARNING" : "ERROR")
        << ". Different space groups in merged and unmerged files:\n"
        << mi.spacegroup_str() << " and " << ui.spacegroup_str() << '\n';
    if (!ok)
      out << "(in the future, this app may recognize compatible space groups\n"
             "and reindex unmerged data if needed; for now, it's on you)\n";
  }

  auto eq = [](double x, double y, double rmsd) { return std::fabs(x - y) < rmsd + 0.02; };
  if(eq(mi.unit_cell.a,     ui.unit_cell.a,     ui.unit_cell_rmsd[0]) &&
     eq(mi.unit_cell.b,     ui.unit_cell.b,     ui.unit_cell_rmsd[1]) &&
     eq(mi.unit_cell.c,     ui.unit_cell.c,     ui.unit_cell_rmsd[2]) &&
     eq(mi.unit_cell.alpha, ui.unit_cell.alpha, ui.unit_cell_rmsd[3]) &&
     eq(mi.unit_cell.beta,  ui.unit_cell.beta,  ui.unit_cell_rmsd[4]) &&
     eq(mi.unit_cell.gamma, ui.unit_cell.gamma, ui.unit_cell_rmsd[5])) {
    out << "The same unit cell parameters.\n";
  } else {
    const UnitCell& mc = mi.unit_cell;
    const UnitCell& uc = ui.unit_cell;
    out << "Unit cell parameters differ:";
    out << "\n    merged: " << mc.a << ' ' << mc.b << ' ' << mc.c << "  "
                            << mc.alpha << ' ' << mc.beta << ' ' << mc.gamma;
    out << "\n  unmerged: " << uc.a << ' ' << uc.b << ' ' << uc.c << "  "
                            << uc.alpha << ' ' << uc.beta << ' ' << uc.gamma;
    out << '\n';
    ok = false;
  }

  size_t ui_size1 = ui.data.size();
  ui.merge_in_place(mi.type);  // it also sorts
  size_t ui_size2 = ui.data.size();
  ui.remove_systematic_absences();
  out << "Unmerged reflections: " << ui_size1 << " (" << ui_size2 << " merged "
      << mi.type_str() << ", " << ui.data.size() << " w/o sysabs)\n";
  mi.switch_to_asu_indices(/*merged=*/true);
  mi.sort();
  size_t mi_size1 = mi.data.size();
  mi.remove_systematic_absences();
  out << "Merged reflections: " << mi_size1 << ' ' << mi.type_str()
      << " (" << mi.data.size() << " w/o sysabs)\n";

  if (mi.staraniso_b.ok()) {
    out << "Taking into account the anisotropy tensor that was used for scaling.\n";
    for (Intensities::Refl& refl : ui.data)
      refl.value *= mi.staraniso_b.scale(refl.hkl, ui.unit_cell);
  }

  // first pass - calculate CC and scale
  gemmi::Correlation corr = calculate_hkl_value_correlation(ui.data, mi.data);
  out << "Corr. coef. of " << corr.n << ' ' << mi.type_str() << " values: "
      << 100 * corr.coefficient() << "%\n";
  double scale = corr.mean_ratio();
  out << "Ratio of compared intensities (merged : unmerged): " << scale << '\n';

  // second pass - check that all reflections match
  double max_weighted_sq_diff = 0.;
  const Intensities::Refl* max_diff_r1 = nullptr;
  const Intensities::Refl* max_diff_r2 = nullptr;
  int differ_count = 0;
  int missing_count = 0;
  auto r1 = ui.data.begin();
  auto r2 = mi.data.begin();
  auto refln_str = [](const Intensities::Refl& r) {
    return r.intensity_label() + std::string(" ") + miller_str(r.hkl);
  };
  while (r1 != ui.data.end() && r2 != mi.data.end()) {
    if (r1->hkl == r2->hkl && r1->isign == r2->isign) {
      if (!relaxed_check) {
        double value1 = scale * r1->value;
        double sigma1 = scale * r1->sigma; // is this approximately correct
        double sq_max = std::max(sq(value1), sq(r2->value));
        double sq_diff = sq(value1 - r2->value);
        // Intensities may happen to be rounded to two decimal places,
        // so if the absolute difference is <0.01 it's OK.
        if (sq_diff > 1e-4 && sq_diff > sq(max_diff) * sq_max) {
          if (differ_count == 0) {
            out << "First difference: " << r1->hkl_label()
                << ' ' << value1 << " vs " << r2->value << '\n';
          }
          ++differ_count;
          double weighted_sq_diff = sq_diff / (sq(sigma1) + sq(r2->sigma));
          if (weighted_sq_diff > max_weighted_sq_diff) {
            max_weighted_sq_diff = weighted_sq_diff;
            max_diff_r1 = &*r1;
            max_diff_r2 = &*r2;
          }
        }
      }
      ++r1;
      ++r2;
    } else if (*r1 < *r2) {
      ++r1;
    } else {
      if (missing_count == 0)
        out << "First missing reflection in unmerged data: " << refln_str(*r1) << '\n';
      ++missing_count;
      ++r2;
    }
  }

  if (differ_count != 0) {
    out << "Most significant difference: " << refln_str(*max_diff_r1) << ' '
        << scale * max_diff_r1->value << " vs " << max_diff_r2->value << '\n';
    out << differ_count << " of " << corr.n << " intensities differ too much (by >"
        << to_str(max_diff * 100) << "%).\n";
    if (differ_count >= 0.001 * corr.n)
      ok = false;
    else
      out << "(less than 0.1% of all intensities -"
          << " probably outlier rejection during merging)\n";
  }
  if (missing_count != 0) {
    out << missing_count << " out of " << mi.data.size()
        << " reflections in the merged file not found in unmerged data\n";
    ok = false;
  }
  if (!relaxed_check) {
    if (differ_count == 0 && missing_count == 0) {
      out << "Intensities match.";
      if (!ok)
        out << " But other problems were found (see above).";
      out << '\n';
    } else {
      out << (ok ? "OK. Intensities almost match.\n"
                 : "ERROR. Intensities do not match.\n");
    }
  }
  return ok;
}

#define WRITE(...) os.write(buf, gf_snprintf(buf, 255, __VA_ARGS__))

// Reorder eigenvalues and change signs of eigenvectors in the STARANISO way.
// It minimises the rotation angle from the basis vectors to the eigenvectors,
// which is equivalent to maximising the trace of the eigenvector matrix.
inline void reorder_staraniso_eigensystem(Mat33& vectors, double (&values)[3]) {
  const int8_t permut[6][3] = {{0,1,2}, {1,2,0}, {2,0,1}, {1,0,2}, {2,1,0}, {0,2,1}};
  const int8_t signs[8][3] = {{1,1,1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1},
                              {-1,-1,-1}, {-1,1,1}, {1,-1,1}, {1,1,-1}};
  double max_trace = -INFINITY;
  int permut_pos = 0, sign_pos = 0;
  bool det_neg = std::signbit(vectors.determinant());
  for (int i = 0; i < 6; ++i) {
    int jbase = det_neg == (i > 2) ? 0 : 4;
    const int8_t (&p)[3] = permut[i];
    for (int j = jbase; j < jbase+4; ++j) {
      double trace = 0.;
      for (int k = 0; k < 3; ++k)
        trace += signs[j][k] * vectors[k][p[k]];
      if (trace > max_trace) {
        max_trace = trace;
        permut_pos = i;
        sign_pos = j;
      }
    }
  }
  const int8_t (&p)[3] = permut[permut_pos];
  double tmp[3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      tmp[j] = signs[sign_pos][j] * vectors.a[i][p[j]];
    std::memcpy(vectors.a[i], tmp, sizeof(tmp));
  }
  for (int j = 0; j < 3; ++j)
    tmp[j] = values[p[j]];
  std::memcpy(values, tmp, sizeof(tmp));
}

inline void write_staraniso_b_in_mmcif(const SMat33<double>& b,
                                       char* buf, std::ostream& os) {
  double eigenvalues[3];
  Mat33 eigenvectors = eigen_decomposition(b, eigenvalues);
  reorder_staraniso_eigensystem(eigenvectors, eigenvalues);
  const char* prefix = "\n_reflns.pdbx_aniso_B_tensor_eigen";
  for (int i = 0; i < 3; ++i) {
    double v = std::fabs(eigenvalues[i]) > 1e-4 ? eigenvalues[i] : 0;
    WRITE("%svalue_%d %.5g", prefix, i+1, v);
    for (int j = 0; j < 3; ++j)
      WRITE("%svector_%d_ortho[%d] %.5g", prefix, i+1, j+1, eigenvectors[j][i]);
  }
  os << '\n';
}

inline void MtzToCif::write_cif(const Mtz& mtz, const Mtz* mtz2,
                                SMat33<double>* staraniso_b, std::ostream& os) {
  if (mtz2 && mtz.is_merged() == mtz2->is_merged())
    fail("If two MTZ files are given, one must be merged and one unmerged,\n"
         "got two ", mtz.is_merged() ? "merged" : "unmerged");
  const Mtz* merged = mtz.is_merged() ? &mtz : mtz2;
  const Mtz* unmerged = mtz.is_merged() ? mtz2 : &mtz;

  char buf[256];
  if (with_comments) {
    os << "# Converted by gemmi-mtz2cif " GEMMI_VERSION "\n";
    for (const Mtz* m : {merged, unmerged})
      if (m) {
        os << "# from " << (m->is_merged() ? "" : "un") << "merged MTZ: "
           << m->source_path << '\n';
        os << "#   title: " << m->title << '\n';
        if (with_history)
          for (size_t i = 0; i != m->history.size(); ++i)
            os << "#   history #" << i << ": " << m->history[i] << '\n';
      }
  }
  os << "data_" << (block_name ? block_name : "mtz");

  os << "\n\n_entry.id " << entry_id << "\n\n";

  write_special_marker_if_requested(os, merged);

  if (unmerged) {
    bool ok = gather_sweep_data(*unmerged);
    std::set<const Mtz::Dataset*> used_datasets;
    std::vector<std::string> crystal_names;
    // Prepare used_datasets, crystal_names and set SweepData::crystal_id.
    for (SweepData& sweep : sweeps)
      if (sweep.dataset) {
        used_datasets.insert(sweep.dataset);
        if (ok) {
          auto it = std::find(crystal_names.begin(), crystal_names.end(),
                              sweep.dataset->crystal_name);
          sweep.crystal_id = int(it - crystal_names.begin()) + 1;
          if (it == crystal_names.end())
            crystal_names.push_back(sweep.dataset->crystal_name);
        }
      }

    if (crystal_names.empty()) {
      os << "_exptl_crystal.id 1\n";
    } else {
      os << "loop_\n_exptl_crystal.id\n";
      for (size_t i = 0; i < crystal_names.size(); ++i)
        os << i+1 << " # " << crystal_names[i] << '\n';
    }

    bool scaled = (unmerged->column_with_label("SCALEUSED") != nullptr);

    os << "\nloop_\n"
          "_diffrn.id\n_diffrn.crystal_id\n_diffrn.details\n";
    for (const SweepData& sweep : sweeps) {
      os << sweep.id << ' ' << sweep.crystal_id << " '";
      if (scaled)
        os << "scaled ";
      os << "unmerged data'\n";
    }

    os << "\nloop_\n"
          "_diffrn_measurement.diffrn_id\n"
          "_diffrn_measurement.details\n";
    for (const SweepData& sweep : sweeps)
      os << sweep.id << " '" << sweep.batch_count << " frames'\n";

    os << "\nloop_\n"
          "_diffrn_radiation.diffrn_id\n"
          "_diffrn_radiation.wavelength_id\n";
    for (const SweepData& sweep : sweeps) {
      os << sweep.id << ' ';
      if (sweep.dataset)
        os << sweep.dataset->id;
      else
        os << '?';
      os << '\n';
    }

    os << "\nloop_\n"
          "_diffrn_radiation_wavelength.id\n"
          "_diffrn_radiation_wavelength.wavelength\n";
    for (const Mtz::Dataset* ds : used_datasets)
      os << ds->id << ' ' << ds->wavelength << '\n';
    os << '\n';

    if (enable_UB) {
      os << "loop_\n"
            "_diffrn_orient_matrix.diffrn_id\n"
            "_diffrn_orient_matrix.UB[1][1]\n"
            "_diffrn_orient_matrix.UB[1][2]\n"
            "_diffrn_orient_matrix.UB[1][3]\n"
            "_diffrn_orient_matrix.UB[2][1]\n"
            "_diffrn_orient_matrix.UB[2][2]\n"
            "_diffrn_orient_matrix.UB[2][3]\n"
            "_diffrn_orient_matrix.UB[3][1]\n"
            "_diffrn_orient_matrix.UB[3][2]\n"
            "_diffrn_orient_matrix.UB[3][3]\n";
      for (const SweepData& sweep : sweeps) {
        Mat33 u = sweep.first_batch->matrix_U();
        Mat33 b = sweep.first_batch->get_cell().calculate_matrix_B();
        Mat33 ub = u.multiply(b);
        WRITE("%d  %#g %#g %#g  %#g %#g %#g  %#g %#g %#g\n", sweep.id,
              ub.a[0][0], ub.a[0][1], ub.a[0][2],
              ub.a[1][0], ub.a[1][1], ub.a[1][2],
              ub.a[2][0], ub.a[2][1], ub.a[2][2]);
      }
      os << '\n';
    }
  } else {
    double w = std::isnan(wavelength) ? get_wavelength(mtz, recipe) : wavelength;
    if (w > 0.)
      os << "_diffrn_radiation.diffrn_id 1\n"
         << "_diffrn_radiation.wavelength_id 1\n"
         << "_diffrn_radiation_wavelength.id 1\n"
         << "_diffrn_radiation_wavelength.wavelength " << to_str(w) << "\n\n";
  }

  if (merged) {
    write_cell_and_symmetry(mtz.get_cell(), nullptr, mtz.spacegroup, buf, os);
  } else {
    double rmsds[6];
    UnitCell cell = unmerged->get_average_cell_from_batch_headers(rmsds);
    write_cell_and_symmetry(cell, rmsds, mtz.spacegroup, buf, os);
  }

  if (write_staraniso_tensor && staraniso_b)
    write_staraniso_b_in_mmcif(*staraniso_b, buf, os);

  if (merged)
    write_main_loop(*merged, buf, os);
  if (unmerged)
    write_main_loop(*unmerged, buf, os);
}

inline void MtzToCif::write_main_loop(const Mtz& mtz, char* buf, std::ostream& os) {
  prepare_recipe(mtz);
  // prepare indices
  std::vector<int> value_indices;  // used for --skip_empty
  std::vector<int> sigma_indices;  // used for status 'x'
  for (const Trans& tr : recipe) {
    if (tr.col_idx < 0)
      continue;
    const Mtz::Column& col = mtz.columns[tr.col_idx];
    if (skip_empty) {
      if (skip_empty_cols.empty() ? col.type != 'H' && col.type != 'I'
                                  : is_in_list(col.label, skip_empty_cols))
        value_indices.push_back(tr.col_idx);
    }
    if (col.type != 'Q' && col.type != 'L' && col.type != 'M')
      sigma_indices.push_back(tr.col_idx);
  }

  bool unmerged = !mtz.is_merged();
  os << "\nloop_\n";
  for (const Trans& tr : recipe) {
    os << (unmerged ? "_diffrn_refln." : "_refln.");
    if (with_comments && tr.col_idx >= 0) {
      WRITE("%-26s # ", tr.tag.c_str());
      const Mtz::Column& col = mtz.columns.at(tr.col_idx);
      const Mtz::Dataset& ds = mtz.dataset(col.dataset_id);
      // dataset is assigned to column only in merged MTZ
      if (unmerged)
        os << col.label;
      else
        WRITE("%-14s from dataset %s", col.label.c_str(), ds.dataset_name.c_str());
    } else {
      os << tr.tag;
    }
    os << '\n';
  }
  int batch_idx = find_column_index("BATCH", mtz);
  std::unordered_map<int, const Mtz::Batch*> batch_by_number;
  for (const Mtz::Batch& b : mtz.batches)
    batch_by_number.emplace(b.number, &b);

  if (free_flag_value < 0) {
    // CCP4 uses flags 0,...N-1 (usually N=20), with default free set 0
    // PHENIX uses 0/1 flags with free set 1
    if (const Trans* tr_status = get_status_translation()) {
      int count = 0;
      for (float val : mtz.columns[tr_status->col_idx])
        if (val == 0.f)
          count++;
      free_flag_value = count < mtz.nreflections / 2 ? 0 : 1;
    }
  }

  auto write_int = [](char* p, int num) {
    //return gf_snprintf(p, 32, "%d", num);
    std::string s = std::to_string(num);
    std::memcpy(p, s.data(), s.size());
    return s.size();
  };

  char* ptr = buf;
  for (int i = 0, idx = 0; i != mtz.nreflections; ++i) {
    const float* row = &mtz.data[i * mtz.columns.size()];
    if (trim > 0) {
      if (row[0] < -trim || row[0] > trim ||
          row[1] < -trim || row[1] > trim ||
          row[2] < -trim || row[2] > trim) {
        continue;
      }
    }
    if (!value_indices.empty())
      if (std::all_of(value_indices.begin(), value_indices.end(),
                      [&](int n) { return std::isnan(row[n]); }))
        continue;
    int batch_number = 0;
    SweepData* sweep = nullptr;
    if (unmerged) {
      if (batch_idx == -1)
        fail("BATCH column not found");
      batch_number = (int) row[batch_idx];
      auto it = batch_by_number.find(batch_number);
      if (it == batch_by_number.end())
        fail("unexpected values in column BATCH");
      const Mtz::Batch& batch = *it->second;
      sweep = &sweeps[sweep_indices.at(batch.number)];
    }
    bool first = true;
    for (const Trans& tr : recipe) {
      if (first)
        first = false;
      else
        *ptr++ = ' ';
      if (ptr - buf > 220) {
        os.write(buf, ptr - buf);
        ptr = buf;
      }
      if (tr.col_idx < 0) {
        switch (tr.col_idx) {
          case Dot: *ptr++ = '.'; break;
          case Qmark: *ptr++ = '?'; break;
          case Counter: ptr += write_int(ptr, ++idx); break;
          case DatasetId: ptr += write_int(ptr, sweep->id); break;
          case Image: ptr += write_int(ptr, batch_number - sweep->offset); break;
        }
      } else {
        float v = row[tr.col_idx];
        if (tr.is_status) {
          char status = 'x';
          if (sigma_indices.empty() ||
              !std::all_of(sigma_indices.begin(), sigma_indices.end(),
                           [&](int n) { return std::isnan(row[n]); }))
            status = int(v) == free_flag_value ? 'f' : 'o';
          *ptr++ = status;
        } else if (std::isnan(v)) {
          // we checked that min_width <= 32
          for (int j = 1; j < tr.min_width; ++j)
            *ptr++ = ' ';
          *ptr++ = '?';
        } else {
#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wformat-nonliteral"
#endif
          ptr += gf_snprintf(ptr, 32, tr.format.c_str(), v);
#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif
        }
      }
    }
    *ptr++ = '\n';
  }
  os.write(buf, ptr - buf);
}

inline void MtzToCif::write_cif_from_xds(const XdsAscii& xds, std::ostream& os) {
  char buf[256];
  if (with_comments) {
    os << "# Converted by gemmi-mtz2cif " GEMMI_VERSION "\n";
    os << "# from scaled unmerged XDS_ASCII: " << xds.source_path << '\n';
  }
  os << "data_" << (block_name ? block_name : "xds");
  os << "\n\n_entry.id " << entry_id << "\n\n";

  write_special_marker_if_requested(os, false);

  os << "_exptl_crystal.id 1\n";

  os << "\nloop_\n_diffrn.id\n_diffrn.crystal_id\n_diffrn.details\n";
  for (const XdsAscii::Iset& iset : xds.isets)
    // We could write iset.input_info as details, but then local paths
    // could end up in the deposition.
    os << iset.id << " 1 ?\n";
  os << '\n';

  os << "loop_\n"
        "_diffrn_measurement.diffrn_id\n"
        "_diffrn_measurement.details\n";
  for (const XdsAscii::Iset& iset : xds.isets) {
    os << iset.id;
    if (iset.frame_count >= 0)
      os << " '" << iset.frame_count << " frames'\n";
    else
      os << " ?\n";
  }
  os << '\n';

  double w_all = wavelength;
  if (std::isnan(w_all) &&
      std::all_of(xds.isets.begin(), xds.isets.end(),
                  [&](const XdsAscii::Iset& i) { return i.wavelength == xds.wavelength; }))
    w_all = xds.wavelength;

  os << "loop_\n"
        "_diffrn_radiation.diffrn_id\n"
        "_diffrn_radiation.wavelength_id\n";
  for (const XdsAscii::Iset& iset : xds.isets)
    os << iset.id << ' ' << (std::isnan(w_all) ? iset.id : 1) << '\n';
  os << '\n';

  os << "loop_\n"
        "_diffrn_radiation_wavelength.id\n"
        "_diffrn_radiation_wavelength.wavelength\n";
  auto number_or_dot = [](double d) { return std::isnan(d) ? "." : to_str(d); };
  if (std::isnan(w_all)) {
    for (const XdsAscii::Iset& iset : xds.isets)
      os << iset.id << ' ' << number_or_dot(iset.wavelength) << '\n';
  } else {
    os << '1' << ' ' << number_or_dot(w_all) << '\n';
  }
  os << '\n';

  const SpaceGroup* sg = find_spacegroup_by_number(xds.spacegroup_number);
  double rmsds[6] = {0., 0., 0., 0., 0., 0.};
  if (xds.isets.size() > 1) {
    double mean[6] = {xds.unit_cell.a, xds.unit_cell.b, xds.unit_cell.c,
                      xds.unit_cell.alpha, xds.unit_cell.beta, xds.unit_cell.gamma};
    int n = 0;
    for (const XdsAscii::Iset& iset : xds.isets) {
      for (int i = 0; i < 6; ++i)
        rmsds[i] += sq(mean[i] - iset.cell_constants[i]) * iset.frame_count;
      n += iset.frame_count;
    }
    for (int i = 0; i < 6; ++i)
      rmsds[i] = std::sqrt(rmsds[i] / n);
  }
  write_cell_and_symmetry(xds.unit_cell, rmsds, sg, buf, os);

  os << "\nloop_"
        "\n_diffrn_refln.diffrn_id"
        "\n_diffrn_refln.id"
        "\n_diffrn_refln.index_h"
        "\n_diffrn_refln.index_k"
        "\n_diffrn_refln.index_l"
        "\n_diffrn_refln.intensity_net"
        "\n_diffrn_refln.intensity_sigma";
  if (xds.oscillation_range != 0.)
    os << "\n_diffrn_refln.pdbx_scan_angle";
  os << "\n_diffrn_refln.pdbx_image_id\n";
  int idx = 0;
  for (const XdsAscii::Refl& refl : xds.data) {
    if (refl.sigma < 0)  // misfit
      continue;
    char* ptr = buf;
    ptr += gf_snprintf(ptr, 128, "%d %d %d %d %d %g %.5g ",
                       refl.iset, ++idx, refl.hkl[0], refl.hkl[1], refl.hkl[2],
                       refl.iobs, refl.sigma);
    if (xds.oscillation_range != 0.) {
      double z = refl.zd - xds.starting_frame + 1;
      double angle = xds.starting_angle + xds.oscillation_range * z;
      ptr += gf_snprintf(ptr, 16, "%.5g ", angle);
    }
    ptr += gf_snprintf(ptr, 16, "%d\n", int(std::ceil(refl.zd)));
    os.write(buf, ptr - buf);
  }
}

inline void MtzToCif::write_cell_and_symmetry(const UnitCell& cell, double* rmsds,
                                              const SpaceGroup* sg,
                                              char* buf, std::ostream& os) const {
  os << "_cell.entry_id " << entry_id << '\n';
  WRITE("_cell.length_a    %8.4f\n", cell.a);
  if (rmsds && rmsds[0] != 0.)
    WRITE("_cell.length_a_esd %7.3f\n", rmsds[0]);
  WRITE("_cell.length_b    %8.4f\n", cell.b);
  if (rmsds && rmsds[1] != 0.)
    WRITE("_cell.length_b_esd %7.3f\n", rmsds[1]);
  WRITE("_cell.length_c    %8.4f\n", cell.c);
  if (rmsds && rmsds[2] != 0.)
    WRITE("_cell.length_c_esd %7.3f\n", rmsds[2]);
  WRITE("_cell.angle_alpha %8.4f\n", cell.alpha);
  if (rmsds && rmsds[3] != 0.)
    WRITE("_cell.angle_alpha_esd %7.3f\n", rmsds[3]);
  WRITE("_cell.angle_beta  %8.4f\n", cell.beta);
  if (rmsds && rmsds[4] != 0.)
    WRITE("_cell.angle_beta_esd %8.3f\n", rmsds[4]);
  WRITE("_cell.angle_gamma %8.4f\n", cell.gamma);
  if (rmsds && rmsds[5] != 0.)
    WRITE("_cell.angle_gamma_esd %7.3f\n", rmsds[5]);
  if (sg) {
    os << "\n_symmetry.entry_id " << entry_id << "\n"
          "_symmetry.space_group_name_H-M '" << sg->hm << "'\n"
          "_symmetry.Int_Tables_number " << sg->number << '\n';
    // could write _symmetry_equiv.pos_as_xyz, but would it be useful?
  }
}

#undef WRITE

} // namespace gemmi
#endif
