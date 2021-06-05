// Copyright 2019 Global Phasing Ltd.
//
// Function read_metadata_from_remarks() that interprets REMARK 3
// and REMARK 200/230/240 filling in Metadata.

#ifndef GEMMI_REMARKS_HPP_
#define GEMMI_REMARKS_HPP_

#include "metadata.hpp"
#include "pdb.hpp"  // for pdb_date_format_to_iso, read_int, ...

namespace gemmi {

namespace pdb_impl {

inline bool is_double(const char* p) {
  while (std::isspace(*p)) ++p;
  if (*p == '-' || *p == '+') ++p;
  while (is_digit(*p)) ++p;
  if (*p == '.') {
    ++p;
    while (is_digit(*++p)) ++p;
  }
  while (std::isspace(*p)) ++p;
  return *p == '\0';
}

inline bool is_tls_item(const std::string& key) {
  return key.size() == 3 &&
    (key[0] == 'T' || key[0] == 'L' || key[0] == 'S') &&
    (key[1] == '1' || key[1] == '2' || key[1] == '3') &&
    (key[2] == '1' || key[2] == '2' || key[2] == '3');
}

// Usually we have one program per line:
//   XDS
//   XDS VERSION NOVEMBER 3, 2014
//   AIMLESS 0.5.17
// but it can also be a list of programs:
//   autoPROC (Version 1.3.0), AIMLESS, STARANISO
//   autoPROC, XDS (VERSION Jan 26, 2018)
// We assume that:
// - the name has only one word (apologies to Queen of Spades,
//   Force Field X, APEX 2 and Insight II).
// - comma not followed by a digit separates programs
// - brackets and the word VERSION are to be removed from version
// Additionally, if version has format: "something (DATE)" where
// the DATE format is either 28-MAR-07 or 28-Mar-2007, then DATE
// is put into _software.date.
inline void add_software(Metadata& meta, SoftwareItem::Classification type,
                         const std::string& name) {
  for (size_t start = 0, end = 0; end != std::string::npos; start = end + 1) {
    end = name.find(',', start);
    while (end != std::string::npos &&
           name[end+1] == ' ' && is_digit(name[end+2]))
      end = name.find(',', end + 1);
    meta.software.emplace_back();
    SoftwareItem& item = meta.software.back();
    item.name = trim_str(name.substr(start, end - start));
    size_t sep = item.name.find(' ');
    if (sep != std::string::npos) {
      size_t ver_start = item.name.find_first_not_of(" (", sep + 1);
      item.version = item.name.substr(ver_start);
      item.name.resize(sep);
      if (!item.version.empty() && item.version.back() == ')') {
        size_t open_br = item.version.find('(');
        if (open_br == std::string::npos) {
          item.version.pop_back();
        } else if (open_br + 11 == item.version.size() ||
                   open_br + 13 == item.version.size()) {
          item.date = pdb_date_format_to_iso(item.version.substr(open_br + 1));
          if (item.date.size() == 10 && item.date[5] != 'x') {
            size_t last = item.version.find_last_not_of(' ', open_br - 1);
            item.version.resize(last + 1);
          } else {
            item.date.clear();
          }
        }
      }
      if (istarts_with(item.version, "version "))
        item.version.erase(0, 8);
    }
    item.classification = type;
    item.pdbx_ordinal = (int) meta.software.size();
  }
}

// REMARK   3   TERM                          COUNT    WEIGHT   FUNCTION.
// REMARK   3    BOND LENGTHS              : 5760   ; 2.000  ; HARMONIC
inline void add_restraint_count_weight(RefinementInfo& ref_info,
                                       const char* key, const char* value) {
  if (*value == 'N') // NULL instead of number
    return;
  ref_info.restr_stats.emplace_back(key);
  RefinementInfo::Restr& restr = ref_info.restr_stats.back();
  const char* endptr;
  restr.count = no_sign_atoi(value, &endptr);
  if (const char* sep = std::strchr(endptr, ';'))
    restr.weight = fast_atof(sep + 1, &endptr);
  if (const char* sep = std::strchr(endptr, ';'))
    restr.function = read_string(sep+1, 50);
}

inline void read_remark3_line(const char* line, Metadata& meta,
                              std::string*& possibly_unfinished_remark3) {
  // Based on:
  // www.wwpdb.org/documentation/file-format-content/format23/remark3.html
  // and analysis of PDB files.
  // In special cases, such as joint X-ray and neutron refinement 5MOO,
  // PDB file can have two REMARK 3 blocks.
  // Generally, after "REMARK   3" we have either a header-like sentance
  // or a key:value pair with a colon, or a continuation of text from the
  // previous line.
  const char* key_start = skip_blank(line + 10);
  const char* colon = std::strchr(key_start, ':');
  const char* key_end = rtrim_cstr(key_start, colon);
  std::string key(key_start, key_end);

  // multi-line continuation requires special handling
  if (possibly_unfinished_remark3) {
    if (key_start > line + 17) {
      *possibly_unfinished_remark3 += ' ';
      possibly_unfinished_remark3->append(key);
      return;
    }
    possibly_unfinished_remark3 = nullptr;
  }

  if (colon) {
    const char* value = skip_blank(colon + 1);
    const char* end = rtrim_cstr(value);
    if (end - value == 4 && std::strncmp(value, "NULL", 4) == 0)
      return;
    if (same_str(key, "PROGRAM"))
      add_software(meta, SoftwareItem::Refinement, std::string(value, end));
    if (meta.refinement.empty())
      return;
    RefinementInfo& ref_info = meta.refinement.back();
    if (same_str(key, "RESOLUTION RANGE HIGH (ANGSTROMS)")) {
      ref_info.resolution_high = fast_atof(value);
    } else if (same_str(key, "RESOLUTION RANGE LOW  (ANGSTROMS)")) {
      ref_info.resolution_low = fast_atof(value);
    } else if (same_str(key, "COMPLETENESS FOR RANGE        (%)")) {
      ref_info.completeness = fast_atof(value);
    } else if (same_str(key, "NUMBER OF REFLECTIONS")) {
      ref_info.reflection_count = std::atoi(value);
    } else if (same_str(key, "CROSS-VALIDATION METHOD")) {
      ref_info.cross_validation_method = std::string(value, end);
    } else if (same_str(key, "FREE R VALUE TEST SET SELECTION")) {
      ref_info.rfree_selection_method = std::string(value, end);
    } else if (same_str(key, "R VALUE     (WORKING + TEST SET)")) {
      ref_info.r_all = fast_atof(value);
    } else if (same_str(key, "R VALUE            (WORKING SET)")) {
      ref_info.r_work = fast_atof(value);
    } else if (same_str(key, "FREE R VALUE")) {
      ref_info.r_free = fast_atof(value);
    } else if (same_str(key, "FREE R VALUE TEST SET COUNT")) {
      ref_info.rfree_set_count = atoi(value);
    } else if (same_str(key, "TOTAL NUMBER OF BINS USED")) {
      ref_info.bin_count = std::atoi(value);
    } else if (same_str(key, "BIN RESOLUTION RANGE HIGH       (A)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().resolution_high = fast_atof(value);
    } else if (same_str(key, "BIN RESOLUTION RANGE LOW        (A)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().resolution_low = fast_atof(value);
    } else if (same_str(key, "BIN COMPLETENESS (WORKING+TEST) (%)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().completeness = fast_atof(value);
    } else if (same_str(key, "REFLECTIONS IN BIN   (WORKING+TEST)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().reflection_count = std::atoi(value);
    } else if (same_str(key, "BIN R VALUE          (WORKING+TEST)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().r_all = fast_atof(value);
    } else if (same_str(key, "BIN R VALUE           (WORKING SET)")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().r_work = fast_atof(value);
    } else if (same_str(key, "BIN FREE R VALUE")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().r_free = fast_atof(value);
    } else if (same_str(key, "BIN FREE R VALUE TEST SET COUNT")) {
      if (!ref_info.bins.empty())
        ref_info.bins.back().rfree_set_count = std::atoi(value);
    } else if (same_str(key, "FROM WILSON PLOT           (A**2)")) {
      // TODO
      // exper.b_wilson = fast_atof(value);
    } else if (same_str(key, "MEAN B VALUE      (OVERALL, A**2)")) {
      ref_info.mean_b = fast_atof(value);
    } else if (same_str(key, "B11 (A**2)")) {
      ref_info.aniso_b[0][0] = fast_atof(value);
    } else if (same_str(key, "B22 (A**2)")) {
      ref_info.aniso_b[1][1] = fast_atof(value);
    } else if (same_str(key, "B33 (A**2)")) {
      ref_info.aniso_b[2][2] = fast_atof(value);
    } else if (same_str(key, "B12 (A**2)")) {
      ref_info.aniso_b[0][1] = fast_atof(value);
    } else if (same_str(key, "B13 (A**2)")) {
      ref_info.aniso_b[0][2] = fast_atof(value);
    } else if (same_str(key, "B23 (A**2)")) {
      ref_info.aniso_b[1][2] = fast_atof(value);
    } else if (same_str(key, "ESD FROM LUZZATI PLOT                    (A)")) {
      ref_info.luzzati_error = fast_atof(value);
    } else if (same_str(key, "DPI (BLOW EQ-10) BASED ON R VALUE        (A)")) {
      ref_info.dpi_blow_r = fast_atof(value);
    } else if (same_str(key, "DPI (BLOW EQ-9) BASED ON FREE R VALUE    (A)")) {
      ref_info.dpi_blow_rfree = fast_atof(value);
    } else if (same_str(key, "DPI (CRUICKSHANK) BASED ON R VALUE       (A)")) {
      ref_info.dpi_cruickshank_r = fast_atof(value);
    } else if (same_str(key, "DPI (CRUICKSHANK) BASED ON FREE R VALUE  (A)")) {
      ref_info.dpi_cruickshank_rfree = fast_atof(value);
    } else if (same_str(key, "CORRELATION COEFFICIENT FO-FC")) {
      ref_info.cc_fo_fc = fast_atof(value);
    } else if (same_str(key, "CORRELATION COEFFICIENT FO-FC FREE")) {
      ref_info.cc_fo_fc_free = fast_atof(value);
    } else if (same_str(key, "BOND LENGTHS")) {
      add_restraint_count_weight(ref_info, "t_bond_d", value);
    } else if (same_str(key, "BOND ANGLES")) {
      add_restraint_count_weight(ref_info, "t_angle_deg", value);
    } else if (same_str(key, "TORSION ANGLES")) {
      add_restraint_count_weight(ref_info, "t_dihedral_angle_d", value);
    } else if (same_str(key, "TRIGONAL CARBON PLANES")) {
      add_restraint_count_weight(ref_info, "t_trig_c_planes", value);
    } else if (same_str(key, "GENERAL PLANES")) {
      add_restraint_count_weight(ref_info, "t_gen_planes", value);
    } else if (same_str(key, "ISOTROPIC THERMAL FACTORS")) {
      add_restraint_count_weight(ref_info, "t_it", value);
    } else if (same_str(key, "BAD NON-BONDED CONTACTS")) {
      add_restraint_count_weight(ref_info, "t_nbd", value);
    } else if (same_str(key, "IMPROPER TORSIONS")) {
      add_restraint_count_weight(ref_info, "t_improper_torsion", value);
    } else if (same_str(key, "CHIRAL IMPROPER TORSION")) {
      add_restraint_count_weight(ref_info, "t_chiral_improper_torsion", value);
    } else if (same_str(key, "SUM OF OCCUPANCIES")) {
      add_restraint_count_weight(ref_info, "t_sum_occupancies", value);
    } else if (same_str(key, "UTILITY DISTANCES")) {
      add_restraint_count_weight(ref_info, "t_utility_distance", value);
    } else if (same_str(key, "UTILITY ANGLES")) {
      add_restraint_count_weight(ref_info, "t_utility_angle", value);
    } else if (same_str(key, "UTILITY TORSION")) {
      add_restraint_count_weight(ref_info, "t_utility_torsion", value);
    } else if (same_str(key, "IDEAL-DIST CONTACT TERM")) {
      add_restraint_count_weight(ref_info, "t_ideal_dist_contact", value);
    } else if (same_str(key, "BOND LENGTHS                       (A)")) {
      impl::find_or_add(ref_info.restr_stats, "t_bond_d").dev_ideal
        = read_double(value, 50);
    } else if (same_str(key, "BOND ANGLES                  (DEGREES)")) {
      impl::find_or_add(ref_info.restr_stats, "t_angle_deg").dev_ideal
        = read_double(value, 50);
    } else if (same_str(key, "PEPTIDE OMEGA TORSION ANGLES (DEGREES)")) {
      impl::find_or_add(ref_info.restr_stats, "t_omega_torsion").dev_ideal
        = read_double(value, 50);
    } else if (same_str(key, "OTHER TORSION ANGLES         (DEGREES)")) {
      impl::find_or_add(ref_info.restr_stats, "t_other_torsion").dev_ideal
        = read_double(value, 50);
    } else if (same_str(key, "TLS GROUP")) {
      ref_info.tls_groups.emplace_back();
      ref_info.tls_groups.back().id = std::string(value, end);
    } else if (same_str(key, "SET") ||
               // "REMARK   3    SELECTION:"            -> TLS
               // "REMARK   3     SELECTION          :" -> NCS
               (same_str(key, "SELECTION") && colon == line + 23)) {
      if (!ref_info.tls_groups.empty()) {
        TlsGroup& group = ref_info.tls_groups.back();
        group.selections.emplace_back();
        group.selections.back().details = std::string(value, end);
        possibly_unfinished_remark3 = &group.selections.back().details;
      }
    } else if (same_str(key, "RESIDUE RANGE")) {
      if (!ref_info.tls_groups.empty() && end > colon+21) {
        TlsGroup& group = ref_info.tls_groups.back();
        group.selections.emplace_back();
        TlsGroup::Selection& sel = group.selections.back();
        sel.chain = read_string(colon+1, 5);
        if (sel.chain == read_string(colon+16, 5)) {
          try {
            sel.res_begin = SeqId(read_string(colon+6, 6));
            sel.res_end = SeqId(read_string(colon+21, 6));
          } catch (std::invalid_argument&) {
            group.selections.pop_back();
          }
        } else {  // unexpected -- TLS group should be in one chain
          group.selections.pop_back();
        }
      }
    } else if (same_str(key, "ORIGIN FOR THE GROUP (A)")) {
      std::vector<std::string> xyz = split_str_multi(std::string(value, end));
      if (ref_info.tls_groups.empty() || xyz.size() != 3)
        return;
      Position& origin = ref_info.tls_groups.back().origin;
      origin.x = fast_atof(xyz[0].c_str());
      origin.y = fast_atof(xyz[1].c_str());
      origin.z = fast_atof(xyz[2].c_str());
    } else if (is_tls_item(key)) {
      if (ref_info.tls_groups.empty())
        return;
      TlsGroup& tls = ref_info.tls_groups.back();
      std::vector<std::string> tokens = split_str_multi(key_start);
      for (size_t i = 0; i + 1 < tokens.size(); i += 2) {
        std::string& k = tokens[i];
        if (k.size() == 4 && k[3] == ':')
          k.resize(3);
        if (is_tls_item(k)) {
          Mat33& m = k[0] == 'T' ? tls.T : k[0] == 'L' ? tls.L : tls.S;
          int x = k[1] - '1';
          int y = k[2] - '1';
          m[x][y] = m[y][x] = fast_atof(tokens[i+1].c_str());
        }
      }
    }
  } else {
    if (same_str(key, "DATA USED IN REFINEMENT.")) {
      meta.refinement.emplace_back();
      meta.refinement.back().id = std::to_string(meta.refinement.size());
    } else if (same_str(key, "FIT IN THE HIGHEST RESOLUTION BIN.")) {
      if (!meta.refinement.empty())
        meta.refinement.back().bins.emplace_back();
    }
  }
}

inline void read_remark_200_230_240(const char* line, Metadata& meta,
                                    std::string*& cryst_desc) {
  // multi-line continuation requires special handling
  if (cryst_desc) {
    if (line[10] == ' ' && line[11] == ' ') {
      const char* start = line + 11;
      cryst_desc->append(start, rtrim_cstr(start) - start);
      return;
    }
    cryst_desc = nullptr;
  }

  const char* key_start = skip_blank(line + 10);
  const char* colon = std::strchr(key_start, ':');
  const char* key_end = rtrim_cstr(key_start, colon);
  std::string key(key_start, key_end);
  if (colon) {
    const char* value = skip_blank(colon + 1);
    const char* end = rtrim_cstr(value);
    if (end - value == 4 && std::strncmp(value, "NULL", 4) == 0)
      return;
    if (same_str(key, "INTENSITY-INTEGRATION SOFTWARE")) {
      add_software(meta, SoftwareItem::DataReduction, std::string(value, end));
    } else if (same_str(key, "DATA SCALING SOFTWARE")) {
      add_software(meta, SoftwareItem::DataScaling, std::string(value, end));
    } else if (same_str(key, "SOFTWARE USED")) {
      add_software(meta, SoftwareItem::Phasing, std::string(value, end));
    } else if (same_str(key, "METHOD USED TO DETERMINE THE STRUCTURE")) {
      meta.solved_by = std::string(value, end);
    } else if (same_str(key, "STARTING MODEL")) {
      meta.starting_model = std::string(value, end);
    } else if (!meta.experiments.empty()) {
      ExperimentInfo& exper = meta.experiments.back();
      DiffractionInfo& diffr = meta.crystals.back().diffractions[0];
      if (same_str(key, "EXPERIMENT TYPE")) {
        exper.method = std::string(value, end);
      } else if (same_str(key, "NUMBER OF CRYSTALS USED")) {
        exper.number_of_crystals = std::atoi(value);
      } else if (same_str(key, "PH")) {
        if (is_double(value))
          meta.crystals.back().ph = fast_atof(value);
        else
          meta.crystals.back().ph_range = std::string(value, end);
      } else if (same_str(key, "DATE OF DATA COLLECTION")) {
        diffr.collection_date = pdb_date_format_to_iso(std::string(value, end));
      } else if (same_str(key, "TEMPERATURE           (KELVIN)")) {
        diffr.temperature = fast_atof(value);
      } else if (same_str(key, "SYNCHROTRON              (Y/N)")) {
        if (*value == 'Y')
          diffr.source = "SYNCHROTRON";
      } else if (same_str(key, "RADIATION SOURCE")) {
        if (same_str(diffr.source, "SYNCHROTRON"))
          diffr.synchrotron = std::string(value, end);
        else
          diffr.source = std::string(value, end);
      } else if (same_str(key, "NEUTRON SOURCE")) {
        diffr.source = std::string(value, end);
      } else if (same_str(key, "BEAMLINE")) {
        diffr.beamline = std::string(value, end);
        if (!diffr.synchrotron.empty() && diffr.source_type.empty())
          diffr.source_type = diffr.synchrotron + " BEAMLINE " + diffr.beamline;
      } else if (same_str(key, "X-RAY GENERATOR MODEL")) {
        diffr.source_type = std::string(value, end);
      } else if (same_str(key, "MONOCHROMATIC OR LAUE    (M/L)")) {
        diffr.mono_or_laue = *value;
      } else if (same_str(key, "WAVELENGTH OR RANGE        (A)")) {
        diffr.wavelengths = std::string(value, end);
      } else if (same_str(key, "MONOCHROMATOR")) {
        diffr.monochromator = std::string(value, end);
      } else if (same_str(key, "OPTICS")) {
        diffr.optics = std::string(value, end);
      } else if (same_str(key, "DETECTOR TYPE")) {
        diffr.detector = std::string(value, end);
      } else if (same_str(key, "DETECTOR MANUFACTURER")) {
        diffr.detector_make = std::string(value, end);
      } else if (same_str(key, "NUMBER OF UNIQUE REFLECTIONS")) {
        exper.unique_reflections = std::atoi(value);
      } else if (same_str(key, "RESOLUTION RANGE HIGH      (A)")) {
        exper.reflections.resolution_high = fast_atof(value);
      } else if (same_str(key, "RESOLUTION RANGE LOW       (A)")) {
        exper.reflections.resolution_low = fast_atof(value);
      } else if (same_str(key, "COMPLETENESS FOR RANGE     (%)")) {
        exper.reflections.completeness = fast_atof(value);
      } else if (same_str(key, "DATA REDUNDANCY")) {
        exper.reflections.redundancy = fast_atof(value);
      } else if (same_str(key, "R MERGE                    (I)")) {
        exper.reflections.r_merge = fast_atof(value);
      } else if (same_str(key, "R SYM                      (I)")) {
        exper.reflections.r_sym = fast_atof(value);
      } else if (same_str(key, "<I/SIGMA(I)> FOR THE DATA SET")) {
        exper.reflections.mean_I_over_sigma = fast_atof(value);
      } else if (same_str(key, "REMARK")) {
        cryst_desc = &meta.crystals.back().description;
        *cryst_desc = std::string(value, end);
      } else if (!exper.shells.empty()) {
        if (same_str(key, "HIGHEST RESOLUTION SHELL, RANGE HIGH (A)")) {
          exper.shells.back().resolution_high = fast_atof(value);
        } else if (same_str(key, "HIGHEST RESOLUTION SHELL, RANGE LOW  (A)")) {
          exper.shells.back().resolution_low = fast_atof(value);
        } else if (same_str(key, "COMPLETENESS FOR SHELL     (%)")) {
          exper.shells.back().completeness = fast_atof(value);
        } else if (same_str(key, "DATA REDUNDANCY IN SHELL")) {
          exper.shells.back().redundancy = fast_atof(value);
        } else if (same_str(key, "R MERGE FOR SHELL          (I)")) {
          exper.shells.back().r_merge = fast_atof(value);
        } else if (same_str(key, "R SYM FOR SHELL            (I)")) {
          exper.shells.back().r_sym = fast_atof(value);
        } else if (same_str(key, "<I/SIGMA(I)> FOR SHELL")) {
          exper.shells.back().mean_I_over_sigma = fast_atof(value);
        }
      }
    }
  } else {
    if (same_str(key, "EXPERIMENTAL DETAILS")) {
      meta.crystals.emplace_back();
      CrystalInfo& c = meta.crystals.back();
      c.id = std::to_string(meta.crystals.size());
      c.diffractions.emplace_back();
      c.diffractions[0].id = c.id;
      meta.experiments.emplace_back();
      meta.experiments.back().diffraction_ids.push_back(c.id);
      if (line[8] == '0' && line[9] == '0')
        c.diffractions[0].scattering_type = "x-ray";
      else if (line[8] == '3' && line[9] == '0')
        c.diffractions[0].scattering_type = "neutron";
      else if (line[8] == '4' && line[9] == '0')
        c.diffractions[0].scattering_type = "electron";
    }
    if (same_str(key, "IN THE HIGHEST RESOLUTION SHELL.")) {
      if (!meta.experiments.empty())
        meta.experiments.back().shells.emplace_back();
    }
  }
}

} // namespace pdb_impl

inline int remark_number(const std::string& remark) {
  if (remark.size() > 11)
    return pdb_impl::read_int(remark.c_str() + 7, 3);
  return 0;
}

inline void read_metadata_from_remarks(Structure& st) {
  std::string* possibly_unfinished_remark3 = nullptr;
  std::string* cr_desc = nullptr;
  for (const std::string& remark : st.raw_remarks)
    switch (remark_number(remark)) {
      case 3:
        pdb_impl::read_remark3_line(remark.c_str(), st.meta,
                                    possibly_unfinished_remark3);
        break;
      case 200:
      case 230:
      case 240:
        pdb_impl::read_remark_200_230_240(remark.c_str(), st.meta, cr_desc);
        break;
      case 300:
        if (!st.meta.remark_300_detail.empty()) {
          st.meta.remark_300_detail += '\n';
          st.meta.remark_300_detail += rtrim_str(remark.substr(11));
        } else if (remark.compare(11, 7, "REMARK:") == 0) {
          st.meta.remark_300_detail = trim_str(remark.substr(18));
        }
        break;
    }
}

// Returns operations corresponding to 1555, 2555, ... N555
inline
std::vector<Op> read_remark_290(const std::vector<std::string>& raw_remarks) {
  std::vector<Op> ops;
  // we only check triplet notation:
  // REMARK 290     NNNMMM   OPERATOR
  // REMARK 290       1555   X,Y,Z
  for (const std::string& remark : raw_remarks)
    if (remark_number(remark) == 290 && remark.size() > 25 &&
        std::memcmp(&remark[10], "     ", 5) == 0 &&
        std::memcmp(&remark[18], "555   ", 6) == 0) {
      if (pdb_impl::read_int(remark.c_str() + 15, 3) != (int)ops.size() + 1)
        fail("Symmetry operators not in order?: " + remark);
      Op op = parse_triplet(pdb_impl::read_string(remark.c_str() + 24, 56));
      ops.push_back(op);
    }
  return ops;
}

} // namespace gemmi
#endif
