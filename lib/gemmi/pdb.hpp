// Copyright 2017 Global Phasing Ltd.
//
// Read PDB file format and store it in Structure.
//
// Based on the format spec:
// https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
// + support for two-character chain IDs (columns 21 and 22)
// + read segment ID (columns 73-76)
// + read hybrid-36 serial numbers (http://cci.lbl.gov/hybrid_36/)
// + hybrid-36 sequence id for sequences longer than 9999 (no such examples)
// + allow for longer REMARK lines (up to 120 characters)

#ifndef GEMMI_PDB_HPP_
#define GEMMI_PDB_HPP_

#include <algorithm>  // for swap
#include <cctype>     // for isalpha
#include <cstdio>     // for FILE, size_t
#include <cstdlib>    // for strtol
#include <cstring>    // for memcpy, strstr, strchr, strcmp
#include <map>        // for map
#include <string>     // for string
#include <vector>     // for vector
#include <unordered_map>

#include "atox.hpp"     // for string_to_int
#include "atof.hpp"     // for fast_from_chars
#include "fail.hpp"     // for fail
#include "fileutil.hpp" // for path_basename, file_open
#include "input.hpp"    // for FileStream
#include "model.hpp"
#include "polyheur.hpp" // for assign_subchains
#include "util.hpp"

namespace gemmi {

namespace pdb_impl {

inline int read_int(const char* p, int field_length) {
  return string_to_int(p, false, field_length);
}

template<int N> int read_base36(const char* p) {
  char zstr[N+1] = {0};
  std::memcpy(zstr, p, N);
  return std::strtol(zstr, nullptr, 36);
}

inline double read_double(const char* p, int field_length) {
  double d = 0.;
  // we don't check for errors here
  fast_from_chars(p, p + field_length, d);
  return d;
}

inline std::string read_string(const char* p, int field_length) {
  // left trim
  while (field_length != 0 && is_space(*p)) {
    ++p;
    --field_length;
  }
  // EOL/EOF ends the string
  for (int i = 0; i < field_length; ++i)
    if (p[i] == '\n' || p[i] == '\r' || p[i] == '\0') {
      field_length = i;
      break;
    }
  // right trim
  while (field_length != 0 && is_space(p[field_length-1]))
    --field_length;
  return std::string(p, field_length);
}

// Compare the first 4 letters of s, ignoring case, with uppercase record.
// Both args must have at least 3+1 chars. ' ' and NUL are equivalent in s.
inline bool is_record_type(const char* s, const char* record) {
  return ialpha4_id(s) == ialpha4_id(record);
}
// for record "TER": "TER ", TER\n, TER\r, TER\t match, TERE, TER1 don't
inline bool is_record_type3(const char* s, const char* record) {
  return (ialpha4_id(s) & ~0xf) == ialpha4_id(record);
}

// The standard charge format is 2+, but some files have +2.
inline signed char read_charge(char digit, char sign) {
  if (sign == ' ' && digit == ' ')  // by far the most common case
    return 0;
  if (sign >= '0' && sign <= '9')
    std::swap(digit, sign);
  if (digit >= '0' && digit <= '9') {
    if (sign != '+' && sign != '-' && sign != '\0' && !is_space(sign))
      fail("Wrong format for charge: " +
           std::string(1, digit) + std::string(1, sign));
    return (digit - '0') * (sign == '-' ? -1 : 1);
  }
  // if we are here the field should be blank, but maybe better not to check
  return 0;
}

inline int read_matrix(Transform& t, char* line, size_t len) {
  if (len < 46)
    return 0;
  char n = line[5] - '0';
  if (n >= 1 && n <= 3) {
    t.mat[n-1][0] = read_double(line+10, 10);
    t.mat[n-1][1] = read_double(line+20, 10);
    t.mat[n-1][2] = read_double(line+30, 10);
    t.vec.at(n-1) = read_double(line+45, 10);
  }
  return n;
}

inline SeqId read_seq_id(const char* str) {
  SeqId seqid;
  if (str[4] != '\r' && str[4] != '\n')
    seqid.icode = str[4];
  // We support hybrid-36 extension, although it is never used in practice
  // as 9999 residues per chain are enough.
  if (str[0] < 'A') {
    for (int i = 4; i != 0; --i, ++str)
      if (!is_space(*str)) {
        seqid.num = read_int(str, i);
        break;
      }
  } else {
    seqid.num = read_base36<4>(str) - 466560 + 10000;
  }
  return seqid;
}

inline ResidueId read_res_id(const char* seq_id, const char* name) {
  return {read_seq_id(seq_id), {}, read_string(name, 3)};
}

inline char read_altloc(char c) { return c == ' ' ? '\0' : c; }

inline int read_serial(const char* ptr) {
  return ptr[0] < 'A' ? read_int(ptr, 5)
                      : read_base36<5>(ptr) - 16796160 + 100000;
}

// "28-MAR-07" -> "2007-03-28"
// (we also accept less standard format "28-Mar-2007" as used by BUSTER)
// We do not check if the date is correct.
// The returned value is one of:
//   DDDD-DD-DD - possibly correct date,
//   DDDD-xx-DD - unrecognized month,
//   empty string - the digits were not there.
inline std::string pdb_date_format_to_iso(const std::string& date) {
  const char months[] = "JAN01FEB02MAR03APR04MAY05JUN06"
                        "JUL07AUG08SEP09OCT10NOV11DEC122222";
  if (date.size() < 9 || !is_digit(date[0]) || !is_digit(date[1]) ||
                         !is_digit(date[7]) || !is_digit(date[8]))
    return std::string();
  std::string iso = "xxxx-xx-xx";
  if (date.size() >= 11 && is_digit(date[9]) && is_digit(date[10])) {
    std::memcpy(&iso[0], &date[7], 4);
  } else {
    std::memcpy(&iso[0], (date[7] > '6' ? "19" : "20"), 2);
    std::memcpy(&iso[2], &date[7], 2);
  }
  char month[4] = {alpha_up(date[3]), alpha_up(date[4]), alpha_up(date[5]), '\0'};
  if (const char* m = std::strstr(months, month))
    std::memcpy(&iso[5], m + 3, 2);
  std::memcpy(&iso[8], &date[0], 2);
  return iso;
}

inline Asu compare_link_symops(const std::string& record) {
  if (record.size() < 72)
    return Asu::Any;  // it could be interpreted as Same
  if (read_string(&record[59], 6) == read_string(&record[66], 6))
    return Asu::Same;
  return Asu::Different;
}

// Atom name and altloc are not provided in the SSBOND record.
// Usually it is SG (cysteine), but other disulfide bonds are also possible.
// If it's not SG, we pick the first sulfur atom in the residue.
inline void complete_ssbond_atom(AtomAddress& ad, const Model& mdl) {
  ad.atom_name = "SG";
  const_CRA cra = mdl.find_cra(ad);
  if (cra.residue && (!cra.atom || cra.atom->element != El::S))
    if (const Atom* a = cra.residue->find_by_element(El::S)) {
      ad.atom_name = a->name;
      ad.altloc = a->altloc;
    }
}

inline
void process_conn(Structure& st, const std::vector<std::string>& conn_records) {
  int disulf_count = 0;
  int covale_count = 0;
  int metalc_count = 0;
  for (const std::string& record : conn_records) {
    if (record[0] == 'S' || record[0] == 's') { // SSBOND
      if (record.length() < 32)
        continue;
      Connection c;
      c.name = "disulf" + std::to_string(++disulf_count);
      c.type = Connection::Disulf;
      const char* r = record.c_str();
      c.partner1.chain_name = read_string(r + 14, 2);
      c.partner1.res_id = read_res_id(r + 17, r + 11);
      c.partner2.chain_name = read_string(r + 28, 2);
      char res_id2[5] = {' ', ' ', ' ', ' ', ' '};
      std::memcpy(res_id2, r + 31, std::min((size_t)5, record.length() - 31));
      c.partner2.res_id = read_res_id(res_id2, r + 25);
      c.asu = compare_link_symops(record);
      if (record.length() > 73)
        c.reported_distance = read_double(r + 73, 5);
      complete_ssbond_atom(c.partner1, st.first_model());
      complete_ssbond_atom(c.partner2, st.models[0]);
      st.connections.emplace_back(c);
    } else if (record[0] == 'L' || record[0] == 'l') { // LINK
      if (record.length() < 57)
        continue;
      Connection c;
      // emulating names used in wwPDB mmCIFs (covaleN and metalcN)
      if (is_metal(find_element(&record[12])) ||
          is_metal(find_element(&record[42]))) {
        c.name = "metalc" + std::to_string(++metalc_count);
        c.type = Connection::MetalC;
      } else {
        c.name = "covale" + std::to_string(++covale_count);
        c.type = Connection::Covale;
      }
      for (int i : {0, 1}) {
        const char* t = record.c_str() + 30 * i;
        AtomAddress& ad = (i == 0 ? c.partner1 : c.partner2);
        ad.chain_name = read_string(t + 20, 2);
        ad.res_id = read_res_id(t + 22, t + 17);
        ad.atom_name = read_string(t + 12, 4);
        ad.altloc = read_altloc(t[16]);
      }
      c.asu = compare_link_symops(record);
      if (record.length() > 73) {
        if (record[4] == 'R')
          c.link_id = read_string(&record[72], 8);
        else
          c.reported_distance = read_double(&record[73], 5);
      }
      st.connections.emplace_back(c);
    } else if (record[0] == 'C' || record[0] == 'c') { // CISPEP
      if (record.length() < 22)
        continue;
      const char* r = record.c_str();
      std::string cname = read_string(r + 14, 2);
      ResidueId rid = read_res_id(r + 17, r + 11);
      for (Model& model : st.models)
        if (Residue* res = model.find_residue(cname, rid))
          res->is_cis = true;
    }
  }
}

template<size_t N>
inline bool same_str(const std::string& s, const char (&literal)[N]) {
  return s.size() == N - 1 && std::strcmp(s.c_str(), literal) == 0;
}

template<typename Input>
Structure read_pdb_from_input(Input&& infile, const std::string& source,
                              const PdbReadOptions& options) {
  using namespace pdb_impl;
  int line_num = 0;
  auto wrong = [&line_num](const std::string& msg) {
    fail("Problem in line " + std::to_string(line_num) + ": " + msg);
  };
  Structure st;
  st.input_format = CoorFormat::Pdb;
  st.name = path_basename(source, {".gz", ".pdb"});
  std::vector<std::string> conn_records;
  st.models.emplace_back("1");
  Model *model = &st.models.back();
  Chain *chain = nullptr;
  Residue *resi = nullptr;
  char line[122] = {0};
  int max_line_length = options.max_line_length;
  if (max_line_length <= 0 || max_line_length > 120)
    max_line_length = 120;
  bool after_ter = false;
  Transform matrix;
  std::unordered_map<ResidueId, int> resmap;
  while (size_t len = copy_line_from_stream(line, max_line_length+1, infile)) {
    ++line_num;
    if (is_record_type(line, "ATOM") || is_record_type(line, "HETATM")) {
      if (len < 66)
        wrong("The line is too short to be correct:\n" + std::string(line));
      std::string chain_name = read_string(line+20, 2);
      ResidueId rid = read_res_id(line+22, line+17);

      if (!chain || chain_name != chain->name) {
        if (!model) {
          // MD trajectories may have frames separated by ENDMDL without MODEL.
          std::string name = std::to_string(st.models.size() + 1);
          if (st.find_model(name))
            wrong("ATOM/HETATM between models");
          st.models.emplace_back(name);
          model = &st.models.back();
        }
        const Chain* prev_part = model->find_chain(chain_name);
        after_ter = prev_part &&
                    prev_part->residues[0].entity_type == EntityType::Polymer;
        model->chains.emplace_back(chain_name);
        chain = &model->chains.back();
        resmap.clear();
        resi = nullptr;
      }
      // Non-standard but widely used 4-character segment identifier.
      // Left-justified, and may include a space in the middle.
      // The segment may be a portion of a chain or a complete chain.
      if (len > 72)
        rid.segment = read_string(line+72, 4);
      if (!resi || !resi->matches(rid)) {
        auto it = resmap.find(rid);
        // In normal PDB files it is fast enough to use
        // resi = chain->find_residue(rid);
        // but in pseudo-PDB files (such as MD files where millions
        // of residues are in the same "chain") it is too slow.
        if (it == resmap.end()) {
          resmap.emplace(rid, (int) chain->residues.size());
          chain->residues.emplace_back(rid);
          resi = &chain->residues.back();

          resi->het_flag = line[0] & ~0x20;
          if (after_ter)
            resi->entity_type = resi->is_water() ? EntityType::Water
                                                 : EntityType::NonPolymer;
        } else {
          resi = &chain->residues[it->second];
        }
      }

      Atom atom;
      atom.serial = read_serial(line+6);
      atom.name = read_string(line+12, 4);
      atom.altloc = read_altloc(line[16]);
      atom.pos.x = read_double(line+30, 8);
      atom.pos.y = read_double(line+38, 8);
      atom.pos.z = read_double(line+46, 8);
      atom.occ = (float) read_double(line+54, 6);
      atom.b_iso = (float) read_double(line+60, 6);
      if (len > 76 && (std::isalpha(line[76]) || std::isalpha(line[77])))
        atom.element = Element(line + 76);
      // Atom names HXXX are ambiguous, but Hg, He, Hf, Ho and Hs (almost)
      // never have 4-character names, so H is assumed.
      else if (alpha_up(line[12]) == 'H' && line[15] != ' ')
        atom.element = El::H;
      else
        atom.element = Element(line + 12);
      atom.charge = (len > 78 ? read_charge(line[78], line[79]) : 0);
      resi->atoms.emplace_back(atom);

    } else if (is_record_type(line, "ANISOU")) {
      if (!model || !chain || !resi || resi->atoms.empty())
        wrong("ANISOU record not directly after ATOM/HETATM.");
      // We assume that ANISOU refers to the last atom.
      // Can it not be the case?
      Atom &atom = resi->atoms.back();
      if (atom.aniso.u11 != 0.)
        wrong("Duplicated ANISOU record or not directly after ATOM/HETATM.");
      atom.aniso.u11 = read_int(line+28, 7) * 1e-4f;
      atom.aniso.u22 = read_int(line+35, 7) * 1e-4f;
      atom.aniso.u33 = read_int(line+42, 7) * 1e-4f;
      atom.aniso.u12 = read_int(line+49, 7) * 1e-4f;
      atom.aniso.u13 = read_int(line+56, 7) * 1e-4f;
      atom.aniso.u23 = read_int(line+63, 7) * 1e-4f;

    } else if (is_record_type(line, "REMARK")) {
      if (line[len-1] == '\n')
        --len;
      if (line[len-1] == '\r')
        --len;
      st.raw_remarks.emplace_back(line, line+len);
      if (len <= 11)
        continue;
      int num = read_int(line + 7, 3);
      // By default, we only look for resolution and REMARK 350.
      // Other parsing of remarks is in interpret_remarks().
      if (num == 2) {
        if (st.resolution == 0.0 && std::strstr(line, "ANGSTROM"))
          st.resolution = read_double(line + 23, 7);
      } else if (num == 3) {
        if (st.resolution == 0.0 &&
            std::strstr(line, "RESOLUTION RANGE HIGH (ANGSTROMS)"))
          if (const char* colon = std::strchr(line + 44, ':'))
            st.resolution = fast_atof(colon + 1);
      } else if (num == 350) {
        const char* colon = std::strchr(line+11, ':');
        if (colon == line+22 && starts_with(line+11, "BIOMOLECULE")) {
          st.assemblies.emplace_back(read_string(line+23, 20));
          continue;
        }
        if (st.assemblies.empty())
          continue;
        Assembly& assembly = st.assemblies.back();
        if (starts_with(line+11, "  BIOMT")) {
          if (read_matrix(matrix, line+13, len-13) == 3)
            if (!assembly.generators.empty()) {
              auto& opers = assembly.generators.back().operators;
              opers.emplace_back();
              opers.back().name = read_string(line+20, 3);
              opers.back().transform = matrix;
              matrix.set_identity();
            }
#define CHECK(cpos, text) (colon == line+cpos && starts_with(line+11, text))
        } else if (CHECK(44, "AUTHOR DETERMINED")) {
          assembly.author_determined = true;
          assembly.oligomeric_details = read_string(line+45, 35);
        } else if (CHECK(51, "SOFTWARE DETERMINED")) {
          assembly.software_determined = true;
          assembly.oligomeric_details = read_string(line+52, 28);
        } else if (CHECK(24, "SOFTWARE USED")) {
          assembly.software_name = read_string(line+25, 55);
        } else if (CHECK(36, "TOTAL BURIED SURFACE AREA")) {
          assembly.absa = read_double(line+37, 12);
        } else if (CHECK(38, "SURFACE AREA OF THE COMPLEX")) {
          assembly.ssa = read_double(line+39, 12);
        } else if (CHECK(40, "CHANGE IN SOLVENT FREE ENERGY")) {
          assembly.more = read_double(line+41, 12);
        } else if (CHECK(40, "APPLY THE FOLLOWING TO CHAINS") ||
                   CHECK(40, "                   AND CHAINS")) {
          if (line[11] == 'A') // first line - APPLY ...
            assembly.generators.emplace_back();
          else if (assembly.generators.empty())
            continue;
          split_str_into_multi(read_string(line+41, 39), ", ",
                               assembly.generators.back().chains);
        }
#undef CHECK
      }

    } else if (is_record_type(line, "CONECT")) {
      // ignore for now

    } else if (is_record_type(line, "SEQRES")) {
      std::string chain_name = read_string(line+10, 2);
      Entity& ent = impl::find_or_add(st.entities, chain_name);
      ent.entity_type = EntityType::Polymer;
      for (int i = 19; i < 68; i += 4) {
        std::string res_name = read_string(line+i, 3);
        if (!res_name.empty())
          ent.full_sequence.emplace_back(res_name);
      }

    } else if (is_record_type(line, "DBREF")) { // DBREF or DBREF1 or DBREF2
      std::string chain_name = read_string(line+11, 2);
      Entity& ent = impl::find_or_add(st.entities, chain_name);
      if (line[5] == ' ' || line[5] == '1')
        ent.dbrefs.emplace_back();
      else if (ent.dbrefs.empty()) // DBREF2 without DBREF1?
        continue;
      Entity::DbRef& dbref = ent.dbrefs.back();
      if (line[5] == ' ' || line[5] == '1') {
        dbref.seq_begin = read_seq_id(line+14);
        dbref.seq_end = read_seq_id(line+20);
        dbref.db_name = read_string(line+26, 6);
        if (line[5] == ' ') {
          dbref.accession_code = read_string(line+33, 8);
          dbref.id_code = read_string(line+42, 12);
          dbref.db_begin.num = read_int(line+55, 5);
          dbref.db_begin.icode = line[60];
          dbref.db_end.num = read_int(line+62, 5);
          dbref.db_end.icode = line[67];
        } else {  // line[5] == '1'
          dbref.id_code = read_string(line+47, 20);
        }
      } else if (line[5] == '2') {
        dbref.accession_code = read_string(line+18, 22);
        dbref.db_begin.num = read_int(line+45, 10);
        dbref.db_end.num = read_int(line+57, 10);
      }
    } else if (is_record_type(line, "HEADER")) {
      if (len > 50)
        st.info["_struct_keywords.pdbx_keywords"] = rtrim_str(std::string(line+10, 40));
      if (len > 59) { // date in PDB has format 28-MAR-07
        std::string date = pdb_date_format_to_iso(std::string(line+50, 9));
        if (!date.empty())
          st.info["_pdbx_database_status.recvd_initial_deposition_date"] = date;
      }
      if (len > 66) {
        std::string entry_id = rtrim_str(std::string(line+62, 4));
        if (!entry_id.empty())
          st.info["_entry.id"] = entry_id;
      }
    } else if (is_record_type(line, "TITLE")) {
      if (len > 10)
        st.info["_struct.title"] += rtrim_str(std::string(line+10, len-10-1));

    } else if (is_record_type(line, "KEYWDS")) {
      if (len > 10)
        st.info["_struct_keywords.text"] += rtrim_str(std::string(line+10, len-10-1));

    } else if (is_record_type(line, "EXPDTA")) {
      if (len > 10)
        st.info["_exptl.method"] += trim_str(std::string(line+10, len-10-1));

    } else if (is_record_type(line, "CRYST1")) {
      if (len > 54)
        st.cell.set(read_double(line+6, 9),
                    read_double(line+15, 9),
                    read_double(line+24, 9),
                    read_double(line+33, 7),
                    read_double(line+40, 7),
                    read_double(line+47, 7));
      if (len > 56)
        st.spacegroup_hm = read_string(line+55, 11);
      if (len > 67) {
        std::string z = read_string(line+66, 4);
        if (!z.empty())
          st.info["_cell.Z_PDB"] = z;
      }
    } else if (is_record_type(line, "MTRIXn")) {
      if (read_matrix(matrix, line, len) == 3) {
        std::string id = read_string(line+7, 3);
        if (matrix.is_identity()) {
          // store only ID that will be used when writing to file
          st.info["_struct_ncs_oper.id"] = id;
        } else {
          bool given = len > 59 && line[59] == '1';
          st.ncs.push_back({id, given, matrix});
          matrix.set_identity();
        }
      }
    } else if (is_record_type(line, "MODEL")) {
      if (model && chain)
        wrong("MODEL without ENDMDL?");
      std::string name = std::to_string(read_int(line+10, 4));
      model = &st.find_or_add_model(name);
      if (!model->chains.empty())
        wrong("duplicate MODEL number: " + name);
      chain = nullptr;

    } else if (is_record_type(line, "ENDMDL")) {
      model = nullptr;
      chain = nullptr;

    } else if (is_record_type3(line, "TER")) { // finishes polymer chains
      // we don't expect more than one TER record in one chain
      if (!chain || after_ter)
        continue;
      if (options.split_chain_on_ter) {
        chain = nullptr;
        // split_chain_on_ter is used for AMBER files that can have TER records
        // in various places. So in such case TER doesn't imply entity_type.
        continue;
      }
      for (Residue& res : chain->residues)
        res.entity_type = EntityType::Polymer;
      after_ter = true;
    } else if (is_record_type(line, "SCALEn")) {
      if (read_matrix(matrix, line, len) == 3) {
        st.cell.set_matrices_from_fract(matrix);
        matrix.set_identity();
      }

    } else if (is_record_type(line, "ORIGX")) {
      st.has_origx = true;
      read_matrix(st.origx, line, len);

    } else if (is_record_type(line, "HELIX")) {
      if (len < 40)
        continue;
      Helix helix;
      helix.start.chain_name = read_string(line+18, 2);
      helix.start.res_id = read_res_id(line+21, line+15);
      helix.end.chain_name = read_string(line+30, 2);
      helix.end.res_id = read_res_id(line+33, line+27);
      helix.set_helix_class_as_int(read_int(line+38, 2));
      if (len > 72)
        helix.length = read_int(line+72, 5);
      st.helices.emplace_back(helix);

    } else if (is_record_type(line, "SHEET")) {
      if (len < 40)
        continue;
      std::string sheet_id = read_string(line+11, 3);
      Sheet& sheet = impl::find_or_add(st.sheets, sheet_id);
      sheet.strands.emplace_back();
      Sheet::Strand& strand = sheet.strands.back();
      strand.start.chain_name = read_string(line+20, 2);
      strand.start.res_id = read_res_id(line+22, line+17);
      strand.end.chain_name = read_string(line+31, 2);
      strand.end.res_id = read_res_id(line+33, line+28);
      strand.sense = read_int(line+38, 2);
      if (len > 67) {
        // the SHEET record has no altloc for atoms of hydrogen bond
        strand.hbond_atom2.atom_name = read_string(line+41, 4);
        strand.hbond_atom2.chain_name = read_string(line+48, 2);
        strand.hbond_atom2.res_id = read_res_id(line+50, line+45);
        strand.hbond_atom1.atom_name = read_string(line+56, 4);
        strand.hbond_atom1.chain_name = read_string(line+63, 2);
        strand.hbond_atom1.res_id = read_res_id(line+65, line+60);
      }

    } else if (is_record_type(line, "SSBOND") ||
               is_record_type(line, "LINK") ||
               is_record_type(line, "CISPEP")) {
      conn_records.emplace_back(line);

    } else if (is_record_type3(line, "END")) {
      break;
    } else if (is_record_type(line, "data")) {
      if (line[4] == '_' && model && model->chains.empty())
        fail("Incorrect file format (perhaps it is cif not pdb?): " + source);
    }
  }

  for (Model& mod : st.models)
    for (Chain& ch : mod.chains)
      if (ch.residues[0].entity_type != EntityType::Unknown)
        assign_subchain_names(ch);

  for (Chain& ch : st.models[0].chains)
    if (Entity* entity = st.get_entity(ch.name))
      if (auto polymer = ch.get_polymer())
        entity->subchains.emplace_back(polymer.subchain_id());

  st.setup_cell_images();

  process_conn(st, conn_records);

  return st;
}

}  // namespace pdb_impl

inline Structure read_pdb_file(const std::string& path,
                               PdbReadOptions options=PdbReadOptions()) {
  auto f = file_open(path.c_str(), "rb");
  return pdb_impl::read_pdb_from_input(FileStream{f.get()}, path, options);
}

inline Structure read_pdb_from_memory(const char* data, size_t size,
                                      const std::string& name,
                                      PdbReadOptions options=PdbReadOptions()) {
  return pdb_impl::read_pdb_from_input(MemoryStream(data, size), name, options);
}

inline Structure read_pdb_string(const std::string& str,
                                 const std::string& name,
                                 PdbReadOptions options=PdbReadOptions()) {
  return read_pdb_from_memory(str.c_str(), str.length(), name, options);
}

// A function for transparent reading of stdin and/or gzipped files.
template<typename T>
inline Structure read_pdb(T&& input, PdbReadOptions options=PdbReadOptions()) {
  if (input.is_stdin())
    return pdb_impl::read_pdb_from_input(FileStream{stdin}, "stdin", options);
  if (input.is_compressed())
    return pdb_impl::read_pdb_from_input(input.get_uncompressing_stream(),
                                         input.path(), options);
  return read_pdb_file(input.path(), options);
}

} // namespace gemmi
#endif
