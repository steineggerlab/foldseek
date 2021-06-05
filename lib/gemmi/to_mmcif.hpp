// Copyright 2017 Global Phasing Ltd.
//
// Create cif::Document (for PDBx/mmCIF file) from Structure.

#ifndef GEMMI_TO_MMCIF_HPP_
#define GEMMI_TO_MMCIF_HPP_

#include <set>
#include "model.hpp"
#include "cifdoc.hpp"

namespace gemmi {

struct MmcifOutputGroups {
  bool atoms:1;
  bool block_name:1;
  bool entry:1;
  bool database_status:1;
  bool cell:1;
  bool symmetry:1;
  bool entity:1;
  bool entity_poly:1;
  bool struct_ref:1;
  bool chem_comp:1;
  bool exptl:1;
  bool diffrn:1;
  bool reflns:1;
  bool refine:1;
  bool title_keywords:1;
  bool ncs:1;
  bool struct_asym:1;
  bool origx:1;
  bool struct_conf:1;
  bool struct_sheet:1;
  bool struct_biol:1;
  bool assembly:1;
  bool conn:1;
  bool cis:1;
  bool scale:1;
  bool atom_type:1;
  bool entity_poly_seq:1;
  bool tls:1;
  bool software:1;
  bool group_pdb:1;  // include _atom_site.group_PDB

  explicit MmcifOutputGroups(bool all)
    : atoms(all), block_name(all), entry(all), database_status(all),
      cell(all), symmetry(all), entity(all),
      entity_poly(false),  // see the comment under "if (groups.entity_poly)"
      struct_ref(all), chem_comp(all), exptl(all), diffrn(all),
      reflns(all), refine(all), title_keywords(all), ncs(all),
      struct_asym(all), origx(all), struct_conf(all), struct_sheet(all),
      struct_biol(all), assembly(all), conn(all), cis(all),
      scale(all), atom_type(all), entity_poly_seq(all), tls(all),
      software(all), group_pdb(false) {}
};

void update_mmcif_block(const Structure& st, cif::Block& block, MmcifOutputGroups groups);
cif::Document make_mmcif_document(const Structure& st);
cif::Block make_mmcif_headers(const Structure& st);
void add_minimal_mmcif_data(const Structure& st, cif::Block& block);

// temporarily we use it in crdrst.cpp
namespace impl {
void write_struct_conn(const Structure& st, cif::Block& block);
}

} // namespace gemmi

#ifdef GEMMI_WRITE_IMPLEMENTATION

#include <cassert>
#include <cmath>  // for isnan
#include <string>
#include <utility>  // std::pair
#include "sprintf.hpp"
#include "enumstr.hpp"  // for entity_type_to_string, polymer_type_to_qstring
#include "calculate.hpp"  // for count_atom_sites

namespace gemmi {

namespace impl {

struct ItemSpan {
  ItemSpan(std::vector<cif::Item>& items, const std::string& cat)
      : items_(items), begin_(0), end_(items.size()) {
    cif::assert_tag(cat);
    while (begin_ != items.size() && !has_prefix(items[begin_], cat))
      ++begin_;
    if (begin_ != end_)
      while (end_ - 1 != begin_ && !has_prefix(items[end_-1], cat))
        --end_;
  }

  // cf. Block::set_pair()
  void set_pair(const std::string& tag, const std::string& value) {
    using namespace gemmi::cif;
    assert_tag(tag);
    for (size_t i = begin_; i != end_; ++i) {
      Item& item = items_[i];
      if (item.type == ItemType::Pair && item.pair[0] == tag) {
        item.pair[1] = value;
        return;
      }
      if (item.type == ItemType::Loop && item.loop.find_tag(tag) != -1) {
        item.set_value(Item(tag, value));
        return;
      }
    }
    items_.emplace(items_.begin() + end_, tag, value);
    ++end_;
  }
private:
  std::vector<cif::Item>& items_;
  size_t begin_, end_;

  static bool has_prefix(const cif::Item& item, const std::string& cat) {
    if (item.type == cif::ItemType::Pair)
      return starts_with(item.pair[0], cat);
    if (item.type == cif::ItemType::Loop)
      return !item.loop.tags.empty() && starts_with(item.loop.tags[0], cat);
    return false;
  }
};

inline std::string pdbx_icode(const SeqId& seqid) {
  return std::string(1, seqid.has_icode() ? seqid.icode : '?');
}
inline std::string pdbx_icode(const ResidueId& rid) {
  return pdbx_icode(rid.seqid);
}

inline std::string subchain_or_dot(const Residue& res) {
  return res.subchain.empty() ? "." : cif::quote(res.subchain);
}

inline std::string number_or_dot(double d) {
  return std::isnan(d) ? "." : to_str(d);
}
inline std::string number_or_qmark(double d) {
  return std::isnan(d) ? "?" : to_str(d);
}

// for use with non-negative Metadata fields that use -1 for N/A
inline std::string int_or_dot(int n) {
  return n == -1 ? "." : std::to_string(n);
}
inline std::string int_or_qmark(int n) {
  return n == -1 ? "?" : std::to_string(n);
}

inline std::string string_or_dot(const std::string& s) {
  return s.empty() ? "." : cif::quote(s);
}
inline std::string string_or_qmark(const std::string& s) {
  return s.empty() ? "?" : cif::quote(s);
}

// Quote chain name or entity id if necessary. It is necessary
// only if the chain name is missing, which was OK in the past.
// Here we use '' rather than . or ?.
inline std::string qchain(const std::string& s) {
  return cif::quote(s);
}


inline void add_cif_atoms(const Structure& st, cif::Block& block, bool use_group_pdb) {
  // atom list
  cif::Loop& atom_loop = block.init_mmcif_loop("_atom_site.", {
      "id",
      "type_symbol",
      "label_atom_id",
      "label_alt_id",
      "label_comp_id",
      "label_asym_id",
      "label_entity_id",
      "label_seq_id",
      "pdbx_PDB_ins_code",
      "Cartn_x",
      "Cartn_y",
      "Cartn_z",
      "occupancy",
      "B_iso_or_equiv",
      "pdbx_formal_charge",
      "auth_seq_id",
      "auth_asym_id",
      "pdbx_PDB_model_num"});
  if (use_group_pdb)
    atom_loop.tags.emplace(atom_loop.tags.begin(), "_atom_site.group_PDB");
  bool has_calc_flag = false;
  bool has_tls_group_id = false;
  for (const Model& model : st.models)
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& atom : res.atoms) {
          if (atom.calc_flag != CalcFlag::NotSet)
            has_calc_flag = true;
          if (atom.tls_group_id >= 0)
            has_tls_group_id = true;
        }
  if (has_calc_flag)
    atom_loop.tags.emplace_back("_atom_site.calc_flag");
  if (has_tls_group_id)
    atom_loop.tags.emplace_back("_atom_site.pdbx_tls_group_id");

  std::vector<std::string>& vv = atom_loop.values;
  vv.reserve(count_atom_sites(st) * atom_loop.tags.size());
  std::vector<std::pair<int, const Atom*>> aniso;
  int serial = 0;
  for (const Model& model : st.models) {
    for (const Chain& chain : model.chains) {
      for (const Residue& res : chain.residues) {
        std::string label_seq_id = res.label_seq.str('.');
        std::string auth_seq_id = res.seqid.num.str();
        std::string entity_id(1, '.');
        if (const Entity* ent = gemmi::find_entity(res.subchain, st.entities))
          entity_id = cif::quote(ent->name);
        for (const Atom& atom : res.atoms) {
          if (use_group_pdb)
            vv.emplace_back(res.het_flag != 'H' ? "ATOM" : "HETATM");
          vv.emplace_back(std::to_string(++serial));
          vv.emplace_back(atom.element.uname());
          vv.emplace_back(cif::quote(atom.name));
          vv.emplace_back(1, atom.altloc_or('.'));
          vv.emplace_back(res.name);
          vv.emplace_back(subchain_or_dot(res));
          vv.emplace_back(entity_id);
          vv.emplace_back(label_seq_id);
          vv.emplace_back(pdbx_icode(res));
          vv.emplace_back(to_str(atom.pos.x));
          vv.emplace_back(to_str(atom.pos.y));
          vv.emplace_back(to_str(atom.pos.z));
          vv.emplace_back(to_str(atom.occ));
          vv.emplace_back(to_str(atom.b_iso));
          vv.emplace_back(atom.charge == 0 ? "?" : std::to_string(atom.charge));
          vv.emplace_back(auth_seq_id);
          vv.emplace_back(impl::qchain(chain.name));
          vv.emplace_back(model.name);
          if (has_calc_flag)
            vv.emplace_back(&".\0d\0c\0dum"[2 * (int) atom.calc_flag]);
          if (has_tls_group_id)
            vv.emplace_back(int_or_qmark(atom.tls_group_id));
          if (atom.aniso.nonzero())
            aniso.emplace_back(serial, &atom);
        }
      }
    }
  }
  if (aniso.empty()) {
    block.find_mmcif_category("_atom_site_anisotrop.").erase();
  } else {
    cif::Loop& aniso_loop = block.init_mmcif_loop("_atom_site_anisotrop.", {
                                    "id", "U[1][1]", "U[2][2]", "U[3][3]",
                                    "U[1][2]", "U[1][3]", "U[2][3]"});
    std::vector<std::string>& aniso_val = aniso_loop.values;
    aniso_val.reserve(aniso_loop.tags.size() * aniso.size());
    for (const auto& a : aniso) {
      aniso_val.emplace_back(std::to_string(a.first));
      aniso_val.emplace_back(to_str(a.second->aniso.u11));
      aniso_val.emplace_back(to_str(a.second->aniso.u22));
      aniso_val.emplace_back(to_str(a.second->aniso.u33));
      aniso_val.emplace_back(to_str(a.second->aniso.u12));
      aniso_val.emplace_back(to_str(a.second->aniso.u13));
      aniso_val.emplace_back(to_str(a.second->aniso.u23));
    }
  }
}

// the names are: monomeric, dimeric, ...meric, 21-meric, 22-meric, ...
int xmeric_to_number(const std::string& oligomeric) {
  static const char names[20][10] = {
    "mono", "di", "tri", "tetra", "penta",
    "hexa", "hepta", "octa", "nona", "deca",
    "undeca", "dodeca", "trideca", "tetradeca", "pentadeca",
    "hexadeca", "heptadeca", "octadeca", "nonadeca", "eicosa"
  };
  size_t len = oligomeric.length();
  const char* p = oligomeric.c_str();
  for (int i = 0; i != 20; ++i)
    if (len == std::strlen(names[i]) + 5 && strncmp(p, names[i], len-5) == 0)
      return i + 1;
  return no_sign_atoi(p);
}

void write_assemblies(const Structure& st, cif::Block& block) {
  block.items.reserve(block.items.size() + 4); // avoid re-allocation
  cif::Loop& a_loop = block.init_mmcif_loop("_pdbx_struct_assembly.",
      {"id", "details", "method_details",
       "oligomeric_details", "oligomeric_count"});
  cif::Loop& prop_loop = block.init_mmcif_loop("_pdbx_struct_assembly_prop.",
      {"biol_id", "type", "value"});
  cif::Loop& gen_loop = block.init_mmcif_loop("_pdbx_struct_assembly_gen.",
      {"assembly_id", "oper_expression", "asym_id_list"});
  cif::Loop& oper_loop = block.init_mmcif_loop("_pdbx_struct_oper_list.",
      {"id", "type",
       "matrix[1][1]", "matrix[1][2]", "matrix[1][3]", "vector[1]",
       "matrix[2][1]", "matrix[2][2]", "matrix[2][3]", "vector[2]",
       "matrix[3][1]", "matrix[3][2]", "matrix[3][3]", "vector[3]"});
  std::vector<const Assembly::Operator*> distinct_oper;
  for (const Assembly& as : st.assemblies) {
    std::string how_defined = "?";
    if (as.author_determined && as.software_determined)
      how_defined = "author_and_software_defined_assembly";
    else if (as.author_determined)
      how_defined = "author_defined_assembly";
    else if (as.software_determined)
      how_defined = "software_defined_assembly";
    else if (as.special_kind == Assembly::SpecialKind::CompleteIcosahedral)
      how_defined = "'complete icosahedral assembly'";
    else if (as.special_kind == Assembly::SpecialKind::RepresentativeHelical)
      how_defined = "'representative helical assembly'";
    else if (as.special_kind == Assembly::SpecialKind::CompletePoint)
      how_defined = "'complete point assembly'";
    std::string oligomer = to_lower(as.oligomeric_details);
    int nmer = as.oligomeric_count != 0 ? as.oligomeric_count
                                        : xmeric_to_number(oligomer);
    // _pdbx_struct_assembly
    a_loop.add_row({as.name,
                    how_defined,
                    impl::string_or_qmark(as.software_name),
                    impl::string_or_qmark(oligomer),
                    nmer == 0 ? "?" : std::to_string(nmer)});

    // _pdbx_struct_assembly_prop
    if (!std::isnan(as.absa))
      prop_loop.add_row({as.name, "'ABSA (A^2)'", to_str(as.absa)});
    if (!std::isnan(as.ssa))
      prop_loop.add_row({as.name, "'SSA (A^2)'", to_str(as.ssa)});
    if (!std::isnan(as.more))
      prop_loop.add_row({as.name, "MORE", to_str(as.more)});

    // _pdbx_struct_assembly_gen and _pdbx_struct_oper_list
    for (const Assembly::Gen& gen : as.generators) {
      std::string subchain_str;
      for (const std::string& name : gen.subchains)
        string_append_sep(subchain_str, ',', name);
      if (subchain_str.empty()) // chain names to subchain names
        for (const Chain& chain : st.models[0].chains)
          if (in_vector(chain.name, gen.chains))
            for (const auto& sub : chain.subchains())
              string_append_sep(subchain_str, ',', sub.front().subchain);
      std::string oper_str;
      for (const Assembly::Operator& oper : gen.operators) {
        size_t k = 0;
        for (; k != distinct_oper.size(); ++k)
          if (distinct_oper[k]->transform.approx(oper.transform, 1e-9))
            break;
        string_append_sep(oper_str, ',', std::to_string(k+1));
        if (k != distinct_oper.size())
          continue;
        distinct_oper.emplace_back(&oper);
        oper_loop.values.emplace_back(std::to_string(k+1));
        if (!oper.type.empty()) {
          oper_loop.values.emplace_back(cif::quote(oper.type));
        } else if (oper.transform.is_identity()) {
          oper_loop.values.emplace_back("'identity operation'");
        } else if (as.author_determined || as.software_determined) {
          oper_loop.values.emplace_back("'crystal symmetry operation'");
        } else {
          oper_loop.values.emplace_back(".");
        }
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j)
            oper_loop.values.emplace_back(to_str(oper.transform.mat[i][j]));
          oper_loop.values.emplace_back(to_str(oper.transform.vec.at(i)));
        }
      }
      gen_loop.add_row({as.name,
                        oper_str.empty() ? "." : oper_str,
                        subchain_str.empty() ? "?" : subchain_str});
    }
  }
}

void write_struct_conn(const Structure& st, cif::Block& block) {
  // example:
  // disulf1 disulf A CYS 3  SG ? 3 ? 1_555 A CYS 18 SG ? 18 ?  1_555 ? 2.045
  std::array<bool,(int)Connection::Type::Unknown+1> type_ids{};
  bool use_ccp4_link_id = false;
  for (const Connection& con : st.connections)
    if (!con.link_id.empty())
      use_ccp4_link_id = true;
  cif::Loop& conn_loop = block.init_mmcif_loop("_struct_conn.",
      {"id", "conn_type_id",
       "ptnr1_auth_asym_id", "ptnr1_label_asym_id", "ptnr1_label_comp_id",
       "ptnr1_label_seq_id", "ptnr1_label_atom_id", "pdbx_ptnr1_label_alt_id",
       "ptnr1_auth_seq_id", "pdbx_ptnr1_PDB_ins_code", "ptnr1_symmetry",
       "ptnr2_auth_asym_id", "ptnr2_label_asym_id", "ptnr2_label_comp_id",
       "ptnr2_label_seq_id", "ptnr2_label_atom_id", "pdbx_ptnr2_label_alt_id",
       "ptnr2_auth_seq_id", "pdbx_ptnr2_PDB_ins_code", "ptnr2_symmetry",
       "details", "pdbx_dist_value"});
  if (use_ccp4_link_id)
    conn_loop.tags.push_back("_struct_conn.ccp4_link_id");
  for (const Connection& con : st.connections) {
    const_CRA cra1 = st.models[0].find_cra(con.partner1, true);
    const_CRA cra2 = st.models[0].find_cra(con.partner2, true);
    if (!cra1.residue || !cra2.residue)
      continue;
    const Atom* at1 = cra1.atom;
    const Atom* at2 = cra2.atom;
    std::string im_pdb_symbol = "?", im_dist_str = "?";
    if (at1 && at2) {
      SymImage im = st.cell.find_nearest_image(at1->pos, at2->pos, con.asu);
      im_pdb_symbol = im.pdb_symbol(true);
      im_dist_str = to_str_prec<4>(im.dist());
    }
    auto& v = conn_loop.values;
    v.emplace_back(impl::string_or_qmark(con.name));        // id
    v.emplace_back(connection_type_to_string(con.type));    // conn_type_id
    v.emplace_back(impl::qchain(con.partner1.chain_name));  // ptnr1_auth_asym_id
    v.emplace_back(subchain_or_dot(*cra1.residue));         // ptnr1_label_asym_id
    v.emplace_back(cra1.residue->name);                     // ptnr1_label_comp_id
    v.emplace_back(cra1.residue->label_seq.str('.'));       // ptnr1_label_seq_id
    v.emplace_back(at1 ? cif::quote(at1->name) : "?");      // ptnr1_label_atom_id
    v.emplace_back(1, at1 ? at1->altloc_or('?') : '?');     // pdbx_ptnr1_label_alt_id
    v.emplace_back(cra1.residue->seqid.num.str());          // ptnr1_auth_seq_id
    v.emplace_back(pdbx_icode(con.partner1.res_id));        // ptnr1_PDB_ins_code
    v.emplace_back("1_555");                                // ptnr1_symmetry
    v.emplace_back(impl::qchain(con.partner2.chain_name));  // ptnr2_auth_asym_id
    v.emplace_back(subchain_or_dot(*cra2.residue));         // ptnr2_label_asym_id
    v.emplace_back(cra2.residue->name);                     // ptnr2_label_comp_id
    v.emplace_back(cra2.residue->label_seq.str('.'));       // ptnr2_label_seq_id
    v.emplace_back(at2 ? cif::quote(at2->name) : "?");      // ptnr2_label_atom_id
    v.emplace_back(1, at2 ? at2->altloc_or('?') : '?');     // pdbx_ptnr2_label_alt_id
    v.emplace_back(cra2.residue->seqid.num.str());          // ptnr2_auth_seq_id
    v.emplace_back(pdbx_icode(con.partner2.res_id));        // ptnr2_PDB_ins_code
    v.emplace_back(im_pdb_symbol);                          // ptnr2_symmetry
    v.emplace_back("?");                                    // details
    v.emplace_back(im_dist_str);                            // pdbx_dist_value
    if (use_ccp4_link_id)
      v.emplace_back(impl::string_or_qmark(con.link_id));   // ccp4_link_id
    type_ids[int(con.type)] = true;
  }

  cif::Loop& type_loop = block.init_mmcif_loop("_struct_conn_type.", {"id"});
  for (int i = 0; i < (int)type_ids.size() - 1; ++i)
    if (type_ids[i])
      type_loop.add_row({connection_type_to_string((Connection::Type)i)});
}

void write_cell_parameters(const UnitCell& cell, ItemSpan& span) {
  span.set_pair("_cell.length_a",    to_str(cell.a));
  span.set_pair("_cell.length_b",    to_str(cell.b));
  span.set_pair("_cell.length_c",    to_str(cell.c));
  span.set_pair("_cell.angle_alpha", to_str(cell.alpha));
  span.set_pair("_cell.angle_beta",  to_str(cell.beta));
  span.set_pair("_cell.angle_gamma", to_str(cell.gamma));
}

void write_ncs_oper(const Structure& st, cif::Block& block) {
  // _struct_ncs_oper (MTRIX)
  if (st.ncs.empty())
    return;
  cif::Loop& ncs_oper = block.init_mmcif_loop("_struct_ncs_oper.",
      {"id", "code",
       "matrix[1][1]", "matrix[1][2]", "matrix[1][3]", "vector[1]",
       "matrix[2][1]", "matrix[2][2]", "matrix[2][3]", "vector[2]",
       "matrix[3][1]", "matrix[3][2]", "matrix[3][3]", "vector[3]"});
  auto add_op = [&ncs_oper](const NcsOp& op) {
    ncs_oper.values.emplace_back(op.id);
    ncs_oper.values.emplace_back(op.given ? "given" : "generate");
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j)
        ncs_oper.values.emplace_back(to_str(op.tr.mat[i][j]));
      ncs_oper.values.emplace_back(to_str(op.tr.vec.at(i)));
    }
  };
  auto identity = st.info.find("_struct_ncs_oper.id");
  if (identity != st.info.end() &&
      !in_vector_f([&](const NcsOp& op) { return op.id == identity->second; }, st.ncs))
    add_op(NcsOp{identity->second, true, {}});
  for (const NcsOp& op : st.ncs)
    add_op(op);
}

bool is_valid_block_name(const std::string& name) {
  return !name.empty() &&
         std::all_of(name.begin(), name.end(), [](char c){ return c >= '!' && c <= '~'; });
}

} // namespace impl

void update_mmcif_block(const Structure& st, cif::Block& block, MmcifOutputGroups groups) {
  using std::to_string;
  if (st.models.empty())
    return;

  if (groups.block_name)
    block.name = impl::is_valid_block_name(st.name) ? st.name : "model";

  auto e_id = st.info.find("_entry.id");
  std::string id = cif::quote(e_id != st.info.end() ? e_id->second : block.name);
  if (groups.entry)
    block.set_pair("_entry.id", id);
  else if (const std::string* val = block.find_value("_entry.id"))
    id = *val;

  if (groups.database_status) {
    auto initial_date =
           st.info.find("_pdbx_database_status.recvd_initial_deposition_date");
    if (initial_date != st.info.end() && !initial_date->second.empty()) {
      impl::ItemSpan span(block.items, "_pdbx_database_status.");
      span.set_pair("_pdbx_database_status.entry_id", id);
      span.set_pair(initial_date->first, initial_date->second);
    }
  }

  if (groups.cell) {
    impl::ItemSpan cell_span(block.items, "_cell.");
    cell_span.set_pair("_cell.entry_id", id);
    impl::write_cell_parameters(st.cell, cell_span);
    auto z_pdb = st.info.find("_cell.Z_PDB");
    if (z_pdb != st.info.end())
      cell_span.set_pair(z_pdb->first, z_pdb->second);
  }

  if (groups.symmetry) {
    impl::ItemSpan span(block.items, "_symmetry.");
    span.set_pair("_symmetry.entry_id", id);
    span.set_pair("_symmetry.space_group_name_H-M",
                   cif::quote(st.spacegroup_hm));
    if (const SpaceGroup* sg = st.find_spacegroup())
      span.set_pair("_symmetry.Int_Tables_number", to_string(sg->number));
  }

  if (groups.entity) {
    cif::Loop& entity_loop = block.init_mmcif_loop("_entity.", {"id", "type"});
    for (const Entity& ent : st.entities)
      entity_loop.add_row({impl::qchain(ent.name),
                           entity_type_to_string(ent.entity_type)});
  }

  if (groups.entity_poly) {
    // In this category we write mandatory tags only.
    // But if such a file is deposited to the wwPDB it causes:
    //  ERROR: Missing '_entity_poly.pdbx_strand_id' item.
    // If pdbx_strand_id was present, OneDep would complain about missing
    // _entity_poly.pdbx_seq_one_letter_code. So just in case, by default
    // MmcifOutputGroups::entity_poly is false.
    cif::Loop& ent_poly_loop = block.init_mmcif_loop("_entity_poly.", {"entity_id", "type"});
    for (const Entity& ent : st.entities)
      if (ent.entity_type == EntityType::Polymer)
        ent_poly_loop.add_row({impl::qchain(ent.name),
                               polymer_type_to_qstring(ent.polymer_type)});
  }

  if (groups.struct_ref) { // _struct_ref, _struct_ref_seq
    block.items.reserve(block.items.size() + 2); // avoid re-allocation
    cif::Loop& ref_loop = block.init_mmcif_loop("_struct_ref.",
                                  {"id", "entity_id", "db_name", "db_code",
                                   "pdbx_db_accession", "pdbx_db_isoform"});
    cif::Loop& seq_loop = block.init_mmcif_loop("_struct_ref_seq.", {
        "align_id", "ref_id", "pdbx_strand_id",
        "seq_align_beg", "seq_align_end",
        "db_align_beg", "db_align_end",
        "pdbx_auth_seq_align_beg", "pdbx_seq_align_beg_ins_code",
        "pdbx_auth_seq_align_end", "pdbx_seq_align_end_ins_code"});
    int counter = 0;
    int counter2 = 0;
    std::map<std::string, std::string> subs = st.models[0].subchain_to_chain();
    for (const Entity& ent : st.entities)
      for (const Entity::DbRef& dbref : ent.dbrefs) {
        ref_loop.add_row({std::to_string(++counter),
                          impl::qchain(ent.name),
                          impl::string_or_dot(dbref.db_name),
                          impl::string_or_dot(dbref.id_code),
                          impl::string_or_qmark(dbref.accession_code),
                          impl::string_or_qmark(dbref.isoform)});
        for (const std::string& subchain : ent.subchains) {
          auto strand_id = subs.find(subchain);
          if (strand_id == subs.end())
            continue;
          // DbRef::label_seq_begin/end (_struct_ref_seq.seq_align_beg/end) is
          // not filled in when reading PDB file, so we check it here.
          Residue::OptionalNum label_begin = dbref.label_seq_begin;
          Residue::OptionalNum label_end = dbref.label_seq_end;
          if (!label_begin || !label_end) {
            ConstResidueSpan span = st.models[0].get_subchain(subchain);
            try {
              label_begin = span.auth_seq_id_to_label(dbref.seq_begin);
              label_end = span.auth_seq_id_to_label(dbref.seq_end);
            } catch (const std::out_of_range&) {}
          }
          seq_loop.add_row({std::to_string(++counter2),
                            std::to_string(counter),
                            strand_id->second,  // pdbx_strand_id
                            label_begin.str(),
                            label_end.str(),
                            dbref.db_begin.num.str(),
                            dbref.db_end.num.str(),
                            dbref.seq_begin.num.str(),
                            impl::pdbx_icode(dbref.seq_begin),
                            dbref.seq_end.num.str(),
                            impl::pdbx_icode(dbref.seq_end)});
        }
      }
  }

  if (groups.chem_comp) {
    std::set<std::string> resnames;
    for (const Model& model : st.models)
      for (const Chain& chain : model.chains)
        for (const Residue& res : chain.residues)
          resnames.insert(res.name);
    for (const Entity& ent : st.entities)
      for (const std::string& item : ent.full_sequence)
        resnames.insert(Entity::first_mon(item));
    cif::Loop& chem_comp_loop = block.init_mmcif_loop("_chem_comp.", {"id", "type"});
    for (const std::string& name : resnames)
      chem_comp_loop.add_row({name, "."});
  }

  if (groups.exptl) {
    // _exptl
    if (!st.meta.experiments.empty()) {
      cif::Loop& loop = block.init_mmcif_loop("_exptl.",
                                              {"entry_id", "method", "crystals_number"});
      for (const ExperimentInfo& exper : st.meta.experiments)
        loop.add_row({id, cif::quote(exper.method),
                      impl::int_or_qmark(exper.number_of_crystals)});
    } else {
      auto exptl_method = st.info.find("_exptl.method");
      if (exptl_method != st.info.end()) {
        cif::Loop& loop = block.init_mmcif_loop("_exptl.", {"entry_id", "method"});
        for (const std::string& m : gemmi::split_str(exptl_method->second, "; "))
          loop.add_row({id, cif::quote(m)});
      }
    }

    // _exptl_crystal
    if (!st.meta.crystals.empty()) {
      cif::Loop& loop = block.init_mmcif_loop("_exptl_crystal.",
                                              {"id", "description"});
      for (const CrystalInfo& cryst : st.meta.crystals)
        loop.add_row({cryst.id, impl::string_or_qmark(cryst.description)});
    }

    // _exptl_crystal_grow
    if (std::any_of(st.meta.crystals.begin(), st.meta.crystals.end(),
          [](const CrystalInfo& c) { return !c.ph_range.empty() || !std::isnan(c.ph); })) {
      cif::Loop& grow_loop = block.init_mmcif_loop("_exptl_crystal_grow.",
                                                   {"crystal_id", "pH", "pdbx_pH_range"});
      for (const CrystalInfo& crystal : st.meta.crystals)
        grow_loop.add_row({cif::quote(crystal.id),
                           impl::number_or_qmark(crystal.ph),
                           impl::string_or_qmark(crystal.ph_range)});
    }
  }

  if (groups.diffrn &&
      std::any_of(st.meta.crystals.begin(), st.meta.crystals.end(),
                  [](const CrystalInfo& c) { return !c.diffractions.empty(); })) {

    cif::Loop& loop = block.init_mmcif_loop("_diffrn.", {"id", "crystal_id", "ambient_temp"});
    for (const CrystalInfo& cryst : st.meta.crystals)
      for (const DiffractionInfo& diffr : cryst.diffractions)
        loop.add_row({diffr.id, cryst.id, impl::number_or_qmark(diffr.temperature)});
    // _diffrn_detector
    cif::Loop& det_loop = block.init_mmcif_loop("_diffrn_detector.",
                                                {"diffrn_id",
                                                 "pdbx_collection_date",
                                                 "detector",
                                                 "type",
                                                 "details"});
    for (const CrystalInfo& cryst : st.meta.crystals)
      for (const DiffractionInfo& diffr : cryst.diffractions)
        det_loop.add_row({diffr.id,
                          impl::string_or_qmark(diffr.collection_date),
                          impl::string_or_qmark(diffr.detector),
                          impl::string_or_qmark(diffr.detector_make),
                          impl::string_or_qmark(diffr.optics)});

    // _diffrn_radiation
    cif::Loop& rad_loop = block.init_mmcif_loop("_diffrn_radiation.",
                                                {"diffrn_id",
                                                 "pdbx_scattering_type",
                                                 "pdbx_monochromatic_or_laue_m_l",
                                                 "monochromator"});
    for (const CrystalInfo& cryst : st.meta.crystals)
      for (const DiffractionInfo& diffr : cryst.diffractions)
        rad_loop.add_row({diffr.id,
                          impl::string_or_qmark(diffr.scattering_type),
                          std::string(1, diffr.mono_or_laue ? diffr.mono_or_laue : '?'),
                          impl::string_or_qmark(diffr.monochromator)});
    // _diffrn_source
    cif::Loop& source_loop = block.init_mmcif_loop("_diffrn_source.",
                                                   {"diffrn_id",
                                                    "source",
                                                    "type",
                                                    "pdbx_synchrotron_site",
                                                    "pdbx_synchrotron_beamline",
                                                    "pdbx_wavelength_list"});
    for (const CrystalInfo& crystal : st.meta.crystals)
      for (const DiffractionInfo& diffr : crystal.diffractions)
        source_loop.add_row({diffr.id,
                             impl::string_or_qmark(diffr.source),
                             impl::string_or_qmark(diffr.source_type),
                             impl::string_or_qmark(diffr.synchrotron),
                             impl::string_or_qmark(diffr.beamline),
                             impl::string_or_qmark(diffr.wavelengths)});
  }

  if (groups.reflns && !st.meta.experiments.empty()) {
    // _reflns
    cif::Loop& loop = block.init_mmcif_loop("_reflns.", {
        "entry_id",
        "pdbx_ordinal",
        "pdbx_diffrn_id",
        "number_obs",
        "d_resolution_high",
        "d_resolution_low",
        "percent_possible_obs",
        "pdbx_redundancy",
        "pdbx_Rmerge_I_obs",
        "pdbx_Rsym_value",
        "pdbx_netI_over_sigmaI",
        /*"B_iso_Wilson_estimate"*/});
    int n = 0;
    for (const ExperimentInfo& exper : st.meta.experiments)
      loop.add_row({id,
                    std::to_string(++n),
                    impl::string_or_dot(join_str(exper.diffraction_ids, ",")),
                    impl::int_or_qmark(exper.unique_reflections),
                    impl::number_or_qmark(exper.reflections.resolution_high),
                    impl::number_or_qmark(exper.reflections.resolution_low),
                    impl::number_or_qmark(exper.reflections.completeness),
                    impl::number_or_qmark(exper.reflections.redundancy),
                    impl::number_or_qmark(exper.reflections.r_merge),
                    impl::number_or_qmark(exper.reflections.r_sym),
                    impl::number_or_qmark(exper.reflections.mean_I_over_sigma),
                    /*impl::number_or_qmark(exper.b_wilson)*/});
    // _reflns_shell
    cif::Loop& shell_loop = block.init_mmcif_loop("_reflns_shell.", {
        "pdbx_ordinal",
        "pdbx_diffrn_id",
        "d_res_high",
        "d_res_low",
        "percent_possible_all",
        "pdbx_redundancy",
        "Rmerge_I_obs",
        "pdbx_Rsym_value",
        "meanI_over_sigI_obs"});
    n = 0;
    for (const ExperimentInfo& exper : st.meta.experiments) {
      std::string diffrn_id =
        impl::string_or_dot(join_str(exper.diffraction_ids, ","));
      for (const ReflectionsInfo& shell : exper.shells)
        shell_loop.add_row({std::to_string(++n),
                            diffrn_id,
                            impl::number_or_qmark(shell.resolution_high),
                            impl::number_or_qmark(shell.resolution_low),
                            impl::number_or_qmark(shell.completeness),
                            impl::number_or_qmark(shell.redundancy),
                            impl::number_or_qmark(shell.r_merge),
                            impl::number_or_qmark(shell.r_sym),
                            impl::number_or_qmark(shell.mean_I_over_sigma)});
    }
  }

  if (groups.refine && !st.meta.refinement.empty()) {
    block.items.reserve(block.items.size() + 4);
    cif::Loop& loop = block.init_mmcif_loop("_refine.", {
        "entry_id",
        "pdbx_refine_id",
        "ls_d_res_high",
        "ls_d_res_low",
        "ls_percent_reflns_obs",
        "ls_number_reflns_obs"});
    cif::Loop& analyze_loop = block.init_mmcif_loop("_refine_analyze.", {
        "entry_id",
        "pdbx_refine_id",
        "Luzzati_coordinate_error_obs"});
    cif::Loop& restr_loop = block.init_mmcif_loop("_refine_ls_restr.", {
        "pdbx_refine_id", "type",
        "number", "weight", "pdbx_restraint_function", "dev_ideal"});
    cif::Loop& shell_loop = block.init_mmcif_loop("_refine_ls_shell.", {
        "pdbx_refine_id",
        "d_res_high",
        "d_res_low",
        "percent_reflns_obs",
        "number_reflns_obs",
        "number_reflns_R_free",
        "R_factor_obs",
        "R_factor_R_work",
        "R_factor_R_free"});
    for (size_t i = 0; i != st.meta.refinement.size(); ++i) {
      const RefinementInfo& ref = st.meta.refinement[i];
      loop.values.push_back(id);
      loop.values.push_back(cif::quote(ref.id));
      loop.values.push_back(impl::number_or_dot(ref.resolution_high));
      loop.values.push_back(impl::number_or_dot(ref.resolution_low));
      loop.values.push_back(impl::number_or_dot(ref.completeness));
      loop.values.push_back(impl::int_or_dot(ref.reflection_count));
      auto add = [&](const std::string& tag, const std::string& val) {
        if (i == 0)
          loop.tags.push_back("_refine." + tag);
        loop.values.push_back(val);
      };
      if (st.meta.has(&RefinementInfo::rfree_set_count))
        add("ls_number_reflns_R_free", impl::int_or_dot(ref.rfree_set_count));
      if (st.meta.has(&RefinementInfo::r_all))
        add("ls_R_factor_obs", impl::number_or_qmark(ref.r_all));
      if (st.meta.has(&RefinementInfo::r_work))
        add("ls_R_factor_R_work", impl::number_or_qmark(ref.r_work));
      if (st.meta.has(&RefinementInfo::r_free))
        add("ls_R_factor_R_free", impl::number_or_qmark(ref.r_free));
      if (st.meta.has(&RefinementInfo::cross_validation_method))
        add("pdbx_ls_cross_valid_method",
            impl::string_or_qmark(ref.cross_validation_method));
      if (st.meta.has(&RefinementInfo::rfree_selection_method))
        add("pdbx_R_Free_selection_details",
            impl::string_or_qmark(ref.rfree_selection_method));
      if (st.meta.has(&RefinementInfo::mean_b))
        add("B_iso_mean", impl::number_or_qmark(ref.mean_b));
      if (st.meta.has(&RefinementInfo::aniso_b)) {
        if (i == 0)
          for (const char* index : {"[1][1]", "[2][2]", "[3][3]",
                                    "[1][2]", "[1][3]", "[2][3]"})
            loop.tags.push_back(std::string("_refine.aniso_B") + index);
        const Mat33& t = ref.aniso_b;
        for (double d : {t[0][0], t[1][1], t[2][2], t[0][1], t[0][2], t[1][2]})
          loop.values.push_back(impl::number_or_qmark(d));
      }
      if (st.meta.has(&RefinementInfo::dpi_blow_r))
        add("pdbx_overall_SU_R_Blow_DPI",
            impl::number_or_qmark(ref.dpi_blow_r));
      if (st.meta.has(&RefinementInfo::dpi_blow_rfree))
        add("pdbx_overall_SU_R_free_Blow_DPI",
            impl::number_or_qmark(ref.dpi_blow_rfree));
      if (st.meta.has(&RefinementInfo::dpi_cruickshank_r))
        add("overall_SU_R_Cruickshank_DPI",
            impl::number_or_qmark(ref.dpi_cruickshank_r));
      if (st.meta.has(&RefinementInfo::dpi_cruickshank_rfree))
        add("pdbx_overall_SU_R_free_Cruickshank_DPI",
            impl::number_or_qmark(ref.dpi_cruickshank_rfree));
      if (st.meta.has(&RefinementInfo::cc_fo_fc))
        add("correlation_coeff_Fo_to_Fc", impl::number_or_qmark(ref.cc_fo_fc));
      if (st.meta.has(&RefinementInfo::cc_fo_fc_free))
        add("correlation_coeff_Fo_to_Fc_free",
            impl::number_or_qmark(ref.cc_fo_fc_free));
      if (!st.meta.solved_by.empty())
        add("pdbx_method_to_determine_struct", impl::string_or_qmark(st.meta.solved_by));
      if (!st.meta.starting_model.empty())
        add("pdbx_starting_model", impl::string_or_qmark(st.meta.starting_model));
      if (!std::isnan(ref.luzzati_error))
        analyze_loop.add_row({id,
                              cif::quote(ref.id),
                              impl::number_or_qmark(ref.luzzati_error)});
      for (const RefinementInfo::Restr& restr : ref.restr_stats)
        restr_loop.add_row({cif::quote(ref.id),
                            cif::quote(restr.name),
                            impl::int_or_qmark(restr.count),
                            impl::number_or_qmark(restr.weight),
                            impl::string_or_qmark(restr.function),
                            impl::number_or_qmark(restr.dev_ideal)});
      for (const BasicRefinementInfo& bin : ref.bins)
        shell_loop.add_row({cif::quote(ref.id),
                            impl::number_or_dot(bin.resolution_high),
                            impl::number_or_qmark(bin.resolution_low),
                            impl::number_or_qmark(bin.completeness),
                            impl::int_or_qmark(bin.reflection_count),
                            impl::int_or_qmark(bin.rfree_set_count),
                            impl::number_or_qmark(bin.r_all),
                            impl::number_or_qmark(bin.r_work),
                            impl::number_or_qmark(bin.r_free)});
    }
    assert(loop.values.size() % loop.tags.size() == 0);
  }

  if (groups.title_keywords) {
    auto title = st.info.find("_struct.title");
    if (title != st.info.end()) {
      impl::ItemSpan span(block.items, "_struct.");
      span.set_pair("_struct.entry_id", id);
      span.set_pair(title->first, cif::quote(title->second));
    }
    auto pdbx_keywords = st.info.find("_struct_keywords.pdbx_keywords");
    auto keywords = st.info.find("_struct_keywords.text");
    impl::ItemSpan span(block.items, "_struct_keywords.");
    if (pdbx_keywords != st.info.end() || keywords != st.info.end())
      span.set_pair("_struct_keywords.entry_id", id);
    if (pdbx_keywords != st.info.end())
      span.set_pair(pdbx_keywords->first, cif::quote(pdbx_keywords->second));
    if (keywords != st.info.end())
      span.set_pair(keywords->first, cif::quote(keywords->second));
  }

  if (groups.ncs)
    impl::write_ncs_oper(st, block);

  if (groups.struct_asym) {
    cif::Loop& asym_loop = block.init_mmcif_loop("_struct_asym.",
                                                 {"id", "entity_id"});
    for (const Chain& chain : st.models[0].chains)
      for (ConstResidueSpan& sub : chain.subchains())
        if (!sub.subchain_id().empty()) {
          const Entity* ent = st.get_entity_of(sub);
          asym_loop.add_row({sub.subchain_id(),
                             (ent ? impl::qchain(ent->name) : "?")});
        }
  }

  if (groups.origx) { // _database_PDB_matrix (ORIGX)
    if (st.has_origx && !st.origx.is_identity()) {
      impl::ItemSpan span(block.items, "_database_PDB_matrix.");
      span.set_pair("_database_PDB_matrix.entry_id", id);
      std::string prefix = "_database_PDB_matrix.origx";
      for (int i = 0; i < 3; ++i) {
        std::string s = "[" + to_string(i+1) + "]";
        span.set_pair(prefix + s + "[1]", to_str(st.origx.mat[i][0]));
        span.set_pair(prefix + s + "[2]", to_str(st.origx.mat[i][1]));
        span.set_pair(prefix + s + "[3]", to_str(st.origx.mat[i][2]));
        span.set_pair(prefix + "_vector" + s, to_str(st.origx.vec.at(i)));
      }
    }
  }

  if (groups.struct_conf && !st.helices.empty()) {
    cif::Loop& struct_conf_loop = block.init_mmcif_loop("_struct_conf.",
        {"conf_type_id", "id",
         "beg_auth_asym_id", "beg_label_asym_id", "beg_label_comp_id",
         "beg_label_seq_id", "beg_auth_seq_id", "pdbx_beg_PDB_ins_code",
         "end_auth_asym_id", "end_label_asym_id", "end_label_comp_id",
         "end_label_seq_id", "end_auth_seq_id", "pdbx_end_PDB_ins_code",
         "pdbx_PDB_helix_class", "pdbx_PDB_helix_length"});
    int count = 0;
    for (const Helix& helix : st.helices) {
      const_CRA cra1 = st.models[0].find_cra(helix.start);
      const_CRA cra2 = st.models[0].find_cra(helix.end);
      if (!cra1.residue || !cra2.residue)
        continue;
      struct_conf_loop.add_row({
        "HELX_P",                                     // conf_type_id
        "H" + std::to_string(++count),                // id
        impl::qchain(cra1.chain->name),               // beg_auth_asym_id
        impl::subchain_or_dot(*cra1.residue),         // beg_label_asym_id
        cra1.residue->name,                           // beg_label_comp_id
        cra1.residue->label_seq.str(),                // beg_label_seq_id
        cra1.residue->seqid.num.str(),                // beg_auth_seq_id
        impl::pdbx_icode(*cra1.residue),              // beg_PDB_ins_code
        impl::qchain(cra2.chain->name),               // end_auth_asym_id
        impl::subchain_or_dot(*cra2.residue),         // end_label_asym_id
        cra2.residue->name,                           // end_label_comp_id
        cra2.residue->label_seq.str(),                // end_label_seq_id
        cra2.residue->seqid.num.str(),                // end_auth_seq_id
        impl::pdbx_icode(*cra2.residue),              // end_PDB_ins_code
        std::to_string((int)helix.pdb_helix_class),   // pdbx_PDB_helix_class
        impl::int_or_qmark(helix.length)              // pdbx_PDB_helix_length
      });
    }
    if (count != 0)
      block.set_pair("_struct_conf_type.id", "HELX_P");
  }

  // _struct_sheet*
  if (groups.struct_sheet && !st.sheets.empty()) {
    cif::Loop& sheet_loop = block.init_mmcif_loop("_struct_sheet.",
                                                  {"id", "number_strands"});
    for (const Sheet& sheet : st.sheets)
      sheet_loop.add_row({impl::string_or_dot(sheet.name),
                          std::to_string(sheet.strands.size())});

    cif::Loop& order_loop = block.init_mmcif_loop("_struct_sheet_order.",
                    {"sheet_id", "range_id_1", "range_id_2", "sense"});
    for (const Sheet& sheet : st.sheets)
      for (size_t i = 1; i < sheet.strands.size(); ++i) {
        const Sheet::Strand& strand = sheet.strands[i];
        if (strand.sense != 0)
          order_loop.add_row({impl::string_or_dot(sheet.name),
                              std::to_string(i), std::to_string(i+1),
                              strand.sense > 0 ? "parallel" : "anti-parallel"});
      }

    cif::Loop& range_loop = block.init_mmcif_loop("_struct_sheet_range.",
        {"sheet_id", "id",
         "beg_auth_asym_id", "beg_label_asym_id", "beg_label_comp_id",
         "beg_label_seq_id", "beg_auth_seq_id", "pdbx_beg_PDB_ins_code",
         "end_auth_asym_id", "end_label_asym_id", "end_label_comp_id",
         "end_label_seq_id", "end_auth_seq_id", "pdbx_end_PDB_ins_code"});
    for (const Sheet& sheet : st.sheets)
      for (size_t i = 0; i < sheet.strands.size(); ++i) {
        const Sheet::Strand& strand = sheet.strands[i];
        const_CRA cra1 = st.models[0].find_cra(strand.start);
        const_CRA cra2 = st.models[0].find_cra(strand.end);
        if (!cra1.residue || !cra2.residue)
          continue;
        range_loop.add_row({
          impl::string_or_dot(sheet.name),            // sheet_id
          std::to_string(i+1),                        // id
          impl::qchain(cra1.chain->name),             // beg_auth_asym_id
          impl::subchain_or_dot(*cra1.residue),       // beg_label_asym_id
          cra1.residue->name,                         // beg_label_comp_id
          cra1.residue->label_seq.str(),              // beg_label_seq_id
          cra1.residue->seqid.num.str(),              // beg_auth_seq_id
          impl::pdbx_icode(*cra1.residue),            // beg_PDB_ins_code
          impl::qchain(cra2.chain->name),             // end_auth_asym_id
          impl::subchain_or_dot(*cra2.residue),       // end_label_asym_id
          cra2.residue->name,                         // end_label_comp_id
          cra2.residue->label_seq.str(),              // end_label_seq_id
          cra2.residue->seqid.num.str(),              // end_auth_seq_id
          impl::pdbx_icode(*cra2.residue)             // end_PDB_ins_code
        });
    }

    cif::Loop& hbond_loop = block.init_mmcif_loop("_pdbx_struct_sheet_hbond.",
        {"sheet_id", "range_id_1", "range_id_2",
         "range_1_auth_asym_id", "range_1_label_asym_id",
         "range_1_label_comp_id", "range_1_label_seq_id", "range_1_auth_seq_id",
         "range_1_PDB_ins_code", "range_1_label_atom_id",
         "range_2_auth_asym_id", "range_2_label_asym_id",
         "range_2_label_comp_id", "range_2_label_seq_id", "range_2_auth_seq_id",
         "range_2_PDB_ins_code", "range_2_label_atom_id"});
    for (const Sheet& sheet : st.sheets)
      for (size_t i = 1; i < sheet.strands.size(); ++i) {
        const Sheet::Strand& strand = sheet.strands[i];
        if (strand.hbond_atom2.atom_name.empty())
          continue;
        // hbond_atomN is not a full atom "address": altloc is missing
        const_CRA cra1 = st.models[0].find_cra(strand.hbond_atom1);
        const_CRA cra2 = st.models[0].find_cra(strand.hbond_atom2);
        if (!cra1.residue || !cra2.residue)
          continue;
        hbond_loop.add_row({
          impl::string_or_dot(sheet.name),            // sheet_id
          std::to_string(i),                          // range_id_1
          std::to_string(i+1),                        // range_id_2
          impl::qchain(cra1.chain->name),             // range_1_auth_asym_id
          impl::subchain_or_dot(*cra1.residue),       // range_1_label_asym_id
          cra1.residue->name,                         // range_1_label_comp_id
          cra1.residue->label_seq.str(),              // range_1_label_seq_id
          cra1.residue->seqid.num.str(),              // range_1_auth_seq_id
          impl::pdbx_icode(*cra1.residue),            // range_1_PDB_ins_code
          cif::quote(strand.hbond_atom1.atom_name),   // range_1_label_atom_id
          impl::qchain(cra2.chain->name),             // range_2_auth_asym_id
          impl::subchain_or_dot(*cra2.residue),       // range_2_label_asym_id
          cra2.residue->name,                         // range_2_label_comp_id
          cra2.residue->label_seq.str(),              // range_2_label_seq_id
          cra2.residue->seqid.num.str(),              // range_2_auth_seq_id
          impl::pdbx_icode(*cra2.residue),            // range_2_PDB_ins_code
          cif::quote(strand.hbond_atom2.atom_name)    // range_2_label_atom_id
        });
    }
  }

  // _pdbx_struct_assembly* and _struct_biol are REMARK 300/350 in PDB
  if (groups.struct_biol && !st.meta.remark_300_detail.empty()) {
    impl::ItemSpan span(block.items, "_struct_biol.");
    span.set_pair("_struct_biol.id", "1");
    span.set_pair("_struct_biol.details", cif::quote(st.meta.remark_300_detail));
  }
  if (groups.assembly && !st.assemblies.empty())
    impl::write_assemblies(st, block);

  if (groups.conn)
    impl::write_struct_conn(st, block);

  if ( groups.cis) {  // _struct_mon_prot_cis
    cif::Loop& prot_cis_loop = block.init_mmcif_loop("_struct_mon_prot_cis.",
        {"pdbx_id", "pdbx_PDB_model_num", "label_asym_id", "label_seq_id",
         "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
         "label_comp_id", "label_alt_id"});
    for (const Model& model : st.models)
      for (const Chain& chain : model.chains)
        for (const Residue& res : chain.residues)
          if (res.is_cis)
            prot_cis_loop.add_row({to_string(prot_cis_loop.length()+1),
                                   model.name, impl::subchain_or_dot(res),
                                   res.label_seq.str(), impl::qchain(chain.name),
                                   res.seqid.num.str(), impl::pdbx_icode(res),
                                   res.name, "."});
  }

  // _atom_sites (SCALE)
  if (groups.scale && (st.has_origx || st.cell.explicit_matrices)) {
    impl::ItemSpan span(block.items, "_atom_sites.");
    span.set_pair("_atom_sites.entry_id", id);
    std::string prefix = "_atom_sites.fract_transf_";
    for (int i = 0; i < 3; ++i) {
      std::string idx = "[" + std::to_string(i + 1) + "]";
      const auto& frac = st.cell.frac;
      span.set_pair(prefix + "matrix" + idx + "[1]", to_str(frac.mat[i][0]));
      span.set_pair(prefix + "matrix" + idx + "[2]", to_str(frac.mat[i][1]));
      span.set_pair(prefix + "matrix" + idx + "[3]", to_str(frac.mat[i][2]));
      span.set_pair(prefix + "vector" + idx, to_str(frac.vec.at(i)));
    }
  }

  // _atom_type
  if (groups.atom_type) {
    std::array<bool, (int)El::END> types{};
    for (const Model& model : st.models)
      for (const Chain& chain : model.chains)
        for (const Residue& res : chain.residues)
          for (const Atom& atom : res.atoms)
            types[atom.element.ordinal()] = true;
    cif::Loop& atom_type_loop = block.init_mmcif_loop("_atom_type.", {"symbol"});
    for (int i = 0; i < (int)El::END; ++i)
      if (types[i])
        atom_type_loop.add_row({Element((El)i).uname()});
  }

  if (groups.entity_poly_seq) {
    // SEQRES from PDB doesn't record microheterogeneity, so if the resulting
    // cif has unknown("?") _entity_poly_seq.num, it cannot be trusted.
    cif::Loop& poly_loop = block.init_mmcif_loop("_entity_poly_seq.",
                                           {"entity_id", "num", "mon_id"});
    for (const Entity& ent : st.entities)
      if (ent.entity_type == EntityType::Polymer)
        for (size_t i = 0; i != ent.full_sequence.size(); ++i) {
          const std::string& mon_ids = ent.full_sequence[i];
          std::string num = std::to_string(i+1);
          size_t start = 0, end;
          while ((end = mon_ids.find(',', start)) != std::string::npos) {
            poly_loop.add_row({impl::qchain(ent.name), num,
                               mon_ids.substr(start, end-start)});
            start = end + 1;
          }
          poly_loop.add_row({impl::qchain(ent.name), num, mon_ids.substr(start)});
        }
  }

  if (groups.atoms)
    impl::add_cif_atoms(st, block, groups.group_pdb);

  if (groups.tls && st.meta.has_tls()) {
    cif::Loop& loop = block.init_mmcif_loop("_pdbx_refine_tls.", {
        "pdbx_refine_id", "id",
        "T[1][1]", "T[2][2]", "T[3][3]", "T[1][2]", "T[1][3]", "T[2][3]",
        "L[1][1]", "L[2][2]", "L[3][3]", "L[1][2]", "L[1][3]", "L[2][3]",
        "S[1][1]", "S[1][2]", "S[1][3]",
        "S[2][1]", "S[2][2]", "S[2][3]",
        "S[3][1]", "S[3][2]", "S[3][3]",
        "origin_x", "origin_y", "origin_z"});
    for (const RefinementInfo& ref : st.meta.refinement)
      for (const TlsGroup& tls : ref.tls_groups) {
        const Mat33& T = tls.T;
        const Mat33& L = tls.L;
        const Mat33& S = tls.S;
        auto q = impl::number_or_qmark;
        loop.add_row({cif::quote(ref.id), impl::string_or_dot(tls.id),
                      q(T[0][0]), q(T[1][1]), q(T[2][2]),
                      q(T[0][1]), q(T[0][2]), q(T[1][2]),
                      q(L[0][0]), q(L[1][1]), q(L[2][2]),
                      q(L[0][1]), q(L[0][2]), q(L[1][2]),
                      q(S[0][0]), q(S[0][1]), q(S[0][2]),
                      q(S[1][0]), q(S[1][1]), q(S[1][2]),
                      q(S[2][0]), q(S[2][1]), q(S[2][2]),
                      q(tls.origin.x), q(tls.origin.y), q(tls.origin.z)});
      }
    cif::Loop& group_loop = block.init_mmcif_loop("_pdbx_refine_tls_group.", {
        "id", "refine_tls_id", "beg_auth_asym_id", "beg_auth_seq_id",
        "end_auth_asym_id", "end_auth_seq_id", "selection_details"});
    int counter = 1;
    for (const RefinementInfo& ref : st.meta.refinement)
      for (const TlsGroup& tls : ref.tls_groups)
        for (const TlsGroup::Selection& sel : tls.selections)
          group_loop.add_row({std::to_string(counter++),
                              impl::string_or_dot(tls.id),
                              impl::string_or_qmark(sel.chain),
                              sel.res_begin.num ? sel.res_begin.str() : "?",
                              impl::string_or_qmark(sel.chain),
                              sel.res_end.num ? sel.res_end.str() : "?",
                              impl::string_or_qmark(sel.details)});
  }

  if (groups.software && !st.meta.software.empty()) {
    cif::Loop& loop = block.init_mmcif_loop("_software.",
                 {"pdbx_ordinal", "classification", "name", "version", "date"});
    for (const SoftwareItem& item : st.meta.software)
      loop.add_row({
          std::to_string(item.pdbx_ordinal),
          cif::quote(software_classification_to_string(item.classification)),
          cif::quote(item.name),
          impl::string_or_dot(item.version),
          impl::string_or_qmark(item.date)});
  }
}

cif::Document make_mmcif_document(const Structure& st) {
  cif::Document doc;
  doc.blocks.resize(1);
  update_mmcif_block(st, doc.blocks[0], MmcifOutputGroups(true));
  return doc;
}

cif::Block make_mmcif_headers(const Structure& st) {
  cif::Block block;
  MmcifOutputGroups groups(true);
  groups.atoms = false;
  update_mmcif_block(st, block, groups);
  return block;
}

void add_minimal_mmcif_data(const Structure& st, cif::Block& block) {
  impl::ItemSpan cell_span(block.items, "_cell.");
  impl::write_cell_parameters(st.cell, cell_span);
  block.set_pair("_symmetry.space_group_name_H-M",
                 cif::quote(st.spacegroup_hm));
  impl::write_ncs_oper(st, block);
  impl::add_cif_atoms(st, block, /*use_group_pdb=*/false);
}

} // namespace gemmi
#endif // GEMMI_WRITE_IMPLEMENTATION

#endif
