// Copyright 2019 Global Phasing Ltd.
//
// Metadata from coordinate files.

#ifndef GEMMI_METADATA_HPP_
#define GEMMI_METADATA_HPP_

#include <algorithm>     // for any_of
#include <string>
#include <vector>
#include "math.hpp"      // for Mat33
#include "unitcell.hpp"  // for Position, Asu
#include "seqid.hpp"     // for SeqId

namespace gemmi {

// corresponds to the mmCIF _software category
struct SoftwareItem {
  enum Classification {
    DataCollection, DataExtraction, DataProcessing, DataReduction,
    DataScaling, ModelBuilding, Phasing, Refinement, Unspecified
  };
  std::string name;
  std::string version;
  std::string date;
  Classification classification = Unspecified;
  int pdbx_ordinal = -1;
};

// Information from REMARK 200/230 is significantly expanded in PDBx/mmCIF.
// These remarks corresponds to data across 12 mmCIF categories
// including categories _exptl, _reflns, _exptl_crystal, _diffrn and others.
// _exptl and _reflns seem to be 1:1. Usually we have one experiment (_exptl),
// except for a joint refinement (e.g. X-ray + neutron data).
// Both crystal (_exptl_crystal) and reflection statistics (_reflns) can
// be associated with multiple diffraction sets (_diffrn).
// But if we use the PDB format, only one diffraction set per method
// can be described.

struct ReflectionsInfo {
  double resolution_high = NAN; // _reflns.d_resolution_high
                                // (or _reflns_shell.d_res_high)
  double resolution_low = NAN;  // _reflns.d_resolution_low
  double completeness = NAN;    // _reflns.percent_possible_obs
  double redundancy = NAN;      // _reflns.pdbx_redundancy
  double r_merge = NAN;         // _reflns.pdbx_Rmerge_I_obs
  double r_sym = NAN;           // _reflns.pdbx_Rsym_value
  double mean_I_over_sigma = NAN; // _reflns.pdbx_netI_over_sigmaI
};

// _exptl has no id, _exptl.method is key item and must be unique
struct ExperimentInfo {
  std::string method;             // _exptl.method
  int number_of_crystals = -1;    // _exptl.crystals_number
  int unique_reflections = -1;    // _reflns.number_obs
  ReflectionsInfo reflections;
  double b_wilson = NAN;          // _reflns.B_iso_Wilson_estimate
  std::vector<ReflectionsInfo> shells;
  std::vector<std::string> diffraction_ids;
};

struct DiffractionInfo {
  std::string id;                // _diffrn.id
  double temperature = NAN;      // _diffrn.ambient_temp
  std::string source;            // _diffrn_source.source
  std::string source_type;       // _diffrn_source.type
  std::string synchrotron;       // _diffrn_source.pdbx_synchrotron_site
  std::string beamline;          // _diffrn_source.pdbx_synchrotron_beamline
  std::string wavelengths;       // _diffrn_source.pdbx_wavelength
  std::string scattering_type;   // _diffrn_radiation.pdbx_scattering_type
  char mono_or_laue = '\0'; // _diffrn_radiation.pdbx_monochromatic_or_laue_m_l
  std::string monochromator;     // _diffrn_radiation.monochromator
  std::string collection_date;   // _diffrn_detector.pdbx_collection_date
  std::string optics;            // _diffrn_detector.details
  std::string detector;          // _diffrn_detector.detector
  std::string detector_make;     // _diffrn_detector.type
};

struct CrystalInfo {
  std::string id;                 // _exptl_crystal.id
  std::string description;        // _exptl_crystal.description
  double ph = NAN;                // _exptl_crystal_grow.pH
  std::string ph_range;           // _exptl_crystal_grow.pdbx_pH_range
  std::vector<DiffractionInfo> diffractions;
};


struct TlsGroup {
  struct Selection {
    std::string chain;
    SeqId res_begin;
    SeqId res_end;
    std::string details;  // _pdbx_refine_tls_group.selection_details
  };
  std::string id;           // _pdbx_refine_tls.id
  std::vector<Selection> selections;
  Position origin;          // _pdbx_refine_tls.origin_x/y/z
  Mat33 T;                  // _pdbx_refine_tls.T[][]
  Mat33 L;                  // _pdbx_refine_tls.L[][]
  Mat33 S;                  // _pdbx_refine_tls.S[][]
};

// RefinementInfo corresponds to REMARK 3.
// BasicRefinementInfo is used for both total and per-bin statistics.
// For per-bin data, each values corresponds to one _refine_ls_shell.* tag.
struct BasicRefinementInfo {
  double resolution_high = NAN;      // _refine.ls_d_res_high
  double resolution_low = NAN;       // _refine.ls_d_res_low
  double completeness = NAN;         // _refine.ls_percent_reflns_obs
  int reflection_count = -1;         // _refine.ls_number_reflns_obs
  int rfree_set_count = -1;          // _refine.ls_number_reflns_R_free
  double r_all = NAN;                // _refine.ls_R_factor_obs
  double r_work = NAN;               // _refine.ls_R_factor_R_work
  double r_free = NAN;               // _refine.ls_R_factor_R_free
};

struct RefinementInfo : BasicRefinementInfo {
  struct Restr {
    std::string name;
    int count = -1;
    double weight = NAN;
    std::string function;
    double dev_ideal = NAN;
    Restr(const std::string& name_) : name(name_) {}
  };
  std::string id;
  std::string cross_validation_method; // _refine.pdbx_ls_cross_valid_method
  std::string rfree_selection_method;  // _refine.pdbx_R_Free_selection_details
  int bin_count = -1;        // _refine_ls_shell.pdbx_total_number_of_bins_used
  std::vector<BasicRefinementInfo> bins;
  double mean_b = NAN;                // _refine.B_iso_mean
  Mat33 aniso_b{NAN};                 // _refine.aniso_B[][]
  double luzzati_error = NAN; // _refine_analyze.Luzzati_coordinate_error_obs
  double dpi_blow_r = NAN;            // _refine.pdbx_overall_SU_R_Blow_DPI
  double dpi_blow_rfree = NAN;        // _refine.pdbx_overall_SU_R_free_Blow_DPI
  double dpi_cruickshank_r = NAN;     // _refine.overall_SU_R_Cruickshank_DPI
  double dpi_cruickshank_rfree = NAN; // _refine.pdbx_overall_SU_R_free_Cruickshank_DPI
  double cc_fo_fc = NAN;              // _refine.correlation_coeff_Fo_to_Fc
  double cc_fo_fc_free = NAN;         // _refine.correlation_coeff_Fo_to_Fc_free
  std::vector<Restr> restr_stats;     // _refine_ls_restr
  std::vector<TlsGroup> tls_groups;   // _pdbx_refine_tls
  std::string remarks;
};


struct Metadata {
  std::vector<ExperimentInfo> experiments;
  std::vector<CrystalInfo> crystals;
  std::vector<RefinementInfo> refinement;
  std::vector<SoftwareItem> software;
  std::string solved_by;       // _refine.pdbx_method_to_determine_struct
  std::string starting_model;  // _refine.pdbx_starting_model
  std::string remark_300_detail; // _struct_biol.details

  bool has(double RefinementInfo::*field) const {
    return std::any_of(refinement.begin(), refinement.end(),
            [&](const RefinementInfo& r) { return !std::isnan(r.*field); });
  }
  bool has(int RefinementInfo::*field) const {
    return std::any_of(refinement.begin(), refinement.end(),
            [&](const RefinementInfo& r) { return r.*field != -1; });
  }
  bool has(std::string RefinementInfo::*field) const {
    return std::any_of(refinement.begin(), refinement.end(),
            [&](const RefinementInfo& r) { return !(r.*field).empty(); });
  }
  bool has(Mat33 RefinementInfo::*field) const {
    return std::any_of(refinement.begin(), refinement.end(),
        [&](const RefinementInfo& r) { return !std::isnan((r.*field)[0][0]); });
  }
  bool has_restr() const {
    return std::any_of(refinement.begin(), refinement.end(),
            [&](const RefinementInfo& r) { return !r.restr_stats.empty(); });
  }
  bool has_tls() const {
    return std::any_of(refinement.begin(), refinement.end(),
            [&](const RefinementInfo& r) { return !r.tls_groups.empty(); });
  }
};


// Entity description.
//
// values corresponding to mmCIF _entity.type
enum class EntityType : unsigned char {
  Unknown,
  Polymer,
  NonPolymer,
  Branched, // introduced in 2020
  // _entity.type macrolide is in PDBx/mmCIF, but no PDB entry uses it
  //Macrolide,
  Water
};

// values corresponding to mmCIF _entity_poly.type
enum class PolymerType : unsigned char {
  Unknown,       // unknown or not applicable
  PeptideL,      // polypeptide(L) in mmCIF (168923 values in the PDB in 2017)
  PeptideD,      // polypeptide(D) (57 values)
  Dna,           // polydeoxyribonucleotide (9905)
  Rna,           // polyribonucleotide (4559)
  DnaRnaHybrid,  // polydeoxyribonucleotide/polyribonucleotide hybrid (156)
  SaccharideD,   // polysaccharide(D) (18)
  SaccharideL,   // polysaccharide(L) (0)
  Pna,           // peptide nucleic acid (2)
  CyclicPseudoPeptide,  // cyclic-pseudo-peptide (1)
  Other,         // other (4)
};

inline bool is_polypeptide(PolymerType pt) {
  return pt == PolymerType::PeptideL || pt == PolymerType::PeptideD;
}

inline bool is_polynucleotide(PolymerType pt) {
  return pt == PolymerType::Dna || pt == PolymerType::Rna ||
         pt == PolymerType::DnaRnaHybrid;
}

struct Entity {
  struct DbRef {
    std::string db_name;
    std::string accession_code;
    std::string id_code;
    std::string isoform;  // pdbx_db_isoform
    SeqId seq_begin, seq_end;
    SeqId db_begin, db_end;
    SeqId::OptionalNum label_seq_begin, label_seq_end;
  };
  std::string name;
  std::vector<std::string> subchains;
  EntityType entity_type = EntityType::Unknown;
  PolymerType polymer_type = PolymerType::Unknown;
  std::vector<DbRef> dbrefs;
  // SEQRES or entity_poly_seq with microheterogeneity as comma-separated names
  std::vector<std::string> full_sequence;

  explicit Entity(std::string name_) noexcept : name(name_) {}
  static std::string first_mon(const std::string& mon_list) {
    return mon_list.substr(0, mon_list.find(','));
  }
};


// A connection. Corresponds to _struct_conn.
// Symmetry operators are not trusted and not stored.
// We assume that the nearest symmetry mate is connected.
struct Connection {
  // in write_struct_conn() we assume that Unknown is at the end
  enum Type { Covale=0, Disulf, Hydrog, MetalC, Unknown };
  std::string name;
  std::string link_id;  // _struct_conn.ccp4_link_id (== _chem_link.id)
  Type type = Unknown;
  Asu asu = Asu::Any;
  AtomAddress partner1, partner2;
  double reported_distance = 0.0;
};


// Secondary structure. PDBx/mmCIF stores helices and sheets separately.

// mmCIF spec defines 32 possible values for _struct_conf.conf_type_id -
// "the type of the conformation of the backbone of the polymer (whether
// protein or nucleic acid)". But as of 2019 only HELX_P is used (not counting
// TURN_P that occurs in only 6 entries). The actual helix type is given
// by numeric value of _struct_conf.pdbx_PDB_helix_class, which corresponds
// to helixClass from the PDB HELIX record. These values are in the range 1-10.
// As of 2019 it's almost only type 1 and 5:
// 3116566 of  1 - right-handed alpha
//      16 of  2 - right-handed omega
//      84 of  3 - right-handed pi
//      79 of  4 - right-handed gamma
// 1063337 of  5 - right-handed 3-10
//      27 of  6 - left-handed alpha
//       5 of  7 - left-handed omega
//       2 of  8 - left-handed gamma
//       8 of  9 - 2-7 ribbon/helix
//      46 of 10 - polyproline
struct Helix {
  enum HelixClass {
    UnknownHelix, RAlpha, ROmega, RPi, RGamma, R310,
    LAlpha, LOmega, LGamma, Helix27, HelixPolyProlineNone
  };
  AtomAddress start, end;
  HelixClass pdb_helix_class = UnknownHelix;
  int length = -1;
  void set_helix_class_as_int(int n) {
    if (n >= 1 && n <= 10)
      pdb_helix_class = static_cast<HelixClass>(n);
  }
};

struct Sheet {
  struct Strand {
    AtomAddress start, end;
    AtomAddress hbond_atom2, hbond_atom1;
    int sense;  // 0 = first strand, 1 = parallel, -1 = anti-parallel.
    std::string name; // optional, _struct_sheet_range.id if from mmCIF
  };
  std::string name;
  std::vector<Strand> strands;
  explicit Sheet(std::string sheet_id) noexcept : name(sheet_id) {}
};


// bioassembly / biomolecule
struct Assembly {
  struct Operator {
    std::string name; // optional
    std::string type; // optional (from mmCIF only)
    Transform transform;
  };
  struct Gen {
    std::vector<std::string> chains;
    std::vector<std::string> subchains;
    std::vector<Operator> operators;
  };
  enum class SpecialKind {
    NA, CompleteIcosahedral, RepresentativeHelical, CompletePoint
  };
  std::string name;
  bool author_determined = false;
  bool software_determined = false;
  SpecialKind special_kind = SpecialKind::NA;
  int oligomeric_count = 0;
  std::string oligomeric_details;
  std::string software_name;
  double absa = NAN; // TOTAL BURIED SURFACE AREA: ... ANGSTROM**2
  double ssa = NAN;  // SURFACE AREA OF THE COMPLEX: ... ANGSTROM**2
  double more = NAN; // CHANGE IN SOLVENT FREE ENERGY: ... KCAL/MOL
  std::vector<Gen> generators;
  Assembly(const std::string& name_) : name(name_) {}
};

} // namespace gemmi
#endif
