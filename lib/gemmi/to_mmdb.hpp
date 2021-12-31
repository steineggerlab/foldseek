// Copyright 2020 Global Phasing Ltd.
//
// copy_to_mmdb(): converts gemmi::Structure to mmdb::Manager.

#ifndef GEMMI_TO_MMDB_HPP_
#define GEMMI_TO_MMDB_HPP_

#include <cstdlib>           // for atoi
#include <cstring>           // for memcpy
#include <gemmi/model.hpp>
#include <gemmi/util.hpp>    // for rtrim_str
#include <mmdb2/mmdb_manager.h>

namespace gemmi {

inline void copy_transform_to_mmdb(const Transform& tr,
                                   mmdb::mat33& mat, mmdb::vect3& vec) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      mat[i][j] = tr.mat[i][j];
    vec[i] = tr.vec.at(i);
  }
}

template<int N>
void strcpy_to_mmdb(char (&dest)[N], const std::string& src) {
  if (src.size() >= N)
    fail("This string is too long: " + src);
  std::memcpy(dest, src.c_str(), src.size() + 1);
}

inline void set_seqid_in_mmdb(int* seqnum, mmdb::InsCode& icode, SeqId seqid) {
  *seqnum = *seqid.num;
  icode[0] = seqid.icode;
  icode[1] = '\0';
}

inline mmdb::Manager* copy_to_mmdb(const Structure& st, mmdb::Manager* manager) {
  for (const std::string& s : st.raw_remarks) {
    std::string line = rtrim_str(s);
    manager->PutPDBString(line.c_str());
  }

  for (int imodel = 0; imodel < (int) st.models.size(); ++imodel) {
    const Model& model = st.models[imodel];
    mmdb::PModel model2 = mmdb::newModel();
    model2->SetMMDBManager(manager, imodel);

    for (const Chain& chain : model.chains) {
      mmdb::PChain chain2 = model2->CreateChain(chain.name.c_str());
      for (const Residue& res : chain.residues) {
        const char icode[2] = {res.seqid.icode, '\0'};
        mmdb::PResidue res2 = chain2->GetResidueCreate(res.name.c_str(),
                                                       *res.seqid.num,
                                                       icode,
                                                       true);
        if (res.is_cis) {
          if (const Residue* next = chain.next_residue(res)) {
            mmdb::CisPep* cispep = new mmdb::CisPep();
            cispep->serNum = model2->GetCisPeps()->Length() + 1;
            strcpy_to_mmdb(cispep->chainID1, chain.name);
            strcpy_to_mmdb(cispep->pep1, res.name);
            set_seqid_in_mmdb(&cispep->seqNum1, cispep->icode1, res.seqid);
            strcpy_to_mmdb(cispep->chainID2, chain.name);
            strcpy_to_mmdb(cispep->pep2, next->name);
            set_seqid_in_mmdb(&cispep->seqNum2, cispep->icode2, next->seqid);
            cispep->modNum = imodel;
            model2->AddCisPep(cispep);
          }
        }
        for (const Atom& atom : res.atoms) {
          mmdb::PAtom atom2 = mmdb::newAtom();
          const char altloc[2] = {atom.altloc, '\0'};
          atom2->SetAtomName(0, atom.serial, atom.padded_name().c_str(),
                             altloc, res.segment.c_str(), atom.element.uname());
          atom2->Het = res.het_flag == 'H';
          atom2->SetCharge(atom.charge);
          atom2->SetCoordinates(atom.pos.x, atom.pos.y, atom.pos.z,
                                atom.occ, atom.b_iso);
          if (atom.aniso.nonzero()) {
            atom2->su11 = atom.aniso.u11;
            atom2->su22 = atom.aniso.u22;
            atom2->su33 = atom.aniso.u33;
            atom2->su12 = atom.aniso.u12;
            atom2->su13 = atom.aniso.u13;
            atom2->su23 = atom.aniso.u23;
            atom2->WhatIsSet |= mmdb::ASET_Anis_tFSigma;
          }
          res2->AddAtom(atom2);
        }
        // TER
        if (res.entity_type == EntityType::Polymer &&
            (&res == &chain.residues.back() ||
             (&res + 1)->entity_type != EntityType::Polymer)) {
          mmdb::PAtom atom2 = mmdb::newAtom();
          atom2->MakeTer();
          atom2->serNum = res.atoms.back().serial + 1;
          res2->AddAtom(atom2);
        }
      }
    }
    int n_model = manager->AddModel(model2);
    manager->GetModel(n_model)->CopyCisPeps(model2);
    manager->PutCell(st.cell.a, st.cell.b, st.cell.c,
                     st.cell.alpha, st.cell.beta, st.cell.gamma, 1);
    manager->SetSpaceGroup(st.spacegroup_hm.c_str());
    mmdb::Cryst* cryst = manager->GetCrystData();
    auto z = st.info.find("_cell.Z_PDB");
    if (z != st.info.end() && !z->second.empty()) {
      cryst->Z = std::atoi(z->second.c_str());
      cryst->WhatIsSet |= mmdb::CSET_ZValue;
    }
    if (st.has_origx && !st.origx.is_identity()) {
      copy_transform_to_mmdb(st.origx, cryst->o, cryst->t);
      cryst->WhatIsSet |= mmdb::CSET_OrigMatrix;
    }
    if (st.cell.explicit_matrices) {
      copy_transform_to_mmdb(st.cell.frac, cryst->s, cryst->u);
      cryst->WhatIsSet |= mmdb::CSET_ScaleMatrix;
    }
    if (!st.ncs.empty()) {
      mmdb::mat33 m;
      mmdb::vect3 v;
      if (st.info.find("_struct_ncs_oper.id") != st.info.end()) {
        copy_transform_to_mmdb(Transform{}, m, v);
        cryst->AddNCSMatrix(m, v, 1);
      }
      for (const NcsOp& op : st.ncs) {
        copy_transform_to_mmdb(op.tr, m, v);
        cryst->AddNCSMatrix(m, v, op.given);
      }
    }
  }
  return manager;
}

} // namespace gemmi
#endif


/*
// Example that uses copy_to_mmdb().
// Reads a coordinate file using gemmi and write it to pdb using mmdb.

#include <gemmi/mmread.hpp>
#include <gemmi/to_mmdb.hpp>

// two arguments expected: input and output paths.
int main (int argc, char** argv) {
  if (argc != 3)
    return 1;

  mmdb::InitMatType();
  mmdb::Manager* manager = new mmdb::Manager();
  try {
    gemmi::Structure st = gemmi::read_structure_file(argv[1]);
    st.merge_chain_parts();
    gemmi::copy_to_mmdb(st, manager);
  } catch(std::runtime_error& e) {
    printf("File reading failed: %s\n", e.what());
    return 1;
  }

  mmdb::ERROR_CODE rc = manager->WritePDBASCII(argv[2]);
  if (rc)
    printf(" ***** ERROR #%i WRITE:\n\n %s\n\n",
           rc, mmdb::GetErrorDescription(rc));
  return 0;
}
*/
