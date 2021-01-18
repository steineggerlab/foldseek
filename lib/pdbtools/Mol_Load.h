///////////////////////////////////////////////////////////////////////////
//  CLePAPS:                                                             //
//  Conformational Letter based Pairwise Alignment of Protein Structure  //
//                                                                       //
//  Authors: S.Wang & WM.Zheng                                           //
///////////////////////////////////////////////////////////////////////////


#pragma once
#include "PDB_To_CLE.h"
#include "Mol_Ori.h"



//====class: Mol_Load====//
//=> input structure (in PDB_Format or XYZ_Format)
class Mol_Load : public PDB_To_CLE,virtual public Mol_Ori
{
public:
	Mol_Load(void);
	~Mol_Load(void);

//=========================== Upload_Input =================================//
public:

//-------- virtual_reloaded_function --------//
public:
	virtual int Upload_Process(string &file,char chain,char chain_bak,int PDBorSEQ,int tag,
		int head_tag,int head_c,int tail_tag,int tail_c,int &moln,
		XYZ *mol,char *AMI,char *CLE,char *ind,PDB_Residue *pdb=0);
};
//====class: Mol_Load====//over
