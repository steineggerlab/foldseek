///////////////////////////////////////////////////////////////////////////
//  CLePAPS:                                                             //
//  Conformational Letter based Pairwise Alignment of Protein Structure  //
//                                                                       //
//  Authors: S.Wang & WM.Zheng                                           //
///////////////////////////////////////////////////////////////////////////


#pragma once
#include <vector>
#include "Confo_Lett.h"
#include "PDB_File.h"


//====class: PDB_To_CLE====//
//=> transfrom PDB_Format to CLE (or to XYZ_Format)
class PDB_To_CLE : virtual public PDB_File ,virtual public Confo_Lett
{
public:
	PDB_To_CLE(void);
	~PDB_To_CLE(void);

//--------- virtual_function ------------//
public:
	virtual void pdb_btb_ori(int head,int totnum,vector <XYZ> &r,vector <char> &CLE);
};
//====class: PDB_To_CLE====//over
