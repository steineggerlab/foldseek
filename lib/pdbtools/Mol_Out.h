#pragma once
#include "PDB_Residue.h" 
#include "PDB_Utility.h"
#include "Computation_Utility.h" 

//====class: Mol_Out====//
//=> output structure
class Mol_Out
{
public:
	Mol_Out(void);
	~Mol_Out(void);

//---- macros ----//
public:
	int OutPrt;       // whether output printf    // (default:no)

//--------- functions ---------//
public:
	//---PART_IV:  output_related
	int Output_XYZ(FILE *fp,int len,XYZ *mol,char *AMI,char *CLE,char *ind,char Chain_ID);
	int Output_PDB(FILE *fp,int len,XYZ *mol,char *AMI,char *ind,char Chain_ID);
	int Output_PDB_I(FILE *fp,int len,XYZ *ca,XYZ *cb,char *AMI,char *ind,char Chain_ID);
	int Output_PDB_II(FILE *fp,int len,XYZ **mol,char *AMI,char *ind,char Chain_ID,int type=1);
	int Output_PDB_III(FILE *fp,int len,PDB_Residue *mol,char Chain_ID,
		int OutType=1,int OutMode=1,int OutGlys=1,int OutLast=1);
	int Output_PDB_III(FILE *fp,int len,vector <PDB_Residue> &mol,char Chain_ID,
		int OutType=1,int OutMode=1,int OutGlys=1,int OutLast=1);
};
