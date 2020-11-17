///////////////////////////////////////////////////////////////////////////
//  CLePAPS:                                                             //
//  Conformational Letter based Pairwise Alignment of Protein Structure  //
//                                                                       //
//  Authors: S.Wang & WM.Zheng                                           //
///////////////////////////////////////////////////////////////////////////
/*
 
Copyright (c)  2007-2008  Institute of Theoretical Physics,Academia Sinica
All Rights Reserved
 
Permission to use, copy, modify and distribute any part of this CLePAPS
software for educational, research and non-profit purposes, without fee,
and without a written agreement is hereby granted, provided that the above
copyright notice, this paragraph and the following three paragraphs appear
in all copies.
 
Those desiring to incorporate this CLePAPS Software into commercial products
or use for commercial purposes should contact Institute of Theoretical Physics,
Academia Sinica,Beijing 100190, China. 
EMAIL: zheng@itp.ac.cn, wangsheng@itp.ac.cn
 
IN NO EVENT SHALL ITP(INSTITUTE OF THEORETICAL PHYSICS) BE LIABLE TO ANY PARTY 
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS CLEPAPS SOFTWARE, EVEN IF ITP HAS 
BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
THE CLEPAPS SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND ITP HAS NO 
OBLIGATION TO PROVIDE MAINTENANCE,SUPPORT, UPDATES, ENHANCEMENTS, OR 
MODIFICATIONS.  ITP MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF 
ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT
THE USE OF THE CLEPAPS SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR
OTHER RIGHTS.
*/

#pragma once
#include "PDB_Residue.h" 
#include "PDB_Utility.h"
#include "Computation_Utility.h" 
#include <vector>


//====class: Mol_Ori====//
//=> input structure (in XYZ_Format)
class Mol_Ori
{
public:
	Mol_Ori(void);
	virtual ~Mol_Ori(void);

//=========================== Database_Input =================================//
public:
	//ind[*]
	//ind[0]   -> chain
	//ind[1-4] -> ResNum
	//ind[5]   -> InsCode
	vector <string> out_map;  // string <13>
	string MOL_ROOT,XYZ_FROM_PDB,DSSP_FROM_PDB;
	int OutPrt;
	int BADorNOT;
	int TERorNOT;
	int name_range;
	//--- memory_limit ---//__130830__//
	int PRE_LOAD;     // default 0; only valid for Mol_File and Mol_Load
	int WARNING_out;  // default 1; printf warning information to stderr
	int MEMORY_LIMIT; // default (PROT_MAX_NUM)

	//--PART_II:  input_real_process
	void range_init(void);
	void range_delete(void);
	int parse_string(string &buf,int &wstag,string &name,string &range);
	int parse_file(string &in,string &file,int wstag);
	int parse_range(const string &in,vector <string> &out_map);
	int check_range_string(const string &in);
	int check_range_char(char in);

	//--PART_III: batch_process
	int parse_file_XYZ(string &in,string &file,int wstag);
	int Input_XYZ_FULL(string &fn,char chain,int PS,int mode,
		int st,char stt,int ed,char edd,int &totnu,
		XYZ *mol,char *AMI,char *CLE,char *ind);
	int XYZ_Input(string &file,const string &range,int form,int &moln,
		XYZ *mol,char *AMI,char *CLE,char *ind,PDB_Residue *pdb=0);	
	int XYZ_Tranform(string &buf,int &lon,char *nam,
		XYZ *mol,char *AMI,char *CLE,char *ind,PDB_Residue *pdb=0);
	int XYZ_Tranform_BATCH(string &list,int &tot,int *lon,char **nam,
		XYZ **mol,char **AMI,char **CLE,char **ind,PDB_Residue **pdb=0);

	//--PART_IV: DSSP_process
	int parse_file_DSSP(string &in,string &file,int wstag);
	int Input_DSSP_FULL(string &fn,char chain,int PS,int mode,
		int st,char stt,int ed,char edd,int &totnu,
		XYZ *mol,char *AMI,char *SSE,char *ind,
		int *acc,int *beta,double *tco,double *kappa,
		double *alpha,double *phi,double *psi);
	int DSSP_Input(string &file,string &range,int form,int &moln,
		XYZ *mol,char *AMI,char *SSE,char *ind,
		int *acc,int *beta,double *tco,double *kappa,double *alpha,double *phi,double *psi);
	int DSSP_Tranform(string &buf,int &lon,char *nam,
		XYZ *mol,char *AMI,char *SSE,char *ind,
		int *acc,int *beta,double *tco,double *kappa,double *alpha,double *phi,double *psi);
	int DSSP_Tranform_BATCH(string &list,int &tot,int *lon,char **nam,
		XYZ **mol,char **AMI,char **CLE,char **ind,
		int **acc,int **beta,double **tco,double **kappa,double **alpha,double **phi,double **psi);


//-------- virtual_reloaded_function --------//
public:
	virtual int Upload_Process(string &file,char chain,char chain_bak,int PDBorSEQ,int tag,
		int head_tag,int head_c,int tail_tag,int tail_c,int &moln,
		XYZ *mol,char *AMI,char *CLE,char *ind,PDB_Residue *pdb=0);
	virtual int Database_Process_XYZ_I(string &name,int wstag);
	virtual int Database_Process_XYZ_II(char chain,char chain_bak,int PDBorSEQ,int tag,
		int head_tag,int head_c,int tail_tag,int tail_c,int &moln,
		XYZ *mol,char *AMI,char *CLE,char *ind);
	virtual int Database_Process_DSSP_I(string &name,int wstag);
	virtual int Database_Process_DSSP_II(char chain,char chain_bak,int PDBorSEQ,int tag,
		int head_tag,int head_c,int tail_tag,int tail_c,int &moln,
		XYZ *mol,char *AMI,char *CLE,char *ind,
		int *acc,int *beta,double *tco,double *kappa,
		double *alpha,double *phi,double *psi);
};
//====class: Mol_Ori====//over
