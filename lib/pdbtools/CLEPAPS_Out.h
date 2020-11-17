#pragma once
#include "Mol_Out.h"
#include "Align_Utility.h"
#include <algorithm>
#include <vector>
using namespace std;

//====class: CLEPAPS====//
//=> CLEPAPS realize
// Using Zheng Wei-Mou's method
class CLEPAPS_Out : virtual public Mol_Out
{
public:
	CLEPAPS_Out(int num=PROT_MAX_NUM);
	~CLEPAPS_Out(void);

//---- variables ----//
public:
	//parameter
	double outcut1;
	double outcut2;
	int outmaxnum;
	//original
	int CH_moln1;
	int CH_moln2;
	XYZ *CH_mol1;
	XYZ *CH_mol2;
	PDB_Residue *CH_pdb1;
	PDB_Residue *CH_pdb2;
	//data
	int *CH_ali1;
	int *CH_ali2;
	int *CH_AFP_Cor;
	XYZ *CH_temp;
	PDB_Residue *CH_pdbt;
	//structure (pointer)
	char *CH_CLE1;
	char *CH_CLE2;
	char *CH_AMI1;
	char *CH_AMI2;
	char *CH_IND1;
	char *CH_IND2;
	//name
	string CH_NAM1;
	string CH_NAM2;

//---- functions ----//
public:
	//[input]
	void CLEPAPS_Input_Func(char *c1,char *c2,char *a1,char *a2,char *ind1,char *ind2,
		PDB_Residue *p1,PDB_Residue *p2,XYZ *m1,XYZ *m2,int n1,int n2,string &nam1,string &nam2);
	//[vice]
	int CLEPAPS_Out_RMSD(XYZ *m1,XYZ *m2,int len,double *rotmat_);
	int CLEPAPS_Out_Script_Simp(Align_Record &align,int *AFP_Cor,int moln1,int moln2);
	int CLEPAPS_Out_Script_Ori(Align_Record &align,int *AFP_Cor,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,double *rotmat_);
	void CLEPAPS_Out_Ali2_To_Ali1(int moln1,int moln2,int *ali1,int *ali2);
	void CLEPAPS_Out_rot_point(XYZ p_old,XYZ &p_new,double *rotmat_);
	void CLEPAPS_Out_rot_mol(XYZ *m_old,XYZ *m_new,int nAtom,double *rotmat_);
	void CLEPAPS_Out_rot_pdb(PDB_Residue *m_old,PDB_Residue *m_new,int nAtom,double *rotmat_);
	//[output]
	void CLEPAPS_Output_Script(FILE *fp,string &name,int *AFP_Cor);        //output mol_script
	void CLEPAPS_Output_MolCA(FILE *fp,Align_Record &align);               //output mol_CA
	void CLEPAPS_Output_Align(FILE *fp,Align_Record &align,int FullOrNot); //output alignment
	void CLEPAPS_Output_MolPDB(FILE *fp,Align_Record &align);              //output mol_CA
	//[main]
	int CLEPAPS_Main_Output(string &name,string &script,vector <Align_Record> &tot,
		int SingleOrNot=1,int FullOrNot=1,int ScriptOrNot=0);
};

