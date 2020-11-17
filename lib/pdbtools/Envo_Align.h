#pragma once
#include "Confo_Lett.h"
#include "TM_align.h"
#include "Bioinfo_Code.h"
#include "Dynamic_Programming.h"

class Envo_Align : virtual public TM_align, virtual public Bioinfo_Code
{
public:
	Envo_Align(int num=PROT_MAX_NUM,int CLESUM=1);
	~Envo_Align(void);
	int Envo_Align_maximal;

//---- variables ----//
public:
	//macro
	int Envo_Align_Do_TMalign;       // whether use TMalign to refine (default: no)
	//normal parameter
	double Envo_Align_Distance;      // the distance cutoff (default: 8.5)
	int Envo_Align_mol_extend;       // determine the initial rotmat (default: 5)
	int Envo_Align_neib_extend;      // determine the addtional neibor (default: 0)
	//DynaProg parameter
	double Envo_Align_AMI_GapOpen;   // the gap open penalty for AMI
	double Envo_Align_AMI_GapExtend; // the gap extend penalty for AMI
	double Envo_Align_CLE_GapOpen;   // the gap open penalty for CLE
	double Envo_Align_CLE_GapExtend; // the gap extend penalty for CLE
	//input
	int Envo_Align_moln1;
	int Envo_Align_moln2;
	XYZ *Envo_Align_mol1;
	XYZ *Envo_Align_mol2;
	int *Envo_Align_ami1;
	int *Envo_Align_ami2;
	int *Envo_Align_cle1;
	int *Envo_Align_cle2;
	//temp
	vector <vector <int> > Envo_Align_Neib1;
	vector <vector <int> > Envo_Align_Neib2;
	double *Envo_Align_DP_sco;  //DynaProg matrix
	XYZ *Envo_Align_tmp1;
	XYZ *Envo_Align_tmp2;
	int *Envo_Align_alit;
	int *Envo_Align_ali2;
	char *Envo_Align_cle;

//---- functions ----//
public:
	//init function
	void Init_Envo_Align(int maxlen);
	void Dele_Envo_Align(void);
	//vice function
	void Envo_Align_Calc_Neib(XYZ *mol,int moln,vector <vector <int> > &Neib,double rcut);
	void Envo_Align_Change_Alignment(vector<pair<int,int> > &alignment,int size2,int *out_ali);
	double Envo_Align_Refine_Alignment(int ii,int jj,int *out_ali);
	//misc function
	double Envo_Align_Calc_AMI(int ii,int jj,int *out_ali); //ami dynaprog
	double Envo_Align_Calc_CLE(int ii,int jj,int *out_ali); //cle dynaprog
	double Envo_Align_Calc_mol(int ii,int jj,int *out_ali); //mol dynaprog
	//main function
	void Input_Envo_Align(XYZ *mol1,XYZ *mol2,int moln1,int moln2,char *ami1,char *ami2);
	void Input_Envo_Align_Total(XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		char *ami1,char *ami2,char *cle1,char *cle2);
	void Envo_Align_Neib_Main(XYZ *mcb1,XYZ *mcb2,int moln1,int moln2);
	void Envo_Align_Calc_Main(vector <vector <double> > &out_mat,int Out_Type);
	//Out_Type: [0] mol, [1] ami, [2] cle
};
