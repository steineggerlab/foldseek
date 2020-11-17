#pragma once
#include "PDB_Utility.h"
#include "CLEFAPS.h"
#include "Bioinfo_Code.h"
#include <algorithm>
#include <vector>
using namespace std;

//====class: CLEPAPS====//
//=> CLEPAPS realize
// Using Zheng Wei-Mou's method
class CLEFAPS_CLE : public CLEFAPS, virtual public Bioinfo_Code
{
public:
	CLEFAPS_CLE(int num=PROT_MAX_NUM,int CLESUM=1);
	~CLEFAPS_CLE(void);

//---- variables ----//
public:
	//macro
	int SFP_Strategy;    //-> using CLE,AMI,GEN or CLE+AMI to generate SFP
	int REF_Strategy;    //-> using refinement strategy
	//parameter
	int SFP_H_Len;
	int SFP_L_Len;
	int SFP_H_Thres;
	int SFP_L_Thres;
	//additional
	int ZM_Upper_Max;    //-> default is 500
	int ZM_Lower_Max;    //-> default is 50
	int ZM_TopL;         //-> retain TopL candidate for TMscore refinement
	double Refine_Cut;   //-> TM_refine cutoff (default: 0.3)
	//structure
	char *IN_CLE1;
	char *IN_CLE2;
	char *IN_AMI1;
	char *IN_AMI2;
	int *IN_DUM1;
	int *IN_DUM2;

//---- functions ----//
public:
	//[input]
	void CLEFAPS_Input_Func(char *c1,char *c2,char *a1,char *a2,
		XYZ *m1,XYZ *m2,int n1,int n2);
	void CLEFAPS_Input_Para(int h_len,int l_len,int h_thres,int l_thres,int Stragety=1);
	//[process]
	void CLEFAPS_Second_Init(void);
	double CLEFAPS_Second_Refine(vector <Align_Record> &tot);
	void CLEFAPS_Second_Final(vector <Align_Record> &tot,int REAL_FAST=0);
	//[main]
	void CLEFAPS_Main_Init(int SFP_Strategy);
	double CLEFAPS_Main_Func(vector <Align_Record> &tot);
	//[mask]
	int Check_Mask(vector <SFP_Record> &in, vector <SFP_Record> &out, int *mas1, int*mas2);
	void Mask_CLE_BreakPoint(char *CLE, int *mask, int len);
};
