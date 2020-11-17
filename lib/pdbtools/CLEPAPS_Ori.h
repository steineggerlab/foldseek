#pragma once
#include "Align_Utility.h"
#include "Computation_Utility.h"
#include "XYZ.h"
#include <vector>
#include <algorithm>
using namespace std;


//====class: CLEPAPS_Ori====//
//=> CLEPAPS Framwork
class CLEPAPS_Ori
{
public:
	CLEPAPS_Ori(void);
	virtual ~CLEPAPS_Ori(void);

//--- macros ---//
public:
	//global
	int HEAD_Chk;          // default 1 (do head check)
	int TAIL_Chk;          // default 1 (do tail check)
	int COMP_Chk;          // default 1 (do comparison check)
	int FAST_Chk;          // default 0 (DO NOT run fast check)
	//partial
	int ZoomOrNot;         // default 1 (do ZoomIn add)
	int RefineOrNot;       // default 1 (do Refinement)
	int NonLinear;         // default 1 (do Kill Nonlinear)
	int FinalOrNot;        // default 1 (do Shrink or Elong)
	int TMSorNOT;          // default 1 (apply TMscore dynamic programming)
//--- parameter ---//
public:
	int CLEP_MIN_LEN;      // default 4    (final cut length)
	int ZOOM_ITER;         // default 3    (Zoom-In number)
	int REFINE_ITER;       // default 10   (Optimization number)
	double FIN_CUT;        // default [5.0]  => final distance cutoff
	double INI_CUT;        // default [10.0] => initial distance cutoff
	int ZM_TopK;           // default 20
	int ZM_TopJ;           // default 50

//---- variables ----//
public:
	//vector temp
	double *CLEP_rotmat;      // CLEP_rotmat	
	double RMSD;              // RMSD
	int LALI;                 // LALI
	Align_Record CLEP_MAX;    // record the best alignment

//---- virtual functions ----//
public:
	//---- input/output functions ------//io
	virtual void CLEPAPS_Input(XYZ *m1,XYZ *m2,int n1,int n2)=0; 
	//---- make_socre virtual functions ----//score
	//given alignment and rotmat, return scores (maybe TMscore,RMSD and LALI, etc)
	virtual double CLEPAPS_Make_Score(Align_Record & align_record)=0;
	//given the align_record, make the secondary score (maybe LALI*TMscore, etc)
	virtual double CLEPAPS_Get_Score(Align_Record &align)=0;
	
	//---- pure virtual functions -----//partial
	//0. Init Other CLEPAPS
	virtual void Other_CLEFAPS_Init(void)=0;
	//1. ZoomIn Add (ret_rot might be updated)
	virtual int ZoomIn_Add(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,double *ret_rot)=0;
	//2 .Refinement (ret_rot might be updated)
	virtual int Refinement(int recur,double FIN_CUT,double *ret_rot)=0;
	//3. Kill NonLinear (ret_rot won't be updated)
	virtual int Kill_NonLinear(double *ret_rot)=0;
	//4. Final Stage (combine Shrink & Elong)
	virtual int Final_Stage(double *ret_rot)=0;
	//5. TMscore refinement (ret_rot might be updated)
	virtual int Final_Refine(double *ret_rot)=0;

	//---- pure virtual functions -----//global
	//0. Init functions
	virtual void Global_CLEPAPS_Init(void)=0;
	//1. Select best pivot
	virtual int Select_Best_Pivot(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,int TopJ)=0;
	//2. Break Check
	virtual int CLEPAPS_Break_Check(int cur,int count,int totnum)=0;
	//3. Get Index
	virtual int CLEPAPS_Get_Index(int cur)=0;
	//4. Head Check
	virtual int CLEPAPS_Head_Check(SFP_Record &SFP)=0;
	//5. Get initial rotmat
	virtual void CLEPAPS_Initial_Rotmat(int index,double *rotmat)=0;
	//6. Tail Check
	virtual int CLEPAPS_Tail_Check(Align_Record &cur)=0;
	//7. Comparison Check
	virtual int CLEPAPS_Comparison_Check(Align_Record &cur,vector <Align_Record> &tot)=0;
	//8. Terminal Process
	virtual void CLEPAPS_Terminal_Process(Align_Record &cur)=0;
	
	
	//---- main functions -----//
	double CLEPAPS_Part(XYZ *m1,XYZ *m2,int n1,int n2,double *rotmat_,vector <SFP_Record> &SFP_Low,
		Align_Record & align_record);
	double CLEPAPS_Full(XYZ *m1,XYZ *m2,int n1,int n2,vector <SFP_Record> &SFP_High,vector <SFP_Record> &SFP_Low,
		vector <Align_Record> & align_records);
};
