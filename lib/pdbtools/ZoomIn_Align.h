#pragma once
#include <vector>
#include "Ali_Ali3.h"
#include "Ali_AFP.h"
#include "TM_align.h"
#include "Align_Utility.h"
using namespace std;

//====class: CLEFAPS====//
//=> CLeFAPS main algorithm
class ZoomIn_Align : public Ali_Ali3, public Ali_AFP, virtual public TM_align
{
public:
	ZoomIn_Align(int num=PROT_MAX_NUM);
	~ZoomIn_Align(void);
	int ZoomIn_Maximal;

//--- macros ---//
public:
	int TMSorNOT;          // default 1 (apply TMscore dynamic programming)
	int HEAD_Chk;          // default 1 (apply head check!!)
//--- parameter ---//
public:
	int ZOOM_ITER;         // default 3    (Zoom-In number)
	int REFINE_ITER;       // default 10   (Optimization number)
	int ANG_CUT;           // default 60   (Vect_Elong's degree cutoff)
	double FIN_CUT;        // default [5.0]  => final distance cutoff
	double INI_CUT;        // default [10.0] => initial distance cutoff
	int ZM_TopK;           // default 20
	int ZM_TopJ;           // default Maximal

//---- variables ----//
public:
	//head fragment
	double *center_mat;     // SFP_High rotmat
	int *head_score;        // SFP_High score
	int *head_index;        // SFP_High index 
	vector <int> valid_idx; // SFP_Low index
	vector <int> align_rec; // record alignment //__110121__//
	//fragment & correspondence	
	int *tali1,*tali2;      // tali1,tali2 (notmal correspondence) [temp]
	int *AFP_Cor;           // AFP_Cor (rigid frag correspondence)  [once]
	int *AFP_Cor_temp;      // AFP_Cor_temp (rigid frag correspondence) [once]
	//inner_iteration
	double *rot_mat;        // rot_mat	
	double *rot_bak;        // rot_bak
	//best record
	int *ZM_BEST_ALI;       // best alignment
	double *ZM_FINMAT;      // best rotmat (mol2 is fixed)
	double *ZM_CURMAT;      // current rotmat (mol2 is fixed)

//---- functions ----//
public:
	//[init function]
	void Init_ZoomIn_Align(int maxlen);
	void Dele_ZoomIn_Align(void);
	void ZM_Input_Mol(XYZ *m1,XYZ *m2,int n1,int n2);
	void CLeFAPS_Init(void);
	//[process minor]
	int ZoomIn_Add(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP);
	int Refinement(int recur,int Range,int ori_sco,double FIN_CUT,double DP_BETA);
	void Elongation(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2);
	void Kill_Bad(double thres,double cutoff,double *rotmat_);
	void Final_Check(double FIN_CUT,double ANG_CUT);
	//[process misc]
	int Select_Single(int cur_SFP,int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP);
	int Select_Best_Pivot(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,int TopJ);
	double T_Similar(double *rotmat_,int *ali1,int *ali2);
	int Comparison_Alignment(int *ali2_cur,int *ali2_comp,int length);
	int Possessed_Check(int ii,int jj,int winlen,int *ali1,int *ali2,int moln1,int moln2);
	int Possessed_Check(int ii,int jj,int winlen,vector <int> &alignment,int moln1,int moln2);
	//[process main]
	double CLeFAPS_Part(XYZ *m1,XYZ *m2,int n1,int n2,double *rotmat_,vector <SFP_Record> &SFP_Low,int *ali2_); //given initial alignment
	double CLeFAPS_Full(XYZ *m1,XYZ *m2,int n1,int n2,vector <SFP_Record> &SFP_High,vector <SFP_Record> &SFP_Low,
		int *ali2_,int *ali2_comp=0);   //given initial alignment
};
