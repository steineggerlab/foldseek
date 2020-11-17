#pragma once
#include <vector>
#include "Ali_Ali3.h"
#include "Ali_AFP.h"
#include "CLEPAPS_Ori.h"
#include "TM_align.h"
using namespace std;

//====class: CLEFAPS====//
//=> CLeFAPS main algorithm
class CLEFAPS : public Ali_Ali3, virtual public Ali_AFP, virtual public CLEPAPS_Ori, virtual public TM_align
{
public:
	CLEFAPS(int num=PROT_MAX_NUM);
	~CLEFAPS(void);
	int CLEF_Maximal;

//---- variables ----//
public:
	//parameter
	int ANG_CUT;            // vector angle threshold
	double Tail_Check;      // tail alignment threshold
	double Comp_Check;      // comparison threshold
	int Comp_Check_Vari;    // comparison partial (default: 0)
	//additional
	double CUR_MaxJ_thres;
	double CUR_MaxJ;        // TopJ maximal
	double CUR_MaxK;        // TopK maximal
	int CUR_INDEX;          // current index
	int INDEX_CUR;          // index current
	double TM_d0;           // d0 to calculate TMscore
	int TM_smaller;         // smaller to calculate TMscore
	//temp
	int *tali1,*tali2;      // tali1 & tali2
	int *AFP_Cor;           // AFP_Cor (rigid frag correspondence)  [once]
	int *AFP_Cor_temp;      // AFP_Cor_temp (rigid frag correspondence) [once]
	double *center_mat;     // SFP_High rotmat
	int *head_score;        // SFP_High score
	int *head_index;        // SFP_High index 
	short *valid_idx;       // SFP_Low index
	vector <int> align_rec; // record alignment //__110121__//

//---- functions ----//
public:
	//[init function]
	void CLEFAPS_Init(int maxlen);
	void CLEFAPS_Dele(void);
	//[vice process]
	double T_Similar(double *rotmat_,int *ali1,int *ali2,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2);
	void Elongation(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2);                             //elongation
	void Shrinking(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2);                             //shrinking
	int Select_Single(int cur_SFP,int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,double *ret_rot);
	double Alg_Comp_Simp(vector <int> &ali_m1,vector <int> &ali_m2,int part=0); //alignment comparison

//---- virtual functions ----//
public:
	//---- input/output functions ------//io
	void CLEPAPS_Input(XYZ *m1,XYZ *m2,int n1,int n2); 
	//---- make_socre virtual functions ----//score
	//given alignment and rotmat, return scores (maybe TMscore,RMSD and LALI, etc)
	double CLEPAPS_Make_Score(Align_Record & align_record);
	//given the align_record, make the secondary score (maybe LALI*TMscore, etc)
	double CLEPAPS_Get_Score(Align_Record &align);
	
	//---- virtual functions -----//partial
	//0. Init Other CLEPAPS
	void Other_CLEFAPS_Init(void);
	//1. ZoomIn Add (ret_rot might be updated)
	int ZoomIn_Add(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,double *ret_rot);
	//2 .Refinement (ret_rot might be updated)
	int Refinement(int recur,double FIN_CUT,double *ret_rot);
	//3. Kill NonLinear (ret_rot won't be updated)
	int Kill_NonLinear(double *ret_rot);
	//4. Final Stage (combine Shrink & Elong)
	int Final_Stage(double *ret_rot);
	//5. TMscore refinement (ret_rot might be updated)
	int Final_Refine(double *ret_rot);

	//---- virtual functions -----//global
	//0. Init functions
	void Global_CLEPAPS_Init(void);
	//1. Select best pivot
	int Select_Best_Pivot(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,int TopJ);
	//2. Break Check
	int CLEPAPS_Break_Check(int cur,int count,int totnum);
	//3. Get Index
	int CLEPAPS_Get_Index(int cur);
	//4. Head Check
	int CLEPAPS_Head_Check(SFP_Record &SFP);
	//5. Get initial rotmat
	void CLEPAPS_Initial_Rotmat(int index,double *rotmat);
	//6. Tail Check
	int CLEPAPS_Tail_Check(Align_Record &cur);
	//7. Comparison Check
	int CLEPAPS_Comparison_Check(Align_Record &cur,vector <Align_Record> &tot);
	//8. Terminal Process
	void CLEPAPS_Terminal_Process(Align_Record &cur);
};
