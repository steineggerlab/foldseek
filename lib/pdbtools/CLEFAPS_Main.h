#pragma once
#include "CLEFAPS_CLE.h"
#include "Envo_Align.h"
#include "TM_Align_Main.h"
using namespace std;

//====class: CLEPAPS====//
//=> CLEPAPS realize
// Using Zheng Wei-Mou's method
class CLEFAPS_Main :  public Envo_Align, public CLEFAPS_CLE, public TM_Align_Main
{
public:
	CLEFAPS_Main(int num=PROT_MAX_NUM,int CLESUM=1);
	~CLEFAPS_Main(void);

	//---------- data structure ------------//
	int Kill_Frag;    //whether kill frag or not
	int Kill_Frag_Para_0;  //-> min_small_len, 25
	int Kill_Frag_Para_1;  //-> gap_len, 25
	int Kill_Frag_Para_2;  //-> nongap_len, 15
	int FM_Kill_Gap;  //whether kill gap or not
	vector <Align_Record> FM_align_tot;

	//----------- functions -----------//
	//kill fragment
	int Alignment_Gap_Process(int *ali2_in,int moln1,int moln2,int *ali2_out,
		int gap_len,int nongap_len);
	//single match kill
	double IN_Calc_Frag_Score(int ii,int jj,int *ali1,int *ali2,
		int moln1,int moln2,XYZ *mol1,XYZ *mol2,double *rotmat_);
	int IN_Ali_To_Cor(int *AFP_Cor,int thres,int moln1,int moln2,int *ali1,int *ali2,
		vector<pair<int,int> > &ws_pair_record);
	void IN_Elongation(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2);
	void FM_Align_Refine(XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali2_in,int *ali2_out,double *rotmat_);
	//gap series kill
	void Insert_Refine_Single(char *seq1,char *seq2,int start1,int start2,int moln1,int moln2,
		double *ws_in,char *out1,char *out2,int HEADorTAIL);
	void Insert_Refine_Total(char *ali1,char *ali2,double *ws_in,int moln1,int moln2,
		char *out1,char *out2,int CutK);
	void Retreive_Alignment(char *out1,char *out2,int *ali2_,int moln2);
	void Transform_Ali(int *ali1,int *ali2,int moln1,int moln2,char *ami1,char *ami2,char *out1,char *out2);
	void Double_Gap_Refine(char *ali1,char *ali2,double *ws_in,int moln1,int moln2,
		char *out1,char *out2,int CutK,double thres);
	double CLEF_Make_Score(Align_Record & align_record,XYZ *mol1,XYZ *mol2,int moln1,int moln2,double *rotmat,int *ali2);
	void CLEF_Single_Match_Kill(XYZ *mol1,XYZ *mol2,int moln1,int moln2,vector <Align_Record> &tot);
	//main function
	void FM_Get_Initial_AMI(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2); //AMI dynamic programming
	void FM_Get_Initial_CLE(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2); //CLE dynamic programming
	void FM_Get_Initial2(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2); //AMI+CLE dynamic programming
	void FM_Get_Initial3(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2); //CLE envo programming
	void Alimeng_Add(vector <Align_Record> &align_tot,int moln2,double tms);
	double FM_Align_Total(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
		int norm_len,double norm_d0,double *MAXSCO=0); //MAXSCO should at least double[8]
	double FM_Align_Normal(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
		int norm_len,double norm_d0,double *MAXSCO=0); //MAXSCO should at least double[8]
	double FM_Align_WithAli(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
		int norm_len,double norm_d0,double *MAXSCO=0); //MAXSCO should at least double[8]
	double FM_Align_Lite(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
		int norm_len,double norm_d0,double *MAXSCO=0); //MAXSCO should at least double[8]
	double FM_Align_Fast(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,double *MAXSCO=0);  //MAXSCO should at least double[8]
};
