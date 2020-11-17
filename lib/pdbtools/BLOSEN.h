#pragma once
#include <vector>
#include "BLOMAPS_Ori.h"
#include "ZoomIn_Align.h"
#include "MultiAlign_Cent.h"
using namespace std;

//--- updated BLOMAPS ---//
//-> using ZoomIn_Align to deal with pairwise
//-> using MultiAlign_Cent to deal with final MSA
//-> the framework is BLOMAPS_Ori
class BLOSEN : public BLOMAPS_Ori, public ZoomIn_Align, public MultiAlign_Cent
{
public:
	BLOSEN(int num=PROT_MAX_NUM);
	~BLOSEN(void);
	int BLOSEN_Maximal;


//--- parameter_set ----//
public:
	double SHAVE_THRES;  //for shaving threshold,      default->0.5
	double COMP_THRES;   //for MSA comparison, default->4
	double CORE_THRES;   //for MSA comparison, default->0.5

//--- data_structure ---//
public:
	//[SFP]
	vector <vector <SFP_Record> > BN_SFP_H;  //for pairwise process (SFP_High)
	vector <vector <SFP_Record> > BN_SFP_L;  //for pairwise process (SFP_Low)
	vector <SFP_Record> TMP_SFP_H;           //for initial-MSA construct (SFP_High)
	//[HSFB]
	int *Fast_Sort_Val;    //just for fast_sort
	int *Fast_Sort_Idx;    //just for fast_sort
	int **Shav_Grid;       //for HSFB shaving
	//[Anchor]
	XYZ **INPUT_MOL_Rot;   //for Anchor-HSFB
	XYZ *BN_tmp1,*BN_tmp2; //for Anchor-HSFB
	int *Pivot_Ali_TTT;    //temporary pivot_ali
	int *Check_ali1;       //for better search
	int *Check_ali2;       //for better search


//--- process_function --//
public:
	//init & input
	void BLOSEN_Create(int num,int len);
	void BLOSEN_Delete(int num,int len);

	//[utility]
	double Compare_Two_MSA(vector < vector <int> > &msa1,vector < vector <int> > &msa2,int len1,int len2,int totnum);

	//[minor]_part_virtual
	virtual double HSFB_Score_Function(int p1,int p2,int ii,int jj,int winlen)=0; //-> HSFB criteria
	virtual void SFP_Generate_Pair(int p1,int p2,vector <SFP_Record> &SFP_H,vector <SFP_Record> &SFP_L)=0;//-> pairwise SFP criteria
	//[minor]_ordinary
	void Pivo_Ancho_Rotate_Single(int pivot,int target,int *In_AFB,int In_Len,double *ret_rot);
	void Pivo_Ancho_Rotate(int pivot,int *In_AFB,int In_Len,double **ret_rot);
	int HSFB_Anchor_List_Single(int pivot,int target,int Sele_Num,int **Sele_AFB,int *Sele_Len,double *ret_rot);

	//[process]_main
	//0.Additional Function
	void Other_Pivot_Function(int pivot);
	//1.HSFB
	void HSFB_Generate_Part(int pivot,int winlen,int &Out_Num,int **Out_AFB,int *Out_Len);
	void HSFB_Shaving(int In_Num,int **In_AFB,int *In_Len,int &Out_Num,int **Out_AFB,int *Out_Len);
	//2.Pivot-select
	int Repre_Mol_Select(int num,int *len,XYZ **mol,char **ami,int *ret_set);
	int Pivot_Mol_Select(int In_Num,int **In_AFB,int *In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,int *ret_set);
	//3.Anchor-HSFB
	int HSFB_Sort_List_Generate(int pivot,int In_Num,int **In_AFB,int *In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,int *ret_set);
	int HSFB_Anchor_List_Generate(int pivot,int *In_AFB,int In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,double cutoff,double **ret_rot);
	//4.pariwise refine
	double Pairwise_Refine(int p1,int p2,double *ini_rot,int *ali2,double *ret_rot);
	double Pairwise_Total(int p1,int p2,int *ali2_comp,int *ali2,double *ret_rot);
	double Refine_HSFB_Generate_Part(int pivot,int **pivot_ali,double cutoff,int *ali_out);
	double Deal_Anchored(int pivot,int target,int color_cut,
		int In_Num,int **In_AFB,int *In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,int *pivot_ali,int *ret_ali,double *ret_rot);
	//5.consensus refine
	double Terminal_Refine_Main(int pivot,double **in_rot,XYZ **in,int *len,int totnum,
		int **alin,int totin,int **aliout,int &totout,double **ret_rot,int &core_len);
};
