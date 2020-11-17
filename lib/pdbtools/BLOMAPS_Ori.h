#pragma once
#include "XYZ.h"
#include "Computation_Utility.h"

//====class: BLOMAPS_Ori====//
//=> BLOMAPS Framwork
class BLOMAPS_Ori
{
public:
	BLOMAPS_Ori(void);
	virtual ~BLOMAPS_Ori(void);

//--- parameter_set ----//
public:
	//macro
	int PRITNForNOT; //printf 
	int UNANCHorNOT; //do unanchor
	int HEADCKorNOT; //do headcheck
	//cutoff parameter
	int HSFB_LEN;    //HSFB length, default 12
	double TMS_CUT;  //TM_score cutoff      -> default 0.5
	double TMS_GOD;  //TM_score good pass   -> default: 0.85
	double ANC_CUT;  //for Anchored_MOL threshold -> default:0.0
	double MID_CUT;  //medium refine cutoff -> default 1.0
	//trim parameter
	int PIVOT_MOL_CUT;   //TopM  -> 25
	int ANCHOR_HSFB_CUT; //TopK  -> 1
	int COLOR1_HSFB_CUT; //TopJ1 -> 5
	int COLOR2_HSFB_CUT; //TopJ2 -> 5

//--- data_structure ---//
public:
	//[ori mol] fixed-> cannot modify
	int INPUT_NUM;       //input number
	int *INPUT_LEN;      //input length
	XYZ **INPUT_MOL;     //input structure
	char **INPUT_AMI;    //input amino acid
	//[BLOCKS]
	int Multi_AFB_Num;   //total number
	int Multi_AFB_Num_;  //shave
	int **Multi_AFB;     //number * position
	int **Multi_AFB_;    //shave
	int *AFB_Len;        //AFB length
	int *AFB_Len_;       //shave
	//[Verti+Horiz]
	int *Mol_Anchor;     //anchored MOL
	int *Pivot_Set;      //candidate pivot
	int *Anchor_Set;     //candidate AFB
	//[temp data] flexi-> can modify
	double **TEMP_ROT;   //temporary rotmat
	double *TEMP_TMS;    //temporary TMscore
	int *Pivot_Ali_Temp; //pivot_aligned temp
	int **Pivot_Ali;     //pivot_aligned corset (All-To-One style)
	int **ALIOUT;        //multi_aligned corset (MSA style)
	int TOTLEN;          //multi_aligned totlen (MSA totlen)
	//[real temp]
	int *Real_Temp_Ali2; //for record real temp
	double *Real_Temp_Rot;//real temp rotmat
	//[final out] fixed-> cannot modify
	int FIN_PIVOT;       //final pivot
	double **FIN_ROT;    //final rotmat
	int **ALIBEST;       //multi_aligned best corset (MSA style)
	int BESTLEN;         //multi_aligned bsst totlen (MSA totlen)
	double FIN_SCORE;    //final Sum-TM
	//---- final record ----//
	vector <vector <vector <double > > > FIN_ROT_REC;
	vector <vector <vector <int > > > FIN_ALI_REC;
	vector <int> FIN_BEST_LEN;
	vector <double> FIN_SCO_REC;
	vector <int> FIN_PIV_REC;
	vector <int> FIN_CORE_LEN;

//--- process_function --//
public:
	//init & input
	void Multi_Ali_Create(int num,int len);
	void Multi_Ali_Delete(int num,int len);
	void Multi_Ali_Load(int totlen,int *inlen,XYZ **mol,char **ami);
	void HSFB_Generate_Full(int winlen,int &Out_Num,int **Out_AFB,int *Out_Len);
	//compare two MSA
	virtual double Compare_Two_MSA(vector < vector <int> > &msa1,vector < vector <int> > &msa2,int len1,int len2,int totnum)=0;
	//[process]_part
	//0.Additional Function
	virtual void Other_Init_Function(void)=0;
	virtual void Other_Pivot_Function(int pivot)=0;
	//1.HSFB
	virtual void HSFB_Generate_Part(int pivot,int winlen,int &Out_Num,int **Out_AFB,int *Out_Len)=0;
	virtual void HSFB_Shaving(int In_Num,int **In_AFB,int *In_Len,int &Out_Num,int **Out_AFB,int *Out_Len)=0;
	//2.Pivot-select
	virtual int Repre_Mol_Select(int num,int *len,XYZ **mol,char **ami,int *ret_set)=0;
	virtual int Pivot_Mol_Select(int In_Num,int **In_AFB,int *In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,int *ret_set)=0;
	//3.Anchor-HSFB
	virtual int HSFB_Sort_List_Generate(int pivot,int In_Num,int **In_AFB,int *In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,int *ret_set)=0;
	virtual int HSFB_Anchor_List_Generate(int pivot,int *In_AFB,int In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,double cutoff,double **ret_rot)=0;
	//4.pariwise refine
	virtual double Pairwise_Refine(int p1,int p2,double *ini_rot,int *ali2,double *ret_rot)=0;
	virtual double Pairwise_Total(int p1,int p2,int *ali2_comp,int *ali2,double *ret_rot)=0;
	virtual double Refine_HSFB_Generate_Part(int pivot,int **pivot_ali,double cutoff,int *ali_out)=0;
	virtual double Deal_Anchored(int pivot,int target,int color_cut,
		int In_Num,int **In_AFB,int *In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,int *pivot_ali,int *ret_ali,double *ret_rot)=0;
	//5.terminal
	virtual double Terminal_Refine_Main(int pivot,double **in_rot,XYZ **in,int *len,int totnum,
		int **alin,int totin,int **aliout,int &totout,double **ret_rot,int &core_len)=0;
	//[main]
	virtual double BLOMAPS_Partial(int num,int *len,XYZ **mol,char **ami,int **aliout,double **ret_rot);
};
