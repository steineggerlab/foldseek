#include "BLOSEN_CLE.h"

BLOSEN_CLE::BLOSEN_CLE(int num,int CLESUM)
:Ali_Ori(num),TM_align(num),
 Bioinfo_Code(CLESUM),BLOSEN(num)
{
	//parameter
	SFP_LEN_L=6;
	SFP_LEN_H=9;
	SFP_THRES_L=0;
	SFP_THRES_H=10;
	SFP_Strategy=4;
}
BLOSEN_CLE::~BLOSEN_CLE(void)
{
}

//------------ create & delete --------------//
void BLOSEN_CLE::BLOSEN_CLE_Create(int num,int len)
{
	NewArray2D(&INPUT_CLE,num,len+1);
	NewArray2D(&INPUT_DUM,num,len);
}
void BLOSEN_CLE::BLOSEN_CLE_Delete(int num,int len)
{
	DeleteArray2D(&INPUT_CLE,num);
	DeleteArray2D(&INPUT_DUM,num);
}


//--------- CLE_Generate ---------//
void BLOSEN_CLE::Input_CLE_Generate(void)
{
	int i;
	Confo_Lett confo_lett;
	for(i=0;i<INPUT_NUM;i++)
	{
		confo_lett.btb_ori(0,0,0,INPUT_LEN[i],INPUT_MOL[i],INPUT_CLE[i]);
		INPUT_CLE[i][INPUT_LEN[i]]='\0';
		AMI_CLE_transform_Ori(INPUT_AMI[i],INPUT_CLE[i],INPUT_DUM[i]);
	}
}
void BLOSEN_CLE::Other_Init_Function(void)
{
	Input_CLE_Generate();
}

//---------------- process ------------//
double BLOSEN_CLE::HSFB_Score_Function(int p1,int p2,int ii,int jj,int winlen)
{
	double score;
	score=Universal_Calc(ii,jj,winlen,INPUT_DUM[p1],INPUT_DUM[p2],INPUT_LEN[p1],INPUT_LEN[p2],SFP_Strategy);
	return score;
}
void BLOSEN_CLE::SFP_Generate_Pair(int p1,int p2,vector <SFP_Record> &SFP_H,vector <SFP_Record> &SFP_L) //-> pairwise SFP criteria
{
	int ww1,ww2;
	ww1=INPUT_LEN[p1];
	ww2=INPUT_LEN[p2];
	int num1,num2;
	SFP_L.resize(ww1*ww2);
	SFP_H.resize(ww1*ww2);
	if(ww1>=ww2)
	{
		Universal_SFP_Seed_II(SFP_L,SFP_H,SFP_LEN_L,SFP_LEN_H,SFP_THRES_L,SFP_THRES_H,
			(SFP_LEN_H-SFP_LEN_L)/2,SFP_LEN_H-SFP_LEN_L/2,SFP_LEN_L/2,SFP_LEN_H-SFP_LEN_L/2,
			INPUT_DUM[p1],INPUT_DUM[p2],ww1,ww2,0,SFP_Strategy,num1,num2);   //use amino acid information
	}
	else
	{
		Universal_SFP_Seed_II(SFP_L,SFP_H,SFP_LEN_L,SFP_LEN_H,SFP_THRES_L,SFP_THRES_H,
			(SFP_LEN_H-SFP_LEN_L)/2,SFP_LEN_H-SFP_LEN_L/2,SFP_LEN_L/2,SFP_LEN_H-SFP_LEN_L/2,
			INPUT_DUM[p2],INPUT_DUM[p1],ww2,ww1,1,SFP_Strategy,num1,num2);   //use amino acid information
	}
	SFP_L.resize(num1);
	SFP_H.resize(num2);
	stable_sort(SFP_L.rbegin(),SFP_L.rend());   //descending sort
	stable_sort(SFP_H.rbegin(),SFP_H.rend());
}

//--------------- parameters --------------//
void BLOSEN_CLE::BLOSEN_CLE_Parameter(int H,int sc,int ac,int C1,int C2,int K,int S)
{
	HSFB_LEN=H;         // 12  -> HSFB length
	PIVOT_MOL_CUT=sc;   // 25  -> pivot cutoff
	ANCHOR_HSFB_CUT=ac; // 1   -> HSFB cutoff
	COLOR1_HSFB_CUT=C1; // 5   -> anchor cutoff
	COLOR2_HSFB_CUT=C2; // 5   -> unanchor cutoff
	ZM_TopK=K;          // 5   -> pairwise cutoff
	CONSEN_STAGE=S;     // 3   -> consensus strategy
}
