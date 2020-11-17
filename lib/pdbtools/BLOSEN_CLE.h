#pragma once
#include <vector>
#include "BLOSEN.h"
#include "Confo_Lett.h"
#include "Bioinfo_Code.h"
using namespace std;

//====class: BLOSEN_CLE====//
//=> BLOSEN + CLE
class BLOSEN_CLE : public BLOSEN, virtual public Bioinfo_Code
{
public:
	BLOSEN_CLE(int num=PROT_MAX_NUM,int CLESUM=1);
	~BLOSEN_CLE(void);

//--- parameter_set ----//
public:
	int SFP_LEN_L;   //for small-length SFP_L -> 6
	int SFP_LEN_H;   //for large-length SFP_H -> 9
	int SFP_THRES_L; //for small-length SFP_L -> 0
	int SFP_THRES_H; //for large-length SFP_H -> 10
	int SFP_Strategy;//consider both AMI+CLE  -> 4

//--- data_structure ---//
public:
	char **INPUT_CLE;    //input CLE
	int **INPUT_DUM;     //input dummy

//--- process_function --//
public:
	//init & input
	void BLOSEN_CLE_Create(int num,int len);
	void BLOSEN_CLE_Delete(int num,int len);
	//[minor]_CLE related
	void Input_CLE_Generate(void);
	void Other_Init_Function(void);
	//[minor]_part_virtual
	double HSFB_Score_Function(int p1,int p2,int ii,int jj,int winlen); //-> HSFB criteria
	void SFP_Generate_Pair(int p1,int p2,vector <SFP_Record> &SFP_H,vector <SFP_Record> &SFP_L); //-> pairwise SFP criteria
	//parameter input
	void BLOSEN_CLE_Parameter(int H,int sc,int ac,int C1,int C2,int K,int S);
};
