#pragma once
#include "XYZ.h"
#include "Computation_Utility.h"



//====class: Hydro_Bond====//
//=> calculate hydrogen bond
class Hydro_Bond
{
public:
	Hydro_Bond(int num=PROT_MAX_NUM);
	~Hydro_Bond(void);
	int HB_MAXIMAL;     //default: PROT_MAX_NUM

//---- data_structure ---//
public:
	//hydrogen_bond
	int HB_moln;        //mol length
	XYZ **HB_mol;       //input mol (N,CA,C,O,CB) [backbone+CB]
	char *HB_ami;       //input ami (important!!proline has no backbone H-atom!) //__100718__//
	char *HB_sse;       //output sse (DSSP)
	int *HB_broke;      //broken chain
	int *HB_simple;     //hydro bond matrix [n*2] (each res should have 2 HB at most)
	int *HB_index;      //beta sheet index  [n*5] (each res should have 4 sheet at most)
	int *HB_bridge;     //beta sheet bridge [n*2]
	int *HB_ladder;     //beta sheet ladder [n*2*4]
	int *HB_ladder_rec; //beta sheet ladder index [n*2]
	int *HB_temp_index; //beta sheet temp index [n]

//---- main_function ---//
public:
	//init
	void HB_Init(int HB_MAXIMAL);
	void HB_Dele(int HB_MAXIMAL);
	//calc Hydro_Bond
	void HB_Input_Mol(XYZ **mol,char *ami,int moln); //intput
	double HB_Calc_Single(XYZ **HB_mol,int i,int j); //calc_single (from i -> j)
	void HB_Calc_Hydro_Bond(XYZ **HB_mol,int HB_moln,int *HB_simple);  //calc_main
	void HB_Calc_Hydro_Bond(void);                   //normal usage
	void HB_Calc_Hydro_Bond(vector <vector <double> > &hb_mat); //-> output hydro_bond matrix //__180520__//
	//calc Secondary_Structure
	//[helix and turn]
	int HB_Broke_Check(int i,int j);                   //check broke
	void HB_Calc_SSE_Helix(int pitch,char h,char t);   //calc_helix(turn)
	//[beta sheet]
	void HB_Calc_SSE_Sheet(char e,char h);             //calc_sheet
	void HB_Calc_SSE_Sheet_Ladder(void);               //calc_sheet(ladder)
	void HB_Calc_SSE_Sheet_Bulge(char ef);             //calc_sheet(bulge)
	void HB_Calc_SSE_Sheet_Lone(char ef,char eb);      //calc_sheet(lone)
	void HB_Calc_SSE_Sheet_Return(char ef,char e);     //calc_sheet(return)
	//[kappa turn]
	void HB_Calc_SSE_Turn(char s);                     //calc_turn
	void HB_Calc_SSE_Lone(char ori,char t,int limit);  //kill_single_helix
	//[main]
	void HB_Calc_SSE(char *sse);                       //calc_main
	char HB_Trans_SSE_Single_A(char c);                //8-digit to 3-digit single, method A
	char HB_Trans_SSE_Single_B(char c);                //8-digit to 3-digit single, method B
	void HB_Trans_SSE(char *in,char *out,int moln,
		int method=0);                             //8-digit to 3-digit //__100820__//
};

