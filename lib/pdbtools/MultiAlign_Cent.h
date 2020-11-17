#pragma once
#include "TM_align.h"

//====class: MultiAlign_Cent====//
//=> input initial multi_ali
//=> output refined multi_ali
class MultiAlign_Cent : virtual public TM_align
{
public:
	MultiAlign_Cent(int tmnum=PROT_MAX_NUM);
	~MultiAlign_Cent(void);

//--- macros ---//
public:
	int TMorLEN;        //default:1, CORE_LEN*TMscore,     (0, TMscore, 2, CORE_LEN)
	int RECorNOT;       //default:1, record initial value  (0, not record)
	int ROTorNOT;       //default:1, final consensus stage (1, not do)
//--- parameter ---//
public:
	int CONSEN_STAGE;   //default:3  -> using total three stragety for refinement
	int REFINE_NUM;     //default:3  -> each stragety, refine how many times
	double MERGE_DIST;  //default:3.0 -> column merge cutoff

//--- data_structure ---//
public:
	//[ori mol]
	int BC_INPUT_NUM;       //input number
	int *BC_INPUT_LEN;      //input length
	XYZ **BC_INPUT_MOL;     //input structure //[this data may be changed!]
	//[consensus]
	XYZ *BLOCEN_Consensus;  //consensus structure
	int *BLOCEN_Consenrec;  //record its original point
	double **BC_Condist;    //distance between consensus
	int *BC_ZeroOne_Cur;    //record current position
	int *BC_ZeroOne_Tmp;    //record current position (temp)
	int *BC_ZeroOne_Len;    //record current length (temp)
	//[refienment]
	XYZ *BC_Tmp1,*BC_Tmp2;  //temp structure
	int **BC_Ali2_Refine;   //use TMalign to refine path
	int **BC_Ali_Temp;      //current path
	//[other temp]
	int *BC_ali1;           //for rotation_version
	int *BC_ali2;           //for rotation_version
	int *BC_alib;           //for rotation_version
	XYZ *BC_best;           //for rotation_version
	//[best record]
	XYZ *BC_Best_Consensus; //record the best consensus
	int BC_Best_Consenlen;  //record the best length
	XYZ **BC_Best_Output;   //superimposed structure set
	int **BC_Ali_Best;      //total best path

//--- process_function --//
public:
	//init
	void MultiAlign_Cent_Create(int num,int len);
	void MultiAlign_Cent_Delete(int num,int len);
	//consensus_vice
	int Single_Ali_To_Multi_Ali(int totnum,int *len,int conlen,int **ali2out,int **aliout);
	int BC_Consen_Refine(XYZ *in,int *rec,int totnum,int totlen,double thres);
	double BC_Calc_SumTM_Given_Mol(XYZ **in,int *len,int totnum,int totlen,int **ali);
	void BC_Get_Initial_Ali2(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2);
	//fixed(for normal refinement)
	int BC_Update_Consen_Given_Ali_Fixed(XYZ **in,int totnum,int totlen,int **ali,XYZ *consen,int *conrec,int cutoff=0);
	int BC_Update_Ali_Gigen_Consen_Fixed(XYZ **in,int *len,int totnum,int conlen,XYZ *consen,int **aliout,int ROTorNOT=0);
	int BC_Return_Conserved_Core(int totnum,int totlen,int **ali);
	void BC_Return_Conserved_Core_Pos(int totnum,int totlen,int **ali,int *core_pos);
	double BC_Return_Conserved_RMSD(int totnum,int totlen,int **ali,XYZ **mol);
	void BC_Return_Conserved_RMSD_Pos(int totnum,int totlen,int **ali,XYZ **mol,double *rmsd_pos);
	//rotation(for advance refinement)
	int BC_Update_Consen_Given_Ali(XYZ **in,int *len,int totnum,int conlen,int *conres,double **conrec,int **aliout,int **ali2out,XYZ *consen);
	double BC_Update_Ali_Gigen_Consen_Single(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *res,double *rec,int *ali2out);
	int BC_Update_Ali_Gigen_Consen(XYZ **in,int *len,int totnum,int conlen,int *conres,double **conrec,int **aliout,XYZ *consen);
	//initial refine
	void BC_Update_Partial_CORE_Single(XYZ **in,int *len,int totnum,int **alin,int totin,int **aliout,int &totout);
	void BC_Update_Partial_CORE(int pivot,XYZ **in,int *len,int totnum,int **alin,int totin,int **aliout,int &totout);
	//main
	double BC_Tota_Update_Function(XYZ **in,int *len,int totnum,int totlen,int **ali); //this could be a virtual function!
	double BC_Tota_Update_Main(XYZ **in,int *len,int totnum,int **alin,int totin,int **aliout,int &totout);
};
