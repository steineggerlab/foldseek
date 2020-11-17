#pragma once
#include "Ali_Ori.h"



//====class: Ali_Ali3====//
//=> one-to-multi correspondence (ali3)
class Ali_Ali3 : virtual public Ali_Ori
{
public:
	Ali_Ali3(int num=PROT_MAX_NUM);
	~Ali_Ali3(void);
	int Ali3_Maximal;

//Ali3:DataStructure//
//	ali3[i][*][*] the i-th index of mol2(smaller)
//	ali3[i][0][0] the total correspondece of i-th index(default:0)
//	ali3[i][0][1] the minimal correspondence of i-th index(default:0,then k)
//  ali3[i][k][*] the k-th correspondence of the i-th index
//	ali3[i][k][0] the k-th correspondence of position in mol1(bigger)
//	ali3[i][k][1] the k-th correspondence of distance in mol1(bigger)(integer)
public:
	int *ali3;             //ali3 universal correspondence
	int ali3_TOT;          //ali3's volumn (default:6)
//pre_sco:DataStructure//
//	pre_sco[i][*][*] the i-th index of mol2(smaller)
//	pre_sco[i][0][*] the null state
//	pre_sco[i][k][*] the k-th state
//	pre_sco[i][x][0] the dynamic programming score
//	pre_sco[i][x][1] the correspondence in mol1
public:
	int *pre_sco;          // ali3 dynamic programming data_structure
	int *pre_temp;         // ali3 dynamic programming data_structure(temp)
public:
	int BAK_GAP;           // default 50   (Ali3_DynaProg Bak_Gap penalty)
	int FOR_GAP;           // default 200  (Ali3_DynaProg For_Gap penalty)
	int SCALE;             // default 100  (Ali3_DynaProg Scaling Factor)


//====== ali3 function_related ======//
public:
	//--create & init--//
	void Ali3_Create(int totlen);                           // ali3_Create
	void Ali3_Delete(void);                                 // ali3_Delete
	void Ali3_Init(int *ali3,int moln2);				    // ali3_Init (complex)
	//---- ali3_kabsch ----//
	double ali3_kabsch(double *rotmat_,double &rmsd,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali3);                                         // ali3_rotation (complex)
	//--main_usage--//
	void Ali3_Add(int ii,int jj,int len,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali3,int *zali1=0,int *zali2=0);               // ali3_Add (complex)
	void Ali3_Cor_Add(int *AFP_Cor,int Range,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali3,int *zali1=0,int *zali2=0);               // ali3_Cor_Add (complex)
	void Ali3_To_Ali(double *rotmat_,int method,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali3,int *ali1,int *ali2);                     // ali3-> ali1&2 (complex)
	//--output_part--//
	void Ali3_Out(FILE *fp,int moln2,int *ali3);            // ali3_Out (complex)


//---- ali3 Dynamic_Programming--//
public:
	//--main_part--//
	//vector_score
	int Ali3_Vector_Score_Forward(int ii,int jj,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int SCALE);                                                 // ali3_Vector_Forward_score (complex)
	int Ali3_Vector_Score_Backward(int ii,int jj,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int SCALE);                                                 // ali3_Vector_Backward_score (complex)
	//ali3_DynaProg
	void Ali3_DynaProg(double *rotmat_,double beta,double d_0,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int FOR_GAP,int BAK_GAP,int SCALE,int *ali3);               // ali3_Dynamic_Programming (complex)
	//ali3_TraceBack
	double Ali3_TraceBack(double *rotmat_,double d_0,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali3,int *ali1,int *ali2);                             // ali3_Trace_Back (complex)
	//--output_part--//
	void Ali3_DynaProg_Out(FILE *fp,int moln2,int *ali3);           // ali3_DP_Out (complex)
};
//====class: Ali_Ali3====//over
