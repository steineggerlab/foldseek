//  CLePAPS:                                                             //
//  Conformational Letter based Pairwise Alignment of Protein Structure  //
//                                                                       //
//  Authors: S.Wang & WM.Zheng                                           //

#pragma once
#include "XYZ.h"
#include "Computation_Utility.h"



//==Definition==//
#ifndef Confo_TriPep
#define Confo_TriPep 3    //peptides   D
#endif
#ifndef Confo_Cluster
#define Confo_Cluster 17  //Confo_Clusters   C 
#endif
//==Definition==//over



//-----extern data_structure ------//
extern int Ori_id[Confo_Cluster];
extern int Ant_id[Confo_Cluster];
extern double Ori_pi[Confo_Cluster];
extern double Ori_sig[Confo_Cluster];
extern double Ori_mid[Confo_Cluster][Confo_TriPep];
extern double Ori_var[Confo_Cluster][Confo_TriPep][Confo_TriPep];



//====class: Confo_Lett====//
//=> Conformational_Letter generating
class Confo_Lett
{
//-----the following is the declaration of functions
public:
	Confo_Lett(void);
	~Confo_Lett(void);

//----------temp_variables----------//
public:
	double cle_x_point[2][3];
	double cle_y_point[2][3];
	double cle_w_point[3];
	double cle_v1[3],cle_v2[3];
	double cle_w1[3],cle_w2[3];
	double cle_x1[3],cle_x2[3];
	double cle_Axis[3],cle_test[3];
	double cle_R_Theta[3][3];
	double cle_R_Thor[3][3];
	double cle_T_Back[3][3];
	double cle_T_Pre[3][3];
	double cle_D_Back[3];
	double cle_D_Pre[3];
	double cle_r1[3][3],cle_r2[3][3],cle_temp[3][3];
	double cle_Z_Axis[3],cle_X_Axis[3];
//---- macro ---//__100225__//
public:
	double CLE_CUT_MIN;  //default:2.5
	double CLE_CUT_MAX;  //default:4.5


//--------------------main_function--------------------//
public:
	//---PART_0: (confo_lett_init)
	void confo_lett_init(void);
	void confo_lett_delete(void);
	//---PART_1:  'XYZ' transform to 'Confo_Angle'(Conformational_Angle)
	double gaussian(int ic,double *xx);
	char btb_stc(double *w_point);
	void profile_stc(double *w_point,double *output);
	void btb_ori(double *dist,double *bend,double *tort,int n,XYZ *r,char *CLE,double **Profile=0);
	//---PART_2: 'Confo_Angle' transform to 'XYZ'
	void Get_Initial(XYZ *in1,XYZ *in2,double r1[3][3],double r2[3][3]);
	void ctc_ori(double *dist,double *bend,double *tort,int n,XYZ *r,XYZ *init);
	//---PART_3: 'Confo_Lett' check
	int Confo_Check_Lett(int code,double *input);
	int Confo_Check_Space(XYZ *bak_mol,XYZ cur_mol,int moln,double cutoff);
};
//====class: Confo_Lett====//over
