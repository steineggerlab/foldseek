#pragma once
#include "XYZ.h"
#include "Computation_Utility.h"

//====class: Confo_Lett====//
//=> Conformational_Letter generating
class XYZ_DBT
{
//-----the following is the declaration of functions
public:
	XYZ_DBT(void);
	~XYZ_DBT(void);

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

//--------------------main_function--------------------//
public:
	//function
	void xyz_to_dbt(double *dist,double *bend,double *tort,int n,XYZ *r);
	void dbt_to_xyz(double *dist,double *bend,double *tort,int n,XYZ *r);
};
