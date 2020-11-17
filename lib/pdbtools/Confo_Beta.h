#pragma once
#include "XYZ.h"
#include "Computation_Utility.h"

//====class: Confo_Beta====//
//=> CA<->CB generating
class Confo_Beta
{
public:
	Confo_Beta(int num=PROT_MAX_NUM);
	~Confo_Beta(void);
	int Confo_Beta_Maximal; 

//--- variable ---//
public:
	//normal_variable
	double beta_v1[3],beta_v2[3];
	double beta_x1[3],beta_x2[3];
	double beta_y[3],beta_ca[3],beta_cb[3];
	double Confo_Beta_rotmat[3][3];      //temp tormat;
	//temp_variable
	double *beta_dist;
	double *beta_bend;
	double *beta_tort;
	//main_variable
	double Confo_Beta_Levit_Angle;  //THETA=0.656
	double Confo_Beta_Levit_Radii;  //DIST=1.53

//--- function ---//
public:
	//init
	void Confo_Beta_Init(int totlen);
	void Confo_Beta_Dele(void);
	//process
	void Confo_Beta_CACB_To_Angle(XYZ *CA,XYZ *CB,int moln,double *bend,double *tort,double *dist);
	void Confo_Beta_Angle_To_CACB(XYZ *CA,XYZ *CB,int moln,double *bend,double *tort,double *dist);
	void Confo_Beta_Levit_To_CACB(XYZ *CA,XYZ *CB,int moln,double *dist);
	//main
	void Recon_Beta_1(XYZ *CA,XYZ *CB,int moln,char *ami); //given CA+CB, and AMI_Dist, update SC
	void Recon_Beta_31(XYZ *CA,XYZ *CB,int moln,char *ami); //given CA, return CB (Levitt)
	void Recon_Beta_32(XYZ *CA,XYZ *CB,int moln,char *ami); //given CA, return SC (Levitt)
	void Recon_Beta_21(XYZ *CA,XYZ *CB,int moln,char *ami,char *cle); //given CA, and CLE_Type, return CB (Liu & Wang)
	void Recon_Beta_22(XYZ *CA,XYZ *CB,int moln,char *ami,char *cle); //given CA, and CLE_Type, return CB (Liu & Wang)	
};
