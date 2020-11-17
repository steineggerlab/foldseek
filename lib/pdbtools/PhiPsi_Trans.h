#pragma once
#include "XYZ_DBT.h"

//====class: PhiPsi_Trans====//
//=> calculate phi_psi from XYZ
class PhiPsi_Trans
{
public:
	PhiPsi_Trans(int num=PROT_MAX_NUM);
	~PhiPsi_Trans(void);
	int PhiPsi_MAXIMAL;

//---- data_structure ---//
public:
	//const structure
	//[dist+bend]
	double trs_dist_CN;
	double trs_dist_NCa;
	double trs_dist_CaC;
	double trs_bend_CaCN;
	double trs_bend_CNCa;
	double trs_bend_NCaC;
	//[cb+co]
	double Confo_Back_CB_Dis;
	double Confo_Back_CB_Ang;  // Ca(i)->CB bend
	double Confo_Back_CB_Tor;  // Ca(i)->CB tort
	double dist_CO;
	double bend_CO;
	//temp structure
	double *PhiPsi_dist;
	double *PhiPsi_bend;
	double *PhiPsi_tort;
	XYZ *PhiPsi_mol;
	//additional
	double PhiPsi_ROT_TOT[3][3];
	XYZ_DBT xyz_to_dbt;

//---- main_function ---//
public:
	//init
	void PhiPsi_Init(int PhiPsi_MAXIMAL);
	void PhiPsi_Dele(void);
	void PhiPsi_Get_Init(int moln);
	//data	
	void PhiPsi_Get_Data(int moln,double *phipsi);
	void PhiPsi_Put_Data(int moln,double *phipsi);
	//recon
	void PhiPsi_Recon_CB_From_NCaC(XYZ N,XYZ Ca,XYZ C,XYZ &Cb);
	void PhiPsi_Recon_CO_From_CaCN(XYZ Ca,XYZ C,XYZ N,XYZ &CO);
	//main
	void PhiPsi_To_BackCB(int moln,char *ami,double *phipsi,XYZ **mol); //phipsi->BackBone+CB
	void BackCB_To_PhiPsi(int moln,char *ami,double *phipsi,XYZ **mol); //BackBone+CB->phipsi
};

