#pragma once
#include <vector>
#include "Computation_Utility.h"
#include "XYZ.h"

//====class: Acc_Surface====//
//=> calculate solvent accessibility
class Acc_Surface
{
public:
	Acc_Surface(int num=PROT_MAX_NUM,int ORDER=2);
	~Acc_Surface(void);

//---- data_structure ---//
public:
	//[process]
	//inner atom data
	double AC_RN;       //atom_N radii   -> 1.65
	double AC_RCA;      //atom_CA radii  -> 1.87
	double AC_RC;       //atom_C radii   -> 1.76
	double AC_RO;       //atom_O radii   -> 1.40
	double AC_RSIDE;    //atom_sidechain -> 1.8
	double AC_WATER;    //atom_water     -> 1.4
	double AC_DATA[5];  //record the above five
	//inner icosahedron
	int icosahedron_order; //order (determine total number) -> default: 2 (4^2)
	int icosahedron_n;     //total number                   -> default: 320 (20*4*4)
	vector <XYZ> icosahedron_p;    //XYZ coordinate
	vector <double> icosahedron_a; //point area
	//acc calculate (res_level)
	XYZ *AC_mol_boxmin;    //boxmin around the residue
	XYZ *AC_mol_boxmax;    //boxmax around the residue
	int *AC_res_neibor_p;  //position of residue neibor
	int AC_res_neibor_n;   //total residue neibor
	//acc calculate (atom_level)
	XYZ AC_atom_center;    //atom center
	double AC_atom_radii;  //atom center radii
	int AC_neibor_n;       //total neibor
	vector <XYZ> AC_neibor_p;      //XYZ coordinate of neibor
	vector <double> AC_neibor_a;   //neibor area

	//[input/output]
	//solvent accessibility structure
	int AC_MAXIMAL;     //maximal length  -> default: PROT_MAX_NUM
	int AC_moln;        //mol length
	XYZ **AC_mol;       //input mol (N,CA,C,O) [backbone]
	int AC_side;        //total sidechain number
	XYZ *AC_sidechain;  //input sidechain molecular (omit sidechain name)
	int *AC_side_rec;   //input sidechain position (record each residue's starting point)
	int *AC_side_num;   //input sidechain number (record each residue's total sidecain)
	//solvent accessibility output
	int *AC_output;     //output acc data //-> Orig_ACC (absolute value)
	int *AC_normal;     //normal acc data //-> Normalized ACC ( betwee 0->100)
	double ACC_BuryCut; //nACC (10.25), for bury
	double ACC_ExpCut;  //nACC (42.9), for exposed

//---- main_function ---//
public:
	//[create]
	void AC_Init(int HB_MAXIMAL);
	void AC_Dele(int HB_MAXIMAL);
	//[init & input]
	void AC_Init_Data(void);  //assign inner atom data
	void AC_Triangle(XYZ x1, XYZ x2, XYZ x3,long level); //this function is recursive
	void AC_Init_Icosahedron(int order); //calculatet the icosahedron, given order
	void AC_Input_Mol(XYZ **mol,int moln,int *side_tot); //intput
	//[process]
	//residue level
	void AC_MinMax(XYZ v,double r,XYZ &vmin,XYZ &vmax);
	void AC_Calc_ResBox(void);
	int AC_Calc_ResNeib(int pos);
	//atomic level
	int AC_InBox(XYZ v,double r,XYZ vmin,XYZ vmax);
	void AC_Calc_AtomNeib_Single(XYZ v,double r,double dist);
	int AC_Calc_AtomNeib(XYZ v,double r);
	double AC_Calc_AtomAcc(void);
	//[main]
	void AC_Calc_SolvAcc(XYZ **mol,char *ami,int moln,
		char *acc,int *side_tot);
};
