#pragma once
#include "Confo_Lett.h"

//====class: Confo_Back====//
//=> CA<->Back generating
class Confo_Back : public Confo_Lett
{
public:
	Confo_Back(int num=PROT_MAX_NUM);
	~Confo_Back(void);
	int Confo_Back_Maximal;

//--- macros ---//
public:
	int Confo_Back_PRINTF;      // whether printf

//--- variable ---//
public:
	//normal dist+bend+angle
	double Confo_Back_CB_Ang;    // Ca(i)-> CB bend
	double Confo_Back_CB_Tor;    // Ca(i)-> CB tort
	double Confo_Back_CB_Dis;    // Ca(i)-> CB dist
	double Confo_Back_CN_Ang;    // Ca(i)-> N angle
	double Confo_Back_CC_Ang;    // Ca(i)-> C angle
	double Confo_Back_CO_Ang;    // Ca(i)-> O angle
	double Confo_Back_CN_Dis;    // Ca(i)-> N dist
	double Confo_Back_CC_Dis;    // Ca(i)-> C dist
	double Confo_Back_CO_Dis;    // Ca(i)-> O dist
	//temp_data
	double Confo_Back_ROT_TOT[3][3]; //3*3
	double Confo_Back_ROT_TMP[3][3]; //3*3
	double Confo_Back_ROT_MAT[3][3]; //3*3
	//additional class //__101125__//
	XYZ *temp_mol;   //pseudo angle to construct mol
	char *temp_cle;
	XYZ *temp_out;
	XYZ *case_mol;   //dealing with special case
	char *case_cle;
	XYZ **case_mcb;

//--- function ---//
public:
	//init
	void Confo_Back_Init(int maxlen);
	void Confo_Back_Dele(void);
	//process
	//[uni]
	void Construct_Mol(XYZ pre,XYZ cur,XYZ nxt,double bend,double tort,double dist,XYZ &out);
	void Construct_Mol_II(XYZ pre,XYZ cur,XYZ nxt,double bend,double tort,double dist,XYZ &out);
	void Construct_CB(XYZ N,XYZ Ca,XYZ C,double bend,double tort,double dist,XYZ &Cb);
	//[special]
	int Recon_Back_2pre(XYZ *mol,char *cle,int moln,XYZ **output,XYZ pre);        //-> deal with moln==2 (use pre)
	int Recon_Back_2nxt(XYZ *mol,char *cle,int moln,XYZ **output,XYZ nxt);        //-> deal with moln==2 (use nxt)
	int Recon_Back_Two(XYZ *mol,char *cle,int moln,XYZ **output,XYZ pre,XYZ nxt); //-> deal with moln==2 (use pre and nxt)
	int Recon_Back_1pre(XYZ *mol,char *cle,int moln,XYZ **output,XYZ pre);        //-> deal with moln==1 (use pre)
	int Recon_Back_1nxt(XYZ *mol,char *cle,int moln,XYZ **output,XYZ nxt);        //-> deal with moln==1 (use nxt)
	int Recon_Back_One(XYZ *mol,char *cle,int moln,XYZ **output,XYZ pre,XYZ nxt); //-> deal with moln==1 (use pre and nxt)
	//[main]
	int Recon_Back(XYZ *mol,char *cle,int moln,XYZ **output);  //input CA -> output Backbone (Method)
	int Recon_Back_Main(XYZ *in,char *cle,int moln,XYZ **out); //main procedure
};
