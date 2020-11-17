#pragma once
#include "XYZ.h"
#include "PDB_Utility.h"
#include "Computation_Utility.h"

//====class: MultiAlign_Out====//
//=> output multiple alignment
class MultiAlign_Out
{
public:
	MultiAlign_Out(void);
	~MultiAlign_Out(void);

//--- data_structure ---//
public:
	//main
	int ***Multi_AFB;       //for script out
	int *Multi_AFB_Record;  //for script out
	int *Real_Block;        //for script out
	//temp
	int *BO_wsrec;
	int *BO_wsbak;
	int *BO_ZeroOne_Cur;

//--- process_function --//
public:
	//init & dele
	void BC_Output_Init(int num,int len);
	void BC_Output_Dele(int num,int len);
	//process
	int BC_Output_PDB(FILE *fp,int len,XYZ *mol,char *AMI,char *ind,char Chain_ID);
	void BC_Ali_To_AFB(int **ali,int totlen,int totnum,int ***AFB_Out,int *AFB_Len,int *Real_Core);
	//output
	void BC_Output_Alignment(FILE *fp,char **in,int totnum,int totlen,int **ali);
	void BC_Output_Superimpose(FILE *fp,XYZ **in,char **ami,int *len,int totnum);
	void BC_Output_RasMol_Script(FILE *fws,int TOT_NUM,int ***Multi_AFB,int *Multi_AFB_Record,int *Real_Block);
	void BC_Output_JMol_Script(FILE *fws,int TOT_NUM,int ***Multi_AFB,int *Multi_AFB_Record,int *Real_Block);
	//main
	//[f3 -> superposed structure rasmol]  -> TYPE=1 for JMol, TYPE=0 for RasMol
	void BC_Output_All(string &f3,XYZ **mol,char **ami,int *len,int totnum,int maxlen,int **ali,int TYPE=1);
};
