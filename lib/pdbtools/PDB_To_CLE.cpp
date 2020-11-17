#include "PDB_To_CLE.h"


//--------------constructor--------------------//
PDB_To_CLE::PDB_To_CLE(void)
{
	//PDB_macro//
	ORIAMI=0;    // whether align ori_AMI
	CaONLY=1;    // only record CA_ATOM
	CbBACK=1;    // consider CB_BACK (output PDB) 
}
PDB_To_CLE::~PDB_To_CLE(void)
{
}

//---------- function ---------//
void PDB_To_CLE::pdb_btb_ori(int head,int totnum,vector <XYZ> &r,vector <char> &CLE)  //__091210__//
{
	int i;
	XYZ *mol=new XYZ[totnum];
	char *cle=new char[totnum];
	for(i=0;i<totnum;i++)mol[i]=r.at(head+i);
	btb_ori(0,0,0,totnum,mol,cle);
	for(i=0;i<totnum;i++)CLE.at(head+i)=cle[i];
	delete [] mol;
	delete [] cle;
}
