#include "Mol_File.h"


//--------------------- Start -----------------//
Mol_File::Mol_File(void)
{
	//init//
	ORIAMI=0;    // whether align ori_AMI
	CaONLY=1;    // only record CA_ATOM 
	CbBACK=1;    // consider CB_BACK (output PDB) 
}
Mol_File::~Mol_File(void)
{
}

//---------Upload-----------//Virtual_Reload//__080830__//================================//
int Mol_File::Upload_Process(string &file,char chain,char chain_bak,int PDBorSEQ,int tag,
	int head_tag,int head_c,int tail_tag,int tail_c,int &moln,XYZ *mol,char *AMI,char *CLE,char *ind,PDB_Residue *pdb)
{
	int moln_;
	int ret_val;
	if(PRE_LOAD==1)
	{
		moln_=0;
		ret_val=Input_XYZ_MINI_II(PDBorSEQ,tag,head_tag,head_c,tail_tag,tail_c,chain,moln_,0,0,0,0,0);
		if(ret_val!=1)return ret_val;      //means partially failed(upload mode)
		Input_XYZ_MINI_II(PDBorSEQ,tag,head_tag,head_c,tail_tag,tail_c,chain,moln,mol,AMI,CLE,ind,pdb);
		return 1;
	}
	
	
	//first try PDB format
	ret_val=PDB_To_XYZ_AMI_CLE(file,chain);           //try chain original
	if(ret_val<0)
	{
		ret_val=PDB_To_XYZ_AMI_CLE(file,chain_bak);   //try chain_bak
		if(ret_val>=1)chain=chain_bak;
	}
	if(ret_val<1)
	{
		return ret_val; //just return!! //__110230__//
		
		//second try XYZ format
		ret_val=Input_XYZ_FULL(file,chain,PDBorSEQ,tag,head_tag,head_c,tail_tag,tail_c,moln,mol,AMI,CLE,ind);
		if(ret_val!=1 && chain_bak!=chain)ret_val=Input_XYZ_FULL(file,chain_bak,PDBorSEQ,tag,head_tag,head_c,tail_tag,tail_c,moln,mol,AMI,CLE,ind);
		if(ret_val!=1)return ret_val;  //means totally failed(upload mode)
	}
	else
	{
		//memory limit
		if(PDB_num_rec>MEMORY_LIMIT)return -12345; //memory limit !!
		moln_=0;
		ret_val=Input_XYZ_MINI_II(PDBorSEQ,tag,head_tag,head_c,tail_tag,tail_c,chain,moln_,0,0,0,0,0);
		if(ret_val!=1)return ret_val;      //means partially failed(upload mode)
		Input_XYZ_MINI_II(PDBorSEQ,tag,head_tag,head_c,tail_tag,tail_c,chain,moln,mol,AMI,CLE,ind,pdb);
	}
	
	//final check
	if(ret_val!=1)return ret_val;      //means partially failed(upload mode)
	return 1;
}




