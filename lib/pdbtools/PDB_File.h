#pragma once 
#include "PDB_Chain.h"
#include "PDB_Residue.h" 
#include "PDB_Utility.h"
#include "Computation_Utility.h" 
#include <vector>
#include <string>
#include <map>

//====class: PDB_File====// 
//=> transfrom PDB_Format to Backbone 
class PDB_File 
{
public: 
	PDB_File(void);
	virtual ~PDB_File(void);
	PDB_File(const PDB_File & file);
	PDB_File & operator =(const PDB_File & file);

//====Macro for log & check==// 
public: 
	//file_out
	FILE *flog;     // ERR log 
	FILE *fpdb;     // PDB out 
	FILE *fxyz;     // XYZ out
	//macro_main
	int LOGOUT;     // whether output log 
	int PRTOUT;     // whether printf ERR 
	int ORIAMI;     // whether record ori_AMI 
	int CaONLY;     // only record CA_ATOM 
	int CbBACK;     // consider CB or ALL_ATOM //__090517__// 
	//macro_vice
	int OUTP_MODE;  // output type             //__090517__// [0: sequential][1: PDB numbering]
	int PROC_MODE;  // process mode            //__090517__// [-2:CA+CB,-1:CA,0:NCaC+CB,+1:ALL] 
	int GLYC_MODE;  // glycine mode            //__090517__// [0: GLY 'CB'='CA'][1: normal GLY]
	int HYDR_MODE;  // hydrogen mode           //__171007__// [1: output hydrogen][0: don't out] 
	//output
	int PDB_OUTPUT; // output PDB file
	int XYZ_OUTPUT; // output XYZ file
	//MODRES
	int MODRES;     // MODRES mapping

//----------- variables --------------//
public:   
	string PDB_root; 
	int PDB_num_rec;
	//data_structure
	vector <int> PDB_int_rec;
	vector <char> PDB_ins_rec;
	vector <char> PDB_tag_rec;
	vector <char> PDB_chn_rec;
	vector <char> PDB_ami_rec;
	vector <char> PDB_cle_rec;
	vector <XYZ> PDB_r_point;
	vector <PDB_Residue> PDB_output;  //-> record total PDB //neo//__090517__//
	//other_structure
	vector <int> TOTAL_CHAIN;  //__100604__//
	vector <char> NAME_CHAIN;  //__100604__//
	vector <string> PDB_record_all;
	//MODRES map
	map <string,string> PDB_MODRES_Map; //__130830__// -> for MODRES mapping, solve UNK problem


//-----the following is the declaration of functions 
public:
	// if chain_id is specified(A, B, C, ... ), coresponding chain is read
	// if chain_id == '_' return default chain
	// if chain_id == '!' return all chains
	int PDB_read_pdb_file(const string & fn,vector<PDB_Chain> &chains,char chain_id=-1);
	// OutMode [-2:CA+CB,-1:CA,0:NCaC+CB,+1:ALL] 
	void PDB_write_pdb_file_single(ostream &os,PDB_Chain &chain,
		int OutType=1,int OutMode=1,int OutGlys=1,int OutHydr=0);
	int PDB_write_pdb_file(ostream &os,vector<PDB_Chain> &chains,char chain_id=-1,
		int OutType=1,int OutMode=1,int OutGlys=1,int OutHydr=0);
	void PDB_write_pdb_file_single(FILE *fp,PDB_Chain &chain,
		int OutType=1,int OutMode=1,int OutGlys=1,int OutHydr=0);
	int PDB_write_pdb_file(FILE *fp,vector<PDB_Chain> &chains,char chain_id=-1,
		int OutType=1,int OutMode=1,int OutGlys=1,int OutHydr=0);
	// output PDB + XYZ
	void PDB_write_pdb_chain(FILE *fp,char *pdbid,char chain,int head,int totnum,
		vector <PDB_Residue> &PDB_output,int OutType=1,int OutMode=1,int OutGlys=1,int OutHydr=0);
	void PDB_write_xyz_chain(FILE *fp,char *pdbid,char chain,int head,int totnum,
		vector <int> &int_,vector <char> &ins_,vector <char> &tag_,vector <char> &chn_,
		vector <char> &ami_,vector <char> &cle_,vector <XYZ> &r_); 

//--------------------main_function--------------------// 
public: 
	//--- range process ----//

	int Input_XYZ_MINI_II(int PS,int mode,int st,char stt,int ed,char edd,char chain,
		int &totnu,XYZ *mol,char *AMI,char *CLE,char *ind,PDB_Residue *pdb);
	//--- For MODRES mapping -----//__130830__//
	int Process_MODRES_Mapping(const string &file,map<string,string> &out);
	void MODRES_Map(string &in,string &out,map <string,string > &ws_mapping);
	//--PART_II:PDB_Process Function
	//single
	int PreProcess_Record_Anis(vector <string> &PDB_record_all); 
	int PreProcess_Record_Alt(vector <string> &PDB_record_all); 
	int PreProcess_Record_Hydro(vector <string> &PDB_record_all,vector <string> &hydro_rec); 
	int PreProcess_Record_Check(vector <string> &PDB_record_all,char ori_c); 
	//total
	// GlyMode [0: glycine will output 'CB'='CA'][1: normal glycine output]
	int PDB_Process_Record(vector <string> &PDB_record_all,char ori_c,PDB_Residue &output); 
	//[check]
	int PDB_Residue_Check(PDB_Residue &output); 
	int PDB_File_Check(char *pdbid,char chain,int head,int totnum,
		vector <char> &PDB_tag_rec,double rate=0.5); 
	//--PART_III:Main        
	int PDB_To_XYZ_AMI_CLE(const string &fn,char ori_id=-1); // fundamental_part 

//--------- virtual_function ------------// 
public: 
	virtual void ori_ami_len_init(void); 
	virtual void pdb_btb_ori(int head,int totnum,vector <XYZ> &r,vector <char> &CLE); 
	virtual void record_ori_ami(char *input); 
	virtual int process_oriami_record(char *pdbid,char chain,int head,int totnum,
		vector <int> &int_,vector <char> &ins_,vector <char> &tag_,vector <char> &ami_); 
	virtual void output_oriami_record(FILE *fp,char *pdbid,char chain,int head,int totnum,
		vector <int> &int_,vector <char> &ins_,vector <char> &tag_,vector <char> &chn_,
		vector <char> &ami_,vector <char> &cle_,vector <XYZ> &r_); 
}; 
//====class: PDB_To_CLE====//over 
