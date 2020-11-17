#include "PDB_File.h" 
#include "Utility.h"
#include <iostream>
using namespace std;

//--------------constructor--------------------// 
PDB_File::PDB_File(void) 
{ 
	//init// 
	PDB_root=".";// output XYZ dir 
	LOGOUT=0;    // whether output LOG 
	PRTOUT=0;    // whether printf ERR 
	ORIAMI=0;    // whether align ori_AMI 
	CaONLY=1;    // only record CA_ATOM (omit Non-Ca atoms!) 
	CbBACK=1;    // consider CB_BACK     //__090517__// 
	//macro
	OUTP_MODE=1; // [0: sequential][1: PDB numbering]
	PROC_MODE=1; // [-2:CA+CB,-1:CA,0:Back+CB,+1:ALL]
	GLYC_MODE=1; // [0: GLY 'CB'='CA'][1: normal GLY]
	HYDR_MODE=0; // [1: output hydrogen][0: don't out]  //__171007__//
	//output
	PDB_OUTPUT=0; // [0: don't output PDB][1: output ]
	XYZ_OUTPUT=0; // [0: don't output XYZ][1: output ]
	//MODRES map
	MODRES=0;     // [0: don't apply MODRES mapping]
} 
PDB_File::~PDB_File(void) 
{ 
}
PDB_File::PDB_File(const PDB_File & file)
{
	//[init]
	//init// 
	PDB_root=".";  // output XYZ dir 
	LOGOUT=0;    // whether output LOG 
	PRTOUT=0;    // whether printf ERR 
	ORIAMI=0;    // whether align ori_AMI 
	CaONLY=1;    // only record CA_ATOM (omit Non-Ca atoms!) 
	CbBACK=1;    // consider CB_BACK     //__090517__// 
	//macro
	GLYC_MODE=1; // [0: GLY 'CB'='CA'][1: normal GLY]
	OUTP_MODE=1; // [0: sequential][1: PDB numbering]
	PROC_MODE=1; // [-2:CA+CB,-1:CA,0:Back+CB,+1:ALL]
	//output
	PDB_OUTPUT=0; // [0: don't output PDB][1: output ]
	XYZ_OUTPUT=0; // [0: don't output XYZ][1: output ]
	//MODRES map
	MODRES=0;     // [0: don't apply MODRES mapping]
	
	//[copy]
	//data
	this->PDB_num_rec=file.PDB_num_rec;
	this->PDB_int_rec=file.PDB_int_rec;
	this->PDB_ins_rec=file.PDB_ins_rec;
	this->PDB_tag_rec=file.PDB_tag_rec;
	this->PDB_chn_rec=file.PDB_chn_rec;
	this->PDB_ami_rec=file.PDB_ami_rec;
	this->PDB_cle_rec=file.PDB_cle_rec;
	this->PDB_r_point=file.PDB_r_point;
	this->PDB_output=file.PDB_output;
	//chain
	this->TOTAL_CHAIN=file.TOTAL_CHAIN; //__100604__//
	this->NAME_CHAIN=file.NAME_CHAIN;   //__100604__//
}
PDB_File &  PDB_File::operator =(const PDB_File & file)
{
	if(this==&file)return *this;
	//data
	this->PDB_num_rec=file.PDB_num_rec;
	this->PDB_int_rec=file.PDB_int_rec;
	this->PDB_ins_rec=file.PDB_ins_rec;
	this->PDB_tag_rec=file.PDB_tag_rec;
	this->PDB_chn_rec=file.PDB_chn_rec;
	this->PDB_ami_rec=file.PDB_ami_rec;
	this->PDB_cle_rec=file.PDB_cle_rec;
	this->PDB_r_point=file.PDB_r_point;
	this->PDB_output=file.PDB_output;
	//chain
	this->TOTAL_CHAIN=file.TOTAL_CHAIN; //__100604__//
	this->NAME_CHAIN=file.NAME_CHAIN;   //__100604__//
	return *this;
}

//============================================= Input&Output_Part ========================================//
// after PDB->XYZ
int PDB_File::Input_XYZ_MINI_II(int PS,int mode,int st,char stt,int ed,char edd,char chain,
	int &totnu,XYZ *mol,char *AMI,char *CLE,char *ind,PDB_Residue *pdb)
{
	int i,ws_i;
	int chain_rec;
	int count;
	int num;
	int rel_num;
	int first,last;
	int ws_exclaim;
	string temp;
	int temp_len;
	char ami,cle;


	//--check--//
	if(totnu<0)return RANGEOVER_ERROR;
	count=totnu;

	//assign_head_tail
	if(stt=='@')first=1;
	else first=0;
	if(edd=='@')last=1;
	else last=0;

	if(mode==1) //head_plus_len
	{
		if(ed==-1)last=1;
		else last=0;
	}

	if(chain=='!')ws_exclaim=1;
	else ws_exclaim=0;

	//--[3]start
	num=0;
	rel_num=0;
	chain_rec=0; //default:FAIL
	if(chain=='_'||chain=='!')chain_rec=1; //chain_record
	for(i=0;i<PDB_num_rec;i++)
	{
		if(ws_exclaim==1)goto next;
		if(chain_rec==0)
		{
			if(chain==PDB_chn_rec.at(i))chain_rec=1; //chain_record
		}

		//[0]chain_check
		if(num==0)
		{
			if(chain!='_'&&chain!=PDB_chn_rec.at(i))continue;
		}
		else
		{
			if(chain!='_'&&chain!=PDB_chn_rec.at(i))break;
		}

		//[1]first_check
		if(first==0)
		{
			rel_num++;
			if(PS==1)  // PDB_numbering
			{
				if(PDB_int_rec.at(i)==st && PDB_ins_rec.at(i)==stt)first=1;
				else continue;
			}
			else       // SEQ_numbering
			{
				if(rel_num==st)
				{
					rel_num--;
					first=1;
				}
				else continue;
			}			
		}

next:
		//[2]record_index
		if(ind!=NULL)
		{
			temp=int2str(PDB_int_rec.at(i));
			temp_len=(int)temp.length();
			if(temp_len>4)return FILE_FORM_ERROR;
			
			for(ws_i=0;ws_i<6;ws_i++)ind[6*count+ws_i]=' ';
			ind[6*count+0]=PDB_chn_rec.at(i);
			for(ws_i=0;ws_i<temp_len;ws_i++)ind[6*count+4-ws_i]=temp[temp_len-ws_i-1];
			ind[6*count+5]=PDB_ins_rec.at(i);
		}


		//[3]input_STC,AMI.XYZ
		if(AMI!=NULL)
		{
			ami=PDB_ami_rec.at(i);
			if(ami<'A' || ami>'Z')return FILE_FORM_ERROR; //error
			AMI[count]=ami;
		}
		if(CLE!=NULL)
		{
			cle=PDB_cle_rec.at(i);
			if(cle<'A' || cle>'R')return FILE_FORM_ERROR; //error
			CLE[count]=cle;
		}
		if(mol!=NULL)mol[count]=PDB_r_point.at(i);
		if(pdb!=NULL)pdb[count]=PDB_output.at(i);


		//[4]step_gain
		count++;
		num++;
		rel_num++;
		if(ws_exclaim==1)continue;


		//[5]last_check
		if(last==0)
		{
			if(mode==0) // head_to_tail
			{
				if(PS==1)  // PDB_numbering
				{
					if(PDB_int_rec.at(i)==ed && PDB_ins_rec.at(i)==edd)
					{
						last=1;
						break;
					}
				}
				else       // SEQ_numbering
				{
					if(rel_num==ed)
					{
						last=1;
						break;
					}
				}
			}
			else       // head_plus_len
			{
				if(num==ed)
				{
					last=1;
					break;
				}
			}
		}
	}
	if(num==0)
	{
		if(chain_rec==0)return CHAINMISS_ERROR;
		else return RESIDOVER_ERROR;
	}

//	if(TERorNOT==1 && CLE!=NULL)
	if(CLE!=NULL)
	{
		CLE[totnu]='R';
		CLE[totnu+1]='R';
		CLE[count-1]='R';
	}
	if(CLE!=NULL)CLE[count]='\0';
	if(AMI!=NULL)AMI[count]='\0';

	totnu=count;
	if(last==0)return 2;
	return 1;
}

//============== CLEF_Neo ==============//__100604__//
//->read
int PDB_File::PDB_read_pdb_file(const string & fn,vector<PDB_Chain> & chains,char chain_id)
{
	//init
	int ret_val;
	ret_val=PDB_To_XYZ_AMI_CLE(fn,chain_id);
	if(ret_val<0)return ret_val;
	//evaluate
	int i,k;
	int length;
	char chainID;
	int count=0;
	chains.resize(ret_val);
	for(i=0;i<ret_val;i++)
	{
		length=TOTAL_CHAIN.at(i);
		chainID=NAME_CHAIN.at(i);
		chains.at(i).initialize_simple(length, chainID);
		for(k=0;k<length;k++)
		{
			chains.at(i).set_residue(k,PDB_output.at(count));
			count++;
		}
	}
	return 0;
}
//->write
//OutType [0: sequential][1: PDB numbering]
//OutMode [-2:CA+CB,-1:CA,0:NCaC+CB,+1:ALL]
//OutGlys [0: GLY 'CB'='CA'][1: normal GLY]
//OutHydr [1: output hydrogen][0: don't out]
//[OS mode]
void PDB_File::PDB_write_pdb_file_single(ostream &os,PDB_Chain &chain,
	int OutType,int OutMode,int OutGlys,int OutHydr)
{
	char fp[100];
	int j,k;
	int length;
	string buf="              ";
	int number;
	char amino;
	const char *dummyaa;
	const char *atomname;
	string pdbind;
	double x,y,z;
	double rfactor,temperature;
	int numb;
	int outtype_atom=1;
	char Real_Chain;
	PDB_Residue PDB;
	XYZ xyz;

	//process
	Real_Chain=chain.get_chain_id();
	if(Real_Chain=='!'||Real_Chain=='_'||Real_Chain==-1)Real_Chain='_';
	length=chain.get_length();
	for(j=0;j<length;j++)
	{
		//init
		chain.get_residue(j, PDB);
		amino=PDB.get_AA();
		dummyaa=One2Three_III(amino);
		PDB.get_PDB_residue_number(pdbind);
		if(Real_Chain=='_')Real_Chain=pdbind[0];
		//backbone_out
		for(k=0;k<4;k++)
		{
			//check
			if(OutMode<0) //only output CA or CB
			{
				if(k!=1)continue;
			}
			//get_backbone
			if(PDB.get_backbone_part_index(k)==0)
			{
				if(PRTOUT==1)fprintf(stderr,"BackBone_Bad_Output[chain;%c][pos:%d][atom:%d]\n!!",Real_Chain,j+1,k+1);
				continue;
			}
			atomname=backbone_atom_name_decode(k);
			PDB.get_backbone_atom(k,xyz, numb, rfactor, temperature);
			x=xyz.X;
			y=xyz.Y;
			z=xyz.Z;
			//output
			if(OutType==1) // PDB numbering
			{
				sprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			else           // sequential numbering
			{
				if(OutType==0)  // partial sequential (atomic number sequential)
				{
					sprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
				else            // total sequential (atomic + residue sequential)
				{
					sprintf(fp,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,j+1,x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
			}
			os<<fp;
		}
		if(OutMode==-1)continue; //only output CA
		//sidechain_out
		number=PDB.get_sidechain_totnum();
		for(k=0;k<number;k++)
		{
			//check
			if(OutMode<=0) //only output CB
			{
				if(k!=0)continue;
			}
			if(OutGlys==1) //glycine won't output CB
			{
				if(amino=='G')continue;
			}
			//get_sidechain
			if(PDB.get_sidechain_part_index(k)==0)
			{
				if(PRTOUT==1)fprintf(stderr,"SideChain_Bad_Output[chain;%c][pos:%d][atom:%d]\n!!",Real_Chain,j+1,k+1);
				continue;
			}
			atomname=sidechain_atom_name_decode(k,amino);
			PDB.get_sidechain_atom(k,xyz, numb, rfactor, temperature);
			x=xyz.X;
			y=xyz.Y;
			z=xyz.Z;
			//output
			if(OutType==1) // PDB numbering
			{
				sprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			else           // sequential numbering
			{
				if(OutType==0)  // partial sequential (atomic number sequential)
				{
					sprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
				else            // total sequential (atomic + residue sequential)
				{
					sprintf(fp,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,j+1,x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
			}
			os<<fp;
		}
		//OutHydr
		if(OutHydr==1)
		{
			vector <XYZ> hydro_out=PDB.hydro_rec;
			int hydr_num=(int)hydro_out.size()<999?(int)hydro_out.size():999;
			for(k=0;k<hydr_num;k++)
			{
				//get hydrogen name
				char hname[5];
				sprintf(hname,"H%-3d",k+1);
				hname[4]='\0';
				string hname_str=hname;
				//get coordinate
				xyz=hydro_out[k];
				x=xyz.X;
				y=xyz.Y;
				z=xyz.Z;
				//get others
				numb=0;
				rfactor=1;
				temperature=20;
				//output
				if(OutType==1) // PDB numbering
				{
					sprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						numb,hname_str.c_str(),dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
				}
				else           // sequential numbering
				{
					if(OutType==0)  // partial sequential (atomic number sequential)
					{
						sprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
							outtype_atom,hname_str.c_str(),dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
						outtype_atom++;
					}
					else            // total sequential (atomic + residue sequential)
					{
						sprintf(fp,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
							outtype_atom,hname_str.c_str(),dummyaa,Real_Chain,j+1,x,y,z,rfactor,temperature,buf.c_str());
						outtype_atom++;
					}
				}
				os<<fp;
			}
		}//end of 'if(OutHydr==1)'
	}
}
int PDB_File::PDB_write_pdb_file(ostream &os,vector<PDB_Chain> &chains,char chain_id,
	int OutType,int OutMode,int OutGlys,int OutHydr)
{
	char fp[100]; //PDB won't exceed 100
	int i;
	int totchain;
	PDB_Chain chain;
	string TER="TER                                                                             ";

	//process
	totchain=(int)chains.size();
	for(i=0;i<totchain;i++)
	{
		if(chain_id=='!'||chain_id==-1)goto begin;
		if(chain_id=='_')goto begin;
		if(chain_id==NAME_CHAIN.at(i))goto begin;
		else continue;
begin:
		chain=chains.at(i);
		PDB_write_pdb_file_single(os,chain,OutType,OutMode,OutGlys,OutHydr);
		sprintf(fp, "%s\n",TER.c_str());
		os<<fp;
		if(chain_id!='!'&&chain_id!=-1)break;
	}
	return 0;
}
//[FP mode]
void PDB_File::PDB_write_pdb_file_single(FILE *fp,PDB_Chain &chain,
	int OutType,int OutMode,int OutGlys,int OutHydr)
{
	int j,k;
	int length;
	string buf="              ";
	int number;
	char amino;
	const char *dummyaa;
	const char *atomname;
	string pdbind;
	double x,y,z;
	double rfactor,temperature;
	int numb;
	int outtype_atom=1;
	char Real_Chain;
	PDB_Residue PDB;
	XYZ xyz;

	//process
	Real_Chain=chain.get_chain_id();
	if(Real_Chain=='!'||Real_Chain=='_'||Real_Chain==-1)Real_Chain='_';
	length=chain.get_length();
	for(j=0;j<length;j++)
	{
		//init
		chain.get_residue(j, PDB);
		amino=PDB.get_AA();
		dummyaa=One2Three_III(amino);
		PDB.get_PDB_residue_number(pdbind);
		if(Real_Chain=='_')Real_Chain=pdbind[0];
		//backbone_out
		for(k=0;k<4;k++)
		{
			//check
			if(OutMode<0) //only output CA or CB
			{
				if(k!=1)continue;
			}
			//get_backbone
			if(PDB.get_backbone_part_index(k)==0)
			{
				if(PRTOUT==1)fprintf(stderr,"BackBone_Bad_Output[chain;%c][pos:%d][atom:%d]\n!!",Real_Chain,j+1,k+1);
				continue;
			}
			atomname=backbone_atom_name_decode(k);
			PDB.get_backbone_atom(k,xyz, numb, rfactor, temperature);
			x=xyz.X;
			y=xyz.Y;
			z=xyz.Z;
			//output
			if(OutType==1) // PDB numbering
			{
				fprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			else           // sequential numbering
			{
				if(OutType==0)  // partial sequential (atomic number sequential)
				{
					fprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
				else            // total sequential (atomic + residue sequential)
				{
					fprintf(fp,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,j+1,x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
			}
		}
		if(OutMode==-1)continue; //only output CA
		//sidechain_out
		number=PDB.get_sidechain_totnum();
		for(k=0;k<number;k++)
		{
			//check
			if(OutMode<=0) //only output CB
			{
				if(k!=0)continue;
			}
			if(OutGlys==1) //glycine won't output CB
			{
				if(amino=='G')continue;
			}
			//get_sidechain
			if(PDB.get_sidechain_part_index(k)==0)
			{
				if(PRTOUT==1)fprintf(stderr,"SideChain_Bad_Output[chain;%c][pos:%d][atom:%d]\n!!",Real_Chain,j+1,k+1);
				continue;
			}
			atomname=sidechain_atom_name_decode(k,amino);
			PDB.get_sidechain_atom(k,xyz, numb, rfactor, temperature);
			x=xyz.X;
			y=xyz.Y;
			z=xyz.Z;
			//output
			if(OutType==1) // PDB numbering
			{
				fprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			else           // sequential numbering
			{
				if(OutType==0)  // partial sequential (atomic number sequential)
				{
					fprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
				else            // total sequential (atomic + residue sequential)
				{
					fprintf(fp,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,j+1,x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
			}
		}
		//OutHydr
		if(OutHydr==1)
		{
			vector <XYZ> hydro_out=PDB.hydro_rec;
			int hydr_num=(int)hydro_out.size()<999?(int)hydro_out.size():999;
			for(k=0;k<hydr_num;k++)
			{
				//get hydrogen name
				char hname[5];
				sprintf(hname,"H%-3d",k+1);
				hname[4]='\0';
				string hname_str=hname;
				//get coordinate
				xyz=hydro_out[k];
				x=xyz.X;
				y=xyz.Y;
				z=xyz.Z;
				//get others
				numb=0;
				rfactor=1;
				temperature=20;
				//output
				if(OutType==1) // PDB numbering
				{
					fprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						numb,hname_str.c_str(),dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
				}
				else           // sequential numbering
				{
					if(OutType==0)  // partial sequential (atomic number sequential)
					{
						fprintf(fp,"ATOM  %5d %4s %3s %6s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
							outtype_atom,hname_str.c_str(),dummyaa,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
						outtype_atom++;
					}
					else            // total sequential (atomic + residue sequential)
					{
						fprintf(fp,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
							outtype_atom,hname_str.c_str(),dummyaa,Real_Chain,j+1,x,y,z,rfactor,temperature,buf.c_str());
						outtype_atom++;
					}
				}
			}
		}//end of 'if(OutHydr==1)'
	}
}
int PDB_File::PDB_write_pdb_file(FILE *fp,vector<PDB_Chain> &chains,
		char chain_id,int OutType,int OutMode,int OutGlys,int OutHydr)
{
	int i;
	int totchain;
	PDB_Chain chain;
	string TER="TER                                                                             ";

	//process
	totchain=(int)chains.size();
	for(i=0;i<totchain;i++)
	{
		if(chain_id=='!'||chain_id==-1)goto begin;
		if(chain_id=='_')goto begin;
		if(chain_id==NAME_CHAIN.at(i))goto begin;
		else continue;
begin:
		chain=chains.at(i);
		PDB_write_pdb_file_single(fp,chain,OutType,OutMode,OutGlys,OutHydr);
		fprintf(fp, "%s\n",TER.c_str());
		if(chain_id!='!'&&chain_id!=-1)break;
	}
	return 0;
}
//==================== normal PDB and XYZ output ==================//
// Output PDB File
void PDB_File::PDB_write_pdb_chain(FILE *fp,char *pdbid,char chain,int head,int totnum,
		vector <PDB_Residue> &PDB_output,int OutType,int OutMode,int OutGlys,int OutHydr)
{
	int k;
	PDB_Chain pdb_chain;
	pdb_chain.initialize_simple(totnum,chain);
	for(k=0;k<totnum;k++)pdb_chain.set_residue(k,PDB_output.at(head+k));
	PDB_write_pdb_file_single(fp,pdb_chain,OutType,OutMode,OutGlys,OutHydr);
}
// Output XYZ File
void PDB_File::PDB_write_xyz_chain(FILE *fp,char *pdbid,char chain,int head,int totnum,
		vector <int> &int_,vector <char> &ins_,vector <char> &tag_,vector <char> &chn_,
		vector <char> &ami_,vector <char> &cle_,vector <XYZ> &r_) 
{ 
	int i; 
	string out="";
	out=out+pdbid+chain;
	int num=-1; 
	fprintf(fp,">%5s %4d %4d\n",out.c_str(),totnum,num); 
	for(i=0;i<totnum;i++) 
	{ 
		fprintf(fp,"%c%c%4d%c %c%c %8.3f%8.3f%8.3f\n",abs(tag_.at(head+i)), 
		chn_.at(head+i),int_.at(head+i),ins_.at(head+i),ami_.at(head+i),cle_.at(head+i),
		r_.at(head+i).X,r_.at(head+i).Y,r_.at(head+i).Z); 
	} 
} 

//================================================PART_II====================================// PDB_Related 
// MODRES Process
int PDB_File::Process_MODRES_Mapping(const string &file,map <string,string > &ws_mapping)
{
	//--- list for mapping ---//
	map<string, string >::iterator iter;
	ws_mapping.clear();
	ifstream fin;
	string buf,temp,str1,str2;
	//read list1 for basic mapping
	fin.open(file.c_str(), ios::in);
	if(fin.fail()!=0)return -1;
	int len;
	char a;
	string key;
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<6)continue;
		temp=buf.substr(0,6);
		if(temp!="MODRES")continue;
		if(len<27)continue;
		str1=buf.substr(12,3);
		str2=buf.substr(24,3);
		//check invalid
		a=Three2One_III(str2.c_str());
		if(a=='X')continue;
		//mapping
		iter = ws_mapping.find(str1);
		if(iter != ws_mapping.end())
		{
			//check duplicated mapping
			key=ws_mapping[str1];
			if(key!=str2)if(PRTOUT==1)fprintf(stderr,"MODRES %s duplicated mapping to %s and %s \n",str1.c_str(),key.c_str(),str2.c_str());
			continue;
		}
		ws_mapping.insert(map < string, string >::value_type(str1, str2));
		count++;
	}
	return count;
}
//MODRES Map
void PDB_File::MODRES_Map(string &in,string &out,map <string,string > &ws_mapping)
{
	out=in;
	map<string, string >::iterator iter;
	iter = ws_mapping.find(in);
	if(iter != ws_mapping.end())out=ws_mapping[in];
}

// PDB_Process_Single
int PDB_File::PreProcess_Record_Anis(vector <string> &PDB_record_all)
{
	int i,k;
	int record_count=(int)PDB_record_all.size();
	char anis[5];
	int count=0;
	for(i=0;i<record_count;i++)
	{
		//kill tail
		for(k=0;k<3;k++)anis[k]=PDB_record_all[i][13+k];
		anis[k]='\0';
		if(strcmp(anis,"OXT")==0)continue;
		//kill header
		for(k=0;k<4;k++)anis[k]=PDB_record_all[i][k];
		anis[k]='\0';
		if(strcmp(anis,"ATOM")==0||strcmp(anis,"HETA")==0)
		{
			PDB_record_all[count]=PDB_record_all[i];
			count++;
		}
	}
	PDB_record_all.resize(count);
	return count;
}
int PDB_File::PreProcess_Record_Alt(vector <string> &PDB_record_all)
{
	int i,k;
	int record_count=(int)PDB_record_all.size();
	char name[5],back[5];
	int count=0;
	for(k=0;k<4;k++)back[k]=' ';
	back[k]='\0';
	for(i=0;i<record_count;i++)
	{
		for(k=0;k<4;k++)name[k]=PDB_record_all[i][12+k];
		name[k]='\0';
		if(strcmp(name,back)!=0)
		{
			PDB_record_all[count]=PDB_record_all[i];
			count++;
		}
		strcpy(back,name);
	}
	PDB_record_all.resize(count);
	return count;
}
int PDB_File::PreProcess_Record_Hydro(vector <string> &PDB_record_all, vector <string> &hydro_rec)
{
	int i;
	int record_count=(int)PDB_record_all.size();
	char h1,h2,hh,ht;
	int count=0;
	hydro_rec.clear();
	for(i=0;i<record_count;i++)
	{
		h1=PDB_record_all[i][12];
		h2=PDB_record_all[i][13];
		if(h1=='H'||h2=='H')hh='H';
		else hh=' ';
		ht=PDB_record_all[i][77];
		if(hh!='H' && ht!='H')
		{
			PDB_record_all[count]=PDB_record_all[i];
			count++;
		}
		else
		{
			hydro_rec.push_back(PDB_record_all[i]);
		}		
	}
	PDB_record_all.resize(count);
	return count;
}
int PDB_File::PreProcess_Record_Check(vector <string> &PDB_record_all,char ori_c)
{
	int i,k;
	int record_count=(int)PDB_record_all.size();
	char a;
	int correct=1;
	int count=0;
	string temp,reco;
	string temp_mod;
	char aachar[4];
	for(i=0;i<record_count;i++)
	{
		for(k=0;k<3;k++)aachar[k]=PDB_record_all[i][17+k];
		aachar[k]='\0';
		temp.assign(aachar);
		//MODRES MAP
		temp_mod=temp;
		if(MODRES==1)MODRES_Map(temp,temp_mod,PDB_MODRES_Map);
		//MODRES MAP over
		a=Three2One_III(temp_mod.c_str());
		if(a==ori_c)
		{
			PDB_record_all[count]=PDB_record_all[i];
			count++;
		}
		else correct=0;
	}
	PDB_record_all.resize(count);
	if(correct==1)return count;
	else return -1*count;
}
// PDB_Process_Full
// GlyMode [0: GLY 'CB'='CA'][1: normal GLY]
int PDB_File::PDB_Process_Record(vector <string> &PDB_record_all,char ori_c,PDB_Residue &output)
{
	int i,k;
	int record_count=(int)PDB_record_all.size();
	int count;
	char atom_serial_buf[6];
	int atom_serial_number;
	char atomchar[4];
	char amino;
	char posx[9];
	char posy[9];
	char posz[9];
	double posxf,posyf,poszf;
	double rfactor,temperature;
	string tempbuf;
	char pdbnum[7];
	int length;
	int back_ret,side_ret;

	//pre_process
	record_count=PreProcess_Record_Anis(PDB_record_all);        //kill "ANIS"
	record_count=PreProcess_Record_Alt(PDB_record_all);         //kill alt
	vector <string> hydro_rec;
	record_count=PreProcess_Record_Hydro(PDB_record_all,hydro_rec); //kill hydrogen
	record_count=PreProcess_Record_Check(PDB_record_all,ori_c); //kill different residue (different to ori_c)
	if(record_count<0)
	{
		if(PRTOUT==1)fprintf(stderr,"RES_MULTICOPY!!\r");
		record_count*=-1;
	}
	if(ori_c=='X')if(PRTOUT==1)fprintf(stderr,"WARNING!! -> RES_ABNORMAL!!\r");
	amino=ori_c;

	//---- process hydrogen -----//__2017_10_07__//
	if(HYDR_MODE==1)
	{
		vector <XYZ> hydro_out;
		for(i=0;i<(int)hydro_rec.size();i++)
		{
			tempbuf=hydro_rec[i];
			//========================= Get XYZ coordinate ===========================//
			for(k=0;k<8;k++)posx[k]=tempbuf[30+k];
			posx[k]='\0';    // get PDB_FILE x coordinate (8)
			for(k=0;k<8;k++)posy[k]=tempbuf[38+k];
			posy[k]='\0';    // get PDB_FILE y coordinate (8)
			for(k=0;k<8;k++)posz[k]=tempbuf[46+k];
			posz[k]='\0';    // get PDB_FILE z coordinate (8)
			//__080228__//apply Neo_Method , check the coordinate
			if(str2dou(posx,posxf)!=1)
			{
				if(PRTOUT==1)fprintf(stderr,"ERROR => HYDROGEN BAD X_POS AT [%d]\n",i);
				return 0; //coordinate error
			}
			if(str2dou(posy,posyf)!=1)
			{
				if(PRTOUT==1)fprintf(stderr,"ERROR => HYDROGEN BAD Y_POS AT [%d]\n",i);
				return 0; //coordinate error
			}
			if(str2dou(posz,poszf)!=1)
			{
				if(PRTOUT==1)fprintf(stderr,"ERROR => HYDROGEN BAD Z_POS AT [%d]\n",i);
				return 0; //coordinate error
			}
			//----> push_back
			XYZ xyz(posxf, posyf, poszf);
			hydro_out.push_back(xyz);
		}
		output.hydro_rec=hydro_out;
	}

	//init
	output.PDB_residue_backbone_initialize(amino); 
	count=AA26_sidechain_size(amino-'A'); 
	if(count!=record_count-4)
	{
		if(amino=='G'&&record_count==4)
		{
		}
		else
		{
			if(PRTOUT==1)fprintf(stderr,"RES_NOTEQUAL!!\r");
		}
	}
	output.PDB_residue_sidechain_initialize(amino); 
	//get PDB_number
	tempbuf=PDB_record_all[0];
	for(i=0;i<6;i++)pdbnum[i]=tempbuf[21+i];
	pdbnum[6]='\0';
	output.set_PDB_residue_number(string(pdbnum));
	//get total record
	for(i=0;i<record_count;i++)
	{
		tempbuf=PDB_record_all[i];
		//get_atom_serial
		for(k=0;k<5;k++)atom_serial_buf[k]=tempbuf[6+k];
		atom_serial_buf[k]='\0';
		if(str2int(atom_serial_buf,atom_serial_number)!=1)
		{
			if(PRTOUT==1)fprintf(stderr,"ERROR => ATOM SERIAL NUMBER AT [%d]\n",i);
			return 0; //coordinate error
		}
		//get_amino_acid
		for(k=0;k<3;k++)atomchar[k]=tempbuf[13+k];
		atomchar[k]='\0';
		//========================= Get XYZ coordinate ===========================//
		for(k=0;k<8;k++)posx[k]=tempbuf[30+k];
		posx[k]='\0';    // get PDB_FILE x coordinate (8)
		for(k=0;k<8;k++)posy[k]=tempbuf[38+k];
		posy[k]='\0';    // get PDB_FILE y coordinate (8)
		for(k=0;k<8;k++)posz[k]=tempbuf[46+k];
		posz[k]='\0';    // get PDB_FILE z coordinate (8)
		//__080228__//apply Neo_Method , check the coordinate
		if(str2dou(posx,posxf)!=1)
		{
			if(PRTOUT==1)fprintf(stderr,"ERROR => BAD X_POS AT [%d]\n",i);
			return 0; //coordinate error
		}
		if(str2dou(posy,posyf)!=1)
		{
			if(PRTOUT==1)fprintf(stderr,"ERROR => BAD Y_POS AT [%d]\n",i);
			return 0; //coordinate error
		}
		if(str2dou(posz,poszf)!=1)
		{
			if(PRTOUT==1)fprintf(stderr,"ERROR => BAD Z_POS AT [%d]\n",i);
			return 0; //coordinate error
		}
		//========================== Get Temprature =========================//__100601__//
		length=(int)tempbuf.length();
		if(length<66)  //missing
		{
			rfactor=1.0;
			temperature=20.0;
		}
		else
		{
			for(k=0;k<6;k++)posx[k]=tempbuf[54+k]; 
			posx[k]='\0';    // get PDB_FILE rfactor (6) 
			for(k=0;k<6;k++)posy[k]=tempbuf[60+k]; 
			posy[k]='\0';    // get PDB_FILE temperature (6) 
			//__080228__//apply Neo_Method , check the coordinate 
			if(str2dou(posx,rfactor)!=1) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"ERROR => BAD RFACTOR AT [%d]\r",i); 
				rfactor=1.0; //coordinate error 
			} 
			if(str2dou(posy,temperature)!=1) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"ERROR => BAD TEMPRATURE AT [%d]\r",i); 
				temperature=20.0; //coordinate error 
			} 
		}
		//backbone
		back_ret=backbone_atom_name_encode(atomchar); 
		side_ret=sidechain_atom_name_encode(atomchar,amino); 
		if((back_ret<0&&side_ret<0)||(back_ret>=0&&side_ret>=0)) //BAD!!
		{
			if(PRTOUT==1)fprintf(stderr,"WARNING!!! BAD_ATOM!!!\r");
		}
		else
		{
			//backbone
			if(back_ret>=0)
			{
				if(output.get_backbone_part_index(back_ret)==1)
				{
					if(PRTOUT==1)fprintf(stderr,"WARNING!! BACKBONE DOUBLE!!\r");
				}
				else
				{
					XYZ xyz(posxf, posyf, poszf);
					output.set_backbone_atom(back_ret,xyz,atom_serial_number, rfactor, temperature);
				}
			}
			//sidechain
			if(side_ret>=0)
			{
				if(output.get_sidechain_part_index(side_ret)==1)
				{
					if(PRTOUT==1)fprintf(stderr,"WARNING!! SIDECHAIN DOUBLE!!\r");
				}
				else
				{
					XYZ xyz(posxf, posyf, poszf);
					output.set_sidechain_atom(side_ret,xyz,atom_serial_number, rfactor, temperature);
				}
			}
		}
		//additional//__091210__//
		if(amino=='G') //glycine pseudo CB(=CA)
		{
			if(back_ret==1)
			{
				if(output.get_sidechain_part_index(0)==1)
				{
					if(PRTOUT==1)fprintf(stderr,"WARNING!! SIDECHAIN DOUBLE!!\r");
				}
				else
				{
					XYZ xyz(posxf, posyf, poszf);
					output.set_sidechain_atom(0,xyz,atom_serial_number, rfactor, temperature);
				}
			}
		}
	}
	return 1;
}
// PDB_Check
int PDB_File::PDB_Residue_Check(PDB_Residue &output) 
{ 
	int ret_val; 
	//CB_check 
	ret_val=output.PDB_residue_CB_check(); 
	if(ret_val!=1)return -10;  //CB error 
	//back_check 
	ret_val=output.PDB_residue_backbone_check(4); 
	if(ret_val!=1)return ret_val;  //backbone error 
	//all_check 
	ret_val=output.PDB_residue_sidechain_check(); 
	if(ret_val!=1)return 0;   //all-atom error 
	//term 
	return 1; 
}
int PDB_File::PDB_File_Check(char *pdbid,char chain,int head,int totnum,vector <char> &PDB_tag_rec,double rate) 
{ 
	int i; 
	int con_num=0; 
	for(i=0;i<totnum;i++)if(PDB_tag_rec.at(head+i)=='x')con_num++; 
	if(1.0*con_num/PDB_num_rec>rate) 
	{ 
		if(PRTOUT==1)fprintf(stderr,"KILLED => PDB_ID[%4s]:CHAIN[%c]\n",pdbid,chain); 
		return 0; 
	} 
	else return 1; 
} 


//====================PDB->CLE(print_out version)============//__071228__// 
int PDB_File::PDB_To_XYZ_AMI_CLE(const string &fn,char ori_id) // given chain 
{
	//--------MODRES process ---------//__130830__//
	PDB_MODRES_Map.clear();
	if(MODRES==1)
	{
		int retv=Process_MODRES_Mapping(fn,PDB_MODRES_Map);
		if(retv==-1)
		{
			if(PRTOUT==1)fprintf(stderr,"PDB File %s Not Found!!!\n",fn.c_str());
			return FILE_LOAD_ERROR;  // no such file;
		}
	}

 
	//---------claimation PART_I----------// I/O 
	ifstream fin; 
	string buf,temp; 
	string wwwtemp;
	char tempbuf[100];  //each line of PDB format won't exceed 100 !! //__110230__//
	int len; 
	//---------claimation PART_II----------// PDB_ID & PDB_Chain 
	int rescount=-1; 
	char pdbid[5]; 
	char ChianID; 
	char ChianID_C=-2; 
	char Process_ID=-2;           // __070806__// 
	//---------claimation PART_III----------// PDB_ATOM 
	char AtomName[5]; 
	char ResName[4];
	string ResName_ori;
	string ResName_mod;
	char ResSeq[5]; 
	char posx[9]; 
	char posy[9]; 
	char posz[9]; 
	//-------- temp record -------//
	int ResSeqI; 
	char iCodeC; 
	double posxf,posyf,poszf; 
	double posxt,posyt,poszt;   // __061106__// temp 
	char ResChar; 
	int ResError=-1; 
	int BrkError=-1; 
	int FirstTime;    // whether the first distance 
	double distance;        // distacne between two continuous res 
	//----- chain record -----//
	int OneChain=-1; //__070806__// 
	int ws_i;  //__070811__// read in buf 
	char tag[5]; 
	string temptag; 
	int FirstRes;  // whether the first CA_res 
	//__080305__// 
	int written_file; 
	//__080306__// 
	int Has_CA; 
	int CA_count; 
	double cax,cay,caz;      
	int ws_danger;   
	//__080307__// 
	int chain_index; 
	int chain_box[63]; 
	//test//__080307__// 
	int ws_ins_totnum=0; 
	int ws_deadend; 
	//neo//__080307__// 
	char ws_ins_bak; 
	int seq_record; 
	char Bak_INS=-1; 
	int Bak_SEQ=-1; 
	int Has_Done; 
	//__test__//__080310__// 
	int ws_DNA_test; 
	char wstag; 
	int wsout; 
	char ws_dead_chain; 
	//--- code last --//
	char code_last=-1; 
	int res_last; 
	int Found_AMI; 
	//neo_mod//__080416__// 
	int ws_exclaim; 
	int wscount; 
	int oo_n; 
	//--- ws's process all --//__090517__//(just_temp) 
	int should_record=0; 
	int ret_val; 
	string TER="TER                                                                             "; 
	//--- ws's process all --//__over__// 


	//===================== Ori_Process =====================//
	written_file=0;    // default[not written file]; 
	wsout=0;           // default[not out] 
	ws_dead_chain=-1;  // default[-1] 
	Found_AMI=0;       // default[not found] 
	ws_exclaim=0;      // default[no exclaim] 
	FirstRes=1;        //__080808__// 
	//--- vector init ---//__110230__//
	//[normal data]
	PDB_int_rec.clear();
	PDB_ins_rec.clear();
	PDB_tag_rec.clear();
	PDB_chn_rec.clear();
	PDB_ami_rec.clear();
	PDB_cle_rec.clear();
	PDB_r_point.clear();
	PDB_output.clear();
	PDB_num_rec=0;     //__070230__//
	//[temp data]
	TOTAL_CHAIN.clear();
	NAME_CHAIN.clear();
	PDB_record_all.clear();

	//======================================================================================//=>[Main Body] 
	/////////-------------- File_Input --------------//////////////=> [PART_0] 
	fin.open(fn.c_str(), ios::in); 
	if (fin.fail()!=0)  
	{ 
		if(PRTOUT==1)fprintf(stderr,"PDB File %s Not Found!!!\n",fn.c_str()); 
		fin.close(); 
		fin.clear(); 
		return FILE_LOAD_ERROR;  // no such file; 
	} 
	//neo//__080416__// 
	if(ori_id=='!') 
	{
		ori_id=-1; 
		ws_exclaim=1; 
	}
	//--- ori_id check ---//__110408__//
	if(ori_id==-1)  //only do when process All_Chains 
	{
		ws_exclaim=1;
		for(ws_i=0;ws_i<63;ws_i++)chain_box[ws_i]=0; //-> '63':  26+26+10=62, 'A-Z'+'a-z'+'0-9'
	}
	else            //assign chain 
	{
		if(ori_id=='_')ori_id=' '; 
		if(ori_id==' ')OneChain=1; //default 
		else OneChain=0; 
		Process_ID=ori_id;  
	}

	/////////-------------- Get_Head --------------////////////////=> [PART_1] 
headpart: 
	debug_info(1, "@PDB_File::PDB_To_XYZ_AMI_CLE. headpart\n");

	if(!getline(fin,buf,'\n'))return PREMATURE_ERROR; 
	getline_end(buf,0x0D); 
	len=(int)buf.length(); 
	if(len!=80) 
	{ 
		if(LOGOUT==1)fprintf(flog,"[L]"); 
		if(len<=4)return DATA_LENG_ERROR; //length <= 3 
		else if(len>80) 
		{ 
			temp=""; 
			temp=buf.substr(0,80); 
			strcpy(tempbuf,temp.c_str()); 
		} 
		else 
		{ 
			int wsi; 
			strcpy(tempbuf,buf.c_str()); 
			for(wsi=len;wsi<80;wsi++)tempbuf[wsi]=' '; 
			tempbuf[wsi]='\0'; 
		} 
	} 
	else strcpy(tempbuf,buf.c_str()); 

	//------------------Judge_Condition[2]------------------//=>break 
	for(ws_i=0;ws_i<4;ws_i++)tag[ws_i]=tempbuf[ws_i]; 
	tag[ws_i]='\0'; 
	temptag.assign(tag);
	if(temptag==("ATOM")||temptag==("HETA")) 
	{ 
		//__RealName__//__NEO__// 
		temp="1pdb"; 
		strcpy(pdbid,temp.c_str()); 
		goto wsmainpart; 
	} 
	else if(temptag==("TITL")||temptag==("REMA")||temptag==("SEQR")) 
	{ 
		//__RealName__//__NEO__// 
		temp="1pdb"; 
		strcpy(pdbid,temp.c_str()); 
		goto chainpart;  
	} 
	else if(temptag==("HEAD")) 
	{ 
		//__RealName__//__OLD__// 
		for(ws_i=0;ws_i<4;ws_i++)pdbid[ws_i]=tempbuf[62+ws_i]; 
		pdbid[ws_i]='\0'; 
		toLowerCase(pdbid); 
		goto chainpart;  
	} 
	else goto headpart;           //restart head part [neo_version!!] 

/////////-------------- Process Chain_Part --------------//////////=> [PART_2] 
chainpart: 
	debug_info(1, "@PDB_File::PDB_To_XYZ_AMI_CLE. chainpart\n");

	if(ORIAMI==1)ori_ami_len_init(); 
	for(;;) 
	{ 
		//------------------Judge_Condition------------------//=>break 
		for(ws_i=0;ws_i<4;ws_i++)tag[ws_i]=tempbuf[ws_i]; 
		tag[ws_i]='\0'; 
		temptag.assign(tag);
		if(ORIAMI==1 && temptag==("SEQR")) 
		{ 
			Found_AMI=1;  //found ori_ami 
			record_ori_ami(tempbuf); 
		} 
		if(temptag==("ATOM")||temptag==("HETA"))goto wsmainpart; 

nextline: 
		//read a new line 
		if(!getline(fin,buf, '\n'))return PREMATURE_ERROR; 
		getline_end(buf,0x0D); 
		len=(int)buf.length(); 
		if(len!=80) 
		{ 
			if(LOGOUT==1)fprintf(flog,"[L]"); 
			if(len<=4)goto nextline; //length <= 3 
			else if(len>80) 
			{ 
				temp=""; 
				temp=buf.substr(0,80); 
				strcpy(tempbuf,temp.c_str()); 
			} 
			else 
			{ 
				int wsi; 
				strcpy(tempbuf,buf.c_str()); 
				for(wsi=len;wsi<80;wsi++)tempbuf[wsi]=' '; 
				tempbuf[wsi]='\0'; 
			} 
		} 
		else strcpy(tempbuf,buf.c_str()); 
	}//end of FOR(::) 


/////////-------------- Process Main_Part --------------///////=> [PART_3] 
wsmainpart:      
	debug_info(1, "@PDB_File::PDB_To_XYZ_AMI_CLE. wsmainpart\n");

	for(;;) 
	{ 
		//-----------------------TOT_PART_I------------// get chain & create ouput_file 
		ChianID=tempbuf[21];     // get PDB_Chain 

		//======================================First Check==================================// 
		for(ws_i=0;ws_i<4;ws_i++)tag[ws_i]=tempbuf[ws_i]; 
		tag[ws_i]='\0';                  
		temptag.assign(tag);
		if(temptag==("CONE")||temptag==("END ")||temptag == ("MAST")||temptag==("ENDM")) 
		{ 
			if(written_file>0)goto wwout;  
			else 
			{ 
				if(PRTOUT==1)fprintf(stderr,"ERROR => FirstCheck EXIT AT PDB_ID[%4s]!!!\n",pdbid);
				if(ori_id!=-1&&ori_id!=' ')return CHAINMISS_ERROR;
				else return FILE_FORM_ERROR; //file_error
			} 
		} 
		if(ws_dead_chain==ChianID)goto wsgetline; 
		else ws_dead_chain=-1; 

		//----------------------------------Has Chain & All Chain-----------------------------// 
		if(ori_id!=-1)             //__[1]__has_chain//__071228__// 
		{ 
			if(OneChain==1) 
			{ 
				Process_ID=ChianID; 
				OneChain=0; 
			} 
			if(ChianID == Process_ID)goto wsbegin; 
		} 
		else                       //__[2]__all chain//__071228__// 
		{ 
			//-----------------------------Final_Check_Condition------------------// 
			if(temptag==("ATOM")||temptag==("HETA")) 
			{ 
				chain_index=CHAIN_to_INT62(ChianID); 
				if(chain_index!=-1 && chain_box[chain_index]==0) 
				{ 
					chain_box[chain_index]=1; 
					goto wsbegin; 
				} 
			} 
		} 

wsgetline: 
		//read a new line 
		if(!getline(fin,buf, '\n')) 
		{ 
			wsout=1; 
			goto wsexit; //premature break 
		} 
		getline_end(buf,0x0D); 
		len=(int)buf.length(); 
		if(len!=80) 
		{ 
			if(LOGOUT==1)fprintf(flog,"[L]"); 
			if(len<=4)goto wsgetline; //length <= 3 
			else if(len>80) 
			{ 
				temp=""; 
				temp=buf.substr(0,80); 
				strcpy(tempbuf,temp.c_str()); 
			} 
			else 
			{ 
				int wsi; 
				strcpy(tempbuf,buf.c_str()); 
				for(wsi=len;wsi<80;wsi++)tempbuf[wsi]=' '; 
				tempbuf[wsi]='\0'; 
			} 
		} 
		else strcpy(tempbuf,buf.c_str()); 
		//__070305__// 
		continue;    

    //======================================Real Begin==================================//[PART_4] 
wsbegin:
 
		//---init temp --// 
		posxt=posyt=poszt=0.0;   
		FirstTime=1; 
		FirstRes=1;  //__070811__// 
		if(ws_exclaim==0)
		{
			PDB_int_rec.clear();
			PDB_ins_rec.clear();
			PDB_tag_rec.clear();
			PDB_chn_rec.clear();
			PDB_ami_rec.clear();
			PDB_cle_rec.clear();
			PDB_r_point.clear();
			PDB_output.clear();
			PDB_num_rec=0;    //__070230__// 
		}
		wscount=0; 
		//__080306__//---process CA_only condition---// 
		Has_CA=1;    //default has CA 
		CA_count=0; 
		cax=cay=caz=0.0;                 
		ws_danger=0; //defualt NOT dangerous 
		//ori_record; 
		ws_ins_bak=-1; 
		seq_record=-99999; 
		//test//__080307__// 
		ws_ins_totnum=0; 
		ws_deadend=0;//default NOT deadend 
		//test//__080308__// 
		Has_Done=0;  //default NOT done 
		//first_record//__090517__// 
		should_record=0; 

		//-----------------------TOT_PART_II------------// get the CA_ATOM of the given_chain 
		for(;;) 
		{ 
			//========================chain_process==============//[TOT_PART_II][part_I] 
			for(ws_i=0;ws_i<4;ws_i++)tag[ws_i]=tempbuf[ws_i]; 
			tag[ws_i]='\0';                  
			temptag.assign(tag);
			if(temptag==("TER ")||temptag==("MAST"))  //terminal
			{ 
wsterandmast: 
				if(Has_CA==0 && CaONLY==0) 
				{ 
					ws_deadend=1; 
					goto wsdeadend; 
wsdeadendback: 
					ws_deadend=0; 
				} 

				//first_record//__090517__// 
				if(CbBACK==1) 
				{ 
					if(should_record==1) 
					{ 
						PDB_Residue pdb_residue;
						char pdb_ami=PDB_ami_rec.at(PDB_num_rec-1);
						ret_val=PDB_Process_Record(PDB_record_all,pdb_ami,pdb_residue); 
						if(ret_val!=1) 
						{ 
							if(PRTOUT==1)fprintf(stderr,"BackBone Error!!\n"); 
							return FILE_FORM_ERROR; //coordinate error 
						} 
						ret_val=PDB_Residue_Check(pdb_residue); 
						if(ret_val<0) 
						{ 
							if(PRTOUT==1)fprintf(stderr,"BackBone Missing!!\r"); 
							if(LOGOUT==1) 
							{ 
								if(ret_val==-10)fprintf(flog,"[CM%4d]",ResSeqI); 
								else fprintf(flog,"[%cM%4d]",-1*ret_val+'0',ResSeqI); 
							} 
						} 
						if(ret_val==0) 
						{ 
							if(PRTOUT==1)fprintf(stderr,"SideChain Missing!!\r"); 
							if(LOGOUT==1)fprintf(flog,"[SM%4d]",ResSeqI); 
						}
						PDB_output.push_back(pdb_residue);
					} 
					PDB_record_all.clear(); 
					should_record=0; 
				} 
				//first_record//__090517__//over 

				if(ori_id!=-1 && FirstRes==0)wsout=1; 
				if(ori_id!=-1 && FirstRes==1) 
				{ 
					OneChain=1; 
					goto wsgetline; 
				} 
				if(ori_id==-1 && FirstRes==1) 
				{ 
					chain_index=CHAIN_to_INT62(ChianID); 
					if(chain_index!=-1)chain_box[chain_index]=0; 
				} 
				goto wsexit; 
			} 
			else if(temptag==("CONE")||temptag==("END ")||temptag == ("MAST")||temptag==("ENDM")) 
			{ 
				wsout=1; 
				goto wsexit; 
			} 
			else if(temptag!=("ATOM")&&temptag!=("HETA"))  goto wsDNAtest; 


			//========================CA_process=================//[TOT_PART_II][part_II] 
			for(ws_i=0;ws_i<4;ws_i++)ResSeq[ws_i]=tempbuf[22+ws_i]; 
			ResSeq[ws_i]='\0';         // get PDB_FILE 7th column (4) 
			iCodeC=tempbuf[26];        // get PDB_FILE 8th column (1) 
			if(str2int(ResSeq,ResSeqI)!=1) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"ERROR => BAD ResSeq AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%4s]!\n",pdbid,ChianID,ResSeq); 
				return STR_TRANS_ERROR; //ResSeqI error 
			} 

			//test//__080307__// ResSeq_Minus TEST 
			if(ResSeqI<-99 || ResSeqI>9999) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"WARNING => ResSeqI ERROR AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%6d]!\n", 
				pdbid,ChianID,ResSeqI); 
			} 
			//test//__080307__// test the range of iCodeC 
			if(iCodeC!=' ' && (iCodeC<'A' || iCodeC>'Z')) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"WARNING => iCodeC ERROR AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%4d]!\n", 
				pdbid,ChianID,ResSeqI); 
			} 


			//__080306__// 
			if(ResSeqI!=seq_record || iCodeC!=ws_ins_bak) 
			{ 
wsdeadend: 
				if(Has_CA==0 && CaONLY==0) 
				{ 
					if(ChianID==' '&&FirstRes==0) 
					{ 
						wsout=1; 
						goto wsexit; 
					} 

					//__080307__//NEO KILL 
					temp.assign(ResName);
					if(temp[0]=='F'||temp==("ACE")||temp==("NH2")) 
					{ 
						if(PRTOUT==1)fprintf(stderr,"WARNING => ResName[XXX] AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%4d]!\n",pdbid,ChianID,seq_record); 
						goto wsfornext; 
					} 
					if(CA_count==0) 
					{ 
						if(PRTOUT==1)fprintf(stderr,"ERROR => FAINT CA AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%4d]!\n",pdbid,ChianID,ResSeqI); 
						return UNEXPECT_ERROR; //faint error                                             
					} 
					if(CA_count>20) 
					{ 
						if(PRTOUT==1)fprintf(stderr,"WARNING => TOO MANY CA[%3d] AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%4d]!\n", 
						CA_count,pdbid,ChianID,ResSeqI); 
					} 

					//swap record 
					Bak_SEQ=ResSeqI; 
					Bak_INS=iCodeC; 
					ResSeqI=seq_record; 
					iCodeC=ws_ins_bak; 

					//goto is DANGEROUS!!!
					ws_danger=1;  // Now is Dangerous!!                                      
					goto wsdanger; 
wsdangerback: 

					//========================= write into XYZ & AMI files ===========================// 
					posxf=cax/CA_count; 
					posyf=cay/CA_count; 
					poszf=caz/CA_count; 

					//__write to STC__// 
					{ 
						XYZ xyz(posxf,posyf,poszf);
						PDB_r_point.push_back(xyz);                                            
					} 
					//__write to STC__//over 

					//__write to char __//
					{
						wstag='x'; 
						PDB_int_rec.push_back(ResSeqI);
						PDB_ins_rec.push_back(iCodeC);
						PDB_tag_rec.push_back(wstag);
						PDB_chn_rec.push_back(ChianID);
						PDB_ami_rec.push_back(ResChar);
						PDB_cle_rec.push_back(' ');
					}
					//__write to char __//over

					//count++//
					PDB_num_rec++; 
					wscount++; 

					//========================= tail check ===========================//                                                             
					rescount++; 
					//__NEO_OUTPUT__//__080226__// 
					if(LOGOUT==1)fprintf(flog,"[X|%4d]",ResSeqI); 

					//restore record 
					ResSeqI=Bak_SEQ; 
					iCodeC=Bak_INS; 
					ws_danger=0;  // Now is NOT Dangerous!! 

					//first_record//__090517__// 
					should_record=1; 
				} 
wsfornext: 
				Has_CA=0; 
				CA_count=0; 
				cax=cay=caz=0.0; 
				Has_Done=0; // default NOT done 
				
				//record_last//__080322__// 
				code_last=ws_ins_bak; 
				res_last=seq_record; 
				
				//record the ResSeq & iCode 
				ws_ins_bak=iCodeC; 
				seq_record=ResSeqI; 

				//first_record//__090517__// 
				if(CbBACK==1) 
				{ 
					if(should_record==1) 
					{
						PDB_Residue pdb_residue;
						char pdb_ami=PDB_ami_rec.at(PDB_num_rec-1);
						ret_val=PDB_Process_Record(PDB_record_all,pdb_ami,pdb_residue); 
						if(ret_val!=1) 
						{ 
							if(PRTOUT==1)fprintf(stderr,"BackBone Error!!\n"); 
							return FILE_FORM_ERROR; //coordinate error 
						}
						ret_val=PDB_Residue_Check(pdb_residue); 
						if(ret_val<0) 
						{ 
							if(PRTOUT==1)fprintf(stderr,"BackBone Missing!!\r"); 
							if(LOGOUT==1) 
							{ 
								if(ret_val==-10)fprintf(flog,"[CM%4d]",ResSeqI); 
								else fprintf(flog,"[%cM%4d]",-1*ret_val+'0',ResSeqI); 
							} 
						} 
						if(ret_val==0) 
						{ 
							if(PRTOUT==1)fprintf(stderr,"SideChain Missing!!\r"); 
							if(LOGOUT==1)fprintf(flog,"[SM%4d]",ResSeqI); 
						}
						PDB_output.push_back(pdb_residue);
					} 
					PDB_record_all.clear();
					should_record=0; 
				} 
				//first_record//__090517__//over 
				if(ws_deadend==1)goto wsdeadendback; 
			} 


			//----- process AtomName ---------//
			debug_info(1, "tempbuf:%s\n", tempbuf);
			for(ws_i=0;ws_i<4;ws_i++)AtomName[ws_i]=tempbuf[12+ws_i]; 
			debug_info(1, "AtomName:%s\n", AtomName);
			AtomName[ws_i]='\0';       // get PDB_FILE 3th column (4) 
			temp.assign(AtomName);
			debug_info(1, "temp:%s\n", temp.c_str());

		//==================== process CA ====================//
		if(temp==(" CA ")) 
		{
			Has_CA=1;  // Has CA_ATOM 
			ResError = 1; // init ResError [default correct] 
			BrkError = 1; // init BrkError [default correct] 

			for(ws_i=0;ws_i<3;ws_i++)ResName[ws_i]=tempbuf[17+ws_i]; 
			ResName[ws_i]='\0';    // get PDB_FILE 5th column (3) 
			ChianID_C=tempbuf[21]; // get PDB_FILE 6th column (1) 

wsdanger: 
			//===First_Res===// 
			if(FirstRes==1) 
			{ 
				//__NEO_OUTPUT__//__080226__// 
				if(LOGOUT==1)
				{
					string wlog="";
					wlog=wlog + pdbid + ChianID;  
					fprintf(flog,"[%5s]=>",wlog.c_str()); 
				}

				//__080305__// 
				rescount=ResSeqI; 
				if(iCodeC!=' ') 
				{ 
					ws_ins_totnum++; 
					//__NEO_OUTPUT__//__080226__// 
					if(LOGOUT==1)fprintf(flog,"[%c_%4d]",iCodeC,ResSeqI); 
				} 
			} 
			else 
			{ 
				//__080304__//=> Alt_Error [Error_2] 
				if(Has_Done==1 && (ResSeqI == seq_record && iCodeC==ws_ins_bak)) 
				{ 
					//__061101__// => Chain error  [Error_5] 
					if(ChianID_C != ChianID)goto wsexit; 

					rescount = ResSeqI; 
					//__NEO_OUTPUT__//__080226__// 
					if(LOGOUT==1)fprintf(flog,"[A|%4d]",ResSeqI);    
                    
					if(ws_danger==1)goto wsdangerback; //__080306__//                                
					goto ws_altins; 
				} 

				//__071224__//=> Ins_Error [Error_3] 
				if(iCodeC!=' ') 
				{ 
					//test//__080307__// 
					ws_ins_totnum++; 

					if(iCodeC==code_last) 
					{ 
						if(ResSeqI != rescount)ResError = 0;  //ResError 
					} 
					else 
					{
						if(ResSeqI != rescount-1)ResError = 0;  //ResError 
						else 
						{ 
							if(code_last==' ') 
							{ 
								if(iCodeC!='A')ResError = 0; 
							} 
							else 
							{ 
								if(abs(iCodeC-code_last)!=1)ResError = 0;
							}
						} 
					} 
					rescount = ResSeqI; 
					//__NEO_OUTPUT__//__080226__// 
					if(LOGOUT==1)fprintf(flog,"[%c_%4d]",iCodeC,ResSeqI);
				}
				// __061106__// =>rescount error [Error_4] 
				else if(ResSeqI != rescount) 
				{
					if(ResSeqI != rescount-1) 
					{
						ResError = 0;    //ResError 
						//__NEO_OUTPUT__//__080226__// 
						if(LOGOUT==1)fprintf(flog,"[M|%4d]",ResSeqI);
					}
					rescount = ResSeqI; 
				}
			}


			//__070811__//=> ResChar_Wrong [Error_1] 
			ResName_ori.assign(ResName);
			//MODRES MAP
			ResName_mod=ResName_ori;
			if(MODRES==1)MODRES_Map(ResName_ori,ResName_mod,PDB_MODRES_Map);
			//MODRES MAP over
			ResChar = Three2One_III(ResName_mod.c_str()); 
			if(ResChar == 'X')             // WARNING =>ResChar error 
			{ 
				//__070811__// 
				if(temptag==("HETA"))
				{
					if(LOGOUT==1)fprintf(flog,"[H|%4d]",ResSeqI);
				}
				else
				{
					if(LOGOUT==1)fprintf(flog,"[W|%4d]",ResSeqI);
				}
			} 


			//__061101__// => Chain error  [Error_5] 
			if(ChianID_C != ChianID)goto wsexit; 


			//========================= tail check ===========================//                     
			if(FirstRes==1)FirstRes=0; 
			if(ws_danger==1) goto wsdangerback; 



			//========================= Get XYZ coordinate ===========================// 
			for(ws_i=0;ws_i<8;ws_i++)posx[ws_i]=tempbuf[30+ws_i]; 
			posx[ws_i]='\0';    // get PDB_FILE x coordinate (8) 
			for(ws_i=0;ws_i<8;ws_i++)posy[ws_i]=tempbuf[38+ws_i]; 
			posy[ws_i]='\0';    // get PDB_FILE y coordinate (8) 
			for(ws_i=0;ws_i<8;ws_i++)posz[ws_i]=tempbuf[46+ws_i]; 
			posz[ws_i]='\0';    // get PDB_FILE z coordinate (8) 


			//__080228__//apply Neo_Method , check the coordinate 
			if(str2dou(posx,posxf)!=1) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"ERROR => BAD X_POS AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%4d]!\n",pdbid,ChianID,ResSeqI); 
				return STR_TRANS_ERROR; //coordinate error 
			} 
			if(str2dou(posy,posyf)!=1) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"ERROR => BAD Y_POS AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%4d]!\n",pdbid,ChianID,ResSeqI); 
				return STR_TRANS_ERROR; //coordinate error 
			} 
			if(str2dou(posz,poszf)!=1) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"ERROR => BAD Z_POS AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%4d]!\n",pdbid,ChianID,ResSeqI); 
				return STR_TRANS_ERROR; //coordinate error 
			} 

			//========================= check whether break ===========================// 
			distance = sqrt( (posxf-posxt)*(posxf-posxt) + (posyf-posyt)*(posyf-posyt) + (poszf-poszt)*(poszf-poszt)); 
			if( ( distance < 2.5 || distance > 4.5) &&  FirstTime == 0) 
			{ 
				BrkError = 0;  //BrkError 
				//__NEO_OUTPUT__//__080226__// 
				if(LOGOUT==1)fprintf(flog,"[D|%4.1f]",distance); 
			}                                        
			posxt =posxf; 
			posyt =posyf; 
			poszt =poszf; 
			if(FirstTime == 1) FirstTime = 0;   

              
			//========================= write into XYZ & AMI files ===========================//     
			//__write to STC__// 
			{ 
				XYZ xyz(posxf,posyf,poszf);
				PDB_r_point.push_back(xyz);                                     
			} 
			//__write to STC__//over 

			//__write to char__//
			{
				if(ResError==0 && BrkError==0)wstag='+';   //Double Error [BrkEr+ResEr] 
				else 
				{ 
					if(ResError==0)wstag='|';              //[ResError] 
					else if(BrkError==0)wstag='-';         //[BrkError] 
					else wstag=' ';                        // Normal 
				}
				PDB_int_rec.push_back(ResSeqI);
				PDB_ins_rec.push_back(iCodeC);
				PDB_tag_rec.push_back(wstag);
				PDB_chn_rec.push_back(ChianID);
				PDB_ami_rec.push_back(ResChar);
				PDB_cle_rec.push_back(' ');
			}
			//__write to char__//over

			//count++//
			PDB_num_rec++; 
			wscount++; 

			//__080308__// 
			Has_Done=1; 

			//first_record//__090517__// 
			should_record=1; 
ws_altins:                               
			rescount++; 
		}// EMD OF IF(temp==(" CA ")) 
		else if(Has_CA==0 && CaONLY==0) 
		{ 
			if(CA_count==0) 
			{ 
				for(ws_i=0;ws_i<3;ws_i++)ResName[ws_i]=tempbuf[17+ws_i]; 
				ResName[ws_i]='\0';    // get PDB_FILE 5th column (3) 
				ChianID_C=tempbuf[21]; // get PDB_FILE 6th column (1) 
				iCodeC=tempbuf[26];    // get PDB_FILE 8th column (1) 
			} 

			//ws_DNA_test// 
			ws_DNA_test=0; 
			for(ws_i=0;ws_i<3;ws_i++)if(ResName[ws_i]==' ')ws_DNA_test++; 
			if(ws_DNA_test!=0) 
			{ 
				if(ChianID==' '&&FirstRes==0) 
				{ 
					wsout=1; 
					goto wsexit; 
				} 

				FirstRes=1; 
				ws_dead_chain=ChianID; 
				goto wsexit; 
			} 
			//water_test// 
			temp.assign(ResName);
			if(temp==("HOH")) 
			{ 
				if(ChianID==' '&&FirstRes==0) 
				{ 
					wsout=1; 
					goto wsexit; 
				} 

				FirstRes=1; 
				ws_dead_chain=ChianID; 
				goto wsexit; 
			} 
			//over// 

			//to dealing with non-CA residue, collect 99 atoms for average coordinate is enough //__110230__//
			if(CA_count>=99)goto wsDNAtest; 
			CA_count++; 

			//========================= Get XYZ coordinate ===========================// 
			for(ws_i=0;ws_i<8;ws_i++)posx[ws_i]=tempbuf[30+ws_i]; 
			posx[ws_i]='\0';    // get PDB_FILE x coordinate (8) 
			for(ws_i=0;ws_i<8;ws_i++)posy[ws_i]=tempbuf[38+ws_i]; 
			posy[ws_i]='\0';    // get PDB_FILE y coordinate (8) 
			for(ws_i=0;ws_i<8;ws_i++)posz[ws_i]=tempbuf[46+ws_i]; 
			posz[ws_i]='\0';    // get PDB_FILE z coordinate (8) 


			//__080228__//apply Neo_Method , check the coordinate 
			if(str2dou(posx,posxf)!=1) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"ERROR => BAD X_POS AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%4d]!\n",pdbid,ChianID,ResSeqI); 
				return STR_TRANS_ERROR; //coordinate error 
			} 
			if(str2dou(posy,posyf)!=1) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"ERROR => BAD Y_POS AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%4d]!\n",pdbid,ChianID,ResSeqI); 
				return STR_TRANS_ERROR; //coordinate error 
			} 
			if(str2dou(posz,poszf)!=1) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"ERROR => BAD Z_POS AT PDB_ID[%4s]:CHAIN[%c]:RES_NUM[%4d]!\n",pdbid,ChianID,ResSeqI); 
				return STR_TRANS_ERROR; //coordinate error 
			} 
			cax+=posxf; 
			cay+=posyf; 
			caz+=poszf; 
		} 

wsDNAtest: 
			//first_record//__090517__// 
			if(CbBACK==1) 
			{
				string ws_record;
				ws_record.assign(tempbuf);
				PDB_record_all.push_back(ws_record);
			} 
			//first_record//__090517__//over 

			//read a new line 
			if(!getline(fin,buf, '\n')) 
			{ 
				wsout=1; 
				goto wsexit; //premature break 
			} 
			getline_end(buf,0x0D); 
			len=(int)buf.length(); 
			if(len!=80) 
			{ 
				if(LOGOUT==1)fprintf(flog,"[L]"); 
				if(len<=3) 
				{ 
					if(ori_id!=-1 && FirstRes==0)wsout=1; 
					if(buf=="END")wsout=1; 
					goto wsterandmast; //length <= 3  //-> this may be "TER"
				} 
				else if(len>80) 
				{ 
					temp=""; 
					temp=buf.substr(0,80); 
					strcpy(tempbuf,temp.c_str()); 
				} 
				else 
				{ 
					int wsi; 
					strcpy(tempbuf,buf.c_str()); 
					for(wsi=len;wsi<80;wsi++)tempbuf[wsi]=' '; 
					tempbuf[wsi]='\0'; 
				} 
			} 
			else strcpy(tempbuf,buf.c_str()); 
		}// END OF for(;;) 

wsexit: 
		if(FirstRes==0) 
		{ 
			if(ws_exclaim==1)oo_n=PDB_num_rec-wscount; 
			else oo_n=0;

			//Debug//__110408__//
			if(wscount==0)
			{
				if(wsout==1)goto wwout;
				else goto check;
			}

			//generate CLE// 
			pdb_btb_ori(oo_n,wscount,PDB_r_point,PDB_cle_rec); 

			//first_record//__090517__// 
			if(CbBACK==1) 
			{ 
				if(should_record==1) 
				{ 
					PDB_Residue pdb_residue;
					char pdb_ami=PDB_ami_rec.at(PDB_num_rec-1);
					ret_val=PDB_Process_Record(PDB_record_all,pdb_ami,pdb_residue); 
					if(ret_val!=1) 
					{ 
						if(PRTOUT==1)fprintf(stderr,"BackBone Error!!\n"); 
						return FILE_FORM_ERROR; //coordinate error 
					} 
					ret_val=PDB_Residue_Check(pdb_residue); 
					if(ret_val<0) 
					{ 
						if(PRTOUT==1)fprintf(stderr,"BackBone Missing!!\r"); 
						if(LOGOUT==1) 
						{ 
							if(ret_val==-10)fprintf(flog,"[CM%4d]",ResSeqI); 
							else fprintf(flog,"[%cM%4d]",-1*ret_val+'0',ResSeqI); 
						} 
					} 
					if(ret_val==0) 
					{ 
						if(PRTOUT==1)fprintf(stderr,"SideChain Missing!!\r"); 
						if(LOGOUT==1)fprintf(flog,"[SM%4d]",ResSeqI); 
					}
					PDB_output.push_back(pdb_residue);
				} 
				PDB_record_all.clear();
				should_record=0; 
			} 
			//first_record//__090517__//over 

			//__NEO_OUTPUT__//__080226__// 
			if(LOGOUT==1)fprintf(flog,"\n"); 

			//output// 
			if(CaONLY==0)if(PDB_File_Check(pdbid,ChianID,oo_n,wscount,PDB_tag_rec)==0)goto check; 
			if(Found_AMI==1)if(process_oriami_record(pdbid,ChianID,oo_n,wscount,
				PDB_int_rec,PDB_ins_rec,PDB_tag_rec,PDB_ami_rec)!=1)goto check; 
//			if(PDBOUT==1) 
			{
				//init output file
				if(written_file==0) 
				{
					string outname;
					string out; 
					string out1; 
					outname=""; 
					outname=outname+pdbid; 
					if(ori_id==-1)out=PDB_root+'/'+outname; 
					else 
					{ 
						if(ori_id==' ')out=PDB_root+'/'+outname+'_'; 
						else out=PDB_root+'/'+outname+ori_id; 
					} 
					if(PDB_OUTPUT==1)
					{
						out1=out+".res";
						fpdb=fopen(out1.c_str(),"wb"); 
					}
					if(XYZ_OUTPUT==1)
					{
						out1=out+".xyz";
						fxyz=fopen(out1.c_str(),"wb"); 
					}
				}
				//output file
				if(written_file>=0)
				{
					if(PDB_OUTPUT==1) //output PDB
					{
						PDB_write_pdb_chain(fpdb,pdbid,ChianID,oo_n,wscount,PDB_output,
							OUTP_MODE,PROC_MODE,GLYC_MODE,HYDR_MODE);
						fprintf(fpdb,"%s\n",TER.c_str()); 
					}
					if(XYZ_OUTPUT==1) //output XYZ
					{ 
						if(Found_AMI==1)output_oriami_record(fxyz,pdbid,ChianID,oo_n,wscount,
							PDB_int_rec,PDB_ins_rec,PDB_tag_rec,PDB_chn_rec,PDB_ami_rec,PDB_cle_rec,PDB_r_point); 
						else PDB_write_xyz_chain(fxyz,pdbid,ChianID,oo_n,wscount,
							PDB_int_rec,PDB_ins_rec,PDB_tag_rec,PDB_chn_rec,PDB_ami_rec,PDB_cle_rec,PDB_r_point); 
					}
				}
			} 

			//__080305__// 
			TOTAL_CHAIN.push_back(wscount); //record chain length
			NAME_CHAIN.push_back(ChianID);  //record chain number
			written_file++;

			//---- init_temp ----//__110408__// 
			{
				FirstRes=1;  //__070811__// 
				wscount=0; 
			}

			//test//__080307__// 
			if(ws_ins_totnum>99) 
			{ 
				if(PRTOUT==1)fprintf(stderr,"WARNING => INS TOTNUM[%4d] EXCEED AT PDB_ID[%4s]:CHAIN[%c]\n", 
				ws_ins_totnum,pdbid,ChianID); 
			} 
		} 
		if(wsout==1)goto wwout; 
	            
check: 
		; 
	}// END OF for(;;) 

wwout: 

	if(written_file==0)return PREMATURE_ERROR; //premature break 
	if(written_file>0)
	{
		if(PDB_OUTPUT==1)fclose(fpdb);
		if(XYZ_OUTPUT==1)fclose(fxyz);
	}
	PDB_MODRES_Map.clear();
	return written_file;  // succeed 
} 


//---------------------- Virtual_Function ----------------// DynamicProgramming_Related 
void PDB_File::ori_ami_len_init(void) 
{ 
} 
void PDB_File::pdb_btb_ori(int head,int totnum,vector <XYZ> &r,vector <char> &CLE)  //__091210__// 
{ 
} 
void PDB_File::record_ori_ami(char *input) 
{ 
} 
int PDB_File::process_oriami_record(char *pdbid,char chain,int head,int totnum,
		vector <int> &int_,vector <char> &ins_,vector <char> &tag_,vector <char> &ami_) 
{ 
	return 1; 
} 
void PDB_File::output_oriami_record(FILE *fp,char *pdbid,char chain,int head,int totnum,
		vector <int> &int_,vector <char> &ins_,vector <char> &tag_,vector <char> &chn_,
		vector <char> &ami_,vector <char> &cle_,vector <XYZ> &r_) 
{ 
	PDB_write_xyz_chain(fp,pdbid,chain,head,totnum,int_,ins_,tag_,chn_,ami_,cle_,r_); 
}  
