#include "Mol_Out.h"

//--------------------- Start -----------------//
Mol_Out::Mol_Out(void)
{
	OutPrt=0;       // whether output printf    // (default:no)
}
Mol_Out::~Mol_Out(void)
{
}

//============================ Output_Related ============================//
int Mol_Out::Output_XYZ(FILE *fp,int len,XYZ *mol,char *AMI,char *CLE,char *ind,char Chain_ID)  // output normal XYZ file
{
	int i,j;
	char aminul='X';
	char clenul=' ';
	char amidum;
	char cledum;
	char dummychn;
	char out[7];

	//judge
	if(mol==NULL)return -1;

	//process
	for(i=0;i<len;i++)
	{
		if(AMI==NULL)amidum=aminul;
		else amidum=AMI[i];
		if(CLE==NULL)cledum=clenul;
		else cledum=CLE[i];
		if(Chain_ID==0)dummychn=' ';
		else dummychn=Chain_ID;
		
		if(ind==NULL)
		{
			fprintf(fp,":%c%4d  %c%c %8.3f%8.3f%8.3f\n",dummychn,i+1,amidum,cledum,mol[i].X,mol[i].Y,mol[i].Z);
		}
		else
		{
			for(j=0;j<6;j++)out[j]=ind[6*i+j];
			out[j]='\0';
			fprintf(fp,":%6s %c%c %8.3f%8.3f%8.3f\n",out,amidum,cledum,mol[i].X,mol[i].Y,mol[i].Z);
		}
	}
	return 1;
}
int Mol_Out::Output_PDB(FILE *fp,int len,XYZ *mol,char *AMI,char *ind,char Chain_ID) //output normal PDB file (CA only)
{
	int i,j;
	string buf;
	const char *dummyaa;
	char dummychn;
	char out[7];
	int atomnum;

	//judge
	if(mol==NULL)return -1;
	//process
	buf="                          ";
	atomnum=1;
	for(i=0;i<len;i++)
	{
		//init
		if(AMI==NULL)dummyaa=One2Three_III(0);
		else dummyaa=One2Three_III(AMI[i]);
		if(Chain_ID==0)dummychn=' ';
		else dummychn=Chain_ID;

		//process
		if(ind==NULL)
		{
			fprintf(fp,"ATOM  %5d  CA  %3s %c%4d    %8.3f%8.3f%8.3f%s\n",atomnum,
			dummyaa,dummychn,i+1,mol[i].X,mol[i].Y,mol[i].Z,buf.c_str());
		}
		else
		{
			for(j=0;j<6;j++)out[j]=ind[6*i+j];
			out[j]='\0';
			fprintf(fp,"ATOM  %5d  CA  %3s %6s   %8.3f%8.3f%8.3f%s\n",atomnum,
			dummyaa,out,mol[i].X,mol[i].Y,mol[i].Z,buf.c_str());
		}
		atomnum++;
	}
	//terminal
	return 1;
}
int Mol_Out::Output_PDB_I(FILE *fp,int len,XYZ *ca,XYZ *cb,char *AMI,char *ind,char Chain_ID) //output CA+CB file
{
	int i,j;
	string buf;
	const char *dummyaa;
	char dummychn;
	char out[7];
	int atomnum;

	//judge
	if(ca==NULL)return -1;
	if(cb==NULL)return -1;
	//process
	buf="                          ";
	atomnum=1;
	for(i=0;i<len;i++)
	{
		//init
		if(AMI==NULL)dummyaa=One2Three_III(0);
		else dummyaa=One2Three_III(AMI[i]);
		if(Chain_ID==0)dummychn=' ';
		else dummychn=Chain_ID;

		//process
		if(ind==NULL)
		{
			fprintf(fp,"ATOM  %5d  CA  %3s %c%4d    %8.3f%8.3f%8.3f%s\n",atomnum+0,
			dummyaa,dummychn,i+1,ca[i].X,ca[i].Y,ca[i].Z,buf.c_str());
			fprintf(fp,"ATOM  %5d  CB  %3s %c%4d    %8.3f%8.3f%8.3f%s\n",atomnum+1,
			dummyaa,dummychn,i+1,cb[i].X,cb[i].Y,cb[i].Z,buf.c_str());
		}
		else
		{
			for(j=0;j<6;j++)out[j]=ind[6*i+j];
			out[j]='\0';
			fprintf(fp,"ATOM  %5d  CA  %3s %6s   %8.3f%8.3f%8.3f%s\n",atomnum+0,
			dummyaa,out,ca[i].X,ca[i].Y,ca[i].Z,buf.c_str());
			fprintf(fp,"ATOM  %5d  CA  %3s %6s   %8.3f%8.3f%8.3f%s\n",atomnum+1,
			dummyaa,out,cb[i].X,cb[i].Y,cb[i].Z,buf.c_str());
		}
		atomnum+=2;
	}
	//terminal
	return 1;
}
int Mol_Out::Output_PDB_II(FILE *fp,int len,XYZ **mol,char *AMI,char *ind,char Chain_ID,int type) //output backbone+CB PDB file
{
	int i,j;
	string buf,ter;
	const char *dummyaa;
	char dummychn;
	char out[7];
	int atomnum;

	//judge
	if(mol==NULL)return -1;
	//process
	buf="                          ";
	atomnum=1;
	for(i=0;i<len;i++)
	{
		//init
		if(AMI==NULL)dummyaa=One2Three_III(0);
		else dummyaa=One2Three_III(AMI[i]);
		if(Chain_ID==0)dummychn=' ';
		else dummychn=Chain_ID;

		//process
		if(ind==NULL)
		{
			fprintf(fp,"ATOM  %5d  N   %3s %c%4d    %8.3f%8.3f%8.3f%s\n",atomnum+0,
			dummyaa,dummychn,i+1,mol[i][0].X,mol[i][0].Y,mol[i][0].Z,buf.c_str());
			fprintf(fp,"ATOM  %5d  CA  %3s %c%4d    %8.3f%8.3f%8.3f%s\n",atomnum+1,
			dummyaa,dummychn,i+1,mol[i][1].X,mol[i][1].Y,mol[i][1].Z,buf.c_str());
			fprintf(fp,"ATOM  %5d  C   %3s %c%4d    %8.3f%8.3f%8.3f%s\n",atomnum+2,
			dummyaa,dummychn,i+1,mol[i][2].X,mol[i][2].Y,mol[i][2].Z,buf.c_str());
			fprintf(fp,"ATOM  %5d  O   %3s %c%4d    %8.3f%8.3f%8.3f%s\n",atomnum+3,
			dummyaa,dummychn,i+1,mol[i][3].X,mol[i][3].Y,mol[i][3].Z,buf.c_str());
			if(type==1)
			{
				fprintf(fp,"ATOM  %5d  CB  %3s %c%4d    %8.3f%8.3f%8.3f%s\n",atomnum+4,
				dummyaa,dummychn,i+1,mol[i][4].X,mol[i][4].Y,mol[i][4].Z,buf.c_str());
			}
		}
		else
		{
			for(j=0;j<6;j++)out[j]=ind[6*i+j];
			out[j]='\0';
			fprintf(fp,"ATOM  %5d  N   %3s %6s   %8.3f%8.3f%8.3f%s\n",atomnum+0,
			dummyaa,out,mol[i][0].X,mol[i][0].Y,mol[i][0].Z,buf.c_str());
			fprintf(fp,"ATOM  %5d  CA  %3s %6s   %8.3f%8.3f%8.3f%s\n",atomnum+1,
			dummyaa,out,mol[i][1].X,mol[i][1].Y,mol[i][1].Z,buf.c_str());
			fprintf(fp,"ATOM  %5d  C   %3s %6s   %8.3f%8.3f%8.3f%s\n",atomnum+2,
			dummyaa,out,mol[i][2].X,mol[i][2].Y,mol[i][2].Z,buf.c_str());
			fprintf(fp,"ATOM  %5d  O   %3s %6s   %8.3f%8.3f%8.3f%s\n",atomnum+3,
			dummyaa,out,mol[i][3].X,mol[i][3].Y,mol[i][3].Z,buf.c_str());
			if(type==1)
			{
				fprintf(fp,"ATOM  %5d  CB  %3s %6s   %8.3f%8.3f%8.3f%s\n",atomnum+4,
				dummyaa,out,mol[i][4].X,mol[i][4].Y,mol[i][4].Z,buf.c_str());
			}
		}
		atomnum+=4;
		if(type==1)atomnum++;
	}
	//terminal
	return 1;
}
//---- neo PDB_Residue output ---//__110230__//
//OutType [0: sequential][1: PDB numbering]
//OutMode [-2:CA+CB,-1:CA,0:NCaC+CB,+1:ALL]
//OutGlys [0: GLY 'CB'='CA'][1: normal GLY]
int Mol_Out::Output_PDB_III(FILE *fp,int len,PDB_Residue *mol,char Chain_ID,
	int OutType,int OutMode,int OutGlys,int OutLast) //output full-atom PDB file
{
	int i,k;
	string TER="TER                                                                             ";
	string buf="              ";
	int number;
	char amino;
	const char *dummyaa;
	const char *atomname;
	string pdbind_;
	string pdbind;
	double x,y,z;
	double rfactor,temperature;
	int numb;
	int outtype_atom=1;
	char Ini_Chain;
	char Real_Chain;
	PDB_Residue PDB;
	XYZ xyz;
	int cur_pos;
	char cur_chain;
	char nxt_chain;
	//ws_new//__110430__//
	char output[100];

	//judge
	if(mol==NULL)return -1;
	Ini_Chain=Chain_ID;
	if(Ini_Chain=='!'||Ini_Chain==-1||Ini_Chain=='_')Ini_Chain='_'; //put the chain from pdbind//__110230__//
	//process
	cur_pos=0;
	for(i=0;i<len;i++)
	{
		//init
		PDB=mol[i];
		amino=PDB.get_AA();
		dummyaa=One2Three_III(amino);
		PDB.get_PDB_residue_number(pdbind_);
		pdbind=pdbind_.substr(1,5);
		if(Ini_Chain=='_')Real_Chain=pdbind_[0];
		else Real_Chain=Ini_Chain;
		cur_chain=pdbind_[0];
		//backbone_out
		number=PDB.get_backbone_totnum();  //this value should be 4
		for(k=0;k<number;k++)
		{
			//check
			if(OutMode<0) //only output CA or CB
			{
				if(k!=1)continue;
			}
			//get_backbone
			if(PDB.get_backbone_part_index(k)==0)
			{
				if(OutPrt==1)fprintf(stderr,"BackBone_Bad_Output[chain;%c][pos:%d][atom:%d]\n!!",Real_Chain,i+1,k+1);
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
				sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			else           // sequential numbering
			{
				if(OutType==0)  // partial sequential (atomic number sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
				else            // total sequential (atomic + residue sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,cur_pos+1,x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
			}
			//real_output//__110430__//
			if(OutLast==1) //-> at column 77, output 'ATOM_NAME' just as column 13
			{
				output[72]=cur_chain;
				output[77]=output[13];
			}
			fprintf(fp,"%s",output);
		}
		if(OutMode==-1)goto termi; //only output CA
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
				if(OutPrt==1)fprintf(stderr,"SideChain_Bad_Output[chain;%c][pos:%d][atom:%d]\n!!",Real_Chain,i+1,k+1);
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
				sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			else           // sequential numbering
			{
				if(OutType==0)  // partial sequential (atomic number sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
				else            // total sequential (atomic + residue sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,cur_pos+1,x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
			}
			//real_output//__110430__//
			if(OutLast==1) //-> at column 77, output 'ATOM_NAME' just as column 13
			{
				output[72]=cur_chain;
				output[77]=output[13];
			}
			fprintf(fp,"%s",output);
		}
termi:
		cur_pos++;
		//terminal_process
		if(Chain_ID=='!'||Chain_ID==-1)
		{
			if(i<len-1)
			{
				PDB_Residue pdb;
				string pdb_ind;
				pdb=mol[i+1];
				pdb.get_PDB_residue_number(pdb_ind);
				nxt_chain=pdb_ind[0];
				if(cur_chain!=nxt_chain)
				{
					fprintf(fp,"%s\n",TER.c_str());
					cur_pos=0;
				}
			}
		}
	}
	//terminal
	return 1;
}
int Mol_Out::Output_PDB_III(FILE *fp,int len,vector <PDB_Residue> &mol,char Chain_ID,
	int OutType,int OutMode,int OutGlys,int OutLast) //output full-atom PDB file
{
	int i,k;
	string TER="TER                                                                             ";
	string buf="              ";
	int number;
	char amino;
	const char *dummyaa;
	const char *atomname;
	string pdbind_;
	string pdbind;
	double x,y,z;
	double rfactor,temperature;
	int numb;
	int outtype_atom=1;
	char Ini_Chain;
	char Real_Chain;
	PDB_Residue PDB;
	XYZ xyz;
	int cur_pos;
	char cur_chain;
	char nxt_chain;
	//ws_new//__110430__//
	char output[100];

	//judge
	if(mol.size()==0)return -1;
	Ini_Chain=Chain_ID;
	if(Ini_Chain=='!'||Ini_Chain==-1||Ini_Chain=='_')Ini_Chain='_'; //put the chain from pdbind//__110230__//
	//process
	cur_pos=0;
	for(i=0;i<len;i++)
	{
		//init
		PDB=mol[i];
		amino=PDB.get_AA();
		dummyaa=One2Three_III(amino);
		PDB.get_PDB_residue_number(pdbind_);
		pdbind=pdbind_.substr(1,5);
		if(Ini_Chain=='_')Real_Chain=pdbind_[0];
		else Real_Chain=Ini_Chain;
		cur_chain=pdbind_[0];
		//backbone_out
		number=PDB.get_backbone_totnum();  //this value should be 4
		for(k=0;k<number;k++)
		{
			//check
			if(OutMode<0) //only output CA or CB
			{
				if(k!=1)continue;
			}
			//get_backbone
			if(PDB.get_backbone_part_index(k)==0)
			{
				if(OutPrt==1)fprintf(stderr,"BackBone_Bad_Output[chain;%c][pos:%d][atom:%d]\n!!",Real_Chain,i+1,k+1);
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
				sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			else           // sequential numbering
			{
				if(OutType==0)  // partial sequential (atomic number sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
				else            // total sequential (atomic + residue sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,cur_pos+1,x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
			}
			//real_output//__110430__//
			if(OutLast==1) //-> at column 77, output 'ATOM_NAME' just as column 13
			{
				output[72]=cur_chain;
				output[77]=output[13];
			}
			fprintf(fp,"%s",output);
		}
		if(OutMode==-1)goto termi; //only output CA
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
				if(OutPrt==1)fprintf(stderr,"SideChain_Bad_Output[chain;%c][pos:%d][atom:%d]\n!!",Real_Chain,i+1,k+1);
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
				sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			else           // sequential numbering
			{
				if(OutType==0)  // partial sequential (atomic number sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
				else            // total sequential (atomic + residue sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,cur_pos+1,x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
			}
			//real_output//__110430__//
			if(OutLast==1) //-> at column 77, output 'ATOM_NAME' just as column 13
			{
				output[72]=cur_chain;
				output[77]=output[13];
			}
			fprintf(fp,"%s",output);
		}
termi:
		cur_pos++;
		//terminal_process
		if(Chain_ID=='!'||Chain_ID==-1)
		{
			if(i<len-1)
			{
				PDB_Residue pdb;
				string pdb_ind;
				pdb=mol[i+1];
				pdb.get_PDB_residue_number(pdb_ind);
				nxt_chain=pdb_ind[0];
				if(cur_chain!=nxt_chain)
				{
					fprintf(fp,"%s\n",TER.c_str());
					cur_pos=0;
				}
			}
		}
	}
	//terminal
	return 1;
}
