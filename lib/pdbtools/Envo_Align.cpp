#include "Envo_Align.h"

//------ constructor -------//
Envo_Align::Envo_Align(int num,int CLESUM)
:TM_align(num),Bioinfo_Code(CLESUM)
{
	//init
	Envo_Align_maximal=num;
	Init_Envo_Align(Envo_Align_maximal);
	//macro
	Envo_Align_Do_TMalign=0;      // whether use TMalign to refine (default: no)
	//normal parameter
	Envo_Align_Distance=8.5;      // the distance cutoff (default: 8.5)
	Envo_Align_mol_extend=4;      // determine the initial rotmat (default: -l,+l)
	Envo_Align_neib_extend=0;     // determine the addtional neibor (default: -k,+k)
	//DynaProg parameter
	Envo_Align_AMI_GapOpen=-10;   // the gap open penalty for AMI
	Envo_Align_AMI_GapExtend=-1;  // the gap extend penalty for AMI
	Envo_Align_CLE_GapOpen=-100;  // the gap open penalty for CLE
	Envo_Align_CLE_GapExtend=-10; // the gap extend penalty for CLE
}
Envo_Align::~Envo_Align(void)
{
	Dele_Envo_Align();
}

//--------- init ---------//
void Envo_Align::Init_Envo_Align(int maxlen)
{
	Envo_Align_DP_sco=new double[maxlen*maxlen];
	//input
	Envo_Align_mol1=new XYZ[maxlen];
	Envo_Align_mol2=new XYZ[maxlen];
	Envo_Align_ami1=new int[maxlen];
	Envo_Align_ami2=new int[maxlen];
	Envo_Align_cle1=new int[maxlen];
	Envo_Align_cle2=new int[maxlen];
	//temp
	Envo_Align_tmp1=new XYZ[maxlen];
	Envo_Align_tmp2=new XYZ[maxlen];
	Envo_Align_alit=new int[maxlen];
	Envo_Align_ali2=new int[maxlen];
	Envo_Align_cle=new char[maxlen+1];
}
void Envo_Align::Dele_Envo_Align(void)
{
	delete [] Envo_Align_DP_sco;
	//input
	delete [] Envo_Align_mol1;
	delete [] Envo_Align_mol2;
	delete [] Envo_Align_ami1;
	delete [] Envo_Align_ami2;
	delete [] Envo_Align_cle1;
	delete [] Envo_Align_cle2;
	//temp
	delete [] Envo_Align_tmp1;
	delete [] Envo_Align_tmp2;
	delete [] Envo_Align_alit;
	delete [] Envo_Align_ali2;
	delete [] Envo_Align_cle;
}

//--------- vice ---------//
void Envo_Align::Envo_Align_Calc_Neib(XYZ *mol,int moln,vector <vector <int> > &Neib,double rcut)
{
	//init
	Neib.clear();
	Neib.resize(moln);
	//process
	int i,j;
	int k;
	double dij;
	double rcut2=rcut*rcut;
	for(i=0;i<moln;i++)
	{
		Neib[i].clear();
		for(j=0;j<moln;j++)Envo_Align_alit[j]=0;
		for(j=0;j<moln;j++)
		{
			dij=mol[i].distance_square(mol[j]);
			if(dij<rcut2)
			{
				Envo_Align_alit[j]=1;
				for(k=1;k<=Envo_Align_neib_extend;k++) //define additional neibors
				{
					if(j+k<moln)Envo_Align_alit[j+k]=1;
					if(j-k>=0)Envo_Align_alit[j-k]=1;
				}
			}
		}
		for(j=0;j<moln;j++)if(Envo_Align_alit[j]==1)Neib[i].push_back(j);
	}
}
void Envo_Align::Envo_Align_Change_Alignment(vector<pair<int,int> > &alignment,int size2,int *out_ali)
{
	int k;
	int wii,wjj;
	int totnum=(int)alignment.size();
	for(k=0;k<size2;k++)out_ali[k]=-1;
	for(k=0;k<totnum;k++)
	{
		wii=alignment[k].first;
		wjj=alignment[k].second;
		if(wii>0&&wjj>0)out_ali[wjj-1]=wii-1;
	}
}
double Envo_Align::Envo_Align_Refine_Alignment(int ii,int jj,int *out_ali)
{
	int i,j;
	int size1,size2;
	size1=(int)Envo_Align_Neib1[ii].size();
	size2=(int)Envo_Align_Neib2[jj].size();
	//collect
	for(i=0;i<size1;i++)Envo_Align_tmp1[i]=Envo_Align_mol1[Envo_Align_Neib1[ii][i]];
	for(j=0;j<size2;j++)Envo_Align_tmp2[j]=Envo_Align_mol2[Envo_Align_Neib2[jj][j]];
	//calc
	TM_Align_Init(size1,size2);
	double TM_CUR=Calc_TM_Align(Envo_Align_tmp1,Envo_Align_tmp2,size1,size2,out_ali,out_ali,
		-1,3,1); //run DynaProg
	//final
	return TM_CUR;
}

//--------- misc ---------//
double Envo_Align::Envo_Align_Calc_AMI(int ii,int jj,int *out_ali)
{
	int i,j;
	int pos1,pos2;
	int size1,size2;
	double score;
	size1=(int)Envo_Align_Neib1[ii].size();
	size2=(int)Envo_Align_Neib2[jj].size();
	//collect
	for(i=0;i<size1;i++)
	{
		int cur_index=i*size2;
		for(j=0;j<size2;j++)
		{
			pos1=Envo_Align_ami1[Envo_Align_Neib1[ii][i]];
			pos2=Envo_Align_ami2[Envo_Align_Neib2[jj][j]];
			score=Ori_BLOSUM[pos1][pos2];
			Envo_Align_DP_sco[cur_index+j]=score;
		}
	}
	//DynaProg
	double ali_sco;
	int ali_len;
	vector<pair<int,int> > alignment;
	ali_len=Advance_Align_Dyna_Prog(size1,size2,Envo_Align_DP_sco,Envo_Align_AMI_GapOpen,Envo_Align_AMI_GapExtend,
		alignment,ali_sco);
	//final
	Envo_Align_Change_Alignment(alignment,size2,out_ali);
	if(Envo_Align_Do_TMalign==0)return ali_sco;
	else return Envo_Align_Refine_Alignment(ii,jj,out_ali);
}
double Envo_Align::Envo_Align_Calc_CLE(int ii,int jj,int *out_ali)
{
	int i,j;
	int pos1,pos2;
	int size1,size2;
	double score;
	size1=(int)Envo_Align_Neib1[ii].size();
	size2=(int)Envo_Align_Neib2[jj].size();
	//collect
	for(i=0;i<size1;i++)
	{
		int cur_index=i*size2;
		for(j=0;j<size2;j++)
		{
			pos1=Envo_Align_cle1[Envo_Align_Neib1[ii][i]];
			pos2=Envo_Align_cle2[Envo_Align_Neib2[jj][j]];
			score=Gen_CLESUM[pos1][pos2];
			Envo_Align_DP_sco[cur_index+j]=score;
		}
	}
	//DynaProg
	double ali_sco;
	int ali_len;
	vector<pair<int,int> > alignment;
	ali_len=Advance_Align_Dyna_Prog(size1,size2,Envo_Align_DP_sco,Envo_Align_CLE_GapOpen,Envo_Align_CLE_GapExtend,
		alignment,ali_sco);
	//final
	Envo_Align_Change_Alignment(alignment,size2,out_ali);
	if(Envo_Align_Do_TMalign==0)return ali_sco;
	else return Envo_Align_Refine_Alignment(ii,jj,out_ali);
}
double Envo_Align::Envo_Align_Calc_mol(int ii,int jj,int *out_ali)
{
	int i,j;
	int size1,size2;
	size1=(int)Envo_Align_Neib1[ii].size();
	size2=(int)Envo_Align_Neib2[jj].size();
	//get initial
	int k;
	int count=0;
	Envo_Align_tmp1[count]=Envo_Align_mol1[ii];
	Envo_Align_tmp2[count]=Envo_Align_mol2[jj];
	count++;
	for(k=1;k<=Envo_Align_mol_extend;k++)
	{
		if(ii+k<Envo_Align_moln1 && jj+k<Envo_Align_moln2)
		{
			Envo_Align_tmp1[count]=Envo_Align_mol1[ii+k];
			Envo_Align_tmp2[count]=Envo_Align_mol2[jj+k];
			count++;
		}
		if(ii-k>=0 && jj-k>=0)
		{
			Envo_Align_tmp1[count]=Envo_Align_mol1[ii-k];
			Envo_Align_tmp2[count]=Envo_Align_mol2[jj-k];
			count++;
		}
	}
	//get rotmat
	double rotmat[12];
	kabsch(Envo_Align_tmp2,Envo_Align_tmp1,count,rotmat);
	//collect
	for(i=0;i<size1;i++)Envo_Align_tmp1[i]=Envo_Align_mol1[Envo_Align_Neib1[ii][i]];
	for(j=0;j<size2;j++)Envo_Align_tmp2[j]=Envo_Align_mol2[Envo_Align_Neib2[jj][j]];
	rot_mol(Envo_Align_tmp1,Envo_Align_tmp1,size1,rotmat);
	//TMalign
	TM_Align_Init(size1,size2);
	TM_Align_Get_Ali(Envo_Align_tmp1,Envo_Align_tmp2,size1,size2,Envo_Align_ali2);
	double TM_CUR=Calc_TM_Align(Envo_Align_tmp1,Envo_Align_tmp2,size1,size2,Envo_Align_ali2,Envo_Align_ali2,
		-1,3,1); //run DynaProg
	//final
	return TM_CUR;
}

//------------------ main_function --------------//
//input original structures
void Envo_Align::Input_Envo_Align(XYZ *mol1,XYZ *mol2,int moln1,int moln2,char *ami1,char *ami2)
{
	//molecular
	EqualArray(Envo_Align_mol1,mol1,moln1);
	EqualArray(Envo_Align_mol2,mol2,moln2);
	Envo_Align_moln1=moln1;
	Envo_Align_moln2=moln2;
	//amino_acid
	AMI_transform((const char*)ami1,Envo_Align_ami1);
	AMI_transform((const char*)ami2,Envo_Align_ami2);
	//confo_lett
	Confo_Lett confo_lett;
	confo_lett.btb_ori(0,0,0,Envo_Align_moln1,Envo_Align_mol1,Envo_Align_cle);
	Envo_Align_cle[Envo_Align_moln1]='\0';
	AMI_CLE_transform((const char*)ami1,(const char*)Envo_Align_cle,Envo_Align_cle1);
	confo_lett.btb_ori(0,0,0,Envo_Align_moln2,Envo_Align_mol2,Envo_Align_cle);
	Envo_Align_cle[Envo_Align_moln2]='\0';
	AMI_CLE_transform((const char*)ami2,(const char*)Envo_Align_cle,Envo_Align_cle2);
}
void Envo_Align::Input_Envo_Align_Total(XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		char *ami1,char *ami2,char *cle1,char *cle2)
{
	//molecular
	EqualArray(Envo_Align_mol1,mol1,moln1);
	EqualArray(Envo_Align_mol2,mol2,moln2);
	Envo_Align_moln1=moln1;
	Envo_Align_moln2=moln2;
	//amino_acid
	AMI_transform((const char*)ami1,Envo_Align_ami1);
	AMI_transform((const char*)ami2,Envo_Align_ami2);
	//confo_lett
	AMI_CLE_transform((const char*)ami1,(const char*)cle1,Envo_Align_cle1);
	AMI_CLE_transform((const char*)ami2,(const char*)cle2,Envo_Align_cle2);
}
//user may assign his own Neib_Definition for this input
void Envo_Align::Envo_Align_Neib_Main(XYZ *mcb1,XYZ *mcb2,int moln1,int moln2)
{
	//neib_calc
	Envo_Align_Calc_Neib(mcb1,moln1,Envo_Align_Neib1,Envo_Align_Distance);
	Envo_Align_Calc_Neib(mcb2,moln2,Envo_Align_Neib2,Envo_Align_Distance);
}
//Out_Type: [0] mol, [1] ami, [2] cle
void Envo_Align::Envo_Align_Calc_Main(vector <vector <double> > &out_mat,int Out_Type)
{
	//init
	int i,j;
	int moln1=Envo_Align_moln1;
	int moln2=Envo_Align_moln2;
	out_mat.resize(moln1);
	for(i=0;i<moln1;i++)out_mat[i].resize(moln2);
	//calc
	for(i=0;i<moln1;i++)
	{
		for(j=0;j<moln2;j++)
		{
			//calc
			if(Out_Type==0) //mol
			{
				out_mat[i][j]=Envo_Align_Calc_mol(i,j,Envo_Align_ali2);
			}
			if(Out_Type==1) //ami
			{
				out_mat[i][j]=Envo_Align_Calc_AMI(i,j,Envo_Align_ali2);
			}
			if(Out_Type==2) //cle
			{
				out_mat[i][j]=Envo_Align_Calc_CLE(i,j,Envo_Align_ali2);
			}
			//printf
//			printf("cur[%4d][%4d]\r",i,j);
		}
	}
}
