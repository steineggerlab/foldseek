#include "CLEPAPS_Out.h"

//================= constructor ================//
CLEPAPS_Out::CLEPAPS_Out(int num)
{
	//vice parameter
	outcut1=1.5;
	outcut2=2.5;
	outmaxnum=10;
	//create
	CH_ali1=new int[num];
	CH_ali2=new int[num];
	CH_temp=new XYZ[num];
	CH_pdbt=new PDB_Residue[num];
	CH_AFP_Cor=new int[4*num];
}
CLEPAPS_Out::~CLEPAPS_Out(void)
{
	//delete
	delete [] CH_ali1;
	delete [] CH_ali2;
	delete [] CH_temp;
	delete [] CH_pdbt;
	delete [] CH_AFP_Cor;
}

//================== input functions ===============//
void CLEPAPS_Out::CLEPAPS_Input_Func(char *c1,char *c2,char *a1,char *a2,char *ind1,char *ind2,
									 PDB_Residue *p1,PDB_Residue *p2,XYZ *m1,XYZ *m2,int n1,int n2,
									 string &nam1,string &nam2)
{
	//copy pointer
	CH_CLE1=c1;
	CH_CLE2=c2;
	CH_AMI1=a1;
	CH_AMI2=a2;
	CH_IND1=ind1;
	CH_IND2=ind2;
	//copy data
	CH_pdb1=p1;
	CH_pdb2=p2;
	CH_mol1=m1;
	CH_mol2=m2;
	CH_moln1=n1;
	CH_moln2=n2;
	//name
	CH_NAM1=nam1;
	CH_NAM2=nam2;
}

//=================== output functions ============//
//[vice]
int CLEPAPS_Out::CLEPAPS_Out_RMSD(XYZ *m1,XYZ *m2,int len,double *rotmat_)
{
	int i;
	double score=0.0;
	int rel_sco;
	XYZ temp;
	if(len==0)return 0;
	for(i=0;i<len;i++)
	{
		CLEPAPS_Out_rot_point(m1[i],temp,rotmat_);
		score+=m2[i].distance_square(temp);
	}
	rel_sco=(int)(100.0*sqrt(1.0*score/len));
	return rel_sco;
}
int CLEPAPS_Out::CLEPAPS_Out_Script_Simp(Align_Record &align,int *AFP_Cor,int moln1,int moln2)
{
	int i;
	int num;
	int ii,jj;
	int count;
	int isFirst;
	int isLast;
	int head1,head2;
	int index;
	int cur_ii;

	//init
	num=0;
	AFP_Cor[0]=0;
	ii=-1;
	jj=-1;
	cur_ii=-1;
	head1=-1;
	head2=-1;
	isLast=0;
	isFirst=1;
	count=INT_MIN_NUM;
	for(i=0;i<moln2;i++)
	{
		cur_ii=align.alignment.at(i);
		if(cur_ii==-1) //purely blank
		{
			if(isFirst==0)
			{
				if(count>0)
				{
					AFP_Cor[0]++;
					index=AFP_Cor[0];
					AFP_Cor[4*index+0] = 0;
					AFP_Cor[4*index+1] = head1;
					AFP_Cor[4*index+2] = head2;
					AFP_Cor[4*index+3] = count;
					num+=count;
				}
				count=0;
				isFirst=1;
			}
			continue;
		}
		if(isFirst==1)
		{
ws_init:
			isFirst=0;
			ii=cur_ii;
			jj=i;
			count=1;
			head1=ii;
			head2=jj;
#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif
			continue;
		}
		if(i==jj+1&&cur_ii==ii+1)
		{
			ii=cur_ii;
			jj=i;
			count++;
			continue;
		}

ws_end:
		if(count>0)
		{
			AFP_Cor[0]++;
			index=AFP_Cor[0];
			AFP_Cor[4*index+0] = 0;
			AFP_Cor[4*index+1] = head1;
			AFP_Cor[4*index+2] = head2;
			AFP_Cor[4*index+3] = count;
			num+=count;
		}
		if(isLast==1)goto end;
		else goto ws_init;
	}
	if(count==INT_MIN_NUM)goto end;
	isLast=1;
	goto ws_end;
end:
	return num;
}
int CLEPAPS_Out::CLEPAPS_Out_Script_Ori(Align_Record &align,int *AFP_Cor,
										XYZ *mol1,XYZ *mol2,int moln1,int moln2,double *rotmat_)
{
	int i;
	int num;
	int ii,jj;
	int count;
	int isFirst;
	int isLast;
	int head1,head2;
	int index;
	int cur_ii;

	//init
	num=0;
	AFP_Cor[0]=0;
	ii=-1;
	jj=-1;
	cur_ii=-1;
	head1=-1;
	head2=-1;
	isLast=0;
	isFirst=1;
	count=INT_MIN_NUM;
	for(i=0;i<moln2;i++)
	{
		cur_ii=align.alignment.at(i);
		if(cur_ii==-1) //purely blank
		{
			if(isFirst==0)
			{
				if(count>0)
				{
					AFP_Cor[0]++;
					index=AFP_Cor[0];
					AFP_Cor[4*index+0] = CLEPAPS_Out_RMSD(mol1+head1,mol2+head2,count,rotmat_);
					AFP_Cor[4*index+1] = head1;
					AFP_Cor[4*index+2] = head2;
					AFP_Cor[4*index+3] = count;
					num+=count;
				}
				count=0;
				isFirst=1;
			}
			continue;
		}
		if(isFirst==1)
		{
ws_init:
			isFirst=0;
			ii=cur_ii;
			jj=i;
			count=1;
			head1=ii;
			head2=jj;
#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif
			continue;
		}
		if(i==jj+1&&cur_ii==ii+1)
		{
			ii=cur_ii;
			jj=i;
			count++;
			continue;
		}

ws_end:
		if(count>0)
		{
			AFP_Cor[0]++;
			index=AFP_Cor[0];
			AFP_Cor[4*index+0] = CLEPAPS_Out_RMSD(mol1+head1,mol2+head2,count,rotmat_);
			AFP_Cor[4*index+1] = head1;
			AFP_Cor[4*index+2] = head2;
			AFP_Cor[4*index+3] = count;
			num+=count;
		}
		if(isLast==1)goto end;
		else goto ws_init;
	}
	if(count==INT_MIN_NUM)goto end;
	isLast=1;
	goto ws_end;
end:
	return num;
}
void CLEPAPS_Out::CLEPAPS_Out_Ali2_To_Ali1(int moln1,int moln2,int *ali1,int *ali2)
{
	int i;
	int ii,jj;
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]>=0)
		{
			ii=ali2[i];
			jj=i;
#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif	
			ali1[ii]=jj;
		}
	}
}
void CLEPAPS_Out::CLEPAPS_Out_rot_point(XYZ p_old,XYZ &p_new,double *rotmat_)
{
	double dx,dy,dz;
	dx=p_old.X;
	dy=p_old.Y;
	dz=p_old.Z;
	p_new.X=dx*rotmat_[0]+dy*rotmat_[1]+dz*rotmat_[2]+rotmat_[9];
	p_new.Y=dx*rotmat_[3]+dy*rotmat_[4]+dz*rotmat_[5]+rotmat_[10];
	p_new.Z=dx*rotmat_[6]+dy*rotmat_[7]+dz*rotmat_[8]+rotmat_[11];
}
void CLEPAPS_Out::CLEPAPS_Out_rot_mol(XYZ *m_old,XYZ *m_new,int nAtom,double *rotmat_)
{
	int l;
	double dx,dy,dz;
	for(l=0;l<nAtom;l++) 
	{
		dx=m_old[l].X; 
		dy=m_old[l].Y;
		dz=m_old[l].Z;
		m_new[l].X=dx*rotmat_[0]+dy*rotmat_[1]+dz*rotmat_[2]+rotmat_[9];
		m_new[l].Y=dx*rotmat_[3]+dy*rotmat_[4]+dz*rotmat_[5]+rotmat_[10];
		m_new[l].Z=dx*rotmat_[6]+dy*rotmat_[7]+dz*rotmat_[8]+rotmat_[11];
	}
}
void CLEPAPS_Out::CLEPAPS_Out_rot_pdb(PDB_Residue *m_old,PDB_Residue *m_new,int nAtom,double *rotmat_)
{
	int i,k;
	int num;
	PDB_Residue pdb;
	XYZ xyz;
	for(i=0;i<nAtom;i++)
	{
		pdb=m_old[i];
		//rotate backbone
		num=pdb.get_backbone_totnum();
		for(k=0;k<num;k++)
		{
			if(pdb.get_backbone_part_index(k)==0)continue;
			pdb.get_backbone_atom(k,xyz);
			CLEPAPS_Out_rot_point(xyz,xyz,rotmat_);
			pdb.set_backbone_atom(k,xyz);
		}
		//rotate sidechain
		num=pdb.get_sidechain_totnum();
		for(k=0;k<num;k++)
		{
			if(pdb.get_sidechain_part_index(k)==0)continue;
			pdb.get_sidechain_atom(k,xyz);
			CLEPAPS_Out_rot_point(xyz,xyz,rotmat_);
			pdb.set_sidechain_atom(k,xyz);
		}
		m_new[i]=pdb;
	}
}
//[function]
void CLEPAPS_Out::CLEPAPS_Output_Script(FILE *fws,string &name,int *AFP_Cor)
{
	//process
	int j,k;
	int jj;
	int totnum;
	int winlen;
	fprintf(fws,"load ./%s\n",name.c_str());
	fprintf(fws,"wireframe off\n");
	fprintf(fws,"backbone 10\n");
	fprintf(fws,"set ambient 20\n");
	fprintf(fws,"set background white\n");
	totnum=AFP_Cor[0];
	for(k=1;k<=totnum;k++)
	{
		for(j=0;j<2;j++)
		{
			//record
			jj=AFP_Cor[k*4+j+1];
			if(jj==-1)continue;
			winlen=AFP_Cor[k*4+3];
			fprintf(fws,"select   %3d-%3d:%c\n", jj+1,jj+winlen,j+'A');
			fprintf(fws,"backbone 100\n");
		}
	}
	fprintf(fws,"select all\n");
}
void CLEPAPS_Out::CLEPAPS_Output_MolCA(FILE *fp,Align_Record &align)
{
	int i;
	double rotmat[12];
	string TER="TER                                                                             ";
	for(i=0;i<12;i++)rotmat[i]=align.rotmat[i];
	CLEPAPS_Out_rot_mol(CH_mol1,CH_temp,CH_moln1,rotmat);
	Output_PDB(fp,CH_moln1,CH_temp,CH_AMI1,0,'A');
	fprintf(fp,"%s\n",TER.c_str());
	Output_PDB(fp,CH_moln2,CH_mol2,CH_AMI2,0,'B');
	fprintf(fp,"%s\n",TER.c_str());
}
void CLEPAPS_Out::CLEPAPS_Output_Align(FILE *fp,Align_Record &align,int FullOrNot)
{
	//init eval
	int moln1=CH_moln1;
	int moln2=CH_moln2;
	XYZ *mol1=CH_mol1;
	XYZ *mol2=CH_mol2;
	XYZ *AFP_tmp1=CH_temp;
	int *ali1=CH_ali1;
	int *ali2=CH_ali2;
	int *AFP_Cor=CH_AFP_Cor;
	double CLEP_rotmat[12];
	//init done
	int i;
	int score;
	int ii,jj,winlen;
	for(i=0;i<moln2;i++)ali2[i]=align.alignment[i];
	for(i=0;i<12;i++)CLEP_rotmat[i]=align.rotmat[i];
	CLEPAPS_Out_Ali2_To_Ali1(moln1,moln2,ali1,ali2);
	CLEPAPS_Out_Script_Ori(align,AFP_Cor,mol1,mol2,moln1,moln2,CLEP_rotmat);
	int totnum=AFP_Cor[0];
	//process
	if(FullOrNot>0) //output single alignment
	{
		int pre_ii,pre_jj;
		int wlen;
		int linear=1; //default: OK
		//check nonlinear
		pre_ii=-1;
		for(i=0;i<moln2;i++)
		{
			ii=ali2[i];
			if(ii==-1)continue;
			if(ii<pre_ii)
			{
				linear=0;
				break;
			}
			pre_ii=ii;
		}
		//output
		if(linear==0) //non-linear
		{
			fprintf(fp,"!non-linear!\n");
			fprintf(fp,"!non-linear!\n");
		}
		else
		{
			//init
			vector<pair<int,int> > alignment;
			alignment.clear();
			//start
			int j;
			pre_ii=0;
			pre_jj=0;
			for(i=1;i<=moln1;i++)
			{
				ii=i;
				jj=ali1[i-1];  //ali1 starts from 0, correspondence also from 0
				if(jj==-1)
				{
					continue;
				}
				else
				{
					jj++;
					//previous_path
					wlen=ii-pre_ii;
					for(j=1;j<wlen;j++)
					{
						pre_ii++;
						alignment.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
					}
					wlen=jj-pre_jj;
					for(j=1;j<wlen;j++)
					{
						pre_jj++;
						alignment.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
					}
					//current_path
					alignment.push_back (pair<int,int>(ii, jj)); //Match
					//update
					pre_ii=ii;
					pre_jj=jj;
				}
			}
			//termi
			pre_ii++;
			for(i=pre_ii;i<=moln1;i++)alignment.push_back (pair<int,int>(i, -pre_jj)); //Ix
			pre_jj++;
			for(i=pre_jj;i<=moln2;i++)alignment.push_back (pair<int,int>(-moln1, i));  //Iy
			int size=(int)alignment.size();
			//output_AMI
			for(i=0;i<size;i++)
			{
				ii=alignment[i].first;
				if(ii<=0)fprintf(fp,"-");
				else fprintf(fp,"%c",CH_AMI1[ii-1]);
			}
			fprintf(fp,"\n");
			for(i=0;i<size;i++)
			{
				jj=alignment[i].second;
				if(jj<=0)fprintf(fp,"-");
				else fprintf(fp,"%c",CH_AMI2[jj-1]);
			}
			fprintf(fp,"\n");
			//output_CLE
			for(i=0;i<size;i++)
			{
				ii=alignment[i].first;
				if(ii<=0)fprintf(fp,"-");
				else fprintf(fp,"%c",CH_CLE1[ii-1]);
			}
			fprintf(fp,"\n");
			for(i=0;i<size;i++)
			{
				jj=alignment[i].second;
				if(jj<=0)fprintf(fp,"-");
				else fprintf(fp,"%c",CH_CLE2[jj-1]);
			}
			fprintf(fp,"\n");
		}
	}
	if(FullOrNot>1) //output rotation matrix
	{
		fprintf(fp,"\n");
		fprintf(fp,"The rigid-body transformation used to superpose the first structure onto the second one is:\n");
		fprintf(fp,"/x'\\   / %9.6f  %9.6f  %9.6f \\   /x\\   / %12.6f \\ \n",CLEP_rotmat[0],CLEP_rotmat[1],CLEP_rotmat[2],CLEP_rotmat[9]);
		fprintf(fp,"|y'| = | %9.6f  %9.6f  %9.6f | * |y| + | %12.6f | \n",CLEP_rotmat[3],CLEP_rotmat[4],CLEP_rotmat[5],CLEP_rotmat[10]);
		fprintf(fp,"\\z'/   \\ %9.6f  %9.6f  %9.6f /   \\z/   \\ %12.6f / \n",CLEP_rotmat[6],CLEP_rotmat[7],CLEP_rotmat[8],CLEP_rotmat[11]);
	}
	if(FullOrNot>2) //output direct information
	{
		//output
		fprintf(fp,"\n");
		fprintf(fp,"TMscore=%5.3f, LALI=%4d, RMSD=%6.3f\n",align.main_sco,align.lali,align.RMSD);
		fprintf(fp,"Frag.  Score  Pos-1  Pos-2  Width \n");
		for(i=1;i<=totnum;i++)
		{
			score=AFP_Cor[i*4+0];
			ii=AFP_Cor[i*4+1];
			jj=AFP_Cor[i*4+2];
			winlen=AFP_Cor[i*4+3];
			fprintf(fp,"%5d  %5d  %5d  %5d  %5d \n",i,score,ii,jj,winlen);
		}
	}
	if(FullOrNot>3) //output full sequence
	{
		//sequence
		CLEPAPS_Out_rot_mol(mol1,AFP_tmp1,moln1,CLEP_rotmat);
		for(i=0;i<moln1;i++)ali1[i]=-1;
		for(i=0;i<moln2;i++)ali2[i]=-1;
		//assign
		string wsline="1   5    10   15   20   25   30   35   40   45   50   55   60   65   70   75   80";
		double wscore2=outcut2*outcut2;
		double wscore1=outcut1*outcut1;
		int k,l;
		XYZ temp;
		double dist2;
		for(k=1;k<=totnum;k++)
		{
			ii=AFP_Cor[k*4+1];
			jj=AFP_Cor[k*4+2];
			winlen=AFP_Cor[k*4+3];
			for(l=0;l<winlen;l++)
			{
				ali1[ii+l]=k;
				ali2[jj+l]=k;
				CLEPAPS_Out_rot_point(mol1[ii+l],temp,CLEP_rotmat);
				dist2=mol2[jj+l].distance_square(temp);
				if(dist2<wscore2)
				{
					ali1[ii+l]+=INT_MAX_NUM;
					ali2[jj+l]+=INT_MAX_NUM;
					if(dist2<wscore1)
					{
						ali1[ii+l]+=INT_MAX_NUM;
						ali2[jj+l]+=INT_MAX_NUM;
					}
				}
			}
		}
		//output_1
		int pos;
		char out1[7];
		char out2[7];
		FILE *fin_log=fp;
		int len=winlen;
		int wscut=80;
		int cur=0;
		int wmax=moln1;
		fprintf(fp,"\n");
		fprintf(fp,"       %s\n",wsline.c_str());
		for(;;)
		{
			len=wmax>=wscut?wscut:wmax;
			//output_STC
			fprintf(fin_log,"[%4d] ",cur+1);
			for(i=0;i<len;i++)
			{
				if(ali1[i+cur]==-1)fprintf(fin_log,"%c",CH_CLE1[i+cur]-'A'+'a');
				else fprintf(fin_log,"%c",CH_CLE1[i+cur]);
			}
			fprintf(fin_log,"\n");
			//output_AMI
			for(i=0;i<6;i++)out1[i]=CH_IND1[6*cur+i];
			out1[i]='\0';
			fprintf(fin_log,"%6s ",out1);
			for(i=0;i<len;i++)
			{
				if(ali1[i+cur]==-1)fprintf(fin_log,"%c",CH_AMI1[i+cur]-'A'+'a');
				else fprintf(fin_log,"%c",CH_AMI1[i+cur]);
			}
			fprintf(fin_log,"\n");
			//output_BLO
			fprintf(fin_log,"[pdb1] ");
			for(i=0;i<len;i++)
			{
				pos=ali1[i+cur];
				if(pos!=-1)
				{
					if(pos>INT_MAX_NUM)
					{
						pos-=INT_MAX_NUM;
						if(pos>INT_MAX_NUM)pos-=INT_MAX_NUM;
					}
					if(pos<=26)fprintf(fin_log,"%c",pos-1+'A');
					else if(pos<=52)fprintf(fin_log,"%c",pos-27+'a');
					else if(pos<=62)fprintf(fin_log,"%c",pos-53+'0');
					else fprintf(fin_log,"-");
				}
				else fprintf(fin_log," ");
			}
			fprintf(fin_log,"\n");
			//output_SIG
			fprintf(fin_log,"       ");
			for(i=0;i<len;i++)
			{
				pos=ali1[i+cur];
				if(pos!=-1)
				{
					if(pos>2*INT_MAX_NUM)fprintf(fin_log,"|");
					else if(pos>INT_MAX_NUM)fprintf(fin_log,":");
					else fprintf(fin_log,".");
				}
				else fprintf(fin_log," ");
			}
			fprintf(fin_log,"\n");
			//next
			if(wmax>=wscut)
			{
				wmax-=wscut;
				cur+=wscut;
			}
			else break;
			fprintf(fin_log,"\n");
		}
		fprintf(fin_log,"\n\n");
		//output_2
		cur=0;
		wmax=moln2;
		fprintf(fin_log,"       %s\n",wsline.c_str());
		for(;;)
		{
			len=wmax>=wscut?wscut:wmax;
			//output_STC
			fprintf(fin_log,"[%4d] ",cur+1);
			for(i=0;i<len;i++)
			{
				if(ali2[i+cur]==-1)fprintf(fin_log,"%c",CH_CLE2[i+cur]-'A'+'a');
				else fprintf(fin_log,"%c",CH_CLE2[i+cur]);
			}
			fprintf(fin_log,"\n");
			//output_AMI
			for(i=0;i<6;i++)out2[i]=CH_IND2[6*cur+i];
			out2[i]='\0';
			fprintf(fin_log,"%6s ",out2);
			for(i=0;i<len;i++)
			{
				if(ali2[i+cur]==-1) fprintf(fin_log,"%c",CH_AMI2[i+cur]-'A'+'a');
				else fprintf(fin_log,"%c",CH_AMI2[i+cur]);
			}
			fprintf(fin_log,"\n");
			//output_BLO
			fprintf(fin_log,"[pdb2] ");
			for(i=0;i<len;i++)
			{
				pos=ali2[i+cur];
				if(pos!=-1)
				{
					if(pos>INT_MAX_NUM)
					{
						pos-=INT_MAX_NUM;
						if(pos>INT_MAX_NUM)pos-=INT_MAX_NUM;
					}
					if(pos<=26)fprintf(fin_log,"%c",pos-1+'A');
					else if(pos<=52)fprintf(fin_log,"%c",pos-27+'a');
					else if(pos<=62)fprintf(fin_log,"%c",pos-53+'0');
					else fprintf(fin_log,"-");
				}
				else fprintf(fin_log," ");
			}
			fprintf(fin_log,"\n");
			//output_SIG
			fprintf(fin_log,"       ");
			for(i=0;i<len;i++)
			{
				pos=ali2[i+cur];
				if(pos!=-1)
				{
					if(pos>2*INT_MAX_NUM)fprintf(fin_log,"|");
					else if(pos>INT_MAX_NUM)fprintf(fin_log,":");
					else fprintf(fin_log,".");
				}
				else fprintf(fin_log," ");
			}
			fprintf(fin_log,"\n");
			//next
			if(wmax>=wscut)
			{
				wmax-=wscut;
				cur+=wscut;
			}
			else break;
			fprintf(fin_log,"\n");
		}
		fprintf(fin_log,"\n\n");
	}
}
void CLEPAPS_Out::CLEPAPS_Output_MolPDB(FILE *fp,Align_Record &align)
{
	int i;
	double rotmat[12];
	string TER="TER                                                                             ";
	for(i=0;i<12;i++)rotmat[i]=align.rotmat[i];
	CLEPAPS_Out_rot_pdb(CH_pdb1,CH_pdbt,CH_moln1,rotmat);
	Output_PDB_III(fp,CH_moln1,CH_pdbt,'A');
	fprintf(fp,"%s\n",TER.c_str());
	Output_PDB_III(fp,CH_moln2,CH_pdb2,'B');
	fprintf(fp,"%s\n",TER.c_str());
}

//================ main functions ============//
//---->main output:
//SingleOrNot -> output single solution or multi solution
//FullOrNot   -> output different level of alignment
//ScriptOrNot -> output RASMOL script
//[name] for main_output, [script] for mol_ca and RASMOL script
int CLEPAPS_Out::CLEPAPS_Main_Output(string &name,string &script,vector <Align_Record> &tot,
	int SingleOrNot,int FullOrNot,int ScriptOrNot)
{
	FILE *fp;
	string pdb_out;
	string ali_out;
	string ca_out;
	string sc_out;
	int size=(int)tot.size();
	if(size==0)
	{
		printf("CLEPAPS FAILED!!\n");
		return -1;
	}
	if(SingleOrNot==1) //single output
	{
		pdb_out=name+".pdb";
		fp=fopen(pdb_out.c_str(),"wb");
		CLEPAPS_Output_MolPDB(fp,tot.at(0));
		fclose(fp);
		ali_out=name+".ali";
		fp=fopen(ali_out.c_str(),"wb");
		fprintf(fp,">%s %d | %s %d\n",CH_NAM1.c_str(),CH_moln1,CH_NAM2.c_str(),CH_moln2);
		CLEPAPS_Output_Align(fp,tot.at(0),FullOrNot);
		fclose(fp);
		if(ScriptOrNot==1) //script output
		{
			ca_out=script+".mol";
			fp=fopen(ca_out.c_str(),"wb");
			CLEPAPS_Output_MolCA(fp,tot.at(0));
			fclose(fp);
			sc_out=script+".scp";
			fp=fopen(sc_out.c_str(),"wb");
			CLEPAPS_Out_Script_Simp(tot.at(0),CH_AFP_Cor,CH_moln1,CH_moln2);
			CLEPAPS_Output_Script(fp,ca_out,CH_AFP_Cor);
			fclose(fp);
		}
	}
	else               //multi output
	{
		int i;
		int totnum=size<outmaxnum?size:outmaxnum;
		char c;
		for(i=0;i<totnum;i++)
		{
			c='0'+i;
			pdb_out="";
			pdb_out=pdb_out+name+"_"+c+".pdb";
			fp=fopen(pdb_out.c_str(),"wb");
			CLEPAPS_Output_MolPDB(fp,tot.at(i));
			fclose(fp);
			ali_out="";
			ali_out=ali_out+name+"_"+c+".ali";
			fp=fopen(ali_out.c_str(),"wb");
			fprintf(fp,">%s %d | %s %d -> solution %d of %d \n",
				CH_NAM1.c_str(),CH_moln1,CH_NAM2.c_str(),CH_moln2,i+1,totnum);
			CLEPAPS_Output_Align(fp,tot.at(i),FullOrNot);
			fclose(fp);
			if(ScriptOrNot==1) //script output
			{
				ca_out="";
				ca_out=ca_out+script+"_"+c+".mol";
				fp=fopen(ca_out.c_str(),"wb");
				CLEPAPS_Output_MolCA(fp,tot.at(i));
				fclose(fp);
				sc_out="";
				sc_out=sc_out+script+"_"+c+".scp";
				fp=fopen(sc_out.c_str(),"wb");
				CLEPAPS_Out_Script_Simp(tot.at(i),CH_AFP_Cor,CH_moln1,CH_moln2);
				CLEPAPS_Output_Script(fp,ca_out,CH_AFP_Cor);
				fclose(fp);
			}
		}
	}
	return 1;
}
