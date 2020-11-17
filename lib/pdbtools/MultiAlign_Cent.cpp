#include "MultiAlign_Cent.h"

//--------------------- Start -----------------//
MultiAlign_Cent::MultiAlign_Cent(int tmnum)
:TM_align(tmnum)
{
	//macros
	TMorLEN=1;   //default: 1 (CORE_LEN*TMscore) -> 0, TMscore; 2, CORE_LEN
	RECorNOT=1;  //default: 1 (record initial value)
	ROTorNOT=1;  //default: 1 (do final consensus stage)
	//parameters
	CONSEN_STAGE=3; //[0], ROT stage; [1], FIX stage; [2], PIX stage
	REFINE_NUM=3;   //how many iterations for refinement
	MERGE_DIST=3.0; //column merge cutoff
}
MultiAlign_Cent::~MultiAlign_Cent(void)
{
}

//------------ create & delete --------------//
void MultiAlign_Cent::MultiAlign_Cent_Create(int num,int len)
{
	//ori_mol
	BC_INPUT_LEN=new int[num];
	NewArray2D(&BC_INPUT_MOL,num,len);
	//consensus
	BLOCEN_Consensus=new XYZ[num*len];
	BLOCEN_Consenrec=new int[num*len];
	NewArray2D(&BC_Condist,num,len);
	BC_ZeroOne_Cur=new int[num];
	BC_ZeroOne_Tmp=new int[num];
	BC_ZeroOne_Len=new int[num];
	//refinement
	BC_Tmp1=new XYZ[num*len];
	BC_Tmp2=new XYZ[num*len];
	NewArray2D(&BC_Ali2_Refine,num,num*len);
	NewArray2D(&BC_Ali_Temp,num,num*len);
	//temp
	BC_ali1=new int[num*len];
	BC_ali2=new int[num*len];
	BC_alib=new int[num*len];
	BC_best=new XYZ[num*len];
	//best
	NewArray2D(&BC_Ali_Best,num,num*len);
	BC_Best_Consensus=new XYZ[num*len];
	NewArray2D(&BC_Best_Output,num,len);
}
void MultiAlign_Cent::MultiAlign_Cent_Delete(int num,int len)
{
	//ori_mol
	delete [] BC_INPUT_LEN;
	DeleteArray2D(&BC_INPUT_MOL,num);
	//consensus
	delete [] BLOCEN_Consensus;
	delete [] BLOCEN_Consenrec;
	DeleteArray2D(&BC_Condist,num);
	delete [] BC_ZeroOne_Cur;
	delete [] BC_ZeroOne_Tmp;
	delete [] BC_ZeroOne_Len;
	//refinement
	delete [] BC_Tmp1;
	delete [] BC_Tmp2;
	DeleteArray2D(&BC_Ali2_Refine,num);
	DeleteArray2D(&BC_Ali_Temp,num);
	//temmp
	delete [] BC_ali1;
	delete [] BC_ali2;
	delete [] BC_alib;
	delete [] BC_best;
	//best
	DeleteArray2D(&BC_Ali_Best,num);
	delete [] BC_Best_Consensus;
	DeleteArray2D(&BC_Best_Output,num);
}

//===================== consensus main process =================//
//[input]
//conlen  -> the one protein's length
//ali2out -> all alignment with one protein 
//[output]
//return -> normal_style MSA's length
//aliout -> normal_style MSA
int MultiAlign_Cent::Single_Ali_To_Multi_Ali(int totnum,int *len,int conlen,int **ali2out,int **aliout)
{
	//once-a-gap,always-a-gap
	int i,j,k,l;
	int cur=0;
	int count;
	int pos;
	int gap;
	for(j=0;j<totnum;j++)BC_ZeroOne_Cur[j]=-1;
	for(i=0;i<conlen;i++)
	{
		//collect gap
		count=0;
		for(j=0;j<totnum;j++)
		{
			pos=ali2out[j][i];
			if(pos!=-1)
			{
				gap=pos-BC_ZeroOne_Cur[j];
				if(gap<=0)
				{
//					int wsiii=0;
					fprintf(stderr,"BAD HERE!!!\n"); //this is IMPOSSIBLE!! (non-linear alignment)
				}
				else if(gap>1)
				{
					//assign
					for(k=0;k<gap-1;k++)
					{
						for(l=0;l<totnum;l++)
						{
							if(l!=j)aliout[l][cur]=0;
							else
							{
								aliout[l][cur]=1;
								BC_ZeroOne_Cur[j]++;
							}
						}
						cur++;
					}
				}
				count++;
			}
		}
		//align
		if(count>0)
		{
			for(j=0;j<totnum;j++)
			{
				pos=ali2out[j][i];
				if(pos!=-1)
				{
					if(count>1)
					{
						aliout[j][cur]=2;
						BC_ZeroOne_Cur[j]++;
					}
					else
					{
						aliout[j][cur]=1;
						BC_ZeroOne_Cur[j]++;
					}
				}
				else aliout[j][cur]=0;
			}
			cur++;
		}
	}
	//terminal//__101120__//
	int terlen;
	for(j=0;j<totnum;j++)
	{
		terlen=len[j];
		pos=BC_ZeroOne_Cur[j];
		for(k=pos;k<terlen-1;k++)
		{
			for(l=0;l<totnum;l++)
			{
				if(l!=j)aliout[l][cur]=0;
				else
				{
					aliout[l][cur]=1;
					BC_ZeroOne_Cur[j]++;
				}
			}
			cur++;
		}
	}
	//terminal_over//
	return cur;
}

//======================== consensus vice process ====================//
//[combine close point]
int MultiAlign_Cent::BC_Consen_Refine(XYZ *in,int *rec,int totnum,int totlen,double thres)
{
	int i,j;
	//second refinement -> kill close
	XYZ temp;
	double dist2;
	double thres2=thres*thres;
	for(i=0;i<totlen;i++)
	{
		if(rec[i]==totnum)continue;
		if(rec[i]<=0)continue;
		for(j=i+1;j<totlen;j++)
		{
			if(rec[j]==totnum)break;
			if(rec[j]<=0)continue;
			dist2=in[i].distance_square(in[j]);
			if(dist2<thres2 && rec[i]+rec[j]<=totnum)
			{
				temp=in[i]*rec[i]+in[j]*rec[j];
				temp/=(rec[i]+rec[j]);
				in[i]=temp;
				rec[i]+=rec[j];
				rec[j]*=-1;
				if(rec[i]==totnum)break;
			}
		}
	}
	int cur=0;
	for(i=0;i<totlen;i++)
	{
		if(rec[i]>0)
		{
			in[cur]=in[i];
			rec[cur]=rec[i];
			cur++;
		}
	}
	return cur;
}
//[calculate SumTM given mol]
double MultiAlign_Cent::BC_Calc_SumTM_Given_Mol(XYZ **in,int *len,int totnum,int totlen,int **ali)
{
	//[temp]
	int pos1,pos2;
	int length;
	//[calc]
	int i,j,k;
	int count;
	int smaller;
	double totscore;
	double score;
	count=0;
	totscore=0.0;
	for(i=0;i<totnum;i++)
	{
		for(j=i+1;j<totnum;j++)
		{
			count++;
			//calc d0
			smaller=len[i]<len[j]?len[i]:len[j];
			if(smaller==0)continue;
			Calc_TM_d0(smaller);
			//collect point
			pos1=0;
			pos2=0;
			length=0;
			for(k=0;k<totlen;k++)
			{
				if(ali[i][k]>0&&ali[j][k]>0) //record
				{
					BC_Tmp1[length]=in[i][pos1];
					BC_Tmp2[length]=in[j][pos2];
					//next
					pos1++;
					pos2++;
					length++;
					
					//judge
					if(pos1>len[i] || pos2>len[j]) //that's impossible!!
					{
//						int wsiii=0;
						fprintf(stderr,"BAD!!!!!!!TM_score_simp!!!!\n");
					}
				}
				else
				{
					if(ali[i][k]>0)pos1++;
					if(ali[j][k]>0)pos2++;

					//judge
					if(pos1>len[i] || pos2>len[j]) //that's impossible!!
					{
//						int wsiii=0;
						fprintf(stderr,"BAD!!!!!!!TM_score_simp!!!!\n");
					}
				}
			}
			//calc TMscore
			score=Calc_TM_Score_Simple(BC_Tmp1,BC_Tmp2,length,d0,d8);
			totscore+=1.0*score/smaller;
		}
	}
	//return
	if(count==0)return 0.0;
	else return 1.0*totscore/count;
}

//=================== main_process ===================//[fixed!!]
//[given multi-aligned structure, return consensus]
int MultiAlign_Cent::BC_Update_Consen_Given_Ali_Fixed(XYZ **in,int totnum,int totlen,int **ali,XYZ *consen,int *conrec,int cutoff)
{
	int i,j;
	int pos;
	int cur;
	int count;
	XYZ temp;

	//first_calculate
	for(j=0;j<totnum;j++)BC_ZeroOne_Cur[j]=0;
	cur=0;
	for(i=0;i<totlen;i++)
	{
		//calc consensus
		temp=0.0;
		count=0;
		for(j=0;j<totnum;j++)
		{
			if(ali[j][i]>0)
			{
				pos=BC_ZeroOne_Cur[j];
				temp+=in[j][pos];
				count++;
				BC_ZeroOne_Cur[j]++;
			}
		}
		//assign
		if(count>cutoff) //__modified at//__110230__//
		{
			temp/=count;
			consen[cur]=temp;
			conrec[cur]=count;
			cur++;
		}
	}
	return cur;
}
//[given consensus and alignment, update multi-structure's alignment]
int MultiAlign_Cent::BC_Update_Ali_Gigen_Consen_Fixed(XYZ **in,int *len,int totnum,int conlen,XYZ *consen,int **aliout,int ROTorNOT)
{
	int i,j;
	int smaller;
	double tms;
	double rms;
	int lali;
	//update
	for(j=0;j<totnum;j++)
	{
		smaller=len[j]<conlen?len[j]:conlen;
		if(smaller==0)
		{
			for(i=0;i<conlen;i++)BC_Ali2_Refine[j][i]=-1;
			continue;
		}
		TM_Align_Get_Ali(in[j],consen,len[j],conlen,BC_Ali2_Refine[j]);
		if(smaller>=3)
		{
//			TM_GAP_TYPE=1;  //advance three-layer four-state DynaProg
//			TM_GAP_STAGE=1; //simple TMalign
//			TM_GAP_TYPE=0;  //normal one-layer DynaProg
//			TM_GAP_STAGE=2; //normal TMalign
//			Calc_TM_Align(in[j],consen,len[j],conlen,BC_Ali2_Refine[j],BC_Ali2_Refine[j]);
			//tmalign_start
			TM_Align_Init(len[j],conlen);
//			TM_GAP_STAGE=2; //normal TMalign
			TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
			Calc_TM_Align(in[j],consen,len[j],conlen,BC_Ali2_Refine[j],BC_Ali2_Refine[j],0,3);    //apply TMscore dynamic programming
			TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
			//tmalign_over
			tms=TM_Align_TM_Score_Simp(in[j],consen,len[j],conlen,BC_Ali2_Refine[j],rms,lali); //use simple version
			if(ROTorNOT==1)rot_mol(in[j],in[j],len[j],finmat);
		}
	}
	return Single_Ali_To_Multi_Ali(totnum,len,conlen,BC_Ali2_Refine,aliout);
}

//---- conserved core calculate ----//
int MultiAlign_Cent::BC_Return_Conserved_Core(int totnum,int totlen,int **ali)
{
	int i,j;
	int correct;
	int len=0;
	for(i=0;i<totlen;i++)
	{
		correct=1;
		for(j=0;j<totnum;j++)
		{
			if(ali[j][i]==0)
			{
				correct=0;
				break;
			}
		}
		if(correct==1)len++;
	}
	return len;
}
void MultiAlign_Cent::BC_Return_Conserved_Core_Pos(int totnum,int totlen,int **ali,int *core_pos)
{
	int i,j;
	for(i=0;i<totlen;i++)
	{
		int number=0;
		for(j=0;j<totnum;j++)
		{
			if(ali[j][i]!=0)number++;
		}
		core_pos[i]=number;
	}
}
//---- conserved rmsd calculate ----//
double MultiAlign_Cent::BC_Return_Conserved_RMSD(int totnum,int totlen,int **ali,XYZ **mol)
{
	int i;
	int j,k;
	int correct;
	int pos1,pos2;
	double dist2;
	double col_rmsd;
	int col_count;
	double ret_rmsd;
	int ret_cur;
	ret_rmsd=0;
	ret_cur=0;
	//init
	for(j=0;j<totnum;j++)BC_ZeroOne_Cur[j]=-1;
	//process
	for(i=0;i<totlen;i++)
	{
		//check
		correct=1;
		for(j=0;j<totnum;j++)
		{
			if(ali[j][i]==0)correct=0;
			else BC_ZeroOne_Cur[j]++;
		}
		if(correct==0)continue;
		//calc
		col_rmsd=0;
		col_count=0;
		for(j=0;j<totnum;j++)
		{
			for(k=j+1;k<totnum;k++)
			{
				pos1=BC_ZeroOne_Cur[j];
				pos2=BC_ZeroOne_Cur[k];
				dist2=mol[j][pos1].distance_square(mol[k][pos2]);
				col_rmsd+=dist2;
				col_count++;
			}
		}
		col_rmsd/=col_count;
		ret_rmsd+=col_rmsd;
		ret_cur++;
	}
	if(ret_cur>0)
	{
		ret_rmsd/=ret_cur;
		ret_rmsd=sqrt(ret_rmsd);
	}
	return ret_rmsd;
}
void MultiAlign_Cent::BC_Return_Conserved_RMSD_Pos(int totnum,int totlen,int **ali,XYZ **mol,double *rmsd_pos)
{
	int i;
	int j,k;
	int pos1,pos2;
	double dist2;
	double col_rmsd;
	int col_count;
	//init
	for(j=0;j<totnum;j++)BC_ZeroOne_Cur[j]=-1;
	//process
	for(i=0;i<totlen;i++)
	{
		//init result
		rmsd_pos[i]=-1;
		//get pos
		for(j=0;j<totnum;j++)
		{
			if(ali[j][i]!=0)BC_ZeroOne_Cur[j]++;
		}
		//calc
		col_rmsd=0;
		col_count=0;
		for(j=0;j<totnum;j++)
		{
			if(ali[j][i]==0)continue;
			for(k=j+1;k<totnum;k++)
			{
				if(ali[k][i]==0)continue;
				pos1=BC_ZeroOne_Cur[j];
				pos2=BC_ZeroOne_Cur[k];
				dist2=mol[j][pos1].distance_square(mol[k][pos2]);
				col_rmsd+=dist2;
				col_count++;
			}
		}
		//final result
		if(col_count==0)continue;
		col_rmsd/=col_count;
		rmsd_pos[i]=sqrt(col_rmsd);
	}
}

//=================== main_process ===================// [final rotation part]
//[given multi-aligned structure, return consensus]
//[input]
//totlen -> normal_style MSA's length
//aliout -> normal_style MSA
//[output]
//return -> the one protein's length (this's the consensus)
//ali2out-> all alignment with one protein 
int MultiAlign_Cent::BC_Update_Consen_Given_Ali(XYZ **in,int *len,int totnum,int totlen,int *conres,double **conrec,int **aliout,int **ali2out,XYZ *consen)
{
	int i,j;
	int pos;
	int cur;
	int count;
	double dist2;
	XYZ temp;

	//init
	for(i=0;i<totnum;i++)for(j=0;j<len[i];j++)conrec[i][j]=INT_MAX_NUM;
	for(i=0;i<totnum;i++)for(j=0;j<totlen;j++)ali2out[i][j]=-1;
	//first_calculate
	for(j=0;j<totnum;j++)BC_ZeroOne_Cur[j]=0;
	cur=0;
	for(i=0;i<totlen;i++)
	{
		//calc consensus
		temp=0.0;
		count=0;
		for(j=0;j<totnum;j++)
		{
			if(aliout[j][i]>0)
			{
				pos=BC_ZeroOne_Cur[j];
				temp+=in[j][pos];
				count++;
				BC_ZeroOne_Cur[j]++;
			}
		}
		//assign
		if(count>0)
		{
			temp/=count;
			consen[cur]=temp;
			conres[cur]=count;
			cur++;
		}
		//record
		if(count>1)
		{
			for(j=0;j<totnum;j++)
			{
				if(aliout[j][i]>0)
				{
					pos=BC_ZeroOne_Cur[j]-1;
					dist2=consen[cur-1].distance_square(in[j][pos]);
					if(conres[cur-1]>1)conrec[j][pos]=dist2;
					ali2out[j][cur-1]=pos;
				}
			}
		}
	}
	return cur;
}
//[given consensus and alignment, update multi-structure's alignment]
//mol2 is fixed (just the consensus)
double MultiAlign_Cent::BC_Update_Ali_Gigen_Consen_Single(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *res,double *rec,int *ali2out)
{
	//real_start
	int i;
	int pos;
	int k;
	double dist2;
	double maxtms;
	double tms;
	int *ali1=BC_ali1;
	int *ali2=BC_ali2;
	int *alib=BC_alib;
	XYZ *best=BC_best;
	double rotmat[12];
//	double epsilu=1.0;

	//initial
	int smaller;
	double ori_d;
	int ii,jj;
	smaller=moln1;
	if(smaller==0)
	{
		for(i=0;i<moln2;i++)ali2out[i]=-1;
		return 0.0;
	}
	Calc_TM_d0(smaller);
	ori_d=d0*d0;
	EqualArray(ali2,ali2out,moln2);
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]>=0)
		{
			ii=ali2[i];
			jj=i;
			ali1[ii]=jj;
		}
	}
	maxtms=0.0;
	for(i=0;i<moln1;i++)
	{
		pos=ali1[i];
		if(pos==-1)continue;
		dist2=rec[i];
		if(res[pos]>1)maxtms+=1.0/(1.0+dist2/ori_d);
	}
	EqualArray(alib,ali2,moln2);
	EqualArray(best,mol1,moln1);

	//find maximal
	for(k=0;k<1;k++)
	{
		//get_update matrix
		TM_Align_Get_Ali(mol1,mol2,moln1,moln2,ali2);
		//ws_method[1] -> using 30 TM iterations
//		TM_GAP_TYPE=1;  //advance three-layer four-state DynaProg
//		TM_GAP_STAGE=1; //simple TMalign
//		TM_GAP_TYPE=0;  //normal one-layer DynaProg
//		TM_GAP_STAGE=2; //normal TMalign
//		Calc_TM_Align(mol1,mol2,moln1,moln2,ali2,ali2); // several TM_rotaion //__101120__//
		//tmalign_start
		TM_Align_Init(moln1,moln2);
//		TM_GAP_STAGE=2; //normal TMalign
		TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
		Calc_TM_Align(mol1,mol2,moln1,moln2,ali2,ali2,0,3); // several TM_rotaion //__101120__//
		TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
		//tmalign_over
		EqualArray(rotmat,TM_rotmat,12);
		//ws_method[2] -> kill those 'bad' positions
		rot_mol(mol1,mol1,moln1,rotmat);
		//tmscore_calc
		for(i=0;i<moln1;i++)ali1[i]=-1;
		for(i=0;i<moln2;i++)
		{
			if(ali2[i]>=0)
			{
				ii=ali2[i];
				jj=i;
				ali1[ii]=jj;
			}
		}
		tms=0.0;
		for(i=0;i<moln1;i++)
		{
			pos=ali1[i];
			if(pos==-1)continue;
			dist2=mol1[i].distance_square(mol2[pos]);
			if(res[pos]>1)tms+=1.0/(1.0+dist2/ori_d);
		}
		//final_judge
		if(tms>maxtms)
		{
			maxtms=tms;
			EqualArray(alib,ali2,moln2);
			EqualArray(best,mol1,moln1);
			k=0;
		}
	}
	EqualArray(ali2out,alib,moln2);
	EqualArray(mol1,best,moln1);

	//return
	maxtms=maxtms/smaller;
	return maxtms;
}
int MultiAlign_Cent::BC_Update_Ali_Gigen_Consen(XYZ **in,int *len,int totnum,int conlen,int *conres,double **conrec,int **aliout,XYZ *consen)
{
	int j;
	double tms;
	//update
	for(j=0;j<totnum;j++)
	{
		//the following procedure only updates the ali2out (consensus_single_correspondence)
		tms=BC_Update_Ali_Gigen_Consen_Single(in[j],consen,len[j],conlen,conres,conrec[j],BC_Ali2_Refine[j]);
	}
	return Single_Ali_To_Multi_Ali(totnum,len,conlen,BC_Ali2_Refine,aliout);
}

//--------------- Update_Partial_CORE -------------//__110115__//
void MultiAlign_Cent::BC_Update_Partial_CORE_Single(XYZ **in,int *len,int totnum,int **alin,int totin,int **aliout,int &totout)
{
	int j,k;
	totout=totin;
	EqualArray2D(aliout,alin,totnum,totin);
	//frist judge
	int member=0;
	for(j=0;j<totnum;j++)
	{
		if(len[j]>0)member++;
	}
	if(member<2)return;
	//real process
	int conlen;
	double cur_score;
	double max_score=-1;
	double ws_epsilu=0.0001;
	int curlen=totin;
	EqualArray2D(BC_Ali_Temp,alin,totnum,curlen);
	for(k=1;k<=1;k++)
	{
		conlen=BC_Update_Consen_Given_Ali_Fixed(in,totnum,curlen,BC_Ali_Temp,BLOCEN_Consensus,BLOCEN_Consenrec,0);
		conlen=BC_Consen_Refine(BLOCEN_Consensus,BLOCEN_Consenrec,totnum,curlen,MERGE_DIST);   //this value [3.0] may be modified...//__101115__//
		curlen=BC_Update_Ali_Gigen_Consen_Fixed(in,len,totnum,conlen,BLOCEN_Consensus,BC_Ali_Temp);
		cur_score=BC_Calc_SumTM_Given_Mol(in,len,totnum,curlen,BC_Ali_Temp);
		if(cur_score>max_score+ws_epsilu)
		{
			//record alignment
			max_score=cur_score;
			totout=curlen;
			EqualArray2D(aliout,BC_Ali_Temp,totnum,curlen);
			k=1;
		}
	}
}
void MultiAlign_Cent::BC_Update_Partial_CORE(int pivot,XYZ **in,int *len,int totnum,int **alin,int totin,int **aliout,int &totout)
{
	int i,j,k;
	//init cur position
	int isFirst=1;
	int last=0;    //current pos
	int curlen=0;  //process length
	int retlen=0;  //return length
	for(j=0;j<totnum;j++)BC_ZeroOne_Tmp[j]=0;
	for(j=0;j<totnum;j++)BC_ZeroOne_Len[j]=0;
	//real process
	int cur=0;
	for(i=0;i<totin;i++)
	{
		if(alin[pivot][i]==0)  //current part should refine!!
		{
			//init
			if(isFirst==1)
			{
				isFirst=0;
				last=i;
				curlen=0;
				retlen=0;
				for(j=0;j<totnum;j++)BC_ZeroOne_Len[j]=0;
			}
			//record
			for(j=0;j<totnum;j++)
			{
				if(alin[j][i]>0)BC_ZeroOne_Len[j]++;
			}
			curlen++;
		}
		else
		{
			//process
			if(isFirst==0)
			{
				//collect
				for(j=0;j<totnum;j++)
				{
					for(k=0;k<BC_ZeroOne_Len[j];k++)BC_INPUT_MOL[j][k]=in[j][BC_ZeroOne_Tmp[j]+k];
					BC_INPUT_LEN[j]=BC_ZeroOne_Len[j];
				}
				for(k=0;k<curlen;k++)
				{
					for(j=0;j<totnum;j++)BC_Ali_Best[j][k]=alin[j][last+k];
				}
				BC_Update_Partial_CORE_Single(BC_INPUT_MOL,BC_INPUT_LEN,totnum,BC_Ali_Best,curlen,BC_Ali_Best,retlen);
				//add_previous
				for(k=0;k<retlen;k++)
				{
					for(j=0;j<totnum;j++)
					{
						aliout[j][cur]=BC_Ali_Best[j][k];
					}
					cur++;
				}
				//init
				isFirst=1;
				for(j=0;j<totnum;j++)
				{
					BC_ZeroOne_Tmp[j]+=BC_ZeroOne_Len[j];
				}
			}
			//add_current
			for(j=0;j<totnum;j++)
			{
				aliout[j][cur]=alin[j][i];
				if(alin[j][i]>0)BC_ZeroOne_Tmp[j]++;
			}
			cur++;
		}
	}
	//last process
	if(isFirst==0)
	{
		//collect
		for(j=0;j<totnum;j++)
		{
			for(k=0;k<BC_ZeroOne_Len[j];k++)BC_INPUT_MOL[j][k]=in[j][BC_ZeroOne_Tmp[j]+k];
			BC_INPUT_LEN[j]=BC_ZeroOne_Len[j];
		}
		for(k=0;k<curlen;k++)
		{
			for(j=0;j<totnum;j++)BC_Ali_Best[j][k]=alin[j][last+k];
		}
		BC_Update_Partial_CORE_Single(BC_INPUT_MOL,BC_INPUT_LEN,totnum,BC_Ali_Best,curlen,BC_Ali_Best,retlen);
		//add_previous
		for(k=0;k<retlen;k++)
		{
			for(j=0;j<totnum;j++)
			{
				aliout[j][cur]=BC_Ali_Best[j][k];
			}
			cur++;
		}
	}
	//termi
	totout=cur;
}




//========================= main_function ===================//
double MultiAlign_Cent::BC_Tota_Update_Function(XYZ **in,int *len,int totnum,int totlen,int **ali) //this is virtual function
{
	double cur_score=BC_Calc_SumTM_Given_Mol(in,len,totnum,totlen,ali);
	int cur_totlen=BC_Return_Conserved_Core(totnum,totlen,ali);
	if(TMorLEN==1)return (cur_score+0.01)*(cur_totlen+1);  //return CORE_LEN*Ave_TMscore
	else if(TMorLEN==0)return cur_score;  //return Ave_TMscore
	else return cur_totlen;               //return CORE_LEN
}
double MultiAlign_Cent::BC_Tota_Update_Main(XYZ **in,int *len,int totnum,int **alin,int totin,int **aliout,int &totout)
{
	int i,j,k;
	int curlen;
	double cur_score;
	int maxlen;
	double max_score;
	int conlen;
	double ws_epsilu=0.0001;

	//init
	BC_INPUT_NUM=totnum;
	EqualArray(BC_INPUT_LEN,len,totnum);
	for(i=0;i<totnum;i++)for(j=0;j<len[i];j++)BC_INPUT_MOL[i][j]=in[i][j];

	//start
	BC_Best_Consenlen=-1;
	max_score=-1;
	maxlen=totin;
	EqualArray2D(BC_Ali_Best,alin,totnum,totin);
	if(RECorNOT==1)max_score=BC_Tota_Update_Function(BC_INPUT_MOL,BC_INPUT_LEN,BC_INPUT_NUM,maxlen,BC_Ali_Best);
	for(i=0;i<totnum;i++)for(j=0;j<len[i];j++)BC_Best_Output[i][j]=BC_INPUT_MOL[i][j];

	//prefixed_version[fix1] -> All the structures are fixed, while update the consensus
	if(CONSEN_STAGE>2)
	{
		curlen=totin;
		for(i=0;i<totnum;i++)for(j=0;j<len[i];j++)BC_INPUT_MOL[i][j]=in[i][j];
		EqualArray2D(BC_Ali_Temp,alin,totnum,curlen);
		for(k=1;k<=REFINE_NUM;k++)
		{
			conlen=BC_Update_Consen_Given_Ali_Fixed(BC_INPUT_MOL,BC_INPUT_NUM,curlen,BC_Ali_Temp,BLOCEN_Consensus,BLOCEN_Consenrec);
			conlen=BC_Consen_Refine(BLOCEN_Consensus,BLOCEN_Consenrec,BC_INPUT_NUM,curlen,MERGE_DIST);   //this value [3.0] may be modified...//__101115__//
			curlen=BC_Update_Ali_Gigen_Consen_Fixed(BC_INPUT_MOL,BC_INPUT_LEN,BC_INPUT_NUM,conlen,BLOCEN_Consensus,BC_Ali_Temp,0); //NOT rotation!!
			cur_score=BC_Tota_Update_Function(BC_INPUT_MOL,BC_INPUT_LEN,BC_INPUT_NUM,curlen,BC_Ali_Temp);
			if(cur_score>max_score+ws_epsilu)
			{
				//record alignment
				max_score=cur_score;
				maxlen=curlen;
				EqualArray2D(BC_Ali_Best,BC_Ali_Temp,totnum,curlen);
				printf("PIX_BEST->%lf\r",max_score);
				//record structure
				BC_Best_Consenlen=conlen;
				EqualArray(BC_Best_Consensus,BLOCEN_Consensus,conlen);
				for(i=0;i<totnum;i++)for(j=0;j<len[i];j++)BC_Best_Output[i][j]=BC_INPUT_MOL[i][j];
				k=1;
			}
			if(BC_Best_Consenlen==-1)
			{
				BC_Best_Consenlen=conlen;
				EqualArray(BC_Best_Consensus,BLOCEN_Consensus,conlen);
			}
		}
	}
	
	//fixed_version[fix2] -> All the structures are rotated against to the consensus, while update the consensus
	if(CONSEN_STAGE>1)
	{
		curlen=totin;
		for(i=0;i<totnum;i++)for(j=0;j<len[i];j++)BC_INPUT_MOL[i][j]=in[i][j];
		EqualArray2D(BC_Ali_Temp,alin,totnum,curlen);
		for(k=1;k<=REFINE_NUM;k++)
		{
			conlen=BC_Update_Consen_Given_Ali_Fixed(BC_INPUT_MOL,BC_INPUT_NUM,curlen,BC_Ali_Temp,BLOCEN_Consensus,BLOCEN_Consenrec);
			conlen=BC_Consen_Refine(BLOCEN_Consensus,BLOCEN_Consenrec,BC_INPUT_NUM,curlen,MERGE_DIST);   //this value [3.0] may be modified...//__101115__//
			curlen=BC_Update_Ali_Gigen_Consen_Fixed(BC_INPUT_MOL,BC_INPUT_LEN,BC_INPUT_NUM,conlen,BLOCEN_Consensus,BC_Ali_Temp,1); //rotation!!
			cur_score=BC_Tota_Update_Function(BC_INPUT_MOL,BC_INPUT_LEN,BC_INPUT_NUM,curlen,BC_Ali_Temp);
			if(cur_score>max_score+ws_epsilu)
			{
				//record alignment
				max_score=cur_score;
				maxlen=curlen;
				EqualArray2D(BC_Ali_Best,BC_Ali_Temp,totnum,curlen);
				printf("FIX_BEST->%lf\r",max_score);
				//record structure
				BC_Best_Consenlen=conlen;
				EqualArray(BC_Best_Consensus,BLOCEN_Consensus,conlen);
				for(i=0;i<totnum;i++)for(j=0;j<len[i];j++)BC_Best_Output[i][j]=BC_INPUT_MOL[i][j];
				k=1;
			}
			if(BC_Best_Consenlen==-1)
			{
				BC_Best_Consenlen=conlen;
				EqualArray(BC_Best_Consensus,BLOCEN_Consensus,conlen);
			}
		}
	}

	//rotation_version[upon the best] -> furthur refine on the best MSA
	if(ROTorNOT==1)
	{
		if(CONSEN_STAGE>0)
		{
			curlen=maxlen;
			for(i=0;i<totnum;i++)for(j=0;j<len[i];j++)BC_INPUT_MOL[i][j]=BC_Best_Output[i][j];
			EqualArray2D(BC_Ali_Temp,BC_Ali_Best,totnum,maxlen);
			for(k=1;k<=REFINE_NUM;k++)
			{
				conlen=BC_Update_Consen_Given_Ali(BC_INPUT_MOL,BC_INPUT_LEN,BC_INPUT_NUM,curlen,BLOCEN_Consenrec,BC_Condist,BC_Ali_Temp,BC_Ali2_Refine,BLOCEN_Consensus);
				curlen=BC_Update_Ali_Gigen_Consen(BC_INPUT_MOL,BC_INPUT_LEN,BC_INPUT_NUM,conlen,BLOCEN_Consenrec,BC_Condist,BC_Ali_Temp,BLOCEN_Consensus);
				cur_score=BC_Tota_Update_Function(BC_INPUT_MOL,BC_INPUT_LEN,BC_INPUT_NUM,curlen,BC_Ali_Temp);
				if(cur_score>max_score+ws_epsilu)
				{
					//record alignment
					max_score=cur_score;
					maxlen=curlen;
					EqualArray2D(BC_Ali_Best,BC_Ali_Temp,totnum,curlen);
					printf("ROT_BEST->%lf\r",max_score);
					//record structure
					BC_Best_Consenlen=conlen;
					EqualArray(BC_Best_Consensus,BLOCEN_Consensus,conlen);
					for(i=0;i<totnum;i++)for(j=0;j<len[i];j++)BC_Best_Output[i][j]=BC_INPUT_MOL[i][j];
					k=1;
				}
				if(BC_Best_Consenlen==-1)
				{
					BC_Best_Consenlen=conlen;
					EqualArray(BC_Best_Consensus,BLOCEN_Consensus,conlen);
				}
			}
		}
	}

	//return
	totout=maxlen;
	for(i=0;i<totnum;i++)for(j=0;j<len[i];j++)BC_INPUT_MOL[i][j]=BC_Best_Output[i][j];
	EqualArray2D(aliout,BC_Ali_Best,BC_INPUT_NUM,maxlen);
	return max_score;
}
