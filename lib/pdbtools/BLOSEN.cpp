#include "BLOSEN.h"

BLOSEN::BLOSEN(int num)
:Ali_Ori(num),TM_align(num),
ZoomIn_Align(num),MultiAlign_Cent(num)
{
	BLOSEN_Maximal=num;
	SHAVE_THRES=0.5; //IMPORTANT!! (grid shave conflict threshold)
	COMP_THRES=4;    //for MSA comparison
	CORE_THRES=0.5;  //for MSA comparison
}
BLOSEN::~BLOSEN(void)
{
}

//------------ create & delete --------------//
void BLOSEN::BLOSEN_Create(int num,int len)
{
	//[SFP]
	BN_SFP_H.resize(num);
	BN_SFP_L.resize(num);
	//[HSFB]
	Fast_Sort_Val=new int[len*num];
	Fast_Sort_Idx=new int[len*num];
	NewArray2D(&Shav_Grid,num,len);
	//[Anchor]
	NewArray2D(&INPUT_MOL_Rot,num,len);
	BN_tmp1=new XYZ[len];
	BN_tmp2=new XYZ[len];
	Pivot_Ali_TTT=new int[len];
	Check_ali1=new int[len];
	Check_ali2=new int[len];
}
void BLOSEN::BLOSEN_Delete(int num,int len)
{
	//[HSFB]
	delete [] Fast_Sort_Val;
	delete [] Fast_Sort_Idx;
	DeleteArray2D(&Shav_Grid,num);
	//[Anchor]
	DeleteArray2D(&INPUT_MOL_Rot,num);
	delete [] BN_tmp1;
	delete [] BN_tmp2;
	delete [] Pivot_Ali_TTT;
	delete [] Check_ali1;
	delete [] Check_ali2;
}

//------------ utility funcion --------//
//Compare Two MSAs
double BLOSEN::Compare_Two_MSA(vector < vector <int> > &msa1,vector < vector <int> > &msa2,
	int len1,int len2,int totnum)
{
	int i,j,k;
	int c1,c2;
	c1=1;
	c2=1;
	double thres=COMP_THRES;
	//create
	double *score_rec=new double[len1*len2];
	memset(score_rec,0,sizeof(double)*len1*len2);
	//assign
	c1=1;
	for(i=0;i<len1;i++)
	{
		//check for core// start
		int number1=0;
		for(k=0;k<totnum;k++)if(msa1[k][i]>0)number1++;
		if(number1<totnum*CORE_THRES)continue;
		c1++;
		//check for core// end

		int cur_index=i*len2;
		c2=1;
		for(j=0;j<len2;j++)
		{
	                //check for core// start
	                int number2=0;
	                for(k=0;k<totnum;k++)if(msa2[k][j]>0)number2++;
	                if(number2<totnum*CORE_THRES)continue;
			c2++;
	                //check for core// end

			//calculate score
			double score=0;
			int dist;
			int count=0;
			for(k=0;k<totnum;k++)
			{
				if(msa1[k][i]>0 && msa2[k][j]>0)
				{
					dist=abs(msa1[k][i]-msa2[k][j]);
					score+=(thres-dist);
					count++;
				}
			}
			if(count!=0)score=10.0*score/count;
			else score=0;
			//assign score
			score_rec[cur_index+j]=score;
		}
	}
	//calculate
	double ali_sco;
	vector<pair<int,int> > WWW_alignment;
	Advance_Align_Dyna_Prog_Double(len1,len2,score_rec,-40,-4,-40,-4,0,0,0,0,WWW_alignment,ali_sco);
	delete [] score_rec;
	//return
	double smaller=1.0*sqrt(1.0*c1*c2);
	return 1.0*ali_sco/smaller;
}


//------------ minor_process ---------//
//given pivot and target plus the AFB, rotate target onto pivot, NOTE that pivot is fixed
void BLOSEN::Pivo_Ancho_Rotate_Single(int pivot,int target,int *In_AFB,int In_Len,double *ret_rot)
{
	int k;
	int pos1,pos2;
	int len;
	//superimpose
	len=In_Len;
	pos1=In_AFB[target];
	pos2=In_AFB[pivot];
	for(k=0;k<len;k++)BN_tmp1[k]=INPUT_MOL[target][pos1+k];
	for(k=0;k<len;k++)BN_tmp2[k]=INPUT_MOL[pivot][pos2+k];
	kabsch(BN_tmp2,BN_tmp1,len,ret_rot);
	rot_mol(INPUT_MOL[target],INPUT_MOL_Rot[target],INPUT_LEN[target],ret_rot);
}
void BLOSEN::Pivo_Ancho_Rotate(int pivot,int *In_AFB,int In_Len,double **ret_rot)
{
	int j;
	//init assign
	EqualArray(INPUT_MOL_Rot[pivot],INPUT_MOL[pivot],INPUT_LEN[pivot]);
	for(int i=0;i<12;i++)ret_rot[pivot][i]=0.0;
	ret_rot[pivot][0]=1.0;
	ret_rot[pivot][4]=1.0;
	ret_rot[pivot][8]=1.0;
	//superimpose
	for(j=0;j<INPUT_NUM;j++)
	{
		if(j==pivot)continue;
		Pivo_Ancho_Rotate_Single(pivot,j,In_AFB,In_Len,ret_rot[j]);
	}
}
//given pivot and target plus one HSFB, calculate all colored fragments
//NOTE!! pivot is fixed, while target is superimposed onto pivot NOW!!
//NOTE!! must first calculate INPUT_MOL_Rot !!
int BLOSEN::HSFB_Anchor_List_Single(int pivot,int target,int Sele_Num,int **Sele_AFB,int *Sele_Len,double *ret_rot)
{
	int l;
	int ii,jj;
	int winlen;
	int totscore;
	//for each Sele_AFB
	TMP_SFP_H.clear();
	for(l=0;l<Sele_Num;l++) 
	{
		//get SFP
		winlen=Sele_Len[l];
		ii=Sele_AFB[l][target];
		jj=Sele_AFB[l][pivot];
		//add to list
		SFP_Record sfp_rec;
		sfp_rec.score=0;
		sfp_rec.ii=ii;
		sfp_rec.jj=jj;
		sfp_rec.winlen=winlen;
		TMP_SFP_H.push_back(sfp_rec);
	}
	//get ali3-DynaProg
	ZM_Input_Mol(INPUT_MOL[target],INPUT_MOL[pivot],INPUT_LEN[target],INPUT_LEN[pivot]);
	EqualArray(rot_mat,ret_rot,12);
	totscore=ZoomIn_Add(1,INI_CUT,INI_CUT,TMP_SFP_H);
	EqualArray(ret_rot,rot_mat,12);
	return totscore;
}

//------------ process ----------//[0.Additional Function]
void BLOSEN::Other_Pivot_Function(int pivot)
{
	//SFP_generate -> WARNING!! Only with this class!!! //__110105__//
	{
		int j;
		int seed=pivot;
		for(j=0;j<INPUT_NUM;j++)
		{
			if(seed==j)
			{
				BN_SFP_H[j].clear();
				BN_SFP_L[j].clear();
				continue;
			}
			SFP_Generate_Pair(j,seed,BN_SFP_H[j],BN_SFP_L[j]);
		}
	}
}

//------------ process ----------//[1.HSFB]
void BLOSEN::HSFB_Generate_Part(int pivot,int winlen,int &Out_Num,int **Out_AFB,int *Out_Len)
{
	int i,j,k;
	int seed=pivot;
	double score;
	double totscore;
	int totcount;
	double wsmax;
	int maxnum;

	//generate
	totcount=Out_Num;
	for(i=0;i<INPUT_LEN[seed]-winlen+1;i++)
	{
		totcount++;
		totscore=0;
		for(j=0;j<INPUT_NUM;j++)
		{
			wsmax=INT_MIN_NUM;
			maxnum=-1;
			if(seed==j)continue;
			for(k=0;k<INPUT_LEN[j]-winlen+1;k++)
			{
				//this should be virtual function!!
				score=HSFB_Score_Function(seed,j,i,k,winlen); //WARNING: score range should between [INT_MIN_NUM,INT_MAX_NUM]
				if(score>wsmax)
				{
					wsmax=score;
					maxnum=k;					
				}
			}
			if(maxnum!=-1) // Has neibors
			{
				Out_AFB[totcount-1][j]=maxnum;
				totscore+=(INT_MAX_NUM+wsmax);
			}
			else
			{
				Out_AFB[totcount-1][j]=-1;
			}
		}
		if(totscore>0)  // Only record the AFB with one more neibors
		{
			Fast_Sort_Val[totcount-1]=(int)totscore;
			Out_AFB[totcount-1][seed]=i;
			Out_Len[totcount-1]=winlen;
		}
		else totcount--;
	}
	Out_Num=totcount;
}
void BLOSEN::HSFB_Shaving(int In_Num,int **In_AFB,int *In_Len,int &Out_Num,int **Out_AFB,int *Out_Len)
{
	int i,j,k;
	int index;
	int count;
	int winlen;
	int ii;
	int cover_num;
	int real_num;
	double cover_rate;


	//sort-> WARNING!! Shaving process MUST follows HSFB_Generate
	fast_sort.fast_sort_1(Fast_Sort_Val,Fast_Sort_Idx,In_Num);
	for(i=0;i<In_Num;i++)
	{
		index=Fast_Sort_Idx[i];
		EqualArray(Out_AFB[i],In_AFB[index],INPUT_NUM);
		Out_Len[i]=In_Len[index];
	}
	EqualArray2D(In_AFB,Out_AFB,In_Num,INPUT_NUM);
	EqualArray(In_Len,Out_Len,In_Num);

	//process
	count=0;
	for(i=0;i<INPUT_NUM;i++)for(j=0;j<INPUT_LEN[i];j++)Shav_Grid[i][j]=-1;
	for(i=0;i<In_Num;i++)
	{
		index=i;
		cover_num=0;
		real_num=0;		
		for(j=0;j<INPUT_NUM;j++)
		{
			ii=In_AFB[index][j];
			if(ii==-1)continue;
			winlen=In_Len[index];
			for(k=0;k<winlen;k++) 
			{
				real_num++;
				if(Shav_Grid[j][ii+k]!=-1)cover_num++;
			}
		}
		cover_rate=1.0*cover_num/real_num;
		if(cover_rate<SHAVE_THRES) // only adding when the overlap-rate is below threshold
		{
			//assign
			for(j=0;j<INPUT_NUM;j++)Out_AFB[count][j]=In_AFB[index][j];
			Out_Len[count]=In_Len[index];
			//cover
			for(j=0;j<INPUT_NUM;j++)
			{
				ii=Out_AFB[count][j];
				if(ii==-1)continue;
				winlen=Out_Len[count];
				for(k=0;k<winlen;k++)Shav_Grid[j][ii+k]=1;
			}
			//termi
			count++;
		}
	}
	Out_Num=count;
}

//------------ process ----------//[2.Pivot-select]
int BLOSEN::Repre_Mol_Select(int num,int *len,XYZ **mol,char **ami,int *ret_set)
{
	int i;
	for(i=0;i<num;i++)ret_set[i]=i;
	return num;
}
int BLOSEN::Pivot_Mol_Select(int In_Num,int **In_AFB,int *In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,int *ret_set)
{
	return 0;
}

//------------ process ----------//[3.Anchor-HSFB]
int BLOSEN::HSFB_Sort_List_Generate(int pivot,int In_Num,int **In_AFB,int *In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,int *ret_set)
{
	int i,j;
	int totscore;
	for(i=0;i<In_Num;i++) //for each In_AFB
	{
		//[1]super_impose
		Pivo_Ancho_Rotate(pivot,In_AFB[i],In_Len[i],TEMP_ROT);
		//[2]collect SFP_H
		totscore=0;
		for(j=0;j<INPUT_NUM;j++)
		{
			if(j==pivot)continue;
			totscore+=HSFB_Anchor_List_Single(pivot,j,Sele_Num,Sele_AFB,Sele_Len,TEMP_ROT[j]);
		}
		Fast_Sort_Val[i]=totscore;
	}
	//fast_sort
	fast_sort.fast_sort_1(Fast_Sort_Val,Fast_Sort_Idx,In_Num);
	EqualArray(ret_set,Fast_Sort_Idx,In_Num);
	return In_Num;
}
int BLOSEN::HSFB_Anchor_List_Generate(int pivot,int *In_AFB,int In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,double cutoff,double **ret_rot)
{
	int j;
	int minlen;
	int count;
	int totnum;
	double thres;
	//init assign
	EqualArray(INPUT_MOL_Rot[pivot],INPUT_MOL[pivot],INPUT_LEN[pivot]);
	for(int i=0;i<12;i++)ret_rot[pivot][i]=0.0;
	ret_rot[pivot][0]=1.0;
	ret_rot[pivot][4]=1.0;
	ret_rot[pivot][8]=1.0;
	Mol_Anchor[pivot]=1; //pivot
	//others process
	totnum=1;
	for(j=0;j<INPUT_NUM;j++)    //for each structure
	{
		if(j==pivot)continue;
		Pivo_Ancho_Rotate_Single(pivot,j,In_AFB,In_Len,ret_rot[j]);
		count=HSFB_Anchor_List_Single(pivot,j,Sele_Num,Sele_AFB,Sele_Len,ret_rot[j]);
		//final judge
		minlen=INPUT_LEN[j]<INPUT_LEN[pivot]?INPUT_LEN[j]:INPUT_LEN[pivot];
		thres=1.0*count/minlen;
		if(thres>cutoff)
		{
			Mol_Anchor[j]=1; //annchored
			totnum++;
		}
		else
		{
			Mol_Anchor[j]=0; //unannchored
		}
	}
	return totnum;
}

//------------ process ----------//[4.pariwise refine]
//-> target is first, pivot is second!!
double BLOSEN::Pairwise_Refine(int p1,int p2,double *ini_rot,int *ali2,double *ret_rot)
{
	double tms=CLeFAPS_Part(INPUT_MOL[p1],INPUT_MOL[p2],INPUT_LEN[p1],INPUT_LEN[p2],ini_rot,BN_SFP_L[p1],ali2);
	EqualArray(ret_rot,ZM_CURMAT,12);
	return tms;
}
double BLOSEN::Pairwise_Total(int p1,int p2,int *ali2_comp,int *ali2,double *ret_rot)
{
	double tms=CLeFAPS_Full(INPUT_MOL[p1],INPUT_MOL[p2],INPUT_LEN[p1],INPUT_LEN[p2],BN_SFP_H[p1],BN_SFP_L[p1],ali2,ali2_comp);
	EqualArray(ret_rot,ZM_FINMAT,12);
	return tms;
}
//-> given pivot_ali and cutoff, return medium aligned set
double BLOSEN::Refine_HSFB_Generate_Part(int pivot,int **pivot_ali,double cutoff,int *ali_out)
{
	int i,j;
	int count;
	int totnum;
	int totlen=0;
	double thres;
	int len=INPUT_LEN[pivot];
	for(i=0;i<len;i++)
	{
		count=0;
		totnum=0;
		for(j=0;j<INPUT_NUM;j++)
		{
			if(Mol_Anchor[j]==0)continue;
			totnum++;
			if(pivot_ali[j][i]!=-1)count++;
		}
		thres=1.0*count/totnum;
		if(thres+0.01>cutoff)
		{
			ali_out[i]=1; //colored
			totlen++;
		}
		else
		{
			ali_out[i]=-1; //uncolored
		}
	}
	return totlen;
}
double BLOSEN::Deal_Anchored(int pivot,int target,int color_cut,
							 int In_Num,int **In_AFB,int *In_Len,int Sele_Num,int **Sele_AFB,int *Sele_Len,
							 int *pivot_ali,int *ret_ali,double *ret_rot)
{
	int i;
	double score;
	double wsmax;
	int match;
	int wsmax_match;
	double tms;
	double wsmax_tms;
	double rotmat[12];
	int pivot_len=INPUT_LEN[pivot];

	//init_score
	tms=TEMP_TMS[target];
	match=Comparison_Alignment(pivot_ali,ret_ali,pivot_len);
	wsmax=(tms+0.01)*(match+1);
	wsmax_tms=tms;
	wsmax_match=match;
	//init_alignment
	EqualArray(Check_ali2,ret_ali,pivot_len);
	Ali2_To_Ali1(INPUT_LEN[target],INPUT_LEN[pivot],Check_ali1,Check_ali2);
	//get higher_score
	int valid;
	int count=0;
	for(i=0;i<In_Num;i++)
	{
		//head_check
		if(count>=color_cut)break;
		if(HEADCKorNOT==1)
		{
			valid=Possessed_Check(In_AFB[i][target],In_AFB[i][pivot],In_Len[i],Check_ali1,Check_ali2,
				INPUT_LEN[target],INPUT_LEN[pivot]);
			if(valid==0)continue;
		}
		//pairwise_refine
		Pivo_Ancho_Rotate_Single(pivot,target,In_AFB[i],In_Len[i],rotmat);
		HSFB_Anchor_List_Single(pivot,target,Sele_Num,Sele_AFB,Sele_Len,rotmat);
		tms=Pairwise_Refine(target,pivot,rotmat,Pivot_Ali_TTT,rotmat);
		//judge
		match=Comparison_Alignment(pivot_ali,Pivot_Ali_TTT,pivot_len);
		score=(tms+0.01)*(match+1);
		if(score>wsmax)
		{
			wsmax=score;
			wsmax_tms=tms;
			wsmax_match=match;
			EqualArray(ret_rot,rotmat,12);
			EqualArray(ret_ali,Pivot_Ali_TTT,pivot_len);
		}
		count++;
	}
	if(wsmax_tms>0.99)wsmax_tms=0.99;
	return wsmax_match+wsmax_tms;
}

//------------ process ----------//[5.consensus refinement]
double BLOSEN::Terminal_Refine_Main(int pivot,double **in_rot,XYZ **in,int *len,int totnum,int **alin,int totin,int **aliout,int &totout,double **ret_rot,int &core_len)
{
	int i;
	for(i=0;i<totnum;i++)rot_mol(in[i],INPUT_MOL_Rot[i],len[i],in_rot[i]);
	totout=Single_Ali_To_Multi_Ali(totnum,len,totin,alin,aliout);
	BC_Update_Partial_CORE(pivot,INPUT_MOL_Rot,len,totnum,aliout,totout,aliout,totout);
	double sumtms=BC_Tota_Update_Main(INPUT_MOL_Rot,len,totnum,aliout,totout,aliout,totout);
	core_len=BC_Return_Conserved_Core(totnum,totout,aliout);
	for(i=0;i<totnum;i++)
	{
		double rmsd;
		rmsd=kabsch(BC_Best_Output[i],in[i],len[i],ret_rot[i]);
	}
	return sumtms;
}



