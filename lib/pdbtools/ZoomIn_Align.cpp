#include "ZoomIn_Align.h"

//------ constructor -------//
ZoomIn_Align::ZoomIn_Align(int num)
:Ali_Ori(num),TM_align(num),
Ali_Ali3(num),Ali_AFP(num)
{
	//init
	ZoomIn_Maximal=num;
	//macros
	TMSorNOT=1;    //apply TM_score dynamic programming
	HEAD_Chk=0;    //apply head check (Full_CLEPAPS only!!)
	//parameter
	ZOOM_ITER=3;
	REFINE_ITER=10;
	ANG_CUT=60;
	FIN_CUT=5.0;   //this can be self-adapted
	INI_CUT=10.0;  //this can be self-adapted
	//TopK & TopJ
	ZM_TopK=20;    //only consider 20
	ZM_TopJ=num;   //using total SFP_H (this should be less than Maximal!!!)
	//create
	Init_ZoomIn_Align(ZoomIn_Maximal);
}
ZoomIn_Align::~ZoomIn_Align(void)
{
	//delete
	Dele_ZoomIn_Align();
}

//--------- init ---------//
void ZoomIn_Align::Init_ZoomIn_Align(int maxlen)
{
	//-- ali3_init --//
	ali1=new int[maxlen];
	ali2=new int[maxlen];
	mol1=new XYZ[maxlen];
	mol2=new XYZ[maxlen];
	//-- head_index --//
	center_mat=new double[12*maxlen];
	head_score=new int[maxlen];
	head_index=new int[maxlen];
	//-- AFP & ali ---//
	tali1=new int[maxlen];
	tali2=new int[maxlen];
	AFP_Cor=new int[4*maxlen];
	AFP_Cor_temp=new int[4*maxlen];
	//-- rot_mat --//
	rot_mat=new double[12];
	rot_bak=new double[12];
	//-- best record --//
	ZM_BEST_ALI=new int[maxlen];
	ZM_CURMAT=new double[12];
	ZM_FINMAT=new double[12];
}
void ZoomIn_Align::Dele_ZoomIn_Align(void)
{
	//-- ali3_dele --//
	delete [] ali1;
	delete [] ali2;
	delete [] mol1;
	delete [] mol2;
	//-- head_index --//
	delete [] center_mat;
	delete [] head_score;
	delete [] head_index;
	//-- AFP & ali ---//
	delete [] tali1;
	delete [] tali2;
	delete [] AFP_Cor;
	delete [] AFP_Cor_temp;
	//-- rot_mat --//
	delete [] rot_mat;
	delete [] rot_bak;
	//-- best record --//
	delete [] ZM_BEST_ALI;
	delete [] ZM_CURMAT;
	delete [] ZM_FINMAT;
}
void ZoomIn_Align::ZM_Input_Mol(XYZ *m1,XYZ *m2,int n1,int n2)
{
	EqualArray(mol1,m1,n1);
	EqualArray(mol2,m2,n2);
	moln1=n1;
	moln2=n2;
	CLeFAPS_Init();
}
void ZoomIn_Align::CLeFAPS_Init(void)
{
	//STEP[0] initialization
	int mol_len=moln1<=moln2?moln1:moln2;   //smaller_length
	double d_0=Calc_TM_d0(mol_len);
	//self_adaptive[1]->RMSD_CUTOFF
	INI_CUT=3*d_0;
	INI_CUT=INI_CUT<MAX_DIST?INI_CUT:MAX_DIST;
	INI_CUT=INI_CUT>MIN_DIST?INI_CUT:MIN_DIST;
	FIN_CUT=d_0;
	FIN_CUT=FIN_CUT<MAX_DIST?FIN_CUT:MAX_DIST;
	FIN_CUT=FIN_CUT>MIN_DIST?FIN_CUT:MIN_DIST;
	//ali_init
	Ali_Init(moln1,moln2,ali1,ali2);
}

//-------- process_minor ---------//
int ZoomIn_Align::ZoomIn_Add(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP)
{
	int i,j;
	int index;
	int ii,jj;
	int winlen;
	int QQ_temp_num;
	int QQ_cur;
	int Mx_K;
	double ret_val;
	double thres,minus;
	int neo_ii,neo_jj,neo_len;
	double score;
	double RMSD;
	double d_0=d0;
	SFP_Record sfp_rec;

	//tot_init
	QQ_cur=0;
	QQ_temp_num=(int)SFP.size();
	valid_idx.resize(QQ_temp_num);
	for(i=0;i<QQ_temp_num;i++)valid_idx[i]=1;
	//rigid_frag_adding
	minus=(INI_CUT-FIN_CUT)/recur;
	for(i=0;i<recur;i++)
	{
		//init
		Ali3_Init(ali3,moln2);
		if(i<recur-1)Mx_K=(int)(QQ_temp_num/pow(2.0,i+1));
		else Mx_K=(int)(QQ_temp_num/pow(2.0,i));
		thres=INI_CUT-(i+1)*minus;
		//check_frag
		for(j=0;j<QQ_cur+Mx_K;j++)
		{
			if(valid_idx[j]==-1)continue;  // -1 means killed
			index=j;
			sfp_rec=SFP.at(index);
			ii=sfp_rec.ii;
			jj=sfp_rec.jj;
			winlen=sfp_rec.winlen;
			calc_frag_dist_thres(ii,jj,winlen,rot_mat,mol1,mol2,moln1,moln2,thres,neo_ii,neo_jj,neo_len,score);
			if(neo_len==0)valid_idx[j]=-1;
			else Ali3_Add(neo_ii,neo_jj,neo_len,rot_mat,mol1,mol2,moln1,moln2,ali3);
		}
		QQ_cur+=Mx_K;
		//rot_mol
		ret_val=ali3_kabsch(rot_mat,RMSD,mol1,mol2,moln1,moln2,ali3); //apply ALI3(weight) rotation
		if(ret_val<0.0 || RMSD<0.0)return -1; //incorrect
	}
	//Ali3_DynaProg
	Ali3_DynaProg(rot_mat,0.0,d_0,mol1,mol2,moln1,moln2,FOR_GAP,BAK_GAP,SCALE,ali3);
	return (int)Ali3_TraceBack(rot_mat,d_0,mol1,mol2,moln1,moln2,ali3,ali1,ali2);
}
int ZoomIn_Align::Refinement(int recur,int Range,int ori_sco,double FIN_CUT,double DP_BETA)
{
	int k;
	double ret_val;
	int score,bak_score;
	int success;
	int count;
	double RMSD;
	double d_0=d0;

	//process
	count=0;
	success=1; //default:success
	bak_score=ori_sco;
	score=ori_sco;
	for(k=0;k<recur;k++)
	{
		//Back
		if(success==1)
		{
			EqualArray(rot_bak,rot_mat,12);
			EqualArray(tali1,ali1,moln1);
			EqualArray(tali2,ali2,moln2);
		}
		if(count==2)return score; //[refine_better]//__081130__//
		//Elongation & Collect
		Elongation(rot_mat,FIN_CUT,INT_MAX_NUM,ali1,ali2);
		Ali_To_Cor(AFP_Cor,1,moln1,moln2,ali1,ali2);
		//Rotation	
		Ali3_Cor_Add(AFP_Cor,0,rot_mat,mol1,mol2,moln1,moln2,ali3);
		ret_val=ali3_kabsch(rot_mat,RMSD,mol1,mol2,moln1,moln2,ali3); //apply ALI3(weight) rotation
		if(ret_val<0.0 || RMSD<0.0)return -1; //incorrect
		//Replace
		Ali3_Cor_Add(AFP_Cor,Range,rot_mat,mol1,mol2,moln1,moln2,ali3);
		Ali3_DynaProg(rot_mat,DP_BETA,d_0,mol1,mol2,moln1,moln2,FOR_GAP,BAK_GAP,SCALE,ali3);
		score=(int)Ali3_TraceBack(rot_mat,d_0,mol1,mol2,moln1,moln2,ali3,ali1,ali2);
		//Check
		if(score>bak_score)
		{
			bak_score=score;
			success=1;
			count=0;
		}
		else
		{
			if(score<bak_score*0.95)
			{
				EqualArray(rot_mat,rot_bak,12);
				EqualArray(ali1,tali1,moln1);
				EqualArray(ali2,tali2,moln2);
				return (bak_score);
			}
			success=0;
			count++;
		}
	}
	return score; //[refine_better]//__081130__//
}
void ZoomIn_Align::Elongation(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2)
{
	int i;
	int ii,jj;
	int index1,index2;
	double Dou_Thres=thres*thres;
	int success;
	XYZ temp1,temp2,temp;
	int cut1,cut2,cut;

	//forward_elongation
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]<0)continue;
		//forward
		index1 = ali2[i];
		index2 = i;

#ifdef DEBUG
	overrange_debug_test(index1,moln1);
	overrange_debug_test(index2,moln2);
	lessrange_debug_test(index1,0);
	lessrange_debug_test(index2,0);
#endif

		for(;;)
		{
			success=0; //default:fail
			index1++;
			index2++;
			if(index1>moln1-1 || index2>moln2-1) break;
			jj=ali1[index1];
			ii=ali2[index2];
			if(jj!=-1 && ii!=-1) break;
			else if(jj!=-1 && ii==-1)
			{
				rot_point(mol1[index1],temp,rotmat_);
				cut1=Ali3_Vector_Score_Forward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				cut2=Ali3_Vector_Score_Forward(index1,jj,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if((mol2[index2].distance_square(temp) < mol2[jj].distance_square(temp)) || cut1>cut2)
				{
					ali2[jj]=-1;
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}
			else if(jj==-1 && ii!=-1)
			{
				rot_point(mol1[index1],temp1,rotmat_);
				rot_point(mol1[ii],temp2,rotmat_);
				cut1=Ali3_Vector_Score_Forward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				cut2=Ali3_Vector_Score_Forward(ii,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if((mol2[index2].distance_square(temp1) < mol2[index2].distance_square(temp2)) || cut1>cut2)
				{
					ali1[ii]=-1;
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}
			else
			{
				rot_point(mol1[index1],temp,rotmat_);
				cut=Ali3_Vector_Score_Forward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if(mol2[index2].distance_square(temp) < Dou_Thres || cut>cutoff)
				{
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}

			if(success==0)break; //failed
			i++;
		}//end of FOR(;;)
	}//end of FOR(i)


	//backward_elongation
	for(i=moln2-1;i>=0;i--)
	{
		if(ali2[i]<0)continue;
		//forward
		index1 = ali2[i];
		index2 = i;

#ifdef DEBUG
	overrange_debug_test(index1,moln1);
	overrange_debug_test(index2,moln2);
	lessrange_debug_test(index1,0);
	lessrange_debug_test(index2,0);
#endif

		for(;;)
		{
			success=0; //default:fail
			index1--;
			index2--;
			if(index1<0 || index2<0) break;
			jj=ali1[index1];
			ii=ali2[index2];
			if(jj!=-1 && ii!=-1) break;
			else if(jj!=-1 && ii==-1)
			{
				rot_point(mol1[index1],temp,rotmat_);
				cut1=Ali3_Vector_Score_Backward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				cut2=Ali3_Vector_Score_Backward(index1,jj,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if((mol2[index2].distance_square(temp) < mol2[jj].distance_square(temp)) || cut1>cut2)
				{
					ali2[jj]=-1;
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}
			else if(jj==-1 && ii!=-1)
			{
				rot_point(mol1[index1],temp1,rotmat_);
				rot_point(mol1[ii],temp2,rotmat_);
				cut1=Ali3_Vector_Score_Backward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				cut2=Ali3_Vector_Score_Backward(ii,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if((mol2[index2].distance_square(temp1) < mol2[index2].distance_square(temp2)) || cut1>cut2)
				{
					ali1[ii]=-1;
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}
			else
			{
				rot_point(mol1[index1],temp,rotmat_);
				cut=Ali3_Vector_Score_Backward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if(mol2[index2].distance_square(temp) < Dou_Thres || cut>cutoff)
				{
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}

			if(success==0)break; //failed
			i--;
		}//end of FOR(;;)
	}//end of FOR(i)
}
void ZoomIn_Align::Kill_Bad(double thres,double cutoff,double *rotmat_)
{
	int i;
	int ii,jj;
	double cut;
	double dist;
	XYZ temp;
	double Dou_Thres=thres*thres;
	double SQURE_DIST=MAX_DIST*MAX_DIST;

	for(i=0;i<moln2;i++)
	{
		ii=ali2[i];
		jj=i;
		if(ii>=0)
		{
			rot_point(mol1[ii],temp,rotmat_);
			dist=mol2[jj].distance_square(temp);
			if(dist>SQURE_DIST)
			{
				ali2[jj]=-1;
				ali1[ii]=-1;
			}
			else if(dist>Dou_Thres)
			{
				cut=Ali3_Vector_Score_Forward(ii,jj,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if(cut<cutoff)
				{
					ali2[jj]=-1;
					ali1[ii]=-1;
				}
			}
		}
	}
}
void ZoomIn_Align::Final_Check(double FIN_CUT,double ANG_CUT)
{
	//--Final_Statis--//
	//ws statis
	double RMSD;
	Elongation(rot_mat,FIN_CUT,INT_MAX_NUM,ali1,ali2);
	Kill_Bad(FIN_CUT,ANG_CUT,rot_mat);
	Ali_To_Cor(AFP_Cor,1,moln1,moln2,ali1,ali2);
	AFP_kabsch(rot_mat,RMSD,mol1,mol2,moln1,moln2,AFP_Cor);      // the only step to generate final rot_mat
}

//---------- process full -------------//
//given a certain SFP, return its' TM_score against other SFPs
int ZoomIn_Align::Select_Single(int cur_SFP,int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP)
{
	int ii,jj,winlen;
	int wsscore;
	double RMSD;
	SFP_Record sfp_rec;
	//init
	sfp_rec=SFP.at(cur_SFP);
	ii=sfp_rec.ii;
	jj=sfp_rec.jj;
	winlen=sfp_rec.winlen;
	//pivot
	if(make_center(ii,jj,winlen,rot_mat,mol1,mol2,moln1,moln2)<0.0)return -1;
	if(ZoomIn_Add(recur,INI_CUT,FIN_CUT,SFP)<0)return -1;
	AFP_Kill_NonLinear(AFP_Cor,AFP_Cor_temp,rot_mat);
	AFP_kabsch(rot_mat,RMSD,mol1,mol2,moln1,moln2,AFP_Cor);
	wsscore=(int)T_Similar(rot_mat,ali1,ali2);
	return wsscore;
}
//given a SFP_list, return its sorted list on the spatial-consistency
int ZoomIn_Align::Select_Best_Pivot(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,int TopJ)
{
	int i;
	int relnum=(int)SFP.size();
	int totnum=TopJ<relnum?TopJ:relnum;
	//real_start
	for(i=0;i<totnum;i++)
	{
		head_score[i]=Select_Single(i,recur,INI_CUT,FIN_CUT,SFP);
		EqualArray(center_mat+12*i,rot_mat,12);
	}
	//fast_sort
	fast_sort.fast_sort_1(head_score,head_index,totnum);
	return totnum;
}
//calculate simple TMscore
double ZoomIn_Align::T_Similar(double *rotmat_,int *ali1,int *ali2)
{
	int i;
	int ii,jj;
	double dist2;
	double ori_d;
	double score;
	double d_0=d0;
	XYZ temp;

	score=0.0;
	temp=0.0;
	ori_d=d_0*d_0;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]==-1)continue;
		if(ali2[i]<0)ii=ali2[i]+INT_MAX_NUM;
		else ii=ali2[i];
		jj=i;
		rot_point(mol1[ii],temp,rotmat_);
		dist2=mol2[jj].distance_square(temp);
		score+=1.0/(1.0+dist2/ori_d);
	}
	return score;
}
//given ali2_current and ali2_compasiron, return match length
int ZoomIn_Align::Comparison_Alignment(int *ali2_cur,int *ali2_comp,int length)
{
	int i;
	int match=0;
	for(i=0;i<length;i++)
	{
		if(ali2_comp[i]>=0&&ali2_cur[i]>=0)match++;
	}
	return match;
}
//given ali1 and ali2, check the SFP is possessed or not
int ZoomIn_Align::Possessed_Check(int ii,int jj,int winlen,int *ali1,int *ali2,int moln1,int moln2)
{
#ifdef DEBUG
	overrange_debug_test(ii+winlen-1,moln1);
	overrange_debug_test(jj+winlen-1,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif
	int k;
	for(k=0;k<winlen;k++)
	{
		if(ali1[ii+k]!=-1&&ali2[jj+k]!=-1)
		{
			if( (ali1[ii+k]!=jj+k) || (ali2[jj+k]!=ii+k))return 1;
		}
	}
	return 0;
}
int ZoomIn_Align::Possessed_Check(int ii,int jj,int winlen,vector <int> &alignment,int moln1,int moln2)
{
#ifdef DEBUG
	overrange_debug_test(ii+winlen-1,moln1);
	overrange_debug_test(jj+winlen-1,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif
	int k;
	for(k=0;k<winlen;k++)
	{
		if(alignment.at((ii+k)*moln2+(jj+k))==0)return 1; //un-possessed!!
	}
	return 0;
}

//==================== main process ====================//-> part ZoomIn
//[input]
//m1 must be initially superimposed onto m2,
//SFP_L -> Similar_Fragment_Pair low-valued list
//[output]
//ali2 -> the final returned alignment
double ZoomIn_Align::CLeFAPS_Part(XYZ *m1,XYZ *m2,int n1,int n2,double *rotmat_,vector <SFP_Record> &SFP_Low,int *ali2_)
{
	//judge
	int i;
	int smaller=n1<n2?n1:n2;
	if(smaller==0)
	{
		for(i=0;i<n2;i++)ali2_[i]=-1;
		return 0.0;
	}
	//init
	ZM_Input_Mol(m1,m2,n1,n2);
	if(rotmat_==0)
	{
		for(int i=0;i<12;i++)rot_mat[i]=0.0;
		rot_mat[0]=1.0;
		rot_mat[4]=1.0;
		rot_mat[8]=1.0;
	}
	else EqualArray(rot_mat,rotmat_,12);
	//ZoomIn stragety
	int ret_val;
	ret_val=ZoomIn_Add(ZOOM_ITER,INI_CUT,FIN_CUT,SFP_Low);
	if(ret_val==-1)return -1;
	ret_val=Refinement(REFINE_ITER,4,ret_val,FIN_CUT,0.5);
	if(ret_val==-1)return -1;
	AFP_Kill_NonLinear(AFP_Cor,AFP_Cor_temp,rot_mat);
	Final_Check(FIN_CUT,ANG_CUT);
	Elongation(rot_mat,FIN_CUT,ANG_CUT,ali1,ali2);
	//TM_align final
	double TM_score;
	double RMSD;
	int LALI;
	if(TMSorNOT==1)
	{
		//TM_init
		TM_Align_Init(n1,n2);
//		TM_GAP_TYPE=0;  //normal one-layer DynaProg
//		TM_GAP_STAGE=2; //normal TMalign
		TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
		Calc_TM_Align(m1,m2,n1,n2,ali2,ali2,0,3);    //apply TMscore dynamic programming
		TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
	}
	TM_score=TM_Align_TM_Score_Simp(m1,m2,n1,n2,ali2,RMSD,LALI); //use simple version
	EqualArray(ZM_CURMAT,finmat,12);
	EqualArray(ali2_,ali2,n2);
	return TM_score;
}

//==================== main process ====================//-> full ZoomIn
//[input]
//SFP_H -> Similar_Fragment_Pair high-valued list
//SFP_L -> Similar_Fragment_Pair low-valued list
//[output]
//ali2 -> the final returned alignment
double ZoomIn_Align::CLeFAPS_Full(XYZ *m1,XYZ *m2,int n1,int n2,vector <SFP_Record> &SFP_High,vector <SFP_Record> &SFP_Low,
									 int *ali2_,int *ali2_comp)
{
	int i,j;
	int ii,jj,winlen;
	int retv;
	int index;
	int totnum;
	double TM_CUR;
	double TM_MAX=-1.0;
	int match=0;
	int max_match=0;
	double score;
	double wsmax=-1.0;
	SFP_Record sfp_rec;
	//judge
	int smaller=n1<n2?n1:n2;
	if(smaller==0)
	{
		for(i=0;i<n2;i++)ali2_[i]=-1;
		return 0.0;
	}
	//init
	ZM_Input_Mol(m1,m2,n1,n2);
	totnum=Select_Best_Pivot(1,INI_CUT,INI_CUT,SFP_High,ZM_TopJ);
	if(HEAD_Chk==1)
	{
		align_rec.resize(n1*n2); //not occupied
		for(i=0;i<n1*n2;i++)align_rec[i]=0;
	}
	//start
	int count=0;
	for(i=0;i<totnum;i++)
	{
		//head_check
		if(count>=ZM_TopK)break;
		index=head_index[i];
		if(HEAD_Chk==1)
		{
			sfp_rec=SFP_High.at(index);
			ii=sfp_rec.ii;
			jj=sfp_rec.jj;
			winlen=sfp_rec.winlen;
			retv=Possessed_Check(ii,jj,winlen,align_rec,n1,n2);
			if(retv==0)continue;
		}
		//CLEFAPS_Part
		EqualArray(rot_mat,center_mat+12*index,12);
		TM_CUR=CLeFAPS_Part(m1,m2,n1,n2,rot_mat,SFP_Low,ZM_BEST_ALI);
		score=TM_CUR;
		if(HEAD_Chk==1)
		{
			for(j=0;j<n2;j++)
			{
				ii=ZM_BEST_ALI[j];
				jj=j;
				if(ii!=-1)align_rec[ii*n2+jj]=1;
			}
		}
		//comparison
		if(ali2_comp!=0)
		{
			match=Comparison_Alignment(ali2,ali2_comp,n2);
			score=(TM_CUR+0.01)*(match+1);
		}
		//comparison_over
		if(score>wsmax)
		{
			wsmax=score;
			max_match=match;
			TM_MAX=TM_CUR;
			EqualArray(ali2_,ZM_BEST_ALI,n2);
			EqualArray(ZM_FINMAT,ZM_CURMAT,12);
		}
		count++;
	}
	return TM_MAX+max_match;
}
