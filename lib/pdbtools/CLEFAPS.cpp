#include "CLEFAPS.h"

//------ constructor -------//
CLEFAPS::CLEFAPS(int num)
:Ali_Ori(num),Ali_AFP(num),TM_align(num),
 Ali_Ali3(num)
{
	//init
	CLEF_Maximal=num;
	CLEFAPS_Init(num);
	CUR_MaxJ_thres=0;
	//parameters
	ANG_CUT=60;         //
	Tail_Check=0.9;     //whether do tail check
	Comp_Check=0.9;     //whether do comparison check
	Comp_Check_Vari=0;  //comparison check (partial count)
}
CLEFAPS::~CLEFAPS(void)
{
	//dele
	CLEFAPS_Dele();
}

//=================== process functions ============//
//--------- init ---------//
void CLEFAPS::CLEFAPS_Init(int maxlen)
{
	//-- ali3_init --//
	ali1=new int[maxlen];
	ali2=new int[maxlen];
	tali1=new int[maxlen];
	tali2=new int[maxlen];
	//-- head_index --//
	valid_idx=new short[maxlen*maxlen];
	center_mat=new double[12*maxlen*10];
	head_score=new int[maxlen*10];
	head_index=new int[maxlen*10];
	//-- AFP & ali ---//
	AFP_Cor=new int[4*maxlen];
	AFP_Cor_temp=new int[4*maxlen];
	//-- CLEP_rotmat --//
	CLEP_rotmat=new double[12];
}
void CLEFAPS::CLEFAPS_Dele(void)
{
	//-- ali3_dele --//
	delete [] ali1;
	delete [] ali2;
	delete [] tali1;
	delete [] tali2;
	//-- head_index --//
	delete [] center_mat;
	delete [] head_score;
	delete [] head_index;
	//-- AFP & ali ---//
	delete [] AFP_Cor;
	delete [] AFP_Cor_temp;
	//-- CLEP_rotmat --//
	delete [] CLEP_rotmat;
}

//--------- vice process ---------//
//calculate simple TMscore
double CLEFAPS::T_Similar(double *rotmat_,int *ali1,int *ali2,
						  XYZ *mol1,XYZ *mol2,int moln1,int moln2)
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
		ii=ali2[i];
		jj=i;
		if(ii==-1)continue;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

		Ali_Cache_Point(mol1,ii,temp,rotmat_);
		dist2=mol2[jj].distance_square(temp);
		score+=1.0/(1.0+dist2/ori_d);
	}
	return score;
}
//----------- Elongation and Shrinking -----------//
void CLEFAPS::Elongation(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
						 XYZ *mol1,XYZ *mol2,int moln1,int moln2)
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
				Ali_Cache_Point(mol1,index1,temp,rotmat_);
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
				Ali_Cache_Point(mol1,index1,temp1,rotmat_);
				Ali_Cache_Point(mol1,ii,temp2,rotmat_);
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
				Ali_Cache_Point(mol1,index1,temp,rotmat_);
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
				Ali_Cache_Point(mol1,index1,temp,rotmat_);
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
				Ali_Cache_Point(mol1,index1,temp1,rotmat_);
				Ali_Cache_Point(mol1,ii,temp2,rotmat_);
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
				Ali_Cache_Point(mol1,index1,temp,rotmat_);
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
void CLEFAPS::Shrinking(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
						XYZ *mol1,XYZ *mol2,int moln1,int moln2)
{
	int i;
	int ii,jj;
	int cut;
	double dist;
	XYZ temp;
	double Dou_Thres=thres*thres;
	double SQURE_DIST=MAX_DIST*MAX_DIST;

	for(i=0;i<moln2;i++)
	{
		ii=ali2[i];
		jj=i;
		if(ii==-1)continue;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

		Ali_Cache_Point(mol1,ii,temp,rotmat_);
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
//------------- Select Best Pivot --------------//
//given a certain SFP, return its' TM_score against other SFPs
int CLEFAPS::Select_Single(int cur_SFP,int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,double *ret_rot)
{
	int ii,jj,winlen;
	double RMSD;
	int wsscore;
	//init
	ii=SFP[cur_SFP].ii;
	jj=SFP[cur_SFP].jj;
	winlen=SFP[cur_SFP].winlen;
	//pivot
	if(make_center(ii,jj,winlen,ret_rot,mol1,mol2,moln1,moln2)<0.0)return -1;
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	if(ZoomIn_Add(recur,INI_CUT,FIN_CUT,SFP,ret_rot)<0)return -1;
	if(NonLinear==1)Kill_NonLinear(ret_rot);
	Final_Rot_Ali(ret_rot,RMSD,mol1,mol2,moln1,moln2,ali1,ali2);
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	wsscore=(int)T_Similar(ret_rot,ali1,ali2,mol1,mol2,moln1,moln2);
	return wsscore;
}
//----------- AlignmentMetricAccuracy ----------------//
double CLEFAPS::Alg_Comp_Simp(vector <int> &ali_m1,vector <int> &ali_m2,int part) //using ali2
{
	int i;
	int ii1,ii2;
	int count=0;
	int parti=0;
	int total=0;
	for(i=0;i<moln2;i++)
	{
		ii1=ali_m1[i];
		ii2=ali_m2[i];
		if(ii1==-1 || ii2==-1)continue;
		if(ii1==ii2)count++;
		else
		{
			if(abs(ii1-ii2)<=part)parti++;
		}
		total++;
	}
	if(total==0)return 0.0;
	else return 1.0*(count+parti)/total;
}

//==================== virtual function ====================//-> part[0]: io
void CLEFAPS::CLEPAPS_Input(XYZ *m1,XYZ *m2,int n1,int n2)
{
	mol1=m1;
	mol2=m2;
	moln1=n1;
	moln2=n2;
	Other_CLEFAPS_Init();
}
//==================== virtual function ====================//-> part[1]: score
//given alignment and rotmat, return scores (maybe TMscore,RMSD and LALI, etc)
double CLEFAPS::CLEPAPS_Make_Score(Align_Record & align_record)
{
	int i;
	int ii,jj;
	double dist2;
	double ori_d;
	double sco;
	double rms;
	double tms;
	int count;
	int lali;
	double rmsd;
	//init equal
	align_record.alignment.resize(moln2);
	for(i=0;i<moln2;i++)align_record.alignment[i]=ali2[i];
	align_record.rotmat.resize(12);
	for(i=0;i<12;i++)align_record.rotmat[i]=CLEP_rotmat[i];
	//calculate
	double d_0=TM_d0;
	XYZ temp=0.0;
	ori_d=d_0*d_0;
	sco=0.0;
	rms=0.0;
	count=0;
	for(i=0;i<moln2;i++)
	{
		ii=ali2[i];
		jj=i;
		if(ii==-1)continue;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

		rot_point(mol1[ii],temp,CLEP_rotmat);
		dist2=mol2[jj].distance_square(temp);
		sco+=1.0/(1.0+dist2/ori_d);
		rms+=dist2;
		count++;
	}
	//evaluate
	rmsd=sqrt(1.0*rms/count);
	lali=count;
	tms=1.0*sco/TM_smaller;
	//final
	align_record.RMSD=rmsd;
	align_record.lali=lali;
	align_record.TMsco=tms;
	align_record.main_sco=tms;
	return tms;
}
//given the align_record, make the secondary score (maybe LALI*TMscore, etc)
double CLEFAPS::CLEPAPS_Get_Score(Align_Record &align)
{
	double tms=align.main_sco;
	return tms;
}

//==================== virtual function ====================//-> part[2]: partial
//0. Init Other CLEFAPS
void CLEFAPS::Other_CLEFAPS_Init(void)
{
	//STEP[0] initialization
	int mol_len=moln1<=moln2?moln1:moln2;   //smaller_length
	double d_0=Calc_TM_d0_Simp(mol_len);
	TM_smaller=mol_len;
	TM_d0=d_0;
	//self_adaptive[1]->RMSD_CUTOFF
	INI_CUT=3*d_0;
	INI_CUT=INI_CUT<MAX_DIST?INI_CUT:MAX_DIST;
	INI_CUT=INI_CUT>MIN_DIST?INI_CUT:MIN_DIST;
	FIN_CUT=d_0;
	FIN_CUT=FIN_CUT<MAX_DIST?FIN_CUT:MAX_DIST;
	FIN_CUT=FIN_CUT>MIN_DIST?FIN_CUT:MIN_DIST;
	//init Maximal
	CUR_MaxJ=0.0;
	CUR_MaxK=0.0;
}
//1. ZoomIn Add (ret_rot might be updated)
int CLEFAPS::ZoomIn_Add(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,double *ret_rot)
{
	int i,j;
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
	//tot_init
	QQ_cur=0;
	QQ_temp_num=(int)SFP.size();
	memset(valid_idx,1,QQ_temp_num*sizeof(short));
	//rigid_frag_adding
	minus=(INI_CUT-FIN_CUT)/recur;
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
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
			ii=SFP[j].ii;
			jj=SFP[j].jj;
			winlen=SFP[j].winlen;
			calc_frag_dist_thres(ii,jj,winlen,ret_rot,mol1,mol2,moln1,moln2,thres,neo_ii,neo_jj,neo_len,score);
			if(neo_len==0)valid_idx[j]=-1;
			else Ali3_Add(neo_ii,neo_jj,neo_len,ret_rot,mol1,mol2,moln1,moln2,ali3);
		}
		QQ_cur+=Mx_K;
		//rot_mol
		ret_val=ali3_kabsch(ret_rot,RMSD,mol1,mol2,moln1,moln2,ali3); //apply ALI3(weight) rotation
		if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
		if(ret_val<0.0 || RMSD<0.0)return -1; //incorrect
	}
	//Ali3_DynaProg
	Ali3_DynaProg(ret_rot,0.0,d_0,mol1,mol2,moln1,moln2,FOR_GAP,BAK_GAP,SCALE,ali3);
	return (int)Ali3_TraceBack(ret_rot,d_0,mol1,mol2,moln1,moln2,ali3,ali1,ali2);
}
//2 .Refinement (ret_rot might be updated)
int CLEFAPS::Refinement(int recur,double FIN_CUT,double *ret_rot)
{
	int k;
	double ret_val;
	int score,bak_score;
	int success;
	int count;
	double rot_bak[12];
	double RMSD;
	int Range=4;
	double DP_BETA=0.5;
	double d_0=d0;

	//process
	count=0;
	success=1; //default:success
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	bak_score=(int)Ali3_TraceBack(ret_rot,d_0,mol1,mol2,moln1,moln2,ali3,ali1,ali2);
	score=bak_score;
	for(k=0;k<recur;k++)
	{
		//Back
		if(success==1)
		{
			memcpy(rot_bak,ret_rot,12*sizeof(double));
			memcpy(tali1,ali1,moln1*sizeof(int));
			memcpy(tali2,ali2,moln2*sizeof(int));
		}
		if(count==2)return score; //[refine_better]//__081130__//
		//Elongation & Collect
		Elongation(ret_rot,FIN_CUT,INT_MAX_NUM,ali1,ali2,mol1,mol2,moln1,moln2);
		Ali_To_Cor(AFP_Cor,1,moln1,moln2,ali1,ali2);
		//Rotation	
		Ali3_Cor_Add(AFP_Cor,0,ret_rot,mol1,mol2,moln1,moln2,ali3);
		ret_val=ali3_kabsch(ret_rot,RMSD,mol1,mol2,moln1,moln2,ali3); //apply ALI3(weight) rotation
		if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
		if(ret_val<0.0 || RMSD<0.0)return -1; //incorrect
		//Replace
		Ali3_Cor_Add(AFP_Cor,Range,ret_rot,mol1,mol2,moln1,moln2,ali3);
		Ali3_DynaProg(ret_rot,DP_BETA,d_0,mol1,mol2,moln1,moln2,FOR_GAP,BAK_GAP,SCALE,ali3);
		score=(int)Ali3_TraceBack(ret_rot,d_0,mol1,mol2,moln1,moln2,ali3,ali1,ali2);
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
				memcpy(ret_rot,rot_bak,12*sizeof(double));
				memcpy(ali1,tali1,moln1*sizeof(int));
				memcpy(ali2,tali2,moln2*sizeof(int));
				if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
				return (bak_score);
			}
			success=0;
			count++;
		}
	}
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	return score; //[refine_better]//__081130__//
}
//3. Kill NonLinear (ret_rot won't be updated)
int CLEFAPS::Kill_NonLinear(double *ret_rot)
{
	AFP_Kill_NonLinear(AFP_Cor,AFP_Cor_temp,ret_rot);
	return 1;
}
//4. Final Stage (combine Shrink & Elong)
int CLEFAPS::Final_Stage(double *ret_rot)
{
	//--Final_Statis--//
	double RMSD;
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	Elongation(ret_rot,FIN_CUT,INT_MAX_NUM,ali1,ali2,mol1,mol2,moln1,moln2);
	Shrinking(ret_rot,FIN_CUT,ANG_CUT,ali1,ali2,mol1,mol2,moln1,moln2);
	Final_Rot_Ali(ret_rot,RMSD,mol1,mol2,moln1,moln2,ali1,ali2);
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	Elongation(ret_rot,FIN_CUT,ANG_CUT,ali1,ali2,mol1,mol2,moln1,moln2);
	return 1;
}
//5. TMscore refinement (ret_rot might be updated)
int CLEFAPS::Final_Refine(double *ret_rot)
{
	Calc_TM_Align(mol1,mol2,moln1,moln2,ali2,ali2);    //apply TMscore dynamic programming
	TM_Align_TM_Score_Simp(mol1,mol2,moln1,moln2,ali2,RMSD,LALI); //use simple version
	EqualArray(ret_rot,finmat,12);
	return 1;
}

//==================== virtual function ====================//-> part[3]: global
//0. Init functions
void CLEFAPS::Global_CLEPAPS_Init(void)
{
	int i;
	align_rec.resize(moln1*moln2);
	for(i=0;i<moln1*moln2;i++)align_rec[i]=0;
}
//1. Select best pivot
int CLEFAPS::Select_Best_Pivot(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,int TopJ)
{
	int i;
	int relnum=(int)SFP.size();
	if(relnum==0)return 0;
	int totnum;
	totnum=TopJ;
	totnum=totnum<CLEF_Maximal?totnum:CLEF_Maximal;
	totnum=totnum<relnum?totnum:relnum;
	//real_start
	for(i=0;i<totnum;i++)
	{
		head_score[i]=Select_Single(i,recur,INI_CUT,FIN_CUT,SFP,CLEP_rotmat);
		EqualArray(center_mat+12*i,CLEP_rotmat,12);
	}
	//fast_sort
	fast_sort.fast_sort_1(head_score,head_index,totnum);
	int smaller=moln1<moln2?moln1:moln2;
	CUR_MaxJ=1.0*head_score[head_index[0]]/smaller;
	//fast_check
	if(FAST_Chk==1)
	{
//		printf("%d\n",FAST_Chk);
//		printf("%lf\n",CUR_MaxJ);
		if(CUR_MaxJ<CUR_MaxJ_thres)return 0;
	}
	return totnum;
}
//2. Break Check
int CLEFAPS::CLEPAPS_Break_Check(int cur,int count,int totnum)
{
	if(count>=ZM_TopK)return 0;
	else return 1;
}
//3. Get Index
int CLEFAPS::CLEPAPS_Get_Index(int cur)
{
	CUR_INDEX=cur;
	INDEX_CUR=head_index[cur];
	return head_index[cur];
}
//4. Head Check
int CLEFAPS::CLEPAPS_Head_Check(SFP_Record &SFP)
{
	int ii=SFP.ii;
	int jj=SFP.jj;
	int winlen=SFP.winlen;

#ifdef DEBUG
	overrange_debug_test(ii+winlen-1,moln1);
	overrange_debug_test(jj+winlen-1,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

	int k;
	for(k=0;k<winlen;k++)
	{
		if(align_rec.at((ii+k)*moln2+(jj+k))==0)return 1; //un-possessed!!
	}
	return 0;
}
//5. Get initial rotmat
void CLEFAPS::CLEPAPS_Initial_Rotmat(int index,double *rotmat)
{
	EqualArray(rotmat,center_mat+12*index,12);
}
//6. Tail Check
int CLEFAPS::CLEPAPS_Tail_Check(Align_Record &cur)
{
	//init_check
	if(cur.score>CLEP_MAX.score)return 1;
	//process
	int i;
	int ii,jj;
	int count=0;
	int match=0;
	for(i=0;i<moln2;i++)
	{
		ii=cur.alignment[i];
		jj=i;
		if(ii==-1)continue;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

		count++;
		if(align_rec.at((ii)*moln2+(jj))==1)match++; //possessed!!
	}
	if(count==0)return 0;
	if(1.0*match/count>Tail_Check)return 0;
	return 1;
}
//7. Comparison Check
int CLEFAPS::CLEPAPS_Comparison_Check(Align_Record &cur,vector <Align_Record> &tot)
{
	//init_check
	if(cur.score>CLEP_MAX.score)return 1;
	//process
	int i;
	double score;
	int totnum=(int)tot.size();
	for(i=0;i<totnum;i++)
	{
		score=Alg_Comp_Simp(cur.alignment,tot[i].alignment,Comp_Check_Vari);
		if(score>Comp_Check)return 0;
	}
	return 1;
}
//8. Terminal Process
void CLEFAPS::CLEPAPS_Terminal_Process(Align_Record &cur)
{
	int i;
	int ii,jj;
	for(i=0;i<moln2;i++)
	{
		ii=cur.alignment[i];
		jj=i;
		if(ii==-1)continue;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

		align_rec.at((ii)*moln2+(jj))=1; //possess!!
	}
	//terminal_add
	cur.index.clear();
	cur.index.push_back(INDEX_CUR);
	cur.index.push_back(CUR_INDEX);
}
