#include "Ali_Ali3.h"


//--------------constructor--------------------//
Ali_Ali3::Ali_Ali3(int num)
:Ali_Ori(num)
{
	Ali3_Maximal=num;
	//init ali3's parameter
	ali3_TOT=6;         // defaul:6
	//init ali3_DP's parameter
	FOR_GAP=50;         // Ali3_DynaProg For_Gap penalty  [default:50]
	BAK_GAP=200;        // Ali3_DynaProg Bak_Gap penalty  [default:200]
	SCALE=100;          // Ali3_DynaProg Scoring scaling  [default:100]
	//create
	Ali3_Create(Ali3_Maximal);
}
Ali_Ali3::~Ali_Ali3(void)
{
	//delete
	Ali3_Delete();
}




//------------------------------------//
//-- Ali3_Universal_Correspondence --//
//----------------------------------//
//----------------//---------------------------------------------------------------------Chapter_II--PART_I
//-Ali3_Related--//
//--------------//
//Ali3:DataStructure//
//	ali3[i][*][*] the i-th index of mol2(smaller)
//	ali3[i][0][0] the total correspondece of i-th index(default:0)
//	ali3[i][0][1] the minimal correspondence of i-th index(default:0,then k)
//  ali3[i][k][*] the k-th correspondence of the i-th index
//	ali3[i][k][0] the k-th correspondence of position in mol1(bigger)
//	ali3[i][k][1] the k-th correspondence of distance in mol1(bigger)(integer)
//pre_sco:DataStructure//
//	pre_sco[i][*][*] the i-th index of mol2(smaller)
//	pre_sco[i][0][*] the null state
//	pre_sco[i][k][*] the k-th state
//	pre_sco[i][x][0] the dynamic programming score
//	pre_sco[i][x][1] the correspondence in mol1
void Ali_Ali3::Ali3_Create(int totlen)
{
	ali3=new int[totlen*(ali3_TOT+1)*2];
	pre_sco=new int[totlen*(ali3_TOT+1)*2];
	pre_temp=new int[ali3_TOT+1];
}
void Ali_Ali3::Ali3_Delete(void)
{
	delete [] ali3;
	delete [] pre_sco;
	delete [] pre_temp;
}


//-------------//
//Ali3_Init
//----------------------------------------//
void Ali_Ali3::Ali3_Init(int *ali3,int moln2)  // ali3_init
{
	int i;
	for(i=0;i<moln2;i++)
	{
		ali3[i*(ali3_TOT+1)*2+0]=0;
		ali3[i*(ali3_TOT+1)*2+1]=0;
	}
}

//-------------// -> mol2 is fixed!!!
//ALI3_Weighted//(ws_modification)
//----------------------------------------//
double Ali_Ali3::ali3_kabsch(double *rotmat_,double &rmsd,
							 XYZ *mol1, XYZ *mol2,int moln1,int moln2,
							 int *ali3)
{
//input the ALI3_Cor(One-to-Multi)[weighted]

	XYZ xc,yc;
	XYZ zero;
	double xtot,ytot;
	int i,j,k; 
	int winlen;
	int ii,jj;
	double weight;
	double weight_all;
	double result;

	/* Find center of each set of coordinates. */ 
	weight_all=0.0;
	xc=0.0;
	yc=0.0;
	for(i=0;i<moln2;i++)
	{
		winlen=ali3[i*(ali3_TOT+1)*2+0];
		int cur_index=i*(ali3_TOT+1)*2;
		for(k=1;k<=winlen;k++)
		{			
			ii=ali3[cur_index+k*2+0];
			jj=i;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif	

			//calculate weight
			weight=ali3[cur_index+k*2+1];
			weight=exp(-1.0*(weight)/(10.0*MAX_DIST*MAX_DIST));
			//[weighted_matrix]!!
			xc+=mol1[ii]*weight;
			yc+=mol2[jj]*weight;
			//add weight
			weight_all+=weight;
		}
	}
	if(IsZero(weight_all))return -1.0;
	xc/=weight_all;
	yc/=weight_all;

	//get xtot & ytot//
	/*
	* Initialise and then fill the r matrix.
	* Note that centre is subtracted at this stage.
	*/  
	zero=0.0;
	xtot=0.0;
	ytot=0.0;
	for(i=0;i<3;i++)for(j=0;j<3;j++)kabsch_r[i][j]=0.0;
	for(i=0;i<moln2;i++)
	{
		winlen=ali3[i*(ali3_TOT+1)*2+0];
		int cur_index=i*(ali3_TOT+1)*2;
		for(k=1;k<=winlen;k++)
		{			
			ii=ali3[cur_index+k*2+0];
			jj=i;

			//calculate weight
			weight=ali3[cur_index+k*2+1];
			weight=exp(-1.0*(weight)/(10.0*MAX_DIST*MAX_DIST));
			//[weighted_matrix]!!
			xtot+=((mol1[ii]-xc).distance_square(zero))*weight;
			ytot+=((mol2[jj]-yc).distance_square(zero))*weight;
			//[weighted_matrix]!!
			kabsch_r[0][0]+=(mol2[jj].X-yc.X)*(mol1[ii].X-xc.X)*weight;
			kabsch_r[1][0]+=(mol2[jj].X-yc.X)*(mol1[ii].Y-xc.Y)*weight;
			kabsch_r[2][0]+=(mol2[jj].X-yc.X)*(mol1[ii].Z-xc.Z)*weight;
			kabsch_r[0][1]+=(mol2[jj].Y-yc.Y)*(mol1[ii].X-xc.X)*weight;
			kabsch_r[1][1]+=(mol2[jj].Y-yc.Y)*(mol1[ii].Y-xc.Y)*weight;
			kabsch_r[2][1]+=(mol2[jj].Y-yc.Y)*(mol1[ii].Z-xc.Z)*weight;
			kabsch_r[0][2]+=(mol2[jj].Z-yc.Z)*(mol1[ii].X-xc.X)*weight;
			kabsch_r[1][2]+=(mol2[jj].Z-yc.Z)*(mol1[ii].Y-xc.Y)*weight;
			kabsch_r[2][2]+=(mol2[jj].Z-yc.Z)*(mol1[ii].Z-xc.Z)*weight;
		}
	}
	zero=xc;
	xc=yc;
	yc=zero;

	//input R,XC,YC to kabsch_basic to get the rot_mat(d)
	result=kabsch_base(kabsch_r,xc,yc,rotmat_);
	rmsd=(xtot+ytot-2*result)/weight_all;
	if(IsZero(rmsd))rmsd=0.0;
	return weight_all;
}

//-------------//
//Ali3_Add -> mol2 is fixed!!!
//----------------------------------------//
void Ali_Ali3::Ali3_Add(int ii,int jj,int len,double *rotmat_,
						XYZ *mol1,XYZ *mol2,int moln1,int moln2,
						int *ali3,int *zali1,int *zali2)
{
#ifdef DEBUG
	overrange_debug_test(ii+len-1,moln1);
	overrange_debug_test(jj+len-1,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

	int k,l;
	int curlen;
	int rellen;
	int ws_score;
	int wsmax,maxnum;
	int ws_index;
	int index;
	int score;
	XYZ temp;

	for(k=0;k<len;k++)
	{
		//head_check
		if(zali1==NULL || zali2==NULL)goto add;		
		if(zali1[ii+k]!=-1||zali2[jj+k]!=-1)continue;

add:
		int cur_index=(jj+k)*(ali3_TOT+1)*2;
		//real_add
		rellen=ali3[cur_index+0];
		for(l=1;l<=rellen;l++)
		{
			if(ali3[cur_index+l*2+0]==(ii+k))goto ws_next;
		}
		Ali_Cache_Point(mol1,ii+k,temp,rotmat_);
		ws_score=(int)(10.0*(mol2[jj+k].distance_square(temp)));
		if(rellen<ali3_TOT)
		{
			ali3[cur_index+0]++;
			ws_index=ali3[cur_index+0];
			ali3[cur_index+ws_index*2+0]=ii+k;
			ali3[cur_index+ws_index*2+1]=ws_score;
			//record_maximal//
			index=ali3[cur_index+1];
			if(index==0)ali3[cur_index+1]=1;
			else
			{
				if(ws_score>ali3[cur_index+index*2+1])ali3[cur_index+1]=ws_index;
			}
		}
		else 
		{
			ws_index=ali3[cur_index+1];
			if(ws_score<ali3[cur_index+ws_index*2+1])
			{
				ali3[cur_index+ws_index*2+0]=ii+k;
				ali3[cur_index+ws_index*2+1]=ws_score;
				//get_maximal
				wsmax=ali3[cur_index+1*2+1];
				maxnum=1;
				curlen=ali3[cur_index+0];
				for(l=2;l<=curlen;l++)
				{
					score=ali3[cur_index+l*2+1];
					if(score>wsmax)
					{
						wsmax=score;
						maxnum=l;
					}
				}
				ali3[cur_index+1]=maxnum;
			}
		}
ws_next:
		;
	}
}

//-------------//
//Ali3_Cor_Add -> mol2 is fixed!!
//----------------------------------------//
void Ali_Ali3::Ali3_Cor_Add(int *AFP_Cor,int Range,double *rotmat_,
							XYZ *mol1,XYZ *mol2,int moln1,int moln2,
							int *ali3,int *zali1,int *zali2)
{
	int i,j;
	int totnum;
	int ii,jj,len;
	int kk,ll;
	int fk,bk;

	//init
	Ali3_Init(ali3,moln2);

	//start
	totnum=AFP_Cor[0];
	for(i=1;i<=totnum;i++)
	{
		if(AFP_Cor[i*4+0]<0)continue;
		ii=AFP_Cor[i*4+1];
		jj=AFP_Cor[i*4+2];
		len=AFP_Cor[i*4+3];

#ifdef DEBUG
	overrange_debug_test(ii+len-1,moln1);
	overrange_debug_test(jj+len-1,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

		for(j=0;j<2*Range+1;j++)
		{
			//[0]center
			if(j==Range)
			{
				for(bk=0;bk<Range;bk++)if(ii-bk<=0||jj-bk<=0)break;
				for(fk=0;fk<Range;fk++)if(ii+len+fk>=moln1||jj+len+fk>=moln2)break;
				Ali3_Add(ii-bk,jj-bk,len+bk+fk,rotmat_,mol1,mol2,moln1,moln2,ali3,zali1,zali2);
				continue;
			}
			//[1]left
			kk=ii-Range+j;
			ll=jj;
			if(kk<0)goto next;
			if(kk+len>moln1)goto next;
			Ali3_Add(kk,ll,len,rotmat_,mol1,mol2,moln1,moln2,ali3,zali1,zali2);
next:			
			//[right]
			kk=ii;
			ll=jj-Range+j;
			if(ll<0)continue;
			if(ll+len>moln2)continue;
			Ali3_Add(kk,ll,len,rotmat_,mol1,mol2,moln1,moln2,ali3,zali1,zali2);
		}
	}
}

//-------------//
//Ali3_Cor_Add -> mol2 is fixed!!
//----------------------------------------//
void Ali_Ali3::Ali3_To_Ali(double *rotmat_,int method,
						   XYZ *mol1,XYZ *mol2,int moln1,int moln2,
						   int *ali3,int *ali1,int *ali2)
{
	int i,l;
	int index;
	int ii,jj;
	int rellen;
	int wsmax,maxnum;
	int ori_score,score;
	XYZ temp;

	//ali_init
	Ali_Init(moln1,moln2,ali1,ali2);
	//ali3->ali
	for(i=0;i<moln2;i++)
	{
		rellen=ali3[i*(ali3_TOT+1)*2+0];
		if(rellen>0)
		{
			int cur_index=i*(ali3_TOT+1)*2;

			//part[1]:pick minimal
			wsmax=ali3[cur_index+1*2+1];
			maxnum=1;
			for(l=2;l<=rellen;l++)
			{
				if(ali3[cur_index+l*2+1]<wsmax)
				{
					wsmax=ali3[cur_index+l*2+1];
					maxnum=l;
				}
			}
			//part[2]:evaluate
			ori_score=ali3[cur_index+maxnum*2+1];
			ii=ali3[cur_index+maxnum*2+0];
			jj=i;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif
			
			//part[3]:check conflict
			if(ali1[ii]==-1)
			{
				ali1[ii]=jj;
				ali2[jj]=ii;
			}
			else
			{
				//method==1 (First_Come_Never_Gone)
				//method!=1 (check conflict)
				if(method==1)continue;
				//real_add
				index=ali1[ii];
				Ali_Cache_Point(mol1,ii,temp,rotmat_);
				score=(int)(10.0*(mol2[index].distance_square(temp)));
				if(ori_score<score)
				{
					ali2[index]=-1;
					ali1[ii]=jj;
					ali2[jj]=ii;
				}
			}
		}//end of IF(ali3[i][0][0]>0)
		else ali2[i]=-1;
	}//end of FOR(i)
}

//-------------//
//Ali3_Out
//----------------------------------------//
void Ali_Ali3::Ali3_Out(FILE *fp,int moln2,int *ali3)
{
	int i,k;
	int totnum;
	int minimal;
	int pos;
	int weight;
	for(i=0;i<moln2;i++)
	{
		totnum=ali3[i*(ali3_TOT+1)*2+0];
		minimal=ali3[i*(ali3_TOT+1)*2+1];
		fprintf(fp,"%4d (%d) [%d] -> ",i+1,totnum,minimal);
		for(k=1;k<=totnum;k++)
		{
			pos=ali3[i*(ali3_TOT+1)*2+k*2+0];
			weight=ali3[i*(ali3_TOT+1)*2+k*2+1];
			fprintf(fp,"[%4d|%4d]",pos+1,weight);
		}
		fprintf(fp,"\n");		
	}
}


//--------------------------------//
//-- Ali3_Dynamic_Programming  --//
//------------------------------//
//----------------//---------------------------------------------------------------------Chapter_II--PART_I
//--- pre_sco ---//
//--------------//
//-------------//
//Ali3_Vector_Score_Forward -> mol2 is fixed!!
//----------------------------------------//
int Ali_Ali3::Ali3_Vector_Score_Forward(int ii,int jj,double *rotmat_,
										XYZ *mol1,XYZ *mol2,int moln1,int moln2,
										int SCALE)
{
#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

	XYZ temp1,temp2;
	double v1[3],v2[3];
	double d1,d2;
	double score;
	if(ii==0||jj==0)return 0;
	else
	{
		Ali_Cache_Point(mol1,ii,temp1,rotmat_);
		Ali_Cache_Point(mol1,ii-1,temp2,rotmat_);
		(temp1-temp2).xyz2double(v1);
		(mol2[jj]-mol2[jj-1]).xyz2double(v2);
		d1=sqrt(dot(v1,v1,3));
        d2=sqrt(dot(v2,v2,3));
		//out the range of C_ALpha's distance
		if(d1>4.5||d2>4.5||d1<2.5||d2<2.5)return -1*SCALE;
		score=dot(v1,v2,3)/d1/d2;
		score=limit(score,-1.0,1.0);
		return (int)(SCALE*score);
	}
}

//-------------//
//Ali3_Vector_Score_Backward -> mol2 is fixed!!
//----------------------------------------//
int Ali_Ali3::Ali3_Vector_Score_Backward(int ii,int jj,double *rotmat_,
										 XYZ *mol1,XYZ *mol2,int moln1,int moln2,
										 int SCALE)
{
#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

	XYZ temp1,temp2;
	double v1[3],v2[3];
	double d1,d2;
	double score;
	if(ii==moln1-1||jj==moln2-1)return 0;
	else
	{
		Ali_Cache_Point(mol1,ii,temp1,rotmat_);
		Ali_Cache_Point(mol1,ii+1,temp2,rotmat_);
		(temp1-temp2).xyz2double(v1);
		(mol2[jj]-mol2[jj+1]).xyz2double(v2);
		d1=sqrt(dot(v1,v1,3));
        d2=sqrt(dot(v2,v2,3));
		//out the range of C_ALpha's distance
		if(d1>4.5||d2>4.5||d1<2.5||d2<2.5)return -1*SCALE;
		score=dot(v1,v2,3)/d1/d2;
		score=limit(score,-1.0,1.0);
		return (int)(SCALE*score);
	}
}

//note:ali3_num must less than 10
//     it's easy to understand,because ali3_num 
//     is the number that mol2[i] correspond to 
//     mol1[j],the near neibor which less than 
//     a certain threshold cannot exceed 10
//-------------//
//Ali3_DynaProg -> mol2 is fixed!!
//----------------------------------------//
void Ali_Ali3::Ali3_DynaProg(double *rotmat_,double beta,double d_0,
							 XYZ *mol1,XYZ *mol2,int moln1,int moln2,
							 int FOR_GAP,int BAK_GAP,int SCALE,int *ali3)
{
	int i,j,k;
	int num,num_bak;
	int pos,pos_bak;
	int score,dlscore,vtscore;
	int index;
	int ii,jj;
	int munus;
	int count;
	int cur;
	int gap_penalty;
	int gap_extend=FOR_GAP/10;
	int extend_cut=(SCALE-FOR_GAP)/gap_extend;
	double dist;
	double ori_d=d_0*d_0;
	XYZ temp;

	//-------Dynamic_Programming--------//
	//init//
	temp=0.0;
	num=ali3[0];
	pre_sco[0]=0;
	pre_sco[1]=0;
	for(k=1;k<=num;k++)
	{
		pre_sco[0*(ali3_TOT+1)*2+k*2+0]=0;
		pre_sco[0*(ali3_TOT+1)*2+k*2+1]=INT_MAX_NUM*ali3[0*(ali3_TOT+1)*2+k*2+0];
	}

	//fill_in//
	for(i=1;i<moln2;i++)
	{
		//state==0 (null state)
		num_bak=ali3[(i-1)*(ali3_TOT+1)*2+0];
		for(j=0;j<=num_bak;j++)pre_temp[j]=pre_sco[(i-1)*(ali3_TOT+1)*2+j*2+0];
		index=Get_Max(pre_temp,num_bak+1);
		pos_bak=pre_sco[(i-1)*(ali3_TOT+1)*2+index*2+1]/INT_MAX_NUM;
		pre_sco[i*(ali3_TOT+1)*2+0]=pre_temp[index]>0?pre_temp[index]:0;
		pre_sco[i*(ali3_TOT+1)*2+1]=index+INT_MAX_NUM*pos_bak;
		//state!=0 (real state)
		num=ali3[i*(ali3_TOT+1)*2+0];
		for(k=1;k<=num;k++)
		{
			int cur_index=(i-1)*(ali3_TOT+1)*2;
			int nxt_index=i*(ali3_TOT+1)*2;

			//get_cur_pos
			jj=i;
			ii=ali3[nxt_index+k*2+0];

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

			//calculate_score
			//[1]dist_score
			Ali_Cache_Point(mol1,ii,temp,rotmat_);
			dist=mol2[jj].distance_square(temp);
			dlscore=(int)(1.0*SCALE/(1.0+dist/ori_d));
			//[2](cl+vt)_score
			if(beta>0.0)vtscore=Ali3_Vector_Score_Forward(ii,jj,rotmat_,mol1,mol2,moln1,moln2,SCALE);
			else vtscore=0;
			score=(int)((1.0-beta)*dlscore+beta*vtscore);  //ws_ali3_neo//__081020__//
			//calculate_trans
			pos=ii;
			for(j=0;j<=num_bak;j++)
			{
				pos_bak=pre_sco[cur_index+j*2+1]/INT_MAX_NUM;
				if(pos==pos_bak+1)
				{
					if(j==0)
					{
						//trace_back
						count=0;
						cur=i-1;
						for(;;)
						{
							if(cur<0)break;
							if(count>extend_cut)break;
							index=pre_sco[cur*(ali3_TOT+1)*2+1]%INT_MAX_NUM;
							if(index!=0)break;
							count++;
							cur--;
						}
						//calc_score
						gap_penalty=FOR_GAP+count*gap_extend;
						if(gap_penalty>SCALE)gap_penalty=SCALE;
						pre_temp[j]=pre_sco[cur_index+j*2+0]+score-gap_penalty;            //default:100
					}
					else pre_temp[j]=pre_sco[cur_index+j*2+0]+score;
				}
				else
				{
					if(pos<pos_bak)pre_temp[j]=pre_sco[cur_index+j*2+0]+score-BAK_GAP;    //default:99999
					else
					{
						munus=abs(pos-pos_bak);
						gap_penalty=FOR_GAP+munus*gap_extend;
						if(gap_penalty>SCALE)gap_penalty=SCALE;
						pre_temp[j]=pre_sco[cur_index+j*2+0]+score-gap_penalty;            //default:100
					}
				}
			}
			//get_max
			index=Get_Max(pre_temp,num_bak+1);
			pre_sco[nxt_index+k*2+0]=pre_temp[index]>0?pre_temp[index]:0;
			pre_sco[nxt_index+k*2+1]=index+INT_MAX_NUM*pos;
		}
	}
}

//-------------//
//Ali3_TraceBack -> mol2 is fixed!!
//----------------------------------------//
double Ali_Ali3::Ali3_TraceBack(double *rotmat_,double d_0,
								XYZ *mol1,XYZ *mol2,int moln1,int moln2,
								int *ali3,int *ali1,int *ali2)
{
	int i,k;
	int num;
	int index;
	int ii,jj,winlen;
	int kk,ll;
	int ww,ss;
	int isFirst;
	int isLast;
	int count;
	double ret_val;
	double dist,dist1,dist2;
	double ori_d=d_0*d_0;
	XYZ temp;
	XYZ temp1,temp2;

	//---Trace_Back---//
	//init//
	ii=-1;
	jj=-1;
	kk=-1;
	ll=-1;
	temp=0.0;
	num=ali3[(moln2-1)*(ali3_TOT+1)*2+0];
	for(k=0;k<=num;k++)pre_temp[k]=pre_sco[(moln2-1)*(ali3_TOT+1)*2+k*2+0];
	index=Get_Max(pre_temp,num+1);
	isFirst=1;
	isLast=0;
	//ali_init
	Ali_Init(moln1,moln2,ali1,ali2);
	//trace_back//
	winlen=0;
	ret_val=0.0;
	for(i=moln2-1;i>=0;i--)
	{
		if(index==0) //null_state
		{
			if(isFirst==1) //no_record,continue
			{
			}
			else  //has_record,break & continue
			{
				//init
				isFirst=1;
				//check & add
				count=0;
				for(k=0;k<winlen;k++)
				{
					//record1&2
					if(ali1[ii+k]!=-1&&ali2[jj+k]!=-1)continue;
					else if(ali1[ii+k]==-1&&ali2[jj+k]==-1)
					{
						ali1[ii+k]=jj+k;
						ali2[jj+k]=ii+k;
						count++;
						//add
						Ali_Cache_Point(mol1,ii+k,temp,rotmat_);
						dist=mol2[jj+k].distance_square(temp);
						ret_val+=1.0/(1.0+dist/ori_d);
					}
					else
					{
						if(ali1[ii+k]!=-1&&ali2[jj+k]==-1)
						{
							ww=ali1[ii+k];
							Ali_Cache_Point(mol1,ii+k,temp1,rotmat_);
							dist1=mol2[ww].distance_square(temp1);
							dist2=mol2[jj+k].distance_square(temp1);
							if(dist1>dist2)
							{
								ali2[ww]=-1;
								ali1[ii+k]=jj+k;
								ali2[jj+k]=ii+k;
								ret_val-=1.0/(1.0+dist1/ori_d);  //minus
								ret_val+=1.0/(1.0+dist2/ori_d);  //add
							}
						}
						if(ali1[ii+k]==-1&&ali2[jj+k]!=-1)
						{
							ss=ali2[jj+k];
							Ali_Cache_Point(mol1,ss,temp1,rotmat_);
							Ali_Cache_Point(mol1,ii+k,temp2,rotmat_);
							dist1=mol2[jj+k].distance_square(temp1);
							dist2=mol2[jj+k].distance_square(temp2);
							if(dist1>dist2)
							{
								ali1[ss]=-1;
								ali2[jj+k]=ii+k;
								ali1[ii+k]=jj+k;
								ret_val-=1.0/(1.0+dist1/ori_d);  //minus
								ret_val+=1.0/(1.0+dist2/ori_d);  //add
							}
						}
					}
				}//end of FOR(k)
			}//end of has_record
		}
		else  //real_state
		{
			if(isFirst==1) //no_record,start
			{
				//init_rec
				isFirst=0;
				ii=ali3[i*(ali3_TOT+1)*2+index*2+0];
				jj=i;
				winlen=1;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif				
			}
			else  //has_record,check
			{
				//check contiguous
				kk=ali3[i*(ali3_TOT+1)*2+index*2+0];
				ll=i;
				if(kk==ii-1 && ll==jj-1)
				{
					ii--;
					jj--;
					winlen++;
				}
				else
				{
					//check & add
ws_last:
					count=0;
					for(k=0;k<winlen;k++)
					{
						//record1&2
						if(ali1[ii+k]!=-1&&ali2[jj+k]!=-1)continue;
						else if(ali1[ii+k]==-1&&ali2[jj+k]==-1)
						{
							ali1[ii+k]=jj+k;
							ali2[jj+k]=ii+k;
							count++;
							//add
							Ali_Cache_Point(mol1,ii+k,temp,rotmat_);
							dist=mol2[jj+k].distance_square(temp);
							ret_val+=1.0/(1.0+dist/ori_d);
						}
						else
						{
							if(ali1[ii+k]!=-1&&ali2[jj+k]==-1)
							{
								ww=ali1[ii+k];
								Ali_Cache_Point(mol1,ii+k,temp1,rotmat_);
								dist1=mol2[ww].distance_square(temp1);
								dist2=mol2[jj+k].distance_square(temp1);
								if(dist1>dist2)
								{
									ali2[ww]=-1;
									ali1[ii+k]=jj+k;
									ali2[jj+k]=ii+k;
									ret_val-=1.0/(1.0+dist1/ori_d);  //minus
									ret_val+=1.0/(1.0+dist2/ori_d);  //add
								}
							}
							if(ali1[ii+k]==-1&&ali2[jj+k]!=-1)
							{
								ss=ali2[jj+k];
								Ali_Cache_Point(mol1,ss,temp1,rotmat_);
								Ali_Cache_Point(mol1,ii+k,temp2,rotmat_);
								dist1=mol2[jj+k].distance_square(temp1);
								dist2=mol2[jj+k].distance_square(temp2);
								if(dist1>dist2)
								{
									ali1[ss]=-1;
									ali2[jj+k]=ii+k;
									ali1[ii+k]=jj+k;
									ret_val-=1.0/(1.0+dist1/ori_d);  //minus
									ret_val+=1.0/(1.0+dist2/ori_d);  //add
								}
							}
						}
					}//end of FOR(k)
					if(isLast==1)goto end;
					//init_rec
					isFirst=0;
					ii=kk;
					jj=ll;
					winlen=1;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif	
				}
			}
		}
		//next
		index=pre_sco[i*(ali3_TOT+1)*2+index*2+1]%INT_MAX_NUM;
	}

	//check & add (last)
	isLast=1;
	goto ws_last;
end:
	return ret_val;
}

//-------------//
//Ali3_DynaProg_Out
//----------------------------------------//
void Ali_Ali3::Ali3_DynaProg_Out(FILE *fp,int moln2,int *ali3)
{
	int i,k;
	int totnum;
	int score;
	int pos,index;

	for(i=0;i<moln2;i++)
	{
		totnum=ali3[i*(ali3_TOT+1)*2+0];
		fprintf(fp,"%4d (%d) -> ",i+1,totnum);
		for(k=0;k<=totnum;k++)
		{
			score=pre_sco[i*(ali3_TOT+1)*2+k*2+0];
			pos=pre_sco[i*(ali3_TOT+1)*2+k*2+1]/INT_MAX_NUM;
			index=pre_sco[i*(ali3_TOT+1)*2+k*2+1]%INT_MAX_NUM;
			fprintf(fp,"[%6d|%4d|%4d]",score,pos+1,index);
		}
		fprintf(fp,"\n");		
	}
}

