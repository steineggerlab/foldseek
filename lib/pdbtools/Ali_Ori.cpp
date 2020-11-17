#include "Ali_Ori.h"


//--------------constructor--------------------//
Ali_Ori::Ali_Ori(int num)
{
	Ali_Maximal=num;
	Ali_molt=new XYZ[Ali_Maximal];
	Ali_cache=new int[Ali_Maximal];
	AFP_tmp1=new XYZ[Ali_Maximal];
	AFP_tmp2=new XYZ[Ali_Maximal];
	ALI_CACHE=0; //no cache
	//-> no mask
	mas1=0;
	mas2=0;
}
Ali_Ori::~Ali_Ori(void)
{
	delete [] Ali_molt;
	delete [] Ali_cache;
	delete [] AFP_tmp1;
	delete [] AFP_tmp2;
}


//========= Cache_Point =========//
void Ali_Ori::Ali_Cache_Point(XYZ *mol,int index,XYZ &ret_point,double *rotmat)
{
	if(ALI_CACHE==1 && Ali_cache[index]==1)ret_point=Ali_molt[index];
	else
	{
		rot_point(mol[index],ret_point,rotmat);
		Ali_molt[index]=ret_point;
		Ali_cache[index]=1;
	}
}

//------------------ ali_ori --------------------//
//-----------//
//ali3_init
//--------------------------//
void Ali_Ori::Ali_Init(int moln1,int moln2,int *ali1,int *ali2)
{
	int i;
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)ali2[i]=-1;
}

//-----------//
//generate aliI from aliJ
//--------------------------//
void Ali_Ori::Ali1_To_Ali2(int moln1,int moln2,int *ali1,int *ali2)
{
	int i;
	int ii,jj;
	for(i=0;i<moln2;i++)ali2[i]=-1;
	for(i=0;i<moln1;i++)
	{
		if(ali1[i]>=0)
		{
			ii=i;
			jj=ali1[i];

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif	

			ali2[jj]=ii;
		}
	}
}
void Ali_Ori::Ali2_To_Ali1(int moln1,int moln2,int *ali1,int *ali2)
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

//----------//
//final_rot -> mol2 is fixed!!
//-----------------------------------------------//
double Ali_Ori::Final_Rot_Ali(double *rotmat_,double &rmsd,
							  XYZ *mol1,XYZ *mol2,int moln1,int moln2,
							  int *ali1,int *ali2)
{
//input the ali_cor(One-to-One)[no-weighted]
	int i;
	int ii,jj;
	int count=0;
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
			AFP_tmp1[count]=mol1[ii];
			AFP_tmp2[count]=mol2[jj];
			count++;
		}
	}
	rmsd=kabsch(AFP_tmp2,AFP_tmp1,count,rotmat_);
	return count;
}

//-------------//
//make_center -> mol2 is fixed!! superimpose mol1 onto mol2
//--------------------------------------------------//
double Ali_Ori::make_center(int ii,int jj,int len,double *rotmat_,
							XYZ *mol1,XYZ *mol2,int moln1,int moln2)
{
	double rmsd;
#ifdef DEBUG
	overrange_debug_test(ii+len-1,moln1);
	overrange_debug_test(jj+len-1,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

	rmsd=kabsch(mol2+jj,mol1+ii,len,rotmat_);
	if(rmsd<0.0)return -1.0;
	return rmsd;
}

//----------------//
//calc_frag_dist -> mol2 is fixed!!
//--------------------------------------------------//
double Ali_Ori::calc_frag_dist(int ii,int jj,int len,double *rotmat_,
							   XYZ *mol1,XYZ *mol2,int moln1,int moln2)
{
#ifdef DEBUG
	overrange_debug_test(ii+len-1,moln1);
	overrange_debug_test(jj+len-1,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

	return calc_dist(mol1+ii,mol2+jj,len,rotmat_);
}

//----------------//
//calc_frag_dist_thres -> mol2 is fixed!!
//--------------------------------------------------//
//input [ii][jj][len][rotmat_],return [max_dist],
//and neo [ii_][jj_][len_] within [thres]
double Ali_Ori::calc_frag_dist_thres(int ii,int jj,int len,double *rotmat_,
									 XYZ *mol1,XYZ *mol2,int moln1,int moln2,
									 double thres,int &neo_ii,int &neo_jj,int &neo_len,double &score)
{
#ifdef DEBUG
	overrange_debug_test(ii+len-1,moln1);
	overrange_debug_test(jj+len-1,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

	int i;
	double dist;
	XYZ temp;
	double max;
	int valid_len;
	double dou_thres=thres*thres;

	neo_len=0;
	valid_len=0;
	score=0.0;
	max=0.0;
	for(i=0;i<len;i++)
	{
		Ali_Cache_Point(mol1,ii+i,temp,rotmat_);
		dist=mol2[jj+i].distance_square(temp);
		if(dist>max)max=dist;
		if(dist<dou_thres)valid_len++;
		else
		{
			if(valid_len>neo_len)
			{
				neo_len=valid_len;
				neo_ii=ii+i-neo_len;
				neo_jj=jj+i-neo_len;
			}
			valid_len=0;
		}
		score+=dist;
	}
	if(valid_len>neo_len)
	{
		neo_len=valid_len;
		neo_ii=ii+i-neo_len;
		neo_jj=jj+i-neo_len;
	}
	return max;
}
