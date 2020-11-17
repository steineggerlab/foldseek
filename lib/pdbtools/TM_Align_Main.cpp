#include "TM_Align_Main.h"


//------ constructor -------//
TM_Align_Main::TM_Align_Main(int num)
:TM_align(num)
{
	TM_ADDITION=1;    //using additional data
	TMM_maximal=num;
	TM_Align_Main_Init(TMM_maximal);
}
TM_Align_Main::~TM_Align_Main(void)
{
	TM_Align_Main_Dele();
}

//--------------- init & dele -----------------//
void TM_Align_Main::TM_Align_Main_Init(int maxlen)
{
	TM_dist=new double[maxlen];
	TM_sse1=new int[maxlen];
	TM_sse2=new int[maxlen];
	TM_frag1=new int[2*maxlen];
	TM_frag2=new int[2*maxlen];
	TM_Main_Ali2=new int[maxlen];
	TM_Main_AliT=new int[maxlen];
	TM_FINMAT=new double[12];
	TM_ttt1=new XYZ[maxlen];
	TM_ttt2=new XYZ[maxlen];
}
void TM_Align_Main::TM_Align_Main_Dele(void)
{
	delete [] TM_dist;
	delete [] TM_sse1;
	delete [] TM_sse2;
	delete [] TM_frag1;
	delete [] TM_frag2;
	delete [] TM_Main_Ali2;
	delete [] TM_Main_AliT;
	delete [] TM_FINMAT;
	delete [] TM_ttt1;
	delete [] TM_ttt2;
}

/*****************************************************************
*     quickly calculate TM-score with given invmap(i) in 3 iterations
*****************************************************************/
double TM_Align_Main::TM_Align_Get_GL(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	int i;
	int pos;
	int count;
	int lali;
	double d002;
	double d002t;
	double ori_d;
	double dist2;
	double score;
	double d1,d2,d3;

	//[1]normal TM_score
	//get correspondece
	count=0;
	for(i=0;i<moln2;i++)
	{
		pos=ali2[i];
		if(pos==-1)continue;
		TM_tmp1[count]=mol1[pos];
		TM_tmp2[count]=mol2[i];
		count++;
	}
	lali=count;
	int smaller=moln1<moln2?moln1:moln2;
	if(TM_CALC==0)Calc_TM_d0(smaller);
	ori_d=d0*d0;
	d002=d00*d00;
	//rotmat
	kabsch(TM_tmp2,TM_tmp1,lali,finmat);
	rot_mol(TM_tmp1,TM_tmp1,lali,finmat);
	//calc score
	score=0.0;
	for(i=0;i<lali;i++)
	{
		dist2=TM_tmp2[i].distance_square(TM_tmp1[i]);
		TM_dist[i]=dist2;
		score+=1.0/(1.0+dist2/ori_d);
	}
	d1=score;
	d002t=d002;

	//[2]second TM_score
	//get correspondece
ws_recur1:
	count=0;
	for(i=0;i<lali;i++)
	{
		dist2=TM_dist[i];
		if(dist2<=d002t)
		{
			tmp1[count]=TM_tmp1[i];
			tmp2[count]=TM_tmp2[i];
			count++;
		}
	}
	if(count<3&&lali>3)
	{
		d002t+=0.5;
		goto ws_recur1;
	}
	//rotmat
	kabsch(tmp2,tmp1,count,finmat);
	rot_mol(TM_tmp1,TM_tmp1,lali,finmat);
	//calc score
	score=0.0;
	for(i=0;i<lali;i++)
	{
		dist2=TM_tmp2[i].distance_square(TM_tmp1[i]);
		TM_dist[i]=dist2;
		score+=1.0/(1.0+dist2/ori_d);
	}
	d2=score;
	d002t=d002+1.0;

	//[3]third TM_score
	//get correspondece
ws_recur2:
	count=0;
	for(i=0;i<lali;i++)
	{
		dist2=TM_dist[i];
		if(dist2<=d002t)
		{
			tmp1[count]=TM_tmp1[i];
			tmp2[count]=TM_tmp2[i];
			count++;
		}
	}
	if(count<3&&lali>3)
	{
		d002t+=0.5;
		goto ws_recur2;
	}
	//rotmat
	kabsch(tmp2,tmp1,count,finmat);
	rot_mol(TM_tmp1,TM_tmp1,lali,finmat);
	//calc score
	score=0.0;
	for(i=0;i<lali;i++)
	{
		dist2=TM_tmp2[i].distance_square(TM_tmp1[i]);
		TM_dist[i]=dist2;
		score+=1.0/(1.0+dist2/ori_d);
	}
	d3=score;

	//[4]final judge
	score=d1;
	if(d2>score)score=d2;
	if(d3>score)score=d3;
	return score;
}

//--------- TM_SSE related ---------//
int TM_Align_Main::TM_Assign_SSE_Single(double d13,double d14,double d15,double d24,double d25,double d35)
{
	//[0]coil
	int make_sec=1;
	double delta;
	//[1]helix
	delta=2.1;
	if(fabs(d15-6.37)<delta)
	{
		if(fabs(d14-5.18)<delta)
		{
			if(fabs(d25-5.18)<delta)
			{
				if(fabs(d13-5.45)<delta)
				{
					if(fabs(d24-5.45)<delta)
					{
						if(fabs(d35-5.45)<delta)
						{
							make_sec=2; //helix
							return make_sec;
						}
					}
				}
			}
		}
	}
	//[2]sheet
	delta=1.42;
	if(fabs(d15-13.0)<delta)
	{
		if(fabs(d14-10.4)<delta)
		{
			if(fabs(d25-10.4)<delta)
			{
				if(fabs(d13-6.1)<delta)
				{
					if(fabs(d24-6.1)<delta)
					{
						if(fabs(d35-6.1)<delta)
						{
							make_sec=4; //sheet
							return make_sec;
						}
					}
				}
			}
		}
	}
	//[3]turn
	if(d15<8.0)
	{
		make_sec=3; //turn
		return make_sec;
	}
	//return
	return make_sec;
}
//1->coil, 2->helix, 3->turn, 4->strand
void TM_Align_Main::TM_Assign_SSE(int *sse,XYZ *mol,int moln)
{
	int i;
	int j1,j2,j3,j4,j5;
	double d13,d14,d15,d24,d25,d35;
	for(i=0;i<moln;i++)
	{
		sse[i]=1;
		j1=i-2;
		j2=i-1;
		j3=i;
		j4=i+1;
		j5=i+2;
		if(j1>=0&&j5<moln)
		{
			d13=mol[j1].distance(mol[j3]);
			d14=mol[j1].distance(mol[j4]);
			d15=mol[j1].distance(mol[j5]);
			d24=mol[j2].distance(mol[j4]);
			d25=mol[j2].distance(mol[j5]);
			d35=mol[j3].distance(mol[j5]);
			sse[i]=TM_Assign_SSE_Single(d13,d14,d15,d24,d25,d35);
		}
	}
}
void TM_Align_Main::TM_Smooth_SSE(int *sse,int moln)
{
	int i;
	int j;
	//--- smooth single ----//
	for(i=2;i<moln-2;i++)  //is a bug??//__Added at//__101205__//
	{
		if(sse[i]==2||sse[i]==4)
		{
			j=sse[i];
			if(sse[i-2]!=j)
			{
				if(sse[i-1]!=j)
				{
					if(sse[i+1]!=j)
					{
						if(sse[i+2]!=j)  //is a bug??//__Added at//__101205__//
						{
							sse[i]=1;
						}
					}
				}
			}
		}
	}
	//--- smooth double ----//
	for(i=0;i<moln-5;i++)  //is a bug??//__Added at//__101205__//
	{
		//helix
		if(sse[i]!=2)
		{
			if(sse[i+1]!=2)
			{
				if(sse[i+2]==2)
				{
					if(sse[i+3]==2)
					{
						if(sse[i+4]!=2)
						{
							if(sse[i+5]!=2)
							{
								sse[i+2]=1;
								sse[i+3]=1;
							}
						}
					}
				}
			}
		}
		//sheet
		if(sse[i]!=4)
		{
			if(sse[i+1]!=4)
			{
				if(sse[i+2]==4)
				{
					if(sse[i+3]==4)
					{
						if(sse[i+4]!=4)
						{
							if(sse[i+5]!=4)
							{
								sse[i+2]=1;
								sse[i+3]=1;
							}
						}
					}
				}
			}
		}
	}
	//--- smooth connect ----//
	for(i=0;i<moln-2;i++)
	{
		//helix
		if(sse[i]==2)
		{
			if(sse[i+1]!=2)
			{
				if(sse[i+2]==2)
				{
					sse[i+1]=2;
				}
			}
		}
		//sheet
		if(sse[i]==4)
		{
			if(sse[i+1]!=4)
			{
				if(sse[i+2]==4)
				{
					sse[i+1]=4;
				}
			}
		}
	}
}

//========================= major function ======================//
/**************************************************************
*     get initial alignment invmap0(i) from gapless threading
**************************************************************/
void TM_Align_Main::TM_Get_Initial1(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	//init
	int i,j;
	int smaller=moln1<moln2?moln1:moln2;
	if(TM_CALC==0)Calc_TM_d0(smaller);
	double d01=d0+1.5;
	if(d01<0.5)d01=0.5;
	double d02=d01*d01;
	int aL=smaller;
	int idel;
	if(TM_ADDITION==1)idel=(int)(aL/2.3); //-> additional version
	else idel=(int)(aL/2.0);
	if(idel<=5)idel=5;
	int n1=-moln2+idel;
	int n2=moln1-idel;
	double GL_max_cut=0.95;
	double GL_maxA=0.0;
	double GL_max=0.0;
	double GL_score;
	//process
	int count;
	int ishift;
	int n_jump;
	if(TM_ADDITION==1)n_jump=5; //-> additional version
	else n_jump=1;
	double n_cuta=0.75;
	for(ishift=n1;ishift<=n2;ishift+=n_jump)  // for (ishift = n1; n_jump < 0 ? ishift >= n2 : ishift <= n2; ishift +=n_jump)
	{
		count=0;
		for(j=0;j<moln2;j++)
		{
			i=j+ishift;
			if(i>=0&&i<moln1)
			{
				TM_DP_ali2[j]=i;

				//additional
				if(TM_ADDITION==1)
				{
					TM_ttt1[count].X=mol1[i].X;
					TM_ttt1[count].Y=mol1[i].Y;
					TM_ttt1[count].Z=mol1[i].Z;
					TM_ttt2[count].X=mol2[j].X;
					TM_ttt2[count].Y=mol2[j].Y;
					TM_ttt2[count].Z=mol2[j].Z;
				}
				//additional//__over__//

				count++;
			}
			else TM_DP_ali2[j]=-1;
		}
		if(count>=idel)
		{
			GL_score=TM_Align_Get_GL(mol1,mol2,moln1,moln2,TM_DP_ali2);
			if(GL_score>GL_max)
			{
				GL_max=GL_score;
				EqualArray(ali2,TM_DP_ali2,moln2);
			}
			//additional
			if(TM_ADDITION==1)
			{
				if(GL_score>GL_max*GL_max_cut)
				{
					double dist2;
					kabsch(TM_ttt2,TM_ttt1,count,finmat);       //the kabsch has ERROR!!!!//__101206__//
					rot_mol(mol1,TM_tmp1,moln1,finmat);
					for(i=0;i<moln1;i++)
					{
						for(j=0;j<moln2;j++)
						{
							dist2=TM_tmp1[i].distance_square(mol2[j]);
							TM_DP_sco[i*moln2+j]=1.0/(1.0+dist2/d02);  // this is NOT d0 !!!
						}
					}
					double gap_open=0.0;            // this is NOT defined in TM_align.f !!! //__Added at//__101205__//
					int align_len1=TM_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,TM_DP_best,gap_open);
					if(align_len1>count*n_cuta && align_len1>idel)
					{
						GL_score=TM_Align_Get_GL(mol1,mol2,moln1,moln2,TM_DP_best);
						if(GL_score>GL_max && GL_score>GL_maxA)
						{
							GL_maxA=GL_score;
							EqualArray(TM_Main_AliT,TM_DP_best,moln2);
							AliT_Rec=1;
						}
					}
				}//end of IF(GL_score>GL_max*GL_max_cut)
			}
			//additional//__over__//
		}//end of IF(count>=idel)
	}
}

/**************************************************************
*     get initial alignment invmap0(i) from secondary structure
**************************************************************/
void TM_Align_Main::TM_Get_Initial2(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	//assign SSE
	TM_Assign_SSE(TM_sse1,mol1,moln1);
	TM_Smooth_SSE(TM_sse1,moln1);
	TM_Assign_SSE(TM_sse2,mol2,moln2);
	TM_Smooth_SSE(TM_sse2,moln2);
	//get alignment path
	int i,j;
	for(i=0;i<moln1;i++)
	{
		int cur_index=i*moln2;
		for(j=0;j<moln2;j++)
		{
			if(TM_sse1[i]==TM_sse2[j])TM_DP_sco[cur_index+j]=1.0;
			else TM_DP_sco[cur_index+j]=0.0;
		}
	}
	//apply dynamic programming
	double gap_open=-1.0;
	TM_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,ali2,gap_open);
}

/**************************************************************
*     get initial alignment invmap0(i) from secondary structure 
*     and previous alignments
**************************************************************/
void TM_Align_Main::TM_Get_Initial3(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	int i,j;
	int pos;
	int count;
	//get correspondece
	count=0;
	for(i=0;i<moln2;i++)
	{
		pos=ali2[i];
		if(pos==-1)continue;
		TM_tmp1[count]=mol1[pos];
		TM_tmp2[count]=mol2[i];
		count++;
	}
	kabsch(TM_tmp2,TM_tmp1,count,finmat);
	rot_mol(mol1,TM_tmp1,moln1,finmat);
	//calculate score matrix
	int smaller=moln1<moln2?moln1:moln2;
	if(TM_CALC==0)Calc_TM_d0(smaller);
	double d01=d0+1.5;
	if(d01<0.5)d01=0.5;
	double d02=d01*d01;
	double dist2;
	for(i=0;i<moln1;i++)
	{
		int cur_index=i*moln2;
		for(j=0;j<moln2;j++)
		{
			dist2=TM_tmp1[i].distance_square(mol2[j]);
			TM_DP_sco[cur_index+j]=1.0/(1.0+dist2/d02);  // this is NOT d0 !!!
			if(TM_sse1[i]==TM_sse2[j])TM_DP_sco[cur_index+j]+=0.5;
		}
	}
	//apply dynamic programming
	double gap_open=-1.0;
	TM_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,ali2,gap_open);
}

/**************************************************************
*     get initial alignment invmap0(i) from fragment gapless threading
**************************************************************/
void TM_Align_Main::TM_Get_Initial4(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	//init
	int fra_min=4;           //  >=4,minimum fragment for search
	int fra_min1=fra_min-1;  //  cutoff for shift, save time
	double dcu0=3.85;

//[1]  Find the maximal continuous fragments -------->
	int mol1max,mol2max;
	int i,j,k;
	int nfr;
	int r_min;
	int Lfr_max;
	double dcu;
	double dis;
	//-> for mol1
	dcu=dcu0;
	r_min=(int)(moln1/3.0);
	if(r_min>fra_min)r_min=fra_min;
ws_recur1:
	nfr=1;                    // number of fragments
	TM_frag1[(nfr-1)*2+0]=0;  // what residue
	TM_frag1[(nfr-1)*2+1]=1;  // length of the fragment
	for(i=1;i<moln1;i++)
	{
		dis=mol1[i-1].distance(mol1[i]);
		if(dis<dcu)TM_frag1[(nfr-1)*2+1]++;
		else
		{
			nfr++;
			TM_frag1[(nfr-1)*2+0]=i;
			TM_frag1[(nfr-1)*2+1]=1;
		}
	}
	Lfr_max=0;
	mol1max=0;
	for(k=0;k<nfr;k++)
	{
		if(Lfr_max<TM_frag1[k*2+1])
		{
			Lfr_max=TM_frag1[k*2+1];
			mol1max=k;
		}
	}
	if(Lfr_max<r_min)
	{
		dcu*=1.1;
		goto ws_recur1;
	}
	//-> for mol2
	dcu=dcu0;
	r_min=(int)(moln2/3.0);
	if(r_min>fra_min)r_min=fra_min;
ws_recur2:
	nfr=1;                    // number of fragments
	TM_frag2[(nfr-1)*2+0]=0;  // what residue
	TM_frag2[(nfr-1)*2+1]=1;  // length of the fragment
	for(i=1;i<moln2;i++)
	{
		dis=mol2[i-1].distance(mol2[i]);
		if(dis<dcu)TM_frag2[(nfr-1)*2+1]++;
		else
		{
			nfr++;
			TM_frag2[(nfr-1)*2+0]=i;
			TM_frag2[(nfr-1)*2+1]=1;
		}
	}
	Lfr_max=0;
	mol2max=0;
	for(k=0;k<nfr;k++)
	{
		if(Lfr_max<TM_frag2[k*2+1])
		{
			Lfr_max=TM_frag2[k*2+1];
			mol2max=k;
		}
	}
	if(Lfr_max<r_min)
	{
		dcu*=1.1;
		goto ws_recur2;
	}

//[2]  select what piece will be used -------->
//   select what piece will be used (this may araise ansysmetry, but
//   only when L1=L2 and Lfr1=Lfr2 and L1 ne Lfr1
//   if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1
	int mark;
	if(TM_frag1[mol1max*2+1]<TM_frag2[mol2max*2+1])mark=1;
	else if(TM_frag1[mol1max*2+1]>TM_frag2[mol2max*2+1])mark=2;
	else
	{
		if(moln1<=moln2)mark=1;
		else mark=2;
	}
	int moln;
	int pos;
	int L_fr;
	//get_correspondence
	if(mark==1)  //using mol1
	{
		moln=moln1;
		pos=TM_frag1[mol1max*2+0];
		L_fr=TM_frag1[mol1max*2+1];
		for(k=0;k<L_fr;k++)TM_DP_best[k]=pos+k;
	}
	else         //using mol2
	{
		moln=moln2;
		pos=TM_frag2[mol2max*2+0];
		L_fr=TM_frag2[mol2max*2+1];
		for(k=0;k<L_fr;k++)TM_DP_best[k]=pos+k;
	}
	//final judge
	int n1,n2;
	if(L_fr==moln)
	{
		n1=(int)(moln*0.1);
		n2=(int)(moln*0.89);
		j=0;
		for(i=n1;i<=n2;i++)
		{
			TM_DP_best[j]=TM_DP_best[n1+j];
			j++;
		}
		L_fr=j;
	}

//[3]  get initial ------------->
	int moln1_,moln2_;
	int aL;
	int idel;
	int ishift;
	int count;
	double GL_max=0.0;
	double GL_score;
	if(mark==1)           //using mol1 as template
	{
		moln1_=L_fr;
		aL=moln1_<moln2?moln1_:moln2;
		idel=(int)(aL/2.5);    // minimum size of considered fragment
		if(idel<=fra_min1)idel=fra_min1;
		n1=-moln2+idel; //shift1
		n2=moln1_-idel; //shift2
		for(ishift=n1;ishift<=n2;ishift++)
		{
			count=0;
			for(j=0;j<moln2;j++)
			{
				i=j+ishift;
				if(i>=0&&i<moln1_)
				{
					count++;
					TM_DP_ali2[j]=TM_DP_best[i];
				}
				else TM_DP_ali2[j]=-1;
			}
			if(count>=idel)
			{
				GL_score=TM_Align_Get_GL(mol1,mol2,moln1,moln2,TM_DP_ali2);
				if(GL_score>GL_max)
				{
					GL_max=GL_score;
					EqualArray(ali2,TM_DP_ali2,moln2);
				}
			}
		}
	}
	else                  //using mol2 as template
	{
		moln2_=L_fr;
		aL=moln1<moln2_?moln1:moln2_;
		idel=(int)(aL/2.5);    // minimum size of considered fragment
		if(idel<=fra_min1)idel=fra_min1;
		n1=-moln2_+idel; //shift1
		n2=moln1-idel; //shift2
		for(ishift=n1;ishift<=n2;ishift++)
		{
			count=0;
			for(j=0;j<moln2;j++)TM_DP_ali2[j]=-1;
			for(j=0;j<moln2_;j++)
			{
				i=j+ishift;
				if(i>=0&&i<moln1)
				{
					count++;
					TM_DP_ali2[TM_DP_best[j]]=i;
				}
			}
			if(count>=idel)
			{
				GL_score=TM_Align_Get_GL(mol1,mol2,moln1,moln2,TM_DP_ali2);
				if(GL_score>GL_max)
				{
					GL_max=GL_score;
					EqualArray(ali2,TM_DP_ali2,moln2);
				}
			}
		}
	}
}

/**************************************************************
*    fifth initial alignement. Using local structure super   *
*    position.                                               *
**************************************************************/
void TM_Align_Main::TM_Get_Initial5(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
//[1]-> setting parameters ************************************
	int n_frag=20;
	if(moln1>250&&moln2>250)n_frag=50;
	else if(moln1>200&&moln2>200)n_frag=40;
	else if(moln1>150&&moln2>150)n_frag=30;
	else n_frag=20;

//[2]-> next stage --------//
	int n_jump=n_frag;
	int smaller=moln1<moln2?moln1:moln2;
	if(TM_CALC==0)Calc_TM_d0(smaller);
	double d01=d0+1.5;
	if(d01<0.5)d01=0.5;
	double GL_max=0.0;
	double GL_score;
	int m1=moln1-n_frag-20;
	int m2=moln2-n_frag-20;

//[3]-> real process -------//
	int k;
	int ii,jj;
	int iii,jjj;
	double gap_open=0.0;            // this is NOT defined in TM_align.f !!! //__Added at//__101205__//
	for(ii=20;ii<=m1;ii+=n_jump)
	{
		for(jj=20;jj<=m2;jj+=n_jump)
		{
			//get correspondece
			for(k=0;k<n_frag;k++)
			{
				iii=ii+k-1;
				jjj=jj+k-1;
				TM_tmp1[k]=mol1[iii];
				TM_tmp2[k]=mol2[jjj];
			}
			kabsch(TM_tmp2,TM_tmp1,n_frag,finmat);
			rot_mol(mol1,TM_tmp1,moln1,finmat);
			//calculate score matrix
//			TM_Align_Get_Matrix_TMs(TM_tmp1,mol2,moln1,moln2,d01,TM_DP_sco); //this is NOT d0!!!
			TM_Align_Get_Matrix(TM_tmp1,mol2,moln1,moln2,d01,TM_DP_sco); //this is NOT d0!!!
			//apply dynamic programming
			TM_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,TM_DP_ali2,gap_open);
			GL_score=TM_Align_Get_GL(mol1,mol2,moln1,moln2,TM_DP_ali2);
			if(GL_score>GL_max)
			{
				GL_max=GL_score;
				EqualArray(ali2,TM_DP_ali2,moln2);
			}            
		}
	}
}
void TM_Align_Main::TM_Get_Initial6(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
//[1]-> setting parameters ************************************
	int n_frag=5;
	int n_jump=n_frag;
	double n_all_max=0.0;
	double n_max_cut=0.7;
	double n_cutt=0.0;

//[2]-> next stage --------//
	int smaller=moln1<moln2?moln1:moln2;
	if(TM_CALC==0)Calc_TM_d0(smaller);
	double d01=d0+1.5;
	if(d01<0.5)d01=0.5;
	double GL_score;
	double GL_maxA=0.0;
	double GL_maxA_cut=0.12;
	int n_frag2=n_frag/2;
	int m1=moln1-n_frag+1;
	int m2=moln2-n_frag+1;

//[3]-> real process -------//
	int k;
	int ii,jj;
	int iii,jjj;
	int secmatch;
	int coilcnt;
	double gap_open=0.0;            // this is NOT defined in TM_align.f !!! //__Added at//__101205__//
	int align_len1;
	for(ii=1;ii<=m1;ii+=n_jump)
	{
		for(jj=1;jj<=m2;jj+=n_jump)
		{
			secmatch=0;
			coilcnt=0;
			for(k=0;k<n_frag;k++)
			{
				iii=ii+k-1;
				jjj=jj+k-1;
				if(TM_sse1[iii]==TM_sse2[jjj])secmatch++;
				if(TM_sse1[iii]==1)coilcnt++;
			}
			if(secmatch>n_frag2 && coilcnt<n_frag2)
			{
				//get correspondece
				for(k=0;k<moln2;k++)TM_DP_ali2[k]=-1;
				for(k=0;k<n_frag;k++)
				{
					iii=ii+k-1;
					jjj=jj+k-1;
					TM_DP_ali2[jjj]=iii;
				}
				GL_score=TM_Align_Get_GL(mol1,mol2,moln1,moln2,TM_DP_ali2);
				if(GL_score>GL_maxA*GL_maxA_cut)
				{
					//get correspondece
					for(k=0;k<n_frag;k++)
					{
						iii=ii+k-1;
						jjj=jj+k-1;
						TM_tmp1[k]=mol1[iii];
						TM_tmp2[k]=mol2[jjj];
					}
					kabsch(TM_tmp2,TM_tmp1,n_frag,finmat);
					rot_mol(mol1,TM_tmp1,moln1,finmat);
					//calculate score matrix
//					TM_Align_Get_Matrix_TMs(TM_tmp1,mol2,moln1,moln2,d01,TM_DP_sco); // this is NOT d0 !!!
					TM_Align_Get_Matrix(TM_tmp1,mol2,moln1,moln2,d01,TM_DP_sco); //this is NOT d0!!!
					//apply dynamic programming 
					align_len1=TM_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,TM_DP_best,gap_open);
					n_cutt=n_all_max*n_max_cut;
					if(align_len1>n_cutt)
					{
						if(align_len1>n_all_max)n_all_max=align_len1;
						GL_score=TM_Align_Get_GL(mol1,mol2,moln1,moln2,TM_DP_best);
						if(GL_score>GL_maxA)
						{
							GL_maxA=GL_score;
							EqualArray(ali2,TM_DP_best,moln2);
						}
					}
				}//end of IF(GL_score>GL_max*GL_max_cut)
			}//end of IF(secmatch>n_frag2 && coilcnt<n_frag2)
		}//end of FOR(jj)
	}//end of FOR(ii)
}

//============================= the following is the main function for TM_align =======================//
//[original, according to Zhang Yang's TMalign.f (2009) + TMalign.f (2010)]
double TM_Align_Main::TM_Align_Total(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,
									 int norm_len,double norm_d0,double *MAXSCO)
{
	int k;
	double TM_CUR;
	double TM_BEST=0.0;

	//total_init
	AliT_Rec=0;  //default:NOT
	AliT_Max=0;  //default:NO
	for(k=0;k<moln2;k++)ali2[k]=-1;
	for(k=0;k<moln2;k++)TM_Main_Ali2[k]=-1;

	//--- get_initial1 ---//gapless threading
//ws_init1:
	TM_Get_Initial1(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=1;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init2;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=1;
	}
	//--- get_initial1.5 ---//gapless threading advance
	if(TM_ADDITION==1 && AliT_Rec==1)
	{
		EqualArray(TM_Main_Ali2,TM_Main_AliT,moln2);
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			EqualArray(ali2,TM_Main_Ali2,moln2);
			EqualArray(TM_FINMAT,TM_rotmat,12);
			AliT_Max=7;
		}
		if(TM_CUR<=TM_BEST*0.2)goto ws_init2;  //check bad...
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1,2); //run DynaProg for only 2 cycles!!
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			EqualArray(ali2,TM_Main_Ali2,moln2);
			EqualArray(TM_FINMAT,TM_rotmat,12);
			AliT_Max=7;
		}
	}
	//--- get_initial2 ---//SSE
ws_init2:
	TM_Get_Initial2(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=2;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init5;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=2;
	}
	//--- get_initial5 ---//local structure superposition(neo)
ws_init5:
	TM_Get_Initial5(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=5;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init6;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1,2); //run DynaProg for only 2 cycles!!
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=5;
	}
	//--- get_initial6 ---//local structure superposition(old)
ws_init6:
	if(TM_ADDITION==1)
	{
		TM_Get_Initial6(mol1,mol2,moln1,moln2,TM_Main_Ali2);
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			EqualArray(ali2,TM_Main_Ali2,moln2);
			EqualArray(TM_FINMAT,TM_rotmat,12);
			AliT_Max=6;
		}
		if(TM_CUR<=TM_BEST*0.2)goto ws_init3;  //check bad...
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1,2); //run DynaProg for only 2 cycles!!
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			EqualArray(ali2,TM_Main_Ali2,moln2);
			EqualArray(TM_FINMAT,TM_rotmat,12);
			AliT_Max=6;
		}
	}
	//--- get_initial3 ---//best+SSE
ws_init3:
	EqualArray(TM_Main_Ali2,ali2,moln2);  //assign the current best path!!
	TM_Get_Initial3(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=3;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init4;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=3;
	}
	//--- get_initial4 ---//
ws_init4:
	TM_Get_Initial4(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=4;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_end;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1,30,0); //run DynaProg, no skip!!!
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=4;
	}
	//========= terminal =========//
ws_end:
	TM_FIN_TMS=TM_Align_TM_Score(mol1,mol2,moln1,moln2,ali2,
		norm_len,norm_d0,TM_FIN_RMS,TM_FIN_LALI,MAXSCO);
	EqualArray(TM_FINMAT,finmat,12);
	return TM_FIN_TMS;
}

//[modified, just change the order of the initial alignment]
double TM_Align_Main::TM_Align_Total_II(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,
										int norm_len,double norm_d0,double *MAXSCO)
{
	int k;
	double TM_CUR;
	double TM_BEST=0.0;

	//total_init
	AliT_Rec=0;  //default:NOT
	AliT_Max=0;  //default:NO
	for(k=0;k<moln2;k++)ali2[k]=-1;
	for(k=0;k<moln2;k++)TM_Main_Ali2[k]=-1;

	//--- get_initial5 ---//local structure superposition(neo)
//ws_init5:
	TM_Get_Initial5(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=5;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init6;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=5;
	}

	//--- get_initial6 ---//local structure superposition(old)
ws_init6:
	if(TM_ADDITION==1)
	{
		TM_Get_Initial6(mol1,mol2,moln1,moln2,TM_Main_Ali2);
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			EqualArray(ali2,TM_Main_Ali2,moln2);
			EqualArray(TM_FINMAT,TM_rotmat,12);
			AliT_Max=6;
		}
		if(TM_CUR<=TM_BEST*0.2)goto ws_init2;  //check bad...
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			EqualArray(ali2,TM_Main_Ali2,moln2);
			EqualArray(TM_FINMAT,TM_rotmat,12);
			AliT_Max=6;
		}
	}

	//--- get_initial2 ---//SSE
ws_init2:
	TM_Get_Initial2(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=2;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init1;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=2;
	}

	//--- get_initial1 ---//gapless threading
ws_init1:
	TM_Get_Initial1(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=1;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init4;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=1;
	}

	//--- get_initial1.5 ---//gapless threading advance
	if(TM_ADDITION==1 && AliT_Rec==1)
	{
		EqualArray(TM_Main_Ali2,TM_Main_AliT,moln2);
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			EqualArray(ali2,TM_Main_Ali2,moln2);
			EqualArray(TM_FINMAT,TM_rotmat,12);
			AliT_Max=7;
		}
		if(TM_CUR<=TM_BEST*0.2)goto ws_init4;  //check bad...
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			EqualArray(ali2,TM_Main_Ali2,moln2);
			EqualArray(TM_FINMAT,TM_rotmat,12);
			AliT_Max=7;
		}
	}

	//--- get_initial4 ---//
ws_init4:
	TM_Get_Initial4(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=4;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init3;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=4;
	}

	//--- get_initial3 ---//best+SSE
ws_init3:
	EqualArray(TM_Main_Ali2,ali2,moln2);  //assign the current best path!!
	TM_Get_Initial3(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=3;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_end;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=3;
	}

	//========= terminal =========//
ws_end:
	TM_FIN_TMS=TM_Align_TM_Score(mol1,mol2,moln1,moln2,ali2,
		norm_len,norm_d0,TM_FIN_RMS,TM_FIN_LALI,MAXSCO);
	EqualArray(TM_FINMAT,finmat,12);
	return TM_FIN_TMS;
}
