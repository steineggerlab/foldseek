#include "TM_score.h"

//------ constructor -------//
TM_score::TM_score(int num)
{
	//init
	TMS_maximal=num;
	Init_TM_Score(TMS_maximal);
	//default
	d0_min=0.5;
	d0=5.0;
	d8=10.0;
	d00=5.0;
	//macros
	TM_CALC=0; //when dealing with a new structure-pair, this value MUST be set to ZERO !!
	TM_CACHE=0; //no cache
}
TM_score::~TM_score(void)
{
	Dele_TM_Score();
}

//--------- init ---------//
void TM_score::Init_TM_Score(int maxlen)
{
	TMs_molt=new XYZ[maxlen];
	TMs_cache=new int[maxlen];
	tmp1=new XYZ[maxlen];
	tmp2=new XYZ[maxlen];
	i_ali=new int[maxlen];
	L_ini=new int[maxlen];
	k_ali=new int[maxlen];
	k_fin=new int[maxlen];
	iL0_simp=new int[maxlen];
	rotmat=new double[12];
	finmat=new double[12];
}
void TM_score::Dele_TM_Score(void)
{
	delete [] TMs_molt;
	delete [] TMs_cache;
	delete [] tmp1;
	delete [] tmp2;
	delete [] i_ali;
	delete [] L_ini;
	delete [] k_ali;
	delete [] k_fin;
	delete [] iL0_simp;
	delete [] rotmat;
	delete [] finmat;
}
void TM_score::Input_TM_Score(XYZ *mol1,XYZ *mol2,int lali,double d0,double d8)
{
	ori1=mol1;
	ori2=mol2;
	ali_orin=lali;
	d0_input=d0;
	d8_input=d8;
}
//========= Cache_Point =========//
void TM_score::TMs_Cache_Point(XYZ *mol,int index,XYZ &ret_point,double *rotmat)
{
	if(TM_CACHE==1 && TMs_cache[index]==1)ret_point=TMs_molt[index];
	else
	{
		rot_point(mol[index],ret_point,rotmat);
		TMs_molt[index]=ret_point;
		TMs_cache[index]=1;
	}
}

//--------- calc d0 -------//
double TM_score::Calc_TM_d0(double len) //init
{
/*
      d0_min=0.5
      if(m_d0_min.eq.1)then
         d0_min=d0_min_input    !for search
      endif
      anseq_min=min(nseq1,nseq2)
      anseq=anseq_min           !length for defining TMscore in search
      d8=1.5*anseq_min**0.3+3.5 !remove pairs with dis>d8 during search & final
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      if(m_d0.eq.1)d0=d0_fix
      d00=d0                    !for quickly calculate TM-score in searching
      if(d00.gt.8)d00=8
      if(d00.lt.4.5)d00=4.5
*/
	double T_Len;
	double T_Len8;
	double T_LLen;
	double T_Min;
	//calculate
	T_Min=0.5;
	T_Len8=1.5*exp(0.3*log(len))+3.5;
	if(len>15)T_Len=1.24*exp(1.0/3.0*log(len-15.0))-1.8;
	else T_Len=T_Min;
	if(T_Len<T_Min)T_Len=T_Min;
	T_LLen=T_Len;
	if(T_LLen>8.0)T_LLen=8.0;
	if(T_LLen<4.5)T_LLen=4.5;
	//evaluate
	d0_min=T_Min;
	d0=T_Len;
	d8=T_Len8;
	d00=T_LLen;
	//return
//	TM_CALC=1;
	return d0;
}
double TM_score::Calc_TM_d0_Simp(double len) //init
{
	double T_Len;
	double T_Min;
	//calculate
	T_Min=0.5;
	if(len>15)T_Len=1.24*exp(1.0/3.0*log(len-15.0))-1.8;
	else T_Len=T_Min;
	if(T_Len<T_Min)T_Len=T_Min;
	//return
	return T_Len;
}

//--------- minor ---------//
//note: mol2 is fixed!! While mol1 is aligned to mol2
int TM_score::Calc_Small_Dist(XYZ *mol1,XYZ *mol2,int lali,double *rotmat_,double d0,int *ali_tmp)
{
	int i;
	double dist2;
	double ori_d0=d0;
	double ori_d;
	XYZ temp;
	int n_cut;

tag_21:
	n_cut=0;
	ori_d=ori_d0*ori_d0;
	for(i=0;i<lali;i++)
	{
		TMs_Cache_Point(mol1,i,temp,rotmat_);
		dist2=mol2[i].distance_square(temp);
		if(dist2<ori_d)
		{
			ali_tmp[n_cut]=i;
			n_cut++;
		}
	}
	if(n_cut<3 && lali>3)
	{
		ori_d0+=0.5;
		goto tag_21;
	}
	return n_cut;
}
double TM_score::Calc_TM_Score_Single(XYZ *mol1,XYZ *mol2,int lali,double *rotmat_,double d0,double d8,int TM8orTM,
									  double *ALLSCO)
{
	int i;
	double dist2;
	double ori_d;
	double ori_d8;
	XYZ temp;
	double score=0.0;
	ori_d=d0*d0;
	ori_d8=d8*d8;

	// additional score //__110230__//
	if(ALLSCO!=NULL) //8-digit!!
	{
		//GDT_score
		ALLSCO[0]=0;   //number of residues<0.5
		ALLSCO[1]=0;   //number of residues<1.0
		ALLSCO[2]=0;   //number of residues<2.0
		ALLSCO[3]=0;   //number of residues<4.0
		ALLSCO[4]=0;   //number of residues<8.0
		//others
		ALLSCO[5]=0;   //maxsub score
		ALLSCO[6]=0;   //TMscore
		ALLSCO[7]=0;   //TMscore10
	}
	// additional score //__110230__//over

	for(i=0;i<lali;i++)
	{
		TMs_Cache_Point(mol1,i,temp,rotmat_);
		dist2=mol2[i].distance_square(temp);

		// additional score //__110230__//
		if(ALLSCO!=NULL)
		{
			//GDT_score
			if(dist2<=64.0)
			{
				ALLSCO[4]++;
				if(dist2<=16.0)
				{
					ALLSCO[3]++;
					if(dist2<=4.0)
					{
						ALLSCO[2]++;
						if(dist2<=1.0)
						{
							ALLSCO[1]++;
							if(dist2<=0.25)
							{
								ALLSCO[0]++;
							}
						}
					}
				}
			}
			//MAXsub_score
			if(dist2<12.25)
			{
				ALLSCO[5]+=1.0/(1.0+dist2/12.25);
			}
			//TM_score
			ALLSCO[6]+=1.0/(1.0+dist2/ori_d);
			if(dist2<100.0)
			{
				ALLSCO[7]+=1.0/(1.0+dist2/ori_d);
			}
		}
		// additional score //__110230__//over

		if(TM8orTM==1)
		{
			if(dist2>ori_d8)continue;
		}
		score+=1.0/(1.0+dist2/ori_d);
	}
	return score;
}
double TM_score::Calc_TM_Score_Simple(XYZ *mol1,XYZ *mol2,int lali,double d0,double d8,int TM8orTM,
									  double *ALLSCO)
{
	int i;
	double dist2;
	double ori_d;
	double ori_d8;
	double score=0.0;
	ori_d=d0*d0;
	ori_d8=d8*d8;

	// additional score //__110230__//
	if(ALLSCO!=NULL) //8-digit!!
	{
		//GDT_score
		ALLSCO[0]=0;   //number of residues<0.5
		ALLSCO[1]=0;   //number of residues<1.0
		ALLSCO[2]=0;   //number of residues<2.0
		ALLSCO[3]=0;   //number of residues<4.0
		ALLSCO[4]=0;   //number of residues<8.0
		//others
		ALLSCO[5]=0;   //maxsub score
		ALLSCO[6]=0;   //TMscore
		ALLSCO[7]=0;   //TMscore10
	}
	// additional score //__110230__//over

	for(i=0;i<lali;i++)
	{
		dist2=mol2[i].distance_square(mol1[i]);

		// additional score //__110230__//
		if(ALLSCO!=NULL)
		{
			//GDT_score
			if(dist2<=64.0)
			{
				ALLSCO[4]++;
				if(dist2<=16.0)
				{
					ALLSCO[3]++;
					if(dist2<=4.0)
					{
						ALLSCO[2]++;
						if(dist2<=1.0)
						{
							ALLSCO[1]++;
							if(dist2<=0.25)
							{
								ALLSCO[0]++;
							}
						}
					}
				}
			}
			//MAXsub_score
			if(dist2<12.25)
			{
				ALLSCO[5]+=1.0/(1.0+dist2/12.25);
			}
			//TM_score
			ALLSCO[6]+=1.0/(1.0+dist2/ori_d);
			if(dist2<100.0)
			{
				ALLSCO[7]+=1.0/(1.0+dist2/ori_d);
			}
		}
		// additional score //__110230__//over

		if(TM8orTM==1)
		{
			if(dist2>ori_d8)continue;
		}
		score+=1.0/(1.0+dist2/ori_d);
	}
	return score;
}

//---------- main -----------//
//mol2 is fixed, mol1 superimpose onto mol2, to make the TMscore best, given the correspondence
double TM_score::Calc_TM_Score(XYZ *mol1,XYZ *mol2,int lali,double d0,double d8,int SEARCHorNOT,int TM8orTM,
							   double *MAXSCO)
{
	//init
	Input_TM_Score(mol1,mol2,lali,d0,d8);
	if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*lali);

	// additional score //__110230__//
	double *ALLSCO=NULL;
	if(MAXSCO!=NULL) //must 8-digit !! (5gdt+1msb+2tms)
	{
		//allsco create
		ALLSCO=new double[8];
		//maximal init
		for(int wwi=0;wwi<8;wwi++)MAXSCO[wwi]=-1;
	}
	// additional score //__110230__//over

	//parameter
	//[1]d0 related
	if(d0_input<0.5)d0_input=0.5;
	d0_serch=d0_input;
	if(d0_serch>8.0)d0_serch=8.0;
	if(d0_serch<4.5)d0_serch=4.5;
	//[2]iteration related
	int i;
	int n_it=20;          //maximum number of iterations
	int n_init;
	int n_init_max=6;     //maximum number of L_init
	int L_ini_min=4;
	if(ali_orin<4)L_ini_min=ali_orin;
	for(i=0;i<n_init_max-1;i++)
	{
		L_ini[i]=(int)(1.0*ali_orin/(exp(1.0*i*log(2.0))));
		if(L_ini[i]<=L_ini_min)
		{
			L_ini[i]=L_ini_min;
			i++;
			goto tag_402;
		}
	}
	L_ini[i]=L_ini_min;
	i++;
tag_402:
	n_init=i;

	//find the maximum score starting from local structures superposition
	double ddd;
	double score;
	double score_max=-1.0;
	int i_init;
	int L_init;
	int iL;
	int iL1;
	int iL_max;
	int ka;
	int k;
	int ws;
	int n_shift;
	int it;
	int neq;
	for(i_init=0;i_init<n_init;i_init++)
	{
		L_init=L_ini[i_init];
		iL_max=ali_orin-L_init+1;

		//this is simplication//__101115__//
		if(SEARCHorNOT==1)
		{
			ws=0;
			for(k=1;k<=iL_max;k+=40)
			{
				ws++;
				iL0_simp[ws-1]=k;
			}
			if(iL0_simp[ws-1]<iL_max)
			{
				ws++;
				iL0_simp[ws-1]=iL_max;
			}
			n_shift=ws;
		}
		else
		{
			ws=0;
			for(k=1;k<=iL_max;k++)
			{
				ws++;
				iL0_simp[ws-1]=k;
			}
			if(iL0_simp[ws-1]<iL_max)
			{
				ws++;
				iL0_simp[ws-1]=iL_max;
			}
			n_shift=ws;
		}
		//this is simplication//__101115__//over

		for(iL=0;iL<n_shift;iL++)    //on aligned residues, [1,nseqA]
		{
			iL1=iL0_simp[iL]-1;
			//get initial superposition -------------------------------->
			for(i=0;i<L_init;i++)
			{
				k=iL1+i;
				tmp1[i]=ori1[k];
				tmp2[i]=ori2[k];
				k_ali[i]=k;
			}
			ka=L_init;
			kabsch(tmp2,tmp1,ka,rotmat);
			if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*lali);
			ddd=d0_serch-1.0;
			score=Calc_TM_Score_Single(ori1,ori2,ali_orin,rotmat,d0_input,d8_input,TM8orTM,ALLSCO);
			n_cut=Calc_Small_Dist(ori1,ori2,ali_orin,rotmat,ddd,i_ali);
			if(score_max<score)
			{
				score_max=score;
				k_alin=ka;
				memcpy(k_fin,k_ali,k_alin*sizeof(int));
				memcpy(finmat,rotmat,12*sizeof(double));
			}

			// additional score //__110230__//
			if(MAXSCO!=NULL)
			{
				for(int wwi=0;wwi<8;wwi++)if(MAXSCO[wwi]<ALLSCO[wwi])MAXSCO[wwi]=ALLSCO[wwi];
			}
			// additional score //__110230__//over
            
			//iteration for extending ---------------------------------->
			ddd=d0_serch+1.0;
			for(it=0;it<n_it;it++)
			{
				for(i=0;i<n_cut;i++)
				{
					k=i_ali[i];
					tmp1[i]=ori1[k];
					tmp2[i]=ori2[k];
					k_ali[i]=k;
				}
				ka=n_cut;
				kabsch(tmp2,tmp1,ka,rotmat);
				if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*lali);
				score=Calc_TM_Score_Single(ori1,ori2,ali_orin,rotmat,d0_input,d8_input,TM8orTM,ALLSCO);
				n_cut=Calc_Small_Dist(ori1,ori2,ali_orin,rotmat,ddd,i_ali);
				if(score_max<score)
				{
					score_max=score;
					k_alin=ka;
					memcpy(k_fin,k_ali,k_alin*sizeof(int));
					memcpy(finmat,rotmat,12*sizeof(double));
				}

				// additional score //__110230__//
				if(MAXSCO!=NULL)
				{
					for(int wwi=0;wwi<8;wwi++)if(MAXSCO[wwi]<ALLSCO[wwi])MAXSCO[wwi]=ALLSCO[wwi];
				}
				// additional score //__110230__//over

				//final judgement -------------------------------------->
				if(it==n_it-1)break;
				if(n_cut==ka)
				{
					neq=0;
					for(i=0;i<n_cut;i++)if(i_ali[i]==k_ali[i])neq++;
					if(n_cut==neq)break;
				}
			}//end of extending
		}//end of aligned residue
	}//end of n_init


	// additional score //__110230__//
	if(MAXSCO!=NULL)delete [] ALLSCO;
	// additional score //__110230__//over

	//final_return
	return score_max;
}
