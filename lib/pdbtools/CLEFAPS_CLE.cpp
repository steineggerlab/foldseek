#include "CLEFAPS_CLE.h"

//================= constructor ================//
CLEFAPS_CLE::CLEFAPS_CLE(int num,int CLESUM)
:Ali_Ori(num),Ali_AFP(num),TM_align(num),
 Bioinfo_Code(CLESUM),CLEFAPS(num)
{
	//init
	IN_DUM1=new int[num];
	IN_DUM2=new int[num];
	//main parameters
	SFP_Strategy=1; //CLE
	REF_Strategy=1; //refine
	SFP_L_Len=6;
	SFP_H_Len=9;
	SFP_L_Thres=0;
	SFP_H_Thres=10;
	//refine parameter
	Refine_Cut=0.3; //only TM_refine 0.3*max
	//upper and lower max
	ZM_Upper_Max=500;
	ZM_Lower_Max=100;
	ZM_TopL=10;
}
CLEFAPS_CLE::~CLEFAPS_CLE(void)
{
	//delete
	delete [] IN_DUM1;
	delete [] IN_DUM2;
}

//=================== input functions ============//
void CLEFAPS_CLE::CLEFAPS_Input_Func(char *c1,char *c2,char *a1,char *a2,
									 XYZ *m1,XYZ *m2,int n1,int n2)
{
	//copy pointer
	IN_CLE1=c1;
	IN_CLE2=c2;
	IN_AMI1=a1;
	IN_AMI2=a2;
	//copy data
	mol1=m1;
	mol2=m2;
	moln1=n1;
	moln2=n2;
	//mask check for CLE
	if(mas1!=0 && mas2!=0)
	{
		//-> mask the break point 
		Mask_CLE_BreakPoint(IN_CLE1,mas1,moln1);
		Mask_CLE_BreakPoint(IN_CLE2,mas2,moln2);
	}
}
void CLEFAPS_CLE::CLEFAPS_Input_Para(int h_len,int l_len,int h_thres,int l_thres,int Stragety)
{
	SFP_H_Len=h_len;
	SFP_L_Len=l_len;
	SFP_H_Thres=h_thres;
	SFP_L_Thres=l_thres;
	SFP_Strategy=Stragety;
}

//---------- process functions -----------//__for Hyper_CLEFAPS Only !!! __//__110408__//
void CLEFAPS_CLE::CLEFAPS_Second_Init(void)
{
	//macros
	TMSorNOT=0;               //no TMscore refine!!
	HEAD_Chk=0;               //no head check!!
	COMP_Chk=0;               //no comparison check!!
	//TopJ init  -> 1st level cut (From SFP_H -> Top_Pivot)
	ZM_TopJ=(int)(moln1*moln2/(SFP_H_Len*SFP_H_Len*6.25));
	if(ZM_TopJ>ZM_Upper_Max)ZM_TopJ=ZM_Upper_Max;  //this must be killed later!! //__110930__//
	if(ZM_TopJ<ZM_Lower_Max)ZM_TopJ=ZM_Lower_Max;  //this could be modified in the future!! //__110930__//
	//TopK init  -> 2nd level cut (From Top_Pivot -> Candidate)
//	ZM_TopK=20;    //this could be modified in the future...
	//TopL init  -> 3rd level cut (From Candidate -> Fin_Alignment)
//	ZM_TopL=10;    //this could be modified in the future...
}
//suppose the input "tot" has already be sorted
double CLEFAPS_CLE::CLEFAPS_Second_Refine(vector <Align_Record> &tot)
{
	//init_check
	int size=(int)tot.size();
	if(size==0)return -1.0;
	//re-add
	vector <Align_Record> tmp;
	tmp.clear();
	int i,j;
	int cursize;
	double score;
	int success;
	int count=0;
	for(i=0;i<size;i++)
	{
		success=1; //defualt:1
		cursize=(int)tmp.size();
		for(j=0;j<cursize;j++)
		{
			score=Alg_Comp_Simp(tot[i].alignment,tmp[j].alignment,Comp_Check_Vari);
			if(score>Comp_Check)
			{
				success=0;
				break;
			}
		}
		if(success==1)
		{
			tot[i].index.push_back(count);
			tmp.push_back(tot[i]);
			count++;
		}
	}
	//re-refine (using TMscore optimization)
	int k;
	int relnum=(int)tmp.size();
	int totnum;
	totnum=ZM_TopL;
	totnum=totnum<relnum?totnum:relnum;
	double wssco;
	double tmsco;
	double wsmax=0.0;
	double tmmax=0.0;
	int wstot=0;
	for(j=0;j<totnum;j++)
	{
		for(k=0;k<moln2;k++)ali2[k]=tmp[j].alignment[k];
		wssco=Calc_TM_Align(mol1,mol2,moln1,moln2,ali2,ali2,0,REFINE_ITER); //apply TMscore dynamic programming
		EqualArray(CLEP_rotmat,TM_rotmat,12);
		tmsco=CLEPAPS_Make_Score(tmp[j]);
		tmp[j].main_sco=wssco;
		if(wssco>wsmax)wsmax=wssco;
		if(tmsco>tmmax)tmmax=tmsco;
		if(wssco<wsmax*Refine_Cut && tmsco<tmmax*Refine_Cut)break;
		wstot++;
	}
	//final-add
	tmp.resize(wstot);
	stable_sort(tmp.rbegin(),tmp.rend());   //descending sort
	tot.clear();
	for(i=0;i<wstot;i++)
	{
		success=1; //defualt:1
		cursize=(int)tot.size();
		for(j=0;j<cursize;j++)
		{
			score=Alg_Comp_Simp(tmp[i].alignment,tot[j].alignment,Comp_Check_Vari);
			if(score>Comp_Check)
			{
				success=0;
				break;
			}
		}
		if(success==1)tot.push_back(tmp[i]);
	}
	//final
	CLEP_MAX=tot[0];
	return CLEP_MAX.main_sco;
}

//final TMscore refine
void CLEFAPS_CLE::CLEFAPS_Second_Final(vector <Align_Record> &tot,int REAL_FAST)
{
	//init_check
	int i,k;
	int size=(int)tot.size();
	if(size==0)return;
	//refine
	double wssco;
	for(i=0;i<size;i++)
	{
		for(k=0;k<moln2;k++)ali2[k]=tot[i].alignment[k];
		if(REAL_FAST==0) //REAL
		{
//			wssco=TM_Align_TM_Score(mol1,mol2,moln1,moln2,ali2,RMSD,LALI);      //using REAL version (this could be modified!!)
			wssco=TM_Align_TM_Score_Simp(mol1,mol2,moln1,moln2,ali2,RMSD,LALI); //using FAST version (this could be modified!!)
		}
		else             //FAST
		{
			wssco=TM_Align_TM_Score_Simp(mol1,mol2,moln1,moln2,ali2,RMSD,LALI); //using FAST version (this could be modified!!)
		}
		EqualArray(CLEP_rotmat,finmat,12);
		CLEPAPS_Make_Score(tot[i]);
		tot[i].main_sco=wssco;
	}
}

//================ main functions ============//
void CLEFAPS_CLE::CLEFAPS_Main_Init(int SFP_Strategy)
{
	//generate dummy
	switch(SFP_Strategy) 
	{ 
		case 1:   //-> CLE
		{
			CLE_transform(IN_CLE1,IN_DUM1);
			CLE_transform(IN_CLE2,IN_DUM2);
			break;
		}
		case 2:   //-> GEN
		{
			AMI_CLE_transform(IN_AMI1,IN_CLE1,IN_DUM1);
			AMI_CLE_transform(IN_AMI2,IN_CLE2,IN_DUM2);
			break;
		}
		case 3:   //-> AMI
		{
			AMI_transform(IN_AMI1,IN_DUM1);
			AMI_transform(IN_AMI2,IN_DUM2);
			break;
		}
		case 4:   //-> CLE+AMI
		{
			AMI_CLE_transform_Ori(IN_AMI1,IN_CLE1,IN_DUM1);
			AMI_CLE_transform_Ori(IN_AMI2,IN_CLE2,IN_DUM2);
			break;
		}
		default:
		{
			printf("BAD SFP_Stratety !! Should be 0,1,2,3 !!\n");
			exit(-1);
		}
	}
}
double CLEFAPS_CLE::CLEFAPS_Main_Func(vector <Align_Record> &tot)
{
	//init
	CLEFAPS_Main_Init(SFP_Strategy);
	//generate SFP
	vector <SFP_Record> SFP_Low;
	vector <SFP_Record> SFP_High;
	SFP_Low.resize(moln1*moln2);
	SFP_High.resize(moln1*moln2);
	int num1,num2;
	//calculate
	if(moln1>=moln2)
	{
		Universal_SFP_Seed_II(SFP_Low,SFP_High,SFP_L_Len,SFP_H_Len,SFP_L_Thres,SFP_H_Thres,
			(SFP_H_Len-SFP_L_Len)/2,SFP_H_Len-SFP_L_Len/2,SFP_L_Len/2,SFP_H_Len-SFP_L_Len/2,
			IN_DUM1,IN_DUM2,moln1,moln2,0,SFP_Strategy,num1,num2);   //use amino acid information
	}
	else
	{
		Universal_SFP_Seed_II(SFP_Low,SFP_High,SFP_L_Len,SFP_H_Len,SFP_L_Thres,SFP_H_Thres,
			(SFP_H_Len-SFP_L_Len)/2,SFP_H_Len-SFP_L_Len/2,SFP_L_Len/2,SFP_H_Len-SFP_L_Len/2,
			IN_DUM2,IN_DUM1,moln2,moln1,1,SFP_Strategy,num1,num2);   //use amino acid information
	}
	SFP_Low.resize(num1);
	SFP_High.resize(num2);
	//mask
	if(mas1!=0 && mas2!=0)
	{
		//-> process SFP_Low
		vector <SFP_Record> SFP_Low_;
		int num1_=Check_Mask(SFP_Low,SFP_Low_,mas1,mas2);
		SFP_Low=SFP_Low_;
		num1=num1_;
		//-> process SFP_High
		vector <SFP_Record> SFP_High_;
		int num2_=Check_Mask(SFP_High,SFP_High_,mas1,mas2);
		SFP_High=SFP_High_;
		num2=num2_;
	}
	stable_sort(SFP_Low.rbegin(),SFP_Low.rend());   //descending sort
	stable_sort(SFP_High.rbegin(),SFP_High.rend());
	//generate alignment
	double maximal;
	CLEFAPS_Second_Init();
	maximal=CLEPAPS_Full(mol1,mol2,moln1,moln2,SFP_High,SFP_Low,tot);
	CUR_MaxK=maximal;
	//second refinement
	if(REF_Strategy==1)maximal=CLEFAPS_Second_Refine(tot);
	return maximal;
}

//----- check mask -----//
int CLEFAPS_CLE::Check_Mask(vector <SFP_Record> &in, vector <SFP_Record> &out, int *mas1, int*mas2)
{
	int i,k;
	int size=(int)in.size();
	int count=0;
	out.clear();
	for(i=0;i<size;i++)
	{
		int ii=in[i].ii;
		int jj=in[i].jj;
		int len=in[i].winlen;
		//check mask
		int neo_ii=-1;
		int neo_jj=-1;
		int neo_len=0;
		int valid_len=0;
		for(k=0;k<len;k++)
		{
			if(mas1[ii+k]==mas2[jj+k])valid_len++;
			else
			{
				if(valid_len>neo_len)
				{
					neo_len=valid_len;
					neo_ii=ii+k-neo_len;
					neo_jj=jj+k-neo_len;
				}
				valid_len=0;
			}
		}
		if(valid_len>neo_len)
		{
			neo_len=valid_len;
			neo_ii=ii+k-neo_len;
			neo_jj=jj+k-neo_len;
		}
		//push_back
		if(neo_len>0)
		{
			SFP_Record neo=in[i];
			neo.ii=neo_ii;
			neo.jj=neo_jj;
			neo.winlen=neo_len;
			out.push_back(neo);
			count++;
		}
	}
	//return
	return count;
}

//----- mask CLE at the break point -----//
//-> example
/*
        RREADFAGSLAAKHHHHHHR
        11111111122222222222
        RREADFAGRRRAKHHHHHHR
*/
void CLEFAPS_CLE::Mask_CLE_BreakPoint(char *CLE, int *mask, int len)
{
	int i;
	for(i=0;i<len-3;i++)
	{
		if(mask[i]!=mask[i+1] && mask[i+1]==mask[i+2])	
		{
			CLE[i]='R';
			CLE[i+1]='R';
			CLE[i+2]='R';
		}
	}
}

