#include "BLOMAPS_Ori.h"

BLOMAPS_Ori::BLOMAPS_Ori(void)
{
	//cutoff parameter
	HSFB_LEN=12;  //IMPORTANT!!
	TMS_CUT=0.5;  //this is not important.. (anchor threshold)
	TMS_GOD=0.85; //this is not important.. (good threshold)
	ANC_CUT=0.0;  //this is not important..
	MID_CUT=1.0;  //this is not important..
	//trim parameter
	PIVOT_MOL_CUT=25;
	ANCHOR_HSFB_CUT=1;
	COLOR1_HSFB_CUT=5;
	COLOR2_HSFB_CUT=5;
	//macros
	PRITNForNOT=0;  //no printf
	UNANCHorNOT=1;  //do unanchor
	HEADCKorNOT=1;  //do headcheck
}
BLOMAPS_Ori::~BLOMAPS_Ori(void)
{
}

//------------ create & delete --------------//
void BLOMAPS_Ori::Multi_Ali_Create(int num,int len)
{
	//input
	INPUT_LEN=new int[num];
	NewArray2D(&INPUT_MOL,num,len);
	NewArray2D(&INPUT_AMI,num,len+1);
	//block
	NewArray2D(&Multi_AFB,len*num,num);
	NewArray2D(&Multi_AFB_,len*num,num);
	AFB_Len=new int[len*num];
	AFB_Len_=new int[len*num];
	//Verti+Horiz
	Mol_Anchor=new int[num];
	Pivot_Set=new int[num];
	Anchor_Set=new int[len*num];
	//temp
	TEMP_TMS=new double[num];
	NewArray2D(&TEMP_ROT,num,12);
	Pivot_Ali_Temp=new int[len];
	NewArray2D(&Pivot_Ali,num,len);
	NewArray2D(&ALIOUT,num,num*len);
	//real_temp
	Real_Temp_Ali2=new int[len];
	Real_Temp_Rot=new double[12];
	//output
	NewArray2D(&FIN_ROT,num,12);
	NewArray2D(&ALIBEST,num,num*len); //this should be LARGE!!
}
void BLOMAPS_Ori::Multi_Ali_Delete(int num,int len)
{
	//input
	delete [] INPUT_LEN;
	DeleteArray2D(&INPUT_MOL,num);
	DeleteArray2D(&INPUT_AMI,num);
	//block
	DeleteArray2D(&Multi_AFB,len*num);
	DeleteArray2D(&Multi_AFB_,len*num);
	delete [] AFB_Len;
	delete [] AFB_Len_;
	//Verti+Horiz
	delete [] Mol_Anchor;
	delete [] Pivot_Set;
	delete [] Anchor_Set;
	//temp
	delete [] TEMP_TMS;
	DeleteArray2D(&TEMP_ROT,num);
	delete [] Pivot_Ali_Temp;
	DeleteArray2D(&Pivot_Ali,num);
	DeleteArray2D(&ALIOUT,num);
	//output
	DeleteArray2D(&FIN_ROT,num);
	DeleteArray2D(&ALIBEST,num);
}

//---------------- process ------------//
void BLOMAPS_Ori::Multi_Ali_Load(int totlen,int *inlen,XYZ **mol,char **ami)
{
	int i,j;
	INPUT_NUM=totlen;
	for(i=0;i<totlen;i++)
	{
		INPUT_LEN[i]=inlen[i];
		for(j=0;j<inlen[i];j++)INPUT_MOL[i][j]=mol[i][j];
		for(j=0;j<inlen[i];j++)INPUT_AMI[i][j]=ami[i][j];
		INPUT_AMI[i][j]='\0';
	}
}
void BLOMAPS_Ori::HSFB_Generate_Full(int winlen,int &Out_Num,int **Out_AFB,int *Out_Len)
{
	int i;
	Out_Num=0;
	for(i=0;i<INPUT_NUM;i++)HSFB_Generate_Part(i,winlen,Out_Num,Out_AFB,Out_Len);
}

//================ main ================//
//-> HSFB_Partial version
double BLOMAPS_Ori::BLOMAPS_Partial(int num,int *len,XYZ **mol,char **ami,int **aliout,double **ret_rot)
{
	//[1]load structure
	int pivot,anchor;
	int pivot_tot,anchor_tot;
	int pivot_cut,anchor_cut;
	int color1_cut,color2_cut;
	Multi_Ali_Load(num,len,mol,ami);
	Other_Init_Function();
	pivot_tot=Repre_Mol_Select(num,len,mol,ami,Pivot_Set); //-> get representitive list
	pivot_cut=pivot_tot<PIVOT_MOL_CUT?pivot_tot:PIVOT_MOL_CUT;
	//--- the following is the real part!! ---//
	//init
	int match;
	double tms;
	double sumtms;
	int core_len;
	//remember, the second mol is fixed during refinement !!!
	FIN_ROT_REC.clear();
	FIN_ALI_REC.clear();
	FIN_BEST_LEN.clear();
	FIN_SCO_REC.clear();
	FIN_PIV_REC.clear();
	int solunum=0;
	for(int pp=0;pp<pivot_cut;pp++)
	{
		//init
		double wsmax=-1.0;
		int maxnum=-1;
		int core_max=-1;
		//[2]generate self HSFB
		pivot=Pivot_Set[pp]; //-> pivot set
		//printf
		if(PRITNForNOT==1)printf("HSFB_Generate[%4d]->\r",pp+1);
		Multi_AFB_Num=0;
		HSFB_Generate_Part(pivot,HSFB_LEN,Multi_AFB_Num,Multi_AFB,AFB_Len); //-> using pivot as template,generate HSFB (more)
		HSFB_Shaving(Multi_AFB_Num,Multi_AFB,AFB_Len,Multi_AFB_Num_,Multi_AFB_,AFB_Len_); //-> shaving HSFB into HSFB' (less)
		anchor_tot=HSFB_Sort_List_Generate(pivot,Multi_AFB_Num_,Multi_AFB_,AFB_Len_,Multi_AFB_Num,Multi_AFB,AFB_Len,Anchor_Set);
		//resort HSFB' by spatial consistency against HSFBs (more)
		anchor_cut=anchor_tot<ANCHOR_HSFB_CUT?anchor_tot:ANCHOR_HSFB_CUT;
		Other_Pivot_Function(pivot);
		//the above procedure only generate sort_list of HSFBs, according to its Sum-of-Pivot score
		//[3]generate initial MSA
		if(PRITNForNOT==1)printf("Initial_MSA[%4d]\r",pp+1);
		int retv;
		for(int aa=0;aa<anchor_cut;aa++)
		{
			//[3-0]get current TEMP_ROT under pivot+anchor
			anchor=Anchor_Set[aa]; //-> anchor set
			retv=HSFB_Anchor_List_Generate(pivot,Multi_AFB_[anchor],AFB_Len_[anchor],Multi_AFB_Num,Multi_AFB,AFB_Len,ANC_CUT,TEMP_ROT); //get initial MSA(more)
			if(retv==1)continue;
			color1_cut=COLOR1_HSFB_CUT;
			color2_cut=COLOR2_HSFB_CUT;
			//color1 is for anchored, while color2 is for unanchored. Also using shaved list to built initial MSA
			//[3-1]anchord refine
			int pivot_len=len[pivot];
			double tms_cutoff;
			double god_cutoff;
			double ave_tms=0.0;
			int ave_count=0;
			for(int i=0;i<num;i++)
			{
				if(PRITNForNOT==1)printf("Anchor_Process[%d][%4d,%4d]\r",pp+1,aa+1,i+1);
				if(i==pivot) //this is pivot
				{
					Mol_Anchor[i]=1; //'1' means pivot
					for(int j=0;j<pivot_len;j++)Pivot_Ali[i][j]=j;
					for(int j=0;j<12;j++)TEMP_ROT[i][j]=0.0;
					TEMP_ROT[i][0]=1.0;
					TEMP_ROT[i][4]=1.0;
					TEMP_ROT[i][8]=1.0;
					TEMP_TMS[i]=0.99+pivot_len;
					continue;
				}
				for(int j=0;j<pivot_len;j++)Pivot_Ali[i][j]=-1; //init !!
				if(Mol_Anchor[i]==1)
				{
					tms=Pairwise_Refine(i,pivot,TEMP_ROT[i],Pivot_Ali[i],TEMP_ROT[i]);
					if(tms>0.99)tms=0.99;
					TEMP_TMS[i]=tms;
					ave_tms+=tms;
					ave_count++;
				}
			}
			ave_tms/=ave_count;
			//[3-2]judge anchord and unanchored
//			tms_cutoff=TMS_CUT>ave_tms?TMS_CUT:ave_tms;
			tms_cutoff=TMS_CUT;
			for(int i=0;i<num;i++)
			{
				//the above procedure is pairwise refine
				if(TEMP_TMS[i]<tms_cutoff)
				{
					Mol_Anchor[i]=0;
					for(int j=0;j<pivot_len;j++)Pivot_Ali[i][j]=-1;
				}
			}
			//[3-3]unanchord refine
//			god_cutoff=0.5*(tms_cutoff+1.0);
			god_cutoff=TMS_GOD;
			if(UNANCHorNOT==1)
			{
				//-> given pivot and Pivot_Ali, return pivot alignment set
				double match_len=Refine_HSFB_Generate_Part(pivot,Pivot_Ali,MID_CUT,Pivot_Ali_Temp); 
				for(int i=0;i<num;i++)
				{
					if(PRITNForNOT==1)printf("UnAnchor_Process[%d][%4d,%4d]\r",pp+1,aa+1,i+1);
					//first save using refined block
					if(i==pivot)continue;
					if(Mol_Anchor[i]==1) //Deal_Anchored
					{
						if(TEMP_TMS[i]<god_cutoff) //if TM_score > 0.85, its' unneccessary to deal with !!
						{
							tms=Deal_Anchored(pivot,i,color1_cut,Multi_AFB_Num_,Multi_AFB_,AFB_Len_,Multi_AFB_Num,Multi_AFB,AFB_Len,
								Pivot_Ali_Temp,Pivot_Ali[i],TEMP_ROT[i]);
							TEMP_TMS[i]=tms;
							Mol_Anchor[i]=2; //'2' means second stage..
						}
						else
						{
							TEMP_TMS[i]+=match_len;
							Mol_Anchor[i]=2; //'2' means second stage..
						}
					}
					else                 //Deal_UnAnchored
					{
						tms=Deal_Anchored(pivot,i,color2_cut,Multi_AFB_Num_,Multi_AFB_,AFB_Len_,Multi_AFB_Num,Multi_AFB,AFB_Len,
							Pivot_Ali_Temp,Pivot_Ali[i],TEMP_ROT[i]);
						TEMP_TMS[i]=tms;
						match=(int)tms;
						tms=tms-match;
						if(tms>tms_cutoff)Mol_Anchor[i]=2; //'2' means second stage success, also means anchored	
					}
					//second save using total pairwise alignment
					if(Mol_Anchor[i]==0)
					{
						tms=Pairwise_Total(i,pivot,Pivot_Ali_Temp,Real_Temp_Ali2,Real_Temp_Rot);
						if(tms>TEMP_TMS[i])
						{
							EqualArray(Pivot_Ali[i],Real_Temp_Ali2,pivot_len);
							EqualArray(TEMP_ROT[i],Real_Temp_Rot,12);
							TEMP_TMS[i]=tms;
							Mol_Anchor[i]=3; //'3' means third stage 'success', also means unanchored
						}
					}
				}
			}
			//[4]final refine
			if(PRITNForNOT==1)printf("Final_Refine[%4d][%4d]\r",pp+1,aa+1);
			//-> the input Pivot_Ali is all-to-one style, the output ALIOUT is MSA style
			sumtms=Terminal_Refine_Main(pivot,TEMP_ROT,mol,len,num,Pivot_Ali,len[pivot],ALIOUT,TOTLEN,TEMP_ROT,core_len);
			//[5]record best
			if(sumtms>wsmax)
			{
				wsmax=sumtms;
				maxnum=pivot;
				core_max=core_len;
				EqualArray2D(FIN_ROT,TEMP_ROT,INPUT_NUM,12);
				EqualArray2D(ALIBEST,ALIOUT,INPUT_NUM,TOTLEN);
				BESTLEN=TOTLEN;
			}
			//printf
//			printf("->cur[%3d][%3d]->obj[%6.2f]\r",pp+1,aa+1,wsmax);
		}
		//termi
		{
//			FIN_ROT_REC.push_back(FIN_ROT);
//			FIN_ALI_REC.push_back(ALIBEST);
			vector <vector <double > > fin_rot_tmp (INPUT_NUM, vector <double >(12));
			vector <vector <int > > fin_ali_tmp (INPUT_NUM, vector <int >(BESTLEN));
			for(int ii=0;ii<INPUT_NUM;ii++)for(int jj=0;jj<12;jj++)fin_rot_tmp[ii][jj]=FIN_ROT[ii][jj];
			for(int ii=0;ii<INPUT_NUM;ii++)for(int jj=0;jj<BESTLEN;jj++)fin_ali_tmp[ii][jj]=ALIBEST[ii][jj];
			//check overlap
			int correct=1;
//			for(int kk=0;kk<(int)FIN_ALI_REC.size();kk++)
//			{
//				double overlap_sco=Compare_Two_MSA(fin_ali_tmp,FIN_ALI_REC[kk],BESTLEN,FIN_BEST_LEN[kk],INPUT_NUM);
//				if(overlap_sco > MSA_COMP_THRES)
//				{
//					correct=0;
//					break;
//				}
//			}
			//final add
			if(correct==1)
			{
				FIN_ROT_REC.push_back(fin_rot_tmp);
				FIN_ALI_REC.push_back(fin_ali_tmp);
				FIN_BEST_LEN.push_back(BESTLEN);
				FIN_SCO_REC.push_back(wsmax);
				FIN_PIV_REC.push_back(maxnum);
				FIN_CORE_LEN.push_back(core_max);
				solunum++;
			}
		}
	}

	//[output]
	if(solunum>0)
	{
		//get max
		double curmax=FIN_SCO_REC[0];
		int curnum=0;
		for(int i=1;i<solunum;i++)
		{
			if(FIN_SCO_REC[i]>curmax)
			{
				curmax=FIN_SCO_REC[i];
				curnum=i;
			}
		}
		//out max
		FIN_PIVOT=curnum;
		for(int ii=0;ii<INPUT_NUM;ii++)for(int jj=0;jj<12;jj++)ret_rot[ii][jj]=FIN_ROT_REC[curnum][ii][jj];
		for(int ii=0;ii<INPUT_NUM;ii++)for(int jj=0;jj<FIN_BEST_LEN[curnum];jj++)aliout[ii][jj]=FIN_ALI_REC[curnum][ii][jj];
//		EqualArray2D(ret_rot,FIN_ROT_REC[curnum],INPUT_NUM,12);
//		EqualArray2D(aliout,FIN_ALI_REC[curnum],INPUT_NUM,FIN_BEST_LEN[curnum]);
		return curmax;
	}
	else return -1;
}
