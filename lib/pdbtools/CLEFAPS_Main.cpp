#include "CLEFAPS_Main.h"


//------ constructor -------//
CLEFAPS_Main::CLEFAPS_Main(int num,int CLESUM)
:TM_align(num),Bioinfo_Code(CLESUM),
Ali_Ori(num),Ali_AFP(num),
CLEFAPS_CLE(num),TM_Align_Main(num)
{
	FM_Kill_Gap=0;   //default: don't Kill Gaps !!
	Kill_Frag=0;     //default: don't Kill Frag !!
	Kill_Frag_Para_0=25;  //-> min_small_len, 25
	Kill_Frag_Para_1=25;  //-> gap_len, 25
	Kill_Frag_Para_2=15;  //-> nongap_len, 15
}
CLEFAPS_Main::~CLEFAPS_Main(void)
{
}


//================= Alignment_Gap_Process ============//__130830__//
//[note]: gap_len=25, nongap_len=15
int CLEFAPS_Main::Alignment_Gap_Process(int *ali2_in,int moln1,int moln2,int *ali2_out,
	int gap_len,int nongap_len)
{
	int i,k;
	int *ali_rec=new int[moln2];
	for(i=0;i<moln2;i++)ali_rec[i]=0; //normal
	for(i=0;i<moln2;i++)if(ali2_in[i]!=-1)ali_rec[i]=1; //occupied

	//[0-1] get z_score for gap_length
//	double gap_mean,gap_vari;
	vector <int> gap_length;
	int first;
	int prev,curr;
	int start1,start2;
	int end1,end2;
	int start=-1;
	int len=0;

	//[0-2] get z_score for afp_length
//	double afp_mean,afp_vari;
	vector <int> afp_length;
	first=1;
	start=0;
	len=0;
	for(i=0;i<moln2;i++)
	{
		//check ali1 gap
		if(i>0 && ali2_in[i]!=-1 && ali2_in[i-1]!=-1)
		{
			prev=ali2_in[i-1];
			curr=ali2_in[i];
			if(curr-prev!=1)
			{
				if(first==0)
				{
					first=1;
					if(len<4)
					{
						for(k=0;k<len;k++)ali2_in[start+k]=-1; //saved !!
					}
					else 
					{
						afp_length.push_back(len);
					}
				}
				if(first==1)
				{
					first=0;
					start=i;
					len=0;
				}
				len++;
				continue;
			}
		}
		//check ali2 gap
		if(ali2_in[i]!=-1)
		{
			if(first==1)
			{
				first=0;
				start=i;
				len=0;
			}
			len++;
		}
		else
		{
			if(first==0)
			{
				first=1;
				if(len<4)
				{
					for(k=0;k<len;k++)ali2_in[start+k]=-1; //saved !!
				}
				else
				{
					afp_length.push_back(len);
				}
			}
		}
	}
	//tail process
	if(first==0)
	{
		first=1;
		if(len<4)
		{
			for(k=0;k<len;k++)ali2_in[start+k]=-1; //saved !!
		}
		else
		{
			afp_length.push_back(len);
		}
	}
	//[0-3] calculate z_score
//	int size;
	int gap_len_real,nongap_len_real;
	gap_len_real=gap_len;
	nongap_len_real=nongap_len;
	
	//[1] fix nongap
	start=-1;
	len=0;
	for(i=0;i<moln2;i++)
	{
		if(ali2_in[i]!=-1) //nongap
		{
			//check previous
			if(i>0 && ali2_in[i-1]!=-1)
			{
				prev=ali2_in[i-1];
				curr=ali2_in[i];
				if(prev!=curr-1) //gap!!!
				{
					start=i;
					len=0;
					continue;
				}
			}
			//process
			if(len>0)
			{
				if(start==-1)
				{
					start1=0;
					start2=0;
				}
				else
				{
					start1=ali2_in[start];
					start2=start;
				}
				end1=ali2_in[i];
				end2=i;
				//check
				int maxnum=(end1-start1)>(end2-start2)?(end1-start1):(end2-start2);
				if(maxnum<gap_len_real)
				{
					for(k=start+1;k<start+1+len;k++)ali_rec[k]=1; //occupied
				}
			}
			//clear
			start=i;
			len=0;
		}
		else  //gap
		{
			len++;
		}
	}
	//tail process
	if(len>0)
	{
		if(start==-1)
		{
			start1=0;
			start2=0;
		}
		else
		{
			start1=ali2_in[start];
			start2=start;
		}
		end1=moln1;
		end2=moln2;
		//check
		int maxnum=(end1-start1)>(end2-start2)?(end1-start1):(end2-start2);
		if(maxnum<gap_len_real)
		{
			for(k=start+1;k<start+1+len;k++)ali_rec[k]=1; //occupied
		}
	}

	//[2] hard fix
	first=1;
	start=0;
	len=0;
	for(i=0;i<moln2;i++)
	{
		//check ali1 gap
		if(i>0 && ali2_in[i]!=-1 && ali2_in[i-1]!=-1)
		{
			prev=ali2_in[i-1];
			curr=ali2_in[i];
			if(curr-prev>=gap_len_real)
			{
				if(first==0)
				{
					first=1;
					if(len>=nongap_len_real)for(k=0;k<len;k++)ali_rec[start+k]=2; //saved !!
				}
				if(first==1)
				{
					first=0;
					start=i;
					len=0;
				}
				len++;
				continue;
			}
		}
		//check ali2 gap
		if(ali_rec[i]==1)
		{
			if(first==1)
			{
				first=0;
				start=i;
				len=0;
			}
			len++;
		}
		else
		{
			if(first==0)
			{
				first=1;
				if(len>=nongap_len_real)for(k=0;k<len;k++)ali_rec[start+k]=2; //saved !!
			}
		}
	}
	//tail process
	if(first==0)
	{
		first=1;
		if(len>=nongap_len_real)for(k=0;k<len;k++)ali_rec[start+k]=2; //saved !!
	}

	//[3] generate new alignment
	int lali=0;
	for(i=0;i<moln2;i++)
	{
		if(ali2_in[i]==-1)ali2_out[i]=-1;
		else
		{
			if(ali_rec[i]==2)
			{
				ali2_out[i]=ali2_in[i];
				lali++;
			}
			else ali2_out[i]=-1;
		}
	}
	
	//delete & return
	delete [] ali_rec;
	return lali;
}


//---------- post_refine ----------//__110830__//
//to kill gaps!!
double CLEFAPS_Main::IN_Calc_Frag_Score(int ii,int jj,int *ali1,int *ali2,
					int moln1,int moln2,XYZ *mol1,XYZ *mol2,double *rotmat_)
{
	//-2,+1
	int left=CLEP_MIN_LEN/2;
	int right=CLEP_MIN_LEN/2;
	if(CLEP_MIN_LEN%2==0)right--;
	//judge
	if(ii-left<0 || jj-left<0)return -1.0;
	if(ii+right>=moln1 || jj+right>=moln2)return -1.0;
	//check
	int k;
	int pos1,pos2;
	for(k=0;k<CLEP_MIN_LEN;k++)
	{
		pos1=ii-left+k;
		pos2=jj-left+k;
		if(ali1[pos1]>=0 || ali2[pos2]>=0)return -1.0;
	}
	//final score
	int i;
	double vec1,vec2,vec3;
	double vect;
	int tot_num;
	int range=1;
	double dist2;
	double ori_d=d0*d0;
	double tms;
	XYZ xyz;
	double ws_sco=0.0;
	for(i=0;i<CLEP_MIN_LEN;i++)
	{
		pos1=ii-left+i;
		pos2=jj-left+i;
		//rot_point
		if(rotmat_!=0)Ali_Cache_Point(mol1,pos1,xyz,rotmat_);
		else xyz=mol1[pos1];
		dist2=mol2[pos2].distance_square(xyz);
		//nearby
		vec1=0.0;
		vec2=0.0;
		for(k=1;k<=range;k++)
		{
			vec1+=TM_Vector_Score_Forward(pos1,pos2,k,rotmat_,mol1,mol2,moln1,moln2);
			vec2+=TM_Vector_Score_Backward(pos1,pos2,k,rotmat_,mol1,mol2,moln1,moln2);
		}
		tot_num=2;
		if(TM_cb1!=0 && TM_cb2!=0)
		{
			vec3=TM_Vector_Score_CB(pos1,pos2,rotmat_,mol1,mol2,moln1,moln2);
			tot_num=3;
		}
		else vec3=0.0;
		vect=1.0*(vec1+vec2)/range;
		vect=1.0*(vect+vec3)/tot_num;
		tms=1.0/(1.0+dist2/ori_d);
		//total score
		if(TMali_Weight==0) //no weight
		{
			if(TM_Vect_Score==0)ws_sco+=tms;  // original TM_score
			else ws_sco+=vect*tms;            // current Vect_score
		}
		else                //has weight
		{
			if(TM_Vect_Score==0)ws_sco+=TMali_Weight[pos1*moln2+pos2]*tms;  // original TM_score
			else ws_sco+=TMali_Weight[pos1*moln2+pos2]*vect*tms;            // current Vect_score
		}
	}
	return ws_sco;
}
int CLEFAPS_Main::IN_Ali_To_Cor(int *AFP_Cor,int thres,
						int moln1,int moln2,int *ali1,int *ali2,
						vector<pair<int,int> > &ws_pair_record) 
{
	int i,k;
	int num;
	int ii,jj;
	int count;
	int isFirst;
	int isLast;
	int type;
	int head1,head2;
	int index;

	//init
	num=0;
	AFP_Cor[0]=0;
	head1=-1;
	head2=-1;
	isLast=0;
	isFirst=1;
	type=0;
	count=INT_MIN_NUM;
	ii=-1;
	jj=-1;
	ws_pair_record.clear();
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]==-1) //purely blank
		{
			if(isFirst==0)
			{
				if(count>=thres) 
				{
					AFP_Cor[0]++;
					index=AFP_Cor[0];
					AFP_Cor[4*index+0] = type;
					AFP_Cor[4*index+1] = head1;
					AFP_Cor[4*index+2] = head2;
					AFP_Cor[4*index+3] = count;
					num+=count;
				}
				else
				{
					for(k=0;k<count;k++)
					{
						ali1[head1+k]=-1;
						ali2[head2+k]=-1;
						ws_pair_record.push_back(pair<int,int>(head1+k,head2+k)); 
					}
				}
				count=0;
				isFirst=1;
			}
			continue;
		}

		if(isFirst==1)
		{
ws_init:
			isFirst=0;
			ii=ali2[i];
			type=1; // >0 mode
			jj=i;
			count=1;
			head1=ii;
			head2=jj;
#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif
			continue;
		}
		if(i==jj+1&&ali2[i]==ii+1)
		{
			ii=ali2[i];
			jj=i;
			count++;
			continue;
		}

ws_end:
		if(count>=thres) 
		{
			AFP_Cor[0]++;
			index=AFP_Cor[0];
			AFP_Cor[4*index+0] = type;
			AFP_Cor[4*index+1] = head1;
			AFP_Cor[4*index+2] = head2;
			AFP_Cor[4*index+3] = count;
			num+=count;
		}
		else
		{
			for(k=0;k<count;k++)
			{
				ali1[head1+k]=-1;
				ali2[head2+k]=-1;
				ws_pair_record.push_back(pair<int,int>(head1+k,head2+k)); 
			}
		}

		if(isLast==1)goto end;
		else goto ws_init;
	}

	if(count==INT_MIN_NUM)goto end;
	isLast=1;
	goto ws_end;
end:
	return num;
}
void CLEFAPS_Main::IN_Elongation(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
						 XYZ *mol1,XYZ *mol2,int moln1,int moln2)
{
	int i;
	int ii,jj;
	int index1,index2;
	double Dou_Thres=thres*thres;
	int success;
	XYZ temp1,temp2,temp;
	int cut1,cut2,cut;

	//---- ws_record ---//
	int *wsali1=new int[moln1];
	int *wsali2=new int[moln2];
	memset(wsali1,0,moln1*sizeof(int));
	memset(wsali2,0,moln2*sizeof(int));
	for(i=0;i<moln2;i++)
	{
		ii=ali2[i];
		jj=i;
		if(ii<0)continue;
		wsali1[ii]=1;
		wsali2[jj]=1;
	}

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
			if(jj>=0 && ii>=0) break;
			else if(jj>=0 && ii<0)
			{
				//check
				if(wsali2[jj]==1)break;
				//calc
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
			else if(jj<0 && ii>=0)
			{
				//check
				if(wsali1[ii]==1)break;
				//calc
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
			if(jj>=0 && ii>=0) break;
			else if(jj>=0 && ii<0)
			{
				//check
				if(wsali2[jj]==1)break;
				//calc
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
			else if(jj<0 && ii>=0)
			{
				//check
				if(wsali1[ii]==1)break;
				//calc
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

	//delete
	delete [] wsali1;
	delete [] wsali2;
}
void CLEFAPS_Main::FM_Align_Refine(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_in,int *ali2_out,double *rotmat_)
{
	int i;
	int pos;
	int *wsali1=new int[moln1];
	int *wsali2=new int[moln2];
	//init
	memset(wsali1,-1,moln1*sizeof(int));
	memset(wsali2,-1,moln2*sizeof(int));
	for(i=0;i<moln2;i++)
	{
		pos=ali2_in[i];
		if(pos!=-1)
		{
			wsali1[pos]=i;
			wsali2[i]=pos;
		}
	}
	//kill gaps
	vector<pair<int,int> > ws_pair_record;
	IN_Ali_To_Cor(AFP_Cor,CLEP_MIN_LEN,moln1,moln2,wsali1,wsali2,ws_pair_record);
	//refine iteration
	IN_Elongation(rotmat_,FIN_CUT,ANG_CUT,wsali1,wsali2,mol1,mol2,moln1,moln2);
	//final add
	int size=(int)ws_pair_record.size();
	Fast_Sort<double> ws_fast_sort_double;
	double *ws_rec=new double[size];
	int *ws_idx=new int[size];
	int ii,jj;
	for(i=0;i<size;i++)
	{
		ii=ws_pair_record[i].first;
		jj=ws_pair_record[i].second;
		ws_rec[i]=IN_Calc_Frag_Score(ii,jj,wsali1,wsali2,moln1,moln2,mol1,mol2,rotmat_);
	}
	ws_fast_sort_double.fast_sort_1(ws_rec,ws_idx,size);
	//final check
	int k;
	int pos1,pos2;
	int index;
	int ws_correct;
	int left=CLEP_MIN_LEN/2;
	int right=CLEP_MIN_LEN/2;
	if(CLEP_MIN_LEN%2==0)right--;
	for(i=0;i<size;i++)
	{
		index=ws_idx[i];
		if(ws_rec[index]<0)break;
		//check
		ws_correct=1;
		ii=ws_pair_record[index].first;
		jj=ws_pair_record[index].second;
		for(k=0;k<CLEP_MIN_LEN;k++)
		{
			pos1=ii-left+k;
			pos2=jj-left+k;
			if(wsali1[pos1]>=0 || wsali2[pos2]>=0)
			{
				ws_correct=0;
				break;
			}
		}
		if(ws_correct==0)continue;
		//add
		for(k=0;k<CLEP_MIN_LEN;k++)
		{
			pos1=ii-left+k;
			pos2=jj-left+k;
			wsali1[pos1]=pos2;
			wsali2[pos2]=pos1;
		}
	}
	//refine iteration
	IN_Elongation(rotmat_,FIN_CUT,ANG_CUT,wsali1,wsali2,mol1,mol2,moln1,moln2);
	//evaluate
	for(i=0;i<moln2;i++)ali2_out[i]=-1;
	for(i=0;i<moln2;i++)
	{
		pos=wsali2[i];
		if(pos!=-1)ali2_out[i]=pos;
	}
	//delete
	delete [] wsali1;
	delete [] wsali2;
	delete [] ws_rec;
	delete [] ws_idx;
}

//==================== Gap_Kill ================//
void CLEFAPS_Main::Retreive_Alignment(char *out1,char *out2,int *ali2_,int moln2)
{
	//final retrieve
	int i;
	int ii,jj;
	for(i=0;i<moln2;i++)ali2_[i]=-1;
	ii=0;
	jj=0;
	for(i=0;i<(int)strlen(out1);i++)
	{
		if( (out1[i]!='-') && (out2[i]!='-') )
		{
			ali2_[jj]=ii;
			ii++;
			jj++;
		}
		else
		{
			if(out1[i]!='-')ii++;
			if(out2[i]!='-')jj++;
		}
	}
}
void CLEFAPS_Main::Transform_Ali(int *ali1,int *ali2,int moln1,int moln2,char *ami1,char *ami2,char *out1,char *out2)
{
	//start
	int i,j;
	int ii,jj;
	int wlen;
	int pre_ii=0;
	int pre_jj=0;
	int count=0;
	for(i=1;i<=moln2;i++)
	{
		jj=i;
		ii=ali2[i-1];  //ali1 starts from 0, correspondence also from 0
		if(ii==-1)
		{
			continue;
		}
		else
		{
			ii++;
			//previous_path
			wlen=ii-pre_ii;
			for(j=1;j<wlen;j++)
			{
				pre_ii++;
				DP_align1[count]=pre_ii;
				DP_align2[count]=-pre_jj;
				count++;
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
				pre_jj++;
				DP_align1[count]=-pre_ii;
				DP_align2[count]=pre_jj;
				count++;
			}
			//current_path
			DP_align1[count]=ii;
			DP_align2[count]=jj;
			count++;
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	//termi
	pre_ii++;
	for(i=pre_ii;i<=moln1;i++)
	{
		DP_align1[count]=i;
		DP_align2[count]=-pre_jj;
		count++;
	}
	pre_jj++;
	for(i=pre_jj;i<=moln2;i++)
	{
		DP_align1[count]=-moln1;
		DP_align2[count]=i;
		count++;
	}
	//final
	vector<pair<int,int> > alignment;
	alignment.resize(count);
	for(i=0;i<count;i++)
	{
		alignment[i].first=DP_align1[i];
		alignment[i].second=DP_align2[i];
	}
	int size=(int)alignment.size();
	//output_AMI
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		if(ii<=0)out1[i]='-';
		else out1[i]=ami1[ii-1];
	}
	out1[i]='\0';
	for(i=0;i<size;i++)
	{
		jj=alignment[i].second;
		if(jj<=0)out2[i]='-';
		else out2[i]=ami2[jj-1];
	}
	out2[i]='\0';
}
void CLEFAPS_Main::Insert_Refine_Single(char *seq1,char *seq2,int start1,int start2,int moln1,int moln2,
							 double *ws_in,char *out1,char *out2,int HEADorTAIL)
{
	int i,j;
	int len1,len2;
	len1=(int)strlen(seq1);
	len2=(int)strlen(seq2);
	//make score
	double *CUR_DP_sco=new double[len1*len2];
	for(i=0;i<len1;i++)
	{
		int cur_index=i*len2;
		int nxt_index=(start1+i)*moln2+start2;
		for(j=0;j<len2;j++)
		{
			CUR_DP_sco[cur_index+j]=ws_in[nxt_index+j];
		}
	}
	//dyna_prog
	vector<pair<int,int> > alignment;
	double ali_sco;
	int ali_len;
	if(HEADorTAIL==-1) //head
	{
		ali_len=Advance_Align_Dyna_Prog_Double(len1,len2,CUR_DP_sco,-11,-1,-11,-1,0,-110,0,-110,
			alignment,ali_sco);
	}
	else if(HEADorTAIL==1) //tail
	{
		ali_len=Advance_Align_Dyna_Prog_Double(len1,len2,CUR_DP_sco,-11,-1,-11,-1,-110,0,-110,0,
			alignment,ali_sco);
	}
	else //midd
	{
		ali_len=Advance_Align_Dyna_Prog_Double(len1,len2,CUR_DP_sco,-11,-1,-11,-1,-110,-110,-110,-110,
			alignment,ali_sco);
	}
	delete [] CUR_DP_sco;
	//retrive
	int ii,jj;
	int totnum=(int)alignment.size();
	for(i=0;i<totnum;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0&&jj>0)
		{
			out1[i]=seq1[ii-1];
			out2[i]=seq2[jj-1];
		}
		else
		{
			if(ii>0)
			{
				out1[i]=seq1[ii-1];
				out2[i]='-';
			}
			if(jj>0)
			{
				out1[i]='-';
				out2[i]=seq2[jj-1];
			}
		}
	}
	out1[i]='\0';
	out2[i]='\0';
}
void CLEFAPS_Main::Insert_Refine_Total(char *ali1,char *ali2,double *ws_in,int moln1,int moln2,
										  char *out1,char *out2,int CutK)
{
	char *seq1=new char[moln1+moln2+1];
	char *seq2=new char[moln1+moln2+1];
	char *wwout1=new char[moln1+moln2+1];
	char *wwout2=new char[moln1+moln2+1];
	int len1=(int)strlen(ali1);
	int len2=(int)strlen(ali2);
	if(len1!=len2)
	{
		delete [] seq1;
		delete [] seq2;
		delete [] wwout1;
		delete [] wwout2;
		fprintf(stderr,"len1 != len2, [%d!=%d]\n",len1,len2);
		exit(-1);
	}
	int len=len1;
	//extract sub_alignment
	int i,k;
	int start1;
	int start2;
	int cur1=0;
	int cur2=0;
//	int first=1;
//	int mislen=0;
	int ii,jj;
	int HEADorTAIL;
	string wsout1="";
	string wsout2="";
	int wsstart,wslen;
	int wsprev=0;
	int ws1,ws2;
	string wspre1,wspre2;
	wspre1=ali1;
	wspre2=ali2;
	string wstemp1,wstemp2;
	//Ix->Iy kill
	for(i=0;i<len-1;i++)
	{
		//--- find Ix->Iy ----//
		if((ali1[i]!='-'&&ali2[i]=='-') && (ali1[i+1]=='-'&&ali2[i+1]!='-'))
		{
			//init
			start1=cur1;
			start2=cur2;
			//forward
			wsstart=i;
			wslen=1;
			ws1=1;
			for(;;)
			{
				if(wsstart-1<0)break;
				if(!(ali1[wsstart-1]!='-'&&ali2[wsstart-1]=='-'))break;
				wsstart--;
				wslen++;
				start1--;
				ws1++;
			}
			//backward
			ws2=0;
			for(;;)
			{
				if(wsstart+wslen>=len)break;
				if(!(ali1[wsstart+wslen]=='-'&&ali2[wsstart+wslen]!='-'))break;
				wslen++;
				i++;
				cur2++;
				ws2++;
			}
			//copy
			strncpy(wwout1,ali1+wsstart,wslen);
			strncpy(wwout2,ali2+wsstart,wslen);
			wwout1[wslen]='\0';
			wwout2[wslen]='\0';
			//judge
			if(ws1<CutK || ws2<CutK)
			{
				//retrieve
				ii=0;
				jj=0;
				for(k=0;k<wslen;k++)
				{
					if(ali1[wsstart+k]!='-')
					{
						seq1[ii]=ali1[wsstart+k];
						ii++;
					}
					if(ali2[wsstart+k]!='-')
					{
						seq2[jj]=ali2[wsstart+k];
						jj++;
					}
				}
				seq1[ii]='\0';
				seq2[jj]='\0';
				if(start1==0&&start2==0)HEADorTAIL=-1;     //head
				else if(wsstart+wslen==len)HEADorTAIL=1; //tail
				else HEADorTAIL=0;
				if(HEADorTAIL!=0)
				{
					Insert_Refine_Single(seq1,seq2,start1,start2,moln1,moln2,
						ws_in,wwout1,wwout2,HEADorTAIL);
				}
			}
			//combine
			wstemp1=wspre1.substr(wsprev,wsstart-wsprev);
			wstemp2=wspre2.substr(wsprev,wsstart-wsprev);
			wsout1+=wstemp1;
			wsout1+=wwout1;
			wsout2+=wstemp2;
			wsout2+=wwout2;
			//final
			wsprev=wsstart+wslen;
		}
		else //original
		{
			if(ali1[i]!='-')cur1++;
			if(ali2[i]!='-')cur2++;
		}
	}
	//final
	wstemp1=wspre1.substr(wsprev,len-wsprev);
	wstemp2=wspre2.substr(wsprev,len-wsprev);
	wsout1+=wstemp1;
	wsout2+=wstemp2;
	strcpy(out1,wsout1.c_str());
	strcpy(out2,wsout2.c_str());
	//delete 
	delete [] seq1;
	delete [] seq2;
	delete [] wwout1;
	delete [] wwout2;
}
void CLEFAPS_Main::Double_Gap_Refine(char *ali1,char *ali2,double *ws_in,int moln1,int moln2,
										char *out1,char *out2,int CutK,double thres)
{
	char *wwout1=new char[moln1+moln2+1];
	char *wwout2=new char[moln1+moln2+1];
	int len1=(int)strlen(ali1);
	int len2=(int)strlen(ali2);
	if(len1!=len2)
	{
		delete [] wwout1;
		delete [] wwout2;
		fprintf(stderr,"len1 != len2, [%d!=%d]\n",len1,len2);
		exit(-1);
	}
	int len=len1;
	//extract sub_alignment
	int i,k;
	int start1;
	int start2;
	int cur1=0;
	int cur2=0;
//	int first=1;
//	int mislen=0;
	string wsout1="";
	string wsout2="";
	int wsstart,wslen;
	int wsprev=0;
	string wspre1,wspre2;
	wspre1=ali1;
	wspre2=ali2;
	string wstemp1,wstemp2;
	int isfound;
	double ws_prev,ws_curr;
	ws_prev=0.0;
	//Ix->Iy kill
	for(i=0;i<len-1;i++)
	{
		//--- find double gap ---//type_1
		if((ali1[i]=='-'|| ali2[i]=='-') && (ali1[i+1]!='-'&&ali2[i+1]!='-'))
		{
			//finding
			start1=cur1;
			start2=cur2;
			wsstart=i;
			wslen=2;
			isfound=0;
			for(;;)
			{
				if(wsstart+wslen>=len)
				{
					isfound=-1;
					break;
				}
				if((ali1[wsstart+wslen]=='-') || (ali2[wsstart+wslen]=='-'))break;
				wslen++;
			}
			if(isfound==-1)break;
			if((ali1[i]=='-'&& ali2[i]!='-') && (ali1[wsstart+wslen]!='-') && (ali2[wsstart+wslen]=='-'))isfound=1; //type_1
			if((ali1[i]!='-'&& ali2[i]=='-') && (ali1[wsstart+wslen]=='-') && (ali2[wsstart+wslen]!='-'))isfound=2; //type_2
			//judge
			if(isfound==0)continue;
			if(wslen>CutK)continue;
			//calc
			if(isfound==1) //type_1
			{
				ws_prev=0.0;
				for(k=1;k<wslen-1;k++)ws_prev+=ws_in[(cur1+k-1)*moln1+cur2+k];
			}
			if(isfound==2) //type_2
			{
				ws_prev=0.0;
				for(k=1;k<wslen-1;k++)ws_prev+=ws_in[(cur1+k)*moln1+cur2+k-1];
			}
			ws_curr=0.0;
			for(k=0;k<wslen-1;k++)ws_curr+=ws_in[(cur1+k)*moln1+cur2+k];
			if(ws_curr<ws_prev*thres)continue;
			//change
			for(k=0;k<wslen-1;k++)
			{
				wwout1[k]=ali1[cur1+k];
				wwout2[k]=ali2[cur2+k];
			}
			wwout1[k]='\0';
			wwout2[k]='\0';
			//record
			wstemp1=wspre1.substr(wsprev,wsstart-wsprev);
			wstemp2=wspre2.substr(wsprev,wsstart-wsprev);
			wsout1+=wstemp1;
			wsout1+=wwout1;
			wsout2+=wstemp2;
			wsout2+=wwout2;
			//final
			i+=wslen;
			wsprev=wsstart+wslen;
		}
	}
	//final
	wstemp1=wspre1.substr(wsprev,len-wsprev);
	wstemp2=wspre2.substr(wsprev,len-wsprev);
	wsout1+=wstemp1;
	wsout2+=wstemp2;
	strcpy(out1,wsout1.c_str());
	strcpy(out2,wsout2.c_str());
	//delete
	delete [] wwout1;
	delete [] wwout2;
}
//-------- Kill_Gap_Total ----------//
double CLEFAPS_Main::CLEF_Make_Score(Align_Record & align_record,
								  XYZ *mol1,XYZ *mol2,int moln1,int moln2,
								  double *rotmat,int *ali2)
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
	for(i=0;i<12;i++)align_record.rotmat[i]=rotmat[i];
	//calculate
	double d_0=d0;
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
		Ali_Cache_Point(mol1,ii,temp,rotmat);
		dist2=mol2[jj].distance_square(temp);
		sco+=1.0/(1.0+dist2/ori_d);
		rms+=dist2;
		count++;
	}
	//evaluate
	rmsd=sqrt(1.0*rms/count);
	lali=count;
	int TM_smaller=moln1<moln2?moln1:moln2;
	tms=1.0*sco/TM_smaller;
	//final
	align_record.RMSD=rmsd;
	align_record.lali=lali;
	align_record.main_sco=tms;
	return tms;
}
void CLEFAPS_Main::CLEF_Single_Match_Kill(XYZ *mol1,XYZ *mol2,int moln1,int moln2,vector <Align_Record> &tot)
{
	int i,k;
	int size=(int)tot.size();
	char *seq1=new char[moln1+moln2+1];
	char *seq2=new char[moln1+moln2+1];
	char *out1=new char[moln1+moln2+1];
	char *out2=new char[moln1+moln2+1];
	int *ali1_=new int[moln1+moln2];
	int *ali2_=new int[moln1+moln2];
	double rotmat[12];
	double tms;
	double wms;
	double rmsd;
	int lali;
	int CutK1=4;
//	int CutK2=8;
	double CutK2_thres=0.9;
	for(i=0;i<size;i++)
	{
		//get data
		for(k=0;k<moln2;k++)ali2_[k]=tot[i].alignment[k];
		for(k=0;k<12;k++)rotmat[k]=tot[i].rotmat[k];
		if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
		if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);
		//kill single match
//		int Vect_Sco_Ori=TM_Vect_Score;
//		TM_Vect_Score=1;            //switch on
		FM_Align_Refine(mol1,mol2,moln1,moln2,ali2_,ali2_,rotmat);
		//kill single gap
		TM_Align_Get_Score(mol1,mol2,moln1,moln2,ali2_); //only consider vector score !!!
		Transform_Ali(ali1_,ali2_,moln1,moln2,IN_AMI1,IN_AMI2,seq1,seq2);
		Insert_Refine_Total(seq1,seq2,TM_DP_sco,moln1,moln2,out1,out2,CutK1);
		//kill double gap
		strcpy(seq1,out1);
		strcpy(seq2,out2);
		Double_Gap_Refine(seq1,seq2,TM_DP_sco,moln1,moln2,out1,out2,CutK1,CutK2_thres);
//		TM_Vect_Score=Vect_Sco_Ori; //switch off
		//final
		Retreive_Alignment(out1,out2,ali2_,moln2);
		tms=TM_Align_TM_Score_Simp(mol1,mol2,moln1,moln2,ali2_,rmsd,lali);
		//kill fragment
		if(Kill_Frag==1)
		{
			vector <int> ali2_bak;
			ali2_bak.resize(moln2);
			for(int wsi=0;wsi<moln2;wsi++)ali2_bak[wsi]=ali2_[wsi];
			int temp_small=moln1<moln2?moln1:moln2;
			int ws_ret=0;
			if(temp_small>Kill_Frag_Para_0)
			{
				// parameter is 25, 25 and 15
				// Kill_Frag_Para_0 -> min_small_len
				// Kill_Frag_Para_1 -> gap_len
				// Kill_Frag_Para_2 -> nongap_len
				ws_ret=Alignment_Gap_Process(ali2_,moln1,moln2,ali2_,Kill_Frag_Para_1,Kill_Frag_Para_2);
			}
			if(ws_ret==0)
			{
				for(int wsi=0;wsi<moln2;wsi++)ali2_[wsi]=ali2_bak[wsi];
			}
		}
		//final assessment
		int TM_DIST_CUT_=TM_DIST_CUT;
		TM_DIST_CUT=0;
		tms=TM_Align_TM_Score_Simp(mol1,mol2,moln1,moln2,ali2_,rmsd,lali);
		TM_DIST_CUT=TM_DIST_CUT_;
		if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
		if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);
		CLEF_Make_Score(tot[i],mol1,mol2,moln1,moln2,finmat,ali2_);
		wms=TM_Align_Get_Score_Simp(mol1,mol2,finmat,moln1,moln2,ali2_);
		tot[i].main_sco=tms*wms;   //-> we consider tms*wms

//printf("tms=%lf , wms=%lf, tms*wms=%lf \n",tms,wms,tms*wms);

	}

	//delete
	delete [] seq1;
	delete [] seq2;
	delete [] out1;
	delete [] out2;
	delete [] ali1_;
	delete [] ali2_;
}

//============================= the following is the main function for TM_align =======================//
void CLEFAPS_Main::FM_Get_Initial_AMI(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	//assign SSE
	int *ws_ami1=new int[moln1+moln2];
	int *ws_ami2=new int[moln1+moln2];
	AMI_transform(IN_AMI1,ws_ami1);
	AMI_transform(IN_AMI2,ws_ami2);
	//get alignment path
	int i,j;
	for(i=0;i<moln1;i++)
	{
		int cur_index=i*moln2;
		for(j=0;j<moln2;j++)
		{
			TM_DP_sco[cur_index+j]=Universal_Calc(i,j,1,ws_ami1,ws_ami2,moln1,moln2,3);
		}
	}
	//apply dynamic programming
	double gap_open=-10.0;
	TM_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,ali2,gap_open);
	//delete
	delete [] ws_ami1;
	delete [] ws_ami2;
}
void CLEFAPS_Main::FM_Get_Initial_CLE(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	//assign SSE
	int *ws_cle1=new int[moln1+moln2];
	int *ws_cle2=new int[moln1+moln2];
	CLE_transform(IN_CLE1,ws_cle1);
	CLE_transform(IN_CLE2,ws_cle2);
	//get alignment path
	int i,j;
	for(i=0;i<moln1;i++)
	{
		int cur_index=i*moln2;
		for(j=0;j<moln2;j++)
		{
			TM_DP_sco[cur_index+j]=Universal_Calc(i,j,1,ws_cle1,ws_cle2,moln1,moln2,1);
		}
	}
	//apply dynamic programming
	double gap_open=-100.0;
	TM_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,ali2,gap_open);
	//delete
	delete [] ws_cle1;
	delete [] ws_cle2;
}
void CLEFAPS_Main::FM_Get_Initial2(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	//assign SSE
	int *ws_ami_cle1=new int[moln1+moln2];
	int *ws_ami_cle2=new int[moln1+moln2];
	AMI_CLE_transform_Ori(IN_AMI1,IN_CLE1,ws_ami_cle1);
	AMI_CLE_transform_Ori(IN_AMI2,IN_CLE2,ws_ami_cle2);
	//get alignment path
	int i,j;
	for(i=0;i<moln1;i++)
	{
		int cur_index=i*moln2;
		for(j=0;j<moln2;j++)
		{
			TM_DP_sco[cur_index+j]=Universal_Calc(i,j,1,ws_ami_cle1,ws_ami_cle2,moln1,moln2,4);
		}
	}
	//apply dynamic programming
	double gap_open=-100.0;
	TM_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,ali2,gap_open);
	//delete
	delete [] ws_ami_cle1;
	delete [] ws_ami_cle2;
}

//-------- Envo_Align -------//(obselete)
void CLEFAPS_Main::FM_Get_Initial3(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	//assign SSE
	Envo_Align_Neib_Main(mol1,mol2,moln1,moln2);  //-> using C-alpha to calculate neibors
	Input_Envo_Align_Total(mol1,mol2,moln1,moln2,IN_AMI1,IN_AMI2,IN_CLE1,IN_CLE2);
	//alignment
	vector <vector <double> > out_mat_cle;
	Envo_Align_Calc_Main(out_mat_cle,2);
	//get alignment path
	int i,j;
	double ali_sco;
	vector<pair<int,int> > alignment;
	for(i=0;i<moln1;i++)for(j=0;j<moln2;j++)TM_DP_sco[i*moln2+j]=out_mat_cle[i][j];
	Advance_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,0.0,0,alignment,ali_sco);
	int ii,jj;
	for(i=0;i<moln2;i++)ali2[i]=-1;
	int totnum=(int)alignment.size();
	for(i=0;i<totnum;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0&&jj>0)ali2[jj-1]=ii-1;
	}
}

//================= FM_Align_Series ================//__110720__//
void CLEFAPS_Main::Alimeng_Add(vector <Align_Record> &align_tot,int moln2,double tms)
{
	Align_Record tmp_ali;
	memcpy(ali2,TM_Main_Ali2,moln2*sizeof(int));
	memcpy(CLEP_rotmat,TM_rotmat,12*sizeof(double));
	CLEPAPS_Make_Score(tmp_ali);
	tmp_ali.main_sco=tms;
	//init
	int j;
	double comp_check;
	int success=1; //defualt:1
	int cursize=(int)align_tot.size();
	for(j=0;j<cursize;j++)
	{
		comp_check=Alg_Comp_Simp(tmp_ali.alignment,align_tot[j].alignment,Comp_Check_Vari);
		if(comp_check>Comp_Check)
		{
			success=0;
			break;
		}
	}
	if(success==1)align_tot.push_back(tmp_ali);
}
double CLEFAPS_Main::FM_Align_Total(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
									int norm_len,double norm_d0,double *MAXSCO)
{
	int k;
	double TM_CUR;
	double TM_BEST=0.0;
	//total_init
	AliT_Rec=0;  //default:NOT
	AliT_Max=0;  //default:NO
	memset(ali2_,-1,moln2*sizeof(int));
	memset(TM_Main_Ali2,-1,moln2*sizeof(int));

	//============ CLEFAPS ============//
	//CLEFAPS
//ws_init0:
	FM_align_tot.clear();
	TM_CUR=CLEFAPS_Main_Func(FM_align_tot);
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		for(k=0;k<moln2;k++)ali2_[k]=CLEP_MAX.alignment[k];
		for(k=0;k<12;k++)TM_FINMAT[k]=CLEP_MAX.rotmat[k];
		AliT_Max=0;
	}
	//CLE DynaProg
	FM_Get_Initial_CLE(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=8;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init5;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=8;
	}
	Alimeng_Add(FM_align_tot,moln2,TM_CUR);

	//--- get_initial5 ---//local structure superposition(neo)
ws_init5:
	TM_Get_Initial5(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=5;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init6;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=5;
	}
	Alimeng_Add(FM_align_tot,moln2,TM_CUR);

	//--- get_initial6 ---//local structure superposition(old)
ws_init6:
	if(TM_ADDITION==1)
	{
		TM_Get_Initial6(mol1,mol2,moln1,moln2,TM_Main_Ali2);
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=6;
		}
		if(TM_CUR<=TM_BEST*0.2)goto ws_init2;  //check bad...
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=6;
		}
		Alimeng_Add(FM_align_tot,moln2,TM_CUR);
	}

	//--- get_initial2 ---//SSE
ws_init2:
	TM_Get_Initial2(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=2;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init1;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=2;
	}
	Alimeng_Add(FM_align_tot,moln2,TM_CUR);

	//--- get_initial1 ---//gapless threading
ws_init1:
	TM_Get_Initial1(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=1;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init4;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=1;
	}
	Alimeng_Add(FM_align_tot,moln2,TM_CUR);

	//--- get_initial1.5 ---//gapless threading advance
	if(TM_ADDITION==1 && AliT_Rec==1)
	{
		EqualArray(TM_Main_Ali2,TM_Main_AliT,moln2);
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=7;
		}
		if(TM_CUR<=TM_BEST*0.2)goto ws_init4;  //check bad...
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=7;
		}
		Alimeng_Add(FM_align_tot,moln2,TM_CUR);
	}

	//--- get_initial4 ---//
ws_init4:
	TM_Get_Initial4(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=4;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init3;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=4;
	}
	Alimeng_Add(FM_align_tot,moln2,TM_CUR);

	//--- get_initial3 ---//best+SSE
ws_init3:
	memcpy(TM_Main_Ali2,ali2_,moln2*sizeof(int));  //assign the current best path!!
	TM_Get_Initial3(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=3;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_end;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=3;
	}
	Alimeng_Add(FM_align_tot,moln2,TM_CUR);

	//========= terminal =========//
ws_end:
	TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
	if(FM_Kill_Gap==1)CLEF_Single_Match_Kill(mol1,mol2,moln1,moln2,FM_align_tot);
	TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
	stable_sort(FM_align_tot.rbegin(),FM_align_tot.rend());   //descending sort
	for(k=0;k<moln2;k++)ali2_[k]=FM_align_tot[0].alignment[k];
	EqualArray(TM_FINMAT,finmat,12);
	TM_BEST=FM_align_tot[0].TMsco;
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);

	//---- final add ----//__2016.05.30__//
	{
		//kill
		vector <Align_Record> FM_align_tot_;
		for(int kk=0;kk<(int)FM_align_tot.size();kk++)
		{
			for(int i=0;i<moln2;i++)TM_Main_Ali2[i]=FM_align_tot[kk].alignment[i];
			for(int i=0;i<12;i++)TM_rotmat[i]=FM_align_tot[kk].rotmat[i];
			TM_CUR=FM_align_tot[kk].main_sco;
			Alimeng_Add(FM_align_tot_,moln2,TM_CUR);
		}
		//add
		FM_align_tot.clear();
		for(int kk=0;kk<(int)FM_align_tot_.size();kk++)
		{
			FM_align_tot.push_back(FM_align_tot_[kk]);
		}
	}

	//return
	return TM_BEST;
}

//=========== CLEFAPS_Normal ==========//__110715__//
double CLEFAPS_Main::FM_Align_Normal(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
									int norm_len,double norm_d0,double *MAXSCO)
{
	int k;
	double TM_CUR;
	double TM_BEST=0.0;

	//total_init
	AliT_Rec=0;  //default:NOT
	AliT_Max=0;  //default:NO
	memset(ali2_,-1,moln2*sizeof(int));
	memset(TM_Main_Ali2,-1,moln2*sizeof(int));

	//============ CLEFAPS ============//
	//CLEFAPS
//ws_init0:
	FM_align_tot.clear();
	TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
	TM_CUR=CLEFAPS_Main_Func(FM_align_tot);
	TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		for(k=0;k<moln2;k++)ali2_[k]=CLEP_MAX.alignment[k];
		for(k=0;k<12;k++)TM_FINMAT[k]=CLEP_MAX.rotmat[k];
		AliT_Max=0;
	}

	//---- ws_debug ---//
	int retval=(int)TM_CUR;
	int tm_origin=TM_ADDITION;
	if(retval==-1)TM_ADDITION=1;

	//--- get_initial5 ---//local structure superposition(neo)
	if(TM_ADDITION==1)
	{
//ws_init5:
		TM_Get_Initial5(mol1,mol2,moln1,moln2,TM_Main_Ali2);
		TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
		TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=5;
		}
		if(TM_CUR<=TM_BEST*Refine_Cut)goto ws_init1;  //check bad...
		TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1,REFINE_ITER); //run DynaProg
		TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=5;
		}
		Alimeng_Add(FM_align_tot,moln2,TM_CUR);


		//--- get_initial1 ---//gapless threading
ws_init1:
		TM_Get_Initial1(mol1,mol2,moln1,moln2,TM_Main_Ali2);
//		TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
//		TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=1;
		}
		if(TM_CUR<=TM_BEST*Refine_Cut)goto ws_init3;  //check bad...
//		TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1,REFINE_ITER); //run DynaProg
//		TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=1;
		}
		Alimeng_Add(FM_align_tot,moln2,TM_CUR);

		//--- get_initial1.5 ---//gapless threading advance
		if(AliT_Rec==1)
		{
			EqualArray(TM_Main_Ali2,TM_Main_AliT,moln2);
//			TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
			TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
//			TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
			if(TM_CUR>TM_BEST)
			{
				TM_BEST=TM_CUR;
				memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
				memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
				AliT_Max=7;
			}
			if(TM_CUR<=TM_BEST*Refine_Cut)goto ws_init3;  //check bad...
//			TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
			TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1,REFINE_ITER); //run DynaProg
//			TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
			if(TM_CUR>TM_BEST)
			{
				TM_BEST=TM_CUR;
				memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
				memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
				AliT_Max=7;
			}
			Alimeng_Add(FM_align_tot,moln2,TM_CUR);
		}

		//--- get_initial3 ---//best+SSE
ws_init3:
		memcpy(TM_Main_Ali2,ali2_,moln2*sizeof(int));  //assign the current best path!!
		TM_Get_Initial3(mol1,mol2,moln1,moln2,TM_Main_Ali2);
		TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
		TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=3;
		}
		if(TM_CUR<=TM_BEST*Refine_Cut)goto ws_end;  //check bad...
		TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1,REFINE_ITER); //run DynaProg
		TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=3;
		}
		Alimeng_Add(FM_align_tot,moln2,TM_CUR);
	}

	//========= terminal =========//
ws_end:
	//---- ws_debug ---//
	if(retval==-1)TM_ADDITION=tm_origin;
	if(FM_align_tot.size()==0)return TM_BEST;
	//---- ws_real_emd -//
	TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
	if(FM_Kill_Gap==1)CLEF_Single_Match_Kill(mol1,mol2,moln1,moln2,FM_align_tot);
	TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
	stable_sort(FM_align_tot.rbegin(),FM_align_tot.rend());   //descending sort
	for(k=0;k<moln2;k++)ali2_[k]=FM_align_tot[0].alignment[k];
	EqualArray(TM_FINMAT,finmat,12);
	TM_BEST=FM_align_tot[0].TMsco;
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);

	//---- final add ----//__2016.05.30__//
	{
		//kill
		vector <Align_Record> FM_align_tot_;
		for(int kk=0;kk<(int)FM_align_tot.size();kk++)
		{
			for(int i=0;i<moln2;i++)TM_Main_Ali2[i]=FM_align_tot[kk].alignment[i];
			for(int i=0;i<12;i++)TM_rotmat[i]=FM_align_tot[kk].rotmat[i];
			TM_CUR=FM_align_tot[kk].main_sco;
			Alimeng_Add(FM_align_tot_,moln2,TM_CUR);
		}
		//add
		FM_align_tot.clear();
		for(int kk=0;kk<(int)FM_align_tot_.size();kk++)
		{
			FM_align_tot.push_back(FM_align_tot_[kk]);
		}
	}

	//return
	return TM_BEST;
}


//=========== CLEFAPS With Alignment ==========//
double CLEFAPS_Main::FM_Align_WithAli(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
									int norm_len,double norm_d0,double *MAXSCO)
{
	int k;
	double TM_CUR;
	double TM_BEST=0.0;
	double tms;
	double rms;
	int lali;

	//total_init
	AliT_Rec=0;  //default:NOT
	AliT_Max=0;  //default:NO

	//============ CLEFAPS ============//
	FM_align_tot.clear();
	{
		memcpy(TM_Main_Ali2,ali2_,moln2*sizeof(int));  //assign the current best path!!
		tms=TM_Align_TM_Score_Simplest(mol1,mol2,moln1,moln2,TM_Main_Ali2,norm_len,norm_d0,rms,lali);
		rot_mol(mol1,TM_tmp1,moln1,finmat);
		memcpy(TM_rotmat,finmat,12*sizeof(double));
		//get bound and dyna_sco
		TM_Align_Get_Ali2_Bound(moln1,moln2,TM_Main_Ali2,TM_bound,TM_bound_neib);
		TM_Align_Get_Matrix_Bound(TM_tmp1,mol2,moln1,moln2,norm_d0,TM_bound,TM_DP_sco);
		TM_CUR=TM_Align_Get_Score_Simp(mol1,mol2,finmat,moln1,moln2,TM_Main_Ali2);  //DeepScore type
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=3;
		}
		TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
		//TM_NOTM=1;
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1,REFINE_ITER); //run DynaProg
		//TM_NOTM=0;
		TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=3;
		}
		Alimeng_Add(FM_align_tot,moln2,TM_CUR);
	}

	//========= terminal =========//
	TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
	if(FM_Kill_Gap==1)CLEF_Single_Match_Kill(mol1,mol2,moln1,moln2,FM_align_tot);
	TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
	stable_sort(FM_align_tot.rbegin(),FM_align_tot.rend());   //descending sort
	for(k=0;k<moln2;k++)ali2_[k]=FM_align_tot[0].alignment[k];
	EqualArray(TM_FINMAT,finmat,12);
	TM_BEST=FM_align_tot[0].TMsco;
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);

	//---- final add ----//__2016.05.30__//
	{
		//kill
		vector <Align_Record> FM_align_tot_;
		for(int kk=0;kk<(int)FM_align_tot.size();kk++)
		{
			for(int i=0;i<moln2;i++)TM_Main_Ali2[i]=FM_align_tot[kk].alignment[i];
			for(int i=0;i<12;i++)TM_rotmat[i]=FM_align_tot[kk].rotmat[i];
			TM_CUR=FM_align_tot[kk].main_sco;
			Alimeng_Add(FM_align_tot_,moln2,TM_CUR);
		}
		//add
		FM_align_tot.clear();
		for(int kk=0;kk<(int)FM_align_tot_.size();kk++)
		{
			FM_align_tot.push_back(FM_align_tot_[kk]);
		}
	}

	//return
	return TM_BEST;
}


//=========== CLEFAPS_Lite ==========//
double CLEFAPS_Main::FM_Align_Lite(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
									int norm_len,double norm_d0,double *MAXSCO)
{
	int k;
	double TM_CUR;
	double TM_BEST=0.0;

	//total_init
	AliT_Rec=0;  //default:NOT
	AliT_Max=0;  //default:NO
	memset(ali2_,-1,moln2*sizeof(int));
	memset(TM_Main_Ali2,-1,moln2*sizeof(int));

	//============ CLEFAPS ============//
	//CLEFAPS
//ws_init0:
	FM_align_tot.clear();
	TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
	TM_CUR=CLEFAPS_Main_Func(FM_align_tot);
	TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		for(k=0;k<moln2;k++)ali2_[k]=CLEP_MAX.alignment[k];
		for(k=0;k<12;k++)TM_FINMAT[k]=CLEP_MAX.rotmat[k];
		AliT_Max=0;
	}

	//----- final process ----//
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);
	return TM_BEST;
}

//=========== CLEFAPS_Fast ==========//
double CLEFAPS_Main::FM_Align_Fast(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,double *MAXSCO)
{
	int k;
	double TM_CUR;
	double TM_BEST=0.0;

	//total_init
	AliT_Rec=0;  //default:NOT
	AliT_Max=0;  //default:NO
	for(k=0;k<moln2;k++)ali2[k]=-1;
	for(k=0;k<moln2;k++)TM_Main_Ali2[k]=-1;
	CLEFAPS_Main_Init(1);

	//============ CLEFAPS ============//
	//CLE DynaProg
	FM_Get_Initial2(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=8;
	}
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1,2); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		EqualArray(ali2,TM_Main_Ali2,moln2);
		EqualArray(TM_FINMAT,TM_rotmat,12);
		AliT_Max=8;
	}
	return TM_CUR;
}
