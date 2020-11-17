#include "Hydro_Bond.h"

//--------------------- Start -----------------//
Hydro_Bond::Hydro_Bond(int num)
{
	HB_MAXIMAL=num;
	HB_Init(HB_MAXIMAL);
}
Hydro_Bond::~Hydro_Bond(void)
{
	HB_Dele(HB_MAXIMAL);
}

//-------------------- init -------------------//
void Hydro_Bond::HB_Init(int HB_MAXIMAL)
{
	NewArray2D(&HB_mol,HB_MAXIMAL,5);
	HB_sse=new char[HB_MAXIMAL+1];
	HB_ami=new char[HB_MAXIMAL+1];
	HB_broke=new int[HB_MAXIMAL];
	HB_simple=new int[HB_MAXIMAL*2];
	HB_index=new int[HB_MAXIMAL*5];
	HB_bridge=new int[HB_MAXIMAL*2];
	HB_ladder=new int[HB_MAXIMAL*2*4];
	HB_ladder_rec=new int[HB_MAXIMAL*2];
	HB_temp_index=new int[HB_MAXIMAL];
}
void Hydro_Bond::HB_Dele(int HB_MAXIMAL)
{
	DeleteArray2D(&HB_mol,HB_MAXIMAL);
	delete [] HB_sse;
	delete [] HB_ami;
	delete [] HB_broke;
	delete [] HB_simple;
	delete [] HB_index;
	delete [] HB_bridge;
	delete [] HB_ladder;
	delete [] HB_ladder_rec;
	delete [] HB_temp_index;
}

//-------------------- calc hydro_bond-------------------//
void Hydro_Bond::HB_Input_Mol(XYZ **mol,char *ami,int moln)
{
	int i;
	XYZ xyz;
	double dist;
	//input
	HB_moln=moln;
	EqualArray2D(HB_mol,mol,moln,4); //only consider N,CA,C,O
	strcpy(HB_ami,ami);

	//----- create H atom ------//
	for(i=0;i<HB_moln-1;i++)
	{
		//---- Modify ----//__100718__//
		//proline has no backbone hydrogen atom!!
		//so we may neglect the dihetral angle calculation!!
		if(HB_ami[i+1]!='P')
		{
			xyz=HB_mol[i][2]-HB_mol[i][3];
			dist=HB_mol[i][2].distance(HB_mol[i][3]);
			xyz/=dist;
			HB_mol[i+1][4]=HB_mol[i+1][0]+xyz;
		}
		else HB_mol[i+1][4]=HB_mol[i+1][0];
	}
	//first H
	HB_mol[0][4]=HB_mol[0][0]+HB_mol[1][4]-HB_mol[1][0];

	//---- calculate broke ----//
	HB_broke[0]=0;
	for(i=0;i<HB_moln-1;i++)
	{
		//copied from DSSP // Distance(chain[*LINK->lchain].c, LINK->resinfo.n) > BREAKDIST==2.5
		dist=HB_mol[i][2].distance_square(HB_mol[i+1][0]);
		if(dist>6.25)
		{
			HB_broke[i+1]=1;
		}
		else HB_broke[i+1]=0;
	}
}
double Hydro_Bond::HB_Calc_Single(XYZ **HB_mol,int i,int j)
{
	double dho,dhc,dnc,dno;
	dho=HB_mol[i][4].distance(HB_mol[j][3]);
	dhc=HB_mol[i][4].distance(HB_mol[j][2]);
	dno=HB_mol[i][0].distance(HB_mol[j][3]);
	dnc=HB_mol[i][0].distance(HB_mol[j][2]);
	return (332.0*0.42*0.2*(1.0/dno+1.0/dhc-1.0/dho-1.0/dnc));
}
void Hydro_Bond::HB_Calc_Hydro_Bond(XYZ **HB_mol,int HB_moln,int *HB_simple)
{
	int i,j;
	double value;
	double max1,max2;
	int imax1,imax2;
	for(i=0;i<HB_moln;i++)
	{
		//calc value
		max1=0;
		max2=0;
		imax1=-1;
		imax2=-1;
		for(j=0;j<HB_moln;j++)
		{
			//judge
			if(j==i||j==i-1)continue;
			if(HB_mol[i][1].distance_square(HB_mol[j][1])>81.0)continue;
			//hydrogen bond
			value=HB_Calc_Single(HB_mol,i,j);
			if(value<max1)
			{
				max2=max1;
				max1=value;
				imax2=imax1;
				imax1=j;				
			}
			else if(value<max2)
			{
				max2=value;
				imax2=j;
			}
		}
		//assign
		if(imax1!=-1&&max1<-0.5005)HB_simple[i*2+0]=imax1;
		else HB_simple[i*2+0]=-1;
		if(imax2!=-1&&max2<-0.5005)HB_simple[i*2+1]=imax2;
		else HB_simple[i*2+1]=-1;		
	}
}
void Hydro_Bond::HB_Calc_Hydro_Bond(void)
{
	HB_Calc_Hydro_Bond(HB_mol,HB_moln,HB_simple);
}

//---- calculate HydroBond matrix -----//__2018.05.20__//
void Hydro_Bond::HB_Calc_Hydro_Bond(vector <vector <double> > &hb_mat)
{
	int i,j;
	double value;
	//create
	hb_mat.resize(HB_moln);
	for(i=0;i<HB_moln;i++)hb_mat[i].resize(HB_moln);
	//calculate
	for(i=0;i<HB_moln;i++)
		for(j=0;j<HB_moln;j++)
		{
			//-> calc
			if(j==i||j==i-1)value=0;
			else value=HB_Calc_Single(HB_mol,i,j);
			//-> push
			hb_mat[i][j]=value;
		}
}


//------------------------ calc secondary structure ----------------//
//[check broken chain!]
int Hydro_Bond::HB_Broke_Check(int i,int j)
{
	for(int k=i+1;k<=j;k++)if(HB_broke[k]==1)return 1;
	return 0;
}
//[calculate helix and hydro-turn]
void Hydro_Bond::HB_Calc_SSE_Helix(int pitch,char h,char t)
{
	int i,k;
	int pos1,pos2;
	int start,end;
	int first;
	int len;
	int rellen;
	start=-1;
	first=1; //default:OK
	len=0;
	for(i=0;i<HB_moln;i++)
	{
		pos1=HB_simple[i*2+0];
		pos2=HB_simple[i*2+1];
		if(pos1!=-1&&i==pos1+pitch) //satisfy
		{
			if(HB_Broke_Check(pos1,i)==1)goto end;
			if(first==1)
			{
				first=0;
				start=pos1+1;
				len=0;
			}
			else len++;
		}
		else if(pos2!=-1&&i==pos2+pitch) //satisfy
		{
			if(HB_Broke_Check(pos2,i)==1)goto end;
			if(first==1)
			{
				first=0;
				start=pos2+1;
				len=0;
			}
			else len++;
		}
		else
		{
end:
			if(first==0)
			{
				first=1;
				end=i-1;
				//check
				rellen=0;
				for(k=start;k<end;k++)
				{
					if(HB_sse[k]==' '||HB_sse[k]==h||HB_sse[k]==t)rellen++;
				}
				//assign
				if(rellen>=pitch)
				{
					for(k=start;k<end;k++)
					{
						if(HB_sse[k]==' '||HB_sse[k]==t)HB_sse[k]=h;
					}
				}
				else
				{
					for(k=start;k<end;k++)
					{
						if(HB_sse[k]==' ')HB_sse[k]=t;
					}
				}
				len=0;
			}
		}
	}
	//terminal
	if(first==0)
	{
		first=1;
		end=i-1;
		//check
		rellen=0;
		for(k=start;k<end;k++)
		{
			if(HB_sse[k]==' '||HB_sse[k]==h||HB_sse[k]==t)rellen++;
		}
		//assign
		if(rellen>=pitch)
		{
			for(k=start;k<end;k++)
			{
				if(HB_sse[k]==' '||HB_sse[k]==t)HB_sse[k]=h;
			}
		}
		else
		{
			for(k=start;k<end;k++)
			{
				if(HB_sse[k]==' ')HB_sse[k]=t;
			}
		}
		len=0;
	}
}

//[calculate helix and sheet-serials]
void Hydro_Bond::HB_Calc_SSE_Sheet(char e,char h)
{
	int i,j;
	int index;
	int p1,p2,q1,q2;
	int type;	
	//init
	for(i=0;i<HB_moln;i++)HB_index[i*5+0]=0;
	HB_bridge[0]=0;
	//calc
	for(i=1;i<HB_moln-1;i++)
	{
		for(j=i+1;j<HB_moln-1;j++)
		{
			type=0;
			//parallel[1]
			p1=HB_simple[(i+1)*2+0];
			p2=HB_simple[(i+1)*2+1];
			q1=HB_simple[j*2+0];
			q2=HB_simple[j*2+1];
			if((p1==j||p2==j)&&(q1==i-1||q2==i-1))
			{
				type=1;
				goto sheet_assign;
			}
			//parallel[2]
			p1=HB_simple[i*2+0];
			p2=HB_simple[i*2+1];
			q1=HB_simple[(j+1)*2+0];
			q2=HB_simple[(j+1)*2+1];
			if((p1==j-1||p2==j-1)&&(q1==i||q2==i))
			{
				type=2;
				goto sheet_assign;
			}
			//anti-parallel[1]
			p1=HB_simple[(i+1)*2+0];
			p2=HB_simple[(i+1)*2+1];
			q1=HB_simple[(j+1)*2+0];
			q2=HB_simple[(j+1)*2+1];
			if((p1==j-1||p2==j-1)&&(q1==i-1||q2==i-1))
			{
				type=-1;
				goto sheet_assign;
			}
			//anti-parallel[2]
			p1=HB_simple[i*2+0];
			p2=HB_simple[i*2+1];
			q1=HB_simple[j*2+0];
			q2=HB_simple[j*2+1];
			if((p1==j||p2==j)&&(q1==i||q2==i))
			{
				type=-2;
				goto sheet_assign;
			}
			
			//assign
sheet_assign:			
			if(type!=0&&abs(i-j)>2)
			{
				if(HB_broke[i+1]==0&&HB_broke[j+1]==0&&HB_broke[i]==0&&HB_broke[j]==0)  //terminal_check//__100725__//
				{
					if(HB_sse[i]!=h)
					{
						HB_index[i*5+0]++;
						index=HB_index[i*5+0];
						if(index==5)continue;
						HB_index[i*5+index]=j;
						HB_sse[i]=e;
					}
					if(HB_sse[j]!=h)
					{
						HB_index[j*5+0]++;
						index=HB_index[j*5+0];
						if(index==5)continue;
						HB_index[j*5+index]=i;
						HB_sse[j]=e;
					}
					HB_bridge[0]++;
					index=HB_bridge[0];
					HB_bridge[index]=abs(type)*INT_MAX_NUM*PROT_MAX_NUM+i*INT_MAX_NUM+j;
					if(type<0)HB_bridge[index]*=-1;
				}
			}
		}
	}
}
void Hydro_Bond::HB_Calc_SSE_Sheet_Ladder(void)
{
	int i,j;
	int ii,jj,len;
	int cur_ii,cur_jj;
	int type1,type2;
	int rtype1,rtype2;
	int totnum;
	int index;
	int value;
	totnum=HB_bridge[0];
	//init
	HB_ladder[0*4+0]=0;
	for(i=0;i<=totnum;i++)HB_temp_index[i]=1; //default:OK
	//calc
	for(i=1;i<=totnum;i++)
	{
		if(HB_temp_index[i]==0)continue;
		value=HB_bridge[i];
		if(value>0)
		{
			rtype1=1;
			type1=value/(PROT_MAX_NUM*INT_MAX_NUM);
			value=value%(PROT_MAX_NUM*INT_MAX_NUM);
			ii=value/INT_MAX_NUM;
			jj=value%INT_MAX_NUM;
		}
		else
		{
			rtype1=-1;
			value*=-1;
			type1=value/(PROT_MAX_NUM*INT_MAX_NUM);
			value=value%(PROT_MAX_NUM*INT_MAX_NUM);
			ii=value/INT_MAX_NUM;
			jj=value%INT_MAX_NUM;
		}
		HB_temp_index[i]=0;
		len=1;
		for(j=i+1;j<=totnum;j++)
		{
			if(HB_temp_index[j]==0)continue;
			value=HB_bridge[j];
			if(value>0)
			{
				rtype2=1;
				type2=value/(PROT_MAX_NUM*INT_MAX_NUM);
				value=value%(PROT_MAX_NUM*INT_MAX_NUM);
				cur_ii=value/INT_MAX_NUM;
				cur_jj=value%INT_MAX_NUM;
			}
			else
			{
				rtype2=-1;
				value*=-1;
				type2=value/(PROT_MAX_NUM*INT_MAX_NUM);
				value=value%(PROT_MAX_NUM*INT_MAX_NUM);
				cur_ii=value/INT_MAX_NUM;
				cur_jj=value%INT_MAX_NUM;
			}
			if(rtype1!=rtype2)continue;
			if(rtype1==1) //parallel
			{
				if(cur_ii==ii+len&&cur_jj==jj+len)
				{
					len++;
					HB_temp_index[j]=0;
				}
				else continue;
			}
			if(rtype1==-1) //anti-parallel
			{
				if(cur_ii==ii+len&&cur_jj==jj-len)
				{
					len++;
					HB_temp_index[j]=0;
				}
				else continue;
			}
		}
		HB_ladder[0*4+0]++;
		index=HB_ladder[0*4+0];
		HB_ladder[index*4+0]=type1*rtype1;
		HB_ladder[index*4+1]=ii;
		HB_ladder[index*4+2]=jj;
		HB_ladder[index*4+3]=len;
	}
}
void Hydro_Bond::HB_Calc_SSE_Sheet_Bulge(char ef)
{
	int i,j,k;
	int type1,type2;
	int ii1,jj1,len1,ii2,jj2,len2;
	int p1,p2,q1,q2;  //p: one chain, q: other chain, 1: pre, 2: nxt
	int totnum;
	int bulge;
	int linkage;
	int link;
	totnum=HB_ladder[0*4+0];
	//init
	linkage=0;
	for(i=0;i<=totnum;i++)HB_ladder_rec[i]=0;
	for(i=1;i<=totnum;i++)
	{
		//get
		type1=HB_ladder[i*4+0];
		ii1=HB_ladder[i*4+1];
		jj1=HB_ladder[i*4+2];
		len1=HB_ladder[i*4+3];
		bulge=0;
		if(HB_ladder_rec[i]==0)
		{
			linkage++;
			link=linkage;
		}
		else link=HB_ladder_rec[i];
		//search
		for(j=i+1;j<=totnum;j++)
		{
			type2=HB_ladder[j*4+0];
			ii2=HB_ladder[j*4+1];
			jj2=HB_ladder[j*4+2];
			len2=HB_ladder[j*4+3];
			//calc
			if(type1*type2<0)continue;
			if(type1>0) //parallel
			{
				p1=ii1+len1-1;
				p2=ii2;
				q1=jj1+len1-1;
				q2=jj2;
			}
			if(type1<0) //anti-parallel
			{
				p1=ii1+len1-1;
				p2=ii2;
				q1=jj1-len1+1;
				q2=jj2;
			}
			//judge
			if(type1>0) //parallel
			{
				if(p1>p2||q1>q2)continue;
				bulge=0;
				if((p2-p1<6&&q2-q1<3)||(p2-p1<3&&q2-q1<6))
				{
					bulge=1;
					break;
				}
			}
			if(type1<0) //anti-parallel
			{
				if(p1>p2||q2>q1)continue;
				bulge=0;
				if((p2-p1<6&&q1-q2<3)||(p2-p1<3&&q1-q2<6))
				{
					bulge=1;
					break;
				}
			}
		}
		//assign
		if(bulge==1)
		{
			if(type1>0) //parallel
			{
				for(k=p1+1;k<p2;k++)HB_sse[k]=ef;
				for(k=q1+1;k<q2;k++)HB_sse[k]=ef;
			}
			if(type1<0) //anti-parallel
			{
				for(k=p1+1;k<p2;k++)HB_sse[k]=ef;
				for(k=q2+1;k<q1;k++)HB_sse[k]=ef;
			}
			HB_ladder_rec[i]=link;
			HB_ladder_rec[j]=link;
		}
	}
}
void Hydro_Bond::HB_Calc_SSE_Sheet_Lone(char ef,char eb)
{
	int i;
	int ii,jj,len;
	int link;
	int pos1,pos2;
	int totnum;
	totnum=HB_ladder[0*4+0];
	for(i=1;i<=totnum;i++)
	{
		//get
		ii=HB_ladder[i*4+1];
		jj=HB_ladder[i*4+2];
		len=HB_ladder[i*4+3];
		link=HB_ladder_rec[i];
		if(link!=0)continue;
		if(len>1)continue;
		//assign
		pos1=HB_index[ii*5+0];
		pos2=HB_index[jj*5+0];
		if(pos1==1&&HB_sse[ii]!=ef)HB_sse[ii]=eb;
		if(pos2==1&&HB_sse[jj]!=ef)HB_sse[jj]=eb;
	}
}
void Hydro_Bond::HB_Calc_SSE_Sheet_Return(char ef,char e)
{
	int i;
	for(i=0;i<HB_moln;i++)
	{
		if(HB_sse[i]==ef)HB_sse[i]=e;
	}
}

//[calculate kappa-turn]
void Hydro_Bond::HB_Calc_SSE_Turn(char s)
{
	int i,j;
	int correct;
	XYZ a,b;
	double aa[3],bb[3];
	double wa,wb;
	double kappa;
	for(i=2;i<HB_moln-2;i++)
	{
		correct=1; //default:OK
		//first_check
		for(j=1;j<5;j++)
		{
			if(HB_broke[i-2+j]==1)
			{
				correct=0;
				break;
			}
		}
		//kappa-turn
		if(HB_sse[i]==' '&&correct==1)
		{
			a=HB_mol[i][1]-HB_mol[i-2][1];
			b=HB_mol[i+2][1]-HB_mol[i][1];
			a.xyz2double(aa);
			b.xyz2double(bb);
			wa=dot(aa,aa);
			wb=dot(bb,bb);
			kappa=dot(aa,bb);
			kappa/=sqrt(1.0*wa*wb);
			if(kappa<0.34202014332567)HB_sse[i]=s;
		}
	}
}
//[kill singler]
void Hydro_Bond::HB_Calc_SSE_Lone(char ori,char cur,int limit)
{
	int i,k;
	int start;
	int first;
	int len;
	start=-1;
	first=1;
	len=0;
	for(i=0;i<HB_moln;i++)
	{
		if(HB_sse[i]==ori)
		{
			if(first==1)
			{
				start=i;
				first=0;
				len=1;
			}
			else len++;
		}
		else
		{
			if(first==0)
			{
				first=1;
				if(len<limit)
				{
					for(k=0;k<len;k++)HB_sse[start+k]=cur;
				}
			}
		}
	}
	//terminal
	if(first==0)
	{
		if(len<limit)
		{
			for(k=0;k<len;k++)HB_sse[start+k]=cur;
		}
	}
}

//[calculate SSE main]
void Hydro_Bond::HB_Calc_SSE(char *sse)
{
	int i;
	//init
	for(i=0;i<HB_moln;i++)HB_sse[i]=' ';
	HB_sse[i]='\0';
	//[helix]
	HB_Calc_SSE_Helix(4,'H','T');
	//[sheet]
	HB_Calc_SSE_Sheet('E','H');
	HB_Calc_SSE_Sheet_Ladder();
	HB_Calc_SSE_Sheet_Bulge('F');
	HB_Calc_SSE_Sheet_Lone('F','B');
	HB_Calc_SSE_Sheet_Return('F','E');
	HB_Calc_SSE_Lone('E','B',2); //new_add//__100725__//
	//[helix/turn]
	HB_Calc_SSE_Helix(3,'G','T');
	HB_Calc_SSE_Helix(5,'I','T');
	//[kappa/turn]
	HB_Calc_SSE_Turn('S');
	//final
	for(i=0;i<HB_moln;i++)if(HB_sse[i]==' ')HB_sse[i]='L';
	strcpy(sse,HB_sse);
}

//------ transform from 8-digit SSE to 3-digit SSE ----------//
/*
Please refer to the following paper:

"Evaluation and improvement of multiple sequence methods for protein secondary structure prediction"
by James A. Cuff and Geoffrey J. Barton
at 1999, on PROTEINS: Structure, Function, and Genetics 34:508-519

Method A: E,B to E; G,H to H; rest to C
Method B: E to E; H to H; rest to C

The prediction results show that Method B could get better performance than Method A
*/
char Hydro_Bond::HB_Trans_SSE_Single_A(char c)
{
	switch(c)
	{
		case 'G':return 'H'; //H
		case 'H':return 'H'; //H
		case 'I':return 'H'; //H
		case 'E':return 'E'; //E
		case 'B':return 'E'; //E
		default:return 'C';
	}
}
char Hydro_Bond::HB_Trans_SSE_Single_B(char c)
{
	switch(c)
	{
		case 'H':return 'H'; //H
		case 'E':return 'E'; //E
		default:return 'C';
	}
}

//----------- HB_Trans_SSE ----------//
void Hydro_Bond::HB_Trans_SSE(char *in,char *out,int moln,int method)
{
	int i;
	if(method==0)for(i=0;i<moln;i++)out[i]=HB_Trans_SSE_Single_A(in[i]);
	else for(i=0;i<moln;i++)out[i]=HB_Trans_SSE_Single_B(in[i]);
}

