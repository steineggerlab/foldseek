#include "MultiAlign_Out.h"

//--------------------- Start -----------------//
MultiAlign_Out::MultiAlign_Out(void)
{
}
MultiAlign_Out::~MultiAlign_Out(void)
{
}

//------------ create & delete --------------//
void MultiAlign_Out::BC_Output_Init(int num,int len)
{
	//main
	NewArray3D(&Multi_AFB,len,num,2);
	Multi_AFB_Record=new int[len];
	Real_Block=new int[len];
	//temp
	BO_wsrec=new int[num];
	BO_wsbak=new int[num];
	BO_ZeroOne_Cur=new int[num];
}
void MultiAlign_Out::BC_Output_Dele(int num,int len)
{
	//main
	DeleteArray3D(&Multi_AFB,len,num);
	delete [] Multi_AFB_Record;
	delete [] Real_Block;
	//temp
	delete [] BO_wsrec;
	delete [] BO_wsbak;
	delete [] BO_ZeroOne_Cur;
}

//--------- process -----------//
int MultiAlign_Out::BC_Output_PDB(FILE *fp,int len,XYZ *mol,char *AMI,char *ind,char Chain_ID) //output normal PDB file (CA only)
{
	int i,j;
	string buf;
	const char *dummyaa;
	char dummychn;
	char out[7];
	int atomnum;

	//judge
	if(mol==NULL)return -1;
	//process
	buf="                          ";
	atomnum=1;
	for(i=0;i<len;i++)
	{
		//init
		if(AMI==NULL)dummyaa=One2Three_III(0);
		else dummyaa=One2Three_III(AMI[i]);
		if(Chain_ID==0)dummychn='_';
		else dummychn=Chain_ID;

		//process
		if(ind==NULL)
		{
			fprintf(fp,"ATOM  %5d  CA  %3s %c%4d    %8.3f%8.3f%8.3f%s\n",atomnum,
			dummyaa,dummychn,i+1,mol[i].X,mol[i].Y,mol[i].Z,buf.c_str());
		}
		else
		{
			for(j=0;j<6;j++)out[j]=ind[6*i+j];
			out[j]='\0';
			fprintf(fp,"ATOM  %5d  CA  %3s %6s   %8.3f%8.3f%8.3f%s\n",atomnum,
			dummyaa,out,mol[i].X,mol[i].Y,mol[i].Z,buf.c_str());
		}
		atomnum++;
	}
	//terminal
	return 1;
}
void MultiAlign_Out::BC_Ali_To_AFB(int **ali,int totlen,int totnum,int ***AFB_Out,int *AFB_Len,int *Real_Core)
{
	int i,j;
	int pre,cur;
	int valid;
	int first;
	int reclen;
	int curlen=0;

	//init
	int *wsrec=BO_wsrec;
	int *wsbak=BO_wsbak;
	//process
	first=1;
	reclen=0;
	for(j=0;j<totnum;j++)wsbak[j]=0;
	pre=0;
	//start
	for(j=0;j<totnum;j++)BO_ZeroOne_Cur[j]=-1;
	for(i=0;i<totlen;i++)
	{
		valid=1; //default:OK
		//collect
		cur=0;
		for(j=0;j<totnum;j++)
		{
			if(ali[j][i]==2) //aligned
			{
				wsrec[j]=1;
				BO_ZeroOne_Cur[j]++;
				cur++;
			}
			else if(ali[j][i]==1) //not aligned
			{
				wsrec[j]=0;
				BO_ZeroOne_Cur[j]++;
			}
			else //gap
			{
				wsrec[j]=0;
			}
		}
		if(cur!=pre)
		{
			valid=0;
			goto next;
		}
		//check
		for(j=0;j<totnum;j++)
		{
			if(wsrec[j]!=wsbak[j])
			{
				valid=0;
				break;
			}
		}
next:
		//judge
		if(valid==0) //neo_appear
		{
			if(first==1)
			{
				if(cur>1)
				{
					//assign neo
					curlen++;
					for(j=0;j<totnum;j++)
					{
						if(wsrec[j]==1)
						{
							AFB_Out[curlen][j][0]=0;
							AFB_Out[curlen][j][1]=BO_ZeroOne_Cur[j];
						}
						else
						{
							AFB_Out[curlen][j][0]=-1;
							AFB_Out[curlen][j][1]=-1;
						}
					}
					if(cur==totnum)Real_Core[curlen]=cur;
					else Real_Core[curlen]=-1*cur;
					first=0;
					reclen=1;
					//record old
					EqualArray(wsbak,wsrec,totnum);
					pre=cur;
				}
			}
			else
			{
				if(reclen>0)
				{
					AFB_Len[curlen]=reclen;
					reclen=0;
				}
				if(cur>1)
				{
					//assign neo
					curlen++;
					for(j=0;j<totnum;j++)
					{
						if(wsrec[j]==1)
						{
							AFB_Out[curlen][j][0]=0;
							AFB_Out[curlen][j][1]=BO_ZeroOne_Cur[j];
						}
						else
						{
							AFB_Out[curlen][j][0]=-1;
							AFB_Out[curlen][j][1]=-1;
						}
					}
					if(cur==totnum)Real_Core[curlen]=cur;
					else Real_Core[curlen]=-1*cur;
					first=0;
					reclen=1;
					//record old
					EqualArray(wsbak,wsrec,totnum);
					pre=cur;
				}
				else
				{
					first=1;
					reclen=0;
					for(j=0;j<totnum;j++)wsbak[j]=0;
				}
			}
		}
		else //continue
		{
			if(cur>1)reclen++;
		}
	}
	//final
	if(reclen>0)AFB_Len[curlen]=reclen;
	AFB_Out[0][0][0]=curlen;
}

//============================= [output_related] ===========================//
void MultiAlign_Out::BC_Output_Alignment(FILE *fp,char **in,int totnum,int totlen,int **ali)
{
	int j,k;
	int tag;
	int count;
	for(j=0;j<totnum;j++)
	{
		count=0;
		for(k=0;k<totlen;k++)
		{
			tag=ali[j][k];
			if(tag==0) //gap
			{
				fprintf(fp,"-");
			}
			else if(tag==1) //small
			{
				fprintf(fp,"%c",in[j][count]-'A'+'a');
				count++;
			}
			else //capital
			{
				fprintf(fp,"%c",in[j][count]);
				count++;
			}
		}
		fprintf(fp,"\n");
	}
}
void MultiAlign_Out::BC_Output_Superimpose(FILE *fp,XYZ **in,char **ami,int *len,int totnum)
{
	int ori;
	int j;
	string TER="TER                                                                             ";
	string END="END                                                                             ";
	//output
	ori=0;
	for(j=0;j<totnum;j++)	
	{
		if(j<26) BC_Output_PDB(fp,len[j],in[j],ami[j],0,'A'+j);
		else if(j<52) BC_Output_PDB(fp,len[j],in[j],ami[j],0,'a'+j-26);
		else if(j<62) BC_Output_PDB(fp,len[j],in[j],ami[j],0,'0'+j-52);
		else BC_Output_PDB(fp,len[j],in[j],ami[j],0,'-');
		fprintf(fp,"%s\n",TER.c_str());  
	}
	fprintf(fp,"%s\n",END.c_str());  
}
void MultiAlign_Out::BC_Output_RasMol_Script(FILE *fws,
	int TOT_NUM,int ***Multi_AFB,int *Multi_AFB_Record,int *Real_Block)
{
	int j,k;
	int jj;
	int totnum;
	int winlen;

	fprintf(fws,"wireframe off\n");
	fprintf(fws,"backbone 10\n");
	fprintf(fws,"set ambient 20\n");
	fprintf(fws,"set background white\n");

	totnum=Multi_AFB[0][0][0];
	for(k=1;k<=totnum;k++)
	{
		for(j=0;j<TOT_NUM;j++)
		{
			//record
			jj=Multi_AFB[k][j][1];
			if(jj==-1)continue;
			//get chain
			char chain;
			{
				int pos=j;
				if(pos>=0 && pos<=25)chain=pos+'A';
				else if(pos>=26 && pos<=51)chain=pos-26+'a';
				else if(pos>=52 && pos<=61)chain=pos-52+'1';
				else chain='0';
			}
			winlen=Multi_AFB_Record[k];
			fprintf(fws,"select   %d-%d:%c\n", jj+1,jj+winlen,j+'A');
			//__070908__//judge whether is Real_Core or not
			if(Real_Block[k]<0)
			{
				fprintf(fws,"backbone off\n");
				fprintf(fws,"strands\n");
			}
			else
			{
				fprintf(fws,"backbone off\n");
				fprintf(fws,"cartoon\n");
			}
		}
	}
	fprintf(fws,"select all\n");
	fprintf(fws,"color chain\n");
}
void MultiAlign_Out::BC_Output_JMol_Script(FILE *fws,
	int TOT_NUM,int ***Multi_AFB,int *Multi_AFB_Record,int *Real_Block)
{
	int j,k;
	int jj;
	int totnum;
	int winlen;

	fprintf(fws,"rotate x, 180\n");
	fprintf(fws,"backbone only\n");
	fprintf(fws,"backbone 10\n");
	
	totnum=Multi_AFB[0][0][0];
	for(k=1;k<=totnum;k++)
	{
		for(j=0;j<TOT_NUM;j++)
		{
			//record
			jj=Multi_AFB[k][j][1];
			if(jj==-1)continue;
			//get chain
			char chain;
			{
				int pos=j;
				if(pos>=0 && pos<=25)chain=pos+'A';
				else if(pos>=26 && pos<=51)chain=pos-26+'a';
				else if(pos>=52 && pos<=61)chain=pos-52+'1';
				else chain='0';
			}
			winlen=Multi_AFB_Record[k];
			fprintf(fws,"select   %d-%d:%c\n", jj+1,jj+winlen,chain);
			//__070908__//judge whether is Real_Core or not
			if(Real_Block[k]<0)
			{
				fprintf(fws,"backbone off\n");
				fprintf(fws,"strands\n");
			}
			else
			{
				fprintf(fws,"backbone off\n");
				fprintf(fws,"cartoon\n");
			}
		}
	}
	fprintf(fws,"select all\n");
	fprintf(fws,"color chain\n");
}

//---------------- main ----------------//
//[input]:
// mol -> final superposed structures
// ami -> amino acied of those structures
// len -> length of those structures
// totnum -> total input structures number
// maxlen -> MSA(Multiple-Structure-Alignment) total length
// ali -> MSA itself (in 0,1,2 style, where '0':gap, '1':unaligned, '2': aligned)
void MultiAlign_Out::BC_Output_All(string &f3,XYZ **mol,char **ami,int *len,int totnum,int maxlen,int **ali,int TYPE)
{
	FILE *fp;
	//[rasmol]
	BC_Ali_To_AFB(ali,maxlen,totnum,Multi_AFB,Multi_AFB_Record,Real_Block);
	fp=fopen(f3.c_str(),"wb");
	if(fp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",f3.c_str());
	}
	else
	{
		if(TYPE==1) //output JMol
		{
			BC_Output_JMol_Script(fp,totnum,Multi_AFB,Multi_AFB_Record,Real_Block);
		}
		else        //output RasMol
		{
			BC_Output_RasMol_Script(fp,totnum,Multi_AFB,Multi_AFB_Record,Real_Block);
		}
		fclose(fp);
	}
}

