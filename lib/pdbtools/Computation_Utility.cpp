#include "Computation_Utility.h"
#include <cassert>
#include <vector>

using namespace std;

//==================================utility==========================//
//---------------------//
//-- System_Related --//
//-------------------//
char DIR_CHAR=0x2F;

//-------------------//
//-- Mathematical --//
//-----------------//
//--------cubic_root-------------//
double cbrt(double x)
{
	if(IsZero(x))return 0.0;
	if(x>0.0)return exp(1.0*log(x)/3);
	if(x<0.0)return -1.0*exp(1.0*log(fabs(x))/3);
	return 0.0;
}
//---------dot_multiply-----------//
double dot(double *x,double *y,int n)
{ 
	int i;  
	double t=0.;  
	for(i=0;i<n;i++)t+=x[i]*y[i]; 
	return(t);
}
double dot(double x[3],double y[3])
{ 
	int i;  
	int n=3;
	double t=0.;  
	for(i=0;i<n;i++)t+=x[i]*y[i]; 
	return(t);
}
//---------cross_multiply--------//
void cross(double *z,double *x,double *y,int n)
{
	int i;
	int cur0,cur1,cur2;
	for(i=0;i<n;i++)
	{
		cur0=(i+0)%n;
		cur1=(i+1)%n;
		cur2=(i+2)%n;
		z[cur0]=x[cur1]*y[cur2]-x[cur2]*y[cur1];
	}
}
void cross(double z[3],double x[3],double y[3])
{
	z[0]=x[1]*y[2]-x[2]*y[1];
	z[1]=x[2]*y[0]-x[0]*y[2];
	z[2]=x[0]*y[1]-x[1]*y[0];
}
//---------check_limit-----------//
double limit(double x,double low,double high)
{
	if(x>=high)return high;
	if(x<=low)return low;
	return x;
}

//----------- Universal_Rotation -----------//
// Set Vector->Axis(x,y,z) , Rotation->Phi(-pi,pi)
void Universal_Rotation(double Axis[3],double Phi,double R[3][3])
{
	double Axis_Value;
	double Lamda,Miu,Viu;

	Axis_Value = sqrt(dot(Axis,Axis,3));
	Lamda = Axis[0] / Axis_Value;  // X rotation_part
	Miu = Axis[1] / Axis_Value;    // Y rotation_part
	Viu = Axis[2] / Axis_Value;    // Z rotation_part

	R[0][0] = cos(Phi) + (1-cos(Phi)) * Lamda * Lamda  ;
	R[0][1] = (1-cos(Phi)) * Lamda * Miu - Viu * sin(Phi);
	R[0][2] = (1-cos(Phi)) * Lamda * Viu + Miu * sin(Phi);

	R[1][0] = (1-cos(Phi)) * Lamda * Miu + Viu * sin(Phi);
	R[1][1] = cos(Phi) + (1-cos(Phi)) * Miu * Miu ;
	R[1][2] = (1-cos(Phi)) * Viu * Miu - Lamda * sin(Phi);

	R[2][0] = (1-cos(Phi)) * Lamda * Viu - Miu * sin(Phi);
	R[2][1] = (1-cos(Phi)) * Viu * Miu + Lamda * sin(Phi);
	R[2][2] = cos(Phi) + (1-cos(Phi)) * Viu * Viu ;
}

//----------- ND_Gaussian_Train -----------//
// input **data, total num, and dimension
// output *mean, **vari
// [option] *wei->weight
void ND_Gaussian_Train(double **input,int num,int dim,double *mean,double **vari,double *wei)
{
	int i,j,k;
	double weight;
	double total;
	//mean_process
	total=0.0;
	for(i=0;i<dim;i++)mean[i]=0.0;
	for(k=0;k<num;k++)
	{
		if(wei==0)weight=1.0;
		else weight=wei[k];
		for(i=0;i<dim;i++)mean[i]+=weight*input[k][i];
		total+=weight;
	}
	for(i=0;i<dim;i++)mean[i]/=total;
	//vari_process
	total=0.0;
	for(i=0;i<dim;i++)for(j=0;j<dim;j++)vari[i][j]=0.0;
	for(k=0;k<num;k++)
	{
		if(wei==0)weight=1.0;
		else weight=wei[k];
		for(i=0;i<dim;i++)for(j=0;j<dim;j++)vari[i][j]+=weight*(input[k][i]-mean[i])*(input[k][j]-mean[j]);
		total+=weight;
	}
	for(i=0;i<dim;i++)for(j=0;j<dim;j++)vari[i][j]/=total;
}

//--------matrix_related--------//
void Matrix_Unit(double **a,int n1,int n2,double unit)
{
	int i,j;
	for(i=0;i<n1;i++)for(j=0;j<n2;j++)a[i][j]=0.0;
	for(i=0;i<n1;i++)a[i][i]=unit;
}
void Matrix_Unit(double a[3][3],double unit)
{
	int i,j;
	int n1=3;
	int n2=3;
	for(i=0;i<n1;i++)for(j=0;j<n2;j++)a[i][j]=0.0;
	for(i=0;i<n1;i++)a[i][i]=unit;
}
void Matrix_Equal(double **a,double **b,int n1,int n2)
{
	int i,j;
	for(i=0;i<n1;i++)for(j=0;j<n2;j++)a[i][j]=b[i][j];
}
void Matrix_Equal(double a[3][3],double b[3][3])
{
	int i,j;
	int n1=3;
	int n2=3;
	for(i=0;i<n1;i++)for(j=0;j<n2;j++)a[i][j]=b[i][j];
}
void Matrix_Trans(double **a,double **b,int n)
{
	int i,j;
	for(i=0;i<n;i++)for(j=0;j<n;j++)a[i][j]=b[j][i];
}
void Matrix_Trans(double a[3][3],double b[3][3])
{
	int i,j;
	int n=3;
	for(i=0;i<n;i++)for(j=0;j<n;j++)a[i][j]=b[j][i];
}

//---------matrix_determinant-------//
//input matrix:a,dim:n
//output determinant
double Matrix_Determinant(double **a,int n)
{
	int i,j,l;
	double **b;
	NewArray2D(&b,n,n);
	double result=1;
	double k;
	for(i=0;i<n;i++)for(j=0;j<n;j++)b[i][j]=a[i][j];
	for(i=0;i<n-1;i++)
	{
		for(j=i+1;j<n;j++)
		{
			k=b[j][i]/b[i][i];
			b[j][i]=0;
			for(l=i+1;l<n;l++)b[j][l]=b[j][l]-k*b[i][l];
		}
	}
	for(i=0;i<n;i++)result*=b[i][i];
	DeleteArray2D(&b,n);
	return result;
}
double Matrix_Determinant(double a[3][3])
{
	double a0,a1,a2,a3,a4,a5,a6,a7,a8;
	a0=a[0][0];
	a1=a[0][1];
	a2=a[0][2];
	a3=a[1][0];
	a4=a[1][1];
	a5=a[1][2];
	a6=a[2][0];
	a7=a[2][1];
	a8=a[2][2];
	double value;
	value=a0*(a4*a8-a7*a5)-a3*(a1*a8-a7*a2)+a6*(a1*a5-a4*a2);
	return value;
}

//---------matrix_multiply-------//
//input matrix:x and y,dim:n
//output matrix:z=x*y
void Matrix_Multiply(double **z,double **x,double **y,int n)
{
	int i,j,k;
	double sum;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			sum=0.0;
			for(k=0;k<n;k++)sum+=x[i][k]*y[k][j];
			z[i][j]=sum;
		}
	}
}
void Matrix_Multiply(double z[3][3],double x[3][3],double y[3][3])
{
	int i,j,k;
	int n=3;
	double sum;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			sum=0.0;
			for(k=0;k<n;k++)sum+=x[i][k]*y[k][j];
			z[i][j]=sum;
		}
	}
}
void Matrix_Addition(double **z,double **x, double **y,int n)
{
	int i,j;
	for(i=0;i<n;i++)for(j=0;j<n;j++)z[i][j]=x[i][j]+y[i][j];
}
void Matrix_Addition(double z[3][3],double x[3][3], double y[3][3])
{
	int i,j;
	int n=3;
	for(i=0;i<n;i++)for(j=0;j<n;j++)z[i][j]=x[i][j]+y[i][j];
}
void Matrix_Subtract(double **z,double **x, double **y,int n)
{
	int i,j;
	for(i=0;i<n;i++)for(j=0;j<n;j++)z[i][j]=x[i][j]-y[i][j];
}
void Matrix_Subtract(double z[3][3],double x[3][3], double y[3][3])
{
	int i,j;
	int n=3;
	for(i=0;i<n;i++)for(j=0;j<n;j++)z[i][j]=x[i][j]-y[i][j];
}

//----- algebraic_complement ----//__091220__//
//input matrix:b,dim:n
//input (ii,jj)
//output matrix:a,
//return sign
int Matrix_Algebraic(double **a,double **b,int n,int ii,int jj)
{
	int i,j;
	int x,y;
	int sign;
	x=0;
	for(i=0;i<n;i++)
	{
		if(i==ii)continue;
		y=0;
		for(j=0;j<n;j++)
		{
			if(j==jj)continue;
			a[x][y]=b[i][j];
			y++;
		}
		x++;
	}
	if((ii+jj)%2==0)sign=1;
	else sign=-1;
	return sign;
}

//---------matrix_cholesky -------//__091210__//
//matrix B is positive-definite symmetry, then B=A*A^T
//where A is a lower-triangular matrix
int Matrix_Cholesky(double **a,double **b,int n)
{
	int i,j,k;
	double sum;
	for(i=0;i<n;i++)for(j=0;j<n;j++)a[i][j]=b[i][j];
	for(i=0;i<n;i++)
	{
		for(j=i;j<n;j++)
		{
			for(sum=a[i][j],k=i-1;k>=0;k--)sum-=a[i][k]*a[j][k];
			if(i==j)
			{
				if(sum<=0.0)return -1;
				a[i][i]=sqrt(sum);
			}
			else a[j][i]=sum/a[i][i];
		}
	}
	for(i=0;i<n;i++)for(j=0;j<i;j++)a[j][i]=0.0;
	return 1;
}

//---------matrix_inverse ----------//__091220__//
//The input matrix is a[0..n-1][0..n-1]. b[0..n-1][0..m-1] is 
//input containing the m right-hand side vectors.
//On output, a is replaced by its matrix inverse, 
//and b is replaced by the corresponding set of
//solution vectors.
int Matrix_Inverse(double **a,double **b,int n)
{
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv;
	int *indxc=new int[n];
	int *indxr=new int[n];
	int *ipiv=new int[n];
	//process
	icol=0;
	irow=0;
	for(j=0;j<n;j++)ipiv[j]=0;
	for(i=0;i<n;i++)
	{
		big=0.0;
		for(j=0;j<n;j++)
			if(ipiv[j]!=1)
				for(k=0;k<n;k++)
				{
					if(ipiv[k]==0)
					{
						if(fabs(a[j][k])>=big)
						{
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if(irow!=icol)
		{
			for(l=0;l<n;l++)SWAP(a[irow][l],a[icol][l]);
			for(l=0;l<n;l++)SWAP(b[irow][l],b[icol][l]);
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if(a[icol][icol]==0.0)
		{
			if(indxc!=0)delete [] indxc;
			if(indxr!=0)delete [] indxr;
			if(ipiv!=0)delete [] ipiv;
			return -1;
		}
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for(l=0;l<n;l++)a[icol][l]*=pivinv;
		for(l=0;l<n;l++)b[icol][l]*=pivinv;
		for(ll=0;ll<n;ll++)
			if(ll!=icol)
			{
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for(l=0;l<n;l++)a[ll][l]-=a[icol][l]*dum;
				for(l=0;l<n;l++)b[ll][l]-=b[icol][l]*dum;
			}
	}
	for(l=n-1;l>=0;l--)
	{
		if(indxr[l]!=indxc[l])
			for(k=0;k<n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	//terminal
	if(indxc!=0)delete [] indxc;
	if(indxr!=0)delete [] indxr;
	if(ipiv!=0)delete [] ipiv;
	return 1;
}


//---------vector_related-------//
void Vector_Unit(double *a,int n,double unit)
{
	int i;
	for(i=0;i<n;i++)a[i]=unit;
}
double Vector_Normalize(double *a,int n)
{
	int i;
	double norm=sqrt(dot(a,a,n));
	if(!IsZero(norm))
	{
		for(i=0;i<n;i++)a[i]/=norm;
		return norm;
	}
	else
	{
		Vector_Unit(a,n,0.0);
		return 0.0;
	}
}
void Vector_Addition(double *z,double *x, double *y,int n)
{
	int i;
	for(i=0;i<n;i++)z[i]=x[i]+y[i];
}
void Vector_Subtract(double *z,double *x, double *y,int n)
{
	int i;
	for(i=0;i<n;i++)z[i]=x[i]-y[i];
}
void Vector_Multiply(double *z,double **x, double *y,int n)
{
	int i,k;
	double sum;
	for(i=0;i<n;i++)
	{
		sum=0.0;
		for(k=0;k<n;k++)sum+=x[i][k]*y[k];
		z[i]=sum;
	}
}
void Vector_Multiply(double z[3],double x[3][3], double y[3])
{
	int i,k;
	int n=3;
	double sum;
	for(i=0;i<n;i++)
	{
		sum=0.0;
		for(k=0;k<n;k++)sum+=x[i][k]*y[k];
		z[i]=sum;
	}
}
void Vector_Multiply_Trans(double *z,double **x, double *y,int n)
{
	int i,k;
	double sum;
	for(i=0;i<n;i++)
	{
		sum=0.0;
		for(k=0;k<n;k++)sum+=x[k][i]*y[k];
		z[i]=sum;
	}
}
void Vector_Multiply_Trans(double z[3],double x[3][3], double y[3])
{
	int i,k;
	int n=3;
	double sum;
	for(i=0;i<n;i++)
	{
		sum=0.0;
		for(k=0;k<n;k++)sum+=x[k][i]*y[k];
		z[i]=sum;
	}
}
void Vector_Dot(double *z,double a,double *y,int n)
{
	int i;
	for(i=0;i<n;i++)z[i]=a*y[i];
}
double Vector_Angle(double *a,double *b,int n)
{
	double aa=sqrt(dot(a,a,n));
	double bb=sqrt(dot(b,b,n));
	if(IsZero(aa)||IsZero(bb))return 0.0;
	double z=dot(a,b,n)/aa/bb;
	z=limit(z,-1.0,1.0);
	return acos(1.0*z);
}
double Dihedral_Angle(double *a,double *b,double *c,int n)
{
	double *v1=new double[n];
	double *v2=new double[n];
	cross(v1,a,b,n);
	cross(v2,b,c,n);
	double angle=Vector_Angle(v1,v2,n);
	if(dot(v1,c,3)<0.0)angle*=-1.0;
	delete [] v1;
	delete [] v2;
	return angle;
}


//----------------------//
//-- Getline_Ending ---//
//--------------------//
void getline_end(string &input,char kill)
{
	int len=(int)input.length();
	if(input[len-1]==kill)input=input.substr(0,len-1);
}
//-----------------//
//-- Debugging ---//
//---------------//
void overrange_debug_test(int cur,int limit)
{
	if(cur>=limit)
	{
		fprintf(stderr,"EXCEED LIMIT=>BREAK!!\n");
		exit(-1);
	}
}
void lessrange_debug_test(int cur,int limit)
{
	if(cur<limit)
	{
		fprintf(stderr,"UNDER LIMIT=>BREAK!!\n");
		exit(-1);
	}
}

//=================== integer/double conversion ====================//
//-------int_2_string------------//
string int2str(int num)
{
	if(num==0)return "0";
	string str="";
	int num_ = num>0?num:-1*num;
	while(num_)
	{
		str=(char)(num_%10+48)+str;
		num_ /=10;
	}
	if (num<0)str= "-" + str;
	return str;
}
//--------string_2_int-----------//
int str2int(const char *str,int &num)
{
	int i,j,len=(int)strlen(str);
	for(i=0;i<len;i++)if(str[i]!=' ')break;
	for(j=len-1;j>=i;j--)if(str[j]!=' ')break;
	if(str[i]=='-')i++;
	for(;i<=j;i++)if(str[i]<'0'||str[i]>'9')return 0; //error
	num=atoi(str);
	return 1; //correct
}
//--------string_2_double-----------//
int str2dou(const char *str,double &num,int limit)
{
	int i,j,len=(int)strlen(str),count=0;
	for(i=0;i<len;i++)if(str[i]!=' ')break;
	for(j=len-1;j>=i;j--)if(str[j]!=' ')break;
	if(str[i]=='-')i++;
	for(;i<=j;i++)
	{
		if(str[i]=='.')break;
		if(str[i]<'0'||str[i]>'9')return 0; //error
	}
	for(i++;i<=j;i++)
	{
		if(count>=limit)break;
		if(str[i]<'0'||str[i]>'9')return 0; //error
		count++;
	}
	num=atof(str);
	return 1; //correct
}

//=================== upper and lower case ====================//
//----------upper_case-----------//
void toUpperCase(char *buffer) 
{  
	for(int i=0;i<(int)strlen(buffer);i++) 
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
void toUpperCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++) 
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
//----------lower_case-----------//
void toLowerCase(char *buffer)
{  
	for(int i=0;i<(int)strlen(buffer);i++) 
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}
void toLowerCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++) 
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}
//--------- base_name -----------//__110830__//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}

//======================= String Process Utility ==================//
int String_Process_Line_Dou(string &line,char sepa,double *out)
{
	int i;
	string temp;
	int first=1;
	int ini=0;
	int wsc=0;
	int len=(int)line.length();
	int count=0;
	//ori
	for(i=0;i<len;i++)
	{
		if(line[i]==sepa)
		{
			if(first==0)
			{
				temp=line.substr(ini,wsc);
				out[count]=atof(temp.c_str());
				count++;
				first=1;
			}
		}
		else
		{
			if(first==1)
			{
				first=0;
				ini=i;
				wsc=0;
			}
			wsc++;
		}
	}
	//end
	if(first==0)
	{
		temp=line.substr(ini,wsc);
		out[count]=atof(temp.c_str());
		count++;
	}
	return count;
}
int String_Process_Line_Int(string &line,char sepa,int *out)
{
	int i;
	string temp;
	int first=1;
	int ini=0;
	int wsc=0;
	int len=(int)line.length();
	int count=0;
	//ori
	for(i=0;i<len;i++)
	{
		if(line[i]==sepa)
		{
			if(first==0)
			{
				temp=line.substr(ini,wsc);
				out[count]=atoi(temp.c_str());
				count++;
				first=1;
			}
		}
		else
		{
			if(first==1)
			{
				first=0;
				ini=i;
				wsc=0;
			}
			wsc++;
		}
	}
	//end
	if(first==0)
	{
		temp=line.substr(ini,wsc);
		out[count]=atoi(temp.c_str());
		count++;
	}
	return count;
}
int String_Process_Line_Str(string &line,char sepa,char **out)
{
	int i;
	string temp;
	int first=1;
	int ini=0;
	int wsc=0;
	int len=(int)line.length();
	int count=0;
	//ori
	for(i=0;i<len;i++)
	{
		if(line[i]==sepa)
		{
			if(first==0)
			{
				temp=line.substr(ini,wsc);
				strcpy(out[count],temp.c_str());
				count++;
				first=1;
			}
		}
		else
		{
			if(first==1)
			{
				first=0;
				ini=i;
				wsc=0;
			}
			wsc++;
		}
	}
	//end
	if(first==0)
	{
		temp=line.substr(ini,wsc);
		strcpy(out[count],temp.c_str());
		count++;
	}
	return count;
}

//==================== Digit Process Utility ====================//
int MI_Digi_Return(int *digit,int pos,int len)
{
	int i;
	int val=1;
	for(i=len-1;i>pos;i--)
	{
		val*=digit[i];
	}
	return val;
}
void MI_Uni_For_Map(int &pos,int *input,int *digit,int len)  //givin *input, return pos
{
	int i;
	int ival=0;
	for(i=0;i<len;i++)
	{
		ival+=input[i]*MI_Digi_Return(digit,i,len);
	}
	pos=ival;
}
void MI_Uni_Bak_Map(int pos,int *output,int *digit,int len) //given pos, return *output
{
	int i;
	int val;
	int temp=pos;
	for(i=0;i<len;i++)
	{
		val=MI_Digi_Return(digit,i,len);
		output[i]=temp/val;
		temp=temp%val;
	}
}


//================= MingFu Utility ==============//
// log(e^x + e^y) = x + log(1 + e^(b-a)), x is the bigger one
double log_add(double x, double y)
{
	if(x < y)
	{
		double t = y;
		y = x;
		x = t;
	}

	return x + log(1.0 + exp(y - x));
}

double log_add(double x1, double x2, double x3, double x4, double x5)
{
	if(x1 >= x2 && x1 >= x3 && x1 >= x4 && x1 >= x5) return x1 + log(1 + exp(x2 - x1) + exp(x3 - x1) + exp(x4 - x1) + exp(x5 - x1));
	if(x2 >= x1 && x2 >= x3 && x2 >= x4 && x2 >= x5) return x2 + log(1 + exp(x1 - x2) + exp(x3 - x2) + exp(x4 - x2) + exp(x5 - x2));
	if(x3 >= x2 && x3 >= x1 && x3 >= x4 && x3 >= x5) return x3 + log(1 + exp(x2 - x3) + exp(x1 - x3) + exp(x4 - x3) + exp(x5 - x3));
	if(x4 >= x2 && x4 >= x3 && x4 >= x1 && x4 >= x5) return x4 + log(1 + exp(x2 - x4) + exp(x3 - x4) + exp(x1 - x4) + exp(x5 - x4));
	if(x5 >= x2 && x5 >= x3 && x5 >= x4 && x5 >= x1) return x5 + log(1 + exp(x2 - x5) + exp(x3 - x5) + exp(x4 - x5) + exp(x1 - x5));
	return 0;
}

double log_add(double x1, double x2, double x3)
{
	if(x1 >= x2 && x1 >= x3) return x1 + log(1 + exp(x2 - x1) + exp(x3 - x1));
	if(x2 >= x1 && x2 >= x3) return x2 + log(1 + exp(x1 - x2) + exp(x3 - x2));
	if(x3 >= x2 && x3 >= x1) return x3 + log(1 + exp(x2 - x3) + exp(x1 - x3));
	return 0;
}

double MAX(double x1, double x2, double x3, double x4, double x5, int & b)
{
	if(x1 >= x2 && x1 >= x3 && x1 >= x4 && x1 >= x5) { b = 0; return x1;}
	if(x2 >= x1 && x2 >= x3 && x2 >= x4 && x2 >= x5) { b = 1; return x2;}
	if(x3 >= x2 && x3 >= x1 && x3 >= x4 && x3 >= x5) { b = 2; return x3;}
	if(x4 >= x2 && x4 >= x3 && x4 >= x1 && x4 >= x5) { b = 3; return x4;}
	if(x5 >= x2 && x5 >= x3 && x5 >= x4 && x5 >= x1) { b = 4; return x5;}
	return 0;
}

double MAX(double x1, double x2, double x3, int & b)
{
	if(x1 >= x2 && x1 >= x3) { b = 0; return x1;}
	if(x2 >= x1 && x2 >= x3) { b = 1; return x2;}
	if(x3 >= x2 && x3 >= x1) { b = 2; return x3;}
	return 0;
}

double MAX(double x1, double x2, int & b)
{
	if(x1 >= x2) { b = 0; return x1;}
	if(x2 >= x1) { b = 1; return x2;}
	return 0;
}

int Matrix_Multiply(vector< vector<double> > & x, const vector< vector<double> > & y)
{
	assert(x.size() == y.size());
	for(int i = 0; i < (int)x.size(); i++)
	{
		assert(x.at(i).size()==y.at(i).size());
		for(int j = 0; j < (int)x.at(i).size(); j++)
		{
			x.at(i).at(j) = x.at(i).at(j) * y.at(i).at(j);
		}
	}
	return 0;
}

int Matrix_Average(vector< vector<double> > & x, const vector< vector<double> > & y, double weight)
{
	assert(x.size() == y.size());
	for(int i = 0; i < (int)x.size(); i++)
	{
		assert(x.at(i).size()==y.at(i).size());
		for(int j = 0; j < (int)x.at(i).size(); j++)
		{
			x.at(i).at(j) = (1 - weight) * x.at(i).at(j) + weight * y.at(i).at(j);
		}
	}
	return 0;
}

int Matrix_Add(vector< vector<double> > & x, const vector< vector<double> > & y)
{
	assert(x.size() == y.size());
	for(int i = 0; i < (int)x.size(); i++)
	{
		assert(x.at(i).size()==y.at(i).size());
		for(int j = 0; j < (int)x.at(i).size(); j++)
		{
			x.at(i).at(j) = x.at(i).at(j) + y.at(i).at(j);
		}
	}
	return 0;
}

