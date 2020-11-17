#pragma once
#include <math.h>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

//----whether define DEBUG or not ------//
//#define DEBUG     //now is DEBUG situation
//-------------------------------------//

//=========== maximal protein length ========//
#ifndef PROT_MAX_NUM
#define PROT_MAX_NUM 5000
#endif
//===========================================//


//=========== ERROR_TAG ============//
#define STR_TRANS_ERROR -400  //this error occurs when input string cannot transform to [int] or [double]
#define FILE_LOAD_ERROR -999  //this error occurs when input file cannot open (file not found)
#define FILE_FORM_ERROR -888  //this error occurs when input file's format is not correct
#define RANGEOVER_ERROR -333  //this error occurs when threr are some range problem (exceed max,or < 0)
#define UPSIDDOWN_ERROR -987  //this error occurs when head value is larger than tail value
#define CHAINMISS_ERROR -432  //this error occurs when chain is not found
#define RESIDOVER_ERROR -234  //this error occurs when residues is over range
#define PREMATURE_ERROR -200  //this error occurs when input file ends at unexpected position
#define DATA_LENG_ERROR -250  //this error occurs when input string's length is not correct
#define USER_TYPE_ERROR -100  //this error occurs when user's input type is not valid
#define PARAMETER_ERROR -150  //this error occurs when user's parameter value is not valid
#define UNSUPPORT_ERROR -500  //this error occurs when user load a Virtual-Function at Base-Class
#define UNEXPECT_ERROR -777   //this error occurs for some unexpected situations
//==================================//


//==Definition==//
//--PI & NULL--//
#ifndef M_PI
#define M_PI          3.14159265358979323846
#endif
#ifndef NULL
#define NULL 0
#endif

//--integer_boundary--//
#ifndef INT_MAX_NUM
#define INT_MAX_NUM 100000
#endif
#ifndef INT_MIN_NUM
#define INT_MIN_NUM -99999
#endif

//---other_definition--//
#define     EQN_EPS     1.e-9
#define	    IsZero(x)	((x) > -EQN_EPS && (x) < EQN_EPS)
//==Definition==//over



//============== New+Delete+Equal Array ==============//
//new
template <class A>
inline void NewArray2D(A *** warray,int Narray1,int Narray2)
{
	*warray=new A*[Narray1];
	for(int i=0;i<Narray1;i++){
	    *(*warray+i)=new A[Narray2];
    }
}
template <class A>
inline void NewArray3D(A **** warray,int Narray1,int Narray2,int Narray3)
{
	*warray=new A**[Narray1];
	for(int i=0;i<Narray1;i++) *(*warray+i)=new A*[Narray2];
	for(int i=0;i<Narray1;i++){
	    for(int j=0;j<Narray2;j++) {
	        *(*(*warray+i)+j)=new A[Narray3];
        }
    }
}
template <class A>
inline void NewArray4D(A ***** warray,int Narray1,int Narray2,int Narray3,int Narray4)
{
	*warray=new A***[Narray1];
	for(int i=0;i<Narray1;i++) *(*warray+i)=new A**[Narray2];
	for(int i=0;i<Narray1;i++)for(int j=0;j<Narray2;j++) *(*(*warray+i)+j)=new A*[Narray3];
	for(int i=0;i<Narray1;i++)for(int j=0;j<Narray2;j++)for(int k=0;k<Narray3;k++) *(*(*(*warray+i)+j)+k)=new A[Narray4];
}
//delete
template <class A>
inline void DeleteArray2D(A *** warray,int Narray)
{
	if((*warray)==NULL)return;
	for(int i=0;i<Narray;i++) if(*(*warray+i)) delete [] *(*warray+i);
	if(Narray) delete [] (*warray);
	(*warray)=NULL;
}
template <class A>
inline void DeleteArray3D(A **** warray,int Narray1,int Narray2)
{
	if((*warray)==NULL)return;
	for(int i=0;i<Narray1;i++) for(int j=0;j<Narray2;j++)  if(*(*(*warray+i)+j)) delete [] *(*(*warray+i)+j);
	for(int i=0;i<Narray1;i++) if(*(*warray+i)) delete [] *(*warray+i);
	if(Narray1) delete [] (*warray);
	(*warray)=NULL;
}
template <class A>
inline void DeleteArray4D(A ***** warray, int Narray1,int Narray2,int Narray3 )
{
	if((*warray)==NULL)return;
	for(int i=0;i<Narray1;i++) for(int j=0;j<Narray2;j++) for(int k=0;k<Narray3;k++) if(*(*(*(*warray+i)+j)+k)) delete []*(*(*(*warray+i)+j)+k);
	for(int i=0;i<Narray1;i++) for(int j=0;j<Narray2;j++) if(*(*(*warray+i)+j)) delete [] *(*(*warray+i)+j);
	for(int i=0;i<Narray1;i++) if(*(*warray+i)) delete [] *(*warray+i);
	if(Narray1) delete [] (*warray);
	(*warray)=NULL;
}
//equal
template <class A> 
void EqualArray(A *out,A *in,int len)
{
	for(int i=0;i<len;i++) *(out+i)=*(in+i);
}
template <class A> 
void EqualArray2D(A **out,A **in,int len1,int len2)
{
	for(int i=0;i<len1;i++)for(int j=0;j<len2;j++) *(*(out+i)+j)=*(*(in+i)+j);
}
template <class A> 
void EqualArray3D(A ***out,A ***in,int len1,int len2,int len3)
{
	for(int i=0;i<len1;i++)for(int j=0;j<len2;j++)for(int k=0;k<len3;k++) *(*(*(out+i)+j)+k)=*(*(*(in+i)+j)+k);
}
template <class A> 
void EqualArray4D(A ****out,A ****in,int len1,int len2,int len3,int len4)
{
	for(int i=0;i<len1;i++)for(int j=0;j<len2;j++)for(int k=0;k<len3;k++)for(int l=0;l<len4;l++) *(*(*(*(out+i)+j)+k)+l)=*(*(*(*(in+i)+j)+k)+l);
}

//============== Maximal_Minimal ========//
template <class A> 
inline int Get_Max(A *input,int totnum,int start=0)
{
	int i;
	int maxnum=start;
	A dist,wsmax=input[start];
	for(i=1;i<totnum;i++)
	{
		dist=input[i+start];
		if(dist>wsmax)
		{
			wsmax=dist;
			maxnum=i+start;
		}
	}
	return maxnum;
}
template <class A>
inline int Get_Min(A *input,int totnum,int start=0)
{
	int i;
	int maxnum=start;
	A  dist,wsmax=input[start];
	for(i=1;i<totnum;i++)
	{
		dist=input[i+start];
		if(dist<wsmax)
		{
			wsmax=dist;
			maxnum=i+start;
		}
	}
	return maxnum;
}

//============ Utility =========//
template<class A>
inline void SWAP(A &a, A &b)
{
	A dum=a;
	a=b;
	b=dum;
}
template<class A>
inline const A SIGN(const A &a, const A &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
template<class A>
inline const A MAX(const A &a, const A &b)
{
	return b > a ? (b) : (a);
}
template<class A>
inline const A MIN(const A &a, const A &b)
{
	return b < a ? (b) : (a);
}



//============= System_Related ==========//
extern char DIR_CHAR;
extern void getline_end(string &input,char kill);
//============= Debug_Related ==========//
extern void overrange_debug_test(int cur,int limit);
extern void lessrange_debug_test(int cur,int limit);


//============== Math_Utility ===========//
extern double cbrt(double x);
extern double dot(double *x,double *y,int n);
extern double dot(double x[3],double y[3]);
extern void cross(double *z,double *x,double *y,int n);
extern void cross(double z[3],double x[3],double y[3]);
extern double limit(double x,double low,double high);
extern void Universal_Rotation(double Axis[3],double Phi,double R[3][3]);
extern void ND_Gaussian_Train(double **input,int num,int dim,double *mean,double **vari,double *wei=0);
//============== Matrix_Related =========//
extern void Matrix_Unit(double **a,int n1,int n2,double unit);
extern void Matrix_Unit(double a[3][3],double unit);
extern void Matrix_Equal(double **a,double **b,int n1,int n2);
extern void Matrix_Equal(double a[3][3],double b[3][3]);
extern void Matrix_Trans(double **a,double **b,int n);
extern void Matrix_Trans(double a[3][3],double b[3][3]);
//============== Matrix_Process =========//
extern double Matrix_Determinant(double **a,int n);                        //|A|
extern double Matrix_Determinant(double a[3][3]);                          //|3|
extern void Matrix_Multiply(double **z,double **x, double **y,int n);      //Z=X*Y
extern void Matrix_Multiply(double z[3][3],double x[3][3],double y[3][3]); //3=3*3
extern void Matrix_Addition(double **z,double **x, double **y,int n);      //Z=X+Y
extern void Matrix_Addition(double z[3][3],double x[3][3], double y[3][3]);//3=3+3
extern void Matrix_Subtract(double **z,double **x, double **y,int n);      //Z=X-Y
extern void Matrix_Subtract(double z[3][3],double x[3][3], double y[3][3]);//3=3-3
extern int Matrix_Algebraic(double **a,double **b,int n,int ii,int jj);    //A=B[ii,jj]
extern int Matrix_Cholesky(double **a,double **b,int n);                   //A*A^T=B
extern int Matrix_Inverse(double **a,double **b,int n);                    //A=A'
//============== Vector_Related =========//
extern void Vector_Unit(double *a,int n,double unit);
extern double Vector_Normalize(double *a,int n);
extern void Vector_Addition(double *z,double *x, double *y,int n);         //z=x+y
extern void Vector_Subtract(double *z,double *x, double *y,int n);         //z=x-y
extern void Vector_Multiply(double *z,double **x, double *y,int n);        //z=R.y
extern void Vector_Multiply(double z[3],double x[3][3], double y[3]);      //z=3.y
extern void Vector_Multiply_Trans(double *z,double **x, double *y,int n);  //z=y.R
extern void Vector_Multiply_Trans(double z[3],double x[3][3], double y[3]);//z=y.3
extern void Vector_Dot(double *z,double a,double *y,int n);                //z=ay
extern double Vector_Angle(double *a,double *b,int n);                     //a^b angle
extern double Dihedral_Angle(double *a,double *b,double *c,int n);         //a^b^c angle (dihedral angle)

//======= Integer/double Conversion =====//
extern string int2str(int num);
extern int str2dou(const char *str,double &num,int limit=15);
extern int str2int(const char *str,int &num);
//============== UpperLowerCase =========//
extern void toUpperCase(char *buffer);
extern void toLowerCase(char *buffer);
extern void toUpperCase(string &buffer);
extern void toLowerCase(string &buffer);
//============== Base and Root name ===//
extern void getBaseName(string &in,string &out,char slash,char dot);
extern void getRootName(string &in,string &out,char slash);
//======= String Process Utility ====//
extern int String_Process_Line_Dou(string &line,char sepa,double *out);
extern int String_Process_Line_Int(string &line,char sepa,int *out);
extern int String_Process_Line_Str(string &line,char sepa,char **out);
//======= Digit Process Utility ====//
extern int MI_Digi_Return(int *digit,int pos,int len);
extern void MI_Uni_For_Map(int &pos,int *input,int *digit,int len);  //givin *input, return pos
extern void MI_Uni_Bak_Map(int pos,int *output,int *digit,int len);  //given pos, return *output


//================= MingFu Utility ==============//
extern double log_add(double x, double y);
extern double log_add(double x1, double x2, double x3);
extern double log_add(double x1, double x2, double x3, double x4, double x5);
extern double MAX(double x1, double x2, double x3, double x4, double x5, int & b);
extern double MAX(double x1, double x2, double x3, int & b);
extern double MAX(double x1, double x2, int & b);
// x <- x * y
extern int Matrix_Add(vector< vector<double> > & x, const vector< vector<double > > & y);
extern int Matrix_Multiply(vector< vector<double> > & x, const vector< vector<double > > & y);
extern int Matrix_Average(vector< vector<double> > & x, const vector< vector<double > > & y, double weight);
