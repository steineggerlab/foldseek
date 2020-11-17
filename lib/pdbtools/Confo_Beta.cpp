#include "Confo_Beta.h"

double CB_AMI_Dist[26]=
{ 1.5,-1.0,2.0,2.5,3.1,3.4,0.0,3.1,2.3,-1.0,3.5,2.6,3.0,2.5,-1.0,1.9,3.1,4.1,1.9,1.9,-1.0,2.0,3.9,-1.0,3.8,-1.0};
// A    B   C   D   E   F   G   H   I    J   K   L   M   N    O   P   Q   R   S   T    U   V   W    X   Y    Z
// 0    1   2   3   4   5   6   7   8    9  10  11  12  14   14  15  16  17  18  19   20  21  22   23  24   25
double CB_CLE_Type[17][2]=
{
	{0.9952,-2.0677},
	{1.0543,-2.3374},
	{1.0045,-2.4657},
	{1.0746,-2.5422},
	{1.1064,-2.5927},
	{1.1015,-2.5619},
	{1.0282,-2.5262},
	{0.9489,-1.9581},
	{0.9733,-2.0332},
	{0.9628,-2.0315},
	{1.0555,-2.1765},
	{1.1752,-2.5056},
	{0.9995,-2.1944},
	{1.1792,-2.2233},
	{1.2871,-2.2675},
	{1.0195,-2.1373},
	{1.1023,-2.2232}
};
// first->bend, second->tort

//--------------constructor--------------------//
Confo_Beta::Confo_Beta(int num)
{
	Confo_Beta_Maximal=num;
	//data_init
	Confo_Beta_Init(Confo_Beta_Maximal);
	//evaluate_init
	Confo_Beta_Levit_Angle=0.656;  //levit CA<->CB angle
	Confo_Beta_Levit_Radii=1.53;   //normal CA<->CB dist
}
Confo_Beta::~Confo_Beta(void)
{
	Confo_Beta_Dele();
}

//------------ init ------------//
void Confo_Beta::Confo_Beta_Init(int totlen)
{
	//temp_structure
	beta_dist=new double[totlen];
	beta_bend=new double[totlen];
	beta_tort=new double[totlen];
}
void Confo_Beta::Confo_Beta_Dele(void)
{
	//temp_structure
	delete [] beta_dist;
	delete [] beta_bend;
	delete [] beta_tort;
}

//--------------function------------//
//given CA+CB -> return dist,bend,tort
void Confo_Beta::Confo_Beta_CACB_To_Angle(XYZ *CA,XYZ *CB,int moln,double *bend,double *tort,double *dist)
{
	int i;
	double r,b,t;
	XYZ xyz;
	//process
	for(i=1;i<moln-1;i++)
	{
		//get_point
		xyz=(CA[i]-CA[i-1]);
		xyz.xyz2double(beta_v1);
		xyz=(CA[i+1]-CA[i]);
		xyz.xyz2double(beta_v2);
		xyz=(CB[i]-CA[i]);
		xyz.xyz2double(beta_y);
		//calc
		r=dot(beta_y,beta_y);
		if(r<1.e-3)
		{
			r=0.0;
			b=0.0;
			t=0.0;
		}
		else
		{
			//calc_bend
			b=Vector_Angle(beta_v1,beta_y,3);
			//calc_tort
			cross(beta_x1,beta_v1,beta_v2);
			cross(beta_x2,beta_v1,beta_y);
			t=Vector_Angle(beta_x1,beta_x2,3);
			if(dot(beta_x1,beta_y)<0.0)t*=-1.0;
		}
		//evaluate
		if(dist!=NULL)dist[i]=r;
		if(bend!=NULL)bend[i]=b;
		if(tort!=NULL)tort[i]=t;
	}
	//assign head_tail
	if(dist!=NULL)dist[0]=-1.0;
	if(bend!=NULL)bend[0]=0.0;
	if(tort!=NULL)tort[0]=0.0;
	if(dist!=NULL)dist[moln-1]=-1.0;
	if(bend!=NULL)bend[moln-1]=0.0;
	if(tort!=NULL)tort[moln-1]=0.0;
}
//given CA and dist,bend,tort -> return CB (Zheng & Liu's method)
void Confo_Beta::Confo_Beta_Angle_To_CACB(XYZ *CA,XYZ *CB,int moln,double *bend,double *tort,double *dist)
{
	int i,j;
	XYZ xyz;
	double beta,radii;
	//process
	for(i=1;i<moln-1;i++)
	{
		//init
		if(dist==NULL)beta=Confo_Beta_Levit_Radii;
		else beta=dist[i];
		if(beta<0.0)
		{
			CB[i]=CA[i];
			continue;
		}
		//calc
		xyz=(CA[i]-CA[i-1]);
		xyz.xyz2double(beta_v1);
		xyz=(CA[i+1]-CA[i]);
		xyz.xyz2double(beta_v2);
		cross(beta_y,beta_v1,beta_v2);
		Universal_Rotation(beta_y,bend[i],Confo_Beta_rotmat);
		Vector_Multiply(beta_x1,Confo_Beta_rotmat,beta_v1);
		Universal_Rotation(beta_v1,tort[i],Confo_Beta_rotmat);
		Vector_Multiply(beta_x2,Confo_Beta_rotmat,beta_x1);
		//evaluate
		radii=sqrt(dot(beta_x2,beta_x2));
		CA[i].xyz2double(beta_ca);
		for(j=0;j<3;j++)beta_cb[j]=beta_ca[j]+beta*beta_x2[j]/radii;
		CB[i].double2xyz(beta_cb);
	}
	//assign head_tail //this might be modified in the future...
	if(moln>1)
	{
		CB[0]=CB[1]-CA[1]+CA[0];
		CB[moln-1]=CB[moln-2]-CA[moln-2]+CA[moln-1];
	}
	else
	{
		CB[0]=CA[0];
		CB[moln-1]=CA[moln-1];
	}
}
//given CA -> return CB (Levitt's method)
void Confo_Beta::Confo_Beta_Levit_To_CACB(XYZ *CA,XYZ *CB,int moln,double *dist)
{
	int i,j;
	XYZ xyz;
	double beta,theta;

	//process
	theta=Confo_Beta_Levit_Angle;
	for(i=1;i<moln-1;i++) //omit head and tail
	{
		//init
		if(dist==NULL)beta=Confo_Beta_Levit_Radii;
		else beta=dist[i];
		if(beta<0.0)
		{
			CB[i]=CA[i];
			continue;
		}
		//calc
		xyz=(CA[i]-CA[i-1]);
		xyz.xyz2double(beta_v1);
		xyz=(CA[i]-CA[i+1]);
		xyz.xyz2double(beta_v2);
		Vector_Addition(beta_x1,beta_v1,beta_v2,3);
		Vector_Normalize(beta_x1,3);
		cross(beta_x2,beta_v1,beta_v2);
		Vector_Normalize(beta_x2,3);
		//evaluate
		CA[i].xyz2double(beta_ca);
		for(j=0;j<3;j++)beta_cb[j]=beta_ca[j]+beta*(cos(theta)*beta_x1[j]+sin(theta)*beta_x2[j]);
		CB[i].double2xyz(beta_cb);
	}
	//assign head_tail //this might be modified in the future...
	if(moln>1)
	{
		CB[0]=CB[1]-CA[1]+CA[0];
		CB[moln-1]=CB[moln-2]-CA[moln-2]+CA[moln-1];
	}
	else
	{
		CB[0]=CA[0];
		CB[moln-1]=CA[moln-1];
	}
}

//===================== Recon_Beta ==================//
//given CA+CB, and AMI_Dist, update SC
void Confo_Beta::Recon_Beta_1(XYZ *CA,XYZ *CB,int moln,char *ami)  //return SC
{
	int i;
	XYZ xyz;
	double radii;
	double multi;
	double v[3];
	for(i=0;i<moln;i++)
	{
		xyz=(CB[i]-CA[i]);
		xyz.xyz2double(v);
		radii=dot(v,v);
		if(radii<1.e-3)continue;
		Vector_Normalize(v,3);
		multi=CB_AMI_Dist[ami[i]-'A'];
		if(multi<0.0)continue;
		Vector_Dot(v,multi,v,3);
		xyz.double2xyz(v);
		CB[i]=CA[i]+xyz;
	}
}
//given CA, return CB (Levitt)
void Confo_Beta::Recon_Beta_31(XYZ *CA,XYZ *CB,int moln,char *ami)  //return CB
{
	int i;
	if(moln<=0)return;
	//init
	for(i=0;i<moln;i++)
	{
		//dist
		if(ami[i]=='G')beta_dist[i]=0.0;
		else beta_dist[i]=1.53;
	}
	//process
	Confo_Beta_Levit_To_CACB(CA,CB,moln,beta_dist);
}
void Confo_Beta::Recon_Beta_32(XYZ *CA,XYZ *CB,int moln,char *ami)  //return SC
{
	int i;
	double value;
	if(moln<=0)return;
	//init
	for(i=0;i<moln;i++)
	{
		//dist
		value=CB_AMI_Dist[ami[i]-'A'];
		if(value<0.0)beta_dist[i]=1.53;
		else beta_dist[i]=value;
	}
	//process
	Confo_Beta_Levit_To_CACB(CA,CB,moln,beta_dist);
}
//given CA, and CLE_Type, return CB (Liu & Wang)
void Confo_Beta::Recon_Beta_21(XYZ *CA,XYZ *CB,int moln,char *ami,char *cle)  //recon CB
{
	int i;
	char code;
	if(moln<=0)return;
	//init
	for(i=0;i<moln;i++)
	{
		//dist
		if(ami[i]=='G')beta_dist[i]=0.0;
		else beta_dist[i]=1.53;
		//bend+tort
		code=cle[i];
		if(code=='R')code='Q'; //use the most random one
		beta_bend[i]=CB_CLE_Type[code-'A'][0];
		beta_tort[i]=CB_CLE_Type[code-'A'][1];
	}	
	//process
	Confo_Beta_Angle_To_CACB(CA,CB,moln,beta_bend,beta_tort,beta_dist);
}
void Confo_Beta::Recon_Beta_22(XYZ *CA,XYZ *CB,int moln,char *ami,char *cle)  //recon SC
{
	int i;
	char code;
	double value;
	if(moln<=0)return;
	//init
	for(i=0;i<moln;i++)
	{
		//dist
		value=CB_AMI_Dist[ami[i]-'A'];
		if(value<0.0)beta_dist[i]=1.53;
		else beta_dist[i]=value;
		//bend+tort
		code=cle[i];
		if(code=='R')code='Q'; //use the most random one
		beta_bend[i]=CB_CLE_Type[code-'A'][0];
		beta_tort[i]=CB_CLE_Type[code-'A'][1];
	}	
	//process
	Confo_Beta_Angle_To_CACB(CA,CB,moln,beta_bend,beta_tort,beta_dist);
}
