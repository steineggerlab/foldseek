#include "XYZ_DBT.h"

//--------------constructor--------------------//
XYZ_DBT::XYZ_DBT(void)
{
	//init_process (determine_first_three)
	cle_X_Axis[0]=1.0;
	cle_X_Axis[1]=0.0;
	cle_X_Axis[2]=0.0;
	cle_Z_Axis[0]=0.0;
	cle_Z_Axis[1]=0.0;
	cle_Z_Axis[2]=1.0;
}
XYZ_DBT::~XYZ_DBT(void)
{
}

//--------------- function ----------------//
void XYZ_DBT::xyz_to_dbt(double *dist,double *bend,double *tort,int n,XYZ *r)
{
	int i,j;
	double x0,x1,y0,y1;
	double b,t;	

	//------------------//
	//----real_start---//[0-2]
	//----------------//
	//calculate_bend
	(r[1]-r[0]).xyz2double(cle_x_point[0]);
	x0=sqrt(dot(cle_x_point[0],cle_x_point[0]));
	(r[2]-r[1]).xyz2double(cle_x_point[1]);	
	x1=sqrt(dot(cle_x_point[1],cle_x_point[1]));
	cross(cle_y_point[0],cle_x_point[0],cle_x_point[1]); 
	y0=sqrt(dot(cle_y_point[0],cle_y_point[0]));
	b=dot(cle_x_point[0],cle_x_point[1])/x0/x1;
	b=limit(b,-1.0,1.0);
	cle_w_point[0]=acos(1.0*b); 
	//evaluate
	if(dist!=NULL)dist[0]=-1.0;
	if(dist!=NULL)dist[1]=x0;
	if(dist!=NULL)dist[2]=x1;
	if(bend!=NULL)bend[0]=-9.9;
	if(bend!=NULL)bend[1]=-9.9;
	if(bend!=NULL)bend[2]=cle_w_point[0];
	if(tort!=NULL)tort[0]=-9.9;
	if(tort!=NULL)tort[1]=-9.9;
	if(tort!=NULL)tort[2]=-9.9;
	//update
	x0=x1;
	for(j=0;j<3;j++)cle_x_point[0][j]=cle_x_point[1][j];

	//--------------------//
	//----real_process---//[3-n]
	//------------------//
	for(i=3;i<n;i++)
	{
		//calculate_tort
		(r[i]-r[i-1]).xyz2double(cle_x_point[1]);
		x1=sqrt(dot(cle_x_point[1],cle_x_point[1]));
		cross(cle_y_point[1],cle_x_point[0],cle_x_point[1]); 
		y1=sqrt(dot(cle_y_point[1],cle_y_point[1]));
		t=dot(cle_y_point[0],cle_y_point[1])/y0/y1;
		t=limit(t,-1.0,1.0);
		t=acos(1.0*t);
		if(dot(cle_y_point[0],cle_x_point[1])<0.0) t*=-1.0;    
		cle_w_point[1]=t;
		//calculate_bend
		b=dot(cle_x_point[0],cle_x_point[1])/x0/x1;
		b=limit(b,-1.0,1.0);
		cle_w_point[2]=acos(1.0*b);
		//evaluate
		if(dist!=NULL)dist[i]=x1;
		if(bend!=NULL)bend[i]=cle_w_point[2];
		if(tort!=NULL)tort[i]=cle_w_point[1];
		//update
		x0=x1;
		y0=y1;
		for(j=0;j<3;j++)cle_x_point[0][j]=cle_x_point[1][j];
		for(j=0;j<3;j++)cle_y_point[0][j]=cle_y_point[1][j];
		cle_w_point[0]=cle_w_point[2];
	}
}
void XYZ_DBT::dbt_to_xyz(double *dist,double *bend,double *tort,int n,XYZ *r)
{
	int i;
	double radi;
	XYZ xyz;
	xyz=0.0;

	//real_process//
	Matrix_Unit(cle_T_Back,1.0);
	for(i=0;i<n;i++)
	{
		radi=dist[i];
		if(i==0) 
		{
			r[i]=0.0;
		} 
		else if(i==1) 
		{
			r[i]=0.0;
			r[i].X=radi;
		} 
		else if(i==2) 
		{
			Universal_Rotation(cle_Z_Axis,bend[i],cle_R_Theta);
			Vector_Dot(cle_D_Back,radi,cle_X_Axis,3);
			Matrix_Multiply(cle_T_Pre,cle_T_Back,cle_R_Theta);
			Vector_Multiply(cle_D_Pre,cle_T_Pre,cle_D_Back);
			Matrix_Equal(cle_T_Back,cle_T_Pre);
			xyz.double2xyz(cle_D_Pre);
			r[i]=r[i-1]+xyz;
		} 
		else 
		{
			Universal_Rotation(cle_Z_Axis,bend[i],cle_R_Theta);
			Universal_Rotation(cle_X_Axis,tort[i],cle_R_Thor);
			Vector_Dot(cle_D_Back,radi,cle_X_Axis,3);
			Matrix_Multiply(cle_temp,cle_R_Thor,cle_R_Theta);
			Matrix_Multiply(cle_T_Pre,cle_T_Back,cle_temp);
			Vector_Multiply(cle_D_Pre,cle_T_Pre,cle_D_Back);
			Matrix_Equal(cle_T_Back,cle_T_Pre);
			xyz.double2xyz(cle_D_Pre);
			r[i]=r[i-1]+xyz;
		}
	}
}
