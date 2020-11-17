#include "Confo_Lett.h"


//--------------constructor--------------------//
Confo_Lett::Confo_Lett(void)
{
	//init_process (determine_first_three)
	cle_X_Axis[0]=1.0;
	cle_X_Axis[1]=0.0;
	cle_X_Axis[2]=0.0;
	cle_Z_Axis[0]=0.0;
	cle_Z_Axis[1]=0.0;
	cle_Z_Axis[2]=1.0;
	//macros
	CLE_CUT_MIN=2.5;  //default
	CLE_CUT_MAX=4.5;  //default
}
Confo_Lett::~Confo_Lett(void)
{
}

//================================================PART_I====================================// STC_Related
//pi:mode weight; ni:mode counts; sig:|\sqrt(2\pi)\sigma|^-1; mid:mean; var:variance

//the main data of the 17 conformational letter
int Ori_id[Confo_Cluster]=
{'I','J','H','K','F','E','C','D','A','B','G','L','M','N','O','P','Q'};
int Ant_id[Confo_Cluster]=
{8,9,6,7,5,4,10,2,0,1,3,11,12,13,14,15,16};

double Ori_pi[Confo_Cluster]=
{0.0819,0.0727,0.1617,0.0594,0.0491,0.1158,0.0749,0.0544,0.0433,0.0387,0.0556,0.0533,0.0371,0.0315,0.0214,0.0318,0.0173};

double Ori_sig[Confo_Cluster]=
{1881.06,1796.72,10424.63,254.11,104.55,108.95,99.85,77.74,202.76,66.29,132.85,40.02,144.28,73.52,247.16,206.13,24.82};

double Ori_mid[Confo_Cluster][Confo_TriPep]=
{
	{1.517 ,  0.833  , 1.517},
	{1.577 ,  1.045  , 1.550},
	{1.550 ,  0.878  , 1.552},
	{1.479 ,  0.699  , 1.426},
	{1.094 , -2.722  , 0.909},
	{1.021 , -2.977  , 0.954},
	{1.014 , -1.878  , 1.138},
	{0.786 , -2.303  , 1.028},
	{1.025 , -2.003  , 1.551},
	{1.057 , -2.940  , 1.344},
	{1.491 ,  2.087  , 1.046},
	{1.396 ,  0.754  , 0.840},
	{1.466 ,  1.643  , 1.440},
	{1.121 ,  0.142  , 1.490},
	{1.537 , -1.891  , 1.484},
	{1.235 , -2.982  , 1.493},
	{0.863 , -0.372  , 1.008}
};
double Ori_var[Confo_Cluster][Confo_TriPep][Confo_TriPep]=
{
	{{275.430 , -28.322 , 106.938},
	 {-28.322 , 84.345  ,-46.140 },
	 {106.938 ,-46.140  ,214.445 }},

	{{314.298 , -10.252 , 37.793 },
	 {-10.252 , 46.021  , -69.977},
	 {37.793  ,-69.977  ,332.804 }},

	{{706.634 , -93.940 , 128.933},
	 {-93.940 , 245.505 ,-171.758},
	 {128.933 , -171.758,786.116 }},

	{{73.820  , -13.740 , 15.478 },
	 {-13.740 ,  21.548 , -25.281},
	 {15.478  ,-25.281  , 75.727 }},

	{{24.130  , 1.909   , -11.179},
 	 {1.909   , 10.921  ,  -8.787},
	 {-11.179 , -8.787  ,  53.040}},

	{{34.340  , 4.171   , -9.288 },
	 {4.171   , 15.243  , -22.456},
	 {-9.288  , -22.456 ,  56.837}},

	{{28.008  , 4.114   , 2.318  },
	 {4.114   ,  6.183  , -5.128 },
	 {2.318   , -5.128  , 69.368 }},

	{{56.205  , 3.764   , -10.830},
	 {3.764   ,  4.168  , -2.148 },
	 {-10.830 ,  -2.148 , 30.057 }},

	{{30.524  , 9.083   ,  6.041 },
	 {9.083   , 8.692   , 5.674  },
	 {6.041   , 5.674   ,228.584 }},

	{{26.872  , 4.632   , 9.545  },
	 {4.632   , 4.865   ,  -4.959},
	 {9.545   , -4.959  ,  54.317}},

	{{163.926 , 0.615   , 1.975  },
	 {0.615   , 3.777   , -3.748 },
	 {1.975   ,-3.748   ,32.287  }},

	{{43.653  ,  2.493  ,  -6.966},
	 {2.493   ,  1.425  ,  -2.864},
	 {-6.966  , -2.864  , 34.454 }},

	{{72.892  , 2.130   , 1.870  },
	 {2.130   , 4.842   ,-7.864  },
	 {1.870   , -7.864  , 72.914 }},

	{{25.283  , 3.210   , 9.908  },
	 {3.210   , 3.113   , 0.872  },
	 {9.908   , 0.872   , 82.959 }},

	{{170.760 , -0.715  , -4.138 },
	 {-0.715  , 3.726   , 3.069  },
	 {-4.138  , 3.069   , 98.705 }},

	{{47.969  , 8.244   , -4.912 },
	 {8.244   , 7.340   , -6.570 },
	 {-4.912  , -6.570  ,155.581 }},

	{{28.408  ,  1.513  , 3.380  },
	 {1.513   , 1.214   , 0.059  },
	 {3.380   , 0.059   , 19.549 }}
};

//===================================== main function body =============================//
//-----------------------//---------------------------------------------------//
//----XYZ->Confo_Lett---//
//---------------------//
double Confo_Lett::gaussian(int ic,double *xx)
{ 
	int j,k;  
	double y[Confo_TriPep],z;  
	for(j=0;j<Confo_TriPep;j++) 
	{
		y[j]=xx[j]-Ori_mid[ic][j];    
		if(j==1||j==3) 
		{    	
			if(y[j]>M_PI) y[j]-=2.0*M_PI; 
			else if(y[j]<-1.0*M_PI) y[j]+=2.0*M_PI;
		}
	}  
	for(z=0.0,j=0;j<Confo_TriPep;j++) for(k=0;k<=j;k++) z+=(j-k?1.0:0.5)*y[j]*Ori_var[ic][j][k]*y[k];
	return(exp(-z)*Ori_sig[ic]);
}

char Confo_Lett::btb_stc(double *w_point)
{
	int k;
	int ip=Confo_Cluster-1;
	double q,qq=-1.0;
	for(k=0;k<Confo_Cluster;k++) 
	{
		if((q=Ori_pi[k]*gaussian(k,w_point))>qq) 
		{
			qq=q;
			ip=k;
		}
	}
	return (Ori_id[ip]);
}
void Confo_Lett::profile_stc(double *w_point,double *output)
{
	int k;
	double q;
	double value=0.0;
	for(k=0;k<Confo_Cluster;k++) 
	{
		q=Ori_pi[k]*gaussian(k,w_point);
		output[Ori_id[k]-'A']=q;
		value+=q;
	}
	for(k=0;k<Confo_Cluster;k++)output[k]/=value;
}

//input XYZ_Type data_structure (*r), output Confo_Angle (*dist,*bend,*tort) and Confo_Lett (*CLE)
void Confo_Lett::btb_ori(double *dist,double *bend,double *tort,int n,XYZ *r,char *CLE,double **Profile)
{
	int i,j,k;
	double x0,x1,y0,y1;
	double b,t;	

	//init check
	if(n<3)
	{
		if(dist!=NULL)for(i=0;i<n;i++)dist[i]=-1.0;
		if(bend!=NULL)for(i=0;i<n;i++)bend[i]=-9.9;
		if(tort!=NULL)for(i=0;i<n;i++)tort[i]=-9.9;
		if(CLE!=NULL)for(i=0;i<n;i++)CLE[i]='R';
		if(Profile!=NULL)for(i=0;i<n;i++)for(k=0;k<Confo_Cluster;k++)Profile[i][k]=-1.0;
		if(n==2)
		{
			(r[1]-r[0]).xyz2double(cle_x_point[0]);
			x0=sqrt(dot(cle_x_point[0],cle_x_point[0]));
			if(dist!=NULL)dist[1]=x0;
		}
		return;
	}
	if(CLE!=NULL)for(i=0;i<n;i++)CLE[i]=' ';

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
	if(dist!=NULL)
	{
		dist[0]=-1.0;
		dist[1]=x0;
		dist[2]=x1;
	}
	if(bend!=NULL)
	{
		bend[0]=-9.9;
		bend[1]=-9.9;
		bend[2]=cle_w_point[0];
	}
	if(tort!=NULL)
	{
		tort[0]=-9.9;
		tort[1]=-9.9;
		tort[2]=-9.9;
	}
	if(CLE!=NULL)
	{
		CLE[0]='R';
		CLE[1]='R';
		CLE[n-1]='R';
		if(x0<CLE_CUT_MIN||x0>CLE_CUT_MAX)for(j=0;j<3;j++)if(0+j<n)CLE[0+j]='R';
		if(x1<CLE_CUT_MIN||x1>CLE_CUT_MAX)for(j=0;j<3;j++)if(1+j<n)CLE[1+j]='R';
	}
	if(Profile!=NULL)
	{
		for(k=0;k<Confo_Cluster;k++)Profile[0][k]=-1.0;
		for(k=0;k<Confo_Cluster;k++)Profile[1][k]=-1.0;
		for(k=0;k<Confo_Cluster;k++)Profile[n-1][k]=-1.0;
		if(x0<CLE_CUT_MIN||x0>CLE_CUT_MAX)for(j=0;j<3;j++)if(0+j<n)for(k=0;k<Confo_Cluster;k++)Profile[0+j][k]=-1.0;
		if(x1<CLE_CUT_MIN||x1>CLE_CUT_MAX)for(j=0;j<3;j++)if(1+j<n)for(k=0;k<Confo_Cluster;k++)Profile[1+j][k]=-1.0;
	}
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
		if(CLE!=NULL)
		{
			if(CLE[i-1]!='R')CLE[i-1]=btb_stc(cle_w_point);
			if(x1<CLE_CUT_MIN||x1>CLE_CUT_MAX)for(j=0;j<3;j++)if(i-1+j<n)CLE[i-1+j]='R';
		}
		if(Profile!=NULL)
		{
			if(Profile[i-1][0]!=-1.0)profile_stc(cle_w_point,Profile[i-1]);
			if(x1<CLE_CUT_MIN||x1>CLE_CUT_MAX)for(j=0;j<3;j++)if(i-1+j<n)for(k=0;k<Confo_Cluster;k++)Profile[i-1+j][k]=-1.0;
		}
		//update
		x0=x1;
		y0=y1;
		for(j=0;j<3;j++)cle_x_point[0][j]=cle_x_point[1][j];
		for(j=0;j<3;j++)cle_y_point[0][j]=cle_y_point[1][j];
		cle_w_point[0]=cle_w_point[2];
	}
}



//-----------------------//---------------------------------------------------//
//----Confo_Lett->XYZ---//
//---------------------//
void Confo_Lett::Get_Initial(XYZ *in1,XYZ *in2,double r1[3][3],double r2[3][3])
{
	int i;
	double phi;
	double thi;


	//--get_first_rotation--//	
	(in1[1]-in1[0]).xyz2double(cle_v1);  //v1 -> ori_coordinate
	(in2[1]-in2[0]).xyz2double(cle_v2);  //v2 -> pos_coordinate
	
	phi=dot(cle_v1,cle_v2)/sqrt(dot(cle_v1,cle_v1)*dot(cle_v2,cle_v2));
	phi=limit(phi,-1.0,1.0);
	phi=acos(1.0*phi);
	cross(cle_Axis,cle_v2,cle_v1);
	Universal_Rotation(cle_Axis,phi,r1);


	//--get_second_rotation--//	
	(in1[2]-in1[0]).xyz2double(cle_w1);  //w1 -> ori_coordinate
	(in2[2]-in2[0]).xyz2double(cle_w2);  //w2 -> pos_coordinate

	Vector_Multiply(cle_test,r1,cle_w2);
	for(i=0;i<3;i++)cle_w2[i]=cle_test[i];
	cross(cle_x1,cle_v1,cle_w2);
	cross(cle_x2,cle_v1,cle_w1);
	cross(cle_test,cle_x1,cle_x2);

	thi=dot(cle_x1,cle_x2)/sqrt(dot(cle_x1,cle_x1)*dot(cle_x2,cle_x2));
	thi=limit(thi,-1.0,1.0);
	thi=acos(1.0*thi);
	if(dot(cle_test,cle_v1)<0.0)thi*=-1.0;
	Universal_Rotation(cle_v1,thi,r2);
}


//input Confo_Angle (*dist,*bend,*tort), output XYZ_Type data_structure (*r)
//dist[0]    : -1.0;
//dist[1]-[n]: normal
//bend[0]-[1]: -9.9
//bend[2]-[n]: normal
//tort[0]-[2]: -9.9
//tort[3]-[n]: normal
void Confo_Lett::ctc_ori(double *dist,double *bend,double *tort,int n,XYZ *r,XYZ *init)
{
	int i;
	double radi;
	XYZ xyz;
	xyz=0.0;
	XYZ ori[3],pos[3];	

	//real_process//
	Matrix_Unit(cle_T_Back,1.0);
	for(i=0;i<n;i++)
	{
		if(dist==NULL)radi=3.8;
		else radi=dist[i];

		if(i==0)r[i]=0.0;
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

	//whether do ini_vector
	if(init!=NULL && n>=3)
	{
		for(i=0;i<3;i++)
		{
			ori[i]=init[i];
			pos[i]=r[i];
		}
		Get_Initial(ori,pos,cle_r1,cle_r2);
		Matrix_Multiply(cle_T_Pre,cle_r2,cle_r1);
		for(i=0;i<n;i++)
		{
			r[i].xyz2double(cle_D_Back);			
			Vector_Multiply(cle_D_Pre,cle_T_Pre,cle_D_Back);
			xyz.double2xyz(cle_D_Pre);
			r[i]=ori[0]+xyz;
		}
	}
}

//-------------- check ------------//__100205__//
//Check Space_Barrier
int Confo_Lett::Confo_Check_Space(XYZ *bak_mol,XYZ cur_mol,int moln,double cutoff)
{
	int i;
	double dist;
	double dou_cutoff=cutoff*cutoff;
	int num;
	num=0;
	for(i=0;i<moln;i++)
	{
		dist=bak_mol[i].distance_square(cur_mol);
		if(dist<dou_cutoff)num++;
	}
	return num;
}
//Check Confo_Lett
int Confo_Lett::Confo_Check_Lett(int code,double *input)
{
	int i,k;
	double ori;
	double val;
	int num;
	int real;
	//test
	real=btb_stc(input)-'A';
	//calc_ori
	k=Ant_id[code];
	ori=Ori_pi[k]*gaussian(k,input);
	//calc_other
	num=0;
	for(i=0;i<Confo_Cluster;i++) 
	{
		if(i==code)continue;
		k=Ant_id[i];
		val=Ori_pi[k]*gaussian(k,input);
		if(val>ori)num++;
	}
	return num;
}
