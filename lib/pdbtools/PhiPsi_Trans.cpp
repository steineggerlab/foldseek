#include "PhiPsi_Trans.h"

//--------------------- Start -----------------//
PhiPsi_Trans::PhiPsi_Trans(int num)
{
	PhiPsi_MAXIMAL=num;
	//backbone
	trs_dist_CN=1.32868493;
	trs_dist_NCa=1.45800149;
	trs_dist_CaC=1.52325821;
	trs_bend_CaCN=1.11352;
	trs_bend_CNCa=1.01753;
	trs_bend_NCaC=1.20079;
	//cb+co
	Confo_Back_CB_Dis=1.52083504;
	Confo_Back_CB_Ang=110.4*M_PI/180.0;    // Ca(i)->CB bend
	Confo_Back_CB_Tor=124.0*M_PI/180.0;    // Ca(i)->CB tort
	dist_CO=1.23101556;
	bend_CO=1.03323;
	//init
	PhiPsi_Init(PhiPsi_MAXIMAL);
	PhiPsi_Get_Init(PhiPsi_MAXIMAL);
}
PhiPsi_Trans::~PhiPsi_Trans(void)
{
	PhiPsi_Dele();
}

//-------------------- init -------------------//
void PhiPsi_Trans::PhiPsi_Init(int PhiPsi_MAXIMAL)
{
	PhiPsi_dist=new double[3*PhiPsi_MAXIMAL];
	PhiPsi_bend=new double[3*PhiPsi_MAXIMAL];
	PhiPsi_tort=new double[3*PhiPsi_MAXIMAL];
	PhiPsi_mol=new XYZ[3*PhiPsi_MAXIMAL];
}
void PhiPsi_Trans::PhiPsi_Dele(void)
{
	delete [] PhiPsi_dist;
	delete [] PhiPsi_bend;
	delete [] PhiPsi_tort;
	delete [] PhiPsi_mol;
}
void PhiPsi_Trans::PhiPsi_Get_Init(int moln)
{
	int i;
	for(i=0;i<moln;i++)
	{
		PhiPsi_dist[i*3+0]=trs_dist_CN;
		PhiPsi_dist[i*3+1]=trs_dist_NCa;
		PhiPsi_dist[i*3+2]=trs_dist_CaC;
		PhiPsi_bend[i*3+0]=trs_bend_CaCN;
		PhiPsi_bend[i*3+1]=trs_bend_CNCa;
		PhiPsi_bend[i*3+2]=trs_bend_NCaC;
	}
	PhiPsi_dist[0]=-9.9;
	PhiPsi_bend[0]=-9.9;
	PhiPsi_bend[1]=-9.9;
}

//-------------------- data -------------------//
void PhiPsi_Trans::PhiPsi_Get_Data(int moln,double *phipsi)
{
	int i;
	double phi,psi,ome;
	int count;
	//init
	count=2;
	for(i=0;i<moln;i++)
	{
		phi=phipsi[i*3+0];
		psi=phipsi[i*3+1];
		ome=phipsi[i*3+2];
		PhiPsi_tort[count+0]=phi;
		PhiPsi_tort[count+1]=psi;
		PhiPsi_tort[count+2]=ome;
		count+=3;
	}
	PhiPsi_tort[0]=-9.9;
	PhiPsi_tort[1]=-9.9;
	PhiPsi_tort[2]=-9.9;
}
void PhiPsi_Trans::PhiPsi_Put_Data(int moln,double *phipsi)
{
	int i;
	double phi,psi,ome;
	int count;
	//init
	count=2;
	for(i=0;i<moln-1;i++)
	{
		phi=PhiPsi_tort[count+0];
		psi=PhiPsi_tort[count+1];
		ome=PhiPsi_tort[count+2];
		phipsi[i*3+0]=phi;
		phipsi[i*3+1]=psi;
		phipsi[i*3+2]=ome;
		count+=3;
	}
	if(i==0)phi=M_PI;
	else phi=PhiPsi_tort[count+0];
	phipsi[i*3+0]=phi;
	phipsi[i*3+1]=M_PI;
	phipsi[i*3+2]=M_PI;
	phipsi[0*3+0]=M_PI;
}

//-------------------- recon -------------------//
void PhiPsi_Trans::PhiPsi_Recon_CB_From_NCaC(XYZ N,XYZ Ca,XYZ C,XYZ &Cb)
{
	XYZ xyz;
	double x[3],y[3];
	double z[3];
	double cb1[3],cb2[3];
	double dist,bend,tort;
	dist=Confo_Back_CB_Dis;
	bend=Confo_Back_CB_Ang;
	tort=Confo_Back_CB_Tor;
	//process
	xyz=C-Ca;
	xyz.xyz2double(x);
	xyz=Ca-N;
	xyz.xyz2double(y);
	cross(z,x,y);
	Universal_Rotation(z,-1.0*bend,PhiPsi_ROT_TOT);
	Vector_Multiply(cb1,PhiPsi_ROT_TOT,x);
	Universal_Rotation(x,-1.0*tort,PhiPsi_ROT_TOT);
	Vector_Multiply(cb2,PhiPsi_ROT_TOT,cb1);
	Vector_Normalize(cb2,3);
	Vector_Dot(cb2,dist,cb2,3);
	Cb.double2xyz(cb2);
	Cb+=Ca;
}
void PhiPsi_Trans::PhiPsi_Recon_CO_From_CaCN(XYZ Ca,XYZ C,XYZ N,XYZ &CO)
{
	XYZ xyz;
	double x[3],y[3];
	double z[3];
	double cb2[3];
	double dist,bend;
	dist=dist_CO;
	bend=1.0*M_PI-bend_CO;
	//process
	xyz=Ca-C;
	xyz.xyz2double(x);
	xyz=C-N;
	xyz.xyz2double(y);
	cross(z,x,y);
	Universal_Rotation(z,bend,PhiPsi_ROT_TOT);
	Vector_Multiply(cb2,PhiPsi_ROT_TOT,x);
	Vector_Normalize(cb2,3);
	Vector_Dot(cb2,dist,cb2,3);
	CO.double2xyz(cb2);
	CO+=C;
}

//-------------------- calc phi_psi main -------------------//
//double *phipsi -> 3*moln
//XYZ **phipsi   -> moln x 5 [N,CA,C,O,CB]
void PhiPsi_Trans::PhiPsi_To_BackCB(int moln,char *ami,double *phipsi,XYZ **mol)
{
	int i;
	XYZ N1,CA,CO,N2;
	PhiPsi_Get_Data(moln,phipsi);
	xyz_to_dbt.dbt_to_xyz(PhiPsi_dist,PhiPsi_bend,PhiPsi_tort,3*moln,PhiPsi_mol);
	for(i=0;i<moln;i++)
	{
		mol[i][0]=PhiPsi_mol[i*3+0]; //N
		mol[i][1]=PhiPsi_mol[i*3+1]; //CA
		mol[i][2]=PhiPsi_mol[i*3+2]; //C
		if(i==moln-1)
		{
			N1=PhiPsi_mol[i*3+0];
			CA=PhiPsi_mol[i*3+1];
			CO=PhiPsi_mol[i*3+2];
			mol[i][3]=mol[i-1][3]-mol[i-1][2]+mol[i][2];          //O
			if(ami[i]=='G')mol[i][4]=mol[i][1];
			else PhiPsi_Recon_CB_From_NCaC(N1,CA,CO,mol[i][4]);   //CB
		}
		else
		{
			N1=PhiPsi_mol[i*3+0];
			CA=PhiPsi_mol[i*3+1];
			CO=PhiPsi_mol[i*3+2];
			N2=PhiPsi_mol[i*3+3];
			PhiPsi_Recon_CO_From_CaCN(CA,CO,N2,mol[i][3]);        //O
			if(ami[i]=='G')mol[i][4]=mol[i][1];
			else PhiPsi_Recon_CB_From_NCaC(N1,CA,CO,mol[i][4]);   //CB
		}
	}
}
void PhiPsi_Trans::BackCB_To_PhiPsi(int moln,char *ami,double *phipsi,XYZ **mol)
{
	int i,j;
	for(i=0;i<moln;i++)for(j=0;j<3;j++)PhiPsi_mol[i*3+j]=mol[i][j];
	xyz_to_dbt.xyz_to_dbt(PhiPsi_dist,PhiPsi_bend,PhiPsi_tort,3*moln,PhiPsi_mol);
	PhiPsi_Put_Data(moln,phipsi);
}
