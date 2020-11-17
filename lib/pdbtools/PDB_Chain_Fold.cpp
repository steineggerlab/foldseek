#include "PDB_Chain_Fold.h"
#include "Utility.h"
#include "PhiPsi_Trans.h"
#include "XYZ_DBT.h"
#include "Confo_Lett.h"
#include "Hydro_Bond.h"
#include <iomanip>
using namespace std;


//-------- PDB_Chain_Fold ---------// 
PDB_Chain_Fold::PDB_Chain_Fold()
:PDB_Chain()
{
	SSE_DIGIT=1; //default 8-digit
}
PDB_Chain_Fold::~PDB_Chain_Fold()
{
}
PDB_Chain_Fold::PDB_Chain_Fold(const PDB_Chain_Fold & chain)
:PDB_Chain(chain)
{
}
PDB_Chain_Fold &  PDB_Chain_Fold::operator =(const PDB_Chain_Fold & chain)
{
	if(this==&chain)return *this;
	PDB_Chain::operator=(chain);
	return *this;
}
PDB_Chain_Fold::PDB_Chain_Fold(const PDB_Chain & chain)
:PDB_Chain(chain)
{
}
PDB_Chain_Fold &  PDB_Chain_Fold::operator =(const PDB_Chain & chain)
{
	if(this==&chain)return *this;
	PDB_Chain::operator=(chain);
	return *this;
}

//====================== virtual functions ===========================//
//[phi_psi -> XYZ]
int PDB_Chain_Fold::refold_from_phi_psi_omega(void)
{
	PDB_Residue residue;
	_Phi_Psi_Omega phi_psi_ome;
	double *torsion_angle;
	XYZ	**decoy;
	char *ami;
	int i;
	int retv;
	int wsbad=0;
	int	n = get_length();
	PhiPsi_Trans phipsi_trans(n);
	string s=get_sequence();
	// initialize bond angles, and torsion angles;
	torsion_angle = new double[ 3*n ];
	ami=new char[n+1];
	NewArray2D(&decoy,n,5);
	// format phi-psi to torsion angles
	torsion_angle[0]=-9.9;
	torsion_angle[1]=-9.9;
	torsion_angle[2]=-9.9;
	for(i=0;i<n-1;i++)
	{
		retv=get_phi_psi_omega(i,phi_psi_ome);
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		torsion_angle[3*i+2] = phi_psi_ome.phi;
		torsion_angle[3*i+3] = phi_psi_ome.psi;
		torsion_angle[3*i+4] = phi_psi_ome.omega;
	}
	retv=get_phi_psi_omega(n-1,phi_psi_ome);
	if(retv<0)
	{
		wsbad=1;
		goto end;
	}
	torsion_angle[3*n-1] = phi_psi_ome.phi;
	//ami
	if((int)s.size()!=n)
	{
		wsbad=1;
		goto end;
	}
	strcpy(ami,s.c_str());
	//fold
	phipsi_trans.PhiPsi_To_BackCB(n,ami,torsion_angle,decoy);
	//set
	for(i=0;i<n;i++) 
	{
		retv=get_residue(i,residue);
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		debug_info( DEBUG_DEFAULT, "before N \n");
		retv=residue.set_backbone_atom( "N  ",decoy[i][0] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		debug_info( DEBUG_DEFAULT, "before CA \n");
		retv=residue.set_backbone_atom( "CA ",decoy[i][1] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		debug_info( DEBUG_DEFAULT, "before C \n");
		retv=residue.set_backbone_atom( "C  ",decoy[i][2] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		debug_info( DEBUG_DEFAULT, "before CA \n");
		retv=residue.set_backbone_atom( "O  ",decoy[i][3] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		debug_info( DEBUG_DEFAULT, "before C \n");
		retv=residue.set_sidechain_atom( "CB ",decoy[i][4] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=set_residue(i,residue); //set
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
	}
	//delete
end:
	delete [] torsion_angle;
	delete [] ami;
	DeleteArray2D(&decoy,n);
	if(wsbad==0)return 0;
	else return -1;
}
//[the_tau -> XYZ]
int PDB_Chain_Fold::refold_from_theta_tau(void)
{
	PDB_Residue residue;
	_Theta_Tau theta_tau;
	XYZ_DBT xyz_dbt;
	double	*bond_length, *torsion_angle, *bond_angle;
	XYZ	*decoy;
	//init
	double Pseudo_Bond_Dis=3.8;
	int i;
	int retv;
	int wsbad=0;
	int n = get_length();
	bond_length = new double[n];
	bond_angle = new double[n];
	torsion_angle = new double[n];
	decoy = new XYZ[n];
	//process
	for(i = 0; i < n; i++)
	{
		retv=get_theta_tau(i,theta_tau);
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		bond_length[i] = Pseudo_Bond_Dis;
		bond_angle[i] = theta_tau.theta;
		torsion_angle[i] = theta_tau.tau;
	}
	xyz_dbt.dbt_to_xyz(bond_length, bond_angle, torsion_angle, n, decoy);
	//record
	for(i = 0; i < n; i++)
	{
		retv=get_residue(i,residue);
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=residue.set_backbone_atom( "CA ",decoy[i] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=set_residue(i,residue); //set
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
	}
	//delete
end:
	delete [] bond_length;
	delete [] bond_angle;
	delete [] torsion_angle;
	delete [] decoy;
	if(wsbad==0)return 0;
	else return -1;
}
//[phi_psi <- XYZ]
int PDB_Chain_Fold::calculate_phi_psi_omega(void)
{
	PDB_Residue residue;
	_Phi_Psi_Omega phi_psi_ome;
	XYZ_DBT xyz_dbt;
	int i;
	int retv;
	int wsbad=0;
	int n = get_length();
	XYZ * r = new XYZ[3 * n];
	double * bond_length = new double[3 * n];
	double * bond_angle = new double[3 * n];
	double * torsion_angle = new double[3 * n];
	//get
	for(i = 0; i < n; i++)
	{
		retv=get_residue(i,residue);
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=residue.get_backbone_atom( "N  ",r[3*i+0] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=residue.get_backbone_atom( "CA ",r[3*i+1] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=residue.get_backbone_atom( "C  ",r[3*i+2] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
	}
	//trans
	xyz_dbt.xyz_to_dbt(bond_length, bond_angle, torsion_angle, 3 * n, r);
	//calc
	for(i = 0; i < n-1; i++)
	{
		retv=get_phi_psi_omega(i,phi_psi_ome);
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		if(i==0)phi_psi_ome.phi=M_PI;
		else phi_psi_ome.phi = torsion_angle[i * 3 + 2];
		phi_psi_ome.psi = torsion_angle[i * 3 + 3];
		phi_psi_ome.omega = torsion_angle[i * 3 + 4];
		retv=set_phi_psi_omega(i,phi_psi_ome); //set
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
	}
	//-> terminal
	double phi;
	if(i==0)phi=M_PI;
	else phi=torsion_angle[3*n - 1];
	retv=get_phi_psi_omega(n-1,phi_psi_ome);
	if(retv<0)
	{
		wsbad=1;
		goto end;
	}
	phi_psi_ome.phi=phi;
	phi_psi_ome.psi=M_PI;
	phi_psi_ome.omega=M_PI;
	retv=set_phi_psi_omega(n-1,phi_psi_ome);   //set
	if(retv<0)
	{
		wsbad=1;
		goto end;
	}
	//delete
end:
	delete [] bond_length;
	delete [] bond_angle;
	delete [] torsion_angle;
	delete [] r;
	if(wsbad==0)return 0;
	else return -1;
}
//[the_tau <- XYZ]
int PDB_Chain_Fold::calculate_theta_tau(void)
{
	PDB_Residue residue;
	_Theta_Tau theta_tau;
	XYZ_DBT xyz_dbt;
	int i;
	int retv;
	int wsbad=0;
	int n = get_length();
	XYZ * r = new XYZ[n];
	double * bond_length = new double[n];
	double * bond_angle = new double[n];
	double * torsion_angle = new double[n];
	//get
	for(i = 0; i < n; i++)
	{
		retv=get_residue(i,residue);
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=residue.get_backbone_atom( "CA ",r[i] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
	}
	//trans
	xyz_dbt.xyz_to_dbt(bond_length, bond_angle, torsion_angle, n, r);
	//calc
	for(i = 0; i < n; i++)
	{
		retv=get_theta_tau(i,theta_tau);
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		if(i<2)theta_tau.theta = M_PI;
		else theta_tau.theta = bond_angle[i];
		if(i<3)theta_tau.tau = M_PI;
		else theta_tau.tau = torsion_angle[i];
		retv=set_theta_tau(i,theta_tau);    //set
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
	}
	//delete
end:
	delete [] bond_length;
	delete [] bond_angle;
	delete [] torsion_angle;
	delete [] r;
	if(wsbad==0)return 0;
	else return -1;
}

//calculate CLE,SSE and SEQ
//CLE
int PDB_Chain_Fold::calculate_CLE(void)
{
	PDB_Residue residue;
	Confo_Lett confo_lett;
	string CLE_;
	int i;
	int retv;
	int wsbad=0;
	int n = get_length();
	XYZ *r=new XYZ[n];
	char *cle=new char[n+1];
	//CA
	for(i=0; i<n; i++)
	{
		retv=get_residue(i,residue);
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=residue.get_backbone_atom( "CA ",r[i] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
	}
	//CLE
	confo_lett.btb_ori(0,0,0,n,r,cle);
	cle[n]='\0';
	CLE_.assign(cle);
	retv=set_CLE(CLE_,0);
	if(retv<0)
	{
		wsbad=1;
		goto end;
	}
end:
	delete [] r;
	delete [] cle;
	if(wsbad==0) return 0;
	else return -1;
}
//SSE
int PDB_Chain_Fold::calculate_SSE(void)
{
	int n = get_length();
	PDB_Residue residue;
	Hydro_Bond hydro_bond(n);
	string s=get_sequence();
	string SSE_;
	XYZ **mola;
	NewArray2D(&mola,n,5);
	char *ami=new char[n+1];
	char *sse=new char[n+1];
	int i;
	int retv;
	int wsbad=0;
	//backbone
	for(i=0;i<n;i++)
	{
		retv=get_residue(i,residue);
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=residue.get_backbone_atom( "N  ",mola[i][0] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=residue.get_backbone_atom( "CA ",mola[i][1] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=residue.get_backbone_atom( "C  ",mola[i][2] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=residue.get_backbone_atom( "O  ",mola[i][3] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
		retv=residue.get_sidechain_atom( "CB ",mola[i][4] );
		if(retv<0)
		{
			wsbad=1;
			goto end;
		}
	}
	//ami
	if((int)s.size()!=n)
	{
		wsbad=1;
		goto end;
	}
	strcpy(ami,s.c_str());
	//SSE
	hydro_bond.HB_Input_Mol(mola,ami,n);
	hydro_bond.HB_Calc_Hydro_Bond();
	hydro_bond.HB_Calc_SSE(sse);                        //-> return 8-digit SSE
	if(SSE_DIGIT==0)hydro_bond.HB_Trans_SSE(sse,sse,n); //-> transfer to 3-digit SSE
	SSE_.assign(sse);
	retv=set_SSE(SSE_,0);
	if(retv<0)
	{
		wsbad=1;
		goto end;
	}
end:
	DeleteArray2D(&mola,n);
	delete [] ami;
	delete [] sse;
	if(wsbad==0)return 0;
	else return -1;
}
