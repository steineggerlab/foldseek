#pragma once
#include "Computation_Utility.h"
#include "PDB_Chain.h"
#include "PhiPsi_Trans.h"
#include "XYZ_DBT.h"
#include "Confo_Lett.h"
#include "Hydro_Bond.h"
#include <vector>
#include <string>
using namespace std;

//====class_PDB_Chain====// 
class PDB_Chain_Fold : virtual public PDB_Chain
{
public:
	PDB_Chain_Fold();
	~PDB_Chain_Fold();
	PDB_Chain_Fold(const PDB_Chain_Fold & chain);
	PDB_Chain_Fold & operator =(const PDB_Chain_Fold & chain);
	PDB_Chain_Fold(const PDB_Chain & chain);
	PDB_Chain_Fold & operator =(const PDB_Chain & chain);

//======== macros ========//
public:
	int SSE_DIGIT;  //[0: 3-digit][1: 8-digit]

//-------- virtual_reloaded_function --------//
public:
	//phi_psi <-> XYZ
	virtual int refold_from_phi_psi_omega(void); //from phi_psi_omega to 3D-coordinates
	virtual int refold_from_theta_tau(void);     //from theta_tau to 3D-coordinates
	virtual int calculate_phi_psi_omega(void);   //3D-coordinates to phi_psi_omega
	virtual int calculate_theta_tau(void);       //3D-coordinates to theta_tau
	//calculate CLE,SSE and SEQ
	virtual int calculate_CLE(void);  //calculate Conformational Letter (17 digit)
	virtual int calculate_SSE(void);  //calculate DSSP Secondary Structure (8 digit)
};
