#pragma once
#include "XYZ.h"
#include "Computation_Utility.h"

//====class: Kabsch====//
//=> Kabsch (1976) superposition method
class Kabsch
{
public:
	Kabsch(void);
	~Kabsch(void);
	Kabsch(const Kabsch & kabsch);
	Kabsch & operator =(const Kabsch &xyz);

//----------temp_variables----------//
public:
	//macros
	int Using_Qcprot_Or_Not; //using qcprot or not //default: 1 (yes!)
	//original_kabsch
	double kabsch_rot[4][3];
	double kabsch_r[3][3];	
	double kabsch_rr[3][3];
	double kabsch_rtrans[3][3];
	double kabsch_atrans[3][3];	
	double kabsch_a[3][3];
	double kabsch_b[3][3];
	double kabsch_u[3][3];
	double kabsch_xc[3];
	double kabsch_yc[3];
	double kabsch_mu[3];
	double kabsch_yy[3];
	double kabsch_c[4];
	double kabsch_roots[3];
	double kabsch_tt[3];
	//expansion_kabsch
	int test_len;
	XYZ test_xc,test_yc;
	double test_X2,test_Y2;
	double test_XY[3][3];

//--------------------main_function--------------------//
public:
	//PART_I:(kabsch_rotate)
	//original
	double kabsch_base(double r[3][3],XYZ xc_,XYZ yc_,double *rotmat_);
	double kabsch_main(double *rotmat_);
	double kabsch(XYZ *strBuf1,XYZ *strBuf2,int lali,double *rotmat_);
	//expansion
	double kabsch_elong(XYZ p1,XYZ p2,double *rotmat_);
	double kabsch_shrink(XYZ p1,XYZ p2,double *rotmat_);

	//PART_I:(cubic_roots)
	int cubic_roots(double c[4],double s[3]);
	int eigen_values_simp(double m[3][3],double values[3]);
	int eigen_values_comp(double r[3][3],double m[3][3],double values[3],double *rotmat_);
	int eigen_values(double m[3][3],double values[3],double vectors[3][3]);
	//PART_II:(rot_mol)
	double calc_dist(XYZ *molA,XYZ *molB,int nAtom,double *rotmat_);
	double calc_dist(XYZ *molA,XYZ *molB,int nAtom);
	void rot_mol(XYZ *m_old,XYZ *m_new,int nAtom,double *rotmat_);
	void rot_point(XYZ p_old,XYZ &p_new,double *rotmat_);


	//==================== the following is qcprot ===========================//__110720__//
	double InnerProduct(double *A, XYZ *coords1, XYZ *coords2, XYZ &cen1, XYZ &cen2, int len, double *weight=0);
	/* Calculate the inner product of two structures.
	If weight array is not NULL, calculate the weighted inner product.

			Input: 
				coords1 -- reference structure
				coords2 -- candidate structure
				len     -- the size of the system
				weight  -- the weight array of size len: set to NULL if not needed
			Output:
				A[9]    -- the inner product matrix
				cen1    -- centroid of the first input
				cen2    -- centroid of the second input
			Return:
					(G1 + G2) * 0.5; used as E0 in function 'FastCalcRMSDAndRotation'
	*/

	int FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, int len, double minScore);
	/* Calculate the RMSD, and/or the optimal rotation matrix.

			Input:
					A[9]    -- the inner product of two structures
					E0      -- (G1 + G2) * 0.5
					len     -- the size of the system
					minScore-- if( minScore > 0 && rmsd < minScore) then calculate only the rmsd; 
							otherwise, calculate both the RMSD & the rotation matrix
			Output:
					rot[9]   -- the rotation matrix in the order of xx, xy, xz, yx, yy, yz, zx, zy, zz
					rmsd     -- the RMSD value
			Return: 
					only the rmsd was calculated if < 0
					both the RMSD & rotational matrix calculated if > 0
	*/

	double CalcRMSDRotationalMatrix(XYZ *coords1, XYZ *coords2, int len, double *rot, double *weight=0);
	/* Calculate the RMSD & rotational matrix.

			Input: 
				coords1 -- reference structure
				coords2 -- candidate structure
				len     -- the size of the system
				weight  -- the weight array of size len; set to NULL if not needed
			Output:
				rot[9]  -- rotation matrix
			Return:
				RMSD value

	*/

};
//====class_Kabsch====//over
