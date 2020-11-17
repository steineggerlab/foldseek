#include "Kabsch.h"
#include "XYZ.h"

//------ constructor -------//
Kabsch::Kabsch(void)
{
	Using_Qcprot_Or_Not=1; //default: using qcprot
}
Kabsch::~Kabsch(void)
{
}
Kabsch::Kabsch(const Kabsch & kabsch)
{
	this->test_len=kabsch.test_len;
	this->test_xc=kabsch.test_xc;
	this->test_yc=kabsch.test_yc;
	this->test_X2=kabsch.test_X2;
	this->test_Y2=kabsch.test_Y2;
	for(int i=0;i<3;i++)for(int j=0;j<3;j++)this->test_XY[i][j]=kabsch.test_XY[i][j];
}
Kabsch & Kabsch::operator =(const Kabsch & kabsch) 
{
	if(this==&kabsch)return *this;
	this->test_len=kabsch.test_len;
	this->test_xc=kabsch.test_xc;
	this->test_yc=kabsch.test_yc;
	this->test_X2=kabsch.test_X2;
	this->test_Y2=kabsch.test_Y2;
	for(int i=0;i<3;i++)for(int j=0;j<3;j++)this->test_XY[i][j]=kabsch.test_XY[i][j];
	return *this;
}

//==================================== the following is the main_part ===========================//
//[0]Kabsch_Basic//(original_version)
double Kabsch::kabsch_base(double r[3][3],XYZ xc_,XYZ yc_,double *rotmat_)
{
/*
 * Fit pairs of points using the method described in:
 *
 *  A discussion of the solution for the best rotation to relate
 *  two sets of vectors', W.Kabsch, Acta Cryst. (1978), A34, 827-828.
 *
 * The method generates the 4x4 matrix rot which will map the
 * coordinates y on to the coordinates x.
 *
 * The matrix is arranged such that rot[3][0], rot[3][1], rot[3][2]
 * are the translational components.
 *
 * The function returns the root mean square distance between
 * x and y.  The coordinates in x and y are not modified.
 */

	int i;
	double result;
	int tag;
	int nroots;


	/* Generate the transpose of r and form rtrans x r */  
	Matrix_Trans(kabsch_rtrans,r);
	Matrix_Multiply(kabsch_rr,r,kabsch_rtrans);

	//eigen_values
	if(rotmat_==NULL)nroots=eigen_values_simp(kabsch_rr,kabsch_mu);
	else
	{
		kabsch_xc[0]=xc_.X;
		kabsch_xc[1]=xc_.Y;
		kabsch_xc[2]=xc_.Z;
		kabsch_yc[0]=yc_.X;
		kabsch_yc[1]=yc_.Y;
		kabsch_yc[2]=yc_.Z;
		nroots=eigen_values_comp(r,kabsch_rr,kabsch_mu,rotmat_);
	}

	//---Get_Mu--//__080815__//
	if(Matrix_Determinant(r)>0.0)tag=1;
	else tag=0;

	result=0.0;
	if(nroots==3)
	{
		for(i=0;i<2;i++)result+=sqrt(fabs(kabsch_mu[i]));
		if(tag==1)result+=sqrt(fabs(kabsch_mu[i]));
		else result-=sqrt(fabs(kabsch_mu[i]));
	}
	else for(i=0;i<nroots;i++)result+=sqrt(fabs(kabsch_mu[i]));
	return result;
}

//[1]standart_ali//
//main
double Kabsch::kabsch_main(double *rotmat_)
{
	XYZ xc,yc;
	XYZ zero;
	double xtot,ytot;
	double rmsd,result;

	//mean
	xc=test_xc;
	yc=test_yc;
	//square
	zero=0.0;
	xtot=test_X2-test_len*test_xc.distance_square(zero);
	ytot=test_Y2-test_len*test_yc.distance_square(zero);
	//vari
	kabsch_r[0][0]=test_XY[0][0]-test_len*test_yc.X*test_xc.X;
	kabsch_r[0][1]=test_XY[0][1]-test_len*test_yc.X*test_xc.Y;
	kabsch_r[0][2]=test_XY[0][2]-test_len*test_yc.X*test_xc.Z;
	kabsch_r[1][0]=test_XY[1][0]-test_len*test_yc.Y*test_xc.X;
	kabsch_r[1][1]=test_XY[1][1]-test_len*test_yc.Y*test_xc.Y;
	kabsch_r[1][2]=test_XY[1][2]-test_len*test_yc.Y*test_xc.Z;
	kabsch_r[2][0]=test_XY[2][0]-test_len*test_yc.Z*test_xc.X;
	kabsch_r[2][1]=test_XY[2][1]-test_len*test_yc.Z*test_xc.Y;
	kabsch_r[2][2]=test_XY[2][2]-test_len*test_yc.Z*test_xc.Z;
	
	//----- calc_rotmat ------//
	result=kabsch_base(kabsch_r,xc,yc,rotmat_);
	rmsd=(xtot+ytot-2*result)/test_len;
	if(IsZero(rmsd))rmsd=0.0;
	return rmsd;
}

//original
double Kabsch::kabsch(XYZ *strBuf1,XYZ *strBuf2,int lali,double *rotmat_) 
{
	//init_check
	if(lali==0)return -1.0;

	//======== ws_using_qcprot ========//__110720__//
	if(Using_Qcprot_Or_Not==1)
	{
		return CalcRMSDRotationalMatrix(strBuf1,strBuf2,lali,rotmat_);
	}
	//======== ws_using_qcprot ========//__110720__//over


	XYZ zero;
	int i,j; 

	/* Find center of each set of coordinates. */ 
	test_len=lali;
	test_xc=0.0;
	test_yc=0.0;
	for(i=0;i<test_len;i++)
	{
		test_xc+=strBuf1[i];
		test_yc+=strBuf2[i];
	}
	test_xc/=test_len;
	test_yc/=test_len;

	//get xtot & ytot//
	zero=0.0;
	test_X2=0.0;
	test_Y2=0.0;
	for(i=0;i<test_len;i++)
	{
		test_X2+=strBuf1[i].distance_square(zero);
		test_Y2+=strBuf2[i].distance_square(zero);
	}

	/*
	* Initialise and then fill the r matrix.
	* Note that centre is subtracted at this stage.
	*/
	for(i=0;i<3;i++)for(j=0;j<3;j++)test_XY[i][j]=0.0;
	for(i=0;i<test_len;i++)
	{
		test_XY[0][0]+=strBuf2[i].X*strBuf1[i].X;
		test_XY[0][1]+=strBuf2[i].X*strBuf1[i].Y;
		test_XY[0][2]+=strBuf2[i].X*strBuf1[i].Z;
		test_XY[1][0]+=strBuf2[i].Y*strBuf1[i].X;
		test_XY[1][1]+=strBuf2[i].Y*strBuf1[i].Y;
		test_XY[1][2]+=strBuf2[i].Y*strBuf1[i].Z;
		test_XY[2][0]+=strBuf2[i].Z*strBuf1[i].X;
		test_XY[2][1]+=strBuf2[i].Z*strBuf1[i].Y;
		test_XY[2][2]+=strBuf2[i].Z*strBuf1[i].Z;
	}

	//input R,XC,YC to kabsch_basic to get the rot_mat(d)
	return kabsch_main(rotmat_);
}

//expansion
double Kabsch::kabsch_elong(XYZ p1,XYZ p2,double *rotmat_)
{
	XYZ zero;

	//======update_value======//
	//update_XC_YC
	test_xc=(test_xc*test_len+p1)/(double)(test_len+1);
	test_yc=(test_yc*test_len+p2)/(double)(test_len+1);
	test_len++;	
	//update_X2_Y2
	zero=0.0;
	test_X2+=p1.distance_square(zero);
	test_Y2+=p2.distance_square(zero);
	//update_XY[][]
	test_XY[0][0]+=p2.X*p1.X;
	test_XY[0][1]+=p2.X*p1.Y;
	test_XY[0][2]+=p2.X*p1.Z;
	test_XY[1][0]+=p2.Y*p1.X;
	test_XY[1][1]+=p2.Y*p1.Y;
	test_XY[1][2]+=p2.Y*p1.Z;
	test_XY[2][0]+=p2.Z*p1.X;
	test_XY[2][1]+=p2.Z*p1.Y;
	test_XY[2][2]+=p2.Z*p1.Z;

	//---- calc_rmsd ----//
	return kabsch_main(rotmat_);
}
double Kabsch::kabsch_shrink(XYZ p1,XYZ p2,double *rotmat_)
{
	XYZ zero;

	//======update_value======//
	//update_XC_YC
	test_xc=(test_xc*test_len-p1)/(double)(test_len-1);
	test_yc=(test_yc*test_len-p2)/(double)(test_len-1);
	test_len--;	
	//update_X2_Y2
	zero=0.0;
	test_X2-=p1.distance_square(zero);
	test_Y2-=p2.distance_square(zero);
	//update_XY[][]
	test_XY[0][0]-=p2.X*p1.X;
	test_XY[0][1]-=p2.X*p1.Y;
	test_XY[0][2]-=p2.X*p1.Z;
	test_XY[1][0]-=p2.Y*p1.X;
	test_XY[1][1]-=p2.Y*p1.Y;
	test_XY[1][2]-=p2.Y*p1.Z;
	test_XY[2][0]-=p2.Z*p1.X;
	test_XY[2][1]-=p2.Z*p1.Y;
	test_XY[2][2]-=p2.Z*p1.Z;

	//---- calc_rmsd ----//
	return kabsch_main(rotmat_);
}



//==================================== the following is the details ===========================//
//cubic_roots
int Kabsch::cubic_roots(double c[4],double s[3])
{
	int i, num;
	double sub;
	double A, B, C;
	double A2, p, q;
	double p3, D;

	/* normal form: x^3 + Ax^2 + Bx + C = 0 */
	A = c[2] / c[3];
	B = c[1] / c[3];
	C = c[0] / c[3];

	/*  substitute x = y - A/3 to eliminate quadric term:
		x^3 +px + q = 0 */

	A2 = A * A;
	p = 1.0 / 3 * (-1.0 / 3 * A2 + B);
	q = 1.0 / 2 * (2.0 / 27 * A * A2 - 1.0 / 3 * A * B + C);

	/* use Cardano's formula */
	p3 = p * p * p;
	D = q * q + p3;

	if (IsZero (D))
	{				
		/* one triple solution */
		if (IsZero (q))
		{
			s[0] = 0;
			num = 1;
		}
		else
		{			
			/* one single and one double solution */
			double u = cbrt (-q);
			s[0] = 2 * u;
			s[1] = -u;
			num = 2;
		}
	}
	else if (D < 0)
	{				
		/* Casus irreducibilis: three real solutions */
		double phi = 1.0 / 3 * acos (-q / sqrt (-p3));
		double t = 2 * sqrt (-p);

		s[0] = t * cos (phi);
		s[1] = -t * cos (phi + M_PI / 3);
		s[2] = -t * cos (phi - M_PI / 3);
		num = 3;
	}
	else
	{				
		/* one real solution */
		double sqrt_D = sqrt (D);
		double u = cbrt (sqrt_D - q);
		double v = -cbrt (sqrt_D + q);

		s[0] = u + v;
		num = 1;
	}

	/* resubstitute */
	sub = 1.0 / 3 * A;
	for (i = 0; i < num; ++i) s[i] -= sub;

	return num;
}

//eigen_values
int Kabsch::eigen_values_simp(double m[3][3],double values[3])
/*
 * Find the eigen values and vectors for the matrix m.
 */
{
	double a1 = m[0][0], b1 = m[0][1], c1 = m[0][2];
	double a2 = m[1][0], b2 = m[1][1], c2 = m[1][2];
	double a3 = m[2][0], b3 = m[2][1], c3 = m[2][2];
	double l;
	int nroots, iroot;
	int i;

	/*
	* Expanding the characteristic equation of a 3x3 matrix...
	* Maple tells us the following is the cubic in l
	*/

	kabsch_c[0] = a1 * (b2 * c3 - b3 * c2) + c1 * (a2 * b3 - a3 * b2) - b1 * (a2 * c3 - a3 * c2);
	kabsch_c[1] = (a1 * (-b2 - c3) - b2 * c3 + b3 * c2 + c1 * a3 + b1 * a2);
	kabsch_c[2] = (a1 + b2 + c3);
	kabsch_c[3] = -1.;

	nroots = cubic_roots (kabsch_c, kabsch_roots);

	/* Degenerate roots are not returned individually. */
	for (i = 0; i < nroots; i++)
	{
		iroot = i > nroots ? nroots : i;
		l = kabsch_roots[iroot];
		values[i] = l;
	}
	return nroots;
}
int Kabsch::eigen_values_comp(double r[3][3],double m[3][3],double values[3],double *rotmat_)
{
	/*
	* Get the eigenvalues and vectors.
	* Reform a[2] as cross product of a[0] and a[1] to ensure
	* right handed system.
	*/
	int i,j;
	int nroots=eigen_values(m,values,kabsch_a);
	cross(kabsch_a[2],kabsch_a[0],kabsch_a[1]);


	/* Transform first two eigenvectors and normalise them. */  
	for(i=0;i<2;i++)
	{
		Vector_Multiply_Trans(kabsch_b[i],r,kabsch_a[i]);
		Vector_Normalize(kabsch_b[i],3);
	}
	
	
	/* Make right handed set. */ 
	cross(kabsch_b[2],kabsch_b[0],kabsch_b[1]);
	
	
	/* Form the rotation matrix. */  
	Matrix_Trans(kabsch_atrans,kabsch_a);
	Matrix_Multiply(kabsch_u,kabsch_atrans,kabsch_b);
	for(i=0;i<3;i++)for(j=0;j<3;j++)kabsch_rot[i][j]=kabsch_u[i][j];  
	
	
	/* Transform offset of y coordinates by the rotation. */  
	Vector_Multiply_Trans(kabsch_yy,kabsch_u,kabsch_yc);
	
	
	/* Build translational component of rot from offsets. */  
	for(i=0;i<3;i++)kabsch_rot[3][i]=-kabsch_yy[i]+kabsch_xc[i];
	
	
	/* Figure out the rms deviation of the fitted coordinates. */  
	rotmat_[0] = kabsch_rot[0][0]; rotmat_[1] = kabsch_rot[1][0]; rotmat_[2] = kabsch_rot[2][0];   
	rotmat_[3] = kabsch_rot[0][1]; rotmat_[4] = kabsch_rot[1][1]; rotmat_[5] = kabsch_rot[2][1];  
	rotmat_[6] = kabsch_rot[0][2]; rotmat_[7] = kabsch_rot[1][2]; rotmat_[8] = kabsch_rot[2][2];   
	rotmat_[9] = kabsch_rot[3][0]; rotmat_[10] = kabsch_rot[3][1]; rotmat_[11] = kabsch_rot[3][2];

	//final
	return nroots;
}
int Kabsch::eigen_values(double m[3][3],double values[3],double vectors[3][3])
/*
 * Find the eigen values and vectors for the matrix m.
 */
{
	double a1 = m[0][0], b1 = m[0][1], c1 = m[0][2];
	double a2 = m[1][0], b2 = m[1][1], c2 = m[1][2];
	double a3 = m[2][0], b3 = m[2][1], c3 = m[2][2];
	double l;
	double x, y, z, norm;	
	int nroots, iroot;
	int i;

	/*
	* Expanding the characteristic equation of a 3x3 matrix...
	* Maple tells us the following is the cubic in l
	*/

	kabsch_c[0] = a1 * (b2 * c3 - b3 * c2) + c1 * (a2 * b3 - a3 * b2) - b1 * (a2 * c3 - a3 * c2);
	kabsch_c[1] = (a1 * (-b2 - c3) - b2 * c3 + b3 * c2 + c1 * a3 + b1 * a2);
	kabsch_c[2] = (a1 + b2 + c3);
	kabsch_c[3] = -1.;

	nroots = cubic_roots (kabsch_c, kabsch_roots);


	/* Degenerate roots are not returned individually. */
	for (i = 0; i < nroots; i++)
	{
		iroot = i > nroots ? nroots : i;
		l = kabsch_roots[iroot];
		values[i] = l;
	}
	for (i = 0; i < nroots; i++)
	{
		/*
		* Find the eigen vectors by solving pairs of the
		* three simultaneous equations.  From `Mathematical Methods
		* in Science and Engineering', Heiding, pg.19
		*
		* Sometimes we get x = y = z = 0.0, so try the other two
		* pairs of equations and hope that one of them gives a solution.
		*/

		iroot = i > nroots ? nroots : i;
		l = kabsch_roots[iroot];
		x = b1 * c2 - (b2 - l) * c1;
		y = -((a1 - l) * c2 - a2 * c1);
		z = ((a1 - l) * (b2 - l) - a2 * b1);

		if (IsZero (x) && IsZero (y) && IsZero (z))	
		{
			x = b1 * (c3 - l) - b3 * c1;
			y = -((a1 - l) * (c3 - l) - a3 * c1);
			z = ((a1 - l) * b3 - a3 * b1);		
			if (IsZero (x) && IsZero (y) && IsZero (z))
			{
				x = (b2 - l) * (c3 - l) - b3 * c2;
				y = -(a2 * (c3 - l) - a3 * c2);
				z = (a2 * b3 - a3 * (b2 - l));
				if (IsZero (x) && IsZero (y) && IsZero (z))
				{
					x=0.0;
					y=0.0;
					z=0.0;
				}
			}	
		}

		norm = sqrt (x * x + y * y + z * z);
		if (!IsZero (norm))
		{
			vectors[i][0] = x / norm;
			vectors[i][1] = y / norm;
			vectors[i][2] = z / norm;
		}
		else
		{
			vectors[i][0] = 0.0;
			vectors[i][1] = 0.0;
			vectors[i][2] = 0.0;
		}
	}
	return nroots;
}

/////////////////////////////////=============transformation_related=================//////////////////////////////////

//cRMS
double Kabsch::calc_dist(XYZ *molA,XYZ *molB,int nAtom,double *rotmat_)
{
	int l;
	double rmsd=0.0, dx, dy, dz, d;
	for(l=0;l<nAtom;l++) 
	{
		dx=molA[l].X; 
		dy=molA[l].Y;
		dz=molA[l].Z;
		d=molB[l].X-(dx*rotmat_[0]+dy*rotmat_[1]+dz*rotmat_[2]+rotmat_[9]);
		rmsd+=d*d;
		d=molB[l].Y-(dx*rotmat_[3]+dy*rotmat_[4]+dz*rotmat_[5]+rotmat_[10]);
		rmsd+=d*d;
		d=molB[l].Z-(dx*rotmat_[6]+dy*rotmat_[7]+dz*rotmat_[8]+rotmat_[11]);
		rmsd+=d*d;
	}
	return(rmsd/nAtom);
}
//dRMS
double Kabsch::calc_dist(XYZ *molA,XYZ *molB,int nAtom)
{
	int i,j;
	int count=0;
	double rmsd=0.0,dist;
	for(i=0;i<nAtom-1;i++)
	{
		for(j=i+1;j<nAtom;j++)
		{
			dist=molA[i].distance(molA[j])-molB[i].distance(molB[j]);
			dist=dist*dist;
			rmsd+=dist;
			count++;
		}
	}
	return(rmsd/count);
}

//rotation_related
void Kabsch::rot_mol(XYZ *m_old,XYZ *m_new,int nAtom,double *rotmat_)
{
	//ws_add
	if(rotmat_==0)
	{
		EqualArray(m_new,m_old,nAtom);
		return;
	}
	//original
	int l;
	double dx,dy,dz;
	for(l=0;l<nAtom;l++) 
	{
		dx=m_old[l].X; 
		dy=m_old[l].Y;
		dz=m_old[l].Z;
		m_new[l].X=dx*rotmat_[0]+dy*rotmat_[1]+dz*rotmat_[2]+rotmat_[9];
		m_new[l].Y=dx*rotmat_[3]+dy*rotmat_[4]+dz*rotmat_[5]+rotmat_[10];
		m_new[l].Z=dx*rotmat_[6]+dy*rotmat_[7]+dz*rotmat_[8]+rotmat_[11];
	}
}
void Kabsch::rot_point(XYZ p_old,XYZ &p_new,double *rotmat_)
{
	//ws_add
	if(rotmat_==0)
	{
		p_new=p_old;
		return;
	}
	//original
	double dx,dy,dz;
	dx=p_old.X;
	dy=p_old.Y;
	dz=p_old.Z;
	p_new.X=dx*rotmat_[0]+dy*rotmat_[1]+dz*rotmat_[2]+rotmat_[9];
	p_new.Y=dx*rotmat_[3]+dy*rotmat_[4]+dz*rotmat_[5]+rotmat_[10];
	p_new.Z=dx*rotmat_[6]+dy*rotmat_[7]+dz*rotmat_[8]+rotmat_[11];
}



//==================== the following is qcprot ===========================//__110720__//
/*******************************************************************************
 *  -/_|:|_|_\- 
 *
 *  File:           qcprot.c
 *  Version:        1.2
 *
 *  Function:       Rapid calculation of the least-squares rotation using a 
 *                  quaternion-based characteristic polynomial and 
 *                  a cofactor matrix
 *
 *  Author(s):      Douglas L. Theobald
 *                  Department of Biochemistry
 *                  MS 009
 *                  Brandeis University
 *                  415 South St
 *                  Waltham, MA  02453
 *                  USA
 *
 *                  dtheobald@brandeis.edu
 *                  
 *                  Pu Liu
 *                  Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
 *                  665 Stockton Drive
 *                  Exton, PA  19341
 *                  USA
 *
 *                  pliu24@its.jnj.com
 * 
 *
 *    If you use this QCP rotation calculation method in a publication, please
 *    reference:
 *
 *      Douglas L. Theobald (2005)
 *      "Rapid calculation of RMSD using a quaternion-based characteristic
 *      polynomial."
 *      Acta Crystallographica A 61(4):478-480.
 *
 *      Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009)
 *      "Fast determination of the optimal rotational matrix for macromolecular 
 *      superpositions."
 *      in press, Journal of Computational Chemistry 
 *
 *
 *  Copyright (c) 2009-2010, Pu Liu and Douglas L. Theobald
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without modification, are permitted 
 *  provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this list of 
 *    conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice, this list 
 *    of conditions and the following disclaimer in the documentation and/or other materials 
 *    provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to 
 *    endorse or promote products derived from this software without specific prior written 
 *    permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 *  Source:         started anew.
 *
 *  Change History:
 *    2009/04/13      Started source
 *    2010/03/28      Modified FastCalcRMSDAndRotation() to handle tiny qsqr
 *                    If trying all rows of the adjoint still gives too small
 *                    qsqr, then just return identity matrix. (DLT)
 *    2010/06/30      Fixed prob in assigning A[9] = 0 in InnerProduct()
 *                    invalid mem access
 *    2011/02/21      Made CenterCoords use weights
 *  
 ******************************************************************************/
double Kabsch::InnerProduct(double *A, XYZ *coords1, XYZ *coords2, XYZ &cen1, XYZ &cen2, int len, double *weight)
{
    int             i;
    double          x1, x2, y1, y2, z1, z2;
    double          G1 = 0.0, G2 = 0.0;
    double          xsum1, ysum1, zsum1, xsum2, ysum2, zsum2, wsum1, wsum2;
    xsum1 = ysum1 = zsum1 = 0.0;
	xsum2 = ysum2 = zsum2 = 0.0;
    A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

	//---- ws_modify ----//
    if (weight != NULL)
    {
	//cent
	wsum1 = 0.0;
	wsum2 = 0.0;
        for (i = 0; i < len; ++i)
        {
			//cent1
			xsum1 += weight[i] * coords1[i].X;
			ysum1 += weight[i] * coords1[i].Y;
			zsum1 += weight[i] * coords1[i].Z;
			wsum1 += weight[i];
			//cent2
			xsum2 += weight[i] * coords2[i].X;
			ysum2 += weight[i] * coords2[i].Y;
			zsum2 += weight[i] * coords2[i].Z;
			wsum2 += weight[i];
        }
        xsum1 /= wsum1;
        ysum1 /= wsum1;
        zsum1 /= wsum1;
        xsum2 /= wsum2;
        ysum2 /= wsum2;
        zsum2 /= wsum2;

		//cross
        for (i = 0; i < len; ++i)
        {
            x1 = coords1[i].X - xsum1;
            y1 = coords1[i].Y - ysum1;
            z1 = coords1[i].Z - zsum1;
            G1 += weight[i] * (x1 * x1 + y1 * y1 + z1 * z1);

            x2 = coords2[i].X - xsum2;
            y2 = coords2[i].Y - ysum2;
            z2 = coords2[i].Z - zsum2;
            G2 += weight[i] * (x2 * x2 + y2 * y2 + z2 * z2);

            A[0] +=  (x1 * x2);
            A[1] +=  (x1 * y2);
            A[2] +=  (x1 * z2);

            A[3] +=  (y1 * x2);
            A[4] +=  (y1 * y2);
            A[5] +=  (y1 * z2);

            A[6] +=  (z1 * x2);
            A[7] +=  (z1 * y2);
            A[8] +=  (z1 * z2);   
        }
    }
    else
    {
		//cent
        for (i = 0; i < len; ++i)
        {
			//cent1
            xsum1 += coords1[i].X;
            ysum1 += coords1[i].Y;
            zsum1 += coords1[i].Z;
			//cent2
            xsum2 += coords2[i].X;
            ysum2 += coords2[i].Y;
            zsum2 += coords2[i].Z;
        }
        xsum1 /= len;
        ysum1 /= len;
        zsum1 /= len;
        xsum2 /= len;
        ysum2 /= len;
        zsum2 /= len;

		//cross
        for (i = 0; i < len; ++i)
        {
            x1 = coords1[i].X - xsum1;
            y1 = coords1[i].Y - ysum1;
            z1 = coords1[i].Z - zsum1;
            G1 += (x1 * x1 + y1 * y1 + z1 * z1);

            x2 = coords2[i].X - xsum2;
            y2 = coords2[i].Y - ysum2;
            z2 = coords2[i].Z - zsum2;
            G2 += (x2 * x2 + y2 * y2 + z2 * z2);

            A[0] +=  (x1 * x2);
            A[1] +=  (x1 * y2);
            A[2] +=  (x1 * z2);

            A[3] +=  (y1 * x2);
            A[4] +=  (y1 * y2);
            A[5] +=  (y1 * z2);

            A[6] +=  (z1 * x2);
            A[7] +=  (z1 * y2);
            A[8] +=  (z1 * z2);  
        }
    }
	cen1.X=xsum1;
	cen1.Y=ysum1;
	cen1.Z=zsum1;
	cen2.X=xsum2;
	cen2.Y=ysum2;
	cen2.Z=zsum2;
    return (G1 + G2) * 0.5;
}
int Kabsch::FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, int len, double minScore)
{
    double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
           SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
           SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
           SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    double C[4];
    int i;
    double mxEigenV; 
    double oldg = 0.0;
    double b, a, delta, rms, qsqr;
    double q1, q2, q3, q4, normq;
    double a11, a12, a13, a14, a21, a22, a23, a24;
    double a31, a32, a33, a34, a41, a42, a43, a44;
    double a2, x2, y2, z2; 
    double xy, az, zx, ay, yz, ax; 
    double a3344_4334, a3244_4234, a3243_4233, a3143_4133,a3144_4134, a3142_4132; 
    double evecprec = 1e-6;
    double evalprec = 1e-14;

    Sxx = A[0]; Sxy = A[1]; Sxz = A[2];
    Syx = A[3]; Syy = A[4]; Syz = A[5];
    Szx = A[6]; Szy = A[7]; Szz = A[8];

    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

    C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz + Szx;
    SyzpSzy = Syz + Szy;
    SxypSyx = Sxy + Syx;
    SyzmSzy = Syz - Szy;
    SxzmSzx = Sxz - Szx;
    SxymSyx = Sxy - Syx;
    SxxpSyy = Sxx + Syy;
    SxxmSyy = Sxx - Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
         + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
         + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
         + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

    mxEigenV = E0;
    for (i = 0; i < 50; ++i)
    {
        oldg = mxEigenV;
        x2 = mxEigenV*mxEigenV;
        b = (x2 + C[2])*mxEigenV;
        a = b + C[1];
        delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a));
        mxEigenV -= delta;
        if (fabs(mxEigenV - oldg) < fabs((evalprec)*mxEigenV))
            break;
    }

	if (i == 50) 
	{
//		fprintf(stderr,"\nMore than %d iterations needed!\n", i);
	}

//  rms = sqrt(2.0 * (E0 - mxEigenV)/len);
	rms = 2.0 * (E0 - mxEigenV)/len; //ws_modified//__110720__//
    (*rmsd) = rms;
	if(rot==0)return (-2); //ws_modified//__110720__//


    if (minScore > 0) 
        if (rms < minScore)
            return (-1); // Don't bother with rotation. 

    a11 = SxxpSyy + Szz-mxEigenV; a12 = SyzmSzy; a13 = - SxzmSzx; a14 = SxymSyx;
    a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx;
    a31 = a13; a32 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy;
    a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV;
    a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44-a42*a34;
    a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33;
    a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42-a41*a32;
    q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
    q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
    q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
    q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

/* The following code tries to calculate another column in the adjoint matrix when the norm of the 
   current column is too small.
   Usually this commented block will never be activated.  To be absolutely safe this should be
   uncommented, but it is most likely unnecessary.  
*/
    if (qsqr < evecprec)
    {
        q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
        q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
        q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
        q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
        qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

        if (qsqr < evecprec)
        {
            double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
            double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
            double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

            q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
            q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
            q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
            q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

            if (qsqr < evecprec)
            {
                q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
                qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;
                
                if (qsqr < evecprec)
                {
                    /* if qsqr is still too small, return the identity matrix. */
                    rot[0] = rot[4] = rot[8] = 1.0;
                    rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0;

                    return(0);
                }
            }
        }
    }

    normq = sqrt(qsqr);
    q1 /= normq;
    q2 /= normq;
    q3 /= normq;
    q4 /= normq;

    a2 = q1 * q1;
    x2 = q2 * q2;
    y2 = q3 * q3;
    z2 = q4 * q4;

    xy = q2 * q3;
    az = q1 * q4;
    zx = q4 * q2;
    ay = q1 * q3;
    yz = q3 * q4;
    ax = q1 * q2;

    rot[0] = a2 + x2 - y2 - z2;
    rot[1] = 2 * (xy + az);
    rot[2] = 2 * (zx - ay);
    rot[3] = 2 * (xy - az);
    rot[4] = a2 - x2 + y2 - z2;
    rot[5] = 2 * (yz + ax);
    rot[6] = 2 * (zx + ay);
    rot[7] = 2 * (yz - ax);
    rot[8] = a2 - x2 - y2 + z2;

    return (1);
}
double Kabsch::CalcRMSDRotationalMatrix(XYZ *coords1, XYZ *coords2, int len, double *rot, double *weight)
{
	double A[9], rmsd;
	XYZ cen1,cen2;
	/* calculate the (weighted) inner product of two structures */
	double E0 = InnerProduct(A, coords1, coords2, cen1, cen2, len, weight);
	/* calculate the RMSD & rotational matrix */
	FastCalcRMSDAndRotation(rot, A, &rmsd, E0, len, -1);

	//ws_get_final_rotmat//__110720__//
	int i,j,k;
	if(rot!=0)
	{
		double rr[3][3];
		double y[3];
		double x[3];
		k=0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				rr[i][j]=rot[k];
				k++;
			}
		}
		cen2.xyz2double(x);
		Vector_Multiply(y,rr,x);
		cen2.double2xyz(y);
		rot[9]=-cen2.X+cen1.X;
		rot[10]=-cen2.Y+cen1.Y;
		rot[11]=-cen2.Z+cen1.Z;
	}

	//final return
	return rmsd;
}


