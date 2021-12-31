// Structural superposition, the QCP method.
//
// This is modified code from qcprot.c from https://theobald.brandeis.edu/qcp/
// The original qcprot.c file startsd with the following info and license:

/*******************************************************************************
 *  -/_|:|_|_\-
 *
 *  File:           qcprot.c
 *  Version:        1.5
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
 *      Journal of Computational Chemistry 31(7):1561-1563.
 *
 *  Copyright (c) 2009-2016 Pu Liu and Douglas L. Theobald
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
 */

#ifndef GEMMI_QCP_HPP_
#define GEMMI_QCP_HPP_

#include <cmath>         // for fabs, sqrt
#include <cstdio>        // for fprintf (it's temporary)
#include "math.hpp"      // for Mat33
#include "unitcell.hpp"  // for Position

namespace gemmi {

struct SupResult {
  double rmsd;
  size_t count;
  Position center1, center2;
  Transform transform;
};

// helper function
inline double qcp_inner_product(Mat33& mat,
                                const Position* pos1, const Position& ctr1,
                                const Position* pos2, const Position& ctr2,
                                size_t len, const double* weight) {
  double G1 = 0.0, G2 = 0.0;
  for (size_t i = 0; i < len; ++i) {
    Position f1 = pos1[i] - ctr1;
    Position f2 = pos2[i] - ctr2;
    double w = (weight != nullptr ? weight[i] : 1.);
    Vec3 v1 = w * f1;
    G1 += v1.dot(f1);
    G2 += w * f2.length_sq();
    mat[0][0] += v1.x * f2.x;
    mat[0][1] += v1.x * f2.y;
    mat[0][2] += v1.x * f2.z;
    mat[1][0] += v1.y * f2.x;
    mat[1][1] += v1.y * f2.y;
    mat[1][2] += v1.y * f2.z;
    mat[2][0] += v1.z * f2.x;
    mat[2][1] += v1.z * f2.y;
    mat[2][2] += v1.z * f2.z;
  }
  return (G1 + G2) * 0.5;
}

// helper function
inline int fast_calc_rmsd_and_rotation(Mat33* rot, const Mat33& A, double *rmsd,
                                       double E0, double len, double min_score) {
  const double evecprec = 1e-6;
  const double evalprec = 1e-11;

  double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
  Sxx = A[0][0]; Sxy = A[0][1]; Sxz = A[0][2];
  Syx = A[1][0]; Syy = A[1][1]; Syz = A[1][2];
  Szx = A[2][0]; Szy = A[2][1]; Szz = A[2][2];

  double Sxx2 = Sxx * Sxx;
  double Syy2 = Syy * Syy;
  double Szz2 = Szz * Szz;

  double Sxy2 = Sxy * Sxy;
  double Syz2 = Syz * Syz;
  double Sxz2 = Sxz * Sxz;

  double Syx2 = Syx * Syx;
  double Szy2 = Szy * Szy;
  double Szx2 = Szx * Szx;

  double SyzSzymSyySzz2 = 2.0 * (Syz*Szy - Syy*Szz);
  double Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

  double C[4];
  C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
  C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

  double SxzpSzx = Sxz + Szx;
  double SyzpSzy = Syz + Szy;
  double SxypSyx = Sxy + Syx;
  double SyzmSzy = Syz - Szy;
  double SxzmSzx = Sxz - Szx;
  double SxymSyx = Sxy - Syx;
  double SxxpSyy = Sxx + Syy;
  double SxxmSyy = Sxx - Syy;
  double Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

  C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
    + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
    + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
    + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
    + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
    + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

  /* Newton-Raphson */
  double mxEigenV = E0;
  int i;
  double oldg = 0.0;
  for (i = 0; i < 50; ++i) {
    oldg = mxEigenV;
    double x2 = mxEigenV * mxEigenV;
    double b = (x2 + C[2]) * mxEigenV;
    double a = b + C[1];
    double delta = (a * mxEigenV + C[0]) / (2.0 * x2 * mxEigenV + b + a);
    mxEigenV -= delta;
    // printf("\n diff[%3d]: %16g %16g %16g", i, mxEigenV - oldg, evalprec*mxEigenV, mxEigenV);
    if (std::fabs(mxEigenV - oldg) < std::fabs(evalprec * mxEigenV))
      break;
  }
  if (i == 50)
    std::fprintf(stderr,"\nMore than %d iterations needed!\n", i);

  // the fabs() is to guard against extremely small, but *negative* numbers
  // due to floating point error
  double rms = std::sqrt(std::fabs(2.0 * (E0 - mxEigenV) / len));
  (*rmsd) = rms;
  // printf("\n\n %16g %16g %16g \n", rms, E0, 2.0 * (E0 - mxEigenV)/len);

  if (rot == nullptr)
    return -1;
  if (min_score > 0)
    if (rms < min_score)
      return -1; // Don't bother with rotation.

  double a11, a12, a13, a14, a21, a22, a23, a24;
  a11 = SxxpSyy + Szz-mxEigenV; a12 = SyzmSzy; a13 = - SxzmSzx; a14 = SxymSyx;
  a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx;
  double a31, a32, a33, a34, a41, a42, a43, a44;
  a31 = a13; a32 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy;
  a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV;
  double a3344_4334, a3244_4234, a3243_4233, a3143_4133,a3144_4134, a3142_4132;
  a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44-a42*a34;
  a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33;
  a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42-a41*a32;
  double q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
  double q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
  double q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
  double q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

  double qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

  /* The following code tries to calculate another column in the adjoint matrix
     when the norm of the current column is too small.
     Usually this block will never be activated.  To be absolutely safe this should be
     uncommented, but it is most likely unnecessary.
     */
  if (qsqr < evecprec) {
    q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
    q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
    q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
    q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
    qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

    if (qsqr < evecprec) {
      double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
      double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
      double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

      q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
      q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
      q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
      q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
      qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

      if (qsqr < evecprec) {
        q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
        q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
        q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
        q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
        qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;

        if (qsqr < evecprec) {
          /* if qsqr is still too small, return the identity matrix. */
          *rot = Mat33();
          return 0;
        }
      }
    }
  }

  double normq = std::sqrt(qsqr);
  q1 /= normq;
  q2 /= normq;
  q3 /= normq;
  q4 /= normq;

  double a2 = q1 * q1;
  double x2 = q2 * q2;
  double y2 = q3 * q3;
  double z2 = q4 * q4;

  double xy = q2 * q3;
  double az = q1 * q4;
  double zx = q4 * q2;
  double ay = q1 * q3;
  double yz = q3 * q4;
  double ax = q1 * q2;

  rot->a[0][0] = a2 + x2 - y2 - z2;
  rot->a[0][1] = 2 * (xy + az);
  rot->a[0][2] = 2 * (zx - ay);
  rot->a[1][0] = 2 * (xy - az);
  rot->a[1][1] = a2 - x2 + y2 - z2;
  rot->a[1][2] = 2 * (yz + ax);
  rot->a[2][0] = 2 * (zx + ay);
  rot->a[2][1] = 2 * (yz - ax);
  rot->a[2][2] = a2 - x2 - y2 + z2;

  return 1;
}

// helper function
inline Position qcp_calculate_center(const Position* pos, size_t len, const double *weight) {
  double wsum = 0.0;
  Position ctr;
  for (size_t i = 0; i < len; ++i) {
    double w = (weight != nullptr ? weight[i] : 1.);
    ctr += w * pos[i];
    wsum += w;
  }
  return ctr / wsum;
}

// Calculate superposition of pos2 onto pos1 -- pos2 is movable.
// Does not perform the superposition, only returns the operation to be used.
inline SupResult superpose_positions(const Position* pos1, const Position* pos2,
                                     size_t len, const double* weight) {
  SupResult result;
  result.count = len;

  /* center the structures -- if precentered you can omit this step */
  result.center1 = qcp_calculate_center(pos1, len, weight);
  result.center2 = qcp_calculate_center(pos2, len, weight);

  double wsum = 0.0;
  if (weight == nullptr)
    wsum = (double) len;
  else
    for (size_t i = 0; i < len; ++i)
      wsum += weight[i];

  Mat33 A(0);
  /* calculate the (weighted) inner product of two structures */
  double E0 = qcp_inner_product(A, pos1, result.center1, pos2, result.center2, len, weight);

  /* calculate the RMSD & rotational matrix */
  fast_calc_rmsd_and_rotation(&result.transform.mat, A, &result.rmsd, E0, wsum, -1);
  result.transform.vec = Vec3(result.center1) - result.transform.mat.multiply(result.center2);

  return result;
}

// Similar to superpose_positions(), but calculates RMSD only.
inline double calculate_rmsd_of_superposed_positions(const Position* pos1,
                                                     const Position* pos2,
                                                     size_t len, const double* weight) {
  double result;

  // center the structures
  Position ctr1 = qcp_calculate_center(pos1, len, weight);
  Position ctr2 = qcp_calculate_center(pos2, len, weight);

  double wsum = 0.0;
  if (weight == nullptr)
    wsum = (double) len;
  else
    for (size_t i = 0; i < len; ++i)
      wsum += weight[i];

  // calculate the (weighted) inner product of two structures
  Mat33 A(0);
  double E0 = qcp_inner_product(A, pos1, ctr1, pos2, ctr2, len, weight);

  // calculate the RMSD
  fast_calc_rmsd_and_rotation(nullptr, A, &result, E0, wsum, -1);
  return result;
}

} // namespace gemmi
#endif
