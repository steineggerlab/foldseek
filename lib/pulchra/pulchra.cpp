/*
 * Excerpt from pulchra, does only backbone reconstruction, to be used as libary!
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pulchra_common.h"
#include "nco_data.h"
#include "pulchra.h"

// distance
double calc_distance(double x1, double y1, double z1,
						  		  double x2, double y2, double z2)
{
  double dx,dy,dz;
  double dist2;

    dx = (x1) - (x2);
    dy = (y1) - (y2);
    dz = (z1) - (z2);
    if (dx || dy || dz ) {
      dist2 = dx*dx+dy*dy+dz*dz;
      return (sqrt(dist2));
    } else
      return 0.0;
}

// r14 chiral distance
double calc_r14(double x1, double y1, double z1,
							 double x2, double y2, double z2,
							 double x3, double y3, double z3,
							 double x4, double y4, double z4)
{
  double r, dx, dy, dz;
  double vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3;
  double hand;

    dx = x4-x1;
    dy = y4-y1;
    dz = z4-z1;

    r = sqrt(dx*dx+dy*dy+dz*dz);

    vx1=x2-x1;
    vy1=y2-y1;
    vz1=z2-z1;
    vx2=x3-x2;
    vy2=y3-y2;
    vz2=z3-z2;
    vx3=x4-x3;
    vy3=y4-y3;
    vz3=z4-z3;

    hand = (vy1*vz2-vy2*vz1)*vx3+
           (vz1*vx2-vz2*vx1)*vy3+
           (vx1*vy2-vx2*vy1)*vz3;

    if (hand<0) r=-r;

  return r;
}

// superimposition of two sets for coordinates + optional transformation of tpoints

double superimpose2(double **coords1, double **coords2, int npoints, double **tpoints, int ntpoints)
{
  double mat_s[3][3], mat_a[3][3], mat_b[3][3], mat_g[3][3];
  double mat_u[3][3], tmp_mat[3][3];
  double val, d, alpha, beta, gamma, x, y, z;
  double cx1, cy1, cz1, cx2, cy2, cz2, tmpx, tmpy, tmpz;
  int i, j, k, n;

    cx1=cy1=cz1=cx2=cy2=cz2=0.;

    for (i=0; i<npoints; i++) {
      cx1+=coords1[i][0];
      cy1+=coords1[i][1];
      cz1+=coords1[i][2];
      cx2+=coords2[i][0];
      cy2+=coords2[i][1];
      cz2+=coords2[i][2];
    }

    cx1/=(double)npoints;
    cy1/=(double)npoints;
    cz1/=(double)npoints;

    cx2/=(double)npoints;
    cy2/=(double)npoints;
    cz2/=(double)npoints;

    for (i=0; i<npoints; i++) {
      coords1[i][0]-=cx1;
      coords1[i][1]-=cy1;
      coords1[i][2]-=cz1;
      coords2[i][0]-=cx2;
      coords2[i][1]-=cy2;
      coords2[i][2]-=cz2;
    }

    for (i=0; i<ntpoints; i++) {
      tpoints[i][0]-=cx2;
      tpoints[i][1]-=cy2;
      tpoints[i][2]-=cz2;
    }

    for (i=0; i<3; i++)
      for (j=0; j<3; j++) {
        if (i==j)
          mat_s[i][j]=mat_a[i][j]=mat_b[i][j]=mat_g[i][j]=1.0;
        else
          mat_s[i][j]=mat_a[i][j]=mat_b[i][j]=mat_g[i][j]=0.0;
        mat_u[i][j]=0.;
      }

    for (n=0; n<npoints; n++) {
      mat_u[0][0]+=coords1[n][0]*coords2[n][0];
      mat_u[0][1]+=coords1[n][0]*coords2[n][1];
      mat_u[0][2]+=coords1[n][0]*coords2[n][2];
      mat_u[1][0]+=coords1[n][1]*coords2[n][0];
      mat_u[1][1]+=coords1[n][1]*coords2[n][1];
      mat_u[1][2]+=coords1[n][1]*coords2[n][2];
      mat_u[2][0]+=coords1[n][2]*coords2[n][0];
      mat_u[2][1]+=coords1[n][2]*coords2[n][1];
      mat_u[2][2]+=coords1[n][2]*coords2[n][2];
    }

    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
        tmp_mat[i][j]=0.;

    do {
      d=mat_u[2][1]-mat_u[1][2];
      if (d==0) alpha=0; else alpha=atan(d/(mat_u[1][1]+mat_u[2][2]));
      if (cos(alpha)*(mat_u[1][1]+mat_u[2][2])+sin(alpha)*(mat_u[2][1]-mat_u[1][2])<0.0)       alpha+=M_PI;
      mat_a[1][1]=mat_a[2][2]=cos(alpha);
      mat_a[2][1]=sin(alpha);
      mat_a[1][2]=-mat_a[2][1];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_a[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_a[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      d=mat_u[0][2]-mat_u[2][0];
      if (d==0) beta=0; else beta=atan(d/(mat_u[0][0]+mat_u[2][2]));
      if (cos(beta)*(mat_u[0][0]+mat_u[2][2])+sin(beta)*(mat_u[0][2]-mat_u[2][0])<0.0) beta+=M_PI;
      mat_b[0][0]=mat_b[2][2]=cos(beta);
      mat_b[0][2]=sin(beta);
      mat_b[2][0]=-mat_b[0][2];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_b[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_b[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      d=mat_u[1][0]-mat_u[0][1];
      if (d==0) gamma=0; else gamma=atan(d/(mat_u[0][0]+mat_u[1][1]));
      if (cos(gamma)*(mat_u[0][0]+mat_u[1][1])+sin(gamma)*(mat_u[1][0]-mat_u[0][1])<0.0)
        gamma+=M_PI;
      mat_g[0][0]=mat_g[1][1]=cos(gamma);
      mat_g[1][0]=sin(gamma);
      mat_g[0][1]=-mat_g[1][0];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_g[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_g[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      val=fabs(alpha)+fabs(beta)+fabs(gamma);
    } while (val>0.001);

    val=0.;
    for (i=0; i<npoints; i++) {
      x=coords2[i][0];
      y=coords2[i][1];
      z=coords2[i][2];
      tmpx=x*mat_s[0][0]+y*mat_s[0][1]+z*mat_s[0][2];
      tmpy=x*mat_s[1][0]+y*mat_s[1][1]+z*mat_s[1][2];
      tmpz=x*mat_s[2][0]+y*mat_s[2][1]+z*mat_s[2][2];
      x=coords1[i][0]-tmpx;
      y=coords1[i][1]-tmpy;
      z=coords1[i][2]-tmpz;
      val+=x*x+y*y+z*z;
    }

    for (i=0; i<ntpoints; i++) {
      x=tpoints[i][0];
      y=tpoints[i][1];
      z=tpoints[i][2];
      tpoints[i][0]=x*mat_s[0][0]+y*mat_s[0][1]+z*mat_s[0][2];
      tpoints[i][1]=x*mat_s[1][0]+y*mat_s[1][1]+z*mat_s[1][2];
      tpoints[i][2]=x*mat_s[2][0]+y*mat_s[2][1]+z*mat_s[2][2];
    }

    for (i=0; i<npoints; i++) {
      coords1[i][0]+=cx1;
      coords1[i][1]+=cy1;
      coords1[i][2]+=cz1;
      coords2[i][0]+=cx2;
      coords2[i][1]+=cy2;
      coords2[i][2]+=cz2;
    }

    for (i=0; i<ntpoints; i++) {
      tpoints[i][0]+=cx1;
      tpoints[i][1]+=cy1;
      tpoints[i][2]+=cz1;
    }

  return sqrt(val/(double)npoints);
}

void prepare_rbins(double ** ca, int ** rbins, int len, double ** cacoords, double ** tmpcoords, double ** tmpstat)
{
  int i, j, k, bin13_1, bin13_2, bin14;
  double x1, y1, z1;
  double x2, y2, z2;
  double x3, y3, z3;
  double x4, y4, z4;
  double r13_1, r13_2, r14;


  // rebuild ends...
  for (i=0,j=0;i<5;i++,j++)
    for (k=0;k<3;k++)
      tmpcoords[j][k] = ca[i][k];
  for (i=2,j=0;i<5;i++,j++)
    for (k=0;k<3;k++)
      cacoords[j][k] = ca[i][k];
  for (i=0,j=0;i<3;i++,j++)
    for (k=0;k<3;k++)
      tmpstat[j][k] = ca[i][k];

  superimpose2(tmpstat,cacoords,3,tmpcoords,5);

  for (i=-2,j=0;i<0;i++,j++)
    for (k=0;k<3;k++)
      ca[i][k] = tmpcoords[j][k];

  for (i=len-5,j=0;i<len;i++,j++)
    for (k=0;k<3;k++)
      tmpcoords[j][k] = ca[i][k];
  for (i=len-5,j=0;i<len-2;i++,j++)
    for (k=0;k<3;k++)
      cacoords[j][k] = ca[i][k];
  for (i=len-3,j=0;i<len;i++,j++)
    for (k=0;k<3;k++)
      tmpstat[j][k] = ca[i][k];

  superimpose2(tmpstat,cacoords,3,tmpcoords,5);

  for (i=len-3,j=0;i<len;i++,j++)
    for (k=0;k<3;k++)
      ca[i+3][k] = tmpcoords[j+3][k];

  for (i=0;i<len+1;i++) {
  	x1 = ca[i-2][0];
  	y1 = ca[i-2][1];
  	z1 = ca[i-2][2];

  	x2 = ca[i-1][0];
  	y2 = ca[i-1][1];
  	z2 = ca[i-1][2];

  	x3 = ca[i][0];
  	y3 = ca[i][1];
  	z3 = ca[i][2];

  	x4 = ca[i+1][0];
  	y4 = ca[i+1][1];
  	z4 = ca[i+1][2];

  	r13_1 = calc_distance(x1, y1, z1, x3, y3, z3);
  	r13_2 = calc_distance(x2, y2, z2, x4, y4, z4);
  	r14 = calc_r14(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

    bin13_1 = (int)((r13_1-4.6)/0.3);
    bin13_2 = (int)((r13_2-4.6)/0.3);
    bin14 = (int)((r14+11.)/0.3);

    if (bin13_1<0) bin13_1=0;
    if (bin13_2<0) bin13_2=0;
    if (bin14<0) bin14=0;
    if (bin13_1>9) bin13_1=9;
    if (bin13_2>9) bin13_2=9;
    if (bin14>73) bin14=73;

    rbins[i][0] = bin13_1;
    rbins[i][1] = bin13_2;
    rbins[i][2] = bin14;
  }
}

// Reconstruct N and C atoms from Ca atoms. Aminoacids only needed due to Prolin.
// Note xca: Ca atom coords., where the five first/last fields are empty!!
void pulchra_rebuild_backbone(double ** xca, double ** n, double ** c, char * aa, int len)
{
  double **cacoords, **tmpcoords, **tmpstat;
  double x1, y1, z1;
  double x2, y2, z2;
  double x3, y3, z3;
  double x4, y4, z4;
  double besthit, hit;
  int bestpos;
  int i, j, k, bin13_1, bin13_2, bin14;

  int ** rbins = (int**)calloc(sizeof(int*)*(len+1),1);
  for (i=0;i<len+1;i++)
    rbins[i] = (int*)calloc(sizeof(int)*3,1);

  cacoords = (double**)calloc(sizeof(double*)*(8),1);
  tmpcoords = (double**)calloc(sizeof(double*)*(8),1);
  tmpstat = (double**)calloc(sizeof(double*)*(8),1);
  for (i=0;i<8;i++) {
    cacoords[i] = (double*)calloc(sizeof(double)*3,1);;
    tmpcoords[i] = (double*)calloc(sizeof(double)*3,1);;
    tmpstat[i] = (double*)calloc(sizeof(double)*3,1);;
  }

  double ** ca = &xca[5];
  prepare_rbins(ca, rbins, len, cacoords, tmpcoords, tmpstat);

 for (i=0;i<len+1;i++) {
 	x1 = ca[i-2][0];
 	y1 = ca[i-2][1];
 	z1 = ca[i-2][2];

 	x2 = ca[i-1][0];
 	y2 = ca[i-1][1];
 	z2 = ca[i-1][2];

 	x3 = ca[i][0];
 	y3 = ca[i][1];
 	z3 = ca[i][2];

 	x4 = ca[i+1][0];
 	y4 = ca[i+1][1];
 	z4 = ca[i+1][2];

 	cacoords[0][0] = x1;
 	cacoords[0][1] = y1;
 	cacoords[0][2] = z1;

 	cacoords[1][0] = x2;
  cacoords[1][1] = y2;
  cacoords[1][2] = z2;

  cacoords[2][0] = x3;
  cacoords[2][1] = y3;
  cacoords[2][2] = z3;

  cacoords[3][0] = x4;
  cacoords[3][1] = y4;
  cacoords[3][2] = z4;

  bin13_1 = rbins[i][0];
  bin13_2 = rbins[i][1];
  bin14 = rbins[i][2];

      if (i>0 && aa[i-1]=='P') {
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat_pro[j].bins[0]-bin13_1)+fabs(nco_stat_pro[j].bins[1]-bin13_2)+0.2*fabs(nco_stat_pro[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat_pro[j].bins[0]>=0 && hit>1e-3);
        for (j=0;j<4;j++) {
         	for (k=0;k<3;k++) {
         		tmpstat[j][k] = nco_stat_pro[bestpos].data[j][k];
       	  }
     	  }
        for (j=0;j<8;j++) {
         	for (k=0;k<3;k++) {
         		tmpcoords[j][k] = nco_stat_pro[bestpos].data[j][k];
       	  }
     	  }
      } else {
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat[j].bins[0]-bin13_1)+fabs(nco_stat[j].bins[1]-bin13_2)+0.2*fabs(nco_stat[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat[j].bins[0]>=0 && hit>1e-3);
        for (j=0;j<4;j++) {
         	for (k=0;k<3;k++) {
         		tmpstat[j][k] = nco_stat[bestpos].data[j][k];
       	  }
     	  }
        for (j=0;j<8;j++) {
         	for (k=0;k<3;k++) {
         		tmpcoords[j][k] = nco_stat[bestpos].data[j][k];
       	  }
     	  }
      }

     	superimpose2(cacoords, tmpstat, 4, tmpcoords, 8);

      if (i>0) {
        c[i-1][0] = tmpcoords[4][0];
        c[i-1][1] = tmpcoords[4][1];
        c[i-1][2] = tmpcoords[4][2];
      }

      if (i<len) {
        n[i][0] = tmpcoords[6][0];
        n[i][1] = tmpcoords[6][1];
        n[i][2] = tmpcoords[6][2];
      }
 }
 
  for (i=0;i<8;i++) {
    free(cacoords[i]);
    free(tmpcoords[i]);
    free(tmpstat[i]);
  }
  free(cacoords);
  free(tmpcoords);
  free(tmpstat);

  for (i=0;i<len+1;i++) {
    free(rbins[i]);
  }
  free(rbins);
}

