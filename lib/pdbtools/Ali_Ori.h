///////////////////////////////////////////////////////////////////////////
//  CLePAPS:                                                             //
//  Conformational Letter based Pairwise Alignment of Protein Structure  //
//                                                                       //
//  Authors: S.Wang & WM.Zheng                                           //
///////////////////////////////////////////////////////////////////////////
/*
 
Copyright (c)  2007-2008  Institute of Theoretical Physics,Academia Sinica
All Rights Reserved
 
Permission to use, copy, modify and distribute any part of this CLePAPS
software for educational, research and non-profit purposes, without fee,
and without a written agreement is hereby granted, provided that the above
copyright notice, this paragraph and the following three paragraphs appear
in all copies.
 
Those desiring to incorporate this CLePAPS Software into commercial products
or use for commercial purposes should contact Institute of Theoretical Physics,
Academia Sinica,Beijing 100190, China. 
EMAIL: zheng@itp.ac.cn, wangsheng@itp.ac.cn
 
IN NO EVENT SHALL ITP(INSTITUTE OF THEORETICAL PHYSICS) BE LIABLE TO ANY PARTY 
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS CLEPAPS SOFTWARE, EVEN IF ITP HAS 
BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
THE CLEPAPS SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND ITP HAS NO 
OBLIGATION TO PROVIDE MAINTENANCE,SUPPORT, UPDATES, ENHANCEMENTS, OR 
MODIFICATIONS.  ITP MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF 
ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT
THE USE OF THE CLEPAPS SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR
OTHER RIGHTS.
*/

#pragma once
#include "Kabsch.h"



//==Definition==//
//--distance_boundary--//
#ifndef MAX_DIST
#define MAX_DIST 15.0
#endif
#ifndef MIN_DIST
#define MIN_DIST 5.0
#endif
//==Definition==//over




//====class: Ali_Ori====//
//=> original one-to-one correspondence (ali1 & ali2)
class Ali_Ori:virtual public Kabsch
{
public:
	Ali_Ori(int num=PROT_MAX_NUM);
	~Ali_Ori(void);
	int Ali_Maximal;

//----- data_structure ----//
public:
	int *ali1,*ali2;
	int *mas1,*mas2; //region mask 
	XYZ *mol1,*mol2;
	int moln1,moln2;
	int ALI_CACHE;   //do cache
	XYZ *Ali_molt;   //cache point for mol1
	int *Ali_cache;  //cache index for mol1
	XYZ *AFP_tmp1;
	XYZ *AFP_tmp2;

//====== function =======//
public:
	//--cache //
	void Ali_Cache_Point(XYZ *mol,int index,XYZ &ret_point,double *rotmat);
	//--create & init--//
	void Ali_Init(int moln1,int moln2,int *ali1,int *ali2);         // ali1 & ali2 init (complex)
	//--generate aliI from aliJ --//
	void Ali1_To_Ali2(int moln1,int moln2,int *ali1,int *ali2);     // from ali1 -> ali2 (complex)
	void Ali2_To_Ali1(int moln1,int moln2,int *ali1,int *ali2);     // from ali2 -> ali1 (complex)
	//--final_rot_ali --//
	double Final_Rot_Ali(double *rotmat_,double &rmsd,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali1,int *ali2);                                       // ali_ori rotation (complex)
	//--make_center--//
	double make_center(int ii,int jj,int len,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2);                   // use a fragment to pivot (complex)
	//--main usage --//
	double calc_frag_dist(int ii,int jj,int len,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2);                   // calculate RMSD (complex)
	double calc_frag_dist_thres(int ii,int jj,int len,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		double thres,int &ii_,int &jj_,int &len_,double &score);    // calculate RMSD [thres] (complex)
};
//====class_Ali_Ori====//over
