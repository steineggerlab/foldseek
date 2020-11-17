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
#include <vector>
#include "Ali_Ori.h"
#include "Fast_Sort.h"



//====class: Ali_AFP====//
class Ali_AFP : virtual public Ali_Ori
{
public:
	Ali_AFP(int num=PROT_MAX_NUM);
	~Ali_AFP(void);
	int Ali_AFP_Maximal;

//----- create & delete ----//
public:
	void AFP_Create_Linear(int IN_MAX_DIM);
	void AFP_Delete_Linear(void);

//--------------------------------//
// DataStructure for AFP_Cor  ---//
//------------------------------//
//AFP_Cor (data_structure)//
//  AFP_Cor[0*4+0] -> totnum
//  AFP_Cor[i*4..] -> the i-th fragment
//  AFP_Cor[i*4+0] -> RMSD of the frag + type ('+' or normal)
//  AFP_Cor[i*4+1] -> start position in mol1 (ii)
//  AFP_Cor[i*4+2] -> start position in mol2 (jj)
//  AFP_Cor[i*4+3] -> fragment length (winlen)

//===== AFP_Cor function_related =====//
public:
	//---- AFP_Cor_kabsch ----//
	double AFP_kabsch(double *rotmat_,double &rmsd,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *AFP_Cor);                                               // AFP_Cor_rotation (complex)
	//--- update_correspondence ---//
	void update_corespond(int *AFP_Cor,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2);                    // update_correspondence (complex)
	void update_corespond_II(int *AFP_Cor,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2);                    // update_correspondence_II (complex)
	//---- AFP_Cor -> ali1&2 ---//
	void AFP_Cor_Add(int *AFP_Cor,double *rotmat_,int method,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali1,int *ali2);                                        // AFP_Cor -> ali1&2 (complex)
	//---- ali1&2 -> AFP_Cor ---//
	int Ali_To_Cor(int *AFP_Cor,int thres,
		int moln1,int moln2,int *ali1,int *ali2);                    // ali1&2 -> AFP_Cor (complex)


//---------------------------------------//
// DataStructure for Frag_DynaProg  ----//
//-------------------------------------//
private:
//AFP_Cor Dynamic_Programming
//------Kill_Irrelevant_Cor(ori)
//ori_cor[0][0] totnum
//ori_cor[0][1] realnum
//ori_cor[k][0] score
//ori_cor[k][1] head
//ori_cor[k][2] len
//ori_cor[k][3] winlen  <+,->
//ori_cor[k][4] correspondence[ii]  // because only mol1 will be nonlinear
	int *dp_ori_cor;  //contains cor at the initial process

//------Kill_Irrelevant_Cor(cur)
//cur_cor[0][0] totnum
//cur_cor[0][1] totlen
//cur_cor[k][0] score
//cur_cor[k][1] correspondence[i]
//cur_cor[k][2] correspondence[j]
//cur_cor[k][3] winlen
//cur_cor[k][4] ori_cores
	int *dp_cur_cor;  //contains cor for dynamicprogramming

	//---temp_structure---//
	vector <int> AFP_matrix;
	vector <int> AFP_prerec;
	int *AFP_temp;
	int *AFP_indx;
	int AFP_max_sco;
	int AFP_max_i,AFP_max_j;

public:
	//---- parameter for dynamic_programming ----//
	int AFP_GAP_CUT_1;   // gap cutoff between two fragment(combine) [solid]
	int AFP_GAP_CUT_2;   // gap cutoff between two fragment(divide)  [changeable]
	int AFP_NOL_CUT;     // the nonlinear fragment cutoff            [parameter]


//====//additional class (fast_sort)
public:
	Fast_Sort<int> fast_sort;

//---- AFP_Cor Dynamic Programming ---//
public:
	void AFP_Kill_NonLinear(int *ori_cor,int *bak_cor,double *rot_mat);
	void AFP_Combine_Cor(int *AFP_Cor,int cutoff);
	void AFP_Create_Cur_Cor(int *AFP_Cor);
	void AFP_Linear_Process(int cutoff);
	int AFP_Trace_Back(int *AFP_Cor,int *AFP_Cor_temp,int cutoff);
};
//====class: Ali_AFP====//over
