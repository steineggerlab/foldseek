#pragma once
#include <assert.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
using namespace std;

//--- utility ---//
#define DYNA_PROG_MAXIMAL 25000000   //5000*5000
#define DYNA_PROG_LENGTH  10000      //5000*2
extern int DP_align1[DYNA_PROG_LENGTH];
extern int DP_align2[DYNA_PROG_LENGTH];

//======= Dynamic Programming =====//
extern int Normal_Align_Dyna_Prog(int n1,int n2,double *score,double GAP_OPEN,double GAP_HEAD,
								  vector<pair<int,int> > & alignment,double &ali_sco);     //normal non-extend DynaProg
extern int Normal_Align_Dyna_Prog_II(int n1,int n2,double *score,double GAP_OPEN,double GAP_HEAD,
								     vector<pair<int,int> > & alignment,double &ali_sco);  //modified non-extend DynaProg (for TM_align)
extern int Advance_Align_Dyna_Prog(int n1,int n2,double *score,double GAP_OPEN,double GAP_EXT,
								   vector<pair<int,int> > & alignment,double &ali_sco);    //three-state DynaProg
extern int Advance_Align_Dyna_Prog_Double(int n1,int n2,double *score,
								   double GAP_OPEN1,double GAP_EXT1,double GAP_OPEN2,double GAP_EXT2,
								   double GAP_HEAD1,double GAP_TAIL1,double GAP_HEAD2,double GAP_TAIL2,
								   vector<pair<int,int> > & alignment,double &ali_sco);    //three-state double DynaProg
extern int Advance_Align_Dyna_Prog_II(int n1,int n2,double *score,double GAP_OPEN,double GAP_EXT,
									  vector<pair<int,int> > & alignment,double &ali_sco); //four-state DynaProg (for TM_align)

//---------- ws_additional dyna_prog -----------//
extern int Normal_Align_Dyna_Prog_Fast(int n1,int n2,double *score,double GAP_OPEN,double GAP_HEAD,
								vector<pair<int,int> > &bound,vector<pair<int,int> > &alignment,double &ali_sco);
