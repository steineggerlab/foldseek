#pragma once
#include <vector>
#include "Dynamic_Programming.h"
#include "TM_score.h"
using namespace std;

//==== class: TM_align =====//
class TM_align : public TM_score
{
public:
	TM_align(int num=PROT_MAX_NUM);
	~TM_align(void);
	int TM_maximal;

//---- parameter ----//
public:
	int TM_NOTM;          //[0,1] -> use TMscore to superimpose or not
	int TM_GAP_TYPE;      //[0,1] -> type
	int TM_GAP_STAGE;     //[1,2] -> using which stage
	double TM_GAP_OPEN;   //TM_GAP_OPEN [-0.6]
	double TM_GAP_EXTEND; //TM_GAP_EXTEND [0.0]
	int TM_DIST_CUT;      //cut those dist > d8 (default: CUT)
	double TM_DistCut;    //TM_DistCut value, default=d8
	int TM_Score_Type;    //[0,1] -> score type (0 for TMscore, 1 for DeepScore)
//---- variables ----//
public:
	double *TM_DP_sco;  //DynaProg matrix
	vector<pair<int,int> > TM_dynaprog; //DynaProg alignment //__110720__//
	vector<pair<int,int> > TM_bound;    //DynaProg bound //__110720__//
	int *TM_ali1;                       //DynaProg Ali3 //__110720__//
	int TM_bound_neib;   //default (4)  //DynaProg bound neibor //__110720__//
	XYZ *TM_tmp1;        //TM_align temp1
	XYZ *TM_tmp2;        //TM_align temp2
	int *TM_DP_ali2;     //DynaProg traceback path (temp)
	int *TM_DP_best;     //DynaProg traceback path (best)
	double *TM_rotmat;   //best rotmat
	//--weight--//__110820__//
	int TM_Vect_Score;      //default is COMBINE,
	int TM_Wei_Score;       //default is NO,
	double *TMali_Weight;   //default is NULL, must be n1*n2!!
	XYZ *TM_cb1;
	XYZ *TM_cb2;

//---- functions ----//
public:
	//init function
	void Init_TM_align(int maxlen);
	void Dele_TM_align(void);
	//part scoring function //__110930__//
	double TM_Align_Get_Score_Part_Wei(XYZ *mol1,XYZ *mol2,double *rotmat_,int moln1,int moln2,int *ali2);
	double TM_Align_Get_Score_Part_Vect(XYZ *mol1,XYZ *mol2,double *rotmat_,int moln1,int moln2,int *ali2);
	double TM_Align_Get_Score_Part_TMsco(XYZ *mol1,XYZ *mol2,double *rotmat_,int moln1,int moln2,int *ali2);
	double TM_Align_Get_Score_Part_Rmsd(XYZ *mol1,XYZ *mol2,double *rotmat_,int moln1,int moln2,int *ali2);
	//primary function
	void TM_Align_Init(int moln1,int moln2);
	int TM_Align_Get_XYZ(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2);
	int TM_Align_Get_CUT(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,double d8,double *rotmat);
	double TM_Align_Get_Score_Simp(XYZ *mol1,XYZ *mol2,double *rotmat_,int moln1,int moln2,int *ali2);
	double TM_Align_Get_Score_Simp_MatchWei(XYZ *mol1,XYZ *mol2,double *rotmat_,
		int moln1,int moln2,int *ali2,vector <double> & MatchWei);
	//minor function
	double TM_Vector_Score_Forward(int ii,int jj,int range,double *rotmat_,XYZ *mol1,XYZ *mol2,int moln1,int moln2);
	double TM_Vector_Score_Backward(int ii,int jj,int range,double *rotmat_,XYZ *mol1,XYZ *mol2,int moln1,int moln2);
	double TM_Vector_Score_CB(int ii,int jj,double *rotmat_,XYZ *mol1,XYZ *mol2,int moln1,int moln2);
	void TM_Align_Get_Matrix(XYZ *mol1,XYZ *mol2,int moln1,int moln2,double d0,double *score);
	void TM_Align_Get_Matrix_TMs(XYZ *mol1,XYZ *mol2,int moln1,int moln2,double d0,double *score);
	double TM_Align_Get_Score(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,int NOTM=0);
	void TM_Align_Get_Ali(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2);
	void TM_Align_Get_Matrix_MatchWei(XYZ *mol1,XYZ *mol2,int moln1,int moln2,double d0,
		vector <double> & MatchWei,vector <double> & RetMatrix);
	//main function
	double TM_Align_TM_Score(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,
		int norm_len,double norm_d0,double &rmsd,int &lali,double *MAXSCO=0); //MAXSCO should at least double[8]
	double TM_Align_TM_Score(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,double &rmsd,int &lali,double *MAXSCO=0);
	double TM_Align_TM_Score_Simp(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,
		int norm_len,double norm_d0,double &rmsd,int &lali,double *MAXSCO=0); //MAXSCO should at least double[8]
	double TM_Align_TM_Score_Simp(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,double &rmsd,int &lali,double *MAXSCO=0);
	double TM_Align_TM_Score_Simplest(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,
		int norm_len,double norm_d0,double &rmsd,int &lali,double *MAXSCO=0); //MAXSCO should at least double[8]
	double TM_Align_TM_Score_Simplest(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,double &rmsd,int &lali,double *MAXSCO=0);
	//bound function //__110730__//
	void TM_Align_Get_Bound(int moln1,int moln2,vector<pair<int,int> > &align,
		vector<pair<int,int> > &bound,int neib);
	void TM_Align_Get_Ali2_Bound(int moln1,int moln2,int *ali2,vector<pair<int,int> > &bound,int neib);
	void TM_Align_Get_Matrix_Bound(XYZ *mol1,XYZ *mol2,int moln1,int moln2,double d0,
		vector<pair<int,int> > &bound,double *score);
	//TMalign main
	int TM_Align_Dyna_Prog(int n1,int n2,double *score,int *ali2,
		double gapopen,double gapextend=0.0,int DP_Type=0);
	double Calc_TM_Align(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_in,int *ali2_out,
		int INI_SKIP=0,int ITER_NUM=30,int SKIPorNOT=1);

};

/*
****************************************************************
*     with invmap(i) calculate TM-score and martix score(i,j) for rotation 
****************************************************************
      subroutine get_score
      PARAMETER(nmax=5000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for TM-score:
            xtm1(n_al)=xa(1,i,0) !for TM-score
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
         endif
      enddo
***   calculate TM-score for the given alignment----------->
      d0_input=d0
      call TMscore8_search(d0_input,n_al,xtm1,ytm1,ztm1,n1,
     &     n_al,xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm) !simplified search engine
      TM=TM*n_al/anseq          !TM-score
***   calculate score matrix score(i,j)------------------>
      do i=1,n_al
         r_2(1,i)=xtm1(i)
         r_2(2,i)=ytm1(i)
         r_2(3,i)=ztm1(i)
      enddo
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
            score(i,j)=1/(1+dd/d0**2)
         enddo
      enddo
      
c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
*/


/*
********************************************************************
*     Dynamic programming for alignment.
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
*     
*     Please note this subroutine is not a correct implementation of 
*     the N-W dynamic programming because the score tracks back only 
*     one layer of the matrix. This code was exploited in TM-align 
*     because it is about 1.5 times faster than a complete N-W code
*     and does not influence much the final structure alignment result.
********************************************************************
      SUBROUTINE DP(NSEQ1,NSEQ2)
      PARAMETER(nmax=5000)
      LOGICAL*1 DIR
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      dimension DIR(0:nmax,0:nmax),VAL(0:nmax,0:nmax)
      common/speed/align_len1
      REAL H,V
      align_len1=0
***   initialize the matrix:
      val(0,0)=0
      do i=1,nseq1
        dir(i,0)=.false.
        val(i,0)=0
      enddo
      do j=1,nseq2
        dir(0,j)=.false.
        val(0,j)=0
        invmap(j)=-1
      enddo

***   decide matrix and path:
      DO j=1,NSEQ2
        DO i=1,NSEQ1
          D=VAL(i-1,j-1)+SCORE(i,j)
          H=VAL(i-1,j)
          if(DIR(i-1,j))H=H+GAP_OPEN
          V=VAL(i,j-1)
          if(DIR(i,j-1))V=V+GAP_OPEN
          
          IF((D.GE.H).AND.(D.GE.V)) THEN
            DIR(I,J)=.true.
            VAL(i,j)=D
          ELSE
            DIR(I,J)=.false.
            if(V.GE.H)then
              val(i,j)=v
            else
              val(i,j)=h
            end if
          ENDIF
        ENDDO
      ENDDO
      
***   extract the alignment:
      i=NSEQ1
      j=NSEQ2
      DO WHILE((i.GT.0).AND.(j.GT.0))
        IF(DIR(i,j))THEN
          invmap(j)=i
          align_len1=align_len1+1
          i=i-1
          j=j-1
        ELSE
          H=VAL(i-1,j)
          if(DIR(i-1,j))H=H+GAP_OPEN
          V=VAL(i,j-1)
          if(DIR(i,j-1))V=V+GAP_OPEN
          IF(V.GE.H) THEN
            j=j-1
          ELSE
            i=i-1
          ENDIF
        ENDIF
      ENDDO
     
c^^^^^^^^^^^^^^^Dynamical programming done ^^^^^^^^^^^^^^^^^^^
      return
      END
*/
