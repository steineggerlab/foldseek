#pragma once
#include "TM_align.h"
using namespace std;

//==== class: TM_align_Main =====//
class TM_Align_Main : virtual public TM_align
{
public:
	TM_Align_Main(int num=PROT_MAX_NUM);
	~TM_Align_Main(void);
	int TMM_maximal;

//---- variables ----//
public:
	int TM_ADDITION;   //macro
	double *TM_dist;   //record temp dist
	int *TM_sse1;      //TM_based sse
	int *TM_sse2;      //...
	int *TM_frag1;     //2*n, [0]->record start pos; [1]->record len
	int *TM_frag2;     //...
	int *TM_Main_Ali2; //temporary for best ali2 record
	int *TM_Main_AliT; //temporary for best ali2 record (temp)
	XYZ *TM_ttt1;      //for temporary
	XYZ *TM_ttt2;      //for temporary
	int AliT_Rec;      //AliT got or not //__110408__//
	int AliT_Max;      //record the best method //__110408__//
	double *TM_FINMAT; //best matrix
	double TM_FIN_TMS; //best TMscore
	double TM_FIN_RMS; //best RMSD
	int TM_FIN_LALI;   //best lali


//---- functions ----//
public:
	//create & delete
	void TM_Align_Main_Init(int maxlen);
	void TM_Align_Main_Dele(void);
	//minor functions
	double TM_Align_Get_GL(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2);
	int TM_Assign_SSE_Single(double d13,double d14,double d15,double d24,double d25,double d35);
	void TM_Assign_SSE(int *sse,XYZ *mol,int moln);
	void TM_Smooth_SSE(int *sse,int moln);
	//major functions
	void TM_Get_Initial1(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2); //gapless threading
	void TM_Get_Initial2(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2); //SSE dynamic programming
	void TM_Get_Initial3(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2); //using previous and SSE alignment
	void TM_Get_Initial4(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2); //no-chain-break gapless threading
	void TM_Get_Initial5(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2); //local structure superposition (neo version)
	void TM_Get_Initial6(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2); //local structure superposition (old version)
	//main function
	//[1] Zhang Yang original TMalign version (2009+2010)
	double TM_Align_Total(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,
		int norm_len,double norm_d0,double *MAXSCO=0); //MAXSCO should at least double[8]
	//[2] modified TMalign, only change the order of initial alignment
	double TM_Align_Total_II(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,
		int norm_len,double norm_d0,double *MAXSCO=0);
};

/*
****************************************************************
*     quickly calculate TM-score with given invmap(i) in 3 iterations
****************************************************************
      subroutine get_GL(GL)
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
      common/d00/d00,d002

      dimension xo1(nmax),yo1(nmax),zo1(nmax)
      dimension xo2(nmax),yo2(nmax),zo2(nmax)
      dimension dis2(nmax)

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
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
            xo1(n_al)=xa(1,i,0)
            yo1(n_al)=xa(2,i,0)
            zo1(n_al)=xa(3,i,0)
            xo2(n_al)=xa(1,j,1)
            yo2(n_al)=xa(2,j,1)
            zo2(n_al)=xa(3,j,1)
         endif
      enddo
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      GL=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         GL=GL+1/(1+dis2(i)/(d0**2))
      enddo
ccc   for next iteration------------->
      d002t=d002
 21   j=0
      do i=1,n_al
         if(dis2(i).le.d002t)then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
         endif
      enddo
      if(j.lt.3.and.n_al.gt.3)then
         d002t=d002t+.5
         goto 21
      endif
      L=j
      call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2
      G2=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         G2=G2+1/(1+dis2(i)/(d0**2))
      enddo
ccc   for next iteration------------->
      d002t=d002+1
 22   j=0
      do i=1,n_al
         if(dis2(i).le.d002t)then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
         endif
      enddo
      if(j.lt.3.and.n_al.gt.3)then
         d002t=d002t+.5
         goto 22
      endif
      L=j
      call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2
      G3=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         G3=G3+1/(1+dis2(i)/(d0**2))
      enddo
      if(G2.gt.GL)GL=G2
      if(G3.gt.GL)GL=G3

c^^^^^^^^^^^^^^^^ GL done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
*/



/*
**************************************************************
*     smooth the secondary structure assignment
**************************************************************
      subroutine smooth
      PARAMETER(nmax=5000)
      common/sec/isec(nmax),jsec(nmax)
      common/length/nseq1,nseq2

***   smooth single -------------->
***   --x-- => -----
      do i=1,nseq1
         if(isec(i).eq.2.or.isec(i).eq.4)then
            j=isec(i)
            if(isec(i-2).ne.j)then
               if(isec(i-1).ne.j)then
                  if(isec(i+1).ne.j)then
                     if(isec(i+1).ne.j)then
                        isec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).eq.2.or.jsec(i).eq.4)then
            j=jsec(i)
            if(jsec(i-2).ne.j)then
               if(jsec(i-1).ne.j)then
                  if(jsec(i+1).ne.j)then
                     if(jsec(i+1).ne.j)then
                        jsec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo

***   smooth double -------------->
***   --xx-- => ------
      do i=1,nseq1
         if(isec(i).ne.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
         if(isec(i+3).eq.2)then
         if(isec(i+4).ne.2)then
         if(isec(i+5).ne.2)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif

         if(isec(i).ne.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
         if(isec(i+3).eq.4)then
         if(isec(i+4).ne.4)then
         if(isec(i+5).ne.4)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).ne.2)then
         if(jsec(i+1).ne.2)then
         if(jsec(i+2).eq.2)then
         if(jsec(i+3).eq.2)then
         if(jsec(i+4).ne.2)then
         if(jsec(i+5).ne.2)then
            jsec(i+2)=1
            jsec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif

         if(jsec(i).ne.4)then
         if(jsec(i+1).ne.4)then
         if(jsec(i+2).eq.4)then
         if(jsec(i+3).eq.4)then
         if(jsec(i+4).ne.4)then
         if(jsec(i+5).ne.4)then
            jsec(i+2)=1
            jsec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif
      enddo

***   connect -------------->
***   x-x => xxx
      do i=1,nseq1
         if(isec(i).eq.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
            isec(i+1)=2
         endif
         endif
         endif

         if(isec(i).eq.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
            isec(i+1)=4
         endif
         endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).eq.2)then
         if(jsec(i+1).ne.2)then
         if(jsec(i+2).eq.2)then
            jsec(i+1)=2
         endif
         endif
         endif

         if(jsec(i).eq.4)then
         if(jsec(i+1).ne.4)then
         if(jsec(i+2).eq.4)then
            jsec(i+1)=4
         endif
         endif
         endif
      enddo

      return
      end

*************************************************************
*     assign secondary structure:
*************************************************************
      function make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
      make_sec=1
      delta=2.1
      if(abs(dis15-6.37).lt.delta)then
         if(abs(dis14-5.18).lt.delta)then
            if(abs(dis25-5.18).lt.delta)then
               if(abs(dis13-5.45).lt.delta)then
                  if(abs(dis24-5.45).lt.delta)then
                     if(abs(dis35-5.45).lt.delta)then
                        make_sec=2 !helix
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      delta=1.42
      if(abs(dis15-13).lt.delta)then
         if(abs(dis14-10.4).lt.delta)then
            if(abs(dis25-10.4).lt.delta)then
               if(abs(dis13-6.1).lt.delta)then
                  if(abs(dis24-6.1).lt.delta)then
                     if(abs(dis35-6.1).lt.delta)then
                        make_sec=4 !strand
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      if(dis15.lt.8)then
         make_sec=3
      endif

      return
      end
*/


