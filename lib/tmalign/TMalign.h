/*
=============================================================
   Implementation of TM-align in C/C++

   This program is written by Jianyi Yang at
   Yang Zhang lab
   And it is updated by Jianjie Wu at
   Yang Zhang lab
   Department of Computational Medicine and Bioinformatics
   University of Michigan
   100 Washtenaw Avenue, Ann Arbor, MI 48109-2218

   Please report bugs and questions to zhng@umich.edu
=============================================================
*/

#include "param_set.h"
#include "NW.h"
#include "Kabsch.h"

float LG_score_soa_sse (float r[16],int nat, float *x1,float *y1,float *z1, float *x2, float *y2, float *z2, float *d,float invd0d0){//coords organized x4,y4,z4
    //no hadds - compatible with SSE2
    float fsum=0;
    __m128 r0 = _mm_load_ps(r);
    __m128 r1 = _mm_load_ps(&(r[4]));
    __m128 r2 = _mm_load_ps(&(r[8]));
    __m128 one=_mm_set_ps(1.0f,1.0f,1.0f,1.0f);
    //4th multiplication unecessary as it is all zeros
    __m128 d0 = _mm_load1_ps(&invd0d0);
    __m128 sum=_mm_setzero_ps();
    int i=0, lower_nat4=(nat/4)*4;
    //four points at a time - otherwise do it scalar
    for( ;i <lower_nat4; i+=4){
        __m128 mx1 = _mm_load_ps(&(x1[i]));
        __m128 my1 = _mm_load_ps(&(y1[i]));
        __m128 mz1 = _mm_load_ps(&(z1[i]));
        __m128 mx2 = _mm_load_ps(&(x2[i]));

        __m128 tx1= _mm_add_ps(_mm_mul_ps(my1,_mm_shuffle_ps(r0,r0,0x55)),_mm_mul_ps(mx1,_mm_shuffle_ps(r0,r0,0x00)));
        __m128 tx2= _mm_add_ps(_mm_mul_ps(mz1,_mm_shuffle_ps(r0,r0,0xAA)),_mm_shuffle_ps(r0,r0,0xFF));
        tx2= _mm_add_ps(tx2,tx1);
        tx2= _mm_sub_ps(tx2,mx2);
        mx2 = _mm_load_ps(&(y2[i]));
        __m128 d1 = _mm_mul_ps(tx2,tx2);

        tx1= _mm_add_ps(_mm_mul_ps(my1,_mm_shuffle_ps(r1,r1,0x55)),_mm_mul_ps(mx1,_mm_shuffle_ps(r1,r1,0x00)));
        tx2= _mm_add_ps(_mm_mul_ps(mz1,_mm_shuffle_ps(r1,r1,0xAA)),_mm_shuffle_ps(r1,r1,0xFF));
        tx2= _mm_add_ps(tx2,tx1);
        tx2= _mm_sub_ps(tx2,mx2);
        mx2= _mm_load_ps(&(z2[i]));
        d1 = _mm_add_ps(d1,_mm_mul_ps(tx2,tx2));

        tx1= _mm_add_ps(_mm_mul_ps(my1,_mm_shuffle_ps(r2,r2,0x55)),_mm_mul_ps(mx1,_mm_shuffle_ps(r2,r2,0x00)));
        tx2= _mm_add_ps(_mm_mul_ps(mz1,_mm_shuffle_ps(r2,r2,0xAA)),_mm_shuffle_ps(r2,r2,0xFF));
        tx2= _mm_add_ps(tx2,tx1);
        tx2= _mm_sub_ps(tx2,mx2);
        d1 = _mm_add_ps(d1,_mm_mul_ps(tx2,tx2));
        _mm_store_ps(&(d[i]),d1); //write out 4 differences
        mx1= _mm_mul_ps(d1,d0);
        mx1= _mm_add_ps(mx1,one);
#ifdef FAST_DIVISION
        mx1= _mm_rcp_ps(mx1);
#else
        mx1=_mm_div_ps(one,mx1);
#endif
        sum=_mm_add_ps(sum,mx1);
    }

    for( ;i <nat; i++){
        float v[4] __attribute__ ((aligned (16))) ={x1[i],y1[i],z1[i],1};
        float w[4] __attribute__ ((aligned (16))) ={x2[i],y2[i],z2[i],1};
        __m128 x1 = _mm_load_ps(v);
        __m128 y1 = _mm_load_ps(w);
        __m128 a1 = _mm_mul_ps(r0,x1);
        __m128 b1 = _mm_mul_ps(r1,x1);

        x1 = _mm_mul_ps(r2,x1);

        a1 = _mm_add_ps(_mm_unpacklo_ps(a1,x1),_mm_unpackhi_ps(a1,x1));
        x1 = _mm_add_ps(_mm_unpacklo_ps(b1,one),_mm_unpackhi_ps(b1,one));
        x1 = _mm_add_ps(_mm_unpacklo_ps(a1,x1),_mm_unpackhi_ps(a1,x1));

        x1 = _mm_sub_ps(x1,y1);
        x1 = _mm_mul_ps(x1,x1);

        x1 = _mm_add_ps(x1, _mm_movehl_ps(x1, x1));
        x1 = _mm_add_ss(x1, _mm_shuffle_ps(x1, x1, 1));

        _mm_store_ss(&(d[i]),x1);
        x1= _mm_mul_ss(x1,d0);
        x1= _mm_add_ss(x1,one);
#ifdef FAST_DIVISION
        x1= _mm_rcp_ps(x1);
#else
        x1=_mm_div_ps(one,x1);
#endif
        sum=_mm_add_ss(sum,x1);
    }
    //not worth writing two versions for SSE/SSE3 - avoid the hadds
    sum=_mm_add_ps(sum,_mm_movehl_ps(sum, sum));
    sum=_mm_add_ss(sum, _mm_shuffle_ps(sum, sum, 1));
    _mm_store_ss(&(fsum),sum);
    return(fsum);
}

int score_fun_soa_sse(int nat, float d0, float d, float *r,float *x1, float *y1, float *z1, float *x2, float *y2, float *z2, int *ialign,int *nalign,float *tm_score){
    //ialign points to atom number
    int k,ncut=0,nchange=0,my_nalign=*nalign,upper_nat= (nat%4)?(nat/4)*4+4 : nat;
    float d2=d*d,invfnat=1.0f/(float)nat;
    float invd0d0=1.0f/(d0*d0);
    float* dist=(float*)mem_align(16,upper_nat*sizeof(float));
    //keep nmin smallest distances && distances < dtmp
    float my_score=LG_score_soa_sse(r,nat,x1,y1,z1,x2,y2,z2,dist,invd0d0);
    for(int k=0;k<nat;k++){
        if(dist[k]<d2) ncut++;
    }

    //adjust d until there are at least 3 the same - rare use another routine for this
    while(ncut <3){
        d2=d*d;
        ncut=0;
        for(k=0;k<nat;k++)
            if(dist[k]<d2) ncut++;
        d+=0.5f;
    }
    ncut=0;
    for(k=0;k<nat;k++){
        if(dist[k]<d2){
            if(ncut < my_nalign && ialign[ncut] == k)ncut++;
            else{
                nchange=1;
                ialign[ncut++]=k;
            }
        }
    }
    free (dist);

    *tm_score=my_score*invfnat;
    if(!nchange)return(0);
    *nalign=ncut;
    return(1);
}

double rmsd_svd(int nat,double *dcoords,double u[3][3], double t[3],bool rmsd_flag){
    const double inv3=1.0/3.0;
    const double dnat=(double)nat;
    const double invdnat=1.0/(double)nat;
    double drms=0,r[9]={0,0,0,0,0,0,0,0,0};
    double s1x=0,s1y=0,s1z=0,s2x=0,s2y=0,s2z=0,ssq=0;
    int m=0;
    for (int i=0;i<nat;i++){
        double c1x=dcoords[m++];
        double c1y=dcoords[m++];
        double c1z=dcoords[m++];
        double c2x=dcoords[m++];
        double c2y=dcoords[m++];
        double c2z=dcoords[m++];
        r[0]+=c1x*c2x; r[1]+=c1x*c2y; r[2]+=c1x*c2z; r[3]+=c1y*c2x; r[4]+=c1y*c2y; r[5]+=c1y*c2z; r[6]+=c1z*c2x; r[7]+=c1z*c2y; r[8]+=c1z*c2z;
        s1x+=c1x;s1y+=c1y;s1z+=c1z;s2x+=c2x;s2y+=c2y;s2z+=c2z;
        if(rmsd_flag) ssq+=c1x*c1x+c1y*c1y+c1z*c1z+c2x*c2x+c2y*c2y+c2z*c2z;
    }
    r[0]-=s1x*s2x*invdnat;
    r[1]-=s1x*s2y*invdnat;
    r[2]-=s1x*s2z*invdnat;
    r[3]-=s1y*s2x*invdnat;
    r[4]-=s1y*s2y*invdnat;
    r[5]-=s1y*s2z*invdnat;
    r[6]-=s1z*s2x*invdnat;
    r[7]-=s1z*s2y*invdnat;
    r[8]-=s1z*s2z*invdnat;
    double det= r[0] * ( (r[4]*r[8]) - (r[5]*r[7]) )- r[1] * ( (r[3]*r[8]) - (r[5]*r[6]) ) + r[2] * ( (r[3]*r[7]) - (r[4]*r[6]) );
    //for symmetric matrix
    //lower triangular matrix rr
    double detsq=det*det;
    //lower triangular matrix rr


    double rr[6]={r[0]*r[0]+ r[1]*r[1]+ r[2]*r[2],
                  r[3]*r[0]+ r[4]*r[1]+ r[5]*r[2],
                  r[3]*r[3]+ r[4]*r[4]+ r[5]*r[5],
                  r[6]*r[0]+ r[7]*r[1]+ r[8]*r[2],
                  r[6]*r[3]+ r[7]*r[4]+ r[8]*r[5],
                  r[6]*r[6]+ r[7]*r[7]+ r[8]*r[8]};

    double spur=((double)(rr[0]+rr[2]+rr[5]))*inv3;
    double cof=((double)(rr[2]*rr[5] - rr[4]*rr[4] + rr[0]*rr[5]- rr[3]*rr[3] + rr[0]*rr[2] - rr[1]*rr[1])) *inv3;
    double e[3] ={spur,spur,spur};
    double h=( spur > 0 )? spur*spur-cof : -1.0;
    if(h>0)
    {
        double g = (spur*cof - detsq)*0.5 - spur*h;
        double sqrth = sqrt(h);
        double d1 = h*h*h - g*g;
        d1= ( d1<0 ) ? atan2(0,-g)*inv3 : atan2(sqrt(d1),-g)*inv3;
        double cth = sqrth * cos(d1);
        double sth = sqrth*SQRT3*sin(d1);
        e[0]+=  cth+cth;
        e[1]+= -cth+sth;
        e[2]+= -cth-sth;
    }
    e[0]=(e[0] < 0) ? 0 : sqrt(e[0]);
    e[1]=(e[1] < 0) ? 0 : sqrt(e[1]);
    e[2]=(e[2] < 0) ? 0 : sqrt(e[2]);

    double d=(det<0)? e[0] + e[1] -e[2] : e[0] + e[1]+e[2];
    if(rmsd_flag){
        double xm=s1x*s1x+s1y*s1y*+s1z*s1z+s2x*s2x+s2y*s2y+s2z*s2z;
        drms=(ssq-xm*dnat-d-d)*invdnat;
        drms=(drms>1e-8)?sqrt(drms) : 0.0f;
    }
    if(rmsd_flag){
        double xm=s1x*s1x+s1y*s1y+s1z*s1z+s2x*s2x+s2y*s2y+s2z*s2z;
        double r2=(ssq-xm*invdnat-d-d)*invdnat;
        drms=(r2>1.0e-8) ? sqrt(r2) : 0.0;
    }
    if(u && t){
        rmatrix(d,r,u);
        t[0] =s2x*invdnat - (u[0][0]*s1x*invdnat + u[1][0]*s1y*invdnat + u[2][0]*s1z*invdnat);
        t[1]= s2y*invdnat - (u[0][1]*s1x*invdnat + u[1][1]*s1y*invdnat + u[2][1]*s1z*invdnat);
        t[2]= s2z*invdnat - (u[0][2]*s1x*invdnat + u[1][2]*s1y*invdnat + u[2][2]*s1z*invdnat);
    }
    return(drms);
}

float tmscore_cpu_soa_sse2(int nat,float *x1, float *y1,float *z1,float *x2,float *y2,float *z2,float bR[3][3], float bt[3],float *rmsd){
    int nalign=0,best_nalign=0;
    int *ialign=new int[nat];
    int *best_align=new int[nat];
    float max_score=-1,rms;
    float r[16] __attribute__ ((aligned (16)));
    //d0
    float d0=1.24*pow((nat-15),(1.0/3.0))-1.8;
    if(d0< 0.5)d0=0.5;
    //d0_search ----->
    float d,d0_search=d0;
    if(d0_search > 8)d0_search=8;
    if(d0_search <4.5)d0_search=4.5;
    //iterative parameters ----->
    int n_it=20;      //maximum number of iterations
    int n_init_max=6; //maximum number of L_init
    int n_init=0;
    int L_ini_min=4;
    int L_ini[6];

    if(nat < 4) L_ini_min=nat;
    int len=nat;
    int divisor=1;
    while(len > L_ini_min && n_init <5){
        L_ini[n_init++]=len;
        divisor*=2;
        len=nat/divisor;
    }
    L_ini[n_init++]=4;
    if (L_ini[n_init-1] > L_ini_min)L_ini[n_init++]=L_ini_min;;

    // find the maximum score starting from local structures superposition
    float score; //TM-score
    for (int seed=0;seed<n_init;seed++)
    {
        //find the initial rotation matrix using the initial seed residues
        int L_init=L_ini[seed];
        for(int istart=0;istart<=nat-L_init;istart++)
        {
            int nchanges=1;
            int nalign=L_init;
            {
                int m=0;
                int n=0;
                for(int i=0;i<nat;i++){
                    if(i>=istart && i<istart+L_init){
                        ialign[n++]=i;
                    }
                }
            }
            if(rmsd && !seed)
                *rmsd=kabsch_quat_soa_sse2(nalign,ialign,x1,y1,z1,x2,y2,z2,r);
            else
                kabsch_quat_soa_sse2(nalign,ialign,x1,y1,z1,x2,y2,z2,r);
            score_fun_soa_sse(nat, d0, d0_search-1,r,x1,y1,z1,x2,y2,z2,ialign,&nalign,&score);
            d=d0_search+1.0f;
            if(score > max_score){
                max_score=score;
                memmove(best_align,ialign,nalign*sizeof(int));
                best_nalign=nalign;
            }
            //extend search from seed
            for (int iter=0;iter<n_it && nchanges;iter++){
                kabsch_quat_soa_sse2(nalign,ialign,x1,y1,z1,x2,y2,z2,r);
                nchanges=score_fun_soa_sse(nat, d0, d,r,x1,y1,z1,x2,y2,z2,ialign,&nalign,&score);
                if(score > max_score){
                    max_score=score;
                    memmove(best_align,ialign,nalign*sizeof(int));
                    best_nalign=nalign;
                }
            }
        }
    }

    //for best frame re-calculate matrix with double precision
    double R[3][3],t[3];
    double *acoords=new double [best_nalign*6];
    for(int k=0;k<best_nalign;k++){
        int i=best_align[k];
        acoords[6*k]  =x1[i];
        acoords[6*k+1]=y1[i];
        acoords[6*k+2]=z1[i];
        acoords[6*k+3]=x2[i];
        acoords[6*k+4]=y2[i];
        acoords[6*k+5]=z2[i];
    }
    rmsd_svd(best_nalign,acoords,R, t,0);
    if(bR){
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                bR[i][j]=R[i][j];
    }
    if(bt){
        for(int i=0;i<3;i++)
            bt[i]=t[i];
    }
    double invd0d0=1.0/(double)(d0*d0);
    double dist;
    float invdnat=1.0/(double)nat;
    double my_score=0;
    for(int k=0;k<nat;k++){
        double u[3];
        double v[3]={x1[k],y1[k],z1[k]};
        double w[3]={x2[k],y2[k],z2[k]};
        u[0]=t[0]+R[0][0]*v[0]+R[1][0]*v[1]+R[2][0]*v[2]-w[0];
        u[1]=t[1]+R[0][1]*v[0]+R[1][1]*v[1]+R[2][1]*v[2]-w[1];
        u[2]=t[2]+R[0][2]*v[0]+R[1][2]*v[1]+R[2][2]*v[2]-w[2];
        dist=u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
        my_score+=1.0/(1.0+dist*invd0d0);
    }

    delete [] best_align;
    delete [] ialign;
    delete [] acoords;
    return(my_score*invdnat);
}


//     1, collect those residues with dis<d;
//     2, calculate TMscore
//int score_fun8( Coordinates &xa, Coordinates &ya, int n_ali, float d, int i_ali[],
//                float *score1, int score_sum_method, const float Lnorm,
//                const float score_d8, const float d0
//                , float * mem
//
//                )
//{
//    float score_sum=0, di;
//    float d_tmp=d*d;
//    float d02=d0*d0;
//    float score_d8_cut = score_d8*score_d8;
//
//    int i, n_cut, inc=0;
//    float *distArray = mem;
//    while(1)
//    {
//        n_cut=0;
//        score_sum=0;
//        simd_float sum = simdf32_set(0);
//        simd_float vscore_d8_cut = simdf32_set(score_d8_cut);
//        simd_float vd02 = simdf32_set(d02);
//        simd_float one = simdf32_set(1.0f);
//        for(i=0; i < n_ali; i+=VECSIZE_FLOAT){
//            //    float d1=xx-yx;
//            //    float d2=xy-yy;
//            //    float d3=xz-yz;
//            //    return (d1*d1 + d2*d2 + d3*d3);
//            simd_float xa_x = simdf32_load(&xa.x[i]);
//            simd_float ya_x = simdf32_load(&ya.x[i]);
//            simd_float xa_y = simdf32_load(&xa.y[i]);
//            simd_float ya_y = simdf32_load(&ya.y[i]);
//            simd_float xa_z = simdf32_load(&xa.z[i]);
//            simd_float ya_z = simdf32_load(&ya.z[i]);
//            ya_x = simdf32_sub(xa_x, ya_x);
//            ya_y = simdf32_sub(xa_y, ya_y);
//            ya_z = simdf32_sub(xa_z, ya_z);
//            ya_x = simdf32_mul(ya_x, ya_x);
//            ya_y = simdf32_mul(ya_y, ya_y);
//            ya_z = simdf32_mul(ya_z, ya_z);
//            simd_float res = simdf32_add(ya_x, ya_y);
//            simd_float di = simdf32_add(res, ya_z);
//            simdf32_store(&distArray[i], di);
//            simd_float di_lt_score_d8 = simdf32_lt(di, vscore_d8_cut);
//            simd_float oneDividedDist = simdf32_div(one, simdf32_add(one, simdf32_div(di,vd02)));
//            sum = simdf32_add(sum, (simd_float)simdi_and((simd_int) di_lt_score_d8, (simd_int) oneDividedDist ));
//        }
//        for(i=0; i < VECSIZE_FLOAT; i++){
//            score_sum+=((float*)&sum)[i];
//        }
//
//        for(i=0; i<n_ali; i++)
//        {
//            di = distArray[i];
//            i_ali[n_cut]=i;
//            n_cut+=(di<d_tmp);
//            //score_sum += (di<=score_d8_cut) ? 1/(1+di/d02) : 0;
//            //else score_sum += 1/(1+di/d02);
//        }
//        //there are not enough feasible pairs, reliefe the threshold
//        if(n_cut<3 && n_ali>3)
//        {
//            inc++;
//            double dinc=(d+inc*0.5);
//            d_tmp = dinc * dinc;
//        }
//        else break;
//    }
//
//    *score1=score_sum/Lnorm;
//    return n_cut;
//}


//     1, collect those residues with dis<d;
//     2, calculate TMscore
int score_fun8( Coordinates &xa, Coordinates &ya, int n_ali, float d, int i_ali[],
                float *score1, int score_sum_method, const float Lnorm,
                const float score_d8, const float d0, float * mem)
{
    float score_sum=0, di;
    float d_tmp=d*d;
    float d02=d0*d0;
    float score_d8_cut = score_d8*score_d8;

    int i, n_cut, inc=0;
    float *distArray = mem;
    int anat= (n_ali%4) ? (n_ali/VECSIZE_FLOAT)*VECSIZE_FLOAT+VECSIZE_FLOAT : n_ali;


    while(1)
    {
        n_cut=0;
        score_sum=0;
        simd_float sum = simdf32_set(0);
        simd_float vscore_d8_cut = simdf32_set(score_d8_cut);
        simd_float vd02 = simdf32_set(d02);
        simd_float one = simdf32_set(1.0f);
        for(i=0; i < n_ali; i+=VECSIZE_FLOAT){
            //    float d1=xx-yx;
            //    float d2=xy-yy;
            //    float d3=xz-yz;
            //    return (d1*d1 + d2*d2 + d3*d3);
            simd_float xa_x = simdf32_load(&xa.x[i]);
            simd_float ya_x = simdf32_load(&ya.x[i]);
            simd_float xa_y = simdf32_load(&xa.y[i]);
            simd_float ya_y = simdf32_load(&ya.y[i]);
            simd_float xa_z = simdf32_load(&xa.z[i]);
            simd_float ya_z = simdf32_load(&ya.z[i]);
            ya_x = simdf32_sub(xa_x, ya_x);
            ya_y = simdf32_sub(xa_y, ya_y);
            ya_z = simdf32_sub(xa_z, ya_z);
            ya_x = simdf32_mul(ya_x, ya_x);
            ya_y = simdf32_mul(ya_y, ya_y);
            ya_z = simdf32_mul(ya_z, ya_z);
            simd_float res = simdf32_add(ya_x, ya_y);
            simd_float di = simdf32_add(res, ya_z);
            simdf32_store(&distArray[i], di);
            simd_float di_lt_score_d8 = simdf32_lt(di, vscore_d8_cut);
            simd_float oneDividedDist = simdf32_div(one, simdf32_add(one, simdf32_div(di,vd02)));
            sum = simdf32_add(sum, (simd_float)simdi_and((simd_int) di_lt_score_d8, (simd_int) oneDividedDist ));
        }
        for(i=0; i < VECSIZE_FLOAT; i++){
            score_sum+=((float*)&sum)[i];
        }

        for(i=0; i<n_ali; i++)
        {
            di = distArray[i];
            i_ali[n_cut]=i;
            n_cut+=(di<d_tmp);
            //score_sum += (di<=score_d8_cut) ? 1/(1+di/d02) : 0;
            //else score_sum += 1/(1+di/d02);
        }
        //there are not enough feasible pairs, reliefe the threshold
        if(n_cut<3 && n_ali>3)
        {
            inc++;
            double dinc=(d+inc*0.5);
            d_tmp = dinc * dinc;
        }
        else break;
    }

    *score1=score_sum/Lnorm;
    return n_cut;
}


//int score_fun8( Coordinates &xa, Coordinates &ya, int n_ali, float d, int i_ali[],
//                float *score1, int score_sum_method, const float Lnorm,
//                const float score_d8, const float d0, float * mem)
//{
//    float score_sum=0, di;
//    float d_tmp=d*d;
//    float d02=d0*d0;
//    float score_d8_cut = score_d8*score_d8;
//
//    int i, n_cut, inc=0;
//
//    while(1)
//    {
//        n_cut=0;
//        score_sum=0;
//        for(i=0; i<n_ali; i++)
//        {
//            di = dist(xa.x[i], xa.y[i], xa.z[i], ya.x[i], ya.y[i], ya.z[i]);
//            if(di<d_tmp)
//            {
//                i_ali[n_cut]=i;
//                n_cut++;
//            }
//            if(score_sum_method==8)
//            {
//                if(di<=score_d8_cut) score_sum += 1/(1+di/d02);
//            }
//            else score_sum += 1/(1+di/d02);
//        }
//        //there are not enough feasible pairs, reliefe the threshold
//        if(n_cut<3 && n_ali>3)
//        {
//            inc++;
//            double dinc=(d+inc*0.5);
//            d_tmp = dinc * dinc;
//        }
//        else break;
//    }
//
//    *score1=score_sum/Lnorm;
//    return n_cut;
//}

int score_fun8_standard(Coordinates &xa, Coordinates &ya, int n_ali, float d,
                        int i_ali[], float *score1, int score_sum_method,
                        float score_d8, float d0)
{
    float score_sum = 0, di;
    float d_tmp = d*d;
    float d02 = d0*d0;
    float score_d8_cut = score_d8*score_d8;

    int i, n_cut, inc = 0;
    while (1)
    {
        n_cut = 0;
        score_sum = 0;
        for (i = 0; i<n_ali; i++)
        {
            di = dist(xa.x[i], xa.y[i], xa.z[i],
                      ya.x[i], ya.y[i], ya.z[i]);
            if (di<d_tmp)
            {
                i_ali[n_cut] = i;
                n_cut++;
            }
            if (score_sum_method == 8)
            {
                if (di <= score_d8_cut) score_sum += 1 / (1 + di / d02);
            }
            else
            {
                score_sum += 1 / (1 + di / d02);
            }
        }
        //there are not enough feasible pairs, reliefe the threshold
        if (n_cut<3 && n_ali>3)
        {
            inc++;
            double dinc = (d + inc*0.5);
            d_tmp = dinc * dinc;
        }
        else break;
    }

    *score1 = score_sum / n_ali;
    return n_cut;
}

bool KabschFast(Coordinates & x,
                 Coordinates & y,
                 int n,
                 int mode,
                 float *rms,
                 float t[3],
                 float u[3][3],
                 float * mem){
    float r[16] __attribute__ ((aligned (16)));
    float rmsdbla = kabsch_quat_soa_avx(n, NULL, x.x, x.y, x.z, y.x, y.y, y.z, r, mem);
    t[0] = r[3];
    t[1] = r[7];
    t[2] = r[11];
    u[0][0] = r[0];
    u[0][1] = r[1];
    u[0][2] = r[2];
    u[1][0] = r[4];
    u[1][1] = r[5];
    u[1][2] = r[6];
    u[2][0] = r[8];
    u[2][1] = r[9];
    u[2][2] = r[10];
    return true;
}

double TMscore8_search(Coordinates &r1, Coordinates &r2,
                       Coordinates &xtm, Coordinates & ytm,
                       Coordinates &xt, int Lali, float t0[3], float u0[3][3], int simplify_step,
                       int score_sum_method, float *Rcomm, float local_d0_search, float Lnorm,
                       float score_d8, float d0, float * mem)
{
    int i, m;
    float score_max, score, rmsd;
    const int kmax=Lali;
    int k_ali[kmax], ka, k;
    float t[3];
    float u[3][3];
    float d;


    //iterative parameters
    int n_it=20;            //maximum number of iterations
    int n_init_max=6; //maximum number of different fragment length
    int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4
    int L_ini_min=4;
    if(Lali<L_ini_min) L_ini_min=Lali;

    int n_init=0, i_init;
    for(i=0; i<n_init_max-1; i++)
    {
        n_init++;
        L_ini[i]=(int) (Lali/pow(2.0, (double) i));
        if(L_ini[i]<=L_ini_min)
        {
            L_ini[i]=L_ini_min;
            break;
        }
    }
    if(i==n_init_max-1)
    {
        n_init++;
        L_ini[i]=L_ini_min;
    }

    score_max=-1;
    //find the maximum score starting from local structures superposition
    int i_ali[kmax], n_cut;
    int L_frag; //fragment length
    int iL_max; //maximum starting postion for the fragment

    for(i_init=0; i_init<n_init; i_init++)
    {
        L_frag=L_ini[i_init];
        iL_max=Lali-L_frag;

        i=0;
        while(1)
        {
            //extract the fragment starting from position i
            ka=0;
            for(k=0; k<L_frag; k++)
            {
                int kk=k+i;
                r1.x[k]=xtm.x[kk];
                r1.y[k]=xtm.y[kk];
                r1.z[k]=xtm.z[kk];

                r2.x[k]=ytm.x[kk];
                r2.y[k]=ytm.y[kk];
                r2.z[k]=ytm.z[kk];

                k_ali[ka]=kk;
                ka++;
            }

            //extract rotation matrix based on the fragment
            //float r[16] __attribute__ ((aligned (16)));
            //float rmsdbla = kabsch_quat_soa_sse2(L_frag, NULL, r1.x, r1.y, r1.z, r2.x, r2.y, r2.z, r);
            //float rmsdbla;
//            float bR[3][3];
//            float bt[3];
            //float TMscore = tmscore_cpu_soa_sse2(L_frag, r1.x, r1.y, r1.z, r2.x, r2.y, r2.z, bR, bt, &rmsdbla);
            KabschFast(r1, r2, L_frag, 1, &rmsd, t, u, mem);
            //Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);//, mem);

            if (simplify_step != 1)
                *Rcomm = 0;
            do_rotation(xtm, xt, Lali, t, u);

            //get subsegment of this fragment
            d = local_d0_search - 1;
            n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score,
                             score_sum_method, Lnorm, score_d8, d0, mem);
            if(score>score_max)
            {
                score_max=score;

                //save the rotation matrix
                for(k=0; k<3; k++)
                {
                    t0[k]=t[k];
                    u0[k][0]=u[k][0];
                    u0[k][1]=u[k][1];
                    u0[k][2]=u[k][2];
                }
            }

            //try to extend the alignment iteratively
            d = local_d0_search + 1;
            for(int it=0; it<n_it; it++)
            {
                ka=0;
                for(k=0; k<n_cut; k++)
                {
                    m=i_ali[k];
                    r1.x[k]=xtm.x[m];
                    r1.y[k]=xtm.y[m];
                    r1.z[k]=xtm.z[m];

                    r2.x[k]=ytm.x[m];
                    r2.y[k]=ytm.y[m];
                    r2.z[k]=ytm.z[m];

                    k_ali[ka]=m;
                    ka++;
                }
                //KabschFast();
                //extract rotation matrix based on the fragment
                //float rmsdbla = kabsch_quat_soa_sse2(L_frag, NULL, r1.x, r1.y, r1.z, r2.x, r2.y, r2.z, r);

                KabschFast(r1, r2, n_cut, 1, &rmsd, t, u, mem);
//                Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);

                do_rotation(xtm, xt, Lali, t, u);
                n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score,
                                 score_sum_method, Lnorm, score_d8, d0, mem);
                if(score>score_max)
                {
                    score_max=score;

                    //save the rotation matrix
                    for(k=0; k<3; k++)
                    {
                        t0[k]=t[k];
                        u0[k][0]=u[k][0];
                        u0[k][1]=u[k][1];
                        u0[k][2]=u[k][2];
                    }
                }

                //check if it converges
                if(n_cut==ka)
                {
                    for(k=0; k<n_cut; k++)
                    {
                        if(i_ali[k]!=k_ali[k]) break;
                    }
                    if(k==n_cut) break;
                }
            } //for iteration

            if(i<iL_max)
            {
                i=i+simplify_step; //shift the fragment
                if(i>iL_max) i=iL_max;  //do this to use the last missed fragment
            }
            else if(i>=iL_max) break;
        }//while(1)
        //end of one fragment
    }//for(i_init
    return score_max;
}


double TMscore8_search_standard(Coordinates &r1, Coordinates &r2,
                                Coordinates &xtm, Coordinates &ytm, Coordinates &xt, int Lali,
                                float t0[3], float u0[3][3], int simplify_step, int score_sum_method,
                                float *Rcomm, float local_d0_search, float score_d8, float d0, float * mem)
{
    int i, m;
    float score_max, score, rmsd;
    const int kmax = Lali;
    int k_ali[kmax], ka, k;
    float t[3];
    float u[3][3];
    float d;

    //iterative parameters
    int n_it = 20;            //maximum number of iterations
    int n_init_max = 6; //maximum number of different fragment length
    int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4
    int L_ini_min = 4;
    if (Lali<L_ini_min) L_ini_min = Lali;

    int n_init = 0, i_init;
    for (i = 0; i<n_init_max - 1; i++)
    {
        n_init++;
        L_ini[i] = (int)(Lali / pow(2.0, (double)i));
        if (L_ini[i] <= L_ini_min)
        {
            L_ini[i] = L_ini_min;
            break;
        }
    }
    if (i == n_init_max - 1)
    {
        n_init++;
        L_ini[i] = L_ini_min;
    }

    score_max = -1;
    //find the maximum score starting from local structures superposition
    int i_ali[kmax], n_cut;
    int L_frag; //fragment length
    int iL_max; //maximum starting postion for the fragment

    for (i_init = 0; i_init<n_init; i_init++)
    {
        L_frag = L_ini[i_init];
        iL_max = Lali - L_frag;

        i = 0;
        while (1)
        {
            //extract the fragment starting from position i
            ka = 0;
            for (k = 0; k<L_frag; k++)
            {
                int kk = k + i;
                r1.x[k] = xtm.x[kk];
                r1.y[k] = xtm.y[kk];
                r1.z[k] = xtm.z[kk];

                r2.x[k] = ytm.x[kk];
                r2.y[k] = ytm.y[kk];
                r2.z[k] = ytm.z[kk];

                k_ali[ka] = kk;
                ka++;
            }
            //extract rotation matrix based on the fragment
            KabschFast(r1, r2, L_frag, 1, &rmsd, t, u, mem);
            if (simplify_step != 1)
                *Rcomm = 0;
            do_rotation(xtm, xt, Lali, t, u);

            //get subsegment of this fragment
            d = local_d0_search - 1;
            n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, &score,
                                        score_sum_method, score_d8, d0);

            if (score>score_max)
            {
                score_max = score;

                //save the rotation matrix
                for (k = 0; k<3; k++)
                {
                    t0[k] = t[k];
                    u0[k][0] = u[k][0];
                    u0[k][1] = u[k][1];
                    u0[k][2] = u[k][2];
                }
            }

            //try to extend the alignment iteratively
            d = local_d0_search + 1;
            for (int it = 0; it<n_it; it++)
            {
                ka = 0;
                for (k = 0; k<n_cut; k++)
                {
                    m = i_ali[k];
                    r1.x[k] = xtm.x[m];
                    r1.y[k] = xtm.y[m];
                    r1.z[k] = xtm.z[m];

                    r2.x[k] = ytm.x[m];
                    r2.y[k] = ytm.y[m];
                    r2.z[k] = ytm.z[m];

                    k_ali[ka] = m;
                    ka++;
                }
                //extract rotation matrix based on the fragment
                KabschFast(r1, r2, n_cut, 1, &rmsd, t, u, mem);
                do_rotation(xtm, xt, Lali, t, u);
                n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, &score,
                                            score_sum_method, score_d8, d0);
                if (score>score_max)
                {
                    score_max = score;

                    //save the rotation matrix
                    for (k = 0; k<3; k++)
                    {
                        t0[k] = t[k];
                        u0[k][0] = u[k][0];
                        u0[k][1] = u[k][1];
                        u0[k][2] = u[k][2];
                    }
                }

                //check if it converges
                if (n_cut == ka)
                {
                    for (k = 0; k<n_cut; k++)
                    {
                        if (i_ali[k] != k_ali[k]) break;
                    }
                    if (k == n_cut) break;
                }
            } //for iteration

            if (i<iL_max)
            {
                i = i + simplify_step; //shift the fragment
                if (i>iL_max) i = iL_max;  //do this to use the last missed fragment
            }
            else if (i >= iL_max) break;
        }//while(1)
        //end of one fragment
    }//for(i_init
    return score_max;
}

//Comprehensive TMscore search engine
// input:   two vector sets: x, y
//          an alignment invmap0[] between x and y
//          simplify_step: 1 or 40 or other integers
//          score_sum_method: 0 for score over all pairs
//                            8 for socre over the pairs with dist<score_d8
// output:  the best rotaion matrix t, u that results in highest TMscore
double detailed_search(Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                       Coordinates &xt, const Coordinates &x, const Coordinates &y, int xlen, int ylen,
                       int invmap0[], float t[3], float u[3][3], int simplify_step,
                       int score_sum_method, float local_d0_search, float Lnorm,
                       float score_d8, float d0, float * mem)
{
    //x is model, y is template, try to superpose onto y
    int i, j, k;
    float tmscore;
    float rmsd;

    k=0;
    for(i=0; i<ylen; i++)
    {
        j=invmap0[i];
        if(j>=0) //aligned
        {
            xtm.x[k]=x.x[j];
            xtm.y[k]=x.y[j];
            xtm.z[k]=x.z[j];

            ytm.x[k]=y.x[i];
            ytm.y[k]=y.y[i];
            ytm.z[k]=y.z[i];
            k++;
        }
    }

    //detailed search 40-->1
    tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u, simplify_step,
                              score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0, mem);
    return tmscore;
}

double detailed_search_standard( Coordinates &r1, Coordinates &r2,
                                 Coordinates &xtm, Coordinates &ytm, Coordinates &xt,
                                 const Coordinates &x, const Coordinates &y,
                                 int xlen, int ylen, int invmap0[], float t[3], float u[3][3],
                                 int simplify_step, int score_sum_method, double local_d0_search,
                                 const bool& bNormalize, float Lnorm, float score_d8, float d0, float * mem)
{
    //x is model, y is template, try to superpose onto y
    int i, j, k;
    float tmscore;
    float rmsd;

    k=0;
    for(i=0; i<ylen; i++)
    {
        j=invmap0[i];
        if(j>=0) //aligned
        {
            xtm.x[k]=x.x[j];
            xtm.y[k]=x.y[j];
            xtm.z[k]=x.z[j];

            ytm.x[k]=y.x[i];
            ytm.y[k]=y.y[i];
            ytm.z[k]=y.z[i];
            k++;
        }
    }

    //detailed search 40-->1
    tmscore = TMscore8_search_standard( r1, r2, xtm, ytm, xt, k, t, u,
                                        simplify_step, score_sum_method, &rmsd, local_d0_search, score_d8, d0, mem);
    if (bNormalize)// "-i", to use standard_TMscore, then bNormalize=true, else bNormalize=false;
        tmscore = tmscore * k / Lnorm;

    return tmscore;
}

//compute the score quickly in three iterations
double get_score_fast( Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                       const Coordinates &x, const Coordinates &y, int xlen, int ylen, int invmap[],
                       float d0, float d0_search, float t[3], float u[3][3], float * mem)
{
    float rms, tmscore, tmscore1, tmscore2;
    int i, j, k;

    k=0;
    for(j=0; j<ylen; j++)
    {
        i=invmap[j];
        if(i>=0)
        {
            r1.x[k]=x.x[i];
            r1.y[k]=x.y[i];
            r1.z[k]=x.z[i];

            r2.x[k]=y.x[j];
            r2.y[k]=y.y[j];
            r2.z[k]=y.z[j];

            xtm.x[k]=x.x[i];
            xtm.y[k]=x.y[i];
            xtm.z[k]=x.z[i];

            ytm.x[k]=y.x[j];
            ytm.y[k]=y.y[j];
            ytm.z[k]=y.z[j];

            k++;
        }
        else if(i!=-1) PrintErrorAndQuit("Wrong map!\n");
    }
    KabschFast(r1, r2, k, 1, &rms, t, u, mem);

    //evaluate score
    double di;
    const int len=k;
    double dis[len];
    double d00=d0_search;
    double d002=d00*d00;
    double d02=d0*d0;

    int n_ali=k;
    float xrot[3];
    tmscore=0;
    for(k=0; k<n_ali; k++)
    {
        transform(t, u, xtm.x[k], xtm.y[k], xtm.z[k], xrot[0], xrot[1], xrot[2]);
        di=dist(xrot[0], xrot[1], xrot[2],
                ytm.x[k], ytm.y[k], ytm.z[k]);
        dis[k]=di;
        tmscore += 1/(1+di/d02);
    }



    //second iteration
    float d002t=d002;
    while(1)
    {
        j=0;
        for(k=0; k<n_ali; k++)
        {
            if(dis[k]<=d002t)
            {
                r1.x[j]=xtm.x[k];
                r1.y[j]=xtm.y[k];
                r1.z[j]=xtm.z[k];

                r2.x[j]=ytm.x[k];
                r2.y[j]=ytm.y[k];
                r2.z[j]=ytm.z[k];

                j++;
            }
        }
        //there are not enough feasible pairs, relieve the threshold
        if(j<3 && n_ali>3) d002t += 0.5;
        else break;
    }

    if(n_ali!=j)
    {
        KabschFast(r1, r2, j, 1, &rms, t, u, mem);
        tmscore1=0;
        for(k=0; k<n_ali; k++)
        {
            transform(t, u, xtm.x[k],  xtm.y[k], xtm.z[k], xrot[0], xrot[1], xrot[2]);
            di=dist(xrot[0], xrot[1], xrot[2],
                    ytm.x[k], ytm.y[k], ytm.z[k]);
            dis[k]=di;
            tmscore1 += 1/(1+di/d02);
        }

        //third iteration
        d002t=d002+1;

        while(1)
        {
            j=0;
            for(k=0; k<n_ali; k++)
            {
                if(dis[k]<=d002t)
                {
                    r1.x[j]=xtm.x[k];
                    r1.y[j]=xtm.y[k];
                    r1.z[j]=xtm.z[k];
                    r2.x[j]=ytm.x[k];
                    r2.y[j]=ytm.y[k];
                    r2.z[j]=ytm.z[k];
                    j++;
                }
            }
            //there are not enough feasible pairs, relieve the threshold
            if(j<3 && n_ali>3) d002t += 0.5;
            else break;
        }

        //evaluate the score
        KabschFast(r1, r2, j, 1, &rms, t, u, mem);
        tmscore2=0;
        for(k=0; k<n_ali; k++)
        {
            transform(t, u, xtm.x[k],  xtm.y[k], xtm.z[k], xrot[0], xrot[1], xrot[2]);
            di=dist(xrot[0], xrot[1], xrot[2],
                    ytm.x[k], ytm.y[k], ytm.z[k]);
            tmscore2 += 1/(1+di/d02);
        }
    }
    else
    {
        tmscore1=tmscore;
        tmscore2=tmscore;
    }

    if(tmscore1>=tmscore) tmscore=tmscore1;
    if(tmscore2>=tmscore) tmscore=tmscore2;
    return tmscore; // no need to normalize this score because it will not be used for latter scoring
}


//perform gapless threading to find the best initial alignment
//input: x, y, xlen, ylen
//output: y2x0 stores the best alignment: e.g.,
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
double get_initial(Coordinates &r1, Coordinates &r2,
                   Coordinates &xtm, Coordinates &ytm,
                   const Coordinates &x, const Coordinates &y, int xlen, int ylen, int *y2x,
                   float d0, float d0_search, const bool fast_opt,
                   float t[3], float u[3][3], float * mem)
{
    int min_len=getmin(xlen, ylen);
    if(min_len<=5) PrintErrorAndQuit("Sequence is too short <=5!\n");

    int min_ali= min_len/2;              //minimum size of considered fragment
    if(min_ali<=5)  min_ali=5;
    int n1, n2;
    n1 = -ylen+min_ali;
    n2 = xlen-min_ali;

    int i, j, k, k_best;
    float tmscore, tmscore_max=-1;

    k_best=n1;
    for(k=n1; k<=n2; k+=(fast_opt)?5:1)
    {
        //get the map
        for(j=0; j<ylen; j++)
        {
            i=j+k;
            if(i>=0 && i<xlen) y2x[j]=i;
            else y2x[j]=-1;
        }

        //evaluate the map quickly in three iterations
        //this is not real tmscore, it is used to evaluate the goodness of the initial alignment
        tmscore=get_score_fast(r1, r2, xtm, ytm,
                               x, y, xlen, ylen, y2x, d0,d0_search, t, u, mem);
        if(tmscore>=tmscore_max)
        {
            tmscore_max=tmscore;
            k_best=k;
        }
    }

    //extract the best map
    k=k_best;
    for(j=0; j<ylen; j++)
    {
        i=j+k;
        if(i>=0 && i<xlen) y2x[j]=i;
        else y2x[j]=-1;
    }

    return tmscore_max;
}

void smooth(int *sec, int len)
{
    int i, j;
    //smooth single  --x-- => -----
    for (i=2; i<len-2; i++)
    {
        if(sec[i]==2 || sec[i]==4)
        {
            j=sec[i];
            if (sec[i-2]!=j && sec[i-1]!=j && sec[i+1]!=j && sec[i+2]!=j)
                sec[i]=1;
        }
    }

    //   smooth double
    //   --xx-- => ------
    for (i=0; i<len-5; i++)
    {
        //helix
        if (sec[i]!=2   && sec[i+1]!=2 && sec[i+2]==2 && sec[i+3]==2 &&
            sec[i+4]!=2 && sec[i+5]!= 2)
        {
            sec[i+2]=1;
            sec[i+3]=1;
        }

        //beta
        if (sec[i]!=4   && sec[i+1]!=4 && sec[i+2]==4 && sec[i+3]==4 &&
            sec[i+4]!=4 && sec[i+5]!= 4)
        {
            sec[i+2]=1;
            sec[i+3]=1;
        }
    }

    //smooth connect
    for (i=0; i<len-2; i++)
    {
        if (sec[i]==2 && sec[i+1]!=2 && sec[i+2]==2) sec[i+1]=2;
        else if(sec[i]==4 && sec[i+1]!=4 && sec[i+2]==4) sec[i+1]=4;
    }

}

int sec_str(float dis13, float dis14, float dis15,
            float dis24, float dis25, float dis35)
{
    int s=1;

    float delta=2.1;
    if (fabsf(dis15-6.37)<delta && fabsf(dis14-5.18)<delta &&
        fabsf(dis25-5.18)<delta && fabsf(dis13-5.45)<delta &&
        fabsf(dis24-5.45)<delta && fabsf(dis35-5.45)<delta)
    {
        s=2; //helix
        return s;
    }

    delta=1.42;
    if (fabsf(dis15-13  )<delta && fabsf(dis14-10.4)<delta &&
        fabsf(dis25-10.4)<delta && fabsf(dis13-6.1 )<delta &&
        fabsf(dis24-6.1 )<delta && fabsf(dis35-6.1 )<delta)
    {
        s=4; //strand
        return s;
    }

    if (dis15 < 8) s=3; //turn
    return s;
}


//1->coil, 2->helix, 3->turn, 4->strand
void make_sec(Coordinates &x, int len, int *sec)
{
    int j1, j2, j3, j4, j5;
    float d13, d14, d15, d24, d25, d35;
    for(int i=0; i<len; i++)
    {
        sec[i]=1;
        j1=i-2;
        j2=i-1;
        j3=i;
        j4=i+1;
        j5=i+2;

        if(j1>=0 && j5<len)
        {
            d13=sqrt(dist(x.x[j1], x.y[j1], x.z[j1], x.x[j3], x.y[j3], x.z[j3]));
            d14=sqrt(dist(x.x[j1], x.y[j1], x.z[j1], x.x[j4], x.y[j4], x.z[j4]));
            d15=sqrt(dist(x.x[j1], x.y[j1], x.z[j1], x.x[j5], x.y[j5], x.z[j5]));
            d24=sqrt(dist(x.x[j2], x.y[j2], x.z[j2], x.x[j4], x.y[j4], x.z[j4]));
            d25=sqrt(dist(x.x[j2], x.y[j2], x.z[j2], x.x[j5], x.y[j5], x.z[j5]));
            d35=sqrt(dist(x.x[j3], x.y[j3], x.z[j3], x.x[j5], x.y[j5], x.z[j5]));
            sec[i]=sec_str(d13, d14, d15, d24, d25, d35);
        }
    }
}


//get initial alignment from secondary structure alignment
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g.,
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
void get_initial_ss( float **score, bool **path, float **val,
                     const int *secx, const int *secy, int xlen, int ylen, int *y2x)
{
    double gap_open=-1.0;
    NWDP_TM(score, path, val, secx, secy, xlen, ylen, gap_open, y2x);
}


// get_initial5 in TMalign fortran, get_intial_local in TMalign c by yangji
//get initial alignment of local structure superposition
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g.,
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
bool get_initial5( Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                   float **score, bool **path, float **val,
                   const Coordinates &x, const Coordinates &y, int xlen, int ylen, int *y2x,
                   float d0, float d0_search, const bool fast_opt, const float D0_MIN, float * mem)
{
    float GL, rmsd;
    float t[3];
    float u[3][3];

    float d01 = d0 + 1.5;
    if (d01 < D0_MIN) d01 = D0_MIN;
    float d02 = d01*d01;

    float GLmax = 0;
    int aL = getmin(xlen, ylen);
    int *invmap = new int[ylen + 1];

    // jump on sequence1-------------->
    int n_jump1 = 0;
    if (xlen > 250)
        n_jump1 = 45;
    else if (xlen > 200)
        n_jump1 = 35;
    else if (xlen > 150)
        n_jump1 = 25;
    else
        n_jump1 = 15;
    if (n_jump1 > (xlen / 3))
        n_jump1 = xlen / 3;

    // jump on sequence2-------------->
    int n_jump2 = 0;
    if (ylen > 250)
        n_jump2 = 45;
    else if (ylen > 200)
        n_jump2 = 35;
    else if (ylen > 150)
        n_jump2 = 25;
    else
        n_jump2 = 15;
    if (n_jump2 > (ylen / 3))
        n_jump2 = ylen / 3;

    // fragment to superimpose-------------->
    int n_frag[2] = { 20, 100 };
    if (n_frag[0] > (aL / 3))
        n_frag[0] = aL / 3;
    if (n_frag[1] > (aL / 2))
        n_frag[1] = aL / 2;

    // start superimpose search-------------->
    if (fast_opt)
    {
        n_jump1*=5;
        n_jump2*=5;
    }
    bool flag = false;
    for (int i_frag = 0; i_frag < 2; i_frag++)
    {
        int m1 = xlen - n_frag[i_frag] + 1;
        int m2 = ylen - n_frag[i_frag] + 1;

        for (int i = 0; i<m1; i = i + n_jump1) //index starts from 0, different from FORTRAN
        {
            for (int j = 0; j<m2; j = j + n_jump2)
            {
                for (int k = 0; k<n_frag[i_frag]; k++) //fragment in y
                {
                    r1.x[k] = x.x[k + i];
                    r1.y[k] = x.y[k + i];
                    r1.z[k] = x.z[k + i];

                    r2.x[k] = y.x[k + j];
                    r2.y[k] = y.y[k + j];
                    r2.z[k] = y.z[k + j];
                }

                // superpose the two structures and rotate it
                KabschFast(r1, r2, n_frag[i_frag], 1, &rmsd, t, u, mem);

                float gap_open = 0.0;
                NWDP_TM(score, path, val,
                        x, y, xlen, ylen, t, u, d02, gap_open, invmap, mem);
                GL = get_score_fast(r1, r2, xtm, ytm, x, y, xlen, ylen,
                                    invmap, d0, d0_search, t, u, mem);
                if (GL>GLmax)
                {
                    GLmax = GL;
                    for (int ii = 0; ii<ylen; ii++) y2x[ii] = invmap[ii];
                    flag = true;
                }
            }
        }
    }

    delete[] invmap;
    return flag;
}

void score_matrix_rmsd_sec( Coordinates &r1,  Coordinates &r2,
                            float **score, const int *secx, const int *secy,
                            const Coordinates &x, const Coordinates &y, int xlen, int ylen,
                            int *y2x, const float D0_MIN, float d0, float * mem)
{
    float t[3], u[3][3];
    float rmsd, dij;
    float d01=d0+1.5;
    if(d01 < D0_MIN) d01=D0_MIN;
    float d02=d01*d01;

    float xx[3];
    int i, k=0;
    for(int j=0; j<ylen; j++)
    {
        i=y2x[j];
        if(i>=0)
        {
            r1.x[k]=x.x[i];
            r1.y[k]=x.y[i];
            r1.z[k]=x.z[i];

            r2.x[k]=y.x[j];
            r2.y[k]=y.y[j];
            r2.z[k]=y.z[j];

            k++;
        }
    }
    KabschFast(r1, r2, k, 1, &rmsd, t, u, mem);


    for(int ii=0; ii<xlen; ii++)
    {
        transform(t, u, x.x[ii], x.y[ii], x.z[ii], xx[0], xx[1], xx[2]);
        for(int jj=0; jj<ylen; jj++)
        {
            dij=dist(xx[0], xx[1], xx[2], y.x[jj], y.y[jj], y.z[jj]);
            if (secx[ii]==secy[jj])
                score[ii+1][jj+1] = 1.0/(1+dij/d02) + 0.5;
            else
                score[ii+1][jj+1] = 1.0/(1+dij/d02);
        }
    }
}


//get initial alignment from secondary structure and previous alignments
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g.,
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
void get_initial_ssplus(Coordinates &r1, Coordinates &r2, float **score, bool **path,
                        float **val, const int *secx, const int *secy, const Coordinates &x, const Coordinates &y,
                        int xlen, int ylen, int *y2x0, int *y2x, const double D0_MIN, double d0, float * mem)
{
    //create score matrix for DP
    score_matrix_rmsd_sec(r1, r2, score, secx, secy, x, y, xlen, ylen,
                          y2x0, D0_MIN,d0, mem);

    float gap_open=-1.0;
    NWDP_TM(score, path, val, xlen, ylen, gap_open, y2x);
}

void find_max_frag(const Coordinates &x, int len, int *start_max,
                   int *end_max, float dcu0, const bool fast_opt)
{
    int r_min, fra_min=4;           //minimum fragment for search
    if (fast_opt) fra_min=8;
    int start;
    int Lfr_max=0;

    r_min= (int) (len*1.0/3.0); //minimum fragment, in case too small protein
    if(r_min > fra_min) r_min=fra_min;

    int inc=0;
    float dcu0_cut=dcu0*dcu0;;
    float dcu_cut=dcu0_cut;

    while(Lfr_max < r_min)
    {
        Lfr_max=0;
        int j=1;    //number of residues at nf-fragment
        start=0;
        for(int i=1; i<len; i++)
        {
            if(dist(x.x[i-1], x.y[i-1], x.z[i-1], x.x[i], x.y[i], x.z[i]) < dcu_cut)
            {
                j++;

                if(i==(len-1))
                {
                    if(j > Lfr_max)
                    {
                        Lfr_max=j;
                        *start_max=start;
                        *end_max=i;
                    }
                    j=1;
                }
            }
            else
            {
                if(j>Lfr_max)
                {
                    Lfr_max=j;
                    *start_max=start;
                    *end_max=i-1;
                }

                j=1;
                start=i;
            }
        }// for i;

        if(Lfr_max < r_min)
        {
            inc++;
            double dinc=pow(1.1, (double) inc) * dcu0;
            dcu_cut= dinc*dinc;
        }
    }//while <;
}

//perform fragment gapless threading to find the best initial alignment
//input: x, y, xlen, ylen
//output: y2x0 stores the best alignment: e.g.,
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
double get_initial_fgt(Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                       const Coordinates &x, const Coordinates &y, int xlen, int ylen,
                       int *y2x, float d0, float d0_search,
                       float dcu0, const bool fast_opt, float t[3], float u[3][3],
                       float * mem)
{
    int fra_min=4;           //minimum fragment for search
    if (fast_opt) fra_min=8;
    int fra_min1=fra_min-1;  //cutoff for shift, save time

    int xstart=0, ystart=0, xend=0, yend=0;

    find_max_frag(x, xlen, &xstart, &xend, dcu0, fast_opt);
    find_max_frag(y, ylen, &ystart, &yend, dcu0, fast_opt);


    int Lx = xend-xstart+1;
    int Ly = yend-ystart+1;
    int *ifr, *y2x_;
    int L_fr=getmin(Lx, Ly);
    ifr= new int[L_fr];
    y2x_= new int[ylen+1];

    //select what piece will be used (this may araise ansysmetry, but
    //only when L1=L2 and Lfr1=Lfr2 and L1 ne Lfr1
    //if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1

    if(Lx<Ly || (Lx==Ly && xlen<=ylen))
    {
        for(int i=0; i<L_fr; i++) ifr[i]=xstart+i;
    }
    else if(Lx>Ly || (Lx==Ly && xlen>ylen))
    {
        for(int i=0; i<L_fr; i++) ifr[i]=ystart+i;
    }


    int L0=getmin(xlen, ylen); //non-redundant to get_initial1
    if(L_fr==L0)
    {
        int n1= (int)(L0*0.1); //my index starts from 0
        int n2= (int)(L0*0.89);

        int j=0;
        for(int i=n1; i<= n2; i++)
        {
            ifr[j]=ifr[i];
            j++;
        }
        L_fr=j;
    }


    //gapless threading for the extracted fragment
    double tmscore, tmscore_max=-1;

    if(Lx<Ly || (Lx==Ly && xlen<=ylen))
    {
        int L1=L_fr;
        int min_len=getmin(L1, ylen);
        int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment
        if(min_ali<=fra_min1)  min_ali=fra_min1;
        int n1, n2;
        n1 = -ylen+min_ali;
        n2 = L1-min_ali;

        int i, j, k;
        for(k=n1; k<=n2; k+=(fast_opt)?3:1)
        {
            //get the map
            for(j=0; j<ylen; j++)
            {
                i=j+k;
                if(i>=0 && i<L1) y2x_[j]=ifr[i];
                else             y2x_[j]=-1;
            }

            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm, x, y, xlen, ylen, y2x_,
                                   d0, d0_search, t, u, mem);

            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }
    }
    else
    {
        int L2=L_fr;
        int min_len=getmin(xlen, L2);
        int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment
        if(min_ali<=fra_min1)  min_ali=fra_min1;
        int n1, n2;
        n1 = -L2+min_ali;
        n2 = xlen-min_ali;

        int i, j, k;

        for(k=n1; k<=n2; k++)
        {
            //get the map
            for(j=0; j<ylen; j++) y2x_[j]=-1;

            for(j=0; j<L2; j++)
            {
                i=j+k;
                if(i>=0 && i<xlen) y2x_[ifr[j]]=i;
            }

            //evaluate the map quickly in three iterations
            tmscore=get_score_fast(r1, r2, xtm, ytm,
                                   x, y, xlen, ylen, y2x_, d0,d0_search, t, u, mem);
            if(tmscore>=tmscore_max)
            {
                tmscore_max=tmscore;
                for(j=0; j<ylen; j++) y2x[j]=y2x_[j];
            }
        }
    }


    delete [] ifr;
    delete [] y2x_;
    return tmscore_max;
}

//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap
double DP_iter(Coordinates &r1, Coordinates &r2,
               Coordinates &xtm, Coordinates &ytm,
               Coordinates &xt, float **score, bool **path, float **val,
               const Coordinates &x, const Coordinates &y, int xlen, int ylen, float t[3], float u[3][3],
               int invmap0[], int g1, int g2, int iteration_max, float local_d0_search,
               float D0_MIN, float Lnorm, float d0, float score_d8, float * mem)
{
    float gap_open[2]={-0.6, 0};
    float rmsd;
    int *invmap=new int[ylen+1];

    int iteration, i, j, k;
    float tmscore, tmscore_max, tmscore_old=0;
    int score_sum_method=8, simplify_step=40;
    tmscore_max=-1;

    //double d01=d0+1.5;
    float d02=d0*d0;
    for(int g=g1; g<g2; g++)
    {
        for(iteration=0; iteration<iteration_max; iteration++)
        {
            NWDP_TM(score, path, val, x, y, xlen, ylen,
                    t, u, d02, gap_open[g], invmap, mem);

            k=0;
            for(j=0; j<ylen; j++)
            {
                i=invmap[j];

                if(i>=0) //aligned
                {
                    xtm.x[k]=x.x[i];
                    xtm.y[k]=x.y[i];
                    xtm.z[k]=x.z[i];

                    ytm.x[k]=y.x[j];
                    ytm.y[k]=y.y[j];
                    ytm.z[k]=y.z[j];
                    k++;
                }
            }

            tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
                                      simplify_step, score_sum_method, &rmsd, local_d0_search,
                                      Lnorm, score_d8, d0, mem);


            if(tmscore>tmscore_max)
            {
                tmscore_max=tmscore;
                for(i=0; i<ylen; i++) invmap0[i]=invmap[i];
            }

            if(iteration>0)
            {
                if(fabs(tmscore_old-tmscore)<0.000001) break;
            }
            tmscore_old=tmscore;
        }// for iteration

    }//for gapopen


    delete []invmap;
    return tmscore_max;
}


void output_superpose(const char *xname, const char *fname_super,
                      float t[3], float u[3][3], const int ter_opt=3)
{
    ifstream fin(xname);
    stringstream buf;
    string line;
    float x[3];  // before transform
    float x1[3]; // after transform
    if (fin.is_open())
    {
        while (fin.good())
        {
            getline(fin, line);
            if (line.compare(0, 6, "ATOM  ")==0 ||
                line.compare(0, 6, "HETATM")==0)
            {
                x[0]=atof(line.substr(30,8).c_str());
                x[1]=atof(line.substr(38,8).c_str());
                x[2]=atof(line.substr(46,8).c_str());
                transform(t, u, x, x1);
                buf<<line.substr(0,30)<<setiosflags(ios::fixed)
                   <<setprecision(3)
                   <<setw(8)<<x1[0] <<setw(8)<<x1[1] <<setw(8)<<x1[2]
                   <<line.substr(54)<<'\n';
            }
            else if (line.size())
                buf<<line<<'\n';
            if (ter_opt>=1 && line.compare(0,3,"END")==0) break;
        }
        fin.close();
    }
    else
    {
        char message[5000];
        sprintf(message, "Can not open file: %s\n", xname);
        PrintErrorAndQuit(message);
    }

    ofstream fp(fname_super);
    fp<<buf.str();
    fp.close();
    buf.str(string()); // clear stream
}

/* extract rotation matrix based on TMscore8 */
void output_rotation_matrix(const char* fname_matrix,
                            const float t[3], const float u[3][3])
{
    fstream fout;
    fout.open(fname_matrix, ios::out | ios::trunc);
    if (fout)// succeed
    {
        fout << "------ The rotation matrix to rotate Chain_1 to Chain_2 ------\n";
        char dest[1000];
        sprintf(dest, "m %18s %14s %14s %14s\n", "t[m]", "u[m][0]", "u[m][1]", "u[m][2]");
        fout << string(dest);
        for (int k = 0; k < 3; k++)
        {
            sprintf(dest, "%d %18.10f %14.10f %14.10f %14.10f\n", k, t[k], u[k][0], u[k][1], u[k][2]);
            fout << string(dest);
        }
        fout << "\nCode for rotating Structure A from (x,y,z) to (X,Y,Z):\n"
                "for(i=0; i<L; i++)\n"
                "{\n"
                "   X[i] = t[0] + u[0][0]*x[i] + u[0][1]*y[i] + u[0][2]*z[i]\n"
                "   Y[i] = t[1] + u[1][0]*x[i] + u[1][1]*y[i] + u[1][2]*z[i]\n"
                "   Z[i] = t[2] + u[2][0]*x[i] + u[2][1]*y[i] + u[2][2]*z[i]\n"
                "}\n";
        fout.close();
    }
    else
        cout << "Open file to output rotation matrix fail.\n";
}

//output the final results
void output_results(
        const char *xname, const char *yname,
        const char *chainID1, const char *chainID2,
        const int xlen, const int ylen, float t[3], float u[3][3],
        const float TM1, const float TM2,
        const float TM3, const float TM4, const float TM5,
        const float rmsd, const float d0_out,
        const char *seqM, const char *seqxA, const char *seqyA, const float Liden,
        const int n_ali8, const int n_ali, const int L_ali,
        const float TM_ali, const float rmsd_ali, const float TM_0,
        const float d0_0, const float d0A, const float d0B,
        const float Lnorm_ass, const float d0_scale,
        const float d0a, const float d0u, const char* fname_matrix,
        const int outfmt_opt, const int ter_opt, const char *fname_super,
        const bool i_opt, const bool I_opt, const bool a_opt,
        const bool u_opt, const bool d_opt)
{
    if (outfmt_opt<=0)
    {
        printf("\nName of Chain_1: %s%s (to be superimposed onto Chain_2)\n",
               xname,chainID1);
        printf("Name of Chain_2: %s%s\n", yname,chainID2);
        printf("Length of Chain_1: %d residues\n", xlen);
        printf("Length of Chain_2: %d residues\n\n", ylen);

        if (i_opt || I_opt)
            printf("User-specified initial alignment: TM/Lali/rmsd = %7.5lf, %4d, %6.3lf\n", TM_ali, L_ali, rmsd_ali);

        printf("Aligned length= %d, RMSD= %6.2f, Seq_ID=n_identical/n_aligned= %4.3f\n", n_ali8, rmsd, Liden/(n_ali8+0.00000001));
        printf("TM-score= %6.5f (if normalized by length of Chain_1, i.e., LN=%d, d0=%.2f)\n", TM2, xlen, d0B);
        printf("TM-score= %6.5f (if normalized by length of Chain_2, i.e., LN=%d, d0=%.2f)\n", TM1, ylen, d0A);

        if (a_opt)
            printf("TM-score= %6.5f (if normalized by average length of two structures, i.e., LN= %.1f, d0= %.2f)\n", TM3, (xlen+ylen)*0.5, d0a);
        if (u_opt)
            printf("TM-score= %6.5f (if normalized by user-specified LN=%.2f and d0=%.2f)\n", TM4, Lnorm_ass, d0u);
        if (d_opt)
            printf("TM-score= %6.5f (if scaled by user-specified d0= %.2f, and LN= %d)\n", TM5, d0_scale, ylen);
        printf("(You should use TM-score normalized by length of the reference protein)\n");

        //output alignment
        printf("\n(\":\" denotes residue pairs of d < %4.1f Angstrom, ", d0_out);
        printf("\".\" denotes other aligned residues)\n");
        printf("%s\n", seqxA);
        printf("%s\n", seqM);
        printf("%s\n", seqyA);
    }
    else if (outfmt_opt==1)
    {
        printf(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
               xname, chainID1, xlen, d0B, Liden/xlen, TM2);
        printf("%s\n", seqxA);
        printf(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
               yname, chainID2, ylen, d0A, Liden/ylen, TM1);
        printf("%s\n", seqyA);

        printf("# Lali=%d\tRMSD=%.2f\tseqID_ali=%.3f\n",
               n_ali8, rmsd, Liden/(n_ali8+0.00000001));

        if (i_opt || I_opt)
            printf("# User-specified initial alignment: TM=%.5lf\tLali=%4d\trmsd=%.3lf\n", TM_ali, L_ali, rmsd_ali);

        if(a_opt)
            printf("# TM-score=%.5f (normalized by average length of two structures: L=%.1f\td0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);

        if(u_opt)
            printf("# TM-score=%.5f (normalized by user-specified L=%.2f\td0=%.2f)\n", TM4, Lnorm_ass, d0u);

        if(d_opt)
            printf("# TM-score=%.5f (scaled by user-specified d0=%.2f\tL=%d)\n", TM5, d0_scale, ylen);

        printf("$$$$\n");
    }
    else if (outfmt_opt==2)
    {
        printf("%s%s\t%s%s\t%.4f\t%.4f\t%.2f\t%.3f\t%4.3f\t%4.3f\t%d\t%d\t%d",
               xname, chainID1, yname, chainID2, TM2, TM1, rmsd,
               Liden/xlen, Liden/ylen, Liden/(n_ali8+0.00000001),
               xlen, ylen, n_ali8);
    }
    cout << endl;

    if (strlen(fname_matrix))
        output_rotation_matrix(fname_matrix, t, u);
    if (strlen(fname_super))
        output_superpose(xname, fname_super, t, u, ter_opt);
}

double standard_TMscore(Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                        Coordinates &xt, Coordinates &x, Coordinates &y, int xlen, int ylen, int invmap[],
                        int& L_ali, float& RMSD, float D0_MIN, float Lnorm, float d0,
                        float d0_search, float score_d8, float t[3], float u[3][3], float * mem)
{
    D0_MIN = 0.5;
    Lnorm = ylen;
    if (Lnorm > 21)
        d0 = (1.24*pow((Lnorm*1.0 - 15), 1.0 / 3) - 1.8);
    else
        d0 = D0_MIN;
    if (d0 < D0_MIN)
        d0 = D0_MIN;
    double d0_input = d0;// Scaled by seq_min

    double tmscore;// collected alined residues from invmap
    int n_al = 0;
    int i;
    for (int j = 0; j<ylen; j++)
    {
        i = invmap[j];
        if (i >= 0)
        {
            xtm.x[n_al] = x.x[i];
            xtm.y[n_al] = x.y[i];
            xtm.z[n_al] = x.z[i];

            ytm.x[n_al] = y.x[j];
            ytm.y[n_al] = y.y[j];
            ytm.z[n_al] = y.z[j];

            r1.x[n_al] = x.x[i];
            r1.y[n_al] = x.y[i];
            r1.z[n_al] = x.z[i];

            r2.x[n_al] = y.x[j];
            r2.y[n_al] = y.y[j];
            r2.z[n_al] = y.z[j];
            n_al++;
        }
        else if (i != -1) PrintErrorAndQuit("Wrong map!\n");
    }
    L_ali = n_al;

    KabschFast(r1, r2, n_al, 0, &RMSD, t, u, mem);
    RMSD = sqrt( RMSD/(1.0*n_al) );

    int temp_simplify_step = 1;
    int temp_score_sum_method = 0;
    d0_search = d0_input;
    float rms = 0.0;
    tmscore = TMscore8_search_standard(r1, r2, xtm, ytm, xt, n_al, t, u,
                                       temp_simplify_step, temp_score_sum_method, &rms, d0_input,
                                       score_d8, d0, mem);
    tmscore = tmscore * n_al / (1.0*Lnorm);

    return tmscore;
}

/* entry function for TMalign */
int TMalign_main(
        const Coordinates &xa, const Coordinates &ya,
        const char *seqx, const char *seqy, const int *secx, const int *secy,
        float t0[3], float u0[3][3],
        float &TM1, float &TM2, float &TM3, float &TM4, float &TM5,
        float &d0_0, float &TM_0,
        float &d0A, float &d0B, float &d0u, float &d0a, float &d0_out,
        string &seqM, string &seqxA, string &seqyA,
        float &rmsd0, int &L_ali, float &Liden,
        float &TM_ali, float &rmsd_ali, int &n_ali, int &n_ali8,
        const int xlen, const int ylen, const float Lnorm_ass,
        const float d0_scale,
        const bool i_opt, const bool I_opt, const bool a_opt,
        const bool u_opt, const bool d_opt, const bool fast_opt, float * mem)
{
    int minlen = min(xlen, ylen);

    float D0_MIN;        //for d0
    float Lnorm;         //normalization length
    float score_d8,d0,d0_search,dcu0;//for TMscore search
    float t[3], u[3][3]; //Kabsch translation vector and rotation matrix
    float **score;     // Input score table for dynamic programming
    bool   **path;     // for dynamic programming
    float **val;       // for dynamic programming
    Coordinates xtm(minlen);
    Coordinates ytm(minlen);
    Coordinates xt(xlen);
    Coordinates r1(minlen);
    Coordinates r2(minlen);


    /***********************/
    /* allocate memory     */
    /***********************/
    NewArray(&score, xlen+1, ylen+1);
    NewArray(&path, xlen+1, ylen+1);
    NewArray(&val, xlen+1, ylen+1);

//    NewArray(&xtm, minlen, 3);
//    NewArray(&ytm, minlen, 3);
//    NewArray(&xt, xlen, 3);
//    NewArray(&r1, minlen, 3);
//    NewArray(&r2, minlen, 3);

    /***********************/
    /*    parameter set    */
    /***********************/
    parameter_set4search(xlen, ylen, D0_MIN, Lnorm,
                         score_d8, d0, d0_search, dcu0);
    int simplify_step    = 40; //for similified search engine
    int score_sum_method = 8;  //for scoring method, whether only sum over pairs with dis<score_d8

    int i;
    int *invmap0         = new int[ylen+1];
    int *invmap          = new int[ylen+1];
    float TM, TMmax=-1;
    for(i=0; i<ylen; i++) invmap0[i]=-1;

    float ddcc=0.4;
    if (Lnorm <= 40) ddcc=0.1;   //Lnorm was setted in parameter_set4search
    float local_d0_search = d0_search;


    /******************************************************/
    /*    get initial alignment with gapless threading    */
    /******************************************************/

    get_initial(r1, r2, xtm, ytm, xa, ya, xlen, ylen, invmap0, d0,
                d0_search, fast_opt, t, u, mem);
    TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap0,
                         t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
                         score_d8, d0, mem);
    if (TM>TMmax) TMmax = TM;
    //run dynamic programing iteratively to find the best alignment
    TM = DP_iter( r1, r2, xtm, ytm, xt, score, path, val, xa, ya,
                  xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30, local_d0_search,
                  D0_MIN, Lnorm, d0, score_d8, mem);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }


    /************************************************************/
    /*    get initial alignment based on secondary structure    */
    /************************************************************/
    get_initial_ss(score, path, val, secx, secy, xlen, ylen, invmap);
    TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap,
                         t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
                         score_d8, d0, mem);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }
    if (TM > TMmax*0.2)
    {
        TM = DP_iter(r1, r2, xtm, ytm, xt, score, path, val, xa, ya,
                     xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30,
                     local_d0_search, D0_MIN, Lnorm, d0, score_d8, mem);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
    }


    /************************************************************/
    /*    get initial alignment based on local superposition    */
    /************************************************************/
    //=initial5 in original TM-align
    if (get_initial5( r1, r2, xtm, ytm, score, path, val, xa, ya,
                      xlen, ylen, invmap, d0, d0_search, fast_opt, D0_MIN, mem))
    {
        TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
                             invmap, t, u, simplify_step, score_sum_method,
                             local_d0_search, Lnorm, score_d8, d0, mem);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
        if (TM > TMmax*ddcc)
        {
            TM = DP_iter(r1, r2, xtm, ytm, xt, score, path, val, xa, ya,
                         xlen, ylen, t, u, invmap, 0, 2, 2, local_d0_search,
                         D0_MIN, Lnorm, d0, score_d8, mem);
            if (TM>TMmax)
            {
                TMmax = TM;
                for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
            }
        }
    }
    else
        cerr << "\n\nWarning: initial alignment from local superposition fail!\n\n" << endl;


    /********************************************************************/
    /* get initial alignment by local superposition+secondary structure */
    /********************************************************************/
    //=initial3 in original TM-align
    get_initial_ssplus(r1, r2, score, path, val, secx, secy, xa, ya,
                       xlen, ylen, invmap0, invmap, D0_MIN, d0, mem);
    TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap,
                         t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
                         score_d8, d0, mem);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }
    if (TM > TMmax*ddcc)
    {
        TM = DP_iter(r1, r2, xtm, ytm, xt, score, path, val, xa, ya,
                     xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30,
                     local_d0_search, D0_MIN, Lnorm, d0, score_d8, mem);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
    }


    /*******************************************************************/
    /*    get initial alignment based on fragment gapless threading    */
    /*******************************************************************/
    //=initial4 in original TM-align
    //TODO
    get_initial_fgt(r1, r2, xtm, ytm, xa, ya, xlen, ylen,
                    invmap, d0, d0_search, dcu0, fast_opt, t, u, mem);
    TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap,
                         t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
                         score_d8, d0, mem);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }
    if (TM > TMmax*ddcc)
    {
        TM = DP_iter(r1, r2, xtm, ytm, xt, score, path, val, xa, ya,
                     xlen, ylen, t, u, invmap, 1, 2, 2, local_d0_search, D0_MIN,
                     Lnorm, d0, score_d8, mem);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
    }

    //*******************************************************************//
    //    The alignment will not be changed any more in the following    //
    //*******************************************************************//
    //check if the initial alignment is generated approately
    bool flag=false;
    for(i=0; i<ylen; i++)
    {
        if(invmap0[i]>=0)
        {
            flag=true;
            break;
        }
    }
    if(!flag)
    {
        cout << "There is no alignment between the two proteins!" << endl;
        cout << "Program stop with no result!" << endl;
        return 1;
    }


    //********************************************************************//
    //    Detailed TMscore search engine --> prepare for final TMscore    //
    //********************************************************************//
    //run detailed TMscore search engine for the best alignment, and
    //extract the best rotation matrix (t, u) for the best alginment
    simplify_step=1;
    if (fast_opt) simplify_step=40;
    score_sum_method=8;
    TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
                                  invmap0, t, u, simplify_step, score_sum_method, local_d0_search,
                                  false, Lnorm, score_d8, d0, mem);

    //select pairs with dis<d8 for final TMscore computation and output alignment
    int k=0;
    int *m1, *m2;
    double d;
    m1=new int[xlen]; //alignd index in x
    m2=new int[ylen]; //alignd index in y
    do_rotation(xa, xt, xlen, t, u);
    k=0;
    for(int j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
            n_ali++;
            d =sqrt( dist(xt.x[i], xt.y[i], xt.z[i], ya.x[j], ya.y[j], ya.z[j]));

            if (d <= score_d8 || (I_opt == true))
            {
                m1[k]=i;
                m2[k]=j;

                xtm.x[k]=xa.x[i];
                xtm.y[k]=xa.y[i];
                xtm.z[k]=xa.z[i];
                ytm.x[k]=ya.x[j];
                ytm.y[k]=ya.y[j];
                ytm.z[k]=ya.z[j];
                r1.x[k] = xt.x[i];
                r1.y[k] = xt.y[i];
                r1.z[k] = xt.z[i];
                r2.x[k] = ya.x[j];
                r2.y[k] = ya.y[j];
                r2.z[k] = ya.z[j];
                k++;
            }
        }
    }
    n_ali8=k;

    KabschFast(r1, r2, n_ali8, 0, &rmsd0, t, u, mem);// rmsd0 is used for final output, only recalculate rmsd0, not t & u
    rmsd0 = sqrt(rmsd0 / n_ali8);


    //****************************************//
    //              Final TMscore             //
    //    Please set parameters for output    //
    //****************************************//
    float rmsd;
    simplify_step=1;
    score_sum_method=0;
    float Lnorm_0=ylen;


    //normalized by length of structure A
    parameter_set4final(Lnorm_0, D0_MIN, Lnorm,
                        score_d8, d0, d0_search, dcu0);
    d0A=d0;
    d0_0=d0A;
    local_d0_search = d0_search;
    TM1 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0, simplify_step,
                          score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0, mem);
    TM_0 = TM1;

    //normalized by length of structure B
    parameter_set4final(xlen+0.0, D0_MIN, Lnorm,
                        score_d8, d0, d0_search, dcu0);
    d0B=d0;
    local_d0_search = d0_search;
    TM2 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t, u, simplify_step,
                          score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0, mem);

    double Lnorm_d0;
    if (a_opt)
    {
        //normalized by average length of structures A, B
        Lnorm_0=(xlen+ylen)*0.5;
        parameter_set4final(Lnorm_0, D0_MIN, Lnorm,
                            score_d8, d0, d0_search, dcu0);
        d0a=d0;
        d0_0=d0a;
        local_d0_search = d0_search;

        TM3 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
                              simplify_step, score_sum_method, &rmsd, local_d0_search, Lnorm,
                              score_d8, d0, mem);
        TM_0=TM3;
    }
    if (u_opt)
    {
        //normalized by user assigned length
        parameter_set4final(Lnorm_ass, D0_MIN, Lnorm,
                            score_d8, d0, d0_search, dcu0);
        d0u=d0;
        d0_0=d0u;
        Lnorm_0=Lnorm_ass;
        local_d0_search = d0_search;
        TM4 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
                              simplify_step, score_sum_method, &rmsd, local_d0_search, Lnorm,
                              score_d8, d0, mem);
        TM_0=TM4;
    }
    if (d_opt)
    {
        //scaled by user assigned d0
        parameter_set4scale(ylen, d0_scale, Lnorm,
                            score_d8, d0, d0_search, dcu0);
        d0_out=d0_scale;
        d0_0=d0_scale;
        //Lnorm_0=ylen;
        Lnorm_d0=Lnorm_0;
        local_d0_search = d0_search;
        TM5 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
                              simplify_step, score_sum_method, &rmsd, local_d0_search, Lnorm,
                              score_d8, d0, mem);
        TM_0=TM5;
    }

    /* derive alignment from superposition */
    int ali_len=xlen+ylen; //maximum length of alignment
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');

    do_rotation(xa, xt, xlen, t, u);

    int kk=0, i_old=0, j_old=0;
    d=0;
    for(int k=0; k<n_ali8; k++)
    {
        for(int i=i_old; i<m1[k]; i++)
        {
            //align x to gap
            seqxA[kk]=seqx[i];
            seqyA[kk]='-';
            seqM[kk]=' ';
            kk++;
        }

        for(int j=j_old; j<m2[k]; j++)
        {
            //align y to gap
            seqxA[kk]='-';
            seqyA[kk]=seqy[j];
            seqM[kk]=' ';
            kk++;
        }

        seqxA[kk]=seqx[m1[k]];
        seqyA[kk]=seqy[m2[k]];
        Liden+=(seqxA[kk]==seqyA[kk]);
        d=sqrt(dist(xt.x[m1[k]], xt.y[m1[k]], xt.z[m1[k]], ya.x[m2[k]], ya.y[m2[k]], ya.z[m2[k]]));

        if(d<d0_out) seqM[kk]=':';
        else         seqM[kk]='.';
        kk++;
        i_old=m1[k]+1;
        j_old=m2[k]+1;
    }

    //tail
    for(int i=i_old; i<xlen; i++)
    {
        //align x to gap
        seqxA[kk]=seqx[i];
        seqyA[kk]='-';
        seqM[kk]=' ';
        kk++;
    }
    for(int j=j_old; j<ylen; j++)
    {
        //align y to gap
        seqxA[kk]='-';
        seqyA[kk]=seqy[j];
        seqM[kk]=' ';
        kk++;
    }
    seqxA=seqxA.substr(0,kk);
    seqyA=seqyA.substr(0,kk);
    seqM =seqM.substr(0,kk);

    /* free memory */
    delete [] invmap0;
    delete [] invmap;
    delete [] m1;
    delete [] m2;
    DeleteArray(&score, xlen+1);
    DeleteArray(&path, xlen+1);
    DeleteArray(&val, xlen+1);

    free(xtm.x);
    free(xtm.y);
    free(xtm.z);

    free(ytm.x);
    free(ytm.y);
    free(ytm.z);

    free(xt.x);
    free(xt.y);
    free(xt.z);

    free(r1.x);
    free(r1.y);
    free(r1.z);

    free(r2.x);
    free(r2.y);
    free(r2.z);


//    DeleteArray(&xtm, minlen);
//    DeleteArray(&ytm, minlen);
//    DeleteArray(&xt, xlen);
//    DeleteArray(&r1, minlen);
//    DeleteArray(&r2, minlen);
    return 0; // zero for no exception
}
