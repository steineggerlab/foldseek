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
#include "Coordinates.h"
#include "simd.h"
#include <string.h>
#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde/x86/avx2.h>

float coords_sum_ssq_xyz_sse2(int nat, float *x, float *y,float *z,float center[3]){
    int lower_nat=(nat/4)*4;
    float invfnat=1.0/(float)nat;
    int i=0;
    float sums[4] __attribute__ ((aligned (16)));
    {
        __m128 sumx = _mm_setzero_ps();
        __m128 sumy = _mm_setzero_ps();
        __m128 sumz = _mm_setzero_ps();
        __m128 ssq  = _mm_setzero_ps();
        for(;i<lower_nat;i+=4){
            __m128 p0 = _mm_load_ps(&(x[i]));
            sumx= _mm_add_ps(sumx,p0);
            ssq= _mm_add_ps(ssq,_mm_mul_ps(p0,p0));
            __m128 p1 = _mm_load_ps(&(y[i]));
            sumy = _mm_add_ps(sumy,p1);
            ssq = _mm_add_ps(ssq,_mm_mul_ps(p1,p1));
            __m128 p2 = _mm_load_ps(&(z[i]));
            sumz = _mm_add_ps(sumz,p2);
            ssq = _mm_add_ps(ssq,_mm_mul_ps(p2,p2));
        }
        __m128 t1  = _mm_add_ps(_mm_unpacklo_ps(sumx,sumz),_mm_unpackhi_ps(sumx,sumz));
        __m128 t2  = _mm_add_ps(_mm_unpacklo_ps(sumy,ssq),_mm_unpackhi_ps(sumy,ssq));
        __m128 sum0= _mm_add_ps(_mm_unpacklo_ps(t1,t2),_mm_unpackhi_ps(t1,t2));
        _mm_store_ps(sums,sum0);
    }

    for(;i<nat;i++){
        sums[0]+=x[i];
        sums[3]+=x[i]*x[i];
        sums[1]+=y[i];
        sums[3]+=y[i]*y[i];
        sums[2]+=z[i];
        sums[3]+=z[i]*z[i];;
    }
    for(int i=0;i<3;i++)
        center[i]=sums[i]*invfnat;
    return(sums[3]-center[0]*sums[0]-center[1]*sums[1]-center[2]*sums[2]);

}

#define SQRT3 1.732050807568877

void R34v4_sse2 (float r[12],float *x,float *Rx){
    __m128 r0 = _mm_load_ps(r);
    __m128 r1 = _mm_load_ps(&(r[4]));
    __m128 r2 = _mm_load_ps(&(r[8]));
    __m128 one= _mm_set_ss(1.0f);
    __m128 mu = _mm_load_ps(x);
    __m128 m0 = _mm_mul_ps(r0,mu);
    __m128 m1 = _mm_mul_ps(r1,mu);
    __m128 m2 = _mm_mul_ps(r2,mu);

    __m128 s1 =_mm_add_ps(_mm_unpacklo_ps(m0,m2),_mm_unpackhi_ps(m0,m2));
    __m128 s2 =_mm_add_ps(_mm_unpacklo_ps(m1,one),_mm_unpackhi_ps(m1,one));
    _mm_store_ps(Rx,_mm_add_ps(_mm_unpacklo_ps(s1,s2),_mm_unpackhi_ps(s1,s2)));
}

template <class T> void rmatrix (T ev,T r[3][3],T u[3][3]){
    //calculate rotation matrix

    T a00=(r[0][0]+r[1][1]+r[2][2]);
    T a01=(r[1][2]-r[2][1]);
    T a02=(r[2][0]-r[0][2]);
    T a03=(r[0][1]-r[1][0]);
    T a11=(r[0][0]-r[1][1]-r[2][2]);
    T a12=(r[0][1]+r[1][0]);
    T a13=(r[2][0]+r[0][2]);
    T a22=(-r[0][0]+r[1][1]-r[2][2]);
    T a23=(r[1][2]+r[2][1]);
    T a33=(-r[0][0]-r[1][1]+r[2][2]);

    //from Theobald
    a00-=ev;a11-=ev;a22-=ev;a33-=ev;
    T a2233_3223 = a22 * a33 - a23 * a23;
    T a1233_3123 = a12 * a33-a13*a23;
    T a1223_3122 = a12 * a23 - a13 * a22;
    T a0232_3022 = a02 * a23-a03*a22;
    T a0233_3023 = a02 * a33 - a03 * a23;
    T a0231_3021 = a02 * a13-a03*a12;

    T q[4]={a11*a2233_3223-a12*a1233_3123+a13*a1223_3122, -a01*a2233_3223+a12*a0233_3023-a13*a0232_3022,a01*a1233_3123-a11*a0233_3023+a13*a0231_3021,-a01*a1223_3122+a11*a0232_3022-a12*a0231_3021};

    T invlen2q=1.0f/(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
    T aj=q[0]*q[0]*invlen2q;
    T xj=q[1]*q[1]*invlen2q;
    T yj=q[2]*q[2]*invlen2q;
    T zj=q[3]*q[3]*invlen2q;
    T  xy = q[1] * q[2]*invlen2q;
    T  az = q[0] * q[3]*invlen2q;
    T  zx = q[3] * q[1]*invlen2q;
    T  ay = q[0] * q[2]*invlen2q;
    T  yz = q[2] * q[3]*invlen2q;
    T  ax = q[0] * q[1]*invlen2q;

    u[0][0]= aj + xj - yj - zj; u[0][1]= 2.0f * (xy + az); u[0][2]= 2.0f * (zx - ay);
    u[1][0]= 2.0f * (xy - az);  u[1][1]=aj - xj + yj - zj; u[1][2]= 2.0f * (yz + ax);
    u[2][0]= 2.0f * (zx + ay),  u[2][1]= 2.0f * (yz - ax); u[2][2]= aj - xj - yj + zj;
}
template <class T> void rmatrix (T ev,T r[9],T u[3][3]){
    //calculate rotation matrix

    T a00=(r[0]+r[4]+r[8]);
    T a01=(r[5]-r[7]);
    T a02=(r[6]-r[2]);
    T a03=(r[1]-r[3]);
    T a11=(r[0]-r[4]-r[8]);
    T a12=(r[1]+r[3]);
    T a13=(r[6]+r[2]);
    T a22=(-r[0]+r[4]-r[8]);
    T a23=(r[5]+r[7]);
    T a33=(-r[0]-r[4]+r[8]);

    //from Theobald
    a00-=ev;a11-=ev;a22-=ev;a33-=ev;
    T a2233_3223 = a22 * a33 - a23 * a23;
    T a1233_3123 = a12 * a33-a13*a23;
    T a1223_3122 = a12 * a23 - a13 * a22;
    T a0232_3022 = a02 * a23-a03*a22;
    T a0233_3023 = a02 * a33 - a03 * a23;
    T a0231_3021 = a02 * a13-a03*a12;

    T q[4]={a11*a2233_3223-a12*a1233_3123+a13*a1223_3122, -a01*a2233_3223+a12*a0233_3023-a13*a0232_3022,a01*a1233_3123-a11*a0233_3023+a13*a0231_3021,-a01*a1223_3122+a11*a0232_3022-a12*a0231_3021};

    T len2q=q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3];
    if(!len2q){
        //return the identity matrix
        u[0][0]= 1.0f; u[0][1]= 0.0f; u[0][2]= 0.0f;
        u[1][0]= 0.0f; u[1][1]= 1.0f; u[1][2]= 0.0f;
        u[2][0]= 0.0f; u[2][1]= 0.0f; u[2][2]= 1.0f;
    }
    else{
        T invlen2q=1.0/len2q;
        T aj=q[0]*q[0]*invlen2q;
        T xj=q[1]*q[1]*invlen2q;
        T yj=q[2]*q[2]*invlen2q;
        T zj=q[3]*q[3]*invlen2q;
        T  xy = q[1] * q[2]*invlen2q;
        T  az = q[0] * q[3]*invlen2q;
        T  zx = q[3] * q[1]*invlen2q;
        T  ay = q[0] * q[2]*invlen2q;
        T  yz = q[2] * q[3]*invlen2q;
        T  ax = q[0] * q[1]*invlen2q;

        u[0][0]= aj + xj - yj - zj; u[0][1]= 2.0 * (xy + az); u[0][2]= 2.0 * (zx - ay);
        u[1][0]= 2.0 * (xy - az);  u[1][1]=aj - xj + yj - zj; u[1][2]= 2.0 * (yz + ax);
        u[2][0]= 2.0 * (zx + ay);  u[2][1]= 2.0 * (yz - ax); u[2][2]= aj - xj - yj + zj;
    }
}

float rmsd_sse2_matrix_xyz (int nat,float *c1x,float *c1y, float *c1z,float *c2x,float *c2y, float *c2z,float center1[3], float center2[3],double ssq,float u[3][3]){
    //c1 has (0,sum(x1),sum(y1),sum(z1));
    //c2 has (0,mean(x1),mean(y2),mean(z2));
    float fnat=(float)nat;
    float c1[4] __attribute__ ((aligned (16))) ={0.0f,center1[0]*fnat,center1[1]*fnat,center1[2]*fnat};
    float c2[4] __attribute__ ((aligned (16))) ={0.0f,center2[0],center2[1],center2[2]};
    float r0[4] __attribute__ ((aligned (16)));
    float r1[4] __attribute__ ((aligned (16)));
    float r2[4] __attribute__ ((aligned (16)));
    float rr[8] __attribute__ ((aligned (16)));
    int nat4=(nat%4)? nat/4+1 : nat/4;
    int padded_nat=nat4*4;
    //SSE3 block - not much faster than SSE2 on Phenom but might be better on Intel
    {
        __m128 mxx = _mm_setzero_ps();
        __m128 mxy = _mm_setzero_ps();
        __m128 mxz = _mm_setzero_ps();
        //these will be reused later - others go into block to be released
        {
            __m128 myx = _mm_setzero_ps();
            __m128 myy = _mm_setzero_ps();
            __m128 myz = _mm_setzero_ps();
            __m128 mzx = _mm_setzero_ps();
            __m128 mzy = _mm_setzero_ps();
            __m128 mzz = _mm_setzero_ps();


            for(int i=0;i<padded_nat;i+=4){
                //load the 4 sets of coords from molecule 2 and then load x,y,z of molecule 1
                __m128 mc2x=_mm_load_ps(&(c2x[i]));
                __m128 mc2y=_mm_load_ps(&(c2y[i]));
                __m128 mc2z=_mm_load_ps(&(c2z[i]));
                __m128 mc1 =_mm_load_ps(&(c1x[i]));

                //generate the 4 sets of products
                mxx=_mm_add_ps(mxx,_mm_mul_ps(mc1,mc2x));
                mxy=_mm_add_ps(mxy,_mm_mul_ps(mc1,mc2y));
                mxz=_mm_add_ps(mxz,_mm_mul_ps(mc1,mc2z));
                mc1=_mm_load_ps(&(c1y[i]));
                myx=_mm_add_ps(myx,_mm_mul_ps(mc1,mc2x));
                myy=_mm_add_ps(myy,_mm_mul_ps(mc1,mc2y));
                myz=_mm_add_ps(myz,_mm_mul_ps(mc1,mc2z));
                mc1=_mm_load_ps(&(c1z[i]));
                mzx=_mm_add_ps(mzx,_mm_mul_ps(mc1,mc2x));
                mzy=_mm_add_ps(mzy,_mm_mul_ps(mc1,mc2y));
                mzz=_mm_add_ps(mzz,_mm_mul_ps(mc1,mc2z));
            }
            //write out the components to the temp variables

            mzy = _mm_add_ps(_mm_unpacklo_ps(mzy,mxy),_mm_unpackhi_ps(mzy,mxy));
            mxy = _mm_add_ps(_mm_unpacklo_ps(mxx,mxz),_mm_unpackhi_ps(mxx,mxz));
            mxx = _mm_add_ps(_mm_unpacklo_ps(mzy,mxy),_mm_unpackhi_ps(mzy,mxy)); //mxx holds sums for zy,xx,xy,xz

            mzy = _mm_add_ps(_mm_unpacklo_ps(mzz,myy),_mm_unpackhi_ps(mzz,myy));
            mxy = _mm_add_ps(_mm_unpacklo_ps(myx,myz),_mm_unpackhi_ps(myx,myz));
            mxy = _mm_add_ps(_mm_unpacklo_ps(mzy,mxy),_mm_unpackhi_ps(mzy,mxy)); //mzy holds sums for zz,yx,yy,yz

            mxz = _mm_add_ps(mzx, _mm_movehl_ps(mzx, mzx));
            mxz = _mm_add_ps(mxz, _mm_unpacklo_ps(mxz,mxz)); //(junk,zx,junk,junk);
            __m128 mask = _mm_set_ss(1);
            mzy = _mm_setzero_ps();
            mask=_mm_cmpeq_ps(mask,mzy);
            mxz=_mm_and_ps(mask,mxz);
            mxz=_mm_movelh_ps(mxz,_mm_unpacklo_ps(mxx,mxy));//(0,zx,zy,zz);
            mxy=_mm_and_ps(mask,mxy);
            mxx=_mm_and_ps(mask,mxx);
        }//end block - only mxx,mxy,mxz remain - correspond to R matrix r0,r1,r2
        __m128 mc1=_mm_load_ps(c1);
        __m128 mc2=_mm_load_ps(c2);
        __m128 mt1=_mm_shuffle_ps(mc1,mc1,0x55);//s1x s1x s1x s1x
        mxx=_mm_sub_ps(mxx,_mm_mul_ps(mt1,mc2));
        mt1=_mm_shuffle_ps(mc1,mc1,0xAA);//s1y s1y s1y s1y
        mxy=_mm_sub_ps(mxy,_mm_mul_ps(mt1,mc2));
        mt1=_mm_shuffle_ps(mc1,mc1,0xFF);//s1y s1y s1y s1y
        mxz=_mm_sub_ps(mxz,_mm_mul_ps(mt1,mc2));
        //write out the matrix
        _mm_store_ps(r0,mxx); _mm_store_ps(r1,mxy); _mm_store_ps(r2,mxz);
        //calculate the determinant using triple product - do addition when calculating rest of dot products
        __m128 mdet = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(mxy, mxy, _MM_SHUFFLE(1,3,2,0)), _mm_shuffle_ps(mxz, mxz, _MM_SHUFFLE(2,1,3,0))),
                                 _mm_mul_ps(_mm_shuffle_ps(mxy, mxy, _MM_SHUFFLE(2,1,3,0)), _mm_shuffle_ps(mxz, mxz, _MM_SHUFFLE(1,3,2,0))));//cross_product
        mdet=_mm_mul_ps(mxx,mdet); //sum to get dot product

        //calculate the necessary 6 dot products - do additions in groups of 4 for sse2
        {
            __m128 mt0=_mm_mul_ps(mxx,mxx);
            __m128 mt1=_mm_mul_ps(mxx,mxy);
            __m128 mt2=_mm_mul_ps(mxy,mxy);
            __m128 mt3=_mm_mul_ps(mxx,mxz);

            mxx = _mm_add_ps(_mm_unpacklo_ps(mt0,mt2),_mm_unpackhi_ps(mt0,mt2));
            mt0 = _mm_add_ps(_mm_unpacklo_ps(mt1,mt3),_mm_unpackhi_ps(mt1,mt3));
            _mm_store_ps(rr,_mm_add_ps(_mm_unpacklo_ps(mxx,mt0),_mm_unpackhi_ps(mxx,mt0)));
        }
        mxx = _mm_mul_ps(mxz,mxy);
        mxy = _mm_mul_ps(mxz,mxz);

        mxx = _mm_add_ps(_mm_unpacklo_ps(mxx,mdet),_mm_unpackhi_ps(mxx,mdet));
        mxy = _mm_add_ps(_mm_unpacklo_ps(mxy,mdet),_mm_unpackhi_ps(mxy,mdet));
        _mm_store_ps(&(rr[4]),_mm_add_ps(_mm_unpacklo_ps(mxx,mxy),_mm_unpackhi_ps(mxx,mxy)));
    }
    //convert to double - can get by with floats except for g and h but the double is faster

    double detsq=rr[6]*rr[6];
    //lower triangular matrix rr
    double inv3=1.0/3.0;
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
    double d=(rr[6]<0)? e[0] + e[1] -e[2] : e[0] + e[1]+e[2];
    double rms=(ssq-d-d)/(double)nat;
    float r[3][3]={{r0[1],r0[2],r0[3]},
                   {r1[1],r1[2],r1[3]},
                   {r2[1],r2[2],r2[3]}};
    rmatrix((float)d,r,u);

    rms=(rms>1e-12)?sqrt(rms) : 0.0f;
    return(rms);
}

float kabsch_quat_soa_sse2(int nat, int *map, float *x1, float *y1,float *z1, float *x2,float *y2, float *z2,float *r){
    int nat4=(nat%4) ? nat/4+1 : nat/4;
    float center1[3],center2[3];
    float u[3][3];
    float *mem = (float*)mem_align(ALIGN_FLOAT,6*nat4*4*sizeof(float));

    memset(mem,0,nat4*24*sizeof(float));
    float* c1x= mem;
    float* c1y= &mem[nat4*4];
    float* c1z= &mem[nat4*8];
    float* c2x= &mem[nat4*12];
    float* c2y= &mem[nat4*16];
    float* c2z= &mem[nat4*20];

    for(int i=0;i<nat;i++){
        int n=(map == NULL) ? i : map[i];
        c1x[i]=x1[n];
        c1y[i]=y1[n];
        c1z[i]=z1[n];
        c2x[i]=x2[n];
        c2y[i]=y2[n];
        c2z[i]=z2[n];
    }
    float ssq=coords_sum_ssq_xyz_sse2(nat,c1x,c1y,c1z,center1);
    ssq+=coords_sum_ssq_xyz_sse2(nat,c2x,c2y,c2z,center2);
    float rms=rmsd_sse2_matrix_xyz(nat,c1x,c1y,c1z,c2x,c2y,c2z,center1,center2,(double)ssq,u);
    float rr[16] __attribute__ ((aligned (16)))=
            {-u[0][0],-u[1][0],-u[2][0],center2[0],
             -u[0][1],-u[1][1],-u[2][1],center2[1],
             -u[0][2],-u[1][2],-u[2][2],center2[2]};
    float w[4]__attribute__ ((aligned (16)));
    float v[4]__attribute__ ((aligned (16)))={center1[0],center1[1],center1[2],1.0f};
    R34v4_sse2(rr,v,w);
    r[0] =u[0][0]; r[1]=u[1][0];r[2] =u[2][0];r[3]=w[0];
    r[4] =u[0][1]; r[5]=u[1][1];r[6] =u[2][1];r[7]=w[1];
    r[8] =u[0][2]; r[9]=u[1][2];r[10]=u[2][2];r[11]=w[2];
    r[12]=0.0f;r[13]=0.0f;r[14]=0.0f;r[15]=1.0f;
    free (mem);
    return(rms);
}

float LG_score_soa_avx (float r[16],int nat, float *x1,float *y1,float *z1, float *x2, float *y2, float *z2, float *d,float invd0d0){//coords organized x4,y4,z4
    float mask_array[56]__attribute__ ((aligned (32)))=
            {1,0,0,0,0,0,0,0,
             1,1,0,0,0,0,0,0,
             1,1,1,0,0,0,0,0,
             1,1,1,1,0,0,0,0,
             1,1,1,1,1,0,0,0,
             1,1,1,1,1,1,0,0,
             1,1,1,1,1,1,1,0};

    float fsum=0;
    int i=0,lower_nat8=(nat/8)*8;
    //arrange r0 in duplicates for easy loading into hi and low word
    __m128 *mr= (__m128*)r;
    __m256 sum= _mm256_setzero_ps();
    {
        __m256 r0 = _mm256_castps128_ps256(mr[0]);
        __m256 r1 = _mm256_castps128_ps256(mr[1]);
        __m256 r2 = _mm256_castps128_ps256(mr[2]);
        r0=_mm256_insertf128_ps(r0,_mm256_castps256_ps128(r0),1);
        r1=_mm256_insertf128_ps(r1,_mm256_castps256_ps128(r1),1);
        r2=_mm256_insertf128_ps(r2,_mm256_castps256_ps128(r2),1);

        //4th multiplication unecessary as it is all zeros
        __m256 d0 = _mm256_broadcast_ss(&invd0d0);
        __m256 one=_mm256_set1_ps(1.0f);
        //8 points at a time then mask out
        for( ;i <lower_nat8; i+=8){
            __m256 mx1 = _mm256_load_ps(&(x1[i]));
            __m256 my1 = _mm256_load_ps(&(y1[i]));
            __m256 mz1 = _mm256_load_ps(&(z1[i]));
            __m256 mx2 = _mm256_load_ps(&(x2[i]));

            __m256 tx1= _mm256_add_ps(_mm256_mul_ps(my1,_mm256_shuffle_ps(r0,r0,0x55)),_mm256_mul_ps(mx1,_mm256_shuffle_ps(r0,r0,0x00)));
            __m256 tx2= _mm256_add_ps(_mm256_mul_ps(mz1,_mm256_shuffle_ps(r0,r0,0xAA)),_mm256_shuffle_ps(r0,r0,0xFF));

            tx2 = _mm256_add_ps(tx2,tx1);
            tx2 = _mm256_sub_ps(tx2,mx2);
            mx2 = _mm256_load_ps(&(y2[i]));
            __m256 d1 = _mm256_mul_ps(tx2,tx2);

            tx1= _mm256_add_ps(_mm256_mul_ps(my1,_mm256_shuffle_ps(r1,r1,0x55)),_mm256_mul_ps(mx1,_mm256_shuffle_ps(r1,r1,0x00)));
            tx2= _mm256_add_ps(_mm256_mul_ps(mz1,_mm256_shuffle_ps(r1,r1,0xAA)),_mm256_shuffle_ps(r1,r1,0xFF));
            tx2= _mm256_add_ps(tx2,tx1);
            tx2= _mm256_sub_ps(tx2,mx2);
            mx2= _mm256_load_ps(&(z2[i]));
            d1 = _mm256_add_ps(d1,_mm256_mul_ps(tx2,tx2));

            tx1= _mm256_add_ps(_mm256_mul_ps(my1,_mm256_shuffle_ps(r2,r2,0x55)),_mm256_mul_ps(mx1,_mm256_shuffle_ps(r2,r2,0x00)));
            tx2= _mm256_add_ps(_mm256_mul_ps(mz1,_mm256_shuffle_ps(r2,r2,0xAA)),_mm256_shuffle_ps(r2,r2,0xFF));
            tx2= _mm256_add_ps(tx2,tx1);
            tx2= _mm256_sub_ps(tx2,mx2);
            d1 = _mm256_add_ps(d1,_mm256_mul_ps(tx2,tx2));
            _mm256_store_ps(&(d[i]),d1); //write out 8 differences
            mx1= _mm256_mul_ps(d1,d0);
            mx1= _mm256_add_ps(mx1,one);
#ifdef FAST_DIVISION
            mx1= _mm256_rcp_ps(mx1);
#else
            mx1= _mm256_div_ps(one,mx1);
#endif
            sum= _mm256_add_ps(sum,mx1);

        }
        for( ;i<nat; i+=8){
            __m256 mask=_mm256_load_ps(&(mask_array[(nat-lower_nat8-1)*8]));
            __m256 mx1 = _mm256_load_ps(&(x1[i]));
            __m256 my1 = _mm256_load_ps(&(y1[i]));
            __m256 mz1 = _mm256_load_ps(&(z1[i]));
            __m256 mx2 = _mm256_load_ps(&(x2[i]));

            __m256 tx1= _mm256_add_ps(_mm256_mul_ps(my1,_mm256_shuffle_ps(r0,r0,0x55)),_mm256_mul_ps(mx1,_mm256_shuffle_ps(r0,r0,0x00)));
            __m256 tx2= _mm256_add_ps(_mm256_mul_ps(mz1,_mm256_shuffle_ps(r0,r0,0xAA)),_mm256_shuffle_ps(r0,r0,0xFF));

            tx2 = _mm256_add_ps(tx2,tx1);
            tx2 = _mm256_sub_ps(tx2,mx2);
            mx2 = _mm256_load_ps(&(y2[i]));
            __m256 d1 = _mm256_mul_ps(tx2,tx2);

            tx1= _mm256_add_ps(_mm256_mul_ps(my1,_mm256_shuffle_ps(r1,r1,0x55)),_mm256_mul_ps(mx1,_mm256_shuffle_ps(r1,r1,0x00)));
            tx2= _mm256_add_ps(_mm256_mul_ps(mz1,_mm256_shuffle_ps(r1,r1,0xAA)),_mm256_shuffle_ps(r1,r1,0xFF));
            tx2= _mm256_add_ps(tx2,tx1);
            tx2= _mm256_sub_ps(tx2,mx2);
            mx2= _mm256_load_ps(&(z2[i]));
            d1 = _mm256_add_ps(d1,_mm256_mul_ps(tx2,tx2));

            tx1= _mm256_add_ps(_mm256_mul_ps(my1,_mm256_shuffle_ps(r2,r2,0x55)),_mm256_mul_ps(mx1,_mm256_shuffle_ps(r2,r2,0x00)));
            tx2= _mm256_add_ps(_mm256_mul_ps(mz1,_mm256_shuffle_ps(r2,r2,0xAA)),_mm256_shuffle_ps(r2,r2,0xFF));
            tx2= _mm256_add_ps(tx2,tx1);
            tx2= _mm256_sub_ps(tx2,mx2);
            d1 = _mm256_add_ps(d1,_mm256_mul_ps(tx2,tx2));

            _mm256_store_ps(&(d[i]),d1); //write out 8 differences
            mx1= _mm256_mul_ps(d1,d0);
            mx1= _mm256_add_ps(mx1,one);

#ifdef FAST_DIVISION
            mx1= _mm256_rcp_ps(mx1);
#else
            mx1= _mm256_div_ps(one,mx1);
#endif
            sum= _mm256_add_ps(sum, _mm256_mul_ps(mx1,mask)); //multiply by mask
        }
    }//end summation

    sum=_mm256_add_ps(sum,_mm256_permute2f128_ps(sum,sum,0x01));
    sum=_mm256_hadd_ps(sum,sum);
    sum=_mm256_hadd_ps(sum,sum);
    _mm_store_ss(&(fsum),_mm256_castps256_ps128(sum));
    return(fsum);
}

void R34v4_sse3 (float r[12],float *x,float *Rx){
    __m128 r0 = _mm_load_ps(r);
    __m128 r1 = _mm_load_ps(&(r[4]));
    __m128 r2 = _mm_load_ps(&(r[8]));
    __m128 one= _mm_set_ss(1.0f);
    __m128 mu = _mm_load_ps(x);
    __m128 m0 = _mm_mul_ps(r0,mu);
    __m128 m1 = _mm_mul_ps(r1,mu);
    __m128 m2 = _mm_mul_ps(r2,mu);
    __m128 s1 = _mm_hadd_ps(m0,m1);
    __m128 s2 = _mm_hadd_ps(m2,one);
    _mm_store_ps(Rx,_mm_hadd_ps(s1,s2));
}

float rmsd_uncentered_avx (int nat,float *c1x,float *c1y, float *c1z,float *c2x,float *c2y, float *c2z, float *rm){

    double u[3][3];
    double invdnat=1.0/(double)nat;
    float fnat=(float)nat;
    float invfnat=1.0f/fnat;
    float c0[8] __attribute__ ((aligned (32)));
    float c1[16] __attribute__ ((aligned (32)));
    float frr[8] __attribute__ ((aligned (32)));
    int upper8 = (nat%8) ? (nat/8)*8+8 : nat;
    {
        //use half registers to do everything at once rather than trying for different products
        //essentially SSE3 recipe but with twice as many registers to save one load pass to center coords
        //the 4 hadds will split the lower and upper halves into a variable with 8 sums
        // A E
        // B F
        // C G
        // D H -> sumA sumB sumC sumD sumE sumF sumG sumH
        //these two vectors/registers will also hold final products after HADDS
        __m256 sxxssq = _mm256_setzero_ps();  // hold cross sums
        __m256 s1xs2x = _mm256_setzero_ps();  // hold sums
        {
            __m256 s1ys2y = _mm256_setzero_ps();
            __m256 s1zs2z = _mm256_setzero_ps();
            __m256 sxysyz = _mm256_setzero_ps();
            __m256 sxzszx = _mm256_setzero_ps();
            __m256 syxszy = _mm256_setzero_ps();
            __m256 syyszz = _mm256_setzero_ps();

            for(int i=0 ; i < upper8; i += 8){
                //load the 4 sets of coords from molecule 2 and then load x,y,z of molecule 1
                __m256 mmc2x = _mm256_loadu_ps(&(c2x[i]));
                __m256 mmc2y = _mm256_loadu_ps(&(c2y[i]));
                __m256 mmc2z = _mm256_loadu_ps(&(c2z[i]));
                __m256 t1    = _mm256_loadu_ps(&(c1x[i])); //do everything necessary with c1x and then use as temp register
                __m256 t2    = _mm256_mul_ps(t1,t1);

                //block with extra temp
                __m256 mmc1y= _mm256_mul_ps(t1,mmc2x); //used as a temp
                s1xs2x      = _mm256_add_ps(s1xs2x,_mm256_permute2f128_ps(t1,mmc2x,0x20));
                s1xs2x      = _mm256_add_ps(s1xs2x,_mm256_permute2f128_ps(t1,mmc2x,0x31));
                sxxssq      = _mm256_add_ps(sxxssq,_mm256_permute2f128_ps(mmc1y,t2,0x20));
                sxxssq      = _mm256_add_ps(sxxssq,_mm256_permute2f128_ps(mmc1y,t2,0x31));
                __m256 mmc1z= _mm256_mul_ps(t1,mmc2y);
                mmc1y       = _mm256_loadu_ps(&(c1y[i]));
                t2          = _mm256_mul_ps(mmc1y,mmc2z);
                sxysyz      = _mm256_add_ps(sxysyz,_mm256_permute2f128_ps(mmc1z,t2,0x20));
                sxysyz      = _mm256_add_ps(sxysyz,_mm256_permute2f128_ps(mmc1z,t2,0x31));
                mmc1z = _mm256_loadu_ps(&(c1z[i]));
                t1          = _mm256_mul_ps(t1,mmc2z); //t1 can be used freely now
                t2          = _mm256_mul_ps(mmc1z,mmc2x);
                sxzszx      = _mm256_add_ps(sxzszx,_mm256_permute2f128_ps(t1,t2,0x20));
                sxzszx     = _mm256_add_ps(sxzszx,_mm256_permute2f128_ps(t1,t2,0x31));

                //finish calculating ssq
                t1=_mm256_mul_ps(mmc2x,mmc2x);
                t2=_mm256_mul_ps(mmc2y,mmc2y);
                t1=_mm256_add_ps(t1,_mm256_mul_ps(mmc2z,mmc2z));
                t2=_mm256_add_ps(t2,_mm256_mul_ps(mmc1y,mmc1y));
                t1=_mm256_add_ps(t1,t2);
                t1=_mm256_add_ps(t1,_mm256_mul_ps(mmc1z,mmc1z));
                sxxssq    = _mm256_add_ps(sxxssq,_mm256_permute2f128_ps(t1,t1,0x08));//pass lower 128 and zero lower  0000 1000
                sxxssq    = _mm256_add_ps(sxxssq,_mm256_permute2f128_ps(t1,t1,0x18));//pass upper 128 and zero lower; 0001 1000

                //do other sums
                s1ys2y     =_mm256_add_ps(s1ys2y,_mm256_permute2f128_ps(mmc1y,mmc2y,0x20));
                s1ys2y     =_mm256_add_ps(s1ys2y,_mm256_permute2f128_ps(mmc1y,mmc2y,0x31));
                s1zs2z     =_mm256_add_ps(s1zs2z,_mm256_permute2f128_ps(mmc1z,mmc2z,0x20));
                s1zs2z     =_mm256_add_ps(s1zs2z,_mm256_permute2f128_ps(mmc1z,mmc2z,0x31));

                //do other cross_sums
                //registers start freeing up
                t1 =  _mm256_mul_ps(mmc1y,mmc2x);
                mmc1y=_mm256_mul_ps(mmc1y,mmc2y);
                t2 =  _mm256_mul_ps(mmc1z,mmc2y);
                mmc1z=_mm256_mul_ps(mmc1z,mmc2z);
                syxszy     = _mm256_add_ps(syxszy,_mm256_permute2f128_ps(t1,t2,0x20));
                syxszy     = _mm256_add_ps(syxszy,_mm256_permute2f128_ps(t1,t2,0x31));
                syyszz     = _mm256_add_ps(syyszz,_mm256_permute2f128_ps(mmc1y,mmc1z,0x20));
                syyszz     = _mm256_add_ps(syyszz,_mm256_permute2f128_ps(mmc1y,mmc1z,0x31));
            }//end loop
            //now do two groups of 4 hadds
            //      a1 a2 a3 a4 b1 b2 b3 b4
            //      c1 c2 c3 c4 d1 d2 d3 d4
            //      a12 a34 c12 c34 b12 b34 d12 d34
            //
            //      e1 e2 e3 e4 f1 f2 f3 f4
            //      g1 g2 g3 g4 h1 h2 h3 h4
            //      e12 e34 g12 g34 f12 f34 h12 h34
            //
            //      a12 a34 c12 c34 b12 b34 d12 d34
            //      e12 e34 g12 g34 f12 f34 h12 h34
            //      A C E G B D F H - even/odd
            s1xs2x =_mm256_hadd_ps(_mm256_hadd_ps(sxxssq,s1xs2x),_mm256_hadd_ps(s1ys2y,s1zs2z));
            sxxssq =_mm256_hadd_ps(_mm256_hadd_ps(sxysyz,sxzszx),_mm256_hadd_ps(syxszy,syyszz));
//   print_m256(s1xs2x);
//   print_m256(sxxssq);
        }//no need for other cross sums anymore

        __m256 normed_sums=_mm256_mul_ps(s1xs2x,_mm256_broadcast_ss(&invfnat));
        _mm256_store_ps(c0,normed_sums);

        //sums  sxx s1x s1y s1z ssq s2x s2y s2z
        //csums sxy sxz syx syy syz szx szy szz
        //      s1x s1x s1y s1y s2z s2x s2y s2z
        //      1   1   2   2   3   1   2   3
        //      01  01  10  10  11  01  10  11  ->
        //      1110 0111 1010 0101 (leftmost number are least significant 2 bits)
        //       E     7    A    5 // A5 E7 in 2 separate permutes to generate the halves and a permute2f to blend them or can load separate mask
        // want low of first word and high of second word   00 11  00 00 = 0x30
        //

        __m256 t1=_mm256_permute2f128_ps(_mm256_permute_ps(normed_sums,0xA5),_mm256_permute_ps(normed_sums,0xE7),0x30);
        __m256 t2=_mm256_permute2f128_ps(s1xs2x,s1xs2x,1);
        // t1   (s1x s1x s1y s1y s2z s2x s2y s2z)*invfnat
        // t2   (ssq s2x s2y s2z sxx s1x s1y s1z)

        //do sxx first
        //second element of t1*s1xs2x has desired element
        //move to first element - rest is junk
        s1xs2x=_mm256_sub_ps(s1xs2x,_mm256_permute_ps(_mm256_mul_ps(t1,t2),0x01));
        s1xs2x=_mm256_permute_ps(s1xs2x,0);
        s1xs2x=_mm256_permute2f128_ps(s1xs2x,s1xs2x,0);
        _mm256_store_ps(c1,s1xs2x);
        //       s2y s2z s2x s2y s1y s1z s1z s1z  <- want this for csums subtraction
        //       2   3   1   2   2   3   3   3
        //      10   11  01  10  10  11  11  11
        //      1111 1110 1001 1110
        //        F   E    9    E

        t2=_mm256_permute2f128_ps(_mm256_permute_ps(t2,0x9E),_mm256_permute_ps(t2,0xFE),0x30);
        sxxssq=_mm256_sub_ps(sxxssq,_mm256_mul_ps(t1,t2));
        _mm256_store_ps(&(c1[8]),sxxssq);
        __m256 zero = _mm256_setzero_ps();
        //put into 3 vectors for calculation of cross_products and determinant
        //sxy sxz syx syy syz szx szy szz
        //        sxy sxz
        //    syx syy
        //syz ...........................
        //
        //
        //sxxssq    sxy sxz syx syy syz szx szy szz
        //t1        sxy sxy sxy sxz syz syz syz szz  permute mask 0001 ->  01 00 00 00 -> 0x40
        //s1xs2x    sxx sxx sxx sxx sxx sxx sxx sxx  blend mask 00000010
        //r0r0      sxy sxx sxy sxz syz syz syz szz
        //r0r0       0  sxx sxy sxz  0   0   0   0
        //r0r0       0  sxx sxy sxz  0  sxx sxy sxz
        //sxxssq    sxy sxz syx syy syz szx szy szz
        //t2        sxy syx syy sxy syz szy szz syz  permute mask 0230   0011 1000 -> 0x38
        //          syz szx szy szz sxy sxz syx syy  swapped sxxssq 0000 0001
        //t1        syz syz syz syz sxy sxy sxy sxy  permute
        //t2        sxy syx syy syz syz szy szz syz  blend mask 00001000 -> 0x08
        //r1r2      sxy syx syy syz syz szx szy szz  low word t2 high word sxxssq 0010 0000
        //r1r2      0   syx syy syz  0  szx szy szz
        t1=_mm256_permute_ps(sxxssq,0x40);
        __m256 r0r0=_mm256_blend_ps(t1,s1xs2x,0x02);
        r0r0 =_mm256_blend_ps(r0r0,zero,0xF1); //mask 11110001
        r0r0 =_mm256_permute2f128_ps(r0r0,r0r0,0);
        t2=_mm256_permute_ps(sxxssq,0x38);
        t1=_mm256_permute_ps(_mm256_permute2f128_ps(sxxssq,sxxssq,0x01),0);
        t2=_mm256_blend_ps(t2,t1,0x08);
        __m256 r1r2 =_mm256_permute2f128_ps(t2,sxxssq,0x30);//00110000
        r1r2=_mm256_blend_ps(r1r2,zero,0x11); //mask 00010001
        //do determinant cross
        //r1r2    0 syx syy syz  0  szx szy szz
        //t1      macromask 1 3 2 0    01 11 10 00  78
        //t2      macromask 2 1 3 0    10 01 11 00  9C
        //        swap
        //  multiply - lower - upper is cross prodduct
        //  multiply by r0r0  subtract two halves and do hadds
        t1=_mm256_permute_ps(r1r2,0x78);
        t2=_mm256_permute_ps(r1r2,0x9C);
        t1=_mm256_mul_ps(t1,_mm256_permute2f128_ps(t2,t2,0x01));
        t1=_mm256_sub_ps(t1,_mm256_permute2f128_ps(t1,t1,0x01));
        __m256 det=_mm256_mul_ps(r0r0,t1);//det sum is in lower half -sum with other dot products
        __m256 r0r1r0r2=_mm256_mul_ps(r0r0,r1r2);
        __m256 r1r1r2r2=_mm256_mul_ps(r1r2,r1r2);
        r0r0=_mm256_mul_ps(r0r0,r0r0);
        r1r2=_mm256_mul_ps(r1r2,_mm256_permute2f128_ps(r1r2,r1r2,0x01));
        r0r0=_mm256_permute2f128_ps(r0r0,r1r2,0x30);//lower half r0r0 upper half r1r2
        //hadds
        //r0r0 r1r2
        //r0r1 r0r2
        //r1r1 r2r2
        // det -det  -> r0r0 r0r1 r1r1 det r1r2 r0r2 r2r2 -det
        _mm256_store_ps(frr,_mm256_hadd_ps(_mm256_hadd_ps(r0r0,r0r1r0r2),_mm256_hadd_ps(r1r1r2r2,det)));


// note that the _MM_SHUFFLE macro reverses the order so that lower integer is first
//                                                   (lower word t1 upperword of t2  -upperword t1 lower word t2)
//  __m256 cross = _mm_sub_ps (_mm_mul_ps(_mm_shuffle_ps(mxy, mxy, _MM_SHUFFLE(1,3,2,0)), _mm_shuffle_ps(mxz, mxz, _MM_SHUFFLE(2,1,3,0))), -
//                             _mm_mul_ps(_mm_shuffle_ps(mxy, mxy, _MM_SHUFFLE(2,1,3,0)), _mm_shuffle_ps(mxz, mxz, _MM_SHUFFLE(1,3,2,0))));//cross_product

    } //end AVX section
    const float *r = &(c1[7]);
    //c0=(ssq s1x s1y s1z sxx s2x s2y s2z)*invfnat
    const float *center1 = &(c0[1]);
    const float *center2 = &(c0[5]);
    double ssq=(double)((c0[4]-c0[1]*c0[1]-c0[2]*c0[2]-c0[3]*c0[3]-c0[5]*c0[5]-c0[6]*c0[6]-c0[7]*c0[7])*fnat);
    //normalize ssq
    //double det= (double)( r[0] * ( (r[4]*r[8]) - (r[5]*r[7]) )- r[1] * ( (r[3]*r[8]) - (r[5]*r[6]) ) + r[2] * ( (r[3]*r[7]) - (r[4]*r[6])));
    double det= (double)frr[3];
    double detsq=det*det;
    double rr[6]={(double) frr[0],(double) frr[1],(double) frr[2],(double) frr[5],(double) frr[4],(double) frr[6]};

    //lower triangular matrix rr
    double inv3=1.0/3.0;
    double spur=(rr[0]+rr[2]+rr[5])*inv3;
    double cof= (rr[2]*rr[5] - rr[4]*rr[4] + rr[0]*rr[5]- rr[3]*rr[3] + rr[0]*rr[2] - rr[1]*rr[1]) *inv3;
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

    double rms=(ssq-d-d)*invdnat;
    rms=(rms>1e-8)?sqrt(rms) : 0.0f;
    double mr[3][3]={{r[0],r[1],r[2]},
                    {r[3],r[4],r[5]},
                    {r[6],r[7],r[8]}};
    rmatrix<double>(d,mr,u);
    float mrr[16] __attribute__ ((aligned (16)))=
            {static_cast<float>(-u[0][0]),static_cast<float>(-u[1][0]),static_cast<float>(-u[2][0]),center2[0],
             static_cast<float>(-u[0][1]),static_cast<float>(-u[1][1]),static_cast<float>(-u[2][1]),center2[1],
             static_cast<float>(-u[0][2]),static_cast<float>(-u[1][2]),static_cast<float>(-u[2][2]),center2[2]};
    float w[4]__attribute__ ((aligned (16)));
    float v[4]__attribute__ ((aligned (16)))={center1[0],center1[1],center1[2],1.0f};
    R34v4_sse3(mrr,v,w);
    rm[0] =u[0][0]; rm[1]=u[1][0];rm[2] =u[2][0];rm[3]=w[0];
    rm[4] =u[0][1]; rm[5]=u[1][1];rm[6] =u[2][1];rm[7]=w[1];
    rm[8] =u[0][2]; rm[9]=u[1][2];rm[10]=u[2][2];rm[11]=w[2];
// rm[12]=0.0f;rm[13]=0.0f;rm[14]=0.0f;rm[15]=1.0f;
    return((float) rms);
}

#define COORDS_BUFFER_SIZE 768

float kabsch_quat_soa_avx(int nat, float *x1, float *y1,float *z1, float *x2,float *y2, float *z2,float *r,float *mem){
    //most of time will be small number of coords
    int upper_nat8= (nat%8)? (nat/8)*8+8: nat;
    int size= 6*upper_nat8*sizeof(float);

    memset(mem,0,size);
    float* c1x= mem;
    float* c1y= &(mem[upper_nat8]);
    float* c1z= &(mem[upper_nat8*2]);
    float* c2x= &(mem[upper_nat8*3]);
    float* c2y= &(mem[upper_nat8*4]);
    float* c2z= &(mem[upper_nat8*5]);
    memcpy(c1x, x1, nat * sizeof(float));
    memcpy(c1y, y1, nat * sizeof(float));
    memcpy(c1z, z1, nat * sizeof(float));
    memcpy(c2x, x2, nat * sizeof(float));
    memcpy(c2y, y2, nat * sizeof(float));
    memcpy(c2z, z2, nat * sizeof(float));

    float rms=rmsd_uncentered_avx(nat,c1x,c1y,c1z,c2x,c2y,c2z,r);
    return(rms);
}


/**************************************************************************
Implemetation of Kabsch algoritm for finding the best rotation matrix
---------------------------------------------------------------------------
x    - x(i,m) are coordinates of atom m in set x            (input)
y    - y(i,m) are coordinates of atom m in set y            (input)
n    - n is number of atom pairs                            (input)
mode  - 0:calculate rms only                                (input)
1:calculate u,t only                                (takes medium)
2:calculate rms,u,t                                 (takes longer)
rms   - sum of w*(ux+t-y)**2 over all atom pairs            (output)
u    - u(i,j) is   rotation  matrix for best superposition  (output)
t    - t(i)   is translation vector for best superposition  (output)
**************************************************************************/
bool Kabsch(Coordinates & x,
            Coordinates & y,
            int n,
            int mode,
            float *rms,
            float t[3],
            float u[3][3])
{
    int i, j, m, m1, l, k;
    double e0, rms1, d, h, g;
    double cth, sth, sqrth, p, det, sigma;
    double xc[3], yc[3];
    double a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
    double sqrt3 = 1.73205080756888, tol = 0.01;
    int ip[] = { 0, 1, 3, 1, 2, 4, 3, 4, 5 };
    int ip2312[] = { 1, 2, 0, 1 };

    int a_failed = 0, b_failed = 0;
    double epsilon = 0.00000001;

    //initializtation
    *rms = 0;
    rms1 = 0;
    e0 = 0;
    double s1[3], s2[3];
    double sx[3], sy[3], sz[3];
    for (i = 0; i < 3; i++)
    {
        s1[i] = 0.0;
        s2[i] = 0.0;

        sx[i] = 0.0;
        sy[i] = 0.0;
        sz[i] = 0.0;
    }

    for (i = 0; i<3; i++)
    {
        xc[i] = 0.0;
        yc[i] = 0.0;
        t[i] = 0.0;
        for (j = 0; j<3; j++)
        {
            u[i][j] = 0.0;
            r[i][j] = 0.0;
            a[i][j] = 0.0;
            if (i == j)
            {
                u[i][j] = 1.0;
                a[i][j] = 1.0;
            }
        }
    }

    if (n<1)
    {
        return false;
    }
    float c1tmp[3], c2tmp[3];
    float s1tmp[3], s2tmp[3];
    float sxtmp[3], sytmp[3], sztmp[3];
    for (i = 0; i < 3; i++)
    {
        s1tmp[i] = 0.0;
        s2tmp[i] = 0.0;
        sxtmp[i] = 0.0;
        sytmp[i] = 0.0;
        sztmp[i] = 0.0;
    }
    //compute centers for vector sets x, y
    for (i = 0; i<n; i++)
    {

        c1tmp[0] = x.x[i];
        c1tmp[1] = x.y[i];
        c1tmp[2] = x.z[i];

        c2tmp[0] = y.x[i];
        c2tmp[1] = y.y[i];
        c2tmp[2] = y.z[i];

        s1tmp[0] += c1tmp[0];
        s1tmp[1] += c1tmp[1];
        s1tmp[2] += c1tmp[2];

        s2tmp[0] += c2tmp[0];
        s2tmp[1] += c2tmp[1];
        s2tmp[2] += c2tmp[2];


        for (j = 0; j < 3; j++)
        {
            sxtmp[j] += c1tmp[0] * c2tmp[j];
            sytmp[j] += c1tmp[1] * c2tmp[j];
            sztmp[j] += c1tmp[2] * c2tmp[j];
        }
    }
    for (i = 0; i < 3; i++)
    {
        s1[i]=s1tmp[i];
        s2[i]=s2tmp[i];
        sx[i]=sxtmp[i];
        sy[i]=sytmp[i];
        sz[i]=sztmp[i];
    }

    for (i = 0; i < 3; i++)
    {
        xc[i] = s1[i] / n;
        yc[i] = s2[i] / n;
    }
    if (mode == 2 || mode == 0)
    {
        for (int mm = 0; mm < n; mm++)
        {
            e0 += (x.x[mm] - xc[0]) * (x.x[mm] - xc[0]) +
                  (y.x[mm] - yc[0]) * (y.x[mm] - yc[0]);
            e0 += (x.y[mm] - xc[1]) * (x.y[mm] - xc[1]) +
                  (y.y[mm] - yc[1]) * (y.y[mm] - yc[1]);
            e0 += (x.z[mm] - xc[2]) * (x.z[mm] - xc[2]) +
                  (y.z[mm] - yc[2]) * (y.z[mm] - yc[2]);
        }
    }

    for (j = 0; j < 3; j++)
    {
        r[j][0] = sx[j] - s1[0] * s2[j] / n;
        r[j][1] = sy[j] - s1[1] * s2[j] / n;
        r[j][2] = sz[j] - s1[2] * s2[j] / n;
    }


    //compute determinat of matrix r
    det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])\
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])\
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    sigma = det;

    //compute tras(r)*r
    m = 0;
    for (j = 0; j<3; j++)
    {
        for (i = 0; i <= j; i++)
        {
            rr[m] = r[0][i] * r[0][j] + r[1][i] * r[1][j] + r[2][i] * r[2][j];
            m++;
        }
    }

    double spur = (rr[0] + rr[2] + rr[5]) / 3.0;
    double cof = (((((rr[2] * rr[5] - rr[4] * rr[4]) + rr[0] * rr[5])\
        - rr[3] * rr[3]) + rr[0] * rr[2]) - rr[1] * rr[1]) / 3.0;
    det = det*det;

    for (i = 0; i<3; i++)
    {
        e[i] = spur;
    }

    if (spur>0)
    {
        d = spur*spur;
        h = d - cof;
        g = (spur*cof - det) / 2.0 - spur*h;

        if (h>0)
        {
            sqrth = sqrt(h);
            d = h*h*h - g*g;
            if (d<0.0) d = 0.0;
            d = atan2(sqrt(d), -g) / 3.0;
            cth = sqrth * cos(d);
            sth = sqrth*sqrt3*sin(d);
            e[0] = (spur + cth) + cth;
            e[1] = (spur - cth) + sth;
            e[2] = (spur - cth) - sth;

            if (mode != 0)
            {//compute a                
                for (l = 0; l<3; l = l + 2)
                {
                    d = e[l];
                    ss[0] = (d - rr[2]) * (d - rr[5]) - rr[4] * rr[4];
                    ss[1] = (d - rr[5]) * rr[1] + rr[3] * rr[4];
                    ss[2] = (d - rr[0]) * (d - rr[5]) - rr[3] * rr[3];
                    ss[3] = (d - rr[2]) * rr[3] + rr[1] * rr[4];
                    ss[4] = (d - rr[0]) * rr[4] + rr[1] * rr[3];
                    ss[5] = (d - rr[0]) * (d - rr[2]) - rr[1] * rr[1];

                    if (fabsf(ss[0]) <= epsilon) ss[0] = 0.0;
                    if (fabsf(ss[1]) <= epsilon) ss[1] = 0.0;
                    if (fabsf(ss[2]) <= epsilon) ss[2] = 0.0;
                    if (fabsf(ss[3]) <= epsilon) ss[3] = 0.0;
                    if (fabsf(ss[4]) <= epsilon) ss[4] = 0.0;
                    if (fabsf(ss[5]) <= epsilon) ss[5] = 0.0;

                    if (fabsf(ss[0]) >= fabs(ss[2]))
                    {
                        j = 0;
                        if (fabsf(ss[0]) < fabs(ss[5]))
                        {
                            j = 2;
                        }
                    }
                    else if (fabsf(ss[2]) >= fabsf(ss[5]))
                    {
                        j = 1;
                    }
                    else
                    {
                        j = 2;
                    }

                    d = 0.0;
                    j = 3 * j;
                    for (i = 0; i<3; i++)
                    {
                        k = ip[i + j];
                        a[i][l] = ss[k];
                        d = d + ss[k] * ss[k];
                    }


                    //if( d > 0.0 ) d = 1.0 / sqrt(d);
                    if (d > epsilon) d = 1.0 / sqrt(d);
                    else d = 0.0;
                    for (i = 0; i<3; i++)
                    {
                        a[i][l] = a[i][l] * d;
                    }
                }//for l

                d = a[0][0] * a[0][2] + a[1][0] * a[1][2] + a[2][0] * a[2][2];
                if ((e[0] - e[1]) >(e[1] - e[2]))
                {
                    m1 = 2;
                    m = 0;
                }
                else
                {
                    m1 = 0;
                    m = 2;
                }
                p = 0;
                for (i = 0; i<3; i++)
                {
                    a[i][m1] = a[i][m1] - d*a[i][m];
                    p = p + a[i][m1] * a[i][m1];
                }
                if (p <= tol)
                {
                    p = 1.0;
                    for (i = 0; i<3; i++)
                    {
                        if (p < fabs(a[i][m]))
                        {
                            continue;
                        }
                        p = fabs(a[i][m]);
                        j = i;
                    }
                    k = ip2312[j];
                    l = ip2312[j + 1];
                    p = sqrt(a[k][m] * a[k][m] + a[l][m] * a[l][m]);
                    if (p > tol)
                    {
                        a[j][m1] = 0.0;
                        a[k][m1] = -a[l][m] / p;
                        a[l][m1] = a[k][m] / p;
                    }
                    else
                    {//goto 40
                        a_failed = 1;
                    }
                }//if p<=tol
                else
                {
                    p = 1.0 / sqrt(p);
                    for (i = 0; i<3; i++)
                    {
                        a[i][m1] = a[i][m1] * p;
                    }
                }//else p<=tol  
                if (a_failed != 1)
                {
                    a[0][1] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
                    a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];
                    a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
                }
            }//if(mode!=0)       
        }//h>0

        //compute b anyway
        if (mode != 0 && a_failed != 1)//a is computed correctly
        {
            //compute b
            for (l = 0; l<2; l++)
            {
                d = 0.0;
                for (i = 0; i<3; i++)
                {
                    b[i][l] = r[i][0] * a[0][l] + r[i][1] * a[1][l] + r[i][2] * a[2][l];
                    d = d + b[i][l] * b[i][l];
                }
                //if( d > 0 ) d = 1.0 / sqrt(d);
                if (d > epsilon) d = 1.0 / sqrt(d);
                else d = 0.0;
                for (i = 0; i<3; i++)
                {
                    b[i][l] = b[i][l] * d;
                }
            }
            d = b[0][0] * b[0][1] + b[1][0] * b[1][1] + b[2][0] * b[2][1];
            p = 0.0;

            for (i = 0; i<3; i++)
            {
                b[i][1] = b[i][1] - d*b[i][0];
                p += b[i][1] * b[i][1];
            }

            if (p <= tol)
            {
                p = 1.0;
                for (i = 0; i<3; i++)
                {
                    if (p<fabs(b[i][0]))
                    {
                        continue;
                    }
                    p = fabs(b[i][0]);
                    j = i;
                }
                k = ip2312[j];
                l = ip2312[j + 1];
                p = sqrt(b[k][0] * b[k][0] + b[l][0] * b[l][0]);
                if (p > tol)
                {
                    b[j][1] = 0.0;
                    b[k][1] = -b[l][0] / p;
                    b[l][1] = b[k][0] / p;
                }
                else
                {
                    //goto 40
                    b_failed = 1;
                }
            }//if( p <= tol )
            else
            {
                p = 1.0 / sqrt(p);
                for (i = 0; i<3; i++)
                {
                    b[i][1] = b[i][1] * p;
                }
            }
            if (b_failed != 1)
            {
                b[0][2] = b[1][0] * b[2][1] - b[1][1] * b[2][0];
                b[1][2] = b[2][0] * b[0][1] - b[2][1] * b[0][0];
                b[2][2] = b[0][0] * b[1][1] - b[0][1] * b[1][0];
                //compute u
                for (i = 0; i<3; i++)
                {
                    for (j = 0; j<3; j++)
                    {
                        u[i][j] = b[i][0] * a[j][0] + b[i][1] * a[j][1]\
                            + b[i][2] * a[j][2];
                    }
                }
            }

            //compute t
            for (i = 0; i<3; i++)
            {
                t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1])\
                    - u[i][2] * xc[2];
            }
        }//if(mode!=0 && a_failed!=1)
    }//spur>0
    else //just compute t and errors
    {
        //compute t
        for (i = 0; i<3; i++)
        {
            t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) - u[i][2] * xc[2];
        }
    }//else spur>0 

    //compute rms
    for (i = 0; i<3; i++)
    {
        if (e[i] < 0) e[i] = 0;
        e[i] = sqrt(e[i]);
    }
    d = e[2];
    if (sigma < 0.0)
    {
        d = -d;
    }
    d = (d + e[1]) + e[0];

    if (mode == 2 || mode == 0)
    {
        rms1 = (e0 - d) - d;
        if (rms1 < 0.0)
            rms1 = 0.0;
    }

    *rms = rms1;

    return true;

}

