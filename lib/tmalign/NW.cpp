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
#include "NW.h"
#include "basic_fun.h"
/*    Please note this fucntion is not a correct implementation of 
*     the N-W dynamic programming because the score tracks back only 
*     one layer of the matrix. This code was exploited in TM-align 
*     because it is about 1.5 times faster than a complete N-W code
*     and does not influence much the final structure alignment result.
*/
void NWDP_TM( float **score, bool **path, float **val,
    int len1, int len2, float gap_open, int j2i[])
{
    //NW dynamic programming for alignment
    //not a standard implementation of NW algorithm
    //Input: score[1:len1, 1:len2], and gap_open
    //Output: j2i[1:len2] \in {1:len1} U {-1}
    //path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical

    int i, j;
    float h, v, d;

    //initialization
    val[0][0]=0;
    for(i=0; i<=len1; i++)
    {
        val[i][0]=0;
        path[i][0]=false; //not from diagonal
    }

    for(j=0; j<=len2; j++)
    {
        val[0][j]=0;
        path[0][j]=false; //not from diagonal
        j2i[j]=-1;    //all are not aligned, only use j2i[1:len2]
    }


    //decide matrix and path
    for(i=1; i<=len1; i++)
    {
        for(j=1; j<=len2; j++)
        {
            d=val[i-1][j-1]+score[i][j]; //diagonal

            //symbol insertion in horizontal (= a gap in vertical)
            h=val[i-1][j];
            h += (path[i-1][j]) ? gap_open : 0.0;

            //symbol insertion in vertical
            v=val[i][j-1];
            v += (path[i][j-1]) ? gap_open : 0.0f;


            path[i][j]= (d>=h && d>=v) ? true : false; //from diagonal
            val[i][j]=(v>=h) ? v : h;
            val[i][j]=(path[i][j] == true) ? d : val[i][j];
        } //for i
    } //for j

    //trace back to extract the alignment
    i=len1;
    j=len2;
    while(i>0 && j>0)
    {
        if(path[i][j]) //from diagonal
        {
            j2i[j-1]=i-1;
            i--;
            j--;
        }
        else
        {
            h=val[i-1][j];
            if(path[i-1][j]) h +=gap_open;

            v=val[i][j-1];
            if(path[i][j-1]) v +=gap_open;

            if(v>=h)
                j--;
            else
                i--;
        }
    }
}


void NWDP_TM( float **score, bool **path, float **val,
    const Coordinates &x, const Coordinates &y, int len1, int len2, float t[3], float u[3][3],
    float d02, float gap_open, int j2i[], float * tmp)
{
    //NW dynamic programming for alignment
    //not a standard implementation of NW algorithm
    //Input: vectors x, y, rotation matrix t, u, scale factor d02, and gap_open
    //Output: j2i[1:len2] \in {1:len1} U {-1}
    //path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical

    int i, j;
    float h, v, d;

    //initialization
    val[0][0]=0;
    for(i=0; i<=len1; i++)
    {
        val[i][0]=0;
        path[i][0]=false; //not from diagonal
    }

    for(j=0; j<=len2; j++)
    {
        val[0][j]=0;
        path[0][j]=false; //not from diagonal
        j2i[j]=-1;    //all are not aligned, only use j2i[1:len2]
    }
    float xx[3], dij;

    simd_float vd02 = simdf32_set(d02);
    simd_float one = simdf32_set(1.0f);
    float * distArray = tmp;
    //decide matrix and path
    for(i=1; i<=len1; i++)
    {
        BasicFunction::transform(t, u, x.x[i-1], x.y[i-1], x.z[i-1], xx[0], xx[1], xx[2]);
        simd_float vXx = simdf32_set(xx[0]);
        simd_float vXy = simdf32_set(xx[1]);
        simd_float vXz = simdf32_set(xx[2]);
        for(j=0; j<=len2; j+=VECSIZE_FLOAT){
            //    float d1=xx-yx;
            //    float d2=xy-yy;
            //    float d3=xz-yz;
            //    return (d1*d1 + d2*d2 + d3*d3);
            simd_float ya_x = simdf32_load(&y.x[j]);
            simd_float ya_y = simdf32_load(&y.y[j]);
            simd_float ya_z = simdf32_load(&y.z[j]);
            ya_x = simdf32_sub(vXx, ya_x);
            ya_y = simdf32_sub(vXy, ya_y);
            ya_z = simdf32_sub(vXz, ya_z);
            ya_x = simdf32_mul(ya_x, ya_x);
            ya_y = simdf32_mul(ya_y, ya_y);
            ya_z = simdf32_mul(ya_z, ya_z);
            simd_float res = simdf32_add(ya_x, ya_y);
            simd_float dij = simdf32_add(res, ya_z);
            simd_float oneDividedDist = simdf32_div(one, simdf32_add(one, simdf32_div(dij,vd02)));
            simdf32_store(&distArray[j], oneDividedDist);
        }

        for(j=1; j<=len2; j++)
        {
            //d=val[i-1][j-1]+score[i][j]; //diagonal
//            dij=dist(xx[0], xx[1], xx[2], y.x[(j-1)], y.y[(j-1)], y.z[(j-1)]);
//            dij = 1.0/(1+dij/d02);
//            d=val[i-1][j-1] + dij;
            d=val[i-1][j-1] + distArray[j-1];

            //symbol insertion in horizontal (= a gap in vertical)
            h=val[i-1][j];
            h += (path[i-1][j]) ? gap_open : 0.0;


            //symbol insertion in vertical
            v=val[i][j-1];
            v += (path[i][j-1]) ? gap_open : 0.0f;

            path[i][j]= (d>=h && d>=v) ? true : false; //from diagonal
            val[i][j]=(v>=h) ? v : h;
            val[i][j]=(path[i][j] == true) ? d : val[i][j];
            
//            const unsigned char mode1 = (curr_sM_G_D_vec[j].H == H) ? M : E;
//            const unsigned char mode2 = (curr_sM_G_D_vec[j].H == curr_sM_G_D_vec[j].F) ? F : E;
//            const unsigned char mode = std::max(mode1, mode2);
//       std::cout << (int) mode << " ";
//            btMatrix[pos/4] |= mode << (pos % 4) * 2;
            
        } //for i
    } //for j
    //trace back to extract the alignment
//    std::cout << val[i-1][j-1] << std::endl;
    i=len1;
    j=len2;
    while(i>0 && j>0)
    {
        if(path[i][j]) //from diagonal
        {
            j2i[j-1]=i-1;
            i--;
            j--;
        }
        else
        {
            h=val[i-1][j];
            if(path[i-1][j]) {
                h +=gap_open;
            }

            v=val[i][j-1];
            if(path[i][j-1]) {
                v +=gap_open;
            }

            if(v>=h) {
                j--;
            } else {
                i--;
            }
        }
    }
}

//+ss
void NWDP_TM(float **score, bool **path, float **val,
             const char *secx, const char *secy, const int len1, const int len2,
             const float gap_open, int j2i[])
{
    //NW dynamic programming for alignment
    //not a standard implementation of NW algorithm
    //Input: secondary structure secx, secy, and gap_open
    //Output: j2i[1:len2] \in {1:len1} U {-1}
    //path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical

    int i, j;
    float h, v, d;

    //initialization
    val[0][0]=0;
    for(i=0; i<=len1; i++)
    {
        val[i][0]=0;
        path[i][0]=false; //not from diagonal
    }

    for(j=0; j<=len2; j++)
    {
        val[0][j]=0;
        path[0][j]=false; //not from diagonal
        j2i[j]=-1;    //all are not aligned, only use j2i[1:len2]
    }

    //decide matrix and path
    for(i=1; i<=len1; i++)
    {
        for(j=1; j<=len2; j++)
        {
            //d=val[i-1][j-1]+score[i][j]; //diagonal
            if(secx[i-1]==secy[j-1])
            {
                d=val[i-1][j-1] + 1.0;
            }
            else
            {
                d=val[i-1][j-1];
            }

            //symbol insertion in horizontal (= a gap in vertical)
            h=val[i-1][j];
            h += (path[i-1][j]) ? gap_open : 0.0f; //aligned in last position


            //symbol insertion in vertical
            v=val[i][j-1];
            v += (path[i][j-1]) ? gap_open : 0.0f; //aligned in last position


            path[i][j]= (d>=h && d>=v) ? true : false; //from diagonal
            val[i][j]=(v>=h) ? v : h;
            val[i][j]=(path[i][j] == true) ? d : val[i][j];
        } //for i
    } //for j
    std::cout << val[i-1][j-1] << std::endl;
    //trace back to extract the alignment
    i=len1;
    j=len2;
    while(i>0 && j>0)
    {
        if(path[i][j]) //from diagonal
        {
            j2i[j-1]=i-1;
            i--;
            j--;
        }
        else
        {
            h=val[i-1][j];
            if(path[i-1][j]) h +=gap_open;

            v=val[i][j-1];
            if(path[i][j-1]) v +=gap_open;

            if(v>=h)
                j--;
            else
                i--;
        }
    }
}
