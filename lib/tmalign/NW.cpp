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
/*    Please note this function is not a correct implementation of
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

