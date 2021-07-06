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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <map>

#include "simd.h"
#include "Coordinates.h"
using namespace std;

class BasicFunction{
public:
    static void PrintErrorAndQuit(string sErrorString)
    {
        cout << sErrorString << endl;
        exit(1);
    }

    static char AAmap(string AA)
    {
        if (AA.compare("ALA")==0) return 'A';
        if (AA.compare("ASX")==0) return 'B';
        if (AA.compare("CYS")==0) return 'C';
        if (AA.compare("ASP")==0) return 'D';
        if (AA.compare("GLU")==0) return 'E';
        if (AA.compare("PHE")==0) return 'F';
        if (AA.compare("GLY")==0) return 'G';
        if (AA.compare("HIS")==0) return 'H';
        if (AA.compare("ILE")==0) return 'I';
        if (AA.compare("LYS")==0) return 'K';
        if (AA.compare("LEU")==0) return 'L';
        if (AA.compare("MET")==0 || AA.compare("MSE")==0) return 'M';
        if (AA.compare("ASN")==0) return 'N';
        if (AA.compare("PYL")==0) return 'O';
        if (AA.compare("PRO")==0) return 'P';
        if (AA.compare("GLN")==0) return 'Q';
        if (AA.compare("ARG")==0) return 'R';
        if (AA.compare("SER")==0) return 'S';
        if (AA.compare("THR")==0) return 'T';
        if (AA.compare("SEC")==0) return 'U';
        if (AA.compare("VAL")==0) return 'V';
        if (AA.compare("TRP")==0) return 'W';
        if (AA.compare("TYR")==0) return 'Y';
        if (AA.compare("GLX")==0) return 'Z';

        if (AA.compare(0,2," D")==0) return tolower(AA[2]);
        if (AA.compare(0,2,"  ")==0) return tolower(AA[2]);
        return 'X';
    }

    static void get_xyz(string line, double *x, double *y, double *z, char *resname, int *no)
    {
        char cstr[50];

        strcpy(cstr, (line.substr(30, 8)).c_str());
        sscanf(cstr, "%lf", x);

        strcpy(cstr, (line.substr(38, 8)).c_str());
        sscanf(cstr, "%lf", y);

        strcpy(cstr, (line.substr(46, 8)).c_str());
        sscanf(cstr, "%lf", z);

        strcpy(cstr, (line.substr(17, 3)).c_str());
        *resname = AAmap(cstr);

        strcpy(cstr, (line.substr(22, 4)).c_str());
        sscanf(cstr, "%d", no);
    }

    static int read_PDB(const vector<string> &PDB_lines, double **a, char *seq, int *resno)
    {
        size_t i;
        for (i=0;i<PDB_lines.size();i++){
            get_xyz(PDB_lines[i], &a[i][0], &a[i][1], &a[i][2], &seq[i], &resno[i]);
        }
        seq[i]='\0';
        return i;
    }

    static float dist(const float xx, const float xy, const float xz,
               const float yx, const float yy, const float yz)
    {
        float d1=xx-yx;
        float d2=xy-yy;
        float d3=xz-yz;
        return (d1*d1 + d2*d2 + d3*d3);
    }

    static double dotProd(const float *a, const float *b)
    {
        return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
    }

    static void transform(float t[3], float u[3][3], const float x[3], float *x1)
    {
        x1[0]=t[0]+dotProd(&u[0][0], x);
        x1[1]=t[1]+dotProd(&u[1][0], x);
        x1[2]=t[2]+dotProd(&u[2][0], x);
    }

    static void transform(float t[3], float u[3][3], const float xx, const float xy, const float xz,
                   float &yx, float & yy, float & yz)
    {
        float xyz[3];
        xyz[0] = xx;
        xyz[1] = xy;
        xyz[2] = xz;
        yx=t[0]+dotProd(&u[0][0], xyz);
        yy=t[1]+dotProd(&u[1][0], xyz);
        yz=t[2]+dotProd(&u[2][0], xyz);
    }

    static void do_rotation(const float *x, float *x1, int len, float t[3], float u[3][3])
    {
        for(int i=0; i<len; i++)
        {
            transform(t, u, &x[i*3], &x1[i*3]);
        }
    }

    static void do_rotation( const Coordinates & x,  Coordinates & y,
                      int len, float t[3], float u[3][3])
    {
        simd_float t0 = simdf32_set(t[0]);
        simd_float t1 = simdf32_set(t[1]);
        simd_float t2 = simdf32_set(t[2]);
//
        simd_float u00 = simdf32_set(u[0][0]);
        simd_float u01 = simdf32_set(u[0][1]);
        simd_float u02 = simdf32_set(u[0][2]);
        simd_float u10 = simdf32_set(u[1][0]);
        simd_float u11 = simdf32_set(u[1][1]);
        simd_float u12 = simdf32_set(u[1][2]);
        simd_float u20 = simdf32_set(u[2][0]);
        simd_float u21 = simdf32_set(u[2][1]);
        simd_float u22 = simdf32_set(u[2][2]);
        for(int i=0; i<len; i+=VECSIZE_FLOAT)
//        for(int i=0; i<len; i++)

        {
//        float xyz[3];
//        xyz[0] = xx;
//        xyz[1] = xy;
//        xyz[2] = xz;
//        yx=t[0]+dotProd(&u[0][0], xyz);
//        yy=t[1]+dotProd(&u[1][0], xyz);
//        yz=t[2]+dotProd(&u[2][0], xyz);
            simd_float x_x = simdf32_load(&x.x[i]);
            simd_float x_y = simdf32_load(&x.y[i]);
            simd_float x_z = simdf32_load(&x.z[i]);
            simd_float xx = simdf32_mul(u00, x_x);
            simd_float yy = simdf32_mul(u01, x_y);
            simd_float zz = simdf32_mul(u02, x_z);
            xx = simdf32_add(xx, yy);
            zz = simdf32_add(xx, zz);
            simdf32_store(&y.x[i], simdf32_add(t0, zz));
            xx = simdf32_mul(u10, x_x);
            yy = simdf32_mul(u11, x_y);
            zz = simdf32_mul(u12, x_z);
            xx = simdf32_add(xx, yy);
            zz = simdf32_add(xx, zz);
            simdf32_store(&y.y[i], simdf32_add(t1, zz));
            xx = simdf32_mul(u20, x_x);
            yy = simdf32_mul(u21, x_y);
            zz = simdf32_mul(u22, x_z);
            xx = simdf32_add(xx, yy);
            zz = simdf32_add(xx, zz);
            simdf32_store(&y.z[i], simdf32_add(t2, zz));

//        transform(t, u, x.x[i], x.y[i], x.z[i], y.x[i], y.y[i], y.z[i]);
        }
    }

};