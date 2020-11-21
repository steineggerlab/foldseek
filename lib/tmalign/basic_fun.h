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
#include <math.h>
#include <time.h>
#include <string.h>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <iomanip>
#include <map>

#include "pstream.h"
#include "simd.h"

using namespace std;


struct Coordinates{
    Coordinates(int size){
        x =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
        y =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
        z =(float*) mem_align(ALIGN_FLOAT, (size+VECSIZE_FLOAT)*sizeof(float));
    }
    Coordinates(){
    }
    float * x;
    float * y;
    float * z;
};


void PrintErrorAndQuit(string sErrorString)
{
    cout << sErrorString << endl;
    exit(1);
}

template <typename T> inline T getmin(const T &a, const T &b)
{
    return b<a?b:a;
}

template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
    *array=new A* [Narray1];
    for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
};

template <class A> void DeleteArray(A *** array, int Narray)
{
    for(int i=0; i<Narray; i++)
        if(*(*array+i)) delete [] *(*array+i);
    if(Narray) delete [] (*array);
    (*array)=NULL;
};

string AAmap(char A)
{
    if (A=='A') return "ALA";
    if (A=='B') return "ASX";
    if (A=='C') return "CYS";
    if (A=='D') return "ASP";
    if (A=='E') return "GLU";
    if (A=='F') return "PHE";
    if (A=='G') return "GLY";
    if (A=='H') return "HIS";
    if (A=='I') return "ILE";
    if (A=='K') return "LYS";
    if (A=='L') return "LEU";
    if (A=='M') return "MET";
    if (A=='N') return "ASN";
    if (A=='O') return "PYL";
    if (A=='P') return "PRO";
    if (A=='Q') return "GLN";
    if (A=='R') return "ARG";
    if (A=='S') return "SER";
    if (A=='T') return "THR";
    if (A=='U') return "SEC";
    if (A=='V') return "VAL";
    if (A=='W') return "TRP";    
    if (A=='Y') return "TYR";
    if (A=='Z') return "GLX";
    if ('a'<=A && A<='z') return "  "+toupper(A);
    return "UNK";
}

char AAmap(string AA)
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

void get_xyz(string line, double *x, double *y, double *z, char *resname, int *no)
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

int get_PDB_lines(std::istream& fin, vector<vector<string> >&PDB_lines,
    vector<string> &chainID_list, vector<string> &resi_vec,
    const int byresi_opt, const int ter_opt=3, const int infmt_opt=0,
    const string atom_opt="auto", const int split_opt=0)
{
    int i=0; // resi
    string line;    
    char chainID=0;
    string resi="";
    bool select_atom=true;
    int model_idx=0;
    vector<string> tmp_str_vec;
    
    int compress_type=0; // uncompressed file


    if (infmt_opt==0) // PDB format
    {
        while (fin.good())
        {
            getline(fin, line);
            if (i > 0)
            {
                if      (ter_opt>=1 && line.compare(0,3,"END")==0) break;
                else if (ter_opt>=3 && line.compare(0,3,"TER")==0) break;
            }
            if (split_opt && line.compare(0,3,"END")==0) chainID=0;
            if (line.compare(0, 6, "ATOM  ")==0 && line.size()>=54 &&
               (line[16]==' ' || line[16]=='A'))
            {
                if (atom_opt=="auto")
                {
                    if (line[17]==' ' && (line[18]=='D'||line[18]==' '))
                         select_atom=(line.compare(12,4," C3'")==0);
                    else select_atom=(line.compare(12,4," CA ")==0);
                }
                else     select_atom=(line.compare(12,4,atom_opt)==0);
                if (select_atom)
                {
                    if (!chainID)
                    {
                        chainID=line[21];
                        model_idx++;
                        stringstream i8_stream;
                        i=0;
                        if (split_opt==2) // split by chain
                        {
                            if (chainID==' ')
                            {
                                if (ter_opt>=1) i8_stream << ":_";
                                else i8_stream<<':'<<model_idx<<":_";
                            }
                            else
                            {
                                if (ter_opt>=1) i8_stream << ':' << chainID;
                                else i8_stream<<':'<<model_idx<<':'<<chainID;
                            }
                            chainID_list.push_back(i8_stream.str());
                        }
                        else if (split_opt==1) // split by model
                        {
                            i8_stream << ':' << model_idx;
                            chainID_list.push_back(i8_stream.str());
                        }
                        PDB_lines.push_back(tmp_str_vec);
                    }
                    else if (ter_opt>=2 && chainID!=line[21]) break;
                    if (split_opt==2 && chainID!=line[21])
                    {
                        chainID=line[21];
                        i=0;
                        stringstream i8_stream;
                        if (chainID==' ')
                        {
                            if (ter_opt>=1) i8_stream << ":_";
                            else i8_stream<<':'<<model_idx<<":_";
                        }
                        else
                        {
                            if (ter_opt>=1) i8_stream << ':' << chainID;
                            else i8_stream<<':'<<model_idx<<':'<<chainID;
                        }
                        chainID_list.push_back(i8_stream.str());
                        PDB_lines.push_back(tmp_str_vec);
                    }

                    if (resi==line.substr(22,5))
                        cerr<<"Warning! Duplicated residue "<<resi<<endl;
                    resi=line.substr(22,5); // including insertion code
                    if (byresi_opt==1) resi_vec.push_back(resi);
                    if (byresi_opt>=2) resi_vec.push_back(resi+line[21]);

                    /* change residue index in line */
                    stringstream i8_stream;
                    i8_stream << setw(4) << i+1;
                    line=line.substr(0,22)+i8_stream.str()+line.substr(26,29);

                    PDB_lines.back().push_back(line);
                    i++;
                }
            }
        }
    }
    else if (infmt_opt==1) // SPICKER format
    {
        int L=0;
        float x,y,z;
        while (fin.good())
        {
            fin   >>L>>x>>y>>z;
            getline(fin, line);
            if (!fin.good()) break;
            model_idx++;
            stringstream i8_stream;
            i8_stream << ':' << model_idx;
            chainID_list.push_back(i8_stream.str());
            PDB_lines.push_back(tmp_str_vec);
            for (i=0;i<L;i++)
            {
                fin   >>x>>y>>z;
                stringstream i8_stream;
                i8_stream<<"ATOM   "<<setw(4)<<i+1<<"  CA  UNK  "<<setw(4)
                    <<i+1<<"    "<<setiosflags(ios::fixed)<<setprecision(3)
                    <<setw(8)<<x<<setw(8)<<y<<setw(8)<<z;
                line=i8_stream.str();
                if (byresi_opt==1) resi_vec.push_back(line.substr(22,5));
                if (byresi_opt>=2) resi_vec.push_back(line.substr(22,5)+' ');
                PDB_lines.back().push_back(line);
            }
            getline(fin, line);
        }
    }
    else if (infmt_opt==2) // xyz format
    {
        int L=0;
        char A;
        while (fin.good())
        {
            getline(fin, line);
            L=atoi(line.c_str());
            getline(fin, line);
            for (i=0;i<line.size();i++)
                if (line[i]==' '||line[i]=='\t') break;
            if (!fin.good()) break;
            chainID_list.push_back(':'+line.substr(0,i));
            PDB_lines.push_back(tmp_str_vec);
            for (i=0;i<L;i++)
            {
                getline(fin, line);
                stringstream i8_stream;
                i8_stream<<"ATOM   "<<setw(4)<<i+1<<"  CA  "
                    <<AAmap(line[0])<<"  "<<setw(4)<<i+1<<"    "
                    <<line.substr(2,8)<<line.substr(11,8)<<line.substr(20,8);
                line=i8_stream.str();
                if (byresi_opt==1) resi_vec.push_back(line.substr(22,5));
                if (byresi_opt>=2) resi_vec.push_back(line.substr(22,5)+' ');
                PDB_lines.back().push_back(line);
            }
        }
    }
    line.clear();
    if (!split_opt) chainID_list.push_back("");
    return PDB_lines.size();
}

/* extract pairwise sequence alignment from residue index vectors,
 * assuming that "sequence" contains two empty strings.
 * return length of alignment, including gap. */
int extract_aln_from_resi(vector<string> &sequence, char *seqx, char *seqy,
    const vector<string> resi_vec1, const vector<string> resi_vec2,
    const int byresi_opt)
{
    int i1=0; // positions in resi_vec1
    int i2=0; // positions in resi_vec2
    int xlen=resi_vec1.size();
    int ylen=resi_vec2.size();
    map<char,int> chainID_map1;
    map<char,int> chainID_map2;
    if (byresi_opt==3)
    {
        vector<char> chainID_vec;
        char chainID;
        int i;
        for (i=0;i<xlen;i++)
        {
            chainID=resi_vec1[i][5];
            if (!chainID_vec.size()|| chainID_vec.back()!=chainID)
            {
                chainID_vec.push_back(chainID);
                chainID_map1[chainID]=chainID_vec.size();
            }
        }
        chainID_vec.clear();
        for (i=0;i<ylen;i++)
        {
            chainID=resi_vec2[i][5];
            if (!chainID_vec.size()|| chainID_vec.back()!=chainID)
            {
                chainID_vec.push_back(chainID);
                chainID_map2[chainID]=chainID_vec.size();
            }
        }
        chainID_vec.clear();
    }
    while(i1<xlen && i2<ylen)
    {
        if ((byresi_opt<=2 && resi_vec1[i1]==resi_vec2[i2]) || (byresi_opt==3
             && resi_vec1[i1].substr(0,5)==resi_vec2[i2].substr(0,5)
             && chainID_map1[resi_vec1[i1][5]]==chainID_map2[resi_vec2[i2][5]]))
        {
            sequence[0]+=seqx[i1++];
            sequence[1]+=seqy[i2++];
        }
        else if (atoi(resi_vec1[i1].substr(0,4).c_str())<=
                 atoi(resi_vec2[i2].substr(0,4).c_str()))
        {
            sequence[0]+=seqx[i1++];
            sequence[1]+='-';
        }
        else
        {
            sequence[0]+='-';
            sequence[1]+=seqy[i2++];
        }
    }
    chainID_map1.clear();
    chainID_map2.clear();
    return sequence[0].size();
}

int read_PDB(const vector<string> &PDB_lines, double **a, char *seq, int *resno)
{
    int i;
    for (i=0;i<PDB_lines.size();i++)
        get_xyz(PDB_lines[i], &a[i][0], &a[i][1], &a[i][2], &seq[i], &resno[i]);
    seq[i]='\0'; 
    return i;
}

float dist(const float xx, const float xy, const float xz,
           const float yx, const float yy, const float yz)
{
    float d1=xx-yx;
    float d2=xy-yy;
    float d3=xz-yz;
    return (d1*d1 + d2*d2 + d3*d3);
}

double dotProd(const float *a, const float *b)
{
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

void transform(float t[3], float u[3][3], const float x[3], float *x1)
{
    x1[0]=t[0]+dotProd(&u[0][0], x);
    x1[1]=t[1]+dotProd(&u[1][0], x);
    x1[2]=t[2]+dotProd(&u[2][0], x);
}

void transform(float t[3], float u[3][3], const float xx, const float xy, const float xz,
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

void do_rotation(const float *x, float *x1, int len, float t[3], float u[3][3])
{
    for(int i=0; i<len; i++)
    {
        transform(t, u, &x[i*3], &x1[i*3]);
    }    
}

void do_rotation( const Coordinates & x,  Coordinates & y,
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

/* strip white space at the begining or end of string */
string Trim(string inputString)
{
    string result = inputString;
    int idxBegin = inputString.find_first_not_of(" \n\r\t");
    int idxEnd = inputString.find_last_not_of(" \n\r\t");
    if (idxBegin >= 0 && idxEnd >= 0)
        result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
    return result;
}
