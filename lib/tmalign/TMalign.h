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
#include "affineneedlemanwunsch.h"
#include <string>

void parameter_set4search(const int xlen, const int ylen,
                          float &D0_MIN, float &Lnorm,
                          float &score_d8, float &d0, float &d0_search, float &dcu0);

void parameter_set4final(const float len, float &D0_MIN, float &Lnorm,
                         float &d0, float &d0_search);

void parameter_set4scale(const int len, const float d_s, float &Lnorm,
                         float &d0, float &d0_search);

//     1, collect those residues with dis<d;
//     2, calculate TMscore
int score_fun8( Coordinates &xa, Coordinates &ya, int n_ali, float d, int i_ali[],
                float *score1, const float Lnorm,
                const float score_d8, const float d0, float * mem);

int score_fun8_standard(Coordinates &xa, Coordinates &ya, int n_ali, float d,
                        int i_ali[], float *score1, int score_sum_method,
                        float score_d8, float d0);


bool KabschFast(Coordinates & x,
                 Coordinates & y,
                 int n,
                 float *rms,
                 float t[3],
                 float u[3][3],
                 float * mem);

double TMscore8_search(Coordinates &r1, Coordinates &r2,
                       Coordinates &xtm, Coordinates & ytm,
                       Coordinates &xt, int Lali, float t0[3], float u0[3][3], int simplify_step,
                       float *Rcomm, float local_d0_search, float Lnorm,
                       float score_d8, float d0, float * mem);


double TMscore8_search_standard(Coordinates &r1, Coordinates &r2,
                                Coordinates &xtm, Coordinates &ytm, Coordinates &xt, int Lali,
                                float t0[3], float u0[3][3], int simplify_step, int score_sum_method,
                                float *Rcomm, float local_d0_search, float score_d8, float d0, float * mem);

//Comprehensive TMscore search engine
// input:   two vector sets: x, y
//          an alignment invmap0[] between x and y
//          simplify_step: 1 or 40 or other integers
//          score_sum_method: 0 for score over all pairs
//                            8 for socre over the pairs with dist<score_d8
// output:  the best rotaion matrix t, u that results in highest TMscore
double detailed_search(Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                       Coordinates &xt, const Coordinates &x, const Coordinates &y, int ylen,
                       int invmap0[], float t[3], float u[3][3], int simplify_step,
                       float local_d0_search, float Lnorm, float score_d8, float d0, float * mem);

double detailed_search_standard( Coordinates &r1, Coordinates &r2,
                                 Coordinates &xtm, Coordinates &ytm, Coordinates &xt,
                                 const Coordinates &x, const Coordinates &y,
                                 int ylen, int invmap0[], float t[3], float u[3][3],
                                 int simplify_step, int score_sum_method, double local_d0_search,
                                 const bool& bNormalize, float Lnorm, float score_d8, float d0, float * mem);
//compute the score quickly in three iterations
double get_score_fast( Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                       const Coordinates &x, const Coordinates &y, int ylen, int invmap[],
                       float d0, float d0_search, float t[3], float u[3][3], float * mem);

int sec_str(float dis13, float dis14, float dis15,
            float dis24, float dis25, float dis35);

//1->coil, 2->helix, 3->turn, 4->strand
void make_sec(Coordinates &x, int len, char *sec);


//get initial alignment from secondary structure alignment
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g.,
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
void get_initial_ss(AffineNeedlemanWunsch *affineNW,
                    const char *secx, const char *secy, int xlen, int ylen, int *y2x);


// get_initial5 in TMalign fortran, get_intial_local in TMalign c by yangji
//get initial alignment of local structure superposition
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g.,
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0
//the jth element in y is aligned to a gap in x if i==-1
bool get_initial5(AffineNeedlemanWunsch *affineNW,
                   Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                   const Coordinates &x, const Coordinates &y, int xlen, int ylen, int *y2x,
                   float d0, float d0_search, const bool fast_opt, const float D0_MIN, float * mem);


void find_max_frag(const Coordinates &x, int len, int *start_max,
                   int *end_max, float dcu0, const bool fast_opt);

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
                       float * mem);

//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap
double DP_iter(AffineNeedlemanWunsch * affineNW, Coordinates &r1, Coordinates &r2,
               Coordinates &xtm, Coordinates &ytm,
               Coordinates &xt, const Coordinates &x, const Coordinates &y,
               int xlen, int ylen, float t[3], float u[3][3],
               int invmap0[], int g1, int g2, int iteration_max, float local_d0_search,
               float Lnorm, float d0, float score_d8, float * mem);

double standard_TMscore(Coordinates &r1, Coordinates &r2, Coordinates &xtm, Coordinates &ytm,
                        Coordinates &xt, Coordinates &x, Coordinates &y, int ylen, int invmap[],
                        int& L_ali, float& RMSD, float D0_MIN, float Lnorm, float d0,
                        float score_d8, float t[3], float u[3][3], float * mem);

/* entry function for TMalign */
int TMalign_main(
        AffineNeedlemanWunsch * affineNW,
        const Coordinates &xa, const Coordinates &ya,
        const char *seqx, const char *seqy, const char *secx, const char *secy,
        float t0[3], float u0[3][3],
        float &TM1, float &TM2, float &TM3, float &TM4, float &TM5,
        float &d0_0, float &TM_0,
        float &d0A, float &d0B, float &d0u, float &d0a, float &d0_out,
        std::string &seqM, std::string &seqxA, std::string &seqyA,
        float &rmsd0, float &Liden, int &n_ali, int &n_ali8,
        const int xlen, const int ylen, const float Lnorm_ass,
        const float d0_scale, const bool I_opt, const bool a_opt,
        const bool u_opt, const bool d_opt, const bool fast_opt, float * mem);
