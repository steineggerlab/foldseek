#ifndef LoLAlign
#define LoLAlign

#include "SubstitutionMatrix.h"
#include "IndexReader.h"
#include "DBReader.h"
#include "Fwbw.h"

#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>

#include "Matcher.h"
#include "tmalign/Coordinates.h"

class lolAlign{
public:
    lolAlign(unsigned int maxSeqLen, bool exact);
    ~lolAlign();

    
    void computeForwardScoreMatrix(
        unsigned char * targetNumAA,
        unsigned char *targetNum3Di,
        int queryLen,
        int targetLen,
        SubstitutionMatrix &subMatAA,
        float** scoreForward
    );
    void addForwardScoreMatrix(
        unsigned char * targetNumAA,
        unsigned char *targetNum3Di,
        int queryLen,
        int targetLen,
        SubstitutionMatrix &subMatAA,
        float** scoreForward
    );
    void calc_gap(int* anchor_query, int* anchor_target, int * gaps,  int queryLen, int targetLen);

    void initQuery(float *x, float *y, float *z, Sequence& qSeqAA, Sequence& qSeq3Di, int queryLen, SubstitutionMatrix &subMatAA, int maxTLen, int md);
    
    Matcher::result_t align(unsigned int dbKey, float *target_x, float *target_y, float *target_z,
                            Sequence& tSeqAA, Sequence& tSeq3Di, int targetLen, SubstitutionMatrix &subMatAA, FwBwAligner* fwbwaln, int md);
    float maxSubArray(float* nums, int numsSize);

    void align_startAnchors(int * anchor_query, int * anchor_target, int max_query, int max_target, int * anchor_length, float** fwbwP, float** G);
    void index_sort(float* nums, int* index, int numsSize);


    void lol_fwbw(float** scoreForward, float** P,
        size_t queryLen, size_t targetLen,size_t assignTargetLen,
        float go, float ge, float T, int length, int blocks, int* gaps);

    void rescaleBlocks(float **matrix, float **scale, size_t rows, size_t length, size_t blocks, size_t targetLen);
    void forwardBackwardSaveBlockMaxLocal(float** S, float** z_init,
                                            float T, float go, float ge,
                                            size_t rows, size_t start, size_t end, size_t memcpy_cols, size_t targetlen,
                                            float** zm, float* zmax, float** zmBlock, float* zeBlock, float* zfBlock, int* gaps);
    void calc_dist_matrix(float *x, float *y, float *z, size_t len, float **d, bool cutoff);
    void reallocate_target(size_t targetL);
    float calc_discore(int * anchor_query, int * anchor_target, int anchor_length);
    void lolmatrix(int *anchor_query, int* anchor_target,int anchor_length, int *gaps, float **d_ij, float **d_kl, float **G, int queryLen, int targetLen, float ** hidden_layer, float * d_dist);
    void lolscore(float* dist, float* d_seq, float* score, int length, float** hidden_layer);
    void set_start_anchor_length(int length) {start_anchor_length = length;}
    

    void lolscore(float* d_dist, float d_seq, float* score, int length, int start, float** hidden_layer);
    void computeDi_score(
        unsigned char * targetNumAA,
        unsigned char *targetNum3Di,
        int anchorLen,
        int* final_anchor_query,
        int* final_anchor_target,
        SubstitutionMatrix &subMatAA,
        float *scoreForward);

private:
    std::string backtrace;
    float * query_x;
    float * query_y;
    float * query_z;
    float * target_x;
    float * target_y;
    float * target_z;
    float ** scoreForward;
    float ** P;
    float ** d_ij;
    float ** d_kl;
    float ** G;
    int ** anchor_query;
    int ** anchor_target;
    int num_sa = 10;
    int SeedNumber = 3;
    float start_anchor_go = -6.0;
    float start_anchor_ge = -3.0;
    float start_anchor_T = 2.0;
    int start_anchor_length = 3;
    float lol_go = -1.5;
    float lol_ge = -0.0;
    float lol_min_p = 0.4;
    float lol_T = 4;
    float** hidden_layer;
    int* sa_index = new int[num_sa];
    float* sa_scores = new float[num_sa];
    int* anchor_length = new int[num_sa];
    int* new_anchor_length = new int[num_sa];
    int* gaps = new int[4]{0, 0, 0, 0};
    float* lol_dist;
    float* lol_seq_dist;
    float* lol_score_vec;
    int* final_anchor_query;
    int* final_anchor_target;
    float QQ_score;

    unsigned char *queryNumAA;
    unsigned char *queryNum3Di;
    float* lol_score_vec_sh;

    float w1[2][3] = {
        {-1.3584513e-04,7.6149112e-01,-8.1348085e-01 },
        {9.9329501e-01 , 5.7029408e-01, 6.0702705e-01}
    };
    float b1[3] = {0.7043129 , 0.374659  , 0.39905924};

    float w2[3] = {-0.776632  ,  0.61055756, 0.5823986};
    float b2 = -0.11200039 +0.5;

    // Load weights and biases into SIMD registers
    simd_float w1_0 = simdf32_set(w1[0][0]); 
    simd_float w1_1 = simdf32_set(w1[0][1]); 
    simd_float w1_2 = simdf32_set(w1[0][2]); 

    simd_float w1_d0 = simdf32_set(w1[1][0]); 
    simd_float w1_d1 = simdf32_set(w1[1][1]); 
    simd_float w1_d2 = simdf32_set(w1[1][2]);

    simd_float b1_0 = simdf32_set(b1[0]); 
    simd_float b1_1 = simdf32_set(b1[1]);
    simd_float b1_2 = simdf32_set(b1[2]); 

    simd_float w2_0 = simdf32_set(w2[0]); 
    simd_float w2_1 = simdf32_set(w2[1]);
    simd_float w2_2 = simdf32_set(w2[2]); 

    simd_float b2_vec = simdf32_set(b2); 

    simd_float zero = simdf32_setzero();
    int min_lolmat_idx;
    int max_lolmat_idx;

    int queryLen;
    const char * querySeq;
    const char * query3diSeq;
    std::string seqM, seqxA, seqyA;// for output alignment

    Coordinates xtm, ytm, xt, r1, r2;
    bool computeExactScore;
    int * invmap;
    const float scoringMatrix3Di[20][20] = {
        // A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T
        {10, -1,  1,  7,  6,  2,  2, -4,  1, -1, -6, -2, -1,  4, -1, -5, -1, -5,  3,  2}, // A
        {-1,  7, -4, -5, -3, -2, -1, -10, -10,  1, -11,  1, -2,  2,  1,  1, -6,  1, -5, -5}, // B
        {1, -4,  0, -3, -1,  0,  1, -4, -5, -6, -5, -4, -3, -3, -2, -6, -2, -6, -2, -2}, // C
        {7, -5, -3, 15,  3, -2,  1, -8, -5, -4, -12, -5, -5,  1, -3, -7, -5, -10, -2,  3}, // D
        {6, -3, -1,  3, 10,  1,  1, -2,  4, -1, -5, -3, -2,  4, -2, -6,  0, -6,  7,  0}, // E
        {2, -2,  0, -2,  1, 10,  7,  3, -3, -5,  3,  1, -2, -1,  5, -1,  7, -5, -1,  3}, // F
        {2, -1,  1,  1,  1,  7, 10,  0, -4, -4, -1,  3, -1,  0,  4, -1,  3, -4, -1,  8}, // G
        {-4, -10, -4, -8, -2,  3,  0, 11, -1, -10, 11, -5, -6, -3,  0, -7,  9, -11, -2, -3}, // H
        {1, -10, -5, -5,  4, -3, -4, -1, 13, -9, -4, -9, -7, -2, -6, -12, -2, -13,  9, -3}, // I
        {-1,  1, -6, -4, -1, -5, -4, -10, -9,  6, -13, -2, -3,  3, -2, -3, -7,  0, -6, -6}, // J
        {-6, -11, -5, -12, -5,  3, -1, 11, -4, -13, 15, -6, -8, -7, -1, -8,  7, -14, -2, -4}, // K
        {-2,  1, -4, -5, -3,  1,  3, -5, -9, -2, -6,  8, -1, -1,  4,  4, -2,  0, -6, -1}, // L
        {-1, -2, -3, -5, -2, -2, -1, -6, -7, -3, -8, -1,  1, -1, -1, -3, -4, -2, -5, -5}, // M
        {4,  2, -3,  1,  4, -1,  0, -3, -2,  3, -7, -1, -1,  7,  0, -2, -3, -1,  1, -1}, // N
        {-1,  1, -2, -3, -2,  5,  4,  0, -6, -2, -1,  4, -1,  0,  8,  4,  3, -1, -3,  1}, // O
        {-5,  1, -6, -7, -6, -1, -1, -7, -12, -3, -8,  4, -3, -2,  4,  6, -4,  0, -9, -5}, // P
        {-1, -6, -2, -5,  0,  7,  3,  9, -2, -7,  7, -2, -4, -3,  3, -4, 11, -8, -2,  0}, // Q
        {-5,  1, -6, -10, -6, -5, -4, -11, -13,  0, -14,  0, -2, -1, -1,  0, -8,  2, -9, -9}, // R
        {3, -5, -2, -2,  7, -1, -1, -2,  9, -6, -2, -6, -5,  1, -3, -9, -2, -9, 11, -2}, // S
        {2, -5, -2,  3,  0,  3,  8, -3, -3, -6, -4, -1, -5, -1,  1, -5,  0, -9, -2, 14}  // T
    };
};

#endif