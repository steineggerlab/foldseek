

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
#include "tmalign/affineneedlemanwunsch.h"
#include "tmalign/Coordinates.h"

class lolAlign{
public:
    lolAlign(unsigned int maxSeqLen, bool exact);
    ~lolAlign();

    struct TMscoreResult{
        TMscoreResult(){
            memset(&u, 0, 3*3*sizeof(float));
            memset(&t, 0, 3*sizeof(float));
            this->tmscore = 0.0;
            this->rmsd = 0.0;
        }

        TMscoreResult(float u[3][3], float t[3], double tmscore, double rmsd) {
            memcpy(this->u,u,3*3*sizeof(float));
            memcpy(this->t,t,3*sizeof(float));
            this->tmscore = tmscore;
            this->rmsd = rmsd;
        }
        float u[3][3];
        float t[3];
        double tmscore;
        double rmsd;
    };
    void computeForwardScoreMatrix(
        unsigned char * targetNumAA,
        unsigned char *targetNum3Di,
        int queryLen,
        int targetLen,
        SubstitutionMatrix &subMatAA,
        SubstitutionMatrix &subMat3Di,
        float** scoreForward
    );
    void calc_gap(int* anchor_query, int* anchor_target, int * gaps,  int queryLen, int targetLen);


    unsigned char* seq2num(const std::string& seq, const unsigned char* aa2num) {
        unsigned char* idx = static_cast<unsigned char*>(malloc(seq.size() * sizeof(unsigned char)));
        for (size_t i = 0; i < seq.size(); ++i) {
            idx[i] = aa2num[static_cast<unsigned char>(seq[i])];
        }
        return idx;
    }


    void initQuery(float *x, float *y, float *z, char * querySeq, char* query3diSeq, int queryLen, int maxTLen, SubstitutionMatrix &subMatAA, SubstitutionMatrix &subMat3Di);
    


    Matcher::result_t align(unsigned int dbKey, float *target_x, float *target_y, float *target_z,
                            char * targetSeq, char* target3diSeq, int targetLen, SubstitutionMatrix &subMatAA, SubstitutionMatrix &subMat3Di, FwBwAligner* fwbwaln);
    float maxSubArray(float* nums, int numsSize);

    void align_startAnchors(int * anchor_query, int * anchor_target, int max_query, int max_target, int * anchor_length, float** P, float** G);
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
        SubstitutionMatrix &subMat3Di,
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
    float lol_go = -5.0;
    float lol_ge = -0.0;
    float lol_min_p = 0.4;
    float lol_T = 3;
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
    char * querySeq;
    char * query3diSeq;
    std::string seqM, seqxA, seqyA;// for output alignment

    Coordinates xtm, ytm, xt, r1, r2;
    bool computeExactScore;
    int * invmap;
};



