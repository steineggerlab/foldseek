

#include "SubstitutionMatrix.h"
#include "IndexReader.h"
#include "DBReader.h"

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
        char* queryNumAA,
        char* queryNum3Di,
        char* targetNumAA,
        char* targetNum3Di,
        int queryLen,
        int targetLen,
        SubstitutionMatrix &subMatAA,
        SubstitutionMatrix &subMat3Di,
        float T,
        float** scoreForward
    );


    unsigned char* seq2num(const std::string& seq, const unsigned char* aa2num) {
        unsigned char* idx = static_cast<unsigned char*>(malloc(seq.size() * sizeof(unsigned char)));
        for (size_t i = 0; i < seq.size(); ++i) {
            idx[i] = aa2num[static_cast<unsigned char>(seq[i])];
        }
        return idx;
    }

    void initQuery(float *x, float *y, float *z, char * querySeq, char* query3diSeq, unsigned int queryLen);
    
    TMscoreResult computeTMscore(float *x, float *y, float *z,
                                 unsigned int targetLen, int qStartPos,
                                 int targetStartPos, const std::string & backtrace,
                                 int normalizationLen);

    Matcher::result_t align(unsigned int dbKey, float *target_x, float *target_y, float *target_z,
                            char * targetSeq, char* target3diSeq, unsigned int targetLen);

    static unsigned int normalization(int mode, unsigned int alignmentLen, unsigned int queryLen, unsigned int targetLen);
    void lol_fwbw(float** scoreForward, float** P,
        size_t queryLen, size_t targetLen,
        float go, float ge, float T, int length, int blocks);

    float** allocateMemory(size_t queryLen, size_t targetLen);
    void rescaleBlocks(float **matrix, float **scale, size_t rows, size_t length, size_t blocks, size_t targetLen);
    void forwardBackwardSaveBlockMaxLocal(float** S, float** z_init,
                                            float T, float go, float ge,
                                            size_t rows, size_t start, size_t end, size_t memcpy_cols, size_t targetlen,
                                            float** zm, float* zmax, float** zmBlock, float* zeBlock, float* zfBlock);

private:




    std::string backtrace;
    float * query_x;
    float * query_y;
    float * query_z;
    float * target_x;
    float * target_y;
    float * target_z;
    char * querySecStruc;
    char * targetSecStruc;
    float ** scoreForward;
    float ** P;
    float *mem;
    float w1[3][2] = {
        {-1.3584513e-04,  9.9329501e-01},
        { 7.6149112e-01,  5.7029408e-01},
        {-8.1348085e-01,  6.0702705e-01}
    };
    float b1[3] = {0.7043129, 0.374659, 0.39905924};

    float w2[3] = {-0.776632, 0.61055756, 0.5823986};
    float b2 = -0.11200039;

    unsigned int queryLen;
    char * querySeq;
    char * query3diSeq;
    std::string seqM, seqxA, seqyA;// for output alignment

    Coordinates xtm, ytm, xt, r1, r2;
    bool computeExactScore;
    int * invmap;
    void lolmatrix(int *anchor, int anchor_length, int *gaps, float **d_ij, float **d_kl, float **G);
    void lolscore(std::vector<float> d_ij, std::vector<float> d_kl, std::vector<float> d_seq, std::vector<float> score); 
    float alignment_lolScore(std::vector<float> d_ij, std::vector<float> d_kl, std::vector<float> anchorpoints, size_t anchor_length);

    void align(char *targetSeq, char * target3diSeq);

};



