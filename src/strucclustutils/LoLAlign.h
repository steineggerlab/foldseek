

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
    lolAlign(unsigned int maxSeqLen, bool tmAlignFast, bool tmScoreOnly, bool exact);
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

    void initQuery(float *x, float *y, float *z, char * querySeq, char* query3diSeq, unsigned int queryLen);
    
    TMscoreResult computeTMscore(float *x, float *y, float *z,
                                 unsigned int targetLen, int qStartPos,
                                 int targetStartPos, const std::string & backtrace,
                                 int normalizationLen);

    Matcher::result_t align(unsigned int dbKey, float *target_x, float *target_y, float *target_z,
                            char * targetSeq, char* target3diSeq, unsigned int targetLen);

    static unsigned int normalization(int mode, unsigned int alignmentLen, unsigned int queryLen, unsigned int targetLen);

private:
    AffineNeedlemanWunsch * affineNW;
    std::string backtrace;
    float * query_x;
    float * query_y;
    float * query_z;
    float * target_x;
    float * target_y;
    float * target_z;
    char * querySecStruc;
    char * targetSecStruc;
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
    bool tmAlignFast;
    Coordinates xtm, ytm, xt, r1, r2;
    bool computeExactScore;
    int * invmap;
    void lolmatrix(int *anchor, int anchor_length, int *gaps, float **d_ij, float **d_kl, float **G);
    void lolscore(std::vector<float> d_ij, std::vector<float> d_kl, std::vector<float> d_seq, std::vector<float> score); 
    float alignment_lolScore(std::vector<float> d_ij, std::vector<float> d_kl, std::vector<float> anchorpoints, size_t anchor_length);

    TMscoreResult computeExactTMscore(float *x, float *y, float *z,
                                      unsigned int targetLen, int qStartPos,
                                      int targetStartPos, const std::string & backtrace,
                                      int normalizationLen);

    TMscoreResult computeAppoximateTMscore(float *x, float *y, float *z,
                                      unsigned int targetLen, int qStartPos,
                                      int targetStartPos, const std::string & backtrace,
                                      int normalizationLen);
};



