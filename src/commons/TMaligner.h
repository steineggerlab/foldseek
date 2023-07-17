//
// Created by Martin Steinegger on 07/07/2022.
//

#ifndef FOLDSEEK_TMALIGNER_H
#define FOLDSEEK_TMALIGNER_H


#include "Matcher.h"
#include "tmalign/affineneedlemanwunsch.h"
#include "tmalign/Coordinates.h"

class TMaligner{
public:
    TMaligner(unsigned int maxSeqLen, bool tmAlignFast, bool tmScoreOnly);
    ~TMaligner();

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
    void initQuery(float * x, float * y, float * z, char * querySeq, unsigned int queryLen);
    TMscoreResult computeTMscore(float *x, float *y, float *z,
                                 unsigned int targetLen, int qStartPos,
                                 int targetStartPos, const std::string & backtrace,
                                 int normalizationLen);
    Matcher::result_t align(unsigned int dbKey, float *target_x, float *target_y, float *target_z,
                            char * targetSeq, unsigned int targetLen, float &TM);

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
    unsigned int queryLen;
    char * querySeq;
    std::string seqM, seqxA, seqyA;// for output alignment
    bool tmAlignFast;
    Coordinates xtm, ytm, xt, r1, r2;
    int * invmap;
};

#endif //FOLDSEEK_TMALIGNER_H
