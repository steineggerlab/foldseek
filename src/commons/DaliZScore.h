#ifndef DALI_Z_H
#define DALI_Z_H

#include <algorithm>
#include <cmath>
#include <string>
#include <iostream> // for debug

class DaliCalculator {
public:
    static constexpr float DALI_THETA_E = 0.2;
    static constexpr float DALI_ALPHA = 400.0;
    static constexpr int DALI_LMAX = 400;
    static constexpr float DALI_C1 = 7.9494;
    static constexpr float DALI_C2 = 0.70852;
    static constexpr float DALI_C3 = 2.5898/10000;
    static constexpr float DALI_C4 = -1.9156/1000000;
    static constexpr float DALI_EPS = 0.000001;

    DaliCalculator(unsigned int maxQueryLength, unsigned int maxTargetLength);
    ~DaliCalculator();

    struct DaliScoreResult {
        DaliScoreResult() {
            daliZScore = 0.0;
        }
        DaliScoreResult(float score) {
            daliZScore = score;
        }
        float daliZScore; 
    };

    float dist(float* arr1, float* arr2);
    void initQuery(unsigned int queryLen, float *qx, float *qy, float *qz);
    void constructAlignHashes(int align_idx, int query_idx, int target_idx);
    DaliScoreResult computeDaliScore(unsigned int targetLen, int qStartPos, int tStartPos, const std::string &backtrace, float *tx, float *ty, float *tz);

private:
    unsigned int queryStart, targetStart, queryLength, targetLength, alignLength;
    unsigned int maxQueryLength, maxTargetLength, maxAlignLength;
    float **query_pos, **target_pos, **query_distmat, **target_distmat;
    float *align_to_query, *align_to_target, *dali;
    std::string cigar; // backtrace
};



#endif