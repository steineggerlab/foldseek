#include "DaliZScore.h"
float dist(float* arr1, float* arr2) {
    float D2 = 0;
    for(int i = 0; i < 3; i++) {
        D2 += (arr1[i] - arr2[i]) * (arr1[i] - arr2[i]);
    }
    return sqrtf(D2);
}

DaliCalculator::DaliCalculator(unsigned int maxQueryLength, unsigned int maxTargetLength)
    : maxQueryLength(maxQueryLength), maxTargetLength(maxTargetLength) {
    maxAlignLength = std::max(maxQueryLength, maxTargetLength);

    query_pos = new float*[maxQueryLength];
    for(unsigned int i = 0; i < maxQueryLength; i++) {
        query_pos[i] = new float[3];
    }
    target_pos = new float*[maxTargetLength];
    for(unsigned int i = 0; i < maxTargetLength; i++) {
        target_pos[i] = new float[3];
    }

    align_to_query = new float[maxAlignLength];
    align_to_target = new float[maxAlignLength];

    query_distmat = new float*[maxQueryLength];
    for(unsigned int i = 0; i < maxQueryLength; i++) {
        query_distmat[i] = new float[i];
    }
    target_distmat = new float*[maxTargetLength];
    for(unsigned int i = 0; i < maxTargetLength; i++) {
        target_distmat[i] = new float[i];
    }

    dali = new float[maxQueryLength];
    for(unsigned int i = 0; i < maxTargetLength; i++) {
        dali[i] = 0.0;
    }
}

DaliCalculator::~DaliCalculator() {
    if(query_pos) {
        for(unsigned int i = 0; i < maxQueryLength; i++) {
            delete[] query_pos[i];
        }
        delete[] query_pos;
    }
    if(target_pos) {
        for(unsigned int i = 0; i < maxTargetLength; i++) {
            delete[] target_pos[i];
        }
        delete[] target_pos;
    }
    if(align_to_query) {
        delete[] align_to_query;
    }
    if(align_to_target) {
        delete[] align_to_target;
    }
    if(query_distmat) {
        for(unsigned int i = 0; i < maxQueryLength; i++) {
            delete[] query_distmat[i];
        }
        delete[] query_distmat;
    }
    if(target_distmat) {
        for(unsigned int i = 0; i < maxTargetLength; i++) {
            delete[] target_distmat[i];
        }
        delete[] target_distmat;
    }
    if(dali) {
        delete[] dali;
    }
}

void DaliCalculator::initQuery(unsigned int queryLen, float *qx, float *qy, float *qz) {
    queryLength = queryLen;
    for(unsigned int i = 0; i < queryLength; i++) {
        query_pos[i][0] = qx[i];
        query_pos[i][1] = qy[i];
        query_pos[i][2] = qz[i];
    }

    for(unsigned int i = 0; i < queryLength; i++) {
        for(unsigned int j = 0; j < queryLength; j++) {
            query_distmat[i][j] = dist(query_pos[i], query_pos[j]);
        }
    }
}

void DaliCalculator::constructAlignHashes(int align_idx, int query_idx, int target_idx) {
    for(std::size_t i = 0; i < cigar.length(); i++) {
        if(cigar[i] == 'M') {
            align_to_query[align_idx] = query_idx;
            align_to_target[align_idx] = target_idx;
            align_idx++; query_idx++; target_idx++;
        } else if(cigar[i] == 'D') {
            target_idx++;
        } else if(cigar[i] == 'I') {
            query_idx++;
        }
    }
    alignLength = align_idx;
}

DaliCalculator::DaliScoreResult DaliCalculator::computeDaliScore(unsigned int targetLen, int qStartPos, int tStartPos, const std::string &backtrace, float *tx, float *ty, float *tz) {
    targetLength = targetLen; 
    queryStart = qStartPos; 
    targetStart = tStartPos; 
    cigar = backtrace;

    for(unsigned int i = 0; i < targetLength; i++) {
        target_pos[i][0] = tx[i];
        target_pos[i][1] = ty[i];
        target_pos[i][2] = tz[i];
    }

    for(unsigned int i = 0; i < targetLength; i++) {
        for(unsigned int j = 0; j < targetLength; j++) {
            target_distmat[i][j] = dist(target_pos[i], target_pos[j]);
        }
    }

    constructAlignHashes(0, qStartPos, tStartPos);

    float daliScore = 0.0;
    for(unsigned int i = 0; i < alignLength; i++) {
        int queryIndex = align_to_query[i];
        int targetIndex = align_to_target[i];
        for(unsigned int j = 0; j < alignLength; j++) {
            if(i == j) {
                dali[i] += DALI_THETA_E;
                continue;
            }
            int queryIndex2 = align_to_query[j];
            int targetIndex2 = align_to_target[j];

            float dist1 = query_distmat[queryIndex][queryIndex2];
            float dist2 = target_distmat[targetIndex][targetIndex2];
            float distAverage = (dist1 + dist2) / 2.0;

            dali[i] += (DALI_THETA_E - fabs(dist1-dist2)/distAverage) * exp(-distAverage*distAverage/DALI_ALPHA);
        }
        daliScore += dali[queryIndex];
        if(dali[queryIndex] < 0) {
            dali[queryIndex] = DALI_EPS;
        }
    }

    float numAtomsEffective = sqrtf(queryLength * targetLength);
    float L = numAtomsEffective;
    if(L > DALI_LMAX) {
        L = DALI_LMAX;
    }
    float mean = DALI_C1 + DALI_C2*L + DALI_C3*L*L + DALI_C4*L*L*L;
    if(numAtomsEffective > DALI_LMAX) {
        mean += numAtomsEffective - L;
    }
    float sigma = 0.50 * mean;

    float daliZScore = (daliScore - mean) / sigma;
    return DaliScoreResult(daliZScore);
}