//
// Created by Charlotte Tumescheit on 2021/12/22.
//

#ifndef FOLDSEEK_NEARESTNEIGHBOURSW_H
#define FOLDSEEK_NEARESTNEIGHBOURSW_H

#include <string>
#include <sstream>
#include <tmalign/TMalign.h>
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "QueryMatcher.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "PareunAlign.h"
#include "simd.h"
#include "structureto3diseqdist.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include <cmath>

class NearestNeighbourSW {
public:

    NearestNeighbourSW(size_t maxSequenceLength, int alphabetSize);
    ~NearestNeighbourSW();





private:

    float * query_x;
    float * query_y;
    float * query_z;
    simd_int* query_profile_word;
    simd_int* query_profile_nn;
    simd_int* query_profile_dist ;

    float * target_x;
    float * target_y;
    float * target_z;

    char * querynn;
    char * targetnn;
    char * querynn_dist;
    char * targetnn_dist;

    simd_int* vHStore;
    simd_int* vHLoad;
    simd_int* vE;
    simd_int* vHmax;
    simd_int* nnScoreVec;
    uint8_t * maxColumn;
    std::pair<float, int> * seqDistList;
    float * seqDistSimd;
};


#endif //FOLDSEEK_NEARESTNEIGHBOURSW_H
