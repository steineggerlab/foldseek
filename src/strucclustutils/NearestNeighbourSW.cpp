//
// Created by Charlotte Tumescheit on 2021/12/22.
//

#include "NearestNeighbourSW.h"
// par.maxSeqLen,subMat.alphabetSize
NearestNeighbourSW::NearestNeighbourSW(size_t maxSequenceLength, int alphabetSize) {
    const int segmentSize = (maxSequenceLength+7)/8;
    query_x = (float*)mem_align(ALIGN_FLOAT, maxSequenceLength * sizeof(float) );
    query_y = (float*)mem_align(ALIGN_FLOAT, maxSequenceLength * sizeof(float) );
    query_z = (float*)mem_align(ALIGN_FLOAT, maxSequenceLength * sizeof(float) );
    query_profile_word = (simd_int*)mem_align(ALIGN_INT, alphabetSize * segmentSize * sizeof(simd_int));
    query_profile_nn  = (simd_int*)mem_align(ALIGN_INT, 6 * segmentSize * sizeof(simd_int));
    query_profile_dist  = (simd_int*)mem_align(ALIGN_INT, 6 * segmentSize * sizeof(simd_int));

    target_x = (float*)mem_align(ALIGN_FLOAT, maxSequenceLength * sizeof(float) );
    target_y = (float*)mem_align(ALIGN_FLOAT, maxSequenceLength * sizeof(float) );
    target_z = (float*)mem_align(ALIGN_FLOAT, maxSequenceLength * sizeof(float) );

    querynn  = (char*)mem_align(ALIGN_INT,maxSequenceLength * sizeof(char) * 8 ); // 8 to number of neighbours
    targetnn = (char*)mem_align(ALIGN_INT, maxSequenceLength * sizeof(char) * 8 );
    querynn_dist  = (char*)mem_align(ALIGN_INT, maxSequenceLength * sizeof(char) * 8 );
    targetnn_dist = (char*)mem_align(ALIGN_INT, maxSequenceLength * sizeof(char) * 8 );

    vHStore = (simd_int*) mem_align(ALIGN_INT, segmentSize * sizeof(simd_int));
    vHLoad  = (simd_int*) mem_align(ALIGN_INT, segmentSize * sizeof(simd_int));
    vE      = (simd_int*) mem_align(ALIGN_INT, segmentSize * sizeof(simd_int));
    vHmax   = (simd_int*) mem_align(ALIGN_INT, segmentSize * sizeof(simd_int));
    nnScoreVec   = (simd_int*) mem_align(ALIGN_INT, segmentSize * sizeof(simd_int));
    maxColumn = new uint8_t[maxSequenceLength*sizeof(uint16_t)];
    seqDistList = new std::pair<float, int>[maxSequenceLength];
    seqDistSimd = (float*)mem_align(ALIGN_FLOAT, maxSequenceLength * sizeof(float));
}

NearestNeighbourSW::~NearestNeighbourSW() {
    free(target_x);
    free(target_y);
    free(target_z);
    free(query_x);
    free(query_y);
    free(query_z);
    free(querynn);
    free(targetnn);
    free(query_profile_word);
    free(query_profile_nn);
    free(query_profile_dist);
    free(vHStore);
    free(vHLoad);
    free(vE);
    free(vHmax);
    free(seqDistSimd);
    free(nnScoreVec);
    free(querynn_dist);
    free(targetnn_dist);

    delete [] maxColumn;
    delete [] seqDistList;
}