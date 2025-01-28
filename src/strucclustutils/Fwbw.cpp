#include "Fwbw.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "SubstitutionMatrix.h"
#include "Matcher.h"
#include "Util.h"
#include "simd.h"
#include "Sequence.h"
#include "Timer.h"
#include "Fwbw.h"

#include <iostream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <numeric>
#include <cmath>
#include <vector>
#include <fstream>


#ifdef OPENMP
#include <omp.h>
#endif

struct States {
    const static uint8_t STOP=0;
    const static uint8_t M=1;
    const static uint8_t I=2;
    const static uint8_t D=3;
};

struct FWBWState {
    const static bool FORWARD = true;
    const static bool BACKWARD = false;
};

inline void calculate_max4(float& max, float& term1, float& term2, float& term3, float& term4, uint8_t& state) {
    max = term1;
    state = States::STOP;
    if (term2 > max) { max = term2; state = States::M; }
    if (term3 > max) { max = term3; state = States::I; }
    if (term4 > max) { max = term4; state = States::D; }
}

inline simd_float simdf32_prefixsum(simd_float a) {
    a = simdf32_add(a, simdi_i2fcast(simdi8_shiftl(simdf_f2icast(a), 4)));
    a = simdf32_add(a, simdi_i2fcast(simdi8_shiftl(simdf_f2icast(a), 8)));

#ifdef AVX2
    a = simdf32_add(a, simdi_i2fcast(simdi8_shiftl(simdf_f2icast(a), 16)));
#endif
    return a;
}

// FwBwAligner constructor for general use case
FwBwAligner::FwBwAligner(size_t length, SubstitutionMatrix &subMat, float gapOpen, float gapExtend, float mact, float temperature, size_t rowsCapacity, size_t colsCapacity)
                : length(length), gapOpen(gapOpen), gapExtend(gapExtend), mact(mact), temperature(temperature), rowsCapacity(rowsCapacity), colsCapacity(colsCapacity) {    
    // ZM
    zm = malloc_matrix<float>(rowsCapacity, colsCapacity);
    // Block
    zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zmBlockPrev = (float *) malloc_simd_float((length+1) * sizeof(float));
    zmBlockCurr = (float *) malloc_simd_float((length+1) * sizeof(float));
    zeBlock = (float *) malloc_simd_float((length+1) * sizeof(float));
    zfBlock = (float *) malloc_simd_float((length+1) * sizeof(float));
    // zInit forward & backward
    zInit = malloc_matrix<float>(3, rowsCapacity);
    // Score Matrix (targetProfile 21xcolsCapacity)
    scoreForwardProfile = malloc_matrix<float>(21, colsCapacity);
    scoreForwardProfile_exp = malloc_matrix<float>(21, colsCapacity);
    scoreBackwardProfile_exp = malloc_matrix<float>(21, colsCapacity);
    btMatrix = malloc_matrix<uint8_t>(rowsCapacity + 1, colsCapacity + 1);

    // btMatrix 
    S_prev = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    S_curr = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    
    // V,J,exp_ge_arr for ZE
    vj = (float *) malloc_simd_float(length * sizeof(float));
    wj = (float *) malloc_simd_float(length * sizeof(float));
    exp_ge_arr = (float *) malloc_simd_float(length * sizeof(float));

    for (size_t i = 0; i < length; ++i) { 
        vj[i] = exp(((length - 1) * gapExtend + gapOpen - i * gapExtend) / temperature);
        wj[i] = exp(((length - 1) * gapExtend - i * gapExtend) / temperature);
    }
    for (size_t i = 0; i < length; ++i) {
        exp_ge_arr[i] = exp((i * gapExtend + gapExtend) / temperature);
    }
    // Gap open and extend
    exp_go = simdf32_set(static_cast<float>(exp(gapOpen / temperature))); 
    exp_ge = simdf32_set(static_cast<float>(exp(gapExtend / temperature)));
    // Blosum matrix
    blosum = malloc_matrix<float>(21, 21);
    for (int i = 0; i < subMat.alphabetSize; ++i) {
        for (int j = 0; j < subMat.alphabetSize; ++j) {
            blosum[i][j] = static_cast<float>(subMat.subMatrix[i][j]);
        }
    }
    // queryAAs and targetAAs
    // queryNum = nullptr;
    // targetNum = nullptr;
}

// FwBwAligner for lolalign
FwBwAligner::FwBwAligner(size_t length, float gapOpen, float gapExtend, float temperature, size_t rowsCapacity, size_t colsCapacity)
                    : length(length), gapOpen(gapOpen), gapExtend(gapExtend), temperature(temperature), rowsCapacity(rowsCapacity), colsCapacity(colsCapacity) {
    
    // scoreForward
    scoreForward = malloc_matrix<float>(rowsCapacity, colsCapacity);
    // ZM
    zm = malloc_matrix<float>(rowsCapacity, colsCapacity);
    // Block
    zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zmBlockPrev = (float *) malloc_simd_float((length+1) * sizeof(float));
    zmBlockCurr = (float *) malloc_simd_float((length+1) * sizeof(float));
    zeBlock = (float *) malloc_simd_float((length+1) * sizeof(float));
    zfBlock = (float *) malloc_simd_float((length+1) * sizeof(float));

    
    // zInit forward & backward
    zInit = malloc_matrix<float>(3, rowsCapacity);
    
    // V,J,exp_ge_arr for ZE
    vj = (float *) malloc_simd_float(length * sizeof(float));
    wj = (float *) malloc_simd_float(length * sizeof(float));
    exp_ge_arr = (float *) malloc_simd_float(length * sizeof(float));

    for (size_t i = 0; i < length; ++i) { 
        vj[i] = exp(((length - 1) * gapExtend + gapOpen - i * gapExtend) / temperature);
        wj[i] = exp(((length - 1) * gapExtend - i * gapExtend) / temperature);
    }
    for (size_t i = 0; i < length; ++i) {
        exp_ge_arr[i] = exp((i * gapExtend + gapExtend) / temperature);
    }
    // Gap open and extend
    exp_go = simdf32_set(static_cast<float>(exp(gapOpen / temperature))); 
    exp_ge = simdf32_set(static_cast<float>(exp(gapExtend / temperature)));

}

FwBwAligner::~FwBwAligner(){
    // matrices used both in lolalign and general cases
    /*free(zm);
    free(zmBlockPrev);
    free(zmBlockCurr);
    free(zeBlock);
    free(zfBlock);
    free(zmFirst);
    free(zeFirst);
    free(zfFirst);
    free(vj);
    free(wj);
    free(zInit);
    free(exp_ge_arr);


    // matrices used in only one case
    if (scoreForward != nullptr) {
        free(scoreForward);
    }
    if (scoreForwardProfile != nullptr) {
        free(scoreForwardProfile);
    }
    if (scoreForwardProfile_exp != nullptr) {
        free(scoreForwardProfile_exp);
    }
    if (scoreBackwardProfile_exp != nullptr) {
        free(scoreBackwardProfile_exp);
    }
    if (btMatrix != nullptr) {
        free(btMatrix);
    }
    if (blosum != nullptr) {
        free(blosum);
    }
    if (S_prev != nullptr) {
        free(S_prev);
    }
    if (S_curr != nullptr) {
        free(S_curr);
    }

    // free(btMatrix);
    // free(blosum);
    // free(S_prev);
    // free(S_curr);
    std::cout << "lol Freeing memory" << std::endl;*/
}


void FwBwAligner::reallocateScoreProfile(size_t newColsCapacity) {
    free(scoreForwardProfile); scoreForwardProfile = malloc_matrix<float>(21, newColsCapacity);
    free(scoreForwardProfile_exp); scoreForwardProfile_exp = malloc_matrix<float>(21, newColsCapacity);
    free(scoreBackwardProfile_exp); scoreBackwardProfile_exp = malloc_matrix<float>(21, newColsCapacity);
}

// lolalign
void FwBwAligner::initScoreMatrix(float** inputScoreMatrix, size_t queryLen, size_t targetLen, int * gaps) {
    tlen = targetLen; qlen = queryLen;
    max_zm = -std::numeric_limits<float>::max(); vMax_zm = simdf32_set(max_zm);
    qlen_padding = ((queryLen + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT) * VECSIZE_FLOAT;
    blocks = (queryLen / length) + (queryLen % length != 0);
    simd_float vTemp = simdf32_set(temperature);
    if (queryLen > colsCapacity) {
        size_t newColsCapacity = ((queryLen + length-1)/length)* length;
        free(scoreForward); scoreForward = malloc_matrix<float>(targetLen, newColsCapacity);
    }
    
    /*for (size_t i=0; i<targetLen; ++i){
        for (size_t j = 0; j < qlen_padding; j += VECSIZE_FLOAT) {
            simd_float vScoreForward = simdf32_loadu(&inputScoreMatrix[i][j]);
            vScoreForward = simdf32_div(vScoreForward, vTemp);
            simdf32_store(&scoreForward[i][j], vScoreForward);
        }
        std::fill(&scoreForward[i][queryLen], &scoreForward[i][qlen_padding], FLT_MIN_EXP);
    }*/
    for (size_t i=0; i<targetLen; ++i){
        for (size_t j = 0; j < queryLen; ++j) {
            float score = inputScoreMatrix[i+gaps[0]][j+gaps[2]]/temperature;
            scoreForward[i][j] = score;
        }
        for(size_t j = colsCapacity; j < qlen_padding; ++j){
            scoreForward[i][j] = FLT_MIN_EXP;
        }
        //std::fill(&scoreForward[i][queryLen], &scoreForward[i][qlen_padding], FLT_MIN_EXP);
    }
    
}

void FwBwAligner::resizeMatrix(size_t newRowsCapacity, size_t newColsCapacity) {
    rowsCapacity = newRowsCapacity;
    colsCapacity = newColsCapacity;
    // blockCapacity = newColsCapacity / length;
    free(scoreForward); scoreForward = malloc_matrix<float>(rowsCapacity, colsCapacity);

    free(zm); zm = malloc_matrix<float>(rowsCapacity, colsCapacity);    
    free(zmFirst); zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zeFirst); zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zfFirst); zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));

    free(zInit); zInit = malloc_matrix<float>(3, rowsCapacity);
    free(S_prev); S_prev = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    free(S_curr); S_curr = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    if (btMatrix != nullptr) {
        free(btMatrix); btMatrix = malloc_matrix<uint8_t>(rowsCapacity + 1, colsCapacity + 1);
    }
    // free(btMatrix); btMatrix = malloc_matrix<uint8_t>(rowsCapacity + 1, colsCapacity + 1);
}

void FwBwAligner::initAlignment(unsigned char* targetAANum, size_t targetLen) {
    targetNum = targetAANum; tlen = targetLen;
    max_zm = -std::numeric_limits<float>::max(); vMax_zm = simdf32_set(max_zm);
    memset(S_prev, 0, (qlen + 1) * sizeof(float));
    memset(S_curr, 0, (qlen + 1) * sizeof(float));

}

void FwBwAligner::initQueryProfile(unsigned char* queryAANum, size_t queryLen) {
    queryNum = queryAANum; qlen = queryLen;
    qlen_padding = ((queryLen + VECSIZE_FLOAT -1) / VECSIZE_FLOAT) * VECSIZE_FLOAT;
    blocks = (queryLen / length) + (queryLen % length != 0);
    if (queryLen > colsCapacity) {
        size_t newColsCapacity = ((queryLen + length-1)/length)* length;
        reallocateScoreProfile(newColsCapacity);
    }
    
    //scoreForwardProfile : 21 * qlen
    for (size_t i=0; i<21; ++i){
        for (size_t j=0; j < queryLen; ++j) {
            float score = blosum[i][queryAANum[j]]/temperature;
            scoreForwardProfile[i][j] = score;
        }
        std::fill(&scoreForwardProfile[i][queryLen], &scoreForwardProfile[i][qlen_padding], FLT_MIN_EXP);
    }
    for (size_t i=0; i<21; ++i){
        for (size_t j=0; j < qlen_padding; j += VECSIZE_FLOAT) {
            simd_float vScoreForward = simdf32_load(&scoreForwardProfile[i][j]);
            vScoreForward = simdf32_exp(vScoreForward);
            simdf32_store(&scoreForwardProfile_exp[i][j], vScoreForward);
        }
        for (size_t j=0; j < queryLen; ++j){
            size_t reverse_j = queryLen - 1 - j;
            scoreBackwardProfile_exp[i][reverse_j] = scoreForwardProfile_exp[i][j];
        }
        //remainder 
        for (size_t j=queryLen; j < qlen_padding; ++j){
            scoreBackwardProfile_exp[i][j] = 0;
        }
    }
}

void FwBwAligner::forward(bool has_Profile) {   // if has_Profile = false, then it is lolalign
    //Init zInit
    for (size_t i = 0 ; i < 3; ++i) {
        std::fill(zInit[i], zInit[i] + rowsCapacity, FLT_MIN_EXP); // rowsCapacity -> tlen
    }  
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, qlen) - start;
        size_t cols = length;
        if (memcpy_cols != length) {
            cols = ((memcpy_cols + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT) * VECSIZE_FLOAT; //padding vecsize_float
        }
        //Init blocks
        memset(zmBlockPrev, 0, (length + 1) * sizeof(float));
        memset(zeBlock, 0, (length + 1) * sizeof(float));
        memset(zfBlock, 0, (length + 1) * sizeof(float));
            
        memcpy(zmFirst + 1, zInit[0], tlen * sizeof(float));
        memcpy(zeFirst + 1, zInit[1], tlen * sizeof(float));
        memcpy(zfFirst + 1, zInit[2], tlen * sizeof(float));

        //Init initial values
        zmBlockPrev[0] = 0; zfBlock[0] = 0; zeBlock[0] = 0;
        zmBlockCurr[0] = exp(zmFirst[1]);
        float ze_i0 = expf(zeFirst[1]);
        float current_max = 0;
        float zmMaxRowBlock = -std::numeric_limits<float>::max();
        float* zmMaxData = (float*)malloc_simd_float(VECSIZE_FLOAT * sizeof(float));
        float log_zmMax = 0;
        simd_float vZmMaxRowBlock;
        for (size_t i = 1; i <= tlen; ++i) {
            simd_float vExpMax = simdf32_set(exp(-current_max));
            simd_float vZeI0 = simdf32_set(ze_i0);
            simd_float vLastPrefixSum = simdf32_setzero(); 
            simd_float vZmax_tmp = simdf32_set(-std::numeric_limits<float>::max());
            // ZM calculation
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_load(&zmBlockPrev[j-1]);
                simd_float vZe = simdf32_load(&zeBlock[j-1]);
                simd_float vZf = simdf32_load(&zfBlock[j-1]);
                simd_float vScoreMatrix;
                if (has_Profile) {
                    vScoreMatrix = simdf32_exp(simdf32_load(&scoreForwardProfile[targetNum[i-1]][start + j - 1]));
                } else {
                    vScoreMatrix = simdf32_exp(simdf32_load(&scoreForward[i-1][start + j - 1]));
                }              
                simd_float vZmCurrUpdate = simdf32_add(simdf32_add(vZmPrev, vZe), simdf32_add(vZf, vExpMax));
                vZmCurrUpdate = simdf32_mul(vZmCurrUpdate, vScoreMatrix);
                vZmax_tmp = simdf32_max(vZmax_tmp, vZmCurrUpdate);
                simdf32_storeu(&zmBlockCurr[j], vZmCurrUpdate);
            }
            simdf32_store(zmMaxData, vZmax_tmp);
            zmMaxRowBlock = *std::max_element(zmMaxData, zmMaxData+VECSIZE_FLOAT); // or hmax
            vZmMaxRowBlock = simdf32_set(zmMaxRowBlock);
            // ZF calculation 
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_loadu(&zmBlockPrev[j]);
                simd_float vZf = simdf32_loadu(&zfBlock[j]);
                simd_float vZfUpdate = simdf32_add(
                                        simdf32_mul(vZmPrev, exp_go),
                                        simdf32_mul(vZf, exp_ge)
                                        );
                vZfUpdate = simdf32_div(vZfUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zfBlock[j], vZfUpdate);
            }
            for (size_t j = 0; j < cols; j += VECSIZE_FLOAT) { 
                simd_float vZmCurr = simdf32_load(&zmBlockCurr[j]);
                simd_float vVj = simdf32_load(&vj[j]);
                simd_float vCumsumZm = simdf32_mul(vZmCurr, vVj);
                vCumsumZm = simdf32_prefixsum(vCumsumZm);
                vCumsumZm = simdf32_add(vCumsumZm, vLastPrefixSum);
                vLastPrefixSum = simdf32_set(vCumsumZm[(VECSIZE_FLOAT - 1)]);
                simd_float vWj = simdf32_load(&wj[j]);
                simd_float vExp_ge_arr = simdf32_load(&exp_ge_arr[j]);
                simd_float vZeUpdate = simdf32_add(
                                        simdf32_div(vCumsumZm, vWj),
                                        simdf32_mul(vZeI0, vExp_ge_arr)
                                        );
                vZeUpdate = simdf32_div(vZeUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zeBlock[j+1], vZeUpdate);
            }

            log_zmMax = log(zmMaxRowBlock);
            current_max += log_zmMax;
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT){
                simd_float vZmCurr = simdf32_loadu(&zmBlockCurr[j]);
                vZmCurr = simdf32_div(vZmCurr, vZmMaxRowBlock);
                simdf32_storeu(&zmBlockCurr[j], vZmCurr);
                vZmCurr = simdf32_add(simdf32_log(vZmCurr), simdf32_set(current_max));
                vMax_zm = simdf32_max(vMax_zm, vZmCurr);
                simdf32_store(&zm[i - 1][start + j-1], vZmCurr);
            }     

            #if defined(AVX512)
                simd_float vNextZinit = _mm512_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
            #elif defined(AVX2)
                simd_float vNextZinit = _mm256_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
            #else // Fallback to SSE
                simd_float vNextZinit = _mm_set_ps(
                    1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
            #endif


            vNextZinit = simdf32_log(vNextZinit);
            vNextZinit = simdf32_add(vNextZinit, simdf32_set(current_max));

            zInit[0][i-1] = vNextZinit[0]; zInit[1][i-1] = vNextZinit[1]; zInit[2][i-1] = vNextZinit[2];
            std::swap(zmBlockCurr, zmBlockPrev);
            //original code
            // zInit[0][i-1] = log(zmBlockCurr[memcpy_cols]) + current_max;
            // zInit[1][i-1] = log(zeBlock[memcpy_cols]) + current_max;
            // zInit[2][i-1] = log(zfBlock[memcpy_cols]) + current_max;
            // zmFirst[i] = exp(zmFirst[i] - log_zmMax); zmBlockPrev[0] = zmFirst[i];
            // zeFirst[i] = exp(zeFirst[i] - log_zmMax); zeBlock[0] = zeFirst[i];
            // zfFirst[i] = exp(zfFirst[i] - current_max); zfBlock[0] = zfFirst[i];
            
            if (i < tlen) {
                zmFirst[i+1] -= current_max;
                zeFirst[i+1] -= current_max;

#if defined(AVX512)
                simd_float vNextFirstExp = _mm512_set_ps(
                    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max, 
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0]; ze_i0 = vNextFirstExp[1];
                zmBlockPrev[0] = vNextFirstExp[2]; zeBlock[0] = vNextFirstExp[3]; zfBlock[0] = vNextFirstExp[4];
#elif defined(AVX2)
                simd_float vNextFirstExp = _mm256_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max, 
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );        
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0]; ze_i0 = vNextFirstExp[1]; 
                zmBlockPrev[0] = vNextFirstExp[2]; zeBlock[0] = vNextFirstExp[3]; zfBlock[0] = vNextFirstExp[4];
#else // Fallback to SSE
                simd_float vNextFirstExp1 = _mm_set_ps(
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );    
                simd_float vNextFirstExp2= _mm_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max    
                );    
                vNextFirstExp1 = simdf32_exp(vNextFirstExp1);
                vNextFirstExp2 = simdf32_exp(vNextFirstExp2);
                zmBlockCurr[0] = vNextFirstExp1[0]; ze_i0 = vNextFirstExp1[1]; 
                zmBlockPrev[0] = vNextFirstExp1[2]; zeBlock[0] = vNextFirstExp1[3];
                zfBlock[0] = vNextFirstExp2[0];
#endif 
            } else{
                zmBlockPrev[0] = exp(zmFirst[i] - log_zmMax);
                zeBlock[0] = exp(zeFirst[i] - log_zmMax); 
                zfBlock[0] = exp(zfFirst[i] - current_max); 
            }
        }  
        free(zmMaxData);  
    }
}

void FwBwAligner::backward(bool has_Profile) {
    //Init zInit
    // Debug(Debug::INFO) << "backward function called with has_Profile: " << has_Profile << '\n';
    size_t vecsize_float = static_cast<size_t>(VECSIZE_FLOAT);
    for (size_t i = 0 ; i < 3; ++i) {
        std::fill(zInit[i], zInit[i] + rowsCapacity, FLT_MIN_EXP); // rowsCapacity -> tlen
    }  
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, qlen) - start;
        size_t cols = length;
        if (memcpy_cols != length) {
            cols = ((memcpy_cols + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT) * VECSIZE_FLOAT; //padding vecsize_float
        }
        //Init blocks
        memset(zmBlockPrev, 0, (length + 1) * sizeof(float));
        memset(zeBlock, 0, (length + 1) * sizeof(float));
        memset(zfBlock, 0, (length + 1) * sizeof(float));

        memcpy(zmFirst + 1, zInit[0], tlen * sizeof(float));
        memcpy(zeFirst + 1, zInit[1], tlen * sizeof(float));
        memcpy(zfFirst + 1, zInit[2], tlen * sizeof(float));

        //Init initial values
        zmBlockPrev[0] = 0; zfBlock[0] = 0; zeBlock[0] = 0;
        zmBlockCurr[0] = exp(zmFirst[1]);
        float ze_i0 = expf(zeFirst[1]);
        float current_max = 0;
        float zmMaxRowBlock = -std::numeric_limits<float>::max();
        float* zmMaxData = (float*)malloc_simd_float(VECSIZE_FLOAT * sizeof(float));
        float log_zmMax = 0;
        simd_float vZmMaxRowBlock;
        for (size_t i = 1; i <= tlen; ++i) {
            simd_float vExpMax = simdf32_set(exp(-current_max));
            simd_float vZeI0 = simdf32_set(ze_i0);
            simd_float vLastPrefixSum = simdf32_setzero(); 
            simd_float vZmax_tmp = simdf32_set(-std::numeric_limits<float>::max());
            
            // ZM calculation
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_load(&zmBlockPrev[j-1]);
                simd_float vZe = simdf32_load(&zeBlock[j-1]);
                simd_float vZf = simdf32_load(&zfBlock[j-1]);
                simd_float vScoreMatrix;
                if (has_Profile) {
                    vScoreMatrix = simdf32_load(&scoreBackwardProfile_exp[targetNum[tlen - i]][start + j - 1]);
                } else {
                    size_t reverse_i = tlen - i;
                    size_t reverse_j = qlen - start - j + 1;
                    simd_float vScoreBackward;
                    if (reverse_j >= vecsize_float) {
                        vScoreBackward = simdf32_loadu(&scoreForward[reverse_i][reverse_j - vecsize_float]);
                    } else {
                        size_t elements_to_fill = reverse_j;
                        vScoreBackward = simdf32_set(FLT_MIN_EXP);
                        // fill from the back
                        for (size_t k = 0; k < elements_to_fill; ++k) {
                            vScoreBackward[VECSIZE_FLOAT - elements_to_fill + k] = scoreForward[reverse_i][k];
                        }
                    }
                    vScoreBackward = simdf32_reverse(vScoreBackward);
                    vScoreMatrix = simdf32_exp(vScoreBackward);
                }
                simd_float vZmCurrUpdate = simdf32_add(simdf32_add(vZmPrev, vZe), simdf32_add(vZf, vExpMax));
                vZmCurrUpdate = simdf32_mul(vZmCurrUpdate, vScoreMatrix);
                vZmax_tmp = simdf32_max(vZmax_tmp, vZmCurrUpdate);
                simdf32_storeu(&zmBlockCurr[j], vZmCurrUpdate);
            }
            simdf32_store(zmMaxData, vZmax_tmp);
            zmMaxRowBlock = *std::max_element(zmMaxData, zmMaxData+VECSIZE_FLOAT); // or hmax
            vZmMaxRowBlock = simdf32_set(zmMaxRowBlock);
            // ZF calculation 
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_loadu(&zmBlockPrev[j]);
                simd_float vZf = simdf32_loadu(&zfBlock[j]);
                simd_float vZfUpdate = simdf32_add(
                                        simdf32_mul(vZmPrev, exp_go),
                                        simdf32_mul(vZf, exp_ge)
                                        );
                vZfUpdate = simdf32_div(vZfUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zfBlock[j], vZfUpdate);
            }
            for (size_t j = 0; j < cols; j += VECSIZE_FLOAT) { 
                simd_float vZmCurr = simdf32_load(&zmBlockCurr[j]);
                simd_float vVj = simdf32_load(&vj[j]);
                simd_float vCumsumZm = simdf32_mul(vZmCurr, vVj);
                vCumsumZm = simdf32_prefixsum(vCumsumZm);
                vCumsumZm = simdf32_add(vCumsumZm, vLastPrefixSum);
                vLastPrefixSum = simdf32_set(vCumsumZm[(VECSIZE_FLOAT - 1)]);
                simd_float vWj = simdf32_load(&wj[j]);
                simd_float vExp_ge_arr = simdf32_load(&exp_ge_arr[j]);
                simd_float vZeUpdate = simdf32_add(
                                        simdf32_div(vCumsumZm, vWj),
                                        simdf32_mul(vZeI0, vExp_ge_arr)
                                        );
                vZeUpdate = simdf32_div(vZeUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zeBlock[j+1], vZeUpdate);
            }

            log_zmMax = log(zmMaxRowBlock);
            current_max += log_zmMax;
            size_t adjusted_memcpycols = memcpy_cols - memcpy_cols % VECSIZE_FLOAT;
            size_t forwardBlockStart = qlen - start;
            for (size_t j = 1; j <= adjusted_memcpycols; j += VECSIZE_FLOAT) {
                forwardBlockStart -= vecsize_float;
                size_t simd_index = forwardBlockStart;
                simd_float vZmCurr = simdf32_loadu(&zmBlockCurr[j]);
                vZmCurr = simdf32_div(vZmCurr, vZmMaxRowBlock);
                simdf32_storeu(&zmBlockCurr[j], vZmCurr);
                vZmCurr = simdf32_add(simdf32_log(vZmCurr), simdf32_set(current_max));

                simd_float vZmForward = simdf32_loadu(&zm[tlen - i][simd_index]);

                simd_float vZmCurr_reverse = simdf32_reverse(vZmCurr);
                simd_float vZmForward_Backward = simdf32_add(vZmForward, vZmCurr_reverse);
                simdf32_storeu(&zm[tlen - i][simd_index], vZmForward_Backward);
            }

            // Handle remainder
            if (memcpy_cols != length) {
                size_t remainder = memcpy_cols % VECSIZE_FLOAT;
                simd_float vZmCurr = simdf32_loadu(&zmBlockCurr[adjusted_memcpycols+1]);
                vZmCurr = simdf32_div(vZmCurr, vZmMaxRowBlock);
                simdf32_storeu(&zmBlockCurr[adjusted_memcpycols+1], vZmCurr);
                vZmCurr = simdf32_add(simdf32_log(vZmCurr), simdf32_set(current_max));
                for (size_t k = 0; k < remainder; ++k) {
                    size_t j_index = remainder - k - 1;
                    zm[tlen - i][j_index] += vZmCurr[k];
                }
            }

#if defined(AVX512)
                simd_float vNextZinit = _mm512_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
#elif defined(AVX2)
                simd_float vNextZinit = _mm256_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
#else // Fallback to SSE
                simd_float vNextZinit = _mm_set_ps(
                    1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
#endif

            vNextZinit = simdf32_log(vNextZinit);
            vNextZinit = simdf32_add(vNextZinit, simdf32_set(current_max));

            zInit[0][i-1] = vNextZinit[0]; zInit[1][i-1] = vNextZinit[1]; zInit[2][i-1] = vNextZinit[2];
            std::swap(zmBlockCurr, zmBlockPrev);
            
            if (i < tlen) {
                zmFirst[i+1] -= current_max;
                zeFirst[i+1] -= current_max;

#if defined(AVX512)
                simd_float vNextFirstExp = _mm512_set_ps(
                    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max, 
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0]; ze_i0 = vNextFirstExp[1];
                zmBlockPrev[0] = vNextFirstExp[2]; zeBlock[0] = vNextFirstExp[3]; zfBlock[0] = vNextFirstExp[4];
#elif defined(AVX2)
                simd_float vNextFirstExp = _mm256_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max, 
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );        
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0]; ze_i0 = vNextFirstExp[1]; 
                zmBlockPrev[0] = vNextFirstExp[2]; zeBlock[0] = vNextFirstExp[3]; zfBlock[0] = vNextFirstExp[4];
#else // Fallback to SSE
                simd_float vNextFirstExp1 = _mm_set_ps(
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );    
                simd_float vNextFirstExp2= _mm_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max     
                );    
                vNextFirstExp1 = simdf32_exp(vNextFirstExp1);
                vNextFirstExp2 = simdf32_exp(vNextFirstExp2);
                zmBlockCurr[0] = vNextFirstExp1[0]; ze_i0 = vNextFirstExp1[1]; 
                zmBlockPrev[0] = vNextFirstExp1[2]; zeBlock[0] = vNextFirstExp1[3];
                zfBlock[0] = vNextFirstExp2[0];
#endif 
            } else{
                zmBlockPrev[0] = exp(zmFirst[i] - log_zmMax);
                zeBlock[0] = exp(zeFirst[i] - log_zmMax); 
                zfBlock[0] = exp(zfFirst[i] - current_max); 
            }
        }  
        free(zmMaxData); 
    }
}
void FwBwAligner::computeProbabilityMatrix(bool has_Profile) {
    //// Run Forward
    forward(has_Profile);
    //Calculate max_zm
    for (size_t i = 0; i < VECSIZE_FLOAT; ++i) {
        max_zm = std::max(max_zm, vMax_zm[i]);
    }
    vMax_zm = simdf32_set(max_zm);

    //Calculate sum_exp
    float sum_exp= 0.0; simd_float vSum_exp = simdf32_setzero();
    for (size_t i = 0; i < tlen; ++i) {
        //Removed remainder handling. Check needed
        for (size_t j = 0; j < qlen_padding; j+= VECSIZE_FLOAT) {
            simd_float vZmForward = simdf32_load(&zm[i][j]);
            vZmForward = simdf32_exp(simdf32_sub(vZmForward, vMax_zm));
            vSum_exp = simdf32_add(vSum_exp, vZmForward);
        }
    }
    sum_exp += simdf32_hadd(vSum_exp);

    //// Backward
    backward(has_Profile);

    float logsumexp_zm = max_zm + log(sum_exp); simd_float vLogsumexp_zm = simdf32_set(logsumexp_zm);
    size_t qLoopCount = qlen / VECSIZE_FLOAT; size_t qLoopEndPos = qLoopCount * VECSIZE_FLOAT;
    P = zm;
    simd_float vMaxP = simdf32_setzero();
    for (size_t i = 0; i < tlen; ++i) {
        // Fill the probability matrix
        //Removed remainder handling. Check needed
        for (size_t j = 0; j < qLoopEndPos; j += VECSIZE_FLOAT) {
            simd_float vZmForward_Backward = simdf32_load(&zm[i][j]);
            simd_float scoreForwardVal;
            if (has_Profile) {
                scoreForwardVal = simdf32_load(&scoreForwardProfile[targetNum[i]][j]);
            } else {
                scoreForwardVal = simdf32_load(&scoreForward[i][j]);
            }
            // simd_float scoreForwardVal = simdf32_load(&scoreForward[targetNum[i]][j]);
            simd_float P_val = simdf32_exp(simdf32_sub(vZmForward_Backward, simdf32_add(scoreForwardVal, vLogsumexp_zm)));
            simdf32_store(&P[i][j], P_val);
            vMaxP = simdf32_max(vMaxP, P_val);
        }    
        for (size_t j = qLoopEndPos; j < qlen; ++j) {
            if (has_Profile) {
                P[i][j] = exp(zm[i][j] - scoreForwardProfile[targetNum[i]][j] - logsumexp_zm);
            } else {
                P[i][j] = exp(zm[i][j] - scoreForward[i][j] - logsumexp_zm);
            }
            // P[i][j] = exp(zm[i][j] - scoreForward[targetNum[i]][j] - logsumexp_zm);
            maxP = std::max(maxP, P[i][j]);
        }
    }

    // Calculate the maximum probability
    for (size_t k = 0; k < VECSIZE_FLOAT; ++k) {
        maxP = std::max(maxP, vMaxP[k]);
    }
}

FwBwAligner::s_align FwBwAligner::computeAlignment() {
    bool has_Profile = true;
    computeProbabilityMatrix(has_Profile);

    // MAC algorithm from HH-suite
    size_t qLoopCount = qlen / VECSIZE_FLOAT; size_t qLoopEndPos = qLoopCount * VECSIZE_FLOAT;
    uint8_t val;
    size_t max_i = 0;
    size_t max_j = 0;
    float term1, term2, term3, term4 = 0.0f;
    float score_MAC = -std::numeric_limits<float>::max();


    simd_float vMact = simdf32_set(mact); 
    simd_float vhalfMact = simdf32_set(0.5f * mact); 
    simd_float vStateStop = simdf32_set(0.0f); 
    simd_float vStateM = simdf32_set(1.0f); 
    simd_float vStateD = simdf32_set(2.0f); 
    simd_float vTerm1, vTerm2, vTerm3, vTerm4;
    simd_float vMax123, vMax1_2;
    simd_float vStateMask_MS, vStateMask_MSI, vStateTmp, vState;
    for (size_t i = 0; i <= tlen; ++i) {
        btMatrix[i][0] = States::STOP;
    }
    for (size_t j = 0; j <= qlen; ++j) {
        btMatrix[0][j] = States::STOP;
    }

    for (size_t i = 1; i <= tlen; ++i) {
        for (size_t j=1; j <= qLoopEndPos; j+= VECSIZE_FLOAT){
            vTerm1 = simdf32_sub(simdf32_loadu(&P[i - 1][j - 1]), vMact);
            vTerm2 = simdf32_add(simdf32_loadu(&S_prev[j - 1]), vTerm1);
            vTerm3 = simdf32_sub(simdf32_loadu(&S_prev[j]), vhalfMact);
            vTerm4 = simdf32_loadu(&S_curr[j - 1]);
            vMax1_2 = simdf32_max(vTerm1, vTerm2);
            vMax123 = simdf32_max(vMax1_2, vTerm3);

            vStateMask_MS = simdf32_gt(vTerm2, vTerm1); // true if Term1 < Term2
            vStateMask_MSI = simdf32_gt(vTerm3, vMax1_2); // Max(Term1, Term2) < Term3
            vStateTmp = simdf32_blendv_ps(vStateStop, vStateM, vStateMask_MS); // Term1 < Term2: choose M
            vState = simdf32_blendv_ps(vStateTmp, vStateD, vStateMask_MSI);

            term4 = vTerm4[0] - 0.5 * mact;
            for (size_t k = 0; k < VECSIZE_FLOAT; ++k) {
                if (term4 > vMax123[k]) {
                    vMax123[k] = term4;
                    vState[k] = States::I;  
                }
                term4 = vMax123[k] - 0.5 * mact;
                btMatrix[i][j+k] = vState[k];
                if (vMax123[k] > score_MAC) {
                    max_i = i;
                    max_j = j + k;
                    score_MAC = vMax123[k];
                }
            }
            // Store the results
            simdf32_storeu(&S_curr[j], vMax123);
            // simdf32_storeu(&btMatrix[i][j], vState);

        }
        for (size_t j = qLoopEndPos+1; j <= qlen; ++j) {
            term1 = P[i - 1][j - 1] - mact;
            term2 = S_prev[j - 1] + P[i - 1][j - 1] - mact;
            term3 = S_prev[j] - 0.5 * mact;
            term4 = S_curr[j - 1] - 0.5 * mact;
            calculate_max4(S_curr[j], term1, term2, term3, term4, val);
            btMatrix[i][j] = val;

            if (S_curr[j] > score_MAC) {
                max_i = i;
                max_j = j;
                score_MAC = S_curr[j];
            }
        }
        std::swap(S_prev, S_curr);
    }

    // traceback 
    s_align result;
    result.cigar = "";
    result.cigar.reserve(qlen + tlen);
    result.qEndPos1 = max_j - 1;
    result.dbEndPos1 = max_i - 1;
    result.score1 = maxP;
    result.score2 = score_MAC;
    while (max_i > 0 && max_j > 0) {
        uint8_t state = btMatrix[max_i][max_j];
        if (state == States::M) {
            --max_i;
            --max_j;
            // Update start positions to the last match
            result.qStartPos1 = max_j;
            result.dbStartPos1 = max_i;
            result.cigar.push_back('M');
        } else if (state == States::I) {
            --max_i;
            result.cigar.push_back('I');
        } else if (state == States::D) {
            --max_j;
            result.cigar.push_back('D');
        } else {
            break;
        }
    }
    while (!result.cigar.empty() && result.cigar.back() != 'M') {
        result.cigar.pop_back();
    }
    result.cigarLen = result.cigar.length();
    std::reverse(result.cigar.begin(), result.cigar.end());
    return result;
}

float** fwbw(float** inputScoreForward, size_t queryLen, size_t targetLen, float gapOpen, float gapExtend, float temperature) {
    size_t length = 16; // AVX2, need to check
    size_t assignSeqLen = VECSIZE_FLOAT * sizeof(float) * 20; 
    float** return_P = malloc_matrix<float>(targetLen, queryLen);
    FwBwAligner subfwbwaligner(length, gapOpen, gapExtend, temperature, assignSeqLen, assignSeqLen);
    //Resizing
    int* gaps = new int[4]{0, 0, 0, 0};
    subfwbwaligner.initScoreMatrix(inputScoreForward, queryLen, targetLen, gaps);
    //Resizing
    if (targetLen > subfwbwaligner.getRowsCapacity() && queryLen > subfwbwaligner.getColsCapacity()) {
        size_t newRowsCapacity = ((targetLen + subfwbwaligner.getBlockLength()-1)/subfwbwaligner.getBlockLength())* subfwbwaligner.getBlockLength();
        size_t newColsCapacity = ((queryLen + subfwbwaligner.getBlockLength()-1)/subfwbwaligner.getBlockLength())* subfwbwaligner.getBlockLength();
        subfwbwaligner.resizeMatrix(newRowsCapacity, newColsCapacity);
    }else if (targetLen > subfwbwaligner.getRowsCapacity()) {
        size_t newRowsCapacity = ((targetLen + subfwbwaligner.getBlockLength()-1)/subfwbwaligner.getBlockLength())* subfwbwaligner.getBlockLength();
        subfwbwaligner.resizeMatrix(newRowsCapacity, subfwbwaligner.getColsCapacity());
    }else if (queryLen > subfwbwaligner.getColsCapacity()) {
        size_t newColsCapacity = ((queryLen + subfwbwaligner.getBlockLength()-1)/subfwbwaligner.getBlockLength())* subfwbwaligner.getBlockLength();
        subfwbwaligner.resizeMatrix(subfwbwaligner.getRowsCapacity(), newColsCapacity);
    }
    
    bool has_Profile = false;
    subfwbwaligner.computeProbabilityMatrix(has_Profile);
    float** P = subfwbwaligner.getZm();
    // copy P to return_P with memcpy
    for (size_t i = 0; i < targetLen; ++i) {
        memcpy(return_P[i], P[i], queryLen * sizeof(float));
    }

    return return_P;  
}