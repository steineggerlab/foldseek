#include "LoLAlign.h"
#include "Fwbw.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "SubstitutionMatrix.h"
#include "Matcher.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Coordinate16.h"

#include <iostream>
#include <algorithm>
#include <limits>
#include <new>
#include <cstdlib>
#include <cfloat>
#include <numeric>
#include <cmath>
#include <vector>
#include <fstream>


#ifdef OPENMP
#include <omp.h>
#endif

LoLAlign::LoLAlign(unsigned int maxSeqLen, bool computeExactScore)
        : xtm(maxSeqLen), ytm(maxSeqLen), xt(maxSeqLen),
          r1(maxSeqLen), r2(maxSeqLen), computeExactScore(computeExactScore)
{
    queryX = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    queryY = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    queryZ = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    targetX = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    targetY = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    targetZ = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    distanceij = malloc_matrix<float>(maxSeqLen, maxSeqLen);
    distancekl = malloc_matrix<float>(maxSeqLen, maxSeqLen);
    lolScoreMatrix = malloc_matrix<float>(maxSeqLen, maxSeqLen);
    probabilityMatrix = malloc_matrix<float>(maxSeqLen, maxSeqLen);
    hiddenLayer = malloc_matrix<float>(maxSeqLen, 3);
    anchorQuery = malloc_matrix<int>(std::max(numStartAnchors, 2*seedNumber), maxSeqLen);
    anchorTarget = malloc_matrix<int>(std::max(numStartAnchors, 2*seedNumber), maxSeqLen);
    lolDist = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    lolSeqDist = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    lolScoreVec = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    lolScoreVecSh = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    
    startAnchorIndex = new int[numStartAnchors];
    startAnchorScores = new float[numStartAnchors];
    anchorLength = new int[numStartAnchors];
    newAnchorLength = new int[numStartAnchors];
    gaps = new int[4]{0, 0, 0, 0};

    finalAnchorQuery = nullptr;
    finalAnchorTarget = nullptr;
    //P = malloc_matrix<float>(maxSeqLen, maxSeqLen);
}



LoLAlign::~LoLAlign()
{
    free(queryX);
    free(queryY);
    free(queryZ);
    free(targetX);
    free(targetY);
    free(targetZ);
    free(distanceij);
    free(distancekl);
    free(probabilityMatrix);
    free(lolScoreMatrix);
    free(anchorQuery);
    free(anchorTarget);
    free(hiddenLayer);
    delete[] startAnchorIndex;
    delete[] startAnchorScores;
    delete[] anchorLength;
    delete[] newAnchorLength;
    delete[] gaps;
    free(lolDist);
    free(lolSeqDist);
    free(lolScoreVec);
    if (finalAnchorQuery != nullptr) {
        delete[] finalAnchorQuery;
    }
    if (finalAnchorTarget != nullptr) {
        delete[] finalAnchorTarget;
    }
    free(lolScoreVecSh);
}

void LoLAlign::calcGap(int* anchorQuery, int* anchorTarget, int* gaps, int queryLen, int targetLen)
{
    gaps[0] = -1;
    int indexQ = gaps[1];
    int indexT = gaps[3];
    while (anchorQuery[indexQ] != 0 || anchorTarget[indexT] != 0) {
        if (anchorQuery[indexQ] != 0 && anchorTarget[indexT] != 0) {
            indexQ++;
            indexT++;
        } else if (anchorQuery[indexQ] == 0) {
            indexQ++;
        } else if (anchorTarget[indexT] == 0) {
            indexT++;
        }
        if (indexQ == queryLen || indexT == targetLen) {
            gaps[0] = -1;
            return;
        }
    }
    gaps[0] = indexQ;
    gaps[2] = indexT;
    while (anchorQuery[indexQ] == 0 || anchorTarget[indexT] == 0) {
        if (anchorQuery[indexQ] == 0 && anchorTarget[indexT] == 0) {
            indexQ++;
            indexT++;
        } else if (anchorQuery[indexQ] == 0) {
            indexQ++;
        } else if (anchorTarget[indexT] == 0) {
            indexT++;
        }
        if (indexQ == queryLen || indexT == targetLen) {
            indexQ = queryLen;
            indexT = targetLen;
            break;
        }
    }
    gaps[1] = indexQ;
    gaps[3] = indexT;
    return;

}
float LoLAlign::maxSubArray(float* nums, int numsSize) {

    float currentMax = nums[0];
    float globalMax = nums[0];

    for (int i = 1; i < numsSize; ++i) {
        currentMax = std::max(nums[i], currentMax + nums[i]);
        globalMax = std::max(globalMax, currentMax);
    }

    return globalMax;
}

void LoLAlign::indexSort(float* nums, int* index, int numsSize) {
    std::sort(index, index + numsSize, [&nums](int i1, int i2) { return nums[i1] < nums[i2]; });
}




Matcher::result_t LoLAlign::align(unsigned int dbKey, float* targetX, float* targetY, float* targetZ,
                                  Sequence& tSeqAA, Sequence& tSeq3Di, int targetLen, SubstitutionMatrix& subMatAA,
                                  FwBwAligner* fwbwaln, int multiDomain)
{
    maxLolMatrixIndex = 0;
    minLolMatrixIndex = queryLen;


    unsigned char *targetNumAA = tSeqAA.numSequence;
    unsigned char *targetNum3Di = tSeq3Di.numSequence;
    const char *targetSeq = tSeqAA.getSeqData();
    computeForwardScoreMatrix(
            targetNumAA,
            targetNum3Di,
            queryLen,
            targetLen,
            subMatAA,
            lolScoreMatrix);

    int maxIndexX = 0;
    int maxIndexY = 0;

    
    for (int i = 0; i < numStartAnchors; i++) {
        startAnchorIndex[i] = i;
        startAnchorScores[i] = 0;
    }

    calcDistMatrix(targetX, targetY, targetZ, targetLen, distancekl, false);

    for (int startAnchor = 0; startAnchor < numStartAnchors; startAnchor++) {
        anchorLength[startAnchor] = 0;
        newAnchorLength[startAnchor] = 0;
        for(int i = 0; i <= maxAnchorLen; i++){
            anchorQuery[startAnchor][i] = 0;
            anchorTarget[startAnchor][i] = 0;

        }

    }


    gaps[0] = 0;
    gaps[1] = queryLen;
    gaps[2] = 0;
    gaps[3] = targetLen;

    fwbwaln->resetParams(startAnchorGo, startAnchorGe, startAnchorT);
    fwbwaln->initScoreMatrix(lolScoreMatrix, gaps);
    fwbwaln->runFwBw<false, 0>();




    for (int startAnchor = 0; startAnchor < numStartAnchors; startAnchor++) {

        maxIndexX = 0;
        maxIndexY = 0;

        float maxScore = 0;
        
        if (startAnchor % 2 == 0) {

            for (int i = startAnchorLength; i < queryLen - startAnchorLength; ++i) {
                for (int j = startAnchorLength; j < targetLen - startAnchorLength; ++j) {
                    float alignProb = fwbwaln->getProbability(i, j);
                    if (alignProb > maxScore) {
                        maxScore = alignProb;
                        maxIndexX = i;
                        maxIndexY = j;
                    }
                }
                if (maxScore >= fwbwaln->maxP){
                    break;
                }

            }
        }
        else{
            for (int i = queryLen - startAnchorLength -1; i > startAnchorLength; --i) {
                for (int j = targetLen - startAnchorLength -1; j > startAnchorLength; --j) {
                    float alignProb = fwbwaln->getProbability(i, j);
                    if (alignProb > maxScore) {
                        maxScore = alignProb;
                        maxIndexX = i;
                        maxIndexY = j;
                    }
                }
                if (maxScore >= fwbwaln->maxP){
                    break;
                }
            }
        }





        int startRow = maxIndexX - std::min(maxIndexX, maxIndexY);
        int startCol = maxIndexY - std::min(maxIndexX, maxIndexY);
        int diagLength = std::min(queryLen - startRow, targetLen - startCol);
        for (int i = 0; i < diagLength; i++) {
            lolScoreVec[i] = lolScoreMatrix[startRow + i][startCol + i];
        }
         
        for (int i = -startAnchorLength; i < startAnchorLength; i++) {
            for (int j = 0; j < diagLength; j++) {
                if (distanceij[maxIndexX+i][startRow + j] > 0) {
                    lolDist[j] = std::abs(distanceij[maxIndexX+i][startRow + j] - distancekl[maxIndexY+i][startCol+j]);
                    lolSeqDist[j] = std::copysign(1.0f, (maxIndexX+i - startRow + j)) * std::log(1 + std::abs((float)(maxIndexX+i - startRow + j)));
                } else {
                    lolDist[j] = -1;
                    lolSeqDist[j] = -1;
                }
            }

            lolScore(lolDist, lolSeqDist, lolScoreVec, diagLength, hiddenLayer);
        }
        startAnchorScores[startAnchor] = maxSubArray(lolScoreVec, diagLength);
        alignStartAnchors(anchorQuery[startAnchor], anchorTarget[startAnchor], maxIndexX, maxIndexY,
                          &newAnchorLength[startAnchor], fwbwaln->getProbabiltyMatrix(), lolScoreMatrix);
        anchorLength[startAnchor] = newAnchorLength[startAnchor];
    }
    
    for (int i = 0; i < queryLen; ++i) {
        std::memset(lolScoreMatrix[i], 0, targetLen * sizeof(lolScoreMatrix[0][0]));
    }
    
    indexSort(startAnchorScores, startAnchorIndex, numStartAnchors);

    gaps[0] = 0;
    gaps[1] = 0;
    gaps[2] = 0;
    gaps[3] = 0;

    fwbwaln->resetParams(lolGo, lolGe, lolT);
    int startAnchor;

    for (int startAnchorIter = 0; startAnchorIter < seedNumber; startAnchorIter++) {
        startAnchor = startAnchorIndex[numStartAnchors - startAnchorIter - 1];
        bool addSeq = false;

        for(int iteration = 0; iteration < 1000; iteration++){
            gaps[0] = 0;
            gaps[1] = 0;
            gaps[2] = 0;
            gaps[3] = 0;


            while((gaps[1] < queryLen && gaps[3] < targetLen)){

                calcGap(anchorQuery[startAnchor], anchorTarget[startAnchor], gaps, queryLen, targetLen);
                

                if (gaps[0] != -1){
                    lolMatrix(anchorQuery[startAnchor], anchorTarget[startAnchor], newAnchorLength[startAnchor], gaps,
                              distanceij, distancekl, lolScoreMatrix, queryLen, targetLen, hiddenLayer, lolDist);
                }else{
                    break;
                }
            }

            for (int i = 0; i < queryLen; i++) {
                if (anchorQuery[startAnchor][i] == 2) {
                    anchorQuery[startAnchor][i] = 1;
                }
            }
            for (int i = 0; i < targetLen; i++) {
                if (anchorTarget[startAnchor][i] == 2) {
                    anchorTarget[startAnchor][i] = 1;
                }
            }
            newAnchorLength[0] = 0;

            gaps[0] = minLolMatrixIndex;
            gaps[1] = 0;
            gaps[2] = 0;
            gaps[3] = 0;
            float maxP = 0;

            if(iteration == 0){
                maxP = 0.5;
            }
            else{
                maxP = lolMinP;
            }
            while (gaps[1] < maxLolMatrixIndex && gaps[3] < targetLen) {
                calcGap(anchorQuery[startAnchor], anchorTarget[startAnchor], gaps, maxLolMatrixIndex, targetLen);
                if(gaps[0] != -1){
                    fwbwaln->initScoreMatrix(lolScoreMatrix, gaps);
                    fwbwaln->runFwBw<false, 0>(); 
                    float** fwbwP = fwbwaln->getProbabiltyMatrix();
                    maxP = std::max(maxP, fwbwaln->maxP);
                    if(fwbwaln->maxP == 0){
                        if (fwbwaln->temperature > 30){
                            fwbwaln->temperature = lolT;
                            break;
                        }
                       
                        fwbwaln->temperature += 2;
                        gaps[0] = 0;
                        gaps[1] = 0;
                        gaps[2] = 0;
                        gaps[3] = 0;
                        continue; 
                        
                    }
                    
                    for (int i = 0; i < gaps[1]-gaps[0]; ++i) {
                        std::copy(fwbwP[i], fwbwP[i] + (gaps[3] - gaps[2]), &probabilityMatrix[i + gaps[0]][gaps[2]]);
                    }
                    
                }
                else{
                    break;
                }         
            }
            fwbwaln->temperature = lolT;
            newAnchorLength[startAnchor] = 0;

            gaps[0] = minLolMatrixIndex;
            gaps[1] = 0;
            gaps[2] = 0;
            gaps[3] = 0;


            while (gaps[1] < maxLolMatrixIndex && gaps[3] < targetLen) {
                calcGap(anchorQuery[startAnchor], anchorTarget[startAnchor], gaps, maxLolMatrixIndex, targetLen);
                if (gaps[0] != -1) {
                    for (int i = gaps[0]; i < gaps[1]; i++) {
                        for (int j = gaps[2]; j < gaps[3]; j++) {
                            if (probabilityMatrix[i][j] > maxP - 0.1) {
                                if (anchorQuery[startAnchor][i] == 0 && anchorTarget[startAnchor][j] == 0) {
                                    anchorQuery[startAnchor][i] = 2;
                                    anchorTarget[startAnchor][j] = 2;
                                    anchorLength[startAnchor] += 1;
                                    newAnchorLength[startAnchor] += 1;
                                    break;
                                }
                            }
                        }
                    }
                } else {
                    break;
                }
            }

            if (newAnchorLength[startAnchor] == 0) {
                if (!addSeq && multiDomain == 1) {
                    addSeq = true;
                    minLolMatrixIndex = 0;
                    maxLolMatrixIndex = queryLen;
                    addForwardScoreMatrix(
                        targetNumAA,
                        targetNum3Di,
                        queryLen,
                        targetLen,
                        subMatAA,
                        lolScoreMatrix);
                    int tempAnchorIdx = startAnchorIndex[numStartAnchors - startAnchorIter - 1 - seedNumber];
                    anchorLength[tempAnchorIdx] = anchorLength[startAnchor];
                    for(int i = 0; i<queryLen; i++){
                        if(anchorQuery[startAnchor][i]!= 0){
                            anchorQuery[tempAnchorIdx][i] = 1;
                        }
                        else{
                            anchorQuery[tempAnchorIdx][i] = 0;
                        }
                    }
                    for(int i = 0; i<targetLen; i++){
                        if(anchorTarget[startAnchor][i]!= 0){
                            anchorTarget[tempAnchorIdx][i] = 1;
                        }
                        else{
                            anchorTarget[tempAnchorIdx][i] = 0;
                        }
                    }
                }
                else{
                    break;
                }

            }
        }

        for (int i = 0; i < queryLen; ++i) {
            std::memset(lolScoreMatrix[i], 0, targetLen * sizeof(lolScoreMatrix[0][0]));
        }

    }

    float maxLolScore = std::numeric_limits<float>::min();
    int maxLolIdx = 0;
    int seedIter = seedNumber;
    if(multiDomain){
        seedIter *= 2;
    }

    
    for (int startAnchorIter = 0; startAnchorIter < seedIter; startAnchorIter++){
        
        startAnchor = startAnchorIndex[numStartAnchors - startAnchorIter -1];
        int startAnchorIndex = 0;
        for(int i = 0; i < queryLen; i++){
            if(anchorQuery[startAnchor][i] != 0){
                finalAnchorQuery[startAnchorIndex] = i;
                startAnchorIndex++;
               
            }
        }
        startAnchorIndex = 0;
        for(int i = 0; i < targetLen; i++){
            if(anchorTarget[startAnchor][i] != 0){
                finalAnchorTarget[startAnchorIndex] = i;
                startAnchorIndex++;
            }
        }


        computeDiScore(targetNumAA, targetNum3Di, anchorLength[startAnchor], finalAnchorQuery, finalAnchorTarget, subMatAA, lolScoreVec);

        for (int i = 0; i < anchorLength[startAnchor]; i++) {
            for (int j = 0;j < anchorLength[startAnchor]; j++) {
                if(distanceij[finalAnchorQuery[i]][finalAnchorQuery[j]] > 0.0){
                    lolDist[j] = std::abs(distanceij[finalAnchorQuery[i]][finalAnchorQuery[j]] 
                                - distancekl[finalAnchorTarget[i]][finalAnchorTarget[j]]);

                    lolSeqDist[j] = std::copysign(1.0f, (finalAnchorQuery[i]-finalAnchorQuery[j])) 
                                    * std::log(1 + std::abs((float)(finalAnchorQuery[i]-finalAnchorQuery[j])));
                }
                else{
                    lolDist[j] = -1;
                    lolSeqDist[j] = -1;
                }
            }
            lolScore(lolDist, lolSeqDist, lolScoreVec, anchorLength[startAnchor], hiddenLayer);

        }

        float totalLolScore = 0;
        for (int i = 0; i < anchorLength[startAnchor]; i++) {
            totalLolScore += lolScoreVec[i];

        }
        if (totalLolScore > maxLolScore){
            maxLolScore = totalLolScore;
            maxLolIdx = startAnchor;
        }
    }
    float seqId = 0.0;

    int startAnchorIndex = 0;
    for(int i = 0; i < queryLen; i++){
        if(anchorQuery[maxLolIdx][i] != 0){
            finalAnchorQuery[startAnchorIndex] = i;
            startAnchorIndex++;            
        }
    }
    startAnchorIndex = 0;
    for(int i = 0; i < targetLen; i++){
        if(anchorTarget[maxLolIdx][i] != 0){
            finalAnchorTarget[startAnchorIndex] = i;
            startAnchorIndex++;
        }
    }
    
    computeDiScore(targetNumAA, targetNum3Di, anchorLength[maxLolIdx], finalAnchorQuery, finalAnchorTarget, subMatAA, lolScoreVec);
    float maxAA3DiScore = 0;
    for (int i = 0; i < anchorLength[maxLolIdx]; i++) {
        maxAA3DiScore += lolScoreVec[i];
    } 

    for(int i = 0; i < queryLen; i++){
        lolScoreVecSh[i] = 0;
        
    }
    for (int i = 0; i < anchorLength[maxLolIdx]; i++) {
        for (int j = 0;j < anchorLength[maxLolIdx]; j++) {
            if(distanceij[finalAnchorQuery[i]][finalAnchorQuery[j]] > 0.0){
                lolDist[j] = 0;
                lolSeqDist[j] = std::copysign(1.0f, (finalAnchorQuery[i]-finalAnchorQuery[j])) * std::log(1 + std::abs((float)(finalAnchorQuery[i]-finalAnchorQuery[j])));
            }
            else{
                lolDist[j] = -1;
                lolSeqDist[j] = -1;
            }
        }
        lolScore(lolDist, lolSeqDist, lolScoreVecSh, anchorLength[maxLolIdx], hiddenLayer);
    }


    for (int i = 0; i < anchorLength[maxLolIdx]; i++) {
        for (int j = 0;j < anchorLength[maxLolIdx]; j++) {
            if(distanceij[finalAnchorQuery[i]][finalAnchorQuery[j]] > 0.0){
                lolDist[j] = std::abs(distanceij[finalAnchorQuery[i]][finalAnchorQuery[j]] - distancekl[finalAnchorTarget[i]][finalAnchorTarget[j]]);
                lolSeqDist[j] = std::copysign(1.0f, (finalAnchorQuery[i]-finalAnchorQuery[j])) * std::log(1 + std::abs((float)(finalAnchorQuery[i]-finalAnchorQuery[j])));
            }
            else{
                lolDist[j] = -1;
                lolSeqDist[j] = -1;
            }
        }
        lolScore(lolDist, lolSeqDist, lolScoreVec, anchorLength[maxLolIdx], hiddenLayer);

    }

    float normalizedLolscoreSelfhit = 0;
    float nanCheck = 0;
    maxLolScore = 0;
    for (int i = 0; i < anchorLength[maxLolIdx]; i++) {
        if(lolScoreVecSh[i] != 0){
            maxLolScore += lolScoreVec[i];
            nanCheck = lolScoreVec[i] / (lolScoreVecSh[i]);
            normalizedLolscoreSelfhit +=  nanCheck == nanCheck ? nanCheck : 0;
        }
        
    }

    std::string backtrace = "";
    int matches = 0;
    int queryCount = 0;
    int targetCount = 0;
    
    while(matches < anchorLength[maxLolIdx]){
        if(anchorQuery[maxLolIdx][queryCount] != 0 && anchorTarget[maxLolIdx][targetCount] != 0){
            backtrace.append(1, 'M');
            matches++;
            if(querySeq[queryCount] == targetSeq[targetCount]){
                seqId += 1.0;
            }
            queryCount++;
            targetCount++;
        }
        else if(anchorTarget[maxLolIdx][targetCount] == 0){
            backtrace.append(1, 'D');
            targetCount++;
        }
        else if(anchorQuery[maxLolIdx][queryCount] == 0){
            backtrace.append(1, 'I');
            queryCount++;
        }

    }

    Matcher::result_t result = Matcher::result_t();
    result.seqId = seqId / (float)anchorLength[maxLolIdx];
    result.qcov = anchorLength[maxLolIdx] / (float)queryLen;
    result.dbcov = anchorLength[maxLolIdx] / (float)targetLen;
    if(multiDomain == 0){
        result.eval = (((maxLolScore +3 * maxAA3DiScore) * normalizedLolscoreSelfhit/anchorLength[maxLolIdx])/qqScore)/ std::pow(queryLen * targetLen, 0.25);
        result.score = (((maxLolScore +3 * maxAA3DiScore) * normalizedLolscoreSelfhit/anchorLength[maxLolIdx]))/ std::pow(queryLen * targetLen, 0.25);
    }
    else{
        result.eval = (((maxLolScore +3 * maxAA3DiScore) * normalizedLolscoreSelfhit/anchorLength[maxLolIdx])/qqScore);
        result.score = (((maxLolScore +3 * maxAA3DiScore) * normalizedLolscoreSelfhit/anchorLength[maxLolIdx]));
    }
    result.dbKey = dbKey;
    result.qStartPos = 0;
    result.dbStartPos = 0;
    result.qEndPos = queryLen - 1;
    result.dbEndPos = targetLen - 1;
    result.qLen = queryLen;
    result.dbLen = targetLen;
    result.alnLength = backtrace.length();
    result.backtrace = Matcher::compressAlignment(backtrace);

    
    // trim backtrace find the first 'M'
    size_t firstM = targetLen;
    int qLeftSkip = 0;
    int tLeftSkip = 0;
    for (size_t i = 0; i < backtrace.size(); i++) {
        if (backtrace[i] == 'M') {
            firstM = i;
            break;
        }
        if (backtrace[i] == 'I') {
            qLeftSkip++;
        } else if (backtrace[i] == 'D') {
            tLeftSkip++;
        }
    }

    result.qStartPos  += qLeftSkip;
    result.dbStartPos += tLeftSkip;
    result.qEndPos = result.qStartPos;
    result.dbEndPos = result.dbStartPos;
    for (size_t i = firstM; i < backtrace.size(); i++) {
        switch (backtrace[i]) {
            case 'M': // match consumes 1 base in query + 1 base in target
                result.qEndPos++;
                result.dbEndPos++;
                break;
            case 'I': // insertion consumes 1 base in query
                result.qEndPos++;
                break;
            case 'D': // deletion consumes 1 base in target
                result.dbEndPos++;
                break;
            default:
                // handle unexpected cigar chars if needed
                break;
        }
    }


    result.qEndPos--;
    result.dbEndPos--;
    result.backtrace = Matcher::compressAlignment(backtrace.substr(firstM));
    result.alnLength = (int)result.backtrace.size();

    return result;

}
void LoLAlign::alignStartAnchors(int* anchorQuery, int* anchorTarget, int maxQuery, int maxTarget,
                                 int* anchorLength, float** fwbwP, float** lolScoreMatrix) {

    for (int i = maxQuery - this->startAnchorLength; i <= maxQuery + this->startAnchorLength; ++i) {
        int maxInd = maxTarget + i - maxQuery;

        anchorQuery[i] = 2;
        anchorTarget[maxInd] = 2;
        fwbwP[i][maxInd] = 0;
        lolScoreMatrix[i][maxInd] = 0;
        *anchorLength = *anchorLength + 1;
    }
}

void LoLAlign::calcDistMatrix(float* x, float* y, float* z, size_t len, float** d, bool cutoff) {
    const float cutoffDistance = 20.0f;
    const float cutoffSq = cutoffDistance * cutoffDistance;

    for (size_t i = 0; i < len; ++i) {
        d[i][i] = 0.0f;
        for (size_t j = i + 1; j < len; ++j) {
            float dx = x[i] - x[j];
            float dy = y[i] - y[j];
            float dz = z[i] - z[j];

            float distSq = dx * dx + dy * dy + dz * dz;

            if (cutoff && distSq > cutoffSq) {
                d[i][j] = 0.0f;
            } else {
                float dist = std::sqrt(distSq);
                d[i][j] = dist;
            }
            d[j][i] = d[i][j];
        }
    }
}

void LoLAlign::computeDiScore(
        unsigned char* targetNumAA,
        unsigned char* targetNum3Di,
        int anchorLen,
        int* finalAnchorQuery,
        int* finalAnchorTarget,
        SubstitutionMatrix& subMatAA,
        float* scoreForward)
{


    for (int i = 0; i < anchorLen; ++i) {
        int q = finalAnchorQuery[i];
        int t = finalAnchorTarget[i];
        scoreForward[i] = static_cast<float>((subMatAA.subMatrix[queryNumAA[q]][targetNumAA[t]] * 1.4)
                                             + (scoringMatrix3Di[queryNum3Di[q]][targetNum3Di[t]] * 2.1));
    }
}

void LoLAlign::initQuery(float* x, float* y, float* z, Sequence& qSeqAA, Sequence& qSeq3Di,
                         int queryLen, SubstitutionMatrix& subMatAA, int maxTLen, int multiDomain)
{
    memcpy(queryX, x, sizeof(float) * queryLen);
    memcpy(queryY, y, sizeof(float) * queryLen);
    memcpy(queryZ, z, sizeof(float) * queryLen);

    this->queryLen = queryLen;
    this->querySeq = qSeqAA.getSeqData();
    this->query3DiSeq = qSeq3Di.getSeqData();
    queryNumAA = qSeqAA.numSequence;
    queryNum3Di = qSeq3Di.numSequence;
    Coordinates queryCaCords;
    queryCaCords.x = queryX;
    queryCaCords.y = queryY;
    queryCaCords.z = queryZ;

    if (finalAnchorQuery != nullptr) {
        delete[] finalAnchorQuery;
    }
    if (finalAnchorTarget != nullptr) {
        delete[] finalAnchorTarget;
    }

    finalAnchorQuery = new int[std::max(queryLen, maxTLen)];
    finalAnchorTarget = new int[std::max(queryLen, maxTLen)];
    maxAnchorLen = std::max(queryLen, maxTLen);
    if(queryLen < 10){
        startAnchorLength = 0;
    }
    else{
        startAnchorLength = 3;
    }
    calcDistMatrix(queryX, queryY, queryZ, queryLen, distanceij, true);
    for (int i = 0; i < queryLen; i++) {
        finalAnchorQuery[i] = i;
    }


    computeDiScore(queryNumAA, queryNum3Di, queryLen, finalAnchorQuery, finalAnchorQuery, subMatAA, lolScoreVec);
    float diScore = 0;
    qqScore = 0;
    for (int i = 0; i < queryLen; i++) {
        diScore += lolScoreVec[i];
    }

    for (int i = 0; i < queryLen; i++) {
        for (int j = 0; j < queryLen; j++) {
            if (distanceij[i][j] > 0) {
                lolDist[j] = 0;
                lolSeqDist[j] = std::copysign(1.0f, (i-j)) * std::log(1 + std::abs((float)(i-j)));
            } else {
                lolDist[j] = -1;
                lolSeqDist[j] = -1;
            }
        }
        lolScore(lolDist, lolSeqDist, lolScoreVec, queryLen, hiddenLayer);
    }



    for (int i = 0; i < queryLen; i++) {
        qqScore += lolScoreVec[i];
    }
    qqScore = (qqScore + 3 * diScore);
    if (multiDomain == 0) {
        qqScore = qqScore / std::pow(queryLen * queryLen, 0.25);
    }
    
    return;
}

void LoLAlign::lolMatrix(int* anchorQuery, int* anchorTarget, int anchorLength,
     int* gaps, float** distanceij, float** distancekl, float** lolScoreMatrix,
      int queryLen, int targetLen, float** hiddenLayer, float* dDist)
{
    int gap0Start = gaps[0];
    int gap0End = gaps[1];
    int gap1Start = gaps[2];
    int gap1End = gaps[3];
    float seqDist = 0.0;
    int anchorQ = 0;
    int anchorT = 0;

    for (int i = 0; i < anchorLength; i++) {
        for (int aq = anchorQ + 1; aq < queryLen; aq++) {
            if (anchorQuery[aq] == 2) {
                anchorQ = aq;
                break;
            }
        }
        for (int at = anchorT + 1; at < targetLen; at++) {
            if (anchorTarget[at] == 2) {
                anchorT = at;
                break;
            }
        }
        

        for (int j = gap0Start; j < gap0End; j++) {
            float dq = distanceij[anchorQ][j];

            if (dq > 0) {
                minLolMatrixIndex = std::min(minLolMatrixIndex, j);
                maxLolMatrixIndex = std::max(maxLolMatrixIndex, j + 1);
                seqDist = std::copysign(1.0f, (anchorQ - j)) * std::log(1 + std::abs((float)(anchorQ - j)));
                for (int l = gap1Start; l < gap1End; l++) {
                    dDist[l - gap1Start] = std::abs(dq - distancekl[anchorT][l]);
                }
                lolScore(dDist, seqDist, lolScoreMatrix[j], gap1End - gap1Start, gap1Start, hiddenLayer);
            }
        }
    }
    return;
}


void LoLAlign::lolScore(float* dDist, float dSeq, float* score, int length, int start, float** hiddenLayer)
{
    simd_float zero = simdf32_setzero();
    
    simd_float seq0 = simdf32_mul(simdf32_set(dSeq), simdf32_set(w1[0][0]));
    simd_float seq1 = simdf32_mul(simdf32_set(dSeq), simdf32_set(w1[0][1]));
    simd_float seq2 = simdf32_mul(simdf32_set(dSeq), simdf32_set(w1[0][2]));
    int i = 0;
    for (; i <= length - VECSIZE_FLOAT; i += VECSIZE_FLOAT) {
        simd_float dDistVec = simdf32_loadu(&dDist[i]);

        simd_float hl0 = simdf32_add(seq0, simdf32_fmadd(dDistVec, simdf32_set(w1[1][0]), simdf32_set(b1[0])));
        simd_float hl1 = simdf32_add(seq1, simdf32_fmadd(dDistVec, simdf32_set(w1[1][1]), simdf32_set(b1[1])));
        simd_float hl2 = simdf32_add(seq2, simdf32_fmadd(dDistVec, simdf32_set(w1[1][2]), simdf32_set(b1[2])));

        hl0 = simdf32_max(hl0, zero);
        hl1 = simdf32_max(hl1, zero);
        hl2 = simdf32_max(hl2, zero);

        simdf32_storeu(&hiddenLayer[i][0], hl0);
        simdf32_storeu(&hiddenLayer[i][1], hl1);
        simdf32_storeu(&hiddenLayer[i][2], hl2);

        simd_float scoreVec = simdf32_loadu(&score[i + start]);
        scoreVec = simdf32_fmadd(hl0, simdf32_set(w2[0]), scoreVec);
        scoreVec = simdf32_fmadd(hl1, simdf32_set(w2[1]), scoreVec);
        scoreVec = simdf32_fmadd(hl2, simdf32_set(w2[2]), scoreVec);
        scoreVec = simdf32_add(scoreVec, simdf32_set(b2));

        simdf32_storeu(&score[i + start], scoreVec);
    }

    for (; i < length; ++i) {
        for (int k = 0; k < 3; ++k) {
            hiddenLayer[i][k] = dSeq * w1[0][k];
            hiddenLayer[i][k] += dDist[i] * w1[1][k];
            hiddenLayer[i][k] += b1[k];
            hiddenLayer[i][k] = std::max(0.0f, hiddenLayer[i][k]);
            score[i + start] += hiddenLayer[i][k] * w2[k];
        }
        score[i + start] += b2;
    }
}



void LoLAlign::lolScore(float* dist, float* dSeq, float* score, int length, float** hiddenLayer)
{

    for (int i = 0; i < length; ++i) {
        if (dist[i] >= 0) {
            for (int k = 0; k < 3; ++k) {
                hiddenLayer[i][k] = dSeq[i] * w1[0][k];
                hiddenLayer[i][k] += dist[i] * w1[1][k];
                hiddenLayer[i][k] += b1[k];
                hiddenLayer[i][k] = std::max(0.0f, hiddenLayer[i][k]);
                score[i] += hiddenLayer[i][k] * w2[k];
            }
            score[i] += b2;
        }
    }
}


void LoLAlign::computeForwardScoreMatrix(
        unsigned char* targetNumAA,
        unsigned char* targetNum3Di,
        int queryLen,
        int targetLen,
        SubstitutionMatrix& subMatAA,
        float** scoreForward)
{
    for (int i = 0; i < queryLen; ++i) {
        for (int j = 0; j < targetLen; ++j) {
            scoreForward[i][j] = static_cast<float>((subMatAA.subMatrix[queryNumAA[i]][targetNumAA[j]] * 1.4) + (scoringMatrix3Di[queryNum3Di[i]][targetNum3Di[j]] * 2.1));
        }
    }
}

void LoLAlign::addForwardScoreMatrix(
    unsigned char* targetNumAA,
    unsigned char* targetNum3Di,
    int queryLen,
    int targetLen,
    SubstitutionMatrix& subMatAA,
    float** scoreForward)
{
    for (int i = 0; i < queryLen; ++i) {
        for (int j = 0; j < targetLen; ++j) {
            scoreForward[i][j] += static_cast<float>((subMatAA.subMatrix[queryNumAA[i]][targetNumAA[j]] * 1.4) + (scoringMatrix3Di[queryNum3Di[i]][targetNum3Di[j]] * 2.1));
        }
    }
}

