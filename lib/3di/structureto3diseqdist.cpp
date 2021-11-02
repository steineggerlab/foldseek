#include <string.h>

#include <iostream>
#include <vector>
#include <math.h>
#include <limits.h>
#include "structureto3diseqdist.h"

StructureTo3diSeqDist::StructureTo3diSeqDist(){}

void StructureTo3diSeqDist::discretizeSeqDistance(std::vector<char> & states,
                                                  std::vector<int> & partnerIdx,
                                                  std::vector<bool> & mask, const size_t len){
    int minDistance;
    char closestState;

    for (size_t i = 0; i < len; i++){
        closestState = Alphabet3diSeqDist::INVALID_STATE;
        if (mask[i]){
            minDistance = INT_MAX;
            int seqDistance = partnerIdx[i] - i;
            for (size_t j = 0; j < Alphabet3diSeqDist::CENTROID_CNT; j++){

                int distToCentroid = abs(Alphabet3diSeqDist::centroids[j] - seqDistance);
                if (distToCentroid < minDistance){
                    closestState = j;
                    minDistance = distToCentroid;
                }
            }
        }
        states[i] = closestState;
    }
}

char * StructureTo3diSeqDist::structure2states(Vec3 * ca, Vec3 * n,
                                        Vec3 * c, Vec3 * cb,
                                        size_t len){
    states.clear();
    partnerIdx.clear();
    mask.clear();

    if(len > states.size()){
        states.resize(len);
        partnerIdx.resize(len);
        mask.resize(len);
    }
    std::fill(partnerIdx.begin(), partnerIdx.begin() + len, -1);

    replaceCBWithVirtualCenter(ca, n, c, cb, len);
    createResidueMask(mask, ca, n, c, len);
    findResiduePartners(partnerIdx, cb, mask, len);
    discretizeSeqDistance(states, partnerIdx,  mask, len);

    return states.data();
}


