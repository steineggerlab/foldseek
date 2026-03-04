#include <string.h>

#include <iostream>
#include <vector>
#include <cmath>
#include "structureto12st.h"

// StructureTo12St

StructureTo12St::StructureTo12St(){
    // Constructor
}

// Override replaceCBWithVirtualCenter to use 12St-specific VIRTUAL_CENTER parameters
void StructureTo12St::replaceCBWithVirtualCenter(Vec3 * ca, Vec3 * n,
                                Vec3 * c, Vec3 * cb, const size_t len){
    // Fix CB positions and create virtual center, but don't modify cb array
    virtualCenter.resize(len);
    for (size_t i = 0; i < len; i++){
        Vec3 cbPos = cb[i];
        if (std::isnan(cb[i].x)){
            cbPos = approxCBetaPosition(ca[i], n[i], c[i]);
            cb[i] = cbPos;  // Fix missing CB positions
        }
        Vec3 virtCenter = calcVirtualCenter(ca[i], cbPos, n[i],
                Alphabet12St::VIRTUAL_CENTER.alpha,
                Alphabet12St::VIRTUAL_CENTER.beta,
                Alphabet12St::VIRTUAL_CENTER.d);
        // Store virtual center in separate array
        virtualCenter[i] = virtCenter;
    }
}
void StructureTo12St::findResiduePartners(std::vector<int> & partnerIdx, std::vector<double> & partnerDistances,
                                          Vec3 * ca, std::vector<Vec3> & virtualCenter, std::vector<bool> & validMask, const size_t n){
    // Pick for each residue the closest neighbour as partner
    // (in terms of distances between their virtual centers).
    //
    // Ignore the first/last and invalid residues.
    // Store the minimum distance for each residue for later use as a feature.
    for(size_t i = 1; i < n - 1; i++){
        if (!validMask[i]){
            continue;  // Skip if this residue is already invalid
        }

        double minDistance = INFINITY;
        for(size_t j = 1; j < n - 1; j++){
            if (i != j && validMask[j]){
                double dist = calcDistanceBetween(virtualCenter[i], ca[j]);  // Distance between virtualCenter and Ca is correct!
                if (dist < minDistance){
                    minDistance = dist;
                    partnerIdx[i] = static_cast<int>(j);
                }
            }
        }
        partnerDistances[i] = minDistance;  // Needed as feature.
        //std::cout << "  residue i=" << i << " -> partnerIdx=" << partnerIdx[i]
        //          << " distance=" << minDistance
        //          << " virtualCenter[i] " << virtualCenter[i].x << ","
        //          << virtualCenter[i].y << ","
        //          << virtualCenter[i].z
        //          << std::endl;

        if (partnerIdx[i] == -1){  // no partner found
            validMask[i] = 0;
        }
    }
}

// Describe interaction of residue i and j
StructureTo12St::Feature StructureTo12St::calcFeatures(Vec3 * ca, int i, int j){
    double features[Alphabet12St::FEATURE_CNT];
    features[0] = copysign(log(fabs(j - i) + 1), j - i );
    features[1] = partnerDistances[i];
    return Feature(features);
}

void StructureTo12St::calcConformationDescriptors(std::vector<Feature> & features, std::vector<int> & partnerIdx,
                                                 Vec3 * ca, std::vector<bool> & mask,
                                                 const size_t len){
    // Note: Could calculate feautres for the first and last residue as well,
    // but it does not really matter.
    for (size_t i = 1; i < len - 1; i++){
        int j = partnerIdx[i];
        if (mask[i] && mask[j]){
            features[i] = calcFeatures(ca, i, j);
        } else {
            mask[i] = 0;
        }
    }

    // The first and last residue are always masked.
    mask[0] = 0;
    mask[len - 1] = 0;
}

void StructureTo12St::predictStates(std::vector<char> & states, std::vector<Feature> & features,
                                     std::vector<bool> & mask, const size_t len){
    // DEBUG: Print features before prediction
    //std::cout << "DEBUG 12St Features (len=" << len << "):" << std::endl;
    //for (size_t i = 0; i < len; i++){
    //    if (mask[i]){
    //        std::cout << "  pos=" << i
    //                  << " feature[0]=" << features[i].f[0]
    //                  << " feature[1]=" << features[i].f[1] << std::endl;
    //    }
    //}

    // Predict states using linear transformation: output = W * features + bias
    // Then take argmax over 12 states
    for (size_t i = 0; i < len; i++){
        if (mask[i]){
            double maxLogit = -INFINITY;
            char maxState = Alphabet12St::INVALID_STATE;

            // Compute logits for each of the 12 states
            for (size_t state = 0; state < Alphabet12St::STATE_CNT; state++){
                double logit = Alphabet12St::LAYER1_BIAS[state];
                for (size_t feat = 0; feat < Alphabet12St::FEATURE_CNT; feat++){
                    logit += Alphabet12St::LAYER1_W[state][feat] * features[i].f[feat];
                }

                // Track the state with maximum logit
                if (logit > maxLogit){
                    maxLogit = logit;
                    maxState = static_cast<char>(state);
                }
            }
            states[i] = maxState;
        } else {
            states[i] = Alphabet12St::INVALID_STATE;
        }
    }
}

char * StructureTo12St::structure2states(Vec3 * ca, Vec3 * n,
                                        Vec3 * c, Vec3 * cb,
                                        size_t len){
    features.clear();
    states.clear();
    partnerIdx.clear();
    partnerDistances.clear();
    mask.clear();
    virtualCenter.clear();

    if(len > features.size()){
        features.resize(len);
        states.resize(len);
        partnerIdx.resize(len);
        partnerDistances.resize(len);
        mask.resize(len);
    }
    std::fill(partnerIdx.begin(), partnerIdx.begin() + len, -1);
    std::fill(partnerDistances.begin(), partnerDistances.begin() + len, INFINITY);

    replaceCBWithVirtualCenter(ca, n, c, cb, len);
    createResidueMask(mask, ca, n, c, len);
    findResiduePartners(partnerIdx, partnerDistances, ca, virtualCenter, mask, len);  // calculates partnerDistances
    calcConformationDescriptors(features, partnerIdx, ca, mask, len);
    predictStates(states, features, mask, len);

    return states.data();
}
