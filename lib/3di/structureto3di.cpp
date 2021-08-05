#include <string.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "structureto3di.h"

Vec3 StructureTo3Di::add(Vec3 a, Vec3 b){
    a.x = a.x + b.x;
    a.y = a.y + b.y;
    a.z = a.z + b.z;
    return a;
}

Vec3 StructureTo3Di::sub(Vec3 a, Vec3 b){
    a.x = a.x - b.x;
    a.y = a.y - b.y;
    a.z = a.z - b.z;
    return a;
}

Vec3 StructureTo3Di::norm(Vec3 a){
    double len = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    a.x = a.x / len;
    a.y = a.y / len;
    a.z = a.z / len;
    return a;
}

Vec3 StructureTo3Di::cross(Vec3 a, Vec3 b){
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x};
}

Vec3 StructureTo3Di::scale(Vec3 a, double f){
    a.x *= f;
    a.y *= f;
    a.z *= f;
    return a;
}

double StructureTo3Di::dot(Vec3 a, Vec3 b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 StructureTo3Di::approxCBetaPosition(Vec3 ca_atom, Vec3 n_atom, Vec3 c_atom){
    // Assumption: CA forms with its four ligands a tetrahedral.
    Vec3 v1 = norm(sub(c_atom, ca_atom));
    Vec3 v2 = norm(sub(n_atom, ca_atom));

    Vec3 b1 = add(v2, scale(v1, 1/3.0));
    Vec3 b2 = cross(v1, b1);

    Vec3 u1 = norm(b1);
    Vec3 u2 = norm(b2);

    //direction from c_alpha to c_beta
    Vec3 v4 = add(
            scale(v1, -1/3.0), scale(sub(scale(u1, -1/2.0),
            scale(u2, sqrt(3)/2.0)), sqrt(8)/3.0)
            );

    return add(ca_atom, scale(v4, Alphabet3Di::DISTANCE_ALPHA_BETA));
}

double StructureTo3Di::degreeToRadians(double degree){
    return (degree / 180) * Alphabet3Di::PI;
}

Vec3 StructureTo3Di::calcVirtualCenter(Vec3 ca, Vec3 cb, Vec3 n,
                                       double alpha, double beta, double d){
    /* Apply Rodrigues rotation formula two times to CB coordinates.
     * First rotate CB by alpha in the plane CA-N-CB around CA.
     * Then the dihedral angle, d is a factor for the distance between CA and CB. */

    alpha = degreeToRadians(alpha);
    beta = degreeToRadians(beta);

    Vec3 v = sub(cb, ca);

    // normal angle (between CA-N and CA-VIRT)
    Vec3 a = sub(cb, ca);
    Vec3 b = sub(n, ca);
    Vec3 k = norm(cross(a, b)); // axis of rotation

    v = add(add(scale(v, cos(alpha)), scale(cross(k, v), sin(alpha))),
            scale(scale(k, dot(k, v)), 1 - cos(alpha)));

    // Dihedral angle (axis: CA-N, CO, VIRT)
    k = norm(sub(n, ca));
    v = add(add(scale(v, cos(beta)), scale(cross(k, v), sin(beta))),
            scale(scale(k, dot(k, v)), 1 - cos(beta)));

    Vec3 virtualCenter = add(ca, scale(v, d));
    return virtualCenter;
}

double StructureTo3Di::calcDistanceBetween(Vec3 & a, Vec3 & b){
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

void StructureTo3Di::findResiduePartners(std::vector<int> & partnerIdx, Vec3 * cb,
                                         std::vector<bool> & validMask, const size_t n){
    // Pick for each residue the closest neighbour as partner
    // (in terms of distances between their virtual centers/C_betas).
    //
    // Ignore the first/last and invalid residues.
    for(size_t i = 1; i < n - 1; i++){
        double minDistance = INFINITY;
        for(size_t j = 1; j < n - 1; j++){
            if (i != j && validMask[j]){
                double dist = calcDistanceBetween(cb[i], cb[j]);
                if (dist < minDistance){
                    minDistance = dist;
                    partnerIdx[i] = static_cast<int>(j);
                }
            }
        }
    }
}

// Describe interaction of residue i and j
StructureTo3Di::Feature StructureTo3Di::calcFeatures(Vec3 * ca, int i, int j){
    Vec3 u1 = norm(sub(ca[i],       ca[i - 1]));
    Vec3 u2 = norm(sub(ca[i + 1],   ca[i]));
    Vec3 u3 = norm(sub(ca[j],       ca[j - 1]));
    Vec3 u4 = norm(sub(ca[j + 1],   ca[j]));
    Vec3 u5 = norm(sub(ca[j],       ca[i]));

    double features[9];
    features[0] = dot(u1, u2);
    features[1] = dot(u3, u4);
    features[2] = dot(u1, u5);
    features[3] = dot(u3, u5);
    features[4] = dot(u1, u4);
    features[5] = dot(u2, u3);
    features[6] = dot(u1, u3);
    features[7] = calcDistanceBetween(ca[i], ca[j]);
    features[8] = copysign(fmin(fabs(j - i), 4), j - i); // clip j-i to [-4, 4]
    return Feature(features);
}

void StructureTo3Di::calcConformationDescriptors(std::vector<Feature> & features, std::vector<int> & partnerIdx,
                                                 Vec3 * ca, std::vector<bool> & mask,
                                                 const size_t len){
    // Calculate for each residue a number of features.
    // Check that all 6 residues are valid and update the mask accordingly.
    //  ^ this is too strict as only the CA atoms are required...
    //  Rows of invalid residues are not initialized!
    std::vector<bool> maskCopy(mask);
    for (size_t i = 0; i < len; i++){  // copy mask
        maskCopy[i] = mask[i];
    }

    for (size_t i = 1; i < len - 1; i++){
        int j = partnerIdx[i];
        if ( maskCopy[i - 1] && maskCopy[i] && maskCopy[i + 1] &&
             maskCopy[j - 1] && maskCopy[j] && maskCopy[j + 1] ){
            features[i] = calcFeatures(ca, i, j);
        } else {
            mask[i] = 0;
        }
    }

    // Resdiues without conformational descriptors are masked false.
    // Before only residues with missing coordinaes were masked.
    mask[0] = 0;
    mask[len - 1] = 0;
}

void StructureTo3Di::discretizeFeatures(std::vector<char> & states, std::vector<Feature> & features,
                                        std::vector<bool> & mask, const size_t len){
    double minDistance;
    char closestState;

    for (size_t i = 0; i < len; i++){
        closestState = Alphabet3Di::INVALID_STATE;
        if (mask[i]){
            minDistance = INFINITY;
            for (size_t j = 0; j < Alphabet3Di::CENTROID_CNT; j++){ // measure squared distance to each centroid
                double sum = 0.0;
                for (size_t k = 0; k < Alphabet3Di::FEATURE_CNT; k++){
                    sum += pow((features[i].f[k] * Alphabet3Di::feature_scaling[k]) - Alphabet3Di::centroids[j][k], 2);
                }
                if (sum < minDistance){
                    closestState = j;
                    minDistance = sum;
                }
            }
        }
        states[i] = closestState;
    }
}

void StructureTo3Di::createResidueMask(std::vector<bool> & validMask,
                                       Vec3 * ca, Vec3 * n, Vec3 * c,
                                       const size_t len){
    for (size_t i = 0; i < len; i++){
        if (isnan(ca[i].x) || isnan(c[i].x) || isnan(n[i].x)){
            validMask[i] = 0;
        } else{
            validMask[i] = 1;
        }
    }
}

char * StructureTo3Di::structure2states(Vec3 * ca, Vec3 * n,
                                        Vec3 * c, Vec3 * cb,
                                        size_t len){
    features.clear();
    states.clear();
    partnerIdx.clear();
    mask.clear();
    if(len > features.size()){
        features.resize(len);
        states.resize(len);
        partnerIdx.resize(len);
        mask.resize(len);
    }
    std::fill(partnerIdx.begin(), partnerIdx.begin() + len, -1);
    // fix CB positions and create virtual center
    for (size_t i = 0; i < len; i++){
        if (isnan(cb[i].x)){
            cb[i] = approxCBetaPosition(ca[i], n[i], c[i]);
        }
        Vec3 virtCenter = calcVirtualCenter(ca[i], cb[i], n[i], 270, 0, 2);
        // replace CB with virtual center
        cb[i] = virtCenter;
    }

    createResidueMask(mask, ca, n, c, len);
    findResiduePartners(partnerIdx, cb, mask, len);
    calcConformationDescriptors(features, partnerIdx, ca, mask, len);
    discretizeFeatures(states, features,  mask, len);

    return states.data();
}

// Profiles3DiFeatureExtractor

void StructureTo3Di::findResiduePartnersProf(std::vector<std::vector<double>> & features,
                                                                std::vector<std::vector<int>> & partnerIdx,
                                                                std::vector<Vec3> & cb,
                                                                std::vector<bool> & validMask, const size_t n){
    // Find for each residue the three closest neighbours.
    // (in terms of distances between their virtual centers/C_betas).
    // Save indices of partners and distance differences.
    //
    // Ignore the first/last and invalid residues.
    for(size_t i = 1; i < n - 1; i++){
        double minDist1 = INFINITY, minDist2 = INFINITY, minDist3 = INFINITY;

        for(size_t j = 1; j < n - 1; j++){
            if (i != j && validMask[j]){
                double dist = calcDistanceBetween(cb[i], cb[j]);
                if (dist < minDist1){
                    minDist3 = minDist2;
                    minDist2 = minDist1;
                    minDist1 = dist;
                    partnerIdx[i][2] = partnerIdx[i][1];
                    partnerIdx[i][1] = partnerIdx[i][0];
                    partnerIdx[i][0] = static_cast<int>(j);
                } else if (dist < minDist2){
                    minDist3 = minDist2;
                    minDist2 = dist;
                    partnerIdx[i][2] = partnerIdx[i][1];
                    partnerIdx[i][1] = static_cast<int>(j);
                } else if (dist < minDist3){
                    minDist3 = dist;
                    partnerIdx[i][2] = static_cast<int>(j);
                }
            }
        }

        features[i][PROFILE_FEATURE_CNT - 2] = minDist2 - minDist1; // distance difference between first and second nearest neighbour
        features[i][PROFILE_FEATURE_CNT - 1] = minDist3 - minDist2;
    }
}

void StructureTo3Di::calcConformationDescriptorsProf(std::vector<std::vector<double>> & features,
                                                                        std::vector<std::vector<int>> & partnerIdx,
                                                                        std::vector<Vec3> & ca, std::vector<bool> & mask,
                                                                        const size_t len){
    // Calculate features for each residue describing its conformation to its three nearest neighbours/interaction partners.
    //
    // features: n_residues x PROFILE_FEATURE_CNT (3 x 9 confirmation descriptors + 2 distances)
    // Note: The two distances (features[i][27], features[i][28]) are calculated in findResiduePartners.
    Feature tmp;
    std::vector<bool> maskCopy(mask);
    for (size_t i = 0; i < len; i++){  // copy mask
        maskCopy[i] = mask[i];
    }

    for (size_t i = 1; i < len - 1; i++){
        // Ensure that every residue has three neighbours
        if (partnerIdx[i][0] == -1 || partnerIdx[i][1] == -1 || partnerIdx[i][2] == -1){
            mask[i] = 0;
            // std::cerr << "Encountered -1 in partnerIdx " << i << "\n";
            continue;
        }
        for (size_t neighbor_idx = 0; neighbor_idx < NEIGHBOUR_CNT; neighbor_idx++){
            int j = partnerIdx[i][neighbor_idx];

            if ( maskCopy[i - 1] && maskCopy[i] && maskCopy[i + 1] &&  // Check validity of atom coords
                 maskCopy[j - 1] && maskCopy[j] && maskCopy[j + 1] ){
                tmp = calcFeatures(ca, i, j);
                for (size_t k = 0; k < Alphabet3Di::FEATURE_CNT; k++){
                    features[i][k + neighbor_idx * Alphabet3Di::FEATURE_CNT] = tmp.f[k] * Alphabet3Di::feature_scaling[k];
                }
            } else {
                mask[i] = 0;
            }
        }
    }

    // Resdiues without conformational descriptors are masked false.
    // Before only residues with missing coordinaes were masked.
    mask[0] = 0;
    mask[len - 1] = 0;

    // Set features for invalid residues NaN (redundant to mask)
    for (size_t i = 0; i < len; i++){
        if (mask[i] == 0){
            for (size_t k = 0; k < PROFILE_FEATURE_CNT; k++){
                features[i][k] = NAN;
            }
        }
    }
}

std::vector<std::vector<double>> StructureTo3Di::extractFeatures(Vec3 * ca, Vec3 * n,Vec3 * c, Vec3 * cb, const size_t len)
{
    //const size_t len = ca.size();
    features.clear();
    partnerIdx.clear();
    mask.clear();
    //if(len > features.size()){
    //    features.resize(len);
    //    partnerIdx.resize(len);
    //    mask.resize(len);
    //}

    //for (size_t i = 0; i < len; i++){
        //features[i].resize(PROFILE_FEATURE_CNT);
        //partnerIdx[i].resize(NEIGHBOUR_CNT);
        //std::fill(partnerIdx[i].begin(), partnerIdx[i].begin() + NEIGHBOUR_CNT, -1);
    //}

    // resizing a 2d vector?
    features.resize(len,std::vector<double>(PROFILE_FEATURE_CNT));
    partnerIdx.resize(len,std::vector<double>(NEIGHBOUR_CNT, -1));

    //replaceCBWithVirtualCenter(ca, n, c, cb, len);
    //createResidueMask(mask, ca, n, c, len);
    StructureTo3Di::findResiduePartners(partnerIdx, cb, mask, len);
    StructureTo3Di::calcConformationDescriptors(features, partnerIdx, ca, mask, len);

    return features;
}

