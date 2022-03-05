#include <string.h>

#include <iostream>
#include <vector>
#include <cmath>
#include "structureto3di.h"
#include "encoder_weights_3di.kerasify.h"

Vec3 StructureTo3DiBase::add(Vec3 a, Vec3 b){
    a.x = a.x + b.x;
    a.y = a.y + b.y;
    a.z = a.z + b.z;
    return a;
}

Vec3 StructureTo3DiBase::sub(Vec3 a, Vec3 b){
    a.x = a.x - b.x;
    a.y = a.y - b.y;
    a.z = a.z - b.z;
    return a;
}

Vec3 StructureTo3DiBase::norm(Vec3 a){
    double len = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    a.x = a.x / len;
    a.y = a.y / len;
    a.z = a.z / len;
    return a;
}

Vec3 StructureTo3DiBase::cross(Vec3 a, Vec3 b){
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x};
}

Vec3 StructureTo3DiBase::scale(Vec3 a, double f){
    a.x *= f;
    a.y *= f;
    a.z *= f;
    return a;
}

double StructureTo3DiBase::dot(Vec3 a, Vec3 b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 StructureTo3DiBase::approxCBetaPosition(Vec3 ca_atom, Vec3 n_atom, Vec3 c_atom){
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

double StructureTo3DiBase::degreeToRadians(double degree){
    return (degree / 180) * Alphabet3Di::PI;
}

Vec3 StructureTo3DiBase::calcVirtualCenter(Vec3 ca, Vec3 cb, Vec3 n,
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

double StructureTo3DiBase::calcDistanceBetween(Vec3 & a, Vec3 & b){
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

void StructureTo3DiBase::replaceCBWithVirtualCenter(Vec3 * ca, Vec3 * n,
                                Vec3 * c, Vec3 * cb, const size_t len){
    // fix CB positions and create virtual center
    for (size_t i = 0; i < len; i++){
        if (std::isnan(cb[i].x)){
            cb[i] = approxCBetaPosition(ca[i], n[i], c[i]);
        }
        Vec3 virtCenter = calcVirtualCenter(ca[i], cb[i], n[i],
                Alphabet3Di::VIRTUAL_CENTER.alpha,
                Alphabet3Di::VIRTUAL_CENTER.beta,
                Alphabet3Di::VIRTUAL_CENTER.d);
        // replace CB with virtual center
        cb[i] = virtCenter;
    }
}

void StructureTo3DiBase::createResidueMask(std::vector<bool> & validMask,
                                           Vec3 * ca, Vec3 * n, Vec3 * c,
                                           const size_t len){
    for (size_t i = 0; i < len; i++){
        if (std::isnan(ca[i].x) || std::isnan(c[i].x) || std::isnan(n[i].x)){
            validMask[i] = 0;
        } else{
            validMask[i] = 1;
        }
    }
}

void StructureTo3DiBase::findResiduePartners(std::vector<int> & partnerIdx, Vec3 * cb,
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
        if (partnerIdx[i] == -1){  // no partner found
            validMask[i] = 0;
        }
    }
}

// StructureTo3Di

StructureTo3Di::StructureTo3Di(){
    encoder.LoadModel(
            std::string((const char *)encoder_weights_3di_kerasify,
                                      encoder_weights_3di_kerasify_len));
}

// Describe interaction of residue i and j
StructureTo3Di::Feature StructureTo3Di::calcFeatures(Vec3 * ca, int i, int j){
    Vec3 u1 = norm(sub(ca[i],       ca[i - 1]));
    Vec3 u2 = norm(sub(ca[i + 1],   ca[i]));
    Vec3 u3 = norm(sub(ca[j],       ca[j - 1]));
    Vec3 u4 = norm(sub(ca[j + 1],   ca[j]));
    Vec3 u5 = norm(sub(ca[j],       ca[i]));

    double features[Alphabet3Di::FEATURE_CNT];
    features[0] = dot(u1, u2);
    features[1] = dot(u3, u4);
    features[2] = dot(u1, u5);
    features[3] = dot(u3, u5);
    features[4] = dot(u1, u4);
    features[5] = dot(u2, u3);
    features[6] = dot(u1, u3);
    features[7] = calcDistanceBetween(ca[i], ca[j]);
    features[8] = copysign(fmin(fabs(j - i), 4), j - i); // clip j-i to [-4, 4]
    features[9] = copysign(log(fabs(j - i) + 1), j - i );
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

void StructureTo3Di::encodeFeatures(std::vector<Embedding> & embeddings, std::vector<Feature> & features,
                                        std::vector<bool> & mask, const size_t len)
{
    for (size_t i = 0; i < len; i++){
        if (mask[i]){
            for (size_t j = 0; j < Alphabet3Di::FEATURE_CNT; j++){
                in.data_[j] = static_cast<float>(features[i].f[j]);
            }
            encoder.Apply(&in, &out);
            for (size_t j = 0; j < Alphabet3Di::EMBEDDING_DIM; j++){
                embeddings[i].f[j] = static_cast<double>(out.data_[j]);
            }
        }
    }
}

void StructureTo3Di::discretizeEmbeddings(std::vector<char> & states, std::vector<Embedding> & embeddings,
                                        std::vector<bool> & mask, const size_t len){
    double minDistance;
    char closestState;

    for (size_t i = 0; i < len; i++){
        closestState = Alphabet3Di::INVALID_STATE;
        if (mask[i]){
            minDistance = INFINITY;
            for (size_t j = 0; j < Alphabet3Di::CENTROID_CNT; j++){ // measure squared distance to each centroid
                double sum = 0.0;
                for (size_t k = 0; k < Alphabet3Di::EMBEDDING_DIM; k++){
                    sum += pow(embeddings[i].f[k] - Alphabet3Di::centroids[j][k], 2);
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

char * StructureTo3Di::structure2states(Vec3 * ca, Vec3 * n,
                                        Vec3 * c, Vec3 * cb,
                                        size_t len){
    features.clear();
    states.clear();
    partnerIdx.clear();
    mask.clear();
    embeddings.clear();
    in = Tensor(Alphabet3Di::FEATURE_CNT);
    out = Tensor(Alphabet3Di::EMBEDDING_DIM);

    if(len > features.size()){
        features.resize(len);
        states.resize(len);
        partnerIdx.resize(len);
        mask.resize(len);
        embeddings.resize(len);
    }
    std::fill(partnerIdx.begin(), partnerIdx.begin() + len, -1);

    replaceCBWithVirtualCenter(ca, n, c, cb, len);
    createResidueMask(mask, ca, n, c, len);
    findResiduePartners(partnerIdx, cb, mask, len);
    calcConformationDescriptors(features, partnerIdx, ca, mask, len);
    encodeFeatures(embeddings, features, mask, len);
    discretizeEmbeddings(states, embeddings,  mask, len);

    return states.data();
}


