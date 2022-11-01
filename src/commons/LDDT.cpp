#include "LDDT.h"
#include <string.h>
#include <algorithm>


static inline float dist(float* arr1, float* arr2) {
    float D2 = 0;
    for(int i = 0; i < 3; i++) {
        D2 += (arr1[i] - arr2[i]) * (arr1[i] - arr2[i]);
    }
    return sqrt(D2);
}

static bool compareByFirstKey(const std::pair<std::tuple<int, int, int>, int>& a, const std::pair<std::tuple<int, int, int>, int>& b) {
    return a.first < b.first;
}

LDDTCalculator::LDDTCalculator(unsigned int maxQueryLength, unsigned int maxTargetLength)
    : maxQueryLength(maxQueryLength), maxTargetLength(maxTargetLength) {
    maxAlignLength = std::max(maxQueryLength, maxTargetLength);
    query_coordinates = new float*[maxQueryLength];
    for(unsigned int i = 0; i < maxQueryLength; i++) {
        query_coordinates[i] = new float[3];
    }
    target_coordinates = new float*[maxTargetLength];
    for(unsigned int i = 0; i < maxTargetLength; i++) {
        target_coordinates[i] = new float[3];
    }    
    dists_to_score = new bool*[maxQueryLength];
    for(unsigned int i = 0; i < maxQueryLength; i++) {
        dists_to_score[i] = new bool[maxQueryLength];
    }

    score = new float*[maxAlignLength];
    for(unsigned int i = 0; i < maxAlignLength; i++) {
        score[i] = new float[maxAlignLength];
    }
    reduce_score = new float[maxAlignLength];
    norm = new float[maxQueryLength];
    query_to_align = new int[maxQueryLength];
    target_to_align = new int[maxTargetLength];
    align_to_query = new int[maxAlignLength];
    align_to_target = new int[maxAlignLength];
}

LDDTCalculator::~LDDTCalculator() {
    if(query_coordinates) {
        for(unsigned int i = 0; i < maxQueryLength; i++) {
            delete[] query_coordinates[i];
        }
        delete[] query_coordinates;
    }
    if(target_coordinates) {
        for(unsigned int i = 0; i < maxTargetLength; i++) {
            delete[] target_coordinates[i];
        }
        delete[] target_coordinates;
    }
    if(dists_to_score) {
        for(unsigned int i = 0; i < maxQueryLength; i++) {
            delete[] dists_to_score[i];
        }
        delete[] dists_to_score;
    }
    if(score) {
        for(unsigned int i = 0; i < maxAlignLength; i++) {
            delete[] score[i];
        }
        delete[] score;
    }
    if(reduce_score) {
        delete[] reduce_score;
    }
    if(norm) {
        delete[] norm;
    }
    if(query_to_align) {
        delete[] query_to_align;
    }
    if(target_to_align) {
        delete[] target_to_align;
    }
    if(align_to_query) {
        delete[] align_to_query;
    }
    if(align_to_target) {
        delete[] align_to_target;
    }
}

void LDDTCalculator::initQuery(unsigned int queryLen, float *qx, float *qy, float *qz) {
    queryLength = queryLen;
    for(unsigned int i = 0; i < queryLength; i++) {
        query_coordinates[i][0] = qx[i];
        query_coordinates[i][1] = qy[i];
        query_coordinates[i][2] = qz[i];
    }
    // Initialize arrays
    for(unsigned int i = 0; i < queryLength; i++) {
        memset(dists_to_score[i], 0, sizeof(bool) * queryLength);
    }

    query_grid = Grid(query_coordinates, queryLength);
    memset(norm, 0, sizeof(float) * queryLength);

    for(unsigned int col = 0; col < queryLength; col++) {
        for (unsigned int row = 0; row < queryLength; row++) {
            float distance = dist(query_coordinates[row], query_coordinates[col]);
            bool isClose = (col != row) && (distance < CUTOFF);
            dists_to_score[col][row] = isClose;
            dists_to_score[row][col] = isClose;
            norm[col] += dists_to_score[row][col];
        }

        if(norm[col] != 0) {
            norm[col] = 1 / norm[col];
        } else {
            norm[col] = INF;
        }
    }

}

LDDTCalculator::LDDTScoreResult LDDTCalculator::computeLDDTScore(unsigned int targetLen, int qStartPos, int tStartPos, const std::string &backtrace, float *tx, float *ty, float *tz) {
    targetLength = targetLen;
    cigar = backtrace;

    for(unsigned int i = 0; i < targetLength; i++) {
        target_coordinates[i][0] = tx[i];
        target_coordinates[i][1] = ty[i];
        target_coordinates[i][2] = tz[i];
    }

    constructAlignHashes(0, qStartPos, tStartPos);
    calculateDistance();
    computeScores();
    return LDDTScoreResult(reduce_score, alignLength);
}

void LDDTCalculator::constructAlignHashes(int align_idx, int query_idx, int target_idx) {
    memset(query_to_align, -1, sizeof(int) * queryLength);
    memset(target_to_align, -1, sizeof(int) * targetLength);
    for(std::size_t i = 0; i < cigar.length(); i++) {
        if(cigar[i] == 'M') {
            align_to_query[align_idx] = query_idx;
            query_to_align[query_idx] = align_idx;
            align_to_target[align_idx] = target_idx;
            target_to_align[target_idx] = align_idx;
            align_idx++; query_idx++; target_idx++;
        } else if(cigar[i] == 'D') {
            target_to_align[target_idx] = -1; // does not align
            target_idx++;
        } else if(cigar[i] == 'I') {
            query_to_align[query_idx] = -1; // does not align
            query_idx++;
        }
    }
    alignLength = align_idx;
}

void LDDTCalculator::calculateDistance() {
    // Directions
    const int DIR = 14;
    int dx[DIR] = {0,1,1, 1,1,1, 1, 1, 1, 1,0,0, 0,0};
    int dy[DIR] = {0,1,1, 1,0,0, 0,-1,-1,-1,1,1, 1,0};
    int dz[DIR] = {0,1,0,-1,1,0,-1, 1, 0,-1,1,0,-1,1};
    typedef std::vector<std::pair<std::tuple<int, int, int>, int>>::const_iterator box_iterator;
    memset(reduce_score, 0, sizeof(float) * alignLength);

    // Iterate through query_grid
    // Update dists_to_score, aligned_dists_to_score, dist_l1
    for(int i = 0; i <= query_grid.num_cells[0]; i++) {
        for(int j = 0; j <= query_grid.num_cells[1]; j++) {
            for(int k = 0; k <= query_grid.num_cells[2]; k++) {
                // dir = 0 (same box)
                std::pair<std::tuple<int, int, int>, int> ref;
                ref.first = std::make_tuple(i, j, k);
                ref.second = 0;

                std::pair<box_iterator, box_iterator> box_members = std::equal_range(query_grid.box.begin(), query_grid.box.end(), ref, compareByFirstKey);
                for(box_iterator it = box_members.first; it != box_members.second; it++) { // it->second contains query_idx corresponding to a "point"
                    int query_idx1 = it->second;
                    int align_idx1 = query_to_align[query_idx1];
                    if(align_idx1 == -1) {
                        continue;
                    }
                    for(box_iterator it2 = it; it2 != box_members.second; it2++) {
                        int query_idx2 = it2->second;
                        int align_idx2 = query_to_align[query_idx2];
                        if(align_idx2 == -1) {
                            continue; // not aligned
                        }
                        if(dists_to_score[query_idx1][query_idx2]) {
                            float distance = dist(query_coordinates[query_idx1], query_coordinates[query_idx2]);
                            float dist_sub = dist(target_coordinates[align_to_target[align_idx1]], target_coordinates[align_to_target[align_idx2]]);
                            float d_l = std::abs(distance - dist_sub);
                            float score = 0.25 * ((d_l < 0.5) + (d_l < 1.0) + (d_l < 2.0) + (d_l < 4.0));
                            reduce_score[align_idx2] += score;
                            reduce_score[align_idx1] += score;
                        }
                    }
                }
                // different box
                for(int dir = 1; dir < DIR; dir++) {
                    std::pair<std::tuple<int, int, int>, int> ref;
                    ref.first = std::make_tuple(i+dx[dir], j+dy[dir], k+dz[dir]);
                    ref.second = 0;
                    std::pair<box_iterator, box_iterator> boxPrime_members = std::equal_range(query_grid.box.begin(), query_grid.box.end(), ref, compareByFirstKey);
                    for(box_iterator it = box_members.first; it != box_members.second; it++) {
                        int query_idx1 = it->second;
                        int align_idx1 = query_to_align[query_idx1];
                        if(align_idx1 == -1) {
                            continue;
                        }
                        for(box_iterator it2 = boxPrime_members.first; it2 != boxPrime_members.second; it2++) {
                            int query_idx2 = it2->second;
                            int align_idx2 = query_to_align[query_idx2];
                            if(align_idx2 == -1) {// not aligned
                                continue;
                            }
                            if(dists_to_score[query_idx1][query_idx2]) {
                                float distance = dist(query_coordinates[query_idx1], query_coordinates[query_idx2]);
                                float dist_sub = dist(target_coordinates[align_to_target[align_idx1]],
                                                      target_coordinates[align_to_target[align_idx2]]);
                                float d_l = std::abs(distance - dist_sub);
                                float score = 0.25 * ((d_l < 0.5) + (d_l < 1.0) + (d_l < 2.0) + (d_l < 4.0));
                                reduce_score[align_idx2] += score;
                                reduce_score[align_idx1] += score;
                            }
                        }
                    }
                }
            }
        }
    }
}

void LDDTCalculator::computeScores() {

    for(unsigned int idx = 0; idx < alignLength; idx++) {
        reduce_score[idx] *= norm[align_to_query[idx]]; // reduce_score[] contains the lddt scores per residue
    }
}
