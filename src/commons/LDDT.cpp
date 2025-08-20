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

LDDTCalculator::LDDTScoreResult LDDTCalculator::computeLDDTScore(unsigned int targetLen, int qStartPos, int tStartPos, const std::string &backtrace,
                                                                 float *tx, float *ty, float *tz) {
    targetLength = targetLen;

    for(unsigned int i = 0; i < targetLength; i++) {
        target_coordinates[i][0] = tx[i];
        target_coordinates[i][1] = ty[i];
        target_coordinates[i][2] = tz[i];
    }

    constructAlignHashes(qStartPos, tStartPos, backtrace);
    calculateLddtScores();
    return LDDTScoreResult(reduce_score, alignLength);
}

void LDDTCalculator::constructAlignHashes(int query_idx, int target_idx, const std::string & cigar) {
    memset(query_to_align, -1, sizeof(int) * queryLength);
    memset(target_to_align, -1, sizeof(int) * targetLength);
    int align_idx = 0;
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

void LDDTCalculator::calculateLddtScores() {
    const int DIR = 14;
    int dx[DIR] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0};
    int dy[DIR] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1, 1, 1, 1, 0};
    int dz[DIR] = {0, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1};
    memset(reduce_score, 0, sizeof(float) * alignLength);
    std::map<std::tuple<int, int, int>, bool> visited_boxes;
    // Create a set of unique boxes containing aligned residues
    for (unsigned int align_idx = 0; align_idx < alignLength; align_idx++) {
        if (align_to_query[align_idx] != -1) {
            std::tuple<int, int, int> box_coord = query_grid.getGridCoordinates(query_coordinates[align_to_query[align_idx]]);
            // If box_coord has been visited already, skip
            if (visited_boxes.find(box_coord) != visited_boxes.end()) {
                continue;
            }
            visited_boxes[box_coord] = true;

            std::pair<size_t, size_t> box_members = query_grid.getBoxMemberRange(box_coord);

            for (size_t i = box_members.first; i < box_members.second; i++) {
                int query_idx1 = query_grid.box[i].second;
                int align_idx1 = query_to_align[query_idx1];
                if (align_idx1 == -1) {
                    continue;
                }

                // Different boxes
                for (int dir = 0; dir < DIR; dir++) {
                    std::tuple<int, int, int> key = std::make_tuple(std::get<0>(box_coord) + dx[dir],
                                    std::get<1>(box_coord) + dy[dir],
                                    std::get<2>(box_coord) + dz[dir]);
                    std::pair<size_t, size_t> boxPrime_members = query_grid.getBoxMemberRange(key);

                    size_t i2 = (dir == 0) ? i + 1 : boxPrime_members.first;
                    for (; i2 < boxPrime_members.second; i2++) {
                        int query_idx2 = query_grid.box[i2].second;;
                        int align_idx2 = query_to_align[query_idx2];
                        if (align_idx2 == -1) {
                            continue;
                        }
                        if (dists_to_score[query_idx1][query_idx2]) {
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

    for(unsigned int idx = 0; idx < alignLength; idx++) {
        reduce_score[idx] *= norm[align_to_query[idx]]; // reduce_score[] contains the lddt scores per residue
    }
}
