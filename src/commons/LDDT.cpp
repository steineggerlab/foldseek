#include "LDDT.h"

float dist(float* arr1, float* arr2) {
    float D2 = 0;
    for(int i = 0; i < 3; i++) {
        D2 += (arr1[i] - arr2[i]) * (arr1[i] - arr2[i]);
    }
    return sqrt(D2);
}

LDDTcalculator::LDDTcalculator(unsigned int maxQueryLength, unsigned int maxTargetLength) {
    unsigned int maxAlignLength = std::max(maxQueryLength, maxTargetLength);
    max_QueryLength = maxQueryLength;
    max_TargetLength = maxTargetLength;
    max_AlignLength = maxAlignLength;

    query_pos = new fptr_t[maxQueryLength];
    for(unsigned int i = 0; i < maxQueryLength; i++) {
        query_pos[i] = new float[3];
    }
    target_pos = new fptr_t[maxTargetLength];
    for(unsigned int i = 0; i < maxTargetLength; i++) {
        target_pos[i] = new float[3];
    }    
    dists_to_score = new fptr_t[maxQueryLength];
    for(unsigned int i = 0; i < maxQueryLength; i++) {
        dists_to_score[i] = new float[maxQueryLength];
    }
    aligned_dists_to_score = new fptr_t[maxAlignLength];
    for(unsigned int i = 0; i < maxAlignLength; i++) {
        aligned_dists_to_score[i] = new float[maxAlignLength];
    }
    dist_l1 = new fptr_t[maxAlignLength];
    for(unsigned int i = 0; i < maxAlignLength; i++) {
        dist_l1[i] = new float[maxAlignLength];
    }
    score = new fptr_t[maxAlignLength];
    for(unsigned int i = 0; i < maxAlignLength; i++) {
        score[i] = new float[maxAlignLength];
    }
    reduce_score = new float[maxAlignLength];
    norm = new float[maxQueryLength];
    norm_aligned = new float[maxAlignLength];
}

LDDTcalculator::~LDDTcalculator() {
    if(query_pos) {
        for(unsigned int i = 0; i < max_QueryLength; i++) {
            delete[] query_pos[i];
        }
        delete[] query_pos;
    }
    if(target_pos) {
        for(unsigned int i = 0; i < max_TargetLength; i++) {
            delete[] target_pos[i];
        }
        delete[] target_pos;
    }
    if(dists_to_score) {
        for(unsigned int i = 0; i < max_QueryLength; i++) {
            delete[] dists_to_score[i];
        }
        delete[] dists_to_score;
    }
    if(aligned_dists_to_score) {
        for(unsigned int i = 0; i < max_AlignLength; i++) {
            delete[] aligned_dists_to_score[i];
        }
        delete[] aligned_dists_to_score;
    }
    if(dist_l1) {
        for(unsigned int i = 0; i < max_AlignLength; i++) {
            delete[] dist_l1[i];
        }
        delete[] dist_l1;
    }
    if(score) {
        for(unsigned int i = 0; i < max_AlignLength; i++) {
            delete[] score[i];
        }
        delete[] score;
    }
    if(reduce_score) delete[] reduce_score;
    if(norm) delete[] norm;
    if(norm_aligned) delete[] norm_aligned;
}

void LDDTcalculator::initQuery(unsigned int queryLen, float *qx, float *qy, float *qz) {
    queryLength = queryLen;
    for(unsigned int i = 0; i < queryLength; i++) {
        query_pos[i][0] = qx[i];
        query_pos[i][1] = qy[i];
        query_pos[i][2] = qz[i];
    }
    query_grid = Grid(query_pos, queryLength);
}

LDDTcalculator::LDDTscoreResult LDDTcalculator::computeLDDTScore(unsigned int targetLen, int qStartPos, int tStartPos, const std::string &backtrace, float *tx, float *ty, float *tz) {
    targetLength = targetLen, queryStart = qStartPos, targetStart = tStartPos, cigar = backtrace;
    for(unsigned int i = 0; i < targetLength; i++) {
        target_pos[i][0] = tx[i];
        target_pos[i][1] = ty[i];
        target_pos[i][2] = tz[i];
    }

    construct_hash_tables_align(0, qStartPos, tStartPos);
    calculate_distance();
    compute_scores();
    return LDDTscoreResult(reduce_score, alignLength);
}

void LDDTcalculator::construct_hash_tables_align(int align_idx, int query_idx, int target_idx) {
    for(int i = 0; i < query_idx; i++) query_to_align[i] = -1;
    for(int i = 0; i < target_idx; i++) target_to_align[i] = -1;
    for(std::size_t i = 0; i < cigar.length(); i++) {
        if(cigar[i] == 'M') {
            align_to_query[align_idx] = query_idx; query_to_align[query_idx] = align_idx;
            align_to_target[align_idx] = target_idx; target_to_align[target_idx] = align_idx;
            align_idx++; query_idx++; target_idx++;
        } else if(cigar[i] == 'D') {
            target_to_align[target_idx] = -1; // does not align
            target_idx++;
        } else if(cigar[i] == 'I') {
            query_to_align[query_idx] = -1; // does not align
            query_idx++;
        }
    }
    for(; query_idx < (int)queryLength; query_idx++) query_to_align[query_idx] = -1;
    for(; target_idx < (int)targetLength; target_idx++) target_to_align[target_idx] = -1;

    alignLength = align_idx;
}

void LDDTcalculator::calculate_distance() {
    // Directions
    const int DIR = 14;
    int dx[DIR] = {0,1,1, 1,1,1, 1, 1, 1, 1,0,0, 0,0};
    int dy[DIR] = {0,1,1, 1,0,0, 0,-1,-1,-1,1,1, 1,0};
    int dz[DIR] = {0,1,0,-1,1,0,-1, 1, 0,-1,1,0,-1,1};
    typedef std::multimap<std::tuple<int, int, int>, int>::iterator box_iterator;

    // Initialize arrays
    for(unsigned int i = 0; i < queryLength; i++) {
        for(unsigned int j = 0; j < queryLength; j++) {
            dists_to_score[i][j] = 0;
        }
    }
    for(unsigned int i = 0; i < alignLength; i++) {
        for(unsigned int j = 0; j < alignLength; j++) {
            aligned_dists_to_score[i][j] = 0;
        }
    }
    for(unsigned int i = 0; i < alignLength; i++) {
        for(unsigned int j = 0; j < alignLength; j++) {
            dist_l1[i][j] = 0;
        }
    }

    // Iterate through query_grid
    // Update dists_to_score, aligned_dists_to_score, dist_l1
    for(int i = 0; i <= query_grid.num_cells[0]; i++) {
        for(int j = 0; j <= query_grid.num_cells[1]; j++) {
            for(int k = 0; k <= query_grid.num_cells[2]; k++) {
                // dir = 0 (same box)
                std::pair<box_iterator, box_iterator> box_members = query_grid.box.equal_range(std::make_tuple(i, j, k));
                for(box_iterator it = box_members.first; it != box_members.second; it++) { // it->second contains query_idx corresponding to a "point"
                    for(box_iterator it2 = it; it2 != box_members.second; it2++) {
                        int query_idx1 = it->second; int query_idx2 = it2->second;
                        if(query_idx1 == query_idx2) continue; // no self-interaction
                        float distance = dist(query_pos[query_idx1], query_pos[query_idx2]);
                        bool isClose = (distance < cutoff);

                        dists_to_score[query_idx1][query_idx2] = isClose;
                        dists_to_score[query_idx2][query_idx1] = isClose;

                        int align_idx1 = query_to_align[query_idx1]; int align_idx2 = query_to_align[query_idx2]; 
                        if(align_idx1 == -1 || align_idx2 == -1) continue; // not aligned
                        aligned_dists_to_score[align_idx1][align_idx2] = isClose;
                        aligned_dists_to_score[align_idx2][align_idx1] = isClose;

                        float dist_sub = dist(target_pos[align_to_target[align_idx1]], target_pos[align_to_target[align_idx2]]);
                        float d_l = std::abs(distance - dist_sub);
                        dist_l1[align_idx1][align_idx2] = d_l;
                        dist_l1[align_idx2][align_idx1] = d_l;
                    }
                }
                // different box
                for(int dir = 1; dir < DIR; dir++) {
                    std::pair<box_iterator, box_iterator> boxPrime_members = query_grid.box.equal_range(std::make_tuple(i+dx[dir], j+dy[dir], k+dz[dir]));
                    for(box_iterator it = box_members.first; it != box_members.second; it++) {
                        for(box_iterator it2 = boxPrime_members.first; it2 != boxPrime_members.second; it2++) {
                            int query_idx1 = it->second; int query_idx2 = it2->second;
                            if(query_idx1 == query_idx2) continue; // no self-interaction
                            float distance = dist(query_pos[query_idx1], query_pos[query_idx2]);
                            bool isClose = (distance < cutoff);

                            dists_to_score[query_idx1][query_idx2] = isClose;
                            dists_to_score[query_idx2][query_idx1] = isClose;

                            int align_idx1 = query_to_align[query_idx1]; int align_idx2 = query_to_align[query_idx2]; 
                            if(align_idx1 == -1 || align_idx2 == -1) continue; // not aligned
                            aligned_dists_to_score[align_idx1][align_idx2] = isClose;
                            aligned_dists_to_score[align_idx2][align_idx1] = isClose;

                            float dist_sub = dist(target_pos[align_to_target[align_idx1]], target_pos[align_to_target[align_idx2]]);
                            float d_l = std::abs(distance - dist_sub);
                            dist_l1[align_idx1][align_idx2] = d_l;
                            dist_l1[align_idx2][align_idx1] = d_l;
                        }
                    }
                }
            }
        }
    }
}

void LDDTcalculator::compute_scores() {
    for(unsigned int i = 0; i < queryLength; i++) {
        norm[i] = 0;
    }
    for(unsigned int i = 0; i < alignLength; i++) {
        norm_aligned[i] = 0;
    }
    for(unsigned int i = 0; i < alignLength; i++) {
        reduce_score[i] = 0;
    }

    // Compute score per residue using the matrices computed above
    for(unsigned int i = 0; i < alignLength; i++) {
        for(unsigned int j = 0; j < alignLength; j++) {
            if(i == j) { // no self-interaction
                score[i][j] = 0;
            } else {
                score[i][j] = 0.25 * ((dist_l1[i][j] < 0.5) +
                                    (dist_l1[i][j] < 1.0) +
                                    (dist_l1[i][j] < 2.0) +
                                    (dist_l1[i][j] < 4.0));
            }
        }
    }

    // Normalize
    for(unsigned int col = 0; col < queryLength; col++) {
        for(unsigned int row = 0; row < queryLength; row++) {
            norm[col] += dists_to_score[row][col];
        }
        if(norm[col] != 0) norm[col] = 1 / norm[col];
        else norm[col] = INF;
    }
    for(unsigned int idx = 0; idx < alignLength; idx++) {
        norm_aligned[idx] = norm[align_to_query[idx]];
    }

    for(unsigned int i = 0; i < alignLength; i++) {
        for(unsigned int j = 0; j < alignLength; j++) {
            score[i][j] *= aligned_dists_to_score[i][j];
            reduce_score[j] += score[i][j];
        }
    }
    for(unsigned int idx = 0; idx < alignLength; idx++) {
        reduce_score[idx] *= norm_aligned[idx]; // reduce_score[] contains the lddt scores per residue
    }
}