#include "LDDT.h"

float dist(float* arr1, float* arr2) {
    float D2 = 0;
    for(int i = 0; i < 3; i++) {
        D2 += (arr1[i] - arr2[i]) * (arr1[i] - arr2[i]);
    }
    return sqrt(D2);
}

void LDDTcalculator::initVariables(unsigned int queryLen, unsigned int targetLen, int qStartPos, int tStartPos, const std::string &backtrace) {
    queryLength = queryLen, targetLength = targetLen, queryStart = qStartPos, targetStart = tStartPos, cigar = backtrace;

    dists_to_score = new float*[queryLength];
    for(int i = 0; i < queryLength; i++) {
        dists_to_score[i] = new float[queryLength];
        for(int j = 0; j < queryLength; j++) {
            dists_to_score[i][j] = 0;
        }
    }
    aligned_dists_to_score = new float*[alignLength];
    for(int i = 0; i < alignLength; i++) {
        aligned_dists_to_score[i] = new float[alignLength];
        for(int j = 0; j < alignLength; j++) {
            aligned_dists_to_score[i][j] = 0;
        }
    }
    dist_l1 = new float*[alignLength];
    for(int i = 0; i < alignLength; i++) {
        dist_l1[i] = new float[alignLength];
        for(int j = 0; j < alignLength; j++) {
            dist_l1[i][j] = 0;
        }
    }
}

LDDTcalculator::LDDTscoreResult LDDTcalculator::computeLDDTScore(float *qx, float *qy, float *qz, float *tx, float *ty, float *tz, int qStartPos, int tStartPos) {
    query_pos = new float*[queryLength];
    for(int i = 0; i < queryLength; i++) {
        query_pos[i] = new float[3]{qx[i], qy[i], qz[i]};
    }
    target_pos = new float*[targetLength];
    for(int i = 0; i < targetLength; i++) {
        target_pos[i] = new float[3]{tx[i], ty[i], tz[i]};
    }

    construct_hash_tables_align(0, qStartPos, tStartPos, 0);
    calculate_distance();
    compute_scores();
    return LDDTscoreResult(reduce_score, alignLength);
}

void LDDTcalculator::construct_hash_tables_align(int align_idx, int query_idx, int target_idx, int cigar_idx) {
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
        cigar_idx = i+1;
    }
    for(; query_idx < queryLength; query_idx++) query_to_align[query_idx] = -1;
    for(; target_idx < targetLength; target_idx++) target_to_align[target_idx] = -1;

    alignLength = align_idx;
}

void LDDTcalculator::calculate_distance() {
    Grid query_grid(query_pos, queryLength);

    // Directions
    const int DIR = 14;
    int dx[DIR] = {0,1,1, 1,1,1, 1, 1, 1, 1,0,0, 0,0};
    int dy[DIR] = {0,1,1, 1,0,0, 0,-1,-1,-1,1,1, 1,0};
    int dz[DIR] = {0,1,0,-1,1,0,-1, 1, 0,-1,1,0,-1,1};
    typedef std::multimap<std::tuple<int, int, int>, int>::iterator box_iterator;

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
    table_t score = new float*[alignLength];
    for(int i = 0; i < alignLength; i++) {
        score[i] = new float[alignLength];
        for(int j = 0; j < alignLength; j++) {
            score[i][j] = 0;
        }
    }

    float* norm = new float[queryLength];
    float* norm_aligned = new float[alignLength]; 
    reduce_score = new float[alignLength];
    for(int i = 0; i < queryLength; i++) {
        norm[i] = 0;
    }
    for(int i = 0; i < alignLength; i++) {
        norm_aligned[i] = reduce_score[i] = 0;
    }

    // Compute score per residue using the matrices computed above
    for(int i = 0; i < alignLength; i++) {
        for(int j = 0; j < alignLength; j++) {
            if(i == j) continue; // no self-interaction
            score[i][j] = 0.25 * ((dist_l1[i][j] < 0.5) +
                                    (dist_l1[i][j] < 1.0) +
                                    (dist_l1[i][j] < 2.0) +
                                    (dist_l1[i][j] < 4.0));
        }
    }

    // Normalize
    for(int col = 0; col < queryLength; col++) {
        for(int row = 0; row < queryLength; row++) {
            norm[col] += dists_to_score[row][col];
        }
        if(norm[col] != 0) norm[col] = 1 / norm[col];
        else norm[col] = INF;
    }
    for(int idx = 0; idx < alignLength; idx++) {
        norm_aligned[idx] = norm[align_to_query[idx]];
    }

    for(int i = 0; i < alignLength; i++) {
        for(int j = 0; j < alignLength; j++) {
            score[i][j] *= aligned_dists_to_score[i][j];
            reduce_score[j] += score[i][j];
        }
    }
    for(int idx = 0; idx < alignLength; idx++) {
        reduce_score[idx] *= norm_aligned[idx]; // reduce_score[] contains the lddt scores per residue
    }

    delete[] norm; delete[] norm_aligned; 
    for(int i = 0; i < alignLength; i++) {
        delete[] score[i];
    }
    delete[] score;
}