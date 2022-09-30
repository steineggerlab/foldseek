#include <iostream>
#include <unordered_map>
#include <cmath>
#include <map>

#ifndef LDDT_H
#define LDDT_H

#define INF 3.40282e+038 // TODO: More elegant way of handling divide-by-zero
const float cutoff = 15.0;
typedef float** table_t;

class LDDTcalculator {
public:
    LDDTcalculator() {}
    LDDTcalculator(int qLen, int tLen) {
        queryLength = qLen, targetLength = tLen;
    }
    ~LDDTcalculator() {
        if(reduce_score) delete[] reduce_score;
        if(query_pos) {
            for(int i = 0; i < queryLength; i++) {
                delete[] query_pos[i];
            }
            delete[] query_pos;
        }
        if(target_pos) {
            for(int i = 0; i < targetLength; i++) {
                delete[] target_pos[i];
            }
            delete[] target_pos;
        }
        if(dists_to_score) {
            for(int i = 0; i < queryLength; i++) {
                delete[] dists_to_score[i];
            }
            delete[] dists_to_score;
        }
        if(aligned_dists_to_score) {
            for(int i = 0; i < alignLength; i++) {
                delete[] aligned_dists_to_score[i];
            }
            delete[] aligned_dists_to_score;
        }
        if(dist_l1) {
            for(int i = 0; i < alignLength; i++) {
                delete[] dist_l1[i];
            }
            delete[] dist_l1;
        }
    }

    struct Grid {
        Grid() {};
        Grid(table_t& m1, int queryLength) {
            for(int i = 0; i < queryLength; i++) {
                for(int dim = 0; dim < 3; dim++) {
                    if(m1[i][dim] < min[dim]) min[dim] = m1[i][dim];
                    if(m1[i][dim] > max[dim]) max[dim] = m1[i][dim];
                }
            }
            for(int i = 0; i < queryLength; i++) {
                int box_coord[3];
                for(int dim = 0; dim < 3; dim++) {
                    box_coord[dim] = (int)((m1[i][dim] - min[dim]) / cutoff);
                }
                box.insert(std::make_pair(std::make_tuple(box_coord[0], box_coord[1], box_coord[2]), i));
            }

            for(int dim = 0; dim < 3; dim++) {
                num_cells[dim] = (int)((max[dim]-min[dim])/cutoff) + 1;
            }
        }

        float min[3] = {INF, INF, INF};
        float max[3] = {-INF, -INF, -INF};
        int num_cells[3];
        std::multimap<std::tuple<int, int, int>, int> box;
    };

    struct LDDTscoreResult {
        LDDTscoreResult() {}
        LDDTscoreResult(float *reduce_score, int alignLength) {
            scoreLength = alignLength;
            if(perCaLddtScore) delete[] perCaLddtScore;
            perCaLddtScore = new float[scoreLength];
            float sum = 0.0;
            for(int i = 0; i < scoreLength; i++) {
                sum += reduce_score[i];
                perCaLddtScore[i] = reduce_score[i];
            }
            avgLddtScore = (double)(sum/(float)scoreLength);
        }

        float *perCaLddtScore = NULL;
        int scoreLength;
        double avgLddtScore;
    };

    void initVariables(unsigned int queryLen, unsigned int targetLen, int qStartPos, int tStartPos, const std::string &backtrace);
    void construct_hash_tables_align(int align_idx, int query_idx, int target_idx, int cigar_idx);
    void calculate_distance();
    void compute_scores();
    LDDTscoreResult computeLDDTScore(float *qx, float *qy, float *qz, float *tx, float *ty, float *tz, int qStartPos, int tStartPos);

private:
    int queryStart, targetStart, queryLength, targetLength, alignLength;
    float *reduce_score = NULL;
    std::unordered_map<int, int> query_to_align, target_to_align, align_to_query, align_to_target;
    std::string cigar; // backtrace
    table_t query_pos = NULL, target_pos = NULL, dists_to_score = NULL, aligned_dists_to_score = NULL, dist_l1 = NULL;
};

#endif