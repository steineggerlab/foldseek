#include <iostream>
#include <unordered_map>
#include <cmath>
#include <map>

#ifndef LDDT_H
#define LDDT_H


typedef float* fptr_t;
typedef float** farrptr_t;

class LDDTcalculator {
public:
    LDDTcalculator(unsigned int maxQueryLength, unsigned int maxTargetLength);
    ~LDDTcalculator();

    struct Grid {
        Grid() {};
        Grid(farrptr_t& m1, unsigned int queryLength) {
            int len = queryLength;
            for(int i = 0; i < len; i++) {
                for(int dim = 0; dim < 3; dim++) {
                    if(m1[i][dim] < min[dim]) min[dim] = m1[i][dim];
                    if(m1[i][dim] > max[dim]) max[dim] = m1[i][dim];
                }
            }
            for(int i = 0; i < (int)queryLength; i++) {
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
        std::multimap<std::tuple<int, int, int>, int> box; // bottleneck
    };

    struct LDDTscoreResult {
        LDDTscoreResult() {
            avgLddtScore = 0.0;
            scoreLength = 0;
        }
        LDDTscoreResult(float *reduce_score, int alignLength) {
            scoreLength = alignLength;
            if(perCaLddtScore) {
                delete[] perCaLddtScore;
            }
            perCaLddtScore = new float[scoreLength];
            float sum = 0.0;
            for(int i = 0; i < scoreLength; i++) {
                sum += reduce_score[i];
                perCaLddtScore[i] = reduce_score[i];
            }
            avgLddtScore = (double)(sum/(float)scoreLength);
        }
        LDDTscoreResult(const LDDTscoreResult& r1) {
            scoreLength = r1.scoreLength;
            avgLddtScore = r1.avgLddtScore;
            perCaLddtScore = new float[scoreLength];
            for(int i = 0; i < scoreLength; i++) {
                perCaLddtScore[i] = r1.perCaLddtScore[i];
            }
        }
        LDDTscoreResult& operator= (const LDDTscoreResult& r1) {
            if(this == &r1) return *this;
            if(perCaLddtScore) {
                delete[] perCaLddtScore;
            }
            scoreLength = r1.scoreLength;
            avgLddtScore = r1.avgLddtScore;
            perCaLddtScore = new float[scoreLength];
            for(int i = 0; i < scoreLength; i++) {
                perCaLddtScore[i] = r1.perCaLddtScore[i];
            }
            return *this;
        }
        ~LDDTscoreResult() {
            if(perCaLddtScore) {
                delete[] perCaLddtScore;
            }
        }

        fptr_t perCaLddtScore = NULL;
        int scoreLength;
        double avgLddtScore;
    };

    void initQuery(unsigned int queryLen, float *qx, float *qy, float *qz);
    void construct_hash_tables_align(int align_idx, int query_idx, int target_idx);
    void calculate_distance();
    void compute_scores();
    LDDTscoreResult computeLDDTScore(unsigned int targetLen, int qStartPos, int tStartPos, const std::string &backtrace, float *tx, float *ty, float *tz);

private:
    static const float cutoff;
    static const float INF;

    unsigned int queryStart, targetStart, queryLength, targetLength, alignLength;
    unsigned int maxQueryLength, maxTargetLength, maxAlignLength;
    fptr_t reduce_score, norm, norm_aligned;
    std::unordered_map<int, int> query_to_align, target_to_align, align_to_query, align_to_target;
    std::string cigar; // backtrace
    farrptr_t query_pos, target_pos, dists_to_score, aligned_dists_to_score, dist_l1, score;
    LDDTcalculator::Grid query_grid;
};

#endif