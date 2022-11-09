// #include <iostream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <map>
#include <limits>

#ifndef LDDT_H
#define LDDT_H

class LDDTCalculator {
public:
    static constexpr float CUTOFF = 15.0;
    static constexpr float INF = std::numeric_limits<float>::infinity();

    LDDTCalculator(unsigned int maxQueryLength, unsigned int maxTargetLength);
    ~LDDTCalculator();

    struct Grid {
        Grid() {};
        Grid(float **& m1, unsigned int queryLength) {
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
                    box_coord[dim] = (int)((m1[i][dim] - min[dim]) / CUTOFF);
                }
                box.insert(std::make_pair(std::make_tuple(box_coord[0], box_coord[1], box_coord[2]), i));
            }

            for(int dim = 0; dim < 3; dim++) {
                num_cells[dim] = (int)((max[dim]-min[dim])/CUTOFF) + 1;
            }
        }

        float min[3] = {INF, INF, INF};
        float max[3] = {-INF, -INF, -INF};
        int num_cells[3];
        std::multimap<std::tuple<int, int, int>, int> box; // bottleneck
    };

    struct LDDTScoreResult {
        LDDTScoreResult() {
            avgLddtScore = 0.0;
            scoreLength = 0;
        }
        LDDTScoreResult(float *reduce_score, int alignLength) {
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
        LDDTScoreResult(const LDDTScoreResult& r1) {
            scoreLength = r1.scoreLength;
            avgLddtScore = r1.avgLddtScore;
            perCaLddtScore = new float[scoreLength];
            for(int i = 0; i < scoreLength; i++) {
                perCaLddtScore[i] = r1.perCaLddtScore[i];
            }
        }
        LDDTScoreResult& operator= (const LDDTScoreResult& r1) {
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
        ~LDDTScoreResult() {
            if(perCaLddtScore) {
                delete[] perCaLddtScore;
            }
        }

        float *perCaLddtScore = NULL;
        int scoreLength;
        double avgLddtScore;
    };

    void initQuery(unsigned int queryLen, float *qx, float *qy, float *qz);
    void constructAlignHashes(int align_idx, int query_idx, int target_idx);
    void calculateDistance();
    void computeScores();
    LDDTScoreResult computeLDDTScore(unsigned int targetLen, int qStartPos, int tStartPos, const std::string &backtrace, float *tx, float *ty, float *tz);

private:
    unsigned int queryLength, targetLength, alignLength;
    unsigned int maxQueryLength, maxTargetLength, maxAlignLength;
    float *reduce_score, *norm;
    int * query_to_align;
    int * target_to_align;
    int * align_to_query;
    int * align_to_target;
    std::string cigar; // backtrace
    float **query_coordinates, **target_coordinates, **score;
    bool **dists_to_score;
    LDDTCalculator::Grid query_grid;
    float dist(float* arr1, float* arr2);
};

#endif
