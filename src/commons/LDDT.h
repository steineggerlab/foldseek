// #include <iostream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <map>
#include <limits>
#include <vector>
#include "FastSort.h"

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
            box.clear();
            box_index.clear();
            for(int i = 0; i < (int)queryLength; i++) {
                int box_coord[3];
                for(int dim = 0; dim < 3; dim++) {
                    box_coord[dim] = (int)((m1[i][dim] - min[dim]) / CUTOFF);
                }
                box.emplace_back(std::make_tuple(box_coord[0], box_coord[1], box_coord[2]), i);
            }
            SORT_SERIAL(box.begin(), box.end());
            // precompute index for box to get range per entry
            // tuple -> (start, end) in box_index
            std::tuple<int, int, int> currBoxKey = box[0].first;
            size_t currBoxStart = 0;
            for(size_t i = 0; i < box.size(); i++){
                if(box[i].first != currBoxKey){
                    box_index.push_back(std::make_pair(currBoxKey, std::make_pair(currBoxStart, i)));
                    currBoxKey = box[i].first;
                    currBoxStart = i;
                }
            }
            box_index.push_back(std::make_pair(currBoxKey, std::make_pair(currBoxStart, box.size())));
            for(int dim = 0; dim < 3; dim++) {
                num_cells[dim] = (int)((max[dim]-min[dim])/CUTOFF) + 1;
            }
        }

        std::tuple<int, int, int> getGridCoordinates(const float *point) {
            int box_coord[3];
            for(int dim = 0; dim < 3; dim++) {
                box_coord[dim] = (int)((point[dim] - min[dim]) / CUTOFF);
            }
            return std::make_tuple(box_coord[0], box_coord[1], box_coord[2]);
        }

        struct PairTupleComparator {
            bool operator()(const std::pair<std::tuple<int, int, int>, std::pair<long unsigned int, long unsigned int>>& lhs,
                            const std::pair<std::tuple<int, int, int>, std::pair<long unsigned int, long unsigned int>>& rhs) const {
                // Compare the tuples first
                if (lhs.first < rhs.first) {
                    return true;
                }
                if (rhs.first < lhs.first) {
                    return false;
                }

                // If the tuples are equal, compare the pairs
                return lhs.second < rhs.second;
            }
        };

        std::pair<size_t, size_t> getBoxMemberRange(std::tuple<int, int, int> &box_coord) {
            auto it = std::lower_bound(box_index.begin(), box_index.end(), std::make_pair(box_coord, std::make_pair(0, 0)), Grid::PairTupleComparator());
            if(it == box_index.end() || it->first != box_coord) {
                return std::make_pair(0, 0);
            }
            return it->second;
        }

        float min[3] = {INF, INF, INF};
        float max[3] = {-INF, -INF, -INF};
        int num_cells[3];
        std::vector<std::pair<std::tuple<int, int, int>, int>> box;
        std::vector<std::pair<std::tuple<int, int, int>, std::pair<size_t, size_t>>> box_index;
    };

    struct LDDTScoreResult {
        LDDTScoreResult() {
            avgLddtScore = 0.0;
            scoreLength = 0;
        }
        LDDTScoreResult(float *reduce_score, int alignLength) {
            if(perCaLddtScore) {
                delete[] perCaLddtScore;
            }
            perCaLddtScore = new float[alignLength];
            float sum = 0.0;
            scoreLength = alignLength;
            for(int i = 0; i < alignLength; i++) {
                if (std::isnan(reduce_score[i])) {
                    scoreLength = scoreLength - 1;
                    perCaLddtScore[i] = 0;
                } else {
                    sum += reduce_score[i];
                    perCaLddtScore[i] = reduce_score[i];
                }
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
    void constructAlignHashes(int query_idx, int target_idx, const std::string & cigar);
    void calculateLddtScores();
    LDDTScoreResult computeLDDTScore(unsigned int targetLen, int qStartPos, int tStartPos, const std::string &backtrace, float *tx, float *ty, float *tz);

private:
    unsigned int queryLength, targetLength, alignLength;
    unsigned int maxQueryLength, maxTargetLength, maxAlignLength;
    float *reduce_score, *norm;
    int * query_to_align;
    int * target_to_align;
    int * align_to_query;
    int * align_to_target;
    float **query_coordinates, **target_coordinates, **score;
    bool **dists_to_score;
    LDDTCalculator::Grid query_grid;
};

#endif
