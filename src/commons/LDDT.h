#include <iostream>
#include <unordered_map>
#include <cmath>
#include <map>

#ifndef LDDT_H
#define LDDT_H

#define INF 3.40282e+038 // TODO: More elegant way of handling divide-by-zero
const float cutoff = 15.0;

class LDDTcalculator {
public:
    LDDTcalculator() {}
    LDDTcalculator(int qLen, int tLen) {
        queryLength = qLen, targetLength = tLen;
    }
    ~LDDTcalculator() {
        if(reduce_score) delete[] reduce_score;
    }

    // TODO: get rid of matrix2D struct
    // Instead, use dynamic float array passed on from structureconvertalis.cpp
    struct matrix2D {
        matrix2D(int dim1 = 0, int dim2 = 0) {
            arr = new float*[dim1]{};
            for(int i = 0; i < dim1; i++) {
                arr[i] = new float[dim2]{};
            }
            xdim = dim1; ydim = dim2;
        }
        ~matrix2D() {
            for(int i = 0; i < xdim; i++) {
                delete[] arr[i];
            }
            delete[] arr;
            arr = nullptr;
            xdim = ydim = 0;
        }
        matrix2D(const matrix2D& m1) { // Copy constructor
            xdim = m1.xdim; ydim = m1.ydim;
            arr = new float*[xdim]{};
            for(int i = 0; i < xdim; i++) {
                arr[i] = new float[ydim]{};
            }
            for(int i = 0; i < xdim; i++) {
                for(int j = 0; j < ydim; j++) {
                    arr[i][j] = m1(i, j);
                }
            }
        } 
        matrix2D& operator= (const matrix2D& m1) { // Assignment operator
            // Self-assignment check
            if(this == &m1) return *this;

            // If data exists, delete it
            if (arr) {
                for(int i = 0; i < xdim; i++) {
                    delete[] arr[i];
                }
                delete[] arr;
                arr = nullptr;
            }

            xdim = m1.xdim; ydim = m1.ydim;
            arr = new float*[xdim]{};
            for(int i = 0; i < xdim; i++) {
                arr[i] = new float[ydim]{};
            }
            for(int i = 0; i < xdim; i++) {
                for(int j = 0; j < ydim; j++) {
                    arr[i][j] = m1(i, j);
                }
            }
            return *this;
        }
        float& operator() (int x_idx, int y_idx) {
            if(x_idx >= xdim || y_idx >= ydim) {
                exit(1); // TODO: Throw error
            }
            return arr[x_idx][y_idx];
        }
        float operator() (int x_idx, int y_idx) const {
            if(x_idx >= xdim || y_idx >= ydim) {
                exit(1); // TODO: Throw error
            }
            return arr[x_idx][y_idx];
        }
        float*& operator() (int x_idx) {
            if(x_idx >= xdim) exit(1);
            return arr[x_idx];
        }
        float* operator() (int x_idx) const {
            if(x_idx >= xdim) exit(1);
            return arr[x_idx];
        }
        bool isCoordinate() const { return ydim == 3; }

        int xdim, ydim;
        float** arr;
    };

    struct Grid {
        Grid() {};
        Grid(matrix2D& m1) {
            int len = m1.xdim;
            for(int i = 0; i < len; i++) {
                for(int dim = 0; dim < 3; dim++) {
                    if(m1(i, dim) < min[dim]) min[dim] = m1(i, dim);
                    if(m1(i, dim) > max[dim]) max[dim] = m1(i, dim);
                }
            }
            for(int i = 0; i < len; i++) {
                int box_coord[3];
                for(int dim = 0; dim < 3; dim++) {
                    box_coord[dim] = (int)((m1(i, dim) - min[dim]) / cutoff);
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

// TODO: encapsulate variables
// private:
    int queryStart, targetStart, queryLength, targetLength, alignLength;
    float *reduce_score = NULL;
    std::unordered_map<int, int> query_to_align, target_to_align, align_to_query, align_to_target;
    std::string cigar; // backtrace
    matrix2D query_pos, target_pos, dists_to_score, aligned_dists_to_score, dist_l1;
};

#endif