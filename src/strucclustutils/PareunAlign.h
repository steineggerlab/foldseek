//
// Created by Martin Steinegger on 1/15/21.
//
//
// Created by Stephanie Kim on 4/28/21.
//

#ifndef STRUCCLUST_PAREUNALIGN_H
#define STRUCCLUST_PAREUNALIGN_H

#include "Alignment.h"
#include "StripedSmithWaterman.h"
#include "tmalign/Coordinates.h"
#include "LocalParameters.h"
#include <vector>
#include <cstring>
#include "structureto3diseqdist.h"
using namespace std;

class PareunAlign  {
private:
    const SubstitutionMatrix * mat;
    short *matrix;

public:

    PareunAlign(unsigned int maxSeqLen, SubstitutionMatrix * subMat): mat(subMat){
        //allocate memory
        matrix = (short*)mem_align(ALIGN_FLOAT,maxSeqLen*maxSeqLen*sizeof(short ));

    }

    ~PareunAlign(){
        //TODO delete things
        free(matrix);
    }
    EvalueComputation * evaluer;

    Matcher::result_t align(Sequence & qSeq, Sequence & tSeq, SubstitutionMatrix * subMat, EvalueComputation evaluer);

    void AlignedResidueIndex(Matcher::result_t  & optAlnResult, int * ires);

    string backtrace2cigar (string backtrace);
};
// todo: clean up
const int DIST_MAT_SIZE = Alphabet3diSeqDist::CENTROID_CNT + 1;
//    A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   X
const short distMat[DIST_MAT_SIZE][DIST_MAT_SIZE] = { {14,5, -19, -50, -50, -50, -50,-50, -50 ,-50, -50, -50, -50, -50, -50, -50, -50, -50, -50,-50, 0},
                              {5, 11,  2, -11, -15, -18, -23, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, 0},
                              {-19,   2,   9,   2,  -7, -13, -16, -23, -50, -24, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                              {-50, -11,   2,   7,   2,  -6, -10, -14, -18, -22, -27, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                              {-50, -15,  -7,   2,   7,   2,  -5,  -9, -12, -17, -22, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                              {-50, -18, -13,  -6,   2,   7,   2,  -4,  -9, -14, -19, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                              {-50, -23, -16, -10,  -5,   2,   8,   4,  -3,  -7, -13, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                              {-50, -50, -23, -14,  -9,  -4,   4,   9,   1,  -3, -10, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                              {-50, -50, -50, -18, -12,  -9,  -3,   1,  11,   2,  -6, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                              {-50, -50, -24, -22, -17, -14,  -7,  -3,   2,  10,  -5, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                              {-50, -50, -50, -27, -22, -19, -13, -10,  -6,  -5,   6, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                              {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50,   8,  -4,  -8, -13, -18, -26, -50, -50, -50,   0},
                              {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50,  -4,   9,   0,  -6, -12, -16, -22, -50, -50,   0},
                              {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50,  -8,   0,   9,   3,  -5, -12, -17, -23, -50,   0},
                              {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -13,  -6,   3,   7,   1,  -8, -12, -19, -50,   0},
                              {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -18, -12,  -5,   1,   6,   1,  -6, -15, -50,   0},
                              {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -26, -16, -12,  -8,   1,   6,   2, -10, -50,   0},
                              {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -22, -17, -12,  -6,   2,   8,   2, -50,   0},
                              {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -23, -19, -15, -10,   2,  10,   4,   0},
                              {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50,   4,  12,   0},
                              { 0,   0,   0,   0,   0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,   0,   0,   0}};

#endif //STRUCCLUST_PAREUNALIGN_H



