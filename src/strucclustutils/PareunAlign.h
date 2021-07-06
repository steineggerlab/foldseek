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

#endif //STRUCCLUST_PAREUNALIGN_H



