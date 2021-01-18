//
// Created by Martin Steinegger on 1/15/21.
//

#ifndef STRUCCLUST_PAREUNALIGN_H
#define STRUCCLUST_PAREUNALIGN_H

#include "Alignment.h"
#include "StripedSmithWaterman.h"
#include "tmalign/Coordinates.h"

class PareunAlign  {
private:
    const SubstitutionMatrix * mat;
//    short *matrix;

public:

    PareunAlign(unsigned int maxSeqLen, SubstitutionMatrix * subMat):mat(subMat){
        //allocate memory
//      matrix = (short*)mem_align(ALIGN_FLOAT,par.maxSeqLen*par.maxSeqLen*sizeof(short ));

    }

    ~PareunAlign(){
        //TODO delete things
//        free(matrix);
    }

    s_align align(char * querySeq, char *  targetSeq, Coordinates * queryCaCords, Coordinates * targetCaCords);
};


#endif //STRUCCLUST_PAREUNALIGN_H
