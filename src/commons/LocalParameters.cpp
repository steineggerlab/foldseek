#include "LocalParameters.h"
#include "mat3di.out.h"


LocalParameters::LocalParameters() :
        Parameters(),
        PARAM_TMSCORE_THRESHOLD(PARAM_TMSCORE_THRESHOLD_ID,"--tmscore-threshold", "TMscore threshold", "accept alignments with a tmsore > thr [0.0,1.0]",typeid(float), (void *) &tmScoreThr, "^0(\\.[0-9]+)?|1(\\.0+)?$")
{
    scoringMatrixFile = "3di.out";
    seedScoringMatrixFile = "3di.out";
    substitutionMatrices.emplace_back("3di.out", mat3di_out, mat3di_out_len);

    tmalign.push_back(&PARAM_MIN_SEQ_ID);
    tmalign.push_back(&PARAM_C);
    tmalign.push_back(&PARAM_COV_MODE);
    tmalign.push_back(&PARAM_MAX_REJECTED);
    tmalign.push_back(&PARAM_MAX_ACCEPT);
    tmalign.push_back(&PARAM_ADD_BACKTRACE);
    tmalign.push_back(&PARAM_INCLUDE_IDENTITY);
    tmalign.push_back(&PARAM_TMSCORE_THRESHOLD);
    tmalign.push_back(&PARAM_THREADS);
    tmalign.push_back(&PARAM_V);
    tmalign.push_back(&PARAM_GAP_OPEN);
    tmalign.push_back(&PARAM_GAP_EXTEND);
    // strucclust
    strucclust = combineList(clust, align);
    strucclust = combineList(strucclust, kmermatcher);
    strucclust.push_back(&PARAM_REMOVE_TMP_FILES);
    strucclust.push_back(&PARAM_RUNNER);
    tmScoreThr = 0.5;
}