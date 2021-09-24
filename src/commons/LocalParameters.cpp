#include "LocalParameters.h"
#include "Command.h"

#include "mat3di.out.h"


const int LocalParameters::DBTYPE_CA_ALPHA = 101;
const int LocalParameters::DBTYPE_TMSCORE = 102;

LocalParameters::LocalParameters() :
        Parameters(),
        PARAM_TMSCORE_THRESHOLD(PARAM_TMSCORE_THRESHOLD_ID,"--tmscore-threshold", "TMscore threshold", "accept alignments with a tmsore > thr [0.0,1.0]",typeid(float), (void *) &tmScoreThr, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        PARAM_ALIGNMENT_TYPE(PARAM_ALIGNMENT_TYPE_ID,"--alignment-type", "Alignment type", "How to compute the alignment:\n0: 3di alignment\n1: TM alignment\n2: pareun alignment\n",typeid(int), (void *) &alignmentType, "^[0-2]{1}$")
{
    scoringMatrixFile = "3di.out";
    seedScoringMatrixFile = "3di.out";
    substitutionMatrices.emplace_back("3di.out", mat3di_out, mat3di_out_len);
    // structurecreatedb
    structurecreatedb.push_back(&PARAM_THREADS);
    structurecreatedb.push_back(&PARAM_WRITE_LOOKUP);
    structurecreatedb.push_back(&PARAM_V);
    // tmalign
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
    // strucclust
    strucclust = combineList(clust, align);
    strucclust = combineList(strucclust, kmermatcher);
    strucclust.push_back(&PARAM_REMOVE_TMP_FILES);
    strucclust.push_back(&PARAM_RUNNER);
    // structuresearchworkflow
    structuresearchworkflow = combineList(align, prefilter);
    structuresearchworkflow = combineList(tmalign, structuresearchworkflow);
    structuresearchworkflow.push_back(&PARAM_ALIGNMENT_TYPE);
    structuresearchworkflow.push_back(&PARAM_REMOVE_TMP_FILES);

    // Setup DbValidation
    easystructuresearchworkflow = combineList(structuresearchworkflow, structurecreatedb);

    alignmentType = ALIGNMENT_TYPE_3DI;
    tmScoreThr = 0.5;

}

std::vector<int> FoldSeekDbValidator::tmscore = {LocalParameters::DBTYPE_TMSCORE};
std::vector<int> FoldSeekDbValidator::cadb = {LocalParameters::DBTYPE_CA_ALPHA};
std::vector<int> FoldSeekDbValidator::flatfileStdinAndFolder = {LocalParameters::DBTYPE_FLATFILE, LocalParameters::DBTYPE_STDIN,LocalParameters::DBTYPE_DIRECTORY};
std::vector<int> FoldSeekDbValidator::flatfileAndFolder = {LocalParameters::DBTYPE_FLATFILE, LocalParameters::DBTYPE_DIRECTORY};