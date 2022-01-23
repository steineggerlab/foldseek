#include "LocalParameters.h"
#include "Command.h"

#include "mat3di.out.h"


const int LocalParameters::DBTYPE_CA_ALPHA = 101;
const int LocalParameters::DBTYPE_TMSCORE = 102;

LocalParameters::LocalParameters() :
        Parameters(),
        PARAM_TMSCORE_THRESHOLD(PARAM_TMSCORE_THRESHOLD_ID,"--tmscore-threshold", "TMscore threshold", "accept alignments with a tmsore > thr [0.0,1.0]",typeid(float), (void *) &tmScoreThr, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        PARAM_ALIGNMENT_TYPE(PARAM_ALIGNMENT_TYPE_ID,"--alignment-type", "Alignment type", "How to compute the alignment:\n0: 3di alignment\n1: TM alignment\n2: 3Di+AA\n",typeid(int), (void *) &alignmentType, "^[0-2]{1}$"),
        PARAM_CHAIN_NAME_MODE(PARAM_CHAIN_NAME_MODE_ID,"--chain-name-mode", "Chain name mode", "Add chain to name:\n0: auto\n1: always add\n",typeid(int), (void *) &chainNameMode, "^[0-1]{1}$"),
        PARAM_TMALIGN_FAST(PARAM_TMALIGN_FAST_ID,"--tmalign-fast", "TMalign fast","turn on fast search in TM-align" ,typeid(int), (void *) &tmAlignFast, "^[0-1]{1}$"),
        PARAM_N_SAMPLE(PARAM_N_SAMPLE_ID, "--n-sample", "Sample size","pick N random sample" ,typeid(int), (void *) &nsample, "^[0-9]{1}[0-9]*$")
{
    scoringMatrixFile = "3di.out";
    seedScoringMatrixFile = "3di.out";
    substitutionMatrices.emplace_back("3di.out", mat3di_out, mat3di_out_len);
    // structurecreatedb
    structurecreatedb.push_back(&PARAM_CHAIN_NAME_MODE);
    structurecreatedb.push_back(&PARAM_WRITE_LOOKUP);
    structurecreatedb.push_back(&PARAM_THREADS);
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
    tmalign.push_back(&PARAM_TMALIGN_FAST);
    tmalign.push_back(&PARAM_PRELOAD_MODE);
    tmalign.push_back(&PARAM_THREADS);
    tmalign.push_back(&PARAM_V);
//    tmalign.push_back(&PARAM_GAP_OPEN);
//    tmalign.push_back(&PARAM_GAP_EXTEND);
    // strucclust
    strucclust = combineList(clust, align);
    strucclust = combineList(strucclust, kmermatcher);
    strucclust.push_back(&PARAM_REMOVE_TMP_FILES);
    strucclust.push_back(&PARAM_RUNNER);
    // structuresearchworkflow
    // structuresearchworkflow
    structuresearchworkflow = combineList(align, prefilter);
    structuresearchworkflow = combineList(tmalign, structuresearchworkflow);
    structuresearchworkflow.push_back(&PARAM_ALIGNMENT_TYPE);
    structuresearchworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
    structuresearchworkflow.push_back(&PARAM_RUNNER);
    structuresearchworkflow.push_back(&PARAM_REUSELATEST);

    // Setup DbValidation
    easystructuresearchworkflow = combineList(structuresearchworkflow, structurecreatedb);
    easystructuresearchworkflow = combineList(easystructuresearchworkflow, convertalignments);

    // Setup DbValidation

    structureclusterworkflow = combineList(prefilter, align);
    structureclusterworkflow = combineList(structureclusterworkflow, rescorediagonal);
    structureclusterworkflow = combineList(structureclusterworkflow, tmalign);
    structureclusterworkflow = combineList(structureclusterworkflow, clust);
    structureclusterworkflow.push_back(&PARAM_CASCADED);
    structureclusterworkflow.push_back(&PARAM_CLUSTER_STEPS);
    //structuresearchworkflow.push_back(&PARAM_CLUSTER_REASSIGN);
    structureclusterworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
    structureclusterworkflow.push_back(&PARAM_REUSELATEST);
    structureclusterworkflow.push_back(&PARAM_RUNNER);
    structureclusterworkflow = combineList(structureclusterworkflow, linclustworkflow);

    databases.push_back(&PARAM_HELP);
    databases.push_back(&PARAM_HELP_LONG);
    databases.push_back(&PARAM_REUSELATEST);
    databases.push_back(&PARAM_REMOVE_TMP_FILES);
    databases.push_back(&PARAM_COMPRESSED);
    databases.push_back(&PARAM_THREADS);
    databases.push_back(&PARAM_V);
    //easystructureclusterworkflow = combineList(structuresearchworkflow, structurecreatedb);
    samplemulambda.push_back(&PARAM_N_SAMPLE);
    samplemulambda.push_back(&PARAM_THREADS);
    samplemulambda.push_back(&PARAM_V);

    alignmentType = ALIGNMENT_TYPE_3DI_AA;
    tmScoreThr = 0.5;
    chainNameMode = 0;
    tmAlignFast = 0;
    nsample = 5000;
}

std::vector<int> FoldSeekDbValidator::tmscore = {LocalParameters::DBTYPE_TMSCORE};
std::vector<int> FoldSeekDbValidator::cadb = {LocalParameters::DBTYPE_CA_ALPHA};
std::vector<int> FoldSeekDbValidator::flatfileStdinAndFolder = {LocalParameters::DBTYPE_FLATFILE, LocalParameters::DBTYPE_STDIN,LocalParameters::DBTYPE_DIRECTORY};
std::vector<int> FoldSeekDbValidator::flatfileAndFolder = {LocalParameters::DBTYPE_FLATFILE, LocalParameters::DBTYPE_DIRECTORY};