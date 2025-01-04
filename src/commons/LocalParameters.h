#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

const int CITATION_FOLDSEEK = CITATION_END;
const int CITATION_FOLDSEEK_MULTIMER = CITATION_FOLDSEEK << 1;
const int CITATION_PROSTT5 = CITATION_FOLDSEEK << 2;

struct FoldSeekDbValidator : public DbValidator {
    static std::vector<int> tmscore;
    static std::vector<int> cadb;
    static std::vector<int> flatfileStdinAndFolder;
    static std::vector<int> flatfileAndFolder;

};

class LocalParameters : public Parameters {
public:
    LocalParameters();
    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initParameterSingleton();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    static const int DBTYPE_CA_ALPHA;
    static const int DBTYPE_TMSCORE;

    static const int ALIGNMENT_TYPE_3DI = 0;
    static const int ALIGNMENT_TYPE_TMALIGN = 1;
    static const int ALIGNMENT_TYPE_3DI_AA = 2;

    static const int TMSCORE_THRESHOLD_MODE_ALIGNMENT = 0;
    static const int TMSCORE_THRESHOLD_MODE_QUERY = 1;
    static const int TMSCORE_THRESHOLD_MODE_TARGET = 2;
    static const int TMSCORE_THRESHOLD_MODE_MIN = 3;

    static const int PREF_MODE_KMER = 0;
    static const int PREF_MODE_UNGAPPED = 1;
    static const int PREF_MODE_EXHAUSTIVE = 2;

    static const int TMALIGN_HIT_ORDER_AVG = 0;
    static const int TMALIGN_HIT_ORDER_QUERY = 1;
    static const int TMALIGN_HIT_ORDER_TARGET = 2;
    static const int TMALIGN_HIT_ORDER_MIN = 3;
    static const int TMALIGN_HIT_ORDER_MAX = 4;

    static const int CHAIN_MODE_AUTO = 0;
    static const int CHAIN_MODE_ADD = 1;

    static const int OUTFMT_QCA = 40;
    static const int OUTFMT_TCA = 41;
    static const int OUTFMT_U = 42;
    static const int OUTFMT_T = 43;
    static const int OUTFMT_ALNTMSCORE = 44;
    static const int OUTFMT_LDDT = 45;
    static const int OUTFMT_LDDT_FULL = 46;
    static const int OUTFMT_RMSD = 47;
    static const int OUTFMT_PROBTP = 48;
    static const int OUTFMT_QTMSCORE = 49;
    static const int OUTFMT_TTMSCORE = 50;
    // for Foldseek-MM
    static const int OUTFMT_QUERY_COMPLEX = 51;
    static const int OUTFMT_TARGET_COMPLEX = 52;
    static const int OUTFMT_Q_COMPLEX_TMSCORE = 53;
    static const int OUTFMT_T_COMPLEX_TMSCORE = 54;
    static const int OUTFMT_ASSIGN_ID = 55;
    static const int OUTFMT_COMPLEX_U = 56;
    static const int OUTFMT_COMPLEX_T = 57;

    static const int DB_EXTRACT_MODE_CHAIN = 0;
    static const int DB_EXTRACT_MODE_INTERFACE = 1;

    static const int COORD_STORE_MODE_CA_FLOAT = 1;
    static const int COORD_STORE_MODE_CA_DIFF  = 2;
    static const int COORD_STORE_MODE_CA_PLAIN_TEXT  = 3;

    static const unsigned int INDEX_DB_CA_KEY_DB1 = 500;
    static const unsigned int INDEX_DB_CA_KEY_DB2 = 502;

    static const int INDEX_EXCLUDE_NONE = 0;
    static const int INDEX_EXCLUDE_KMER_INDEX = 1 << 0;
    static const int INDEX_EXCLUDE_CA = 1 << 1;

    // convert2pdb
    static const int PDB_OUTPUT_MODE_MULTIMODEL = 0;
    static const int PDB_OUTPUT_MODE_SINGLECHAIN = 1;
    static const int PDB_OUTPUT_MODE_COMPLEX = 2;

    // filter mode
    // static const int FILTER_MODE_INTERFACE  = 0;
    // static const int FILTER_MODE_CONFORMATION = 1;
    // static const int FILTER_MODE_LOOSE = 2;

    // TODO
    static const unsigned int FORMAT_ALIGNMENT_PDB_SUPERPOSED = 5;
    std::vector<MMseqsParameter *> strucclust;
    std::vector<MMseqsParameter *> tmalign;
    std::vector<MMseqsParameter *> structurealign;
    std::vector<MMseqsParameter *> structurerescorediagonal;
    std::vector<MMseqsParameter *> structuresearchworkflow;
    std::vector<MMseqsParameter *> structureclusterworkflow;
    std::vector<MMseqsParameter *> databases;
    std::vector<MMseqsParameter *> samplemulambda;
    std::vector<MMseqsParameter *> easystructuresearchworkflow;
    std::vector<MMseqsParameter *> easystructureclusterworkflow;
    std::vector<MMseqsParameter *> structurecreatedb;
    std::vector<MMseqsParameter *> compressca;
    std::vector<MMseqsParameter *> scoremultimer;
    std::vector<MMseqsParameter *> filtermultimer;
    std::vector<MMseqsParameter *> multimerclusterworkflow;
    std::vector<MMseqsParameter *> easymultimerclusterworkflow;
    std::vector<MMseqsParameter *> multimersearchworkflow;
    std::vector<MMseqsParameter *> easymultimersearchworkflow;
    std::vector<MMseqsParameter *> createmultimerreport;
    std::vector<MMseqsParameter *> expandmultimer;
    std::vector<MMseqsParameter *> convert2pdb;

    PARAMETER(PARAM_PREF_MODE)
    PARAMETER(PARAM_TMSCORE_THRESHOLD)
    PARAMETER(PARAM_TMSCORE_THRESHOLD_MODE)
    PARAMETER(PARAM_TMALIGN_HIT_ORDER)
    PARAMETER(PARAM_LDDT_THRESHOLD)
    PARAMETER(PARAM_SORT_BY_STRUCTURE_BITS)
    PARAMETER(PARAM_MASK_BFACTOR_THRESHOLD)
    PARAMETER(PARAM_ALIGNMENT_TYPE)
    PARAMETER(PARAM_CHAIN_NAME_MODE)
    PARAMETER(PARAM_WRITE_MAPPING)
    PARAMETER(PARAM_TMALIGN_FAST)
    PARAMETER(PARAM_EXACT_TMSCORE)
    PARAMETER(PARAM_N_SAMPLE)
    PARAMETER(PARAM_COORD_STORE_MODE)
    PARAMETER(PARAM_MIN_ASSIGNED_CHAINS_THRESHOLD)
    PARAMETER(PARAM_MONOMER_INCLUDE_MODE)
    PARAMETER(PARAM_CLUSTER_SEARCH)
    PARAMETER(PARAM_FILE_INCLUDE)
    PARAMETER(PARAM_FILE_EXCLUDE)
    PARAMETER(PARAM_INDEX_EXCLUDE)
    PARAMETER(PARAM_MULTIMER_REPORT_MODE)
    PARAMETER(PARAM_MULTIMER_REPORT_MODE_BC_COMPAT)
    PARAMETER(PARAM_EXPAND_MULTIMER_EVALUE)
    PARAMETER(PARAM_EXPAND_MULTIMER_EVALUE_BC_COMPAT)
    PARAMETER(PARAM_INPUT_FORMAT)
    PARAMETER(PARAM_PDB_OUTPUT_MODE)
    PARAMETER(PARAM_PROSTT5_MODEL)
    PARAMETER(PARAM_GPU)
    PARAMETER(PARAM_DB_EXTRACTION_MODE)
    PARAMETER(PARAM_DISTANCE_THRESHOLD)
    PARAMETER(PARAM_MULTIMER_TM_THRESHOLD)
    PARAMETER(PARAM_CHAIN_TM_THRESHOLD)
    PARAMETER(PARAM_INTERFACE_LDDT_THRESHOLD)

    int prefMode;
    float tmScoreThr;
    int tmScoreThrMode;
    int tmAlignHitOrder;
    float lddtThr;
    int sortByStructureBits;
    float maskBfactorThreshold;
    int alignmentType;
    int chainNameMode;
    bool writeMapping;
    int tmAlignFast;
    int exactTMscore;
    int nsample;
    int coordStoreMode;
    float minAssignedChainsThreshold;
    int monomerIncludeMode;
    int clusterSearch;
    std::string fileInclude;
    std::string fileExclude;
    int indexExclude;
    int multimerReportMode;
    double eValueThrExpandMultimer;
    int inputFormat;
    int pdbOutputMode;
    float filtMultimerTmThr;
    float filtChainTmThr;
    float filtInterfaceLddtThr;
    std::string prostt5Model;
    int gpu;
    int dbExtractionMode;
    float distanceThreshold;
    int prostt5SplitLength;

    static std::vector<int> getOutputFormat(int formatMode, const std::string &outformat, bool &needSequences, bool &needBacktrace, bool &needFullHeaders,
                                            bool &needLookup, bool &needSource, bool &needTaxonomyMapping, bool &needTaxonomy, bool &needQCa, bool &needTCa, bool &needTMaligner,
                                            bool &needLDDT);


private:
    LocalParameters(LocalParameters const&);
    void operator=(LocalParameters const&);
};
#endif
