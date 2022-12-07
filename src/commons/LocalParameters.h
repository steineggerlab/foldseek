#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

const int CITATION_FOLDSEEK = CITATION_END;

struct FoldSeekDbValidator : public DbValidator {
    static std::vector<int> tmscore;
    static std::vector<int> cadb;
    static std::vector<int> flatfileStdinAndFolder;
    static std::vector<int> flatfileAndFolder;

};

class LocalParameters : public Parameters {
public:
    static void initInstance() {
        new LocalParameters;
    }
    LocalParameters();
    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    static const int DBTYPE_CA_ALPHA;
    static const int DBTYPE_CA_ALPHA_F16;
    static const int DBTYPE_TMSCORE;

    static const int ALIGNMENT_TYPE_3DI = 0;
    static const int ALIGNMENT_TYPE_TMALIGN = 1;
    static const int ALIGNMENT_TYPE_3DI_AA = 2;

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

    static const int COORD_STORE_MODE_CA_FLOAT = 1;
    static const int COORD_STORE_MODE_CA_HALF  = 2;

    static const unsigned int INDEX_DB_CA_KEY = 500;

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
    PARAMETER(PARAM_TMSCORE_THRESHOLD)
    PARAMETER(PARAM_LDDT_THRESHOLD)
    PARAMETER(PARAM_SORT_BY_STRUCTURE_BITS)
    PARAMETER(PARAM_MASK_BFACTOR_THRESHOLD)
    PARAMETER(PARAM_ALIGNMENT_TYPE)
    PARAMETER(PARAM_CHAIN_NAME_MODE)
    PARAMETER(PARAM_TMALIGN_FAST)
    PARAMETER(PARAM_N_SAMPLE)
    PARAMETER(PARAM_COORD_STORE_MODE)

    float tmScoreThr;
    float lddtThr;
    int sortByStructureBits;
    float maskBfactorThreshold;
    int alignmentType;
    int chainNameMode;
    int tmAlignFast;
    int nsample;
    int coordStoreMode;

    static std::vector<int> getOutputFormat(int formatMode, const std::string &outformat, bool &needSequences, bool &needBacktrace, bool &needFullHeaders,
                                            bool &needLookup, bool &needSource, bool &needTaxonomyMapping, bool &needTaxonomy, bool &needCa, bool &needTMaligner, bool &needLDDT);


private:

    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);

};
#endif
