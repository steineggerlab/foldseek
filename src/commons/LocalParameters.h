#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

struct FoldSeekDbValidator : public DbValidator {
    static std::vector<int> tmscore;
    static std::vector<int> cadb;
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
    static const int DBTYPE_TMSCORE;
    std::vector<MMseqsParameter *> strucclust;
    std::vector<MMseqsParameter *> tmalign;
    std::vector<MMseqsParameter *> structuresearchworkflow;
    PARAMETER(PARAM_TMSCORE_THRESHOLD)
    float tmScoreThr;
private:

    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
