#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

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

    static const int DBTYPE_TOSION_SEQUENCE = 20;
    static const int DBTYPE_CA_ALPHA = 21;

    std::vector<MMseqsParameter *> strucclust;
    std::vector<MMseqsParameter *> tmalign;

    PARAMETER(PARAM_TMSCORE_THRESHOLD)
    float tmScoreThr;
private:

    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
