#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

class LocalParameters : public Parameters {
public:
    static void initInstance() {
        new LocalParameters;
    }

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
    LocalParameters() :
            Parameters(),
            PARAM_TMSCORE_THRESHOLD(PARAM_TMSCORE_THRESHOLD_ID,"--tmscore-threshold", "TMscore threshold", "accept alignments with a tmsore > thr [0.0,1.0]",typeid(float), (void *) &tmScoreThr, "^0(\\.[0-9]+)?|1(\\.0+)?$")
    {

        tmalign.push_back(&PARAM_MIN_SEQ_ID);
        tmalign.push_back(&PARAM_C);
        tmalign.push_back(&PARAM_COV_MODE);
        tmalign.push_back(&PARAM_INCLUDE_IDENTITY);
        tmalign.push_back(&PARAM_TMSCORE_THRESHOLD);
        tmalign.push_back(&PARAM_THREADS);
        tmalign.push_back(&PARAM_V);

        // strucclust
        strucclust = combineList(clust, align);
        strucclust = combineList(strucclust, kmermatcher);
        strucclust.push_back(&PARAM_REMOVE_TMP_FILES);
        strucclust.push_back(&PARAM_RUNNER);

        tmScoreThr = 0.5;
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
