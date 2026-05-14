#pragma once
#include <vector>
#include <stddef.h>
#include <cstring>
#include "kerasify/keras_model.h"
#include "structureto3di.h"  // Only for Vec3 struct

namespace Alphabet12St{
    static const size_t STATE_CNT = 12;
    static const char INVALID_STATE = 6; // H
    const double DISTANCE_ALPHA_BETA = 1.5336;
    const double PI = 3.14159265359;
    static const size_t FEATURE_CNT = 2;
    static const struct {
        double alpha, beta, d;
    } VIRTUAL_CENTER = { 322, 48, 6.84 };

    // Linear layer parameters for state prediction
    const double LAYER1_BIAS[STATE_CNT] = {
        -22.109913,  -34.763283,  -12.11722,    5.844174,
         27.193422,    7.8374314,  -0.10719311,  2.1946049,
         26.863173,    0.23405483, -27.650888,  -13.840539
    };

    const double LAYER1_W[STATE_CNT][FEATURE_CNT] = {
        { 11.980139,   -3.2107675  },
        { 10.814536,    0.6751775  },
        {  5.8660164,   0.76846164 },
        {  5.8301463,  -2.9750438  },
        { -0.08006129, -2.771098   },
        { -0.11654578,  0.9788566  },
        { -1.4215535,   1.9951683  },
        { -2.1661475,   1.795795   },
        { -2.9870105,  -3.2082853  },
        { -4.2153015,   1.4156913  },
        {-11.160501,    1.2672309  },
        {-12.911917,   -3.346186   }
    };
}

class StructureTo12St : public StructureTo3DiBase{
public:

    StructureTo12St();
    ~StructureTo12St(){};
    char * structure2states(Vec3 * ca, Vec3 * n,
                            Vec3 * c, Vec3 * cb,
                            size_t len);

    struct Feature{
        double f[Alphabet12St::FEATURE_CNT];
        Feature(){
            memset(f, 0, sizeof(double) * Alphabet12St::FEATURE_CNT);
        }
        Feature(double param_f[Alphabet12St::FEATURE_CNT]){
            memcpy(f, param_f, sizeof(double) * Alphabet12St::FEATURE_CNT);
        }
    };

    std::vector<Feature> getFeatures(){
        return features;
    }

private:
    // store for the class
    std::vector<Feature> features;
    std::vector<char> states;
    std::vector<int> partnerIdx;
    std::vector<double> partnerDistances;  // minimum distance to partner for each residue
    std::vector<bool> mask;
    // virtualCenter is inherited from StructureTo3DiBase

    // Override base class methods with 12St-specific implementations
    void replaceCBWithVirtualCenter(Vec3 * ca, Vec3 * n, Vec3 * c, Vec3 * cb, const size_t len);

    // 12St-specific methods
    void findResiduePartners(std::vector<int> & partnerIdx, std::vector<double> & partnerDistances,
                             Vec3 * ca, std::vector<Vec3> & virtualCenter, std::vector<bool> & validMask, const size_t len);

    Feature calcFeatures(Vec3 * ca, int i, int j);

    void calcConformationDescriptors(std::vector<Feature> & features, std::vector<int> & partnerIdx,
                                     Vec3 * ca, std::vector<bool> & mask, const size_t len);

    void predictStates(std::vector<char> & states, std::vector<Feature> & features,
                       std::vector<bool> & mask, const size_t len);
};
