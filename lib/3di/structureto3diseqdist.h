#pragma once
#include <vector>
#include <stddef.h>
#include <cstring>
#include "structureto3di.h"

namespace Alphabet3diSeqDist{
    static const size_t CENTROID_CNT = 20;
    static const char INVALID_STATE = CENTROID_CNT;
    const int centroids[CENTROID_CNT] =
        {-284,-147,-83,-52,-33,-21,-13,-7,-4,-3,-1,1,3,7,13,24,40,68,123,250};
}

class StructureTo3diSeqDist : StructureTo3DiBase {
public:

    StructureTo3diSeqDist();
    ~StructureTo3diSeqDist(){};
    char * structure2states(Vec3 * ca, Vec3 * n,
                            Vec3 * c, Vec3 * cb,
                            size_t len);

private:

    // store for the class
    std::vector<char> states;
    std::vector<int> partnerIdx;
    std::vector<bool> mask;

    void discretizeSeqDistance(std::vector<char> & states, std::vector<int> & partnerIdx,
                            std::vector<bool> & mask, const size_t len);
};


