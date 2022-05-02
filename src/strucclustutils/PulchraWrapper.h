
#ifndef FOLDSEEK_PUCHRAWRAPPER_H
#define FOLDSEEK_PUCHRAWRAPPER_H
#include <vector>
#include "structureto3di.h" // for Vec3
#include "pulchra.h"

class PulchraWrapper {
public:
    PulchraWrapper();
    ~PulchraWrapper();
    void rebuildBackbone(Vec3 * ca, Vec3 * n, Vec3 * c, char * ami, size_t chainLen);

private:
    // Ca coords array must be longer for pulchra: offset + chain_len + offest
    const size_t offset = 5; 

    std::vector<double*> xcaVec;
    std::vector<double*> nVec;
    std::vector<double*> cVec;
};


#endif //FOLDSEEK_PULCHRAWRAPPER_H
