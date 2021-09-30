#pragma once
#include <vector>
#include <stddef.h>
#include <cstring>
#include "kerasify/keras_model.h"

namespace Alphabet3Di{
    static const size_t CENTROID_CNT = 16;
    static const char INVALID_STATE = CENTROID_CNT; // 'X'
    const double DISTANCE_ALPHA_BETA = 1.5336;
    const double PI = 3.14159265359;
    static const size_t FEATURE_CNT = 10;
    static const size_t EMBEDDING_DIM = 2;
    static const struct {
        double alpha, beta, d;
    } VIRTUAL_CENTER = { 270, 0, 2 };

    const double centroids[CENTROID_CNT][EMBEDDING_DIM] = {
        { -3.1738,  -0.6845},
        { -2.2047,  -0.6065},
        { -2.2506,  -1.7851},
        {  1.3359,   0.2059},
        {  0.7572,  -0.0786},
        { -0.2827,  -0.0022},
        {  1.0574,   0.9981},
        {  2.0778,   0.6344},
        {  1.3023,  -0.8703},
        {  2.7799,   1.0694},
        { -2.3758,   0.1531},
        {  0.1120,  -1.0275},
        { -0.7017,   0.9251},
        { -4.0351,  -1.0812},
        { -1.1916,  -0.1608},
        { -2.5749,   0.9265},
    };
}

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

class StructureTo3DiBase{
protected:

    struct Feature{
        double f[Alphabet3Di::FEATURE_CNT];
        Feature(){}
        Feature(double param_f[Alphabet3Di::FEATURE_CNT]){
            memcpy(f, param_f, sizeof(double) * Alphabet3Di::FEATURE_CNT);
        }
    };

    struct Embedding{
        double f[Alphabet3Di::EMBEDDING_DIM];
        Embedding(){}
        Embedding(double param_f[Alphabet3Di::EMBEDDING_DIM]){
            memcpy(f, param_f, sizeof(double) * Alphabet3Di::EMBEDDING_DIM);
        }
    };

    Vec3 add(Vec3 a, Vec3 b);
    Vec3 sub(Vec3 a, Vec3 b);
    Vec3 norm(Vec3 a);
    Vec3 cross(Vec3 a, Vec3 b);
    Vec3 scale(Vec3 a, double f);
    double dot(Vec3 a, Vec3 b);

    // Distance in Angstroem between CA and CB
    Vec3 approxCBetaPosition(Vec3 ca_atom, Vec3 n_atom, Vec3 c_atom);

    double degreeToRadians(double degree);

    Vec3 calcVirtualCenter(Vec3 ca, Vec3 cb, Vec3 n, double alpha, double beta, double d);

    double calcDistanceBetween(Vec3 & a, Vec3 & b);

    void replaceCBWithVirtualCenter(Vec3 * ca, Vec3 * n,
                                    Vec3 * c, Vec3 * cb, const size_t len);

    void createResidueMask(std::vector<bool> & validMask, Vec3 * ca, Vec3 * n, Vec3 * c, const size_t len);

    // Describe interaction of residue i and j
    Feature calcFeatures(Vec3 * ca, int i, int j);

};

class StructureTo3Di : StructureTo3DiBase{
public:

    StructureTo3Di();
    ~StructureTo3Di(){};
    char * structure2states(Vec3 * ca, Vec3 * n,
                            Vec3 * c, Vec3 * cb,
                            size_t len);

private:
    // Encoding
    KerasModel encoder;
    Tensor in;
    Tensor out;

    // store for the class
    std::vector<Feature> features;
    std::vector<Embedding> embeddings;
    std::vector<char> states;
    std::vector<int> partnerIdx;
    std::vector<bool> mask;

    // find closest member for every c beta atom
    void findResiduePartners(std::vector<int> & partnerIdx, Vec3 * cb,
                             std::vector<bool> & validMask, const size_t len);

    void calcConformationDescriptors(std::vector<Feature> & features, std::vector<int> & partnerIdx,
                                     Vec3 * ca, std::vector<bool> & mask, const size_t len);

    void encodeFeatures(std::vector<Embedding> & embeddings, std::vector<Feature> & features,
                                        std::vector<bool> & mask, const size_t len);

    void discretizeEmbeddings(std::vector<char> & states, std::vector<Embedding> & embeddings,
                            std::vector<bool> & mask, const size_t len);
};


