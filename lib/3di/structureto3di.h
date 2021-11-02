#pragma once
#include <vector>
#include <stddef.h>
#include <cstring>
#include "kerasify/keras_model.h"

namespace Alphabet3Di{
    static const size_t CENTROID_CNT = 20;
    static const char INVALID_STATE = 2; // assign invalid residues to coil state
    const double DISTANCE_ALPHA_BETA = 1.5336;
    const double PI = 3.14159265359;
    static const size_t FEATURE_CNT = 10;
    static const size_t EMBEDDING_DIM = 2;
    static const struct {
        double alpha, beta, d;
    } VIRTUAL_CENTER = { 270, 0, 2 };

    const double centroids[CENTROID_CNT][EMBEDDING_DIM] = {
        { -1.0729,  -0.3600},
        { -0.1356,  -1.8914},
        {  0.4948,  -0.4205},
        { -0.9874,   0.8128},
        { -1.6621,  -0.4259},
        {  2.1394,   0.0486},
        {  1.5558,  -0.1503},
        {  2.9179,   1.1437},
        { -2.8814,   0.9956},
        { -1.1400,  -2.0068},
        {  3.2025,   1.7356},
        {  1.7769,  -1.3037},
        {  0.6901,  -1.2554},
        { -1.1061,  -1.3397},
        {  2.1495,  -0.8030},
        {  2.3060,  -1.4988},
        {  2.5522,   0.6046},
        {  0.7786,  -2.1660},
        { -2.3030,   0.3813},
        {  1.0290,   0.8772},
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

    // find closest member for every c beta atom
    void findResiduePartners(std::vector<int> & partnerIdx, Vec3 * cb,
                             std::vector<bool> & validMask, const size_t len);

};

class StructureTo3Di : StructureTo3DiBase{
public:

    StructureTo3Di();
    ~StructureTo3Di(){};
    char * structure2states(Vec3 * ca, Vec3 * n,
                            Vec3 * c, Vec3 * cb,
                            size_t len);

protected:

    struct Embedding{
        double f[Alphabet3Di::EMBEDDING_DIM];
        Embedding(){}
        Embedding(double param_f[Alphabet3Di::EMBEDDING_DIM]){
            memcpy(f, param_f, sizeof(double) * Alphabet3Di::EMBEDDING_DIM);
        }
    };

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

    // Describe interaction of residue i and j
    Feature calcFeatures(Vec3 * ca, int i, int j);

    void calcConformationDescriptors(std::vector<Feature> & features, std::vector<int> & partnerIdx,
                                     Vec3 * ca, std::vector<bool> & mask, const size_t len);

    void encodeFeatures(std::vector<Embedding> & embeddings, std::vector<Feature> & features,
                                        std::vector<bool> & mask, const size_t len);

    void discretizeEmbeddings(std::vector<char> & states, std::vector<Embedding> & embeddings,
                            std::vector<bool> & mask, const size_t len);
};


