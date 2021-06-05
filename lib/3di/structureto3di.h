#pragma once
#include <vector>

namespace Alphabet3Di{
    static const size_t CENTROID_CNT = 16;
    static const char INVALID_STATE = CENTROID_CNT; // 'X'
    const double DISTANCE_ALPHA_BETA = 1.5336;
    const double PI = 3.14159265359;
    static const size_t FEATURE_CNT = 9;

    const double centroids[CENTROID_CNT][FEATURE_CNT] = {
            {1.707370,1.819577,-0.249517,-0.766728,0.884688,0.782864,1.108104,3.567158,1.212435},
            {2.121513,1.686284,1.174320,-1.252905,-0.831693,-0.901668,-1.148154,4.246405,1.211326},
            {0.301294,0.192883,-0.345992,-0.488502,-0.462889,0.702640,1.086615,3.307813,-1.020397},
            {0.159370,0.441201,0.490513,0.391775,0.596406,-0.586385,1.133961,3.289736,1.017771},
            {0.255064,0.223530,1.062004,-1.145884,0.325438,0.248910,-0.730873,4.501290,-1.212687},
            {0.145179,0.271506,-1.473581,-0.120250,1.447031,-0.793825,0.125083,2.393994,-0.306216},
            {1.940468,1.524076,0.811490,1.391718,0.885128,1.393001,0.996144,2.670570,0.232042},
            {0.398424,2.057004,0.863547,-1.062999,-0.412677,-0.230368,-0.482623,4.224249,-1.202583},
            {2.105008,1.939213,-1.208254,-0.854899,1.306241,1.098909,1.024468,2.715030,-0.578725},
            {0.414251,0.308056,0.165866,1.599132,-0.750401,1.460213,0.189711,2.386347,0.305817},
            {0.410397,2.065107,0.940106,-0.872285,-0.182263,-0.551512,-0.226457,4.247471,1.208209},
            {2.076012,0.468851,0.523333,-0.935681,-0.382627,0.208152,0.219874,3.883404,-1.166472},
            {2.066769,1.962041,1.163406,-1.341045,-0.895239,-0.930776,-1.256450,4.317824,-1.212391},
            {0.408314,0.264871,1.148792,-1.267738,0.146937,-0.140438,-1.001876,4.912282,1.212673},
            {0.453411,0.318703,0.670807,-0.993227,0.395643,0.777184,-0.095377,3.729683,1.212699},
            {0.376471,1.950554,-1.403345,-0.848239,1.402843,0.538040,0.914844,2.463481,-0.359628},
    };

    const double feature_scaling[FEATURE_CNT] =
            {3.702326,3.691894,1.479059,1.603245,1.449159,1.462474,1.694236,0.627700,0.303178};
}

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

class StructureTo3Di{
public:

    StructureTo3Di(){};
    ~StructureTo3Di(){};
    char * structure2states(std::vector<Vec3> & ca, std::vector<Vec3> & n,
                            std::vector<Vec3> & c, std::vector<Vec3> & cb);

private:
    struct Feature{
        double f[Alphabet3Di::FEATURE_CNT];
        Feature(){}
        Feature(double param_f[Alphabet3Di::FEATURE_CNT]){
            memcpy(f, param_f, sizeof(double) * Alphabet3Di::FEATURE_CNT);
        }
    };
    // store for the class
    std::vector<Feature> features;
    std::vector<char> states;
    std::vector<int> partnerIdx;
    std::vector<bool> mask;


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

    // find closest member for every c beta atom
    void findResiduePartners(std::vector<int> & partnerIdx, std::vector<Vec3> & cb,
                             std::vector<bool> & validMask, const size_t len);

    // Describe interaction of residue i and j
    Feature calcFeatures(std::vector<Vec3> & ca, int i, int j);

    void calcConformationDescriptors(std::vector<Feature> & features, std::vector<int> & partnerIdx,
                                     std::vector<Vec3> & ca, std::vector<bool> & mask, const size_t len);

    void discretizeFeatures(std::vector<char> & states, std::vector<Feature> & features,
                            std::vector<bool> & mask, const size_t len);

    void createResidueMask(std::vector<bool> & validMask, std::vector<Vec3> & ca,
                           std::vector<Vec3> & n, std::vector<Vec3> & c, const size_t len);

};