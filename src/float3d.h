#pragma once
#include <cmath>

struct float3d {
    float3d() : x(0), y(0), z(0) {}
    float3d(float x, float y, float z) : x(x), y(y), z(z) {};
    float x;
    float y;
    float z;
};

/**
 * @brief Return the cross product of two vectors
 *
 * @param v1 A 3d vector of float
 * @param v2 A 3d vector of float
 * @return std::vector<float>
 */
static inline float3d crossProduct(float3d v1, float3d v2) {
    float x = v1.y * v2.z - v2.y * v1.z;
    float y = v1.z * v2.x - v2.z * v1.x;
    float z = v1.x * v2.y - v2.x * v1.y;
    return { x, y, z };
}

/**
 * @brief Return the norm of given vector
 *
 * @param v A 3d vector of float
 * @return float
 */
static inline float norm(float3d v) {
    return sqrt(pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2));
}

static inline float getCosineTheta(float3d v1, float3d v2) {
    // TODO: Check the length of v1, v2 to be 3
    // Calculate inner product of two vectors
    float inner_product = (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
    float v1_size = pow(v1.x, 2) + pow(v1.y, 2) + pow(v1.z, 2);
    float v2_size = pow(v2.x, 2) + pow(v2.y, 2) + pow(v2.z, 2);
    float output = inner_product / sqrt(v1_size * v2_size);
    return output;
}

static inline float distance(float3d atm1, float3d atm2) {
    float output = 0.0;
    output = sqrt(
        (pow(atm1.x - atm2.x, 2) +
            pow(atm1.y - atm2.y, 2) +
            pow(atm1.z - atm2.z, 2))
    );
    return output;
}

static inline float angle(float3d atm1, float3d atm2, float3d atm3) {
    float3d d1{
        (atm1.x - atm2.x), (atm1.y - atm2.y), (atm1.z - atm2.z)
    };
    float3d d2{
        (atm3.x - atm2.x), (atm3.y - atm2.y), (atm3.z - atm2.z)
    };
    float cos_theta = getCosineTheta(d1, d2);
    float theta = acos(cos_theta) * 180.0 / M_PI;
    return theta;
}
