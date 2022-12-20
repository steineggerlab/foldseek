#ifndef COORDINATE16_H
#define COORDINATE16_H

#include "LocalParameters.h"

#include "simd.h"
#include "simde/x86/avx2.h"
#include "simde/x86/f16c.h"

#include <vector>

class Coordinate16 {
public:
    Coordinate16(int type) : type(type) {}

    float* read(const char* mem, size_t chainLength) {
        if (type == LocalParameters::DBTYPE_CA_ALPHA) {
            return (float*) mem;
        }
        if (type == LocalParameters::DBTYPE_CA_ALPHA_DIFF) {
            buffer.reserve(chainLength * 3);
            int32_t diffSum = 0;
            int32_t start;
            const char* data = mem;
            memcpy(&start, data, sizeof(int32_t));
            data += sizeof(int32_t);
            buffer[0] = start / 1000.0f;
            int16_t intDiff = 0;
            for (size_t i = 1; i < chainLength; ++i) {
                memcpy(&intDiff, data, sizeof(int16_t));
                data += sizeof(int16_t);
                diffSum += intDiff;
                buffer[i] = (start + diffSum) / 1000.0f;
            }
            diffSum = 0;
            memcpy(&start, data, sizeof(int32_t));
            data += sizeof(int32_t);
            buffer[chainLength] = start / 1000.0f;
            for (size_t i = chainLength + 1; i < 2 * chainLength; ++i) {
                memcpy(&intDiff, data, sizeof(int16_t));
                data += sizeof(int16_t);
                diffSum += intDiff;
                buffer[i] = (start + diffSum) / 1000.0f;
            }
            diffSum = 0;
            memcpy(&start, data, sizeof(int32_t));
            data += sizeof(int32_t);
            buffer[2 * chainLength] = start / 1000.0f;
            for (size_t i = 2 * chainLength + 1; i < 3 * chainLength; ++i) {
                memcpy(&intDiff, data, sizeof(int16_t));
                data += sizeof(int16_t);
                diffSum += intDiff;
                buffer[i] = (start + diffSum) / 1000.0f;
            }
            return buffer.data();
        }
        if (type == LocalParameters::DBTYPE_CA_ALPHA_F16) {
            size_t xyzLength = chainLength * 3;
            size_t simdXyzLength = xyzLength - (chainLength % 8);
            buffer.reserve(chainLength * 3);

            const int16_t* mem16 = (const int16_t*) mem;
            for (size_t i = 0; i < simdXyzLength; i += VECSIZE_FLOAT) {
                __m128i res = _mm_loadu_epi16((const __m128i*)(mem16 + i));
    #ifdef AVX2
                __m256 res2 = _mm256_cvtph_ps(res);
                _mm256_storeu_ps(buffer.data() + i, res2);
    #else
                __m128 res2 = _mm_cvtph_ps(res);
                _mm_storeu_ps(buffer.data() + i, res2);
    #endif
            }
            for (size_t i = simdXyzLength; i < xyzLength; i++) {
                buffer[i] = _mm_cvtss_f32(_mm_cvtph_ps(_mm_set1_epi16(*(mem16 + i))));
            }
            return buffer.data();
        }
        return NULL;
    }

private:
    std::vector<float> buffer;
    int type;
};

#endif