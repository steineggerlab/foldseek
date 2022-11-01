#ifndef COORDINATE16_H
#define COORDINATE16_H

#include "simd.h"
#include "simde/x86/avx2.h"
#include "simde/x86/f16c.h"

#include <vector>

class Coordinate16 {
public:
    void read(const char* mem, size_t chainLength) {
        size_t xyzLength = chainLength * 3;
        size_t simdXyzLength = xyzLength - (chainLength % 8);
        buffer.reserve(chainLength * 3);

        const int16_t* mem16 = (const int16_t*) mem;
        for (size_t i = 0; i < simdXyzLength; i += 8) {
            __m128i res = _mm_loadu_epi16((const __m128i*)(mem16 + i));
            __m256 res2 = _mm256_cvtph_ps(res);
            _mm256_storeu_ps(buffer.data() + i, res2);
        }
        for (size_t i = simdXyzLength; i < xyzLength; i++) {
            buffer[i] = _mm_cvtss_f32(_mm_cvtph_ps(_mm_set1_epi16(*(mem16 + i))));
        }
    }

    float* getBuffer() {
        return buffer.data();
    }

private:
    std::vector<float> buffer;
};

#endif