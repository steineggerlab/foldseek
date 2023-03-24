#ifndef COORDINATE16_H
#define COORDINATE16_H

#include "LocalParameters.h"
#include <vector>

class Coordinate16 {
public:
    Coordinate16() : buffer(NULL), bufferSize(0) {}
    ~Coordinate16(){
        if(buffer != NULL){
            free(buffer);
        }
    }
    float* read(const char* mem, size_t chainLength, size_t entryLength) {
        if (entryLength >= (chainLength * 3) * sizeof(float)) {
            return (float*) mem;
        }
        if(bufferSize < (chainLength * 3)){
            buffer = (float *)realloc(buffer, (chainLength * 3) * sizeof(float));
            bufferSize = (chainLength * 3);
        }
        const char* data = mem;
        int32_t diffSum = 0;
        int32_t start;
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
        return buffer;
    }

    template <typename T>
    static bool convertToDiff16(size_t len, T* data, int16_t* out, int stride = 3) {
        int32_t last = (int)(data[0] * 1000);
        memcpy(out, &last, sizeof(int32_t));
        for (size_t i = 1; i < len; ++i) {
            int32_t curr = (int32_t)(data[i * stride] * 1000);
            int16_t diff;
            bool overflow = __builtin_sub_overflow(curr, last, &diff);
            if (overflow == 1) {
                return true;
            }
            out[i + 1] = diff;
            last = curr;
        }
        return false;
    }

private:
    float * buffer;
    size_t bufferSize;
};

#endif