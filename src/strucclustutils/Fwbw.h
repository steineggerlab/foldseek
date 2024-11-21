#ifndef FWBW
#define FWBW

#include "SubstitutionMatrix.h"
#include "IndexReader.h"
#include "DBReader.h"

#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <string>

class FwBwAligner {
public:
    typedef struct {
        float score1;
        float score2;
        int32_t dbStartPos1;
        int32_t dbEndPos1;
        int32_t	qStartPos1;
        int32_t qEndPos1;
        int32_t ref_end2;
        float qCov;
        float tCov;
        std::string cigar;
        double evalue;
        int identicalAACnt;
        int32_t cigarLen;
        int word;
    } s_align;

    FwBwAligner(size_t maxQueryLen, size_t maxTargetLen, size_t length, size_t blocks);
    ~FwBwAligner();

    s_align align(
        const std::string& querySeqAA,
        const std::string& querySeq3Di,
        const std::string& targetSeqAA,
        const std::string& targetSeq3Di,
        size_t queryLen,
        size_t targetLen,
        size_t length,
        size_t blocks,
        SubstitutionMatrix & subMatAA,
        SubstitutionMatrix & subMat3Di,
        float mact,
        float go,
        float ge,
        float T
    );

    void computeForwardScoreMatrix(const unsigned char* queryNum, const unsigned char* targetNum,
                                   unsigned int queryLen, unsigned int targetLen,
                                   float** blosum, float T, float ** scoreForward);
    void computeForwardScoreMatrix(
        const unsigned char* queryNumAA,
        const unsigned char* queryNum3Di,
        const unsigned char* targetNumAA,
        const unsigned char* targetNum3Di,
        unsigned int queryLen,
        unsigned int targetLen,
        SubstitutionMatrix &subMatAA,
        SubstitutionMatrix &subMat3Di,
        float T,
        float** scoreForward
    );

    unsigned char* seq2num(const std::string& seq, const unsigned char* aa2num) {
        unsigned char* idx = static_cast<unsigned char*>(malloc(seq.size() * sizeof(unsigned char)));
        for (size_t i = 0; i < seq.size(); ++i) {
            idx[i] = aa2num[static_cast<unsigned char>(seq[i])];
        }
        return idx;
    }

    template <typename T>
    T** malloc_matrix(int dim1, int dim2) {
// #define ALIGN_FLOAT 16
#define ICEIL(x_int, fac_int) ((x_int + fac_int - 1) / fac_int) * fac_int
        // Compute mem sizes rounded up to nearest multiple of ALIGN_FLOAT
        size_t size_pointer_array = ICEIL(dim1 * sizeof(T*), ALIGN_FLOAT);
        size_t dim2_padded = ICEIL(dim2 * sizeof(T), ALIGN_FLOAT) / sizeof(T);

        T** matrix = (T**) mem_align(ALIGN_FLOAT, size_pointer_array + dim1 * dim2_padded * sizeof(T));
        if (matrix == NULL)
            return matrix;

        T* ptr = (T*) (matrix + (size_pointer_array / sizeof(T*)));
        for (int i = 0; i < dim1; ++i) {
            matrix[i] = ptr;
            ptr += dim2_padded;
        }
#undef ICEIL
        return matrix;
    }

    inline void* mem_align(size_t boundary, size_t size) {
        void* pointer;
        if (posix_memalign(&pointer, boundary, size) != 0) {
#define MEM_ALIGN_ERROR "mem_align could not allocate memory.\n"
            fwrite(MEM_ALIGN_ERROR, sizeof(MEM_ALIGN_ERROR), 1, stderr);
#undef MEM_ALIGN_ERROR
            exit(3);
        }
        return pointer;
    }

private:

    float** zmForward;
    float** zeForward;
    float** zfForward;
    float** zmBackward;
    float** zeBackward;
    float** zfBackward;
    float** scoreForward;
    float** scoreBackward;
    float** P;
    // uint8_t* btMatrix;
    uint8_t** btMatrix;


    float** zmaxBlocksMaxForward;
    float** zmaxBlocksMaxBackward;

    float* zmaxForward;
    float* zmaxBackward;
    float *vj;
    float *wj;
    // float** mat3di;
    float** blosum;

    // float gapOpen;
    // float gapExtend;
    // float T;

    // const SubstitutionMatrix & subMat3Di;
    // const SubstitutionMatrix & subMat;

    void forwardBackwardSaveBlockMaxLocal(float** S, float** z_init, float* vj, float* wj,
                                          float T, float go, float ge,
                                          size_t rows, size_t start, size_t end, size_t memcpy_cols,
            // output is zm, ze, ze, zmax
                                          float** zm, float** ze, float** zf, float* zmax);

    void rescaleBlocks(float **matrix, float **scale, size_t rows, size_t length, size_t blocks, size_t targetLen);

};


#endif //FWBW_H