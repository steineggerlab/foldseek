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

    FwBwAligner(size_t length, SubstitutionMatrix &subMat, float gapOpen, float gapExtend, float mact, float T, size_t rowsCapacity, size_t colsCapacity);
    FwBwAligner(size_t length, float gapOpen, float gapExtend, float temperature, size_t rowsCapacity, size_t colsCapacity);
    ~FwBwAligner();

    s_align computeAlignment();
    size_t getRowsCapacity() const { return rowsCapacity; }
    size_t getColsCapacity() const { return colsCapacity; }
    size_t getBlockLength() const { return length; }
    void resizeMatrix(size_t newRowsCapacity, size_t newColsCapacity);

    void initQueryProfile(unsigned char* queryNum, size_t queryLen);
    void initAlignment(unsigned char* targetNum, size_t targetLen);
    float** getZm() {return zm;}
    void computeProbabilityMatrix(bool has_Profile);
    void initScoreMatrix(float** inputScoreMatrix, size_t queryLen, size_t targetLen, int * gaps);
    void setParams(float go, float ge, float t, size_t l);

    unsigned char* queryNum = nullptr;
    unsigned char* targetNum = nullptr;
    float** zm;
    float** scoreForward = nullptr;
    float temperature;
    float maxP = 0;


private:

    float** P;
    float* zmFirst;
    float* zeFirst;
    float* zfFirst;
    float* zmBlockPrev;
    float* zmBlockCurr;
    float* zeBlock;
    float* zfBlock;
    float* vj;
    float* wj;
    float** zInit;
    float* exp_ge_arr;

    
    float** scoreForwardProfile = nullptr;
    float** scoreForwardProfile_exp = nullptr;
    float** scoreBackwardProfile_exp = nullptr;
    
    uint8_t** btMatrix = nullptr;
    float** blosum = nullptr;
    float* S_prev = nullptr;
    float* S_curr = nullptr;    

    // float** zmaxBlocksMaxForward;
    // float** zmaxBlocksMaxBackward;

    size_t length;
    // const SubstitutionMatrix & subMat;
    // static SubstitutionMatrix & defaultSubMat;
    float gapOpen;
    float gapExtend;
    float mact;
    
    size_t rowsCapacity;
    size_t colsCapacity;

    size_t blockCapacity;
    size_t blocks;
    size_t qlen;
    size_t tlen;


    simd_float exp_go;
    simd_float exp_ge;
    float max_zm;
    
    
    simd_float vMax_zm;
    size_t qlen_padding;
    void forward(bool has_Profile);
    void backward(bool has_Profile);
    // void forwardBackwardSaveBlockMaxLocal(bool isForward, size_t start, size_t memcpy_cols);
    void reallocateScoreProfile(size_t newColsCapacity);
};


#endif //FWBW_H