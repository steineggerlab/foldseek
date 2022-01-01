//
// Created by Martin Steinegger on 11/26/20.
//

#ifndef STRUCCLUST_AFFINENEEDLEMANWUNSCH_H
#define STRUCCLUST_AFFINENEEDLEMANWUNSCH_H
#include "simd.h"
#include <iostream>
#include <iostream>
#include <unistd.h>
#include <float.h>

class AffineNeedlemanWunsch{

public:

    AffineNeedlemanWunsch(int maxLen, int profileRange);
    ~AffineNeedlemanWunsch();

    static const unsigned int SUBMAT = 0;
    static const unsigned int XYZ = 1;
    static const unsigned int XYZ_SS = 2;

    typedef struct matrix {
        const char * name;
        const int *matrix;
        const int *mapper;
        int size;
        int max;
        int min;
    } matrix_t;

    typedef struct cigar_ {
        uint32_t *seq;
        int len;
        int beg_query;
        int beg_ref;
    } cigar_t;


    typedef struct profile_data {
        void * score1;
        void * score2;
    } profile_data_t;

    typedef struct profile {
        const char *s1;
        int s1Len;
        const matrix_t *matrix;
        struct profile_data profile32;
        void (*free)(void * profile);
        int stop;
    } profile_t;

    profile_t * profile_create(
            const char * s1, const int s1Len,
            const matrix_t *matrix);

    profile_t * profile_xyz_create(
            const char * s1, const int s1Len,
            const float* x, const float* y, const float* z);

    profile_t * profile_xyz_ss_create(
            const char * s1, const int s1Len, const matrix_t *matrix,
            const float* x, const float* y, const float* z);

    typedef struct result_extra_trace {
        void * trace_table;    /* DP table of traceback */

    } result_extra_trace_t;

    typedef struct result {
        float score;
        int end_query;
        int end_ref;
        int flag;       /* bit field for various flags */
        /* union of pointers to extra result data based on the flag */
        union {
            result_extra_trace_t *trace;
        };
    } result_t;

    typedef struct alignment {
        alignment(float score,int start_query, int end_query,
                  int start_target, int end_target,
                  int cigar_len, uint32_t * cigar):
                score(score), start_query(start_query), end_query(end_query),
                start_target(start_target), end_target(end_target),
                cigar_len(cigar_len), cigar(cigar){}
        float score;
        int start_query;
        int end_query;
        int start_target;
        int end_target;
        int cigar_len;
        uint32_t *cigar;
    } alignment_t;

    alignment_t alignXYZ(AffineNeedlemanWunsch::profile_t *profile,
                         long queryLen, long targetLen,
                         const float * targetX, const float * targetY, const float * targetZ,
                         const float d02, float t[3], float u[3][3], float gapopen, float gapextend,
                         int * invmap);


    alignment_t alignXYZ_SS(AffineNeedlemanWunsch::profile_t *profile,
                            long queryLen, long targetLen,
                            const float * targetX, const float * targetY, const float * targetZ,
                            const char * target_ss, const float d02, float t[3],
                            float u[3][3], float gapopen, float gapextend,
                            int * invmap);

    alignment_t align(AffineNeedlemanWunsch::profile_t *profile, long queryLen,
                      const unsigned char *target, long targetLen,
                      float gapopen, float gapextend, int * invmap);

    static char cigar_decode_op(uint32_t cigar_int) {
#define BAM_CIGAR_STR "MIDNSHP=XB"
        return (cigar_int & 0xfU) > 9 ? 'M': BAM_CIGAR_STR[cigar_int & 0xfU];
#undef BAM_CIGAR_STR
    }

    static uint32_t cigar_decode_len(uint32_t cigar_int) {
#define BAM_CIGAR_SHIFT 4u
        return cigar_int >> BAM_CIGAR_SHIFT;
#undef BAM_CIGAR_SHIFT
    }

private:

    simd_float* pvHStore;
    simd_float* pvHLoad;
    simd_float* pvE;
    simd_float* pvEaStore;
    simd_float* pvEaLoad;
    simd_float* pvHT;
    simd_float* traceTableRow;
    simd_float* nextTraceTableRow;
    float * boundary;
    uint32_t * cigarBuffer;
    uint32_t * reverCigarBuffer;
    result_t * result;
    profile_t * profile;
    simd_float* vProfile1;
    simd_float* vProfile2;

    typedef union simd_float_32 {
        simd_float m;
        float v[VECSIZE_FLOAT];
    } simd_float_32_t;

    static inline float _mm_extract_ps_rpl(simd_float a, const int imm) {
        simd_float_32_t A;
        A.m = a;
        return A.v[imm];
    }

    result_t* result_new();
    result_t* result_new_trace(const int a8bit, const int b, const size_t alignment, const size_t size);

    template<unsigned int T>
    result_t* stripedAlign(
            const profile_t * profile,
            const char * s2, const int s2Len,
            const float * targetX, const float * targetY, const float * targetZ,
            const float d02, float t[3], float u[3][3],
            const float open, const float gap);


    cigar_t cigar_striped_32 (
            int lena,
            int lenb,
            result_t *result,
            int * j2i);


    uint32_t* reverse_uint32_t(const uint32_t *s, size_t length);
    uint32_t cigar_encode(uint32_t length, char op_letter);

};

#endif //STRUCCLUST_AFFINENEEDLEMANWUNSCH_H
