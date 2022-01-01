#include "affineneedlemanwunsch.h"

#include <iostream>
#include <iostream>
#include <unistd.h>
#include <float.h>
#include <cstring>


AffineNeedlemanWunsch::AffineNeedlemanWunsch(int maxLen, int profileRange){
    const int32_t segWidth = VECSIZE_FLOAT; /* number of values in vector unit */
    const int32_t segLen = (maxLen + segWidth - 1) / segWidth;
    const int32_t segWidth8Bit = segWidth * 4; /* number of values in vector unit */
    const int32_t segLen8Bit = (maxLen + segWidth8Bit - 1) / segWidth8Bit;

    pvHStore          =  (simd_float*) mem_align(ALIGN_INT, segLen*sizeof(simd_float));
    pvHLoad           =  (simd_float*) mem_align(ALIGN_INT, segLen*sizeof(simd_float));
    pvE               =  (simd_float*) mem_align(ALIGN_INT, segLen*sizeof(simd_float));
    pvEaStore         =  (simd_float*) mem_align(ALIGN_INT, segLen*sizeof(simd_float));
    pvEaLoad          =  (simd_float*) mem_align(ALIGN_INT, segLen*sizeof(simd_float));
    pvHT              =  (simd_float*) mem_align(ALIGN_INT, segLen*sizeof(simd_float));
    traceTableRow     =  (simd_float*) mem_align(ALIGN_INT, (segLen+VECSIZE_FLOAT-1)*sizeof(simd_float));
    nextTraceTableRow =  (simd_float*) mem_align(ALIGN_INT, (segLen+VECSIZE_FLOAT-1)*sizeof(simd_float));
    boundary          =  (float*) mem_align(ALIGN_INT, (maxLen+1)*sizeof(float));
    cigarBuffer       =  (uint32_t *) malloc(sizeof(uint32_t)*(maxLen+maxLen));
    reverCigarBuffer  =  (uint32_t *) malloc(sizeof(uint32_t)*(maxLen+maxLen));
    result = result_new_trace(segLen8Bit, maxLen, ALIGN_INT, sizeof(simd_float));
    profile = (profile_t*)malloc(sizeof(profile_t));
    vProfile1 =  (simd_float*) mem_align(ALIGN_INT, profileRange * segLen * sizeof(simd_float));
    vProfile2 =  (simd_float*) mem_align(ALIGN_INT, profileRange * segLen*sizeof(simd_float));
}

AffineNeedlemanWunsch::~AffineNeedlemanWunsch(){
    free(pvHStore);
    free(pvHLoad);
    free(pvE);
    free(pvEaStore);
    free(pvEaLoad);
    free(pvHT);
    free(traceTableRow);
    free(nextTraceTableRow);
    free(boundary);
    free(cigarBuffer);
    free(reverCigarBuffer);
    free(result->trace->trace_table);
    free(result->trace);
    free(result);
    free(profile);
    free(vProfile1);
    free(vProfile2);
}

AffineNeedlemanWunsch::alignment_t AffineNeedlemanWunsch::alignXYZ_SS(
                        AffineNeedlemanWunsch::profile_t *profile,
                        long queryLen, long targetLen,
                        const float * targetX, const float * targetY, const float * targetZ,
                        const char  *target_ss, const float d02, float t[3],
                        float u[3][3], float gapopen, float gapextend,
                        int * invmap) {
    // fill matrix
    result_t *result = stripedAlign<XYZ_SS>(
            profile, target_ss, targetLen,
                    targetX, targetY, targetZ,
                    d02, t, u, gapopen, gapextend);
    //std::cout << result->score << std::endl;
    // compute backtrace
    cigar_t cigar = cigar_striped_32(
            queryLen,
            targetLen,
            result,
            invmap);
    return alignment_t(result->score, cigar.beg_query,result->end_query,
                       cigar.beg_ref, result->end_ref,
                       cigar.len, cigar.seq);
}

AffineNeedlemanWunsch::alignment_t AffineNeedlemanWunsch::alignXYZ(AffineNeedlemanWunsch::profile_t *profile,
                                                                   long queryLen, long targetLen,
                                                                   const float * targetX, const float * targetY, const float * targetZ,
                                                                   const float d02, float t[3], float u[3][3],
                                                                   float gapopen, float gapextend, int * invmap) {
    // fill matrix
    result_t *result = stripedAlign<XYZ>(
            profile, NULL, targetLen,
            targetX, targetY, targetZ,
            d02, t, u, gapopen, gapextend);
    //std::cout << result->score << std::endl;
    // compute backtrace
    cigar_t cigar = cigar_striped_32(
            queryLen,
            targetLen,
            result,
            invmap);
    return alignment_t(result->score, cigar.beg_query,result->end_query,
                       cigar.beg_ref, result->end_ref,
                       cigar.len, cigar.seq);
}



AffineNeedlemanWunsch::alignment_t AffineNeedlemanWunsch::align(AffineNeedlemanWunsch::profile_t *profile, long queryLen,
                                                                const unsigned char *target, long targetLen,
                                                                float gapopen, float gapextend, int * invmap) {
    // fill matrix
    result_t *result = stripedAlign<SUBMAT>(
            profile, (const char *) target, targetLen,
                    NULL, NULL, NULL, 0.0, NULL, NULL,
                    gapopen, gapextend);
    //std::cout << result->score << std::endl;

    // compute backtrace
    cigar_t cigar = cigar_striped_32(
            queryLen,
            targetLen,
            result,
            invmap);
    return alignment_t(result->score, cigar.beg_query,result->end_query,
                       cigar.beg_ref, result->end_ref,
                       cigar.len, cigar.seq);
}


typedef union simd_float_32 {
    simd_float m;
    float v[4];
} simd_float_32_t;

static inline float _mm_extract_ps_rpl(simd_float a, const int imm) {
    simd_float_32_t A;
    A.m = a;
    return A.v[imm];
}



AffineNeedlemanWunsch::profile_t * AffineNeedlemanWunsch::profile_xyz_ss_create(
        const char * s1, const int s1Len, const matrix_t *matrix,
        const float * x,const float * y, const float * z){
    int32_t i = 0;
    int32_t j = 0;
    int32_t segNum = 0;
    const int32_t segWidth = VECSIZE_FLOAT; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    int32_t index = 0;

    profile->s1 = s1;
    profile->s1Len = s1Len;
    profile->matrix = matrix;
    profile->profile32.score1 = NULL;
    profile->stop = INT32_MAX;
    for (i=0; i<segLen; ++i) {
        simd_float_32_t vX;
        simd_float_32_t vY;
        simd_float_32_t vZ;
        j = i;
        for (segNum=0; segNum<segWidth; ++segNum) {
            vX.v[segNum] = j >= s1Len ? FLT_MIN : x[j];
            vY.v[segNum] = j >= s1Len ? FLT_MIN : y[j];
            vZ.v[segNum] = j >= s1Len ? FLT_MIN : z[j];
            j += segLen;
        }
        simdf32_store((float*)&vProfile1[index], vX.m);
        simdf32_store((float*)&vProfile1[index + 1], vY.m);
        simdf32_store((float*)&vProfile1[index + 2], vZ.m);
        index+=3;
    }
    index = 0;
    for (int k=0; k < matrix->size; ++k) {
        for (i=0; i<segLen; ++i) {
            simd_float_32_t t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t.v[segNum] = (float) j >= s1Len ? 0.0f : (k == matrix->mapper[(unsigned char)s1[j]]) ? 0.5 : 0.0;
                j += segLen;
            }
            simdf32_store((float*)&vProfile2[index], t.m);
            ++index;
        }
    }

    profile->profile32.score1 = vProfile1;
    profile->profile32.score2 = vProfile2;
    return profile;
}



AffineNeedlemanWunsch::profile_t * AffineNeedlemanWunsch::profile_xyz_create(
        const char * s1, const int s1Len,
        const float * x,const float * y, const float * z)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t segNum = 0;
    const int32_t segWidth = VECSIZE_FLOAT; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    int32_t index = 0;

    profile->s1 = s1;
    profile->s1Len = s1Len;
    profile->matrix = NULL;
    profile->profile32.score1 = NULL;
    profile->stop = INT32_MAX;
    for (i=0; i<segLen; ++i) {
        simd_float_32_t vX;
        simd_float_32_t vY;
        simd_float_32_t vZ;
        j = i;
        for (segNum=0; segNum<segWidth; ++segNum) {
            vX.v[segNum] = j >= s1Len ? FLT_MIN : x[j];
            vY.v[segNum] = j >= s1Len ? FLT_MIN : y[j];
            vZ.v[segNum] = j >= s1Len ? FLT_MIN : z[j];
            j += segLen;
        }
        simdf32_store((float*)&vProfile1[index], vX.m);
        simdf32_store((float*)&vProfile1[index + 1], vY.m);
        simdf32_store((float*)&vProfile1[index + 2], vZ.m);
        index+=3;
    }

    profile->profile32.score1 = vProfile1;
    return profile;
}


AffineNeedlemanWunsch::profile_t * AffineNeedlemanWunsch::profile_create(
        const char * s1, const int s1Len,
        const matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = VECSIZE_FLOAT; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    int32_t index = 0;

    profile->s1 = s1;
    profile->s1Len = s1Len;
    profile->matrix = matrix;
    profile->profile32.score1 = NULL;
    profile->stop = INT32_MAX;

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simd_float_32_t t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t.v[segNum] = (float) j >= s1Len ? 0.0f : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            simdf32_store((float*)&vProfile1[index], t.m);
            ++index;
        }
    }

    profile->profile32.score1 = vProfile1;
    return profile;
}

AffineNeedlemanWunsch::result_t* AffineNeedlemanWunsch::result_new()
{
    /* declare all variables */
    result_t *result = NULL;
    result = (result_t*)malloc(sizeof(result_t));
    result->score = 0;
    result->end_query = 0;
    result->end_ref = 0;
    result->flag = 0;
    return result;
}

AffineNeedlemanWunsch::result_t* AffineNeedlemanWunsch::result_new_trace(const int a8bit, const int b,
                                                                         const size_t alignment, const size_t size)
{
    /* declare all variables */
    result_t *result = NULL;
    /* allocate struct to hold memory */
    result = result_new();
    result->trace = (result_extra_trace_t*)malloc(sizeof(result_extra_trace_t));
    result->trace->trace_table = mem_align(alignment, size*a8bit*b);
    return result;
}

#define ZERO_MASK 120 /* all bits set except the first three */
#define E_MASK 103 /* all bits set except the E bits */
#define F_MASK 31 /* all bits set except the F bits */
#define ZERO   0
#define INS    1
#define DEL    2
#define DIAG   4
#define DIAG_E 8
#define INS_E  16
#define DIAG_F 32
#define DEL_F  64



template <unsigned int T>
AffineNeedlemanWunsch::result_t* AffineNeedlemanWunsch::stripedAlign(
        const AffineNeedlemanWunsch::profile_t * profile,
        const char * s2, const int s2Len,
        const float * targetX, const float * targetY, const float * targetZ,
        const float d02, float t[3], float u[3][3], const float open, const float gap)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    const int s1Len = profile->s1Len;
    int32_t end_query = s1Len-1;
    int32_t end_ref = s2Len-1;
    const matrix_t *matrix = profile->matrix;
    const int32_t segWidth = VECSIZE_FLOAT; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t segWidth8Bit = segWidth * 4; /* number of values in vector unit */
    const int32_t segLen8Bit = (s1Len + segWidth8Bit - 1) / segWidth8Bit;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    const simd_float* vProfile = (simd_float*)profile->profile32.score1;
    const simd_float* vProfile2 = (simd_float*) profile->profile32.score2;

    simd_float vGapO = simdf32_set(open);
    simd_float vGapE = simdf32_set(gap);
    simd_float vNegInf = simdf32_set(FLT_MIN);
    simd_float one = simdf32_set(1.0f);
    simd_float vd02 = simdf32_set(d02);

    float score = FLT_MIN;

    simd_float vTIns   = (simd_float) simdi32_set(INS);
    simd_float vTDel   = (simd_float) simdi32_set(DEL);
    simd_float vTDiag  = (simd_float) simdi32_set(DIAG);
    simd_float vTDiagE = (simd_float) simdi32_set(DIAG_E);
    simd_float vTInsE  = (simd_float) simdi32_set(INS_E);
    simd_float vTDiagF = (simd_float) simdi32_set(DIAG_F);
    simd_float vTDelF  = (simd_float) simdi32_set(DEL_F);
    simd_float vTMask  = (simd_float) simdi32_set(ZERO_MASK);
    simd_float vFTMask = (simd_float) simdi32_set(F_MASK);

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            int32_t segNum = 0;
            simd_float_32_t h;
            simd_float_32_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < FLT_MIN ? FLT_MIN : tmp;
                tmp = tmp - open;
                e.v[segNum] = tmp < FLT_MIN ? FLT_MIN : tmp;
            }
            simdf32_store((float*)&pvHStore[index], h.m);
            simdf32_store((float*)&pvE[index], e.m);
            simdf32_store((float*)&pvEaStore[index], e.m);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < FLT_MIN ? FLT_MIN : tmp;
        }
    }

    for (i=0; i<segLen; ++i) {
        simdf32_store((float*)(nextTraceTableRow+i), (simd_float) vTDiagE);
        //arr_store((simd_float*)result->trace->trace_table, (simd_float)vTDiagE, i, segLen, 0);
    }

    simd_float vEF_opn = simdf32_set(0.0);
    simd_float vF_ext = simdf32_set(0.0);
    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simd_float vE;
        simd_float vE_ext;
        simd_float vF;
        simd_float vFa;
        simd_float vFa_ext;
        simd_float vH;
        simd_float vH_dag;
        simd_float vX;
        simd_float vY;
        simd_float vZ;

        const simd_float* vP = NULL;
        const simd_float* vP2 = NULL;

        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        vF = vNegInf;

        /* load final segment of pvHStore and shift left by 4 bytes */
        vH = simdf32_load((float*)&pvHStore[segLen - 1]);
        vH = (simd_float) simdi8_shiftl((simd_int)vH, 4);

        /* insert upper boundary condition */
        vH = (simd_float) simdi32_insert((simd_int)vH,(int) boundary[j], 0);

        /* Correct part of the vProfile1 */
        if(T == XYZ || T == XYZ_SS){
            float xx[3];
            //Transform
            xx[0]=t[0]+(u[0][0] * targetX[j] + u[0][1] * targetY[j] + u[0][2] * targetZ[j]);
            xx[1]=t[1]+(u[1][0] * targetX[j] + u[1][1] * targetY[j] + u[1][2] * targetZ[j]);
            xx[2]=t[2]+(u[2][0] * targetX[j] + u[2][1] * targetY[j] + u[2][2] * targetZ[j]);;
            vX = simdf32_set(xx[0]);
            vY = simdf32_set(xx[1]);
            vZ = simdf32_set(xx[2]);
            vP = vProfile;
            if(T == XYZ_SS){
                vP2 = vProfile2 + matrix->mapper[(unsigned char)s2[j]] * segLen;
            }
        }else {
            vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;
        }

#define SWAP(A,B) { simd_float* tmp = A; A = B; B = tmp; }
        /* Swap the 2 H buffers. */
        SWAP(traceTableRow, nextTraceTableRow)
        SWAP(pvHLoad, pvHStore)
        SWAP(pvEaLoad, pvEaStore)
#undef SWAP
        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vE = simdf32_load((float*)(pvE + i));

            /* Get max from vH, vE and vF. */
            if(T == XYZ||T == XYZ_SS) {
                simd_float queryX = simdf32_load((float*)(vP + i*3));
                simd_float queryY = simdf32_load((float*)(vP + i*3+1));
                simd_float queryZ = simdf32_load((float*)(vP + i*3+2));
                queryX = simdf32_sub(queryX, vX);
                queryY = simdf32_sub(queryY, vY);
                queryZ = simdf32_sub(queryZ, vZ);
                queryX = simdf32_mul(queryX, queryX);
                queryY = simdf32_mul(queryY, queryY);
                queryZ = simdf32_mul(queryZ, queryZ);
                simd_float res = simdf32_add(queryX, queryY);
                simd_float dij = simdf32_add(res, queryZ);
                simd_float oneDividedDist = simdf32_div(one, simdf32_add(one, simdf32_div(dij,vd02)));
                vH_dag = simdf32_add(vH, oneDividedDist);
                if(T == XYZ_SS){
                    vH_dag = simdf32_add(vH_dag, simdf32_load((float*)(vP2 + i)));
                }
            }else{
                vH_dag = simdf32_add(vH, simdf32_load((float*)(vP + i)));
            }
            vH = simdf32_max(vH_dag, vE);
            vH = simdf32_max(vH, vF);
            /* Save vH values. */
            simdf32_store((float*)(pvHStore + i), vH);

            {
                //simd_float vTAll = arr_load((simd_float*)result->trace->trace_table, i, segLen, j);
                simd_float vTAll = simdf32_load((float*)(traceTableRow + i));
                simd_float case1 = simdf32_eq(vH, vH_dag);
                simd_float case2 = simdf32_eq(vH, vF);
                simd_int TInsTDel = simdi8_blend((simd_int)vTIns, (simd_int)vTDel, (simd_int)case2);
                simd_float vT = (simd_float) simdi8_blend(TInsTDel,(simd_int)vTDiag, (simd_int)case1);
                simdf32_store((float*)(pvHT + i), vT);
                vT = simdf32_or(vT, vTAll);
                simdf32_store((float*)(traceTableRow+i), (simd_float) vT);
                //arr_store((simd_float*)result->trace->trace_table, vT, i, segLen, j);
            }

            vEF_opn = simdf32_sub(vH, vGapO);

            /* Update vE value. */
            vE_ext = simdf32_sub(vE, vGapE);
            vE = simdf32_max(vEF_opn, vE_ext);
            simdf32_store((float*)(pvE + i), vE);
            {
                simd_float vEa = simdf32_load((float*)(pvEaLoad + i));
                simd_float vEa_ext = simdf32_sub(vEa, vGapE);
                vEa = simdf32_max(vEF_opn, vEa_ext);
                simdf32_store((float*)(pvEaStore + i), vEa);
                if (j+1<s2Len) {
                    simd_float cond = simdf32_gt(vEF_opn, vEa_ext);
                    simd_float vT = (simd_float) simdi8_blend((simd_int) vTInsE, (simd_int) vTDiagE, (simd_int) cond);
                    //arr_store((simd_float*)result->trace->trace_table, vT, i, segLen, j+1);
                    simdf32_store((float*)(nextTraceTableRow+i), vT);
                }
            }

            /* Update vF value. */
            vF_ext = simdf32_sub(vF, vGapE);
            vF = simdf32_max(vEF_opn, vF_ext);
            if (i+1<segLen) {
                //simd_float vTAll = arr_load((simd_float*)result->trace->trace_table, i+1, segLen, j);
                simd_float vTAll = simdf32_load((float*)(traceTableRow + i + 1));
                simd_float cond = simdf32_gt(vEF_opn, vF_ext);
                simd_float vT = (simd_float) simdi8_blend((simd_int)vTDelF, (simd_int)vTDiagF, (simd_int)cond);
                vT = simdf32_or(vT, vTAll);
                simdf32_store((float*)(traceTableRow+i+1), (simd_float) vT);
                //arr_store((simd_float*)result->trace->trace_table, vT, i+1, segLen, j);
            }

            /* Load the next vH. */
            vH = simdf32_load((float*)(pvHLoad + i));
        }


        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        vFa_ext = vF_ext;
        vFa = vF;
        for (k=0; k<segWidth; ++k) {
            int64_t tmp = boundary[j+1]-open;
            float tmp2 = tmp < FLT_MIN ? FLT_MIN : tmp;
            simd_float vHp = simdf32_load((float*)&pvHLoad[segLen - 1]);
            vHp = (simd_float) simdi8_shiftl((simd_int)vHp, 4);
            vHp = (simd_float) simdi32_insert((simd_int)vHp, (int)boundary[j], 0);
            vEF_opn = (simd_float) simdi8_shiftl((simd_int)vEF_opn, 4);
            vEF_opn = (simd_float) simdi32_insert((simd_int)vEF_opn,(int) tmp2, 0);
            vF_ext = (simd_float) simdi8_shiftl((simd_int)vF_ext, 4);
            vF_ext = (simd_float) simdi32_insert((simd_int)vF_ext, (int)FLT_MIN, 0);
            vF = (simd_float) simdi8_shiftl((simd_int)vF, 4);
            vF = (simd_float) simdi32_insert((simd_int)vF, (int)tmp2, 0);
            vFa_ext = (simd_float) simdi8_shiftl((simd_int)vFa_ext, 4);
            vFa_ext = (simd_float) simdi32_insert((simd_int)vFa_ext,(int) FLT_MIN, 0);
            vFa = (simd_float) simdi8_shiftl((simd_int)vFa, 4);
            vFa = (simd_float) simdi32_insert((simd_int)vFa, (int)tmp2, 0);
            for (i=0; i<segLen; ++i) {
                vH = simdf32_load((float*)(pvHStore + i));
                vH = simdf32_max(vH, vF);
                simdf32_store((float*)(pvHStore + i), vH);

                {
                    simd_float vTAll;
                    simd_float vT;
                    simd_float case1;
                    simd_float case2;
                    simd_float cond;
                    vHp = simdf32_add(vHp, simdf32_load((float*)(vP + i)));
                    case1 = simdf32_eq(vH, vHp);
                    case2 = simdf32_eq(vH, vF);
                    cond = simdf32_andnot(case1,case2);
//                    vTAll = arr_load((simd_float*)result->trace->trace_table, i, segLen, j);
                    vTAll = simdf32_load((float*)(traceTableRow + i));
                    vT = simdf32_load((float*)(pvHT + i));
                    vT = (simd_float) simdi8_blend((simd_int) vT, (simd_int) vTDel, (simd_int) cond);
                    simdf32_store((float*)(pvHT + i), vT);
                    vTAll = simdf32_and(vTAll, vTMask);
                    vTAll = simdf32_or(vTAll, vT);
                    simdf32_store((float*)(traceTableRow+i), (simd_float) vTAll);
                    //arr_store((simd_float*)result->trace->trace_table, vTAll, i, segLen, j);
                }
                /* Update vF value. */
                {
                    //simd_float vTAll = arr_load((simd_float*)result->trace->trace_table, i, segLen, j);
                    simd_float vTAll = simdf32_load((float*)(traceTableRow + i));
                    simd_float cond = simdf32_gt(vEF_opn, vFa_ext);
                    simd_float vT = (simd_float) simdi8_blend((simd_int) vTDelF, (simd_int) vTDiagF, (simd_int) cond);
                    vTAll = simdf32_and(vTAll, vFTMask);
                    vTAll = simdf32_or(vTAll, vT);
                    simdf32_store((float*)(traceTableRow+i), (simd_float) vTAll);
                    //arr_store((simd_float*)result->trace->trace_table, vTAll, i, segLen, j);
                }
                vEF_opn = simdf32_sub(vH, vGapO);
                vF_ext = simdf32_sub(vF, vGapE);
                {
                    simd_float vEa = simdf32_load((float*)(pvEaLoad + i));
                    simd_float vEa_ext = simdf32_sub(vEa, vGapE);
                    vEa = simdf32_max(vEF_opn, vEa_ext);
                    simdf32_store((float*)(pvEaStore + i), vEa);
                    if (j+1<s2Len) {
                        simd_float cond = simdf32_gt(vEF_opn, vEa_ext);
                        simd_float vT = (simd_float) simdi8_blend((simd_int) vTInsE, (simd_int) vTDiagE, (simd_int) cond);
                        //arr_store((simd_float*)result->trace->trace_table, vT, i, segLen, j+1);
                        simdf32_store((float*)(nextTraceTableRow+i), vT);
                    }
                }
                if (! simdi8_movemask(
                        (simd_int)simdf32_or(
                                simdf32_gt(vF_ext, vEF_opn),
                                simdf32_eq(vF_ext, vEF_opn))))
                    goto end;
                /*vF = _mm_max_epi32_rpl(vEF_opn, vF_ext);*/
                vF = vF_ext;
                vFa_ext = simdf32_sub(vFa, vGapE);
                vFa = simdf32_max(vEF_opn, vFa_ext);
                vHp = simdf32_load((float*)(pvHLoad + i));
            }
        }

        end:
        {
//            for (i=0; i<segLen; i++) {
//                arr_store((simd_float*)result->trace->trace_table, *(traceTableRow+i), i, segLen, j);
//            }
            for (i=0; i<segLen; i+=4) {
                simd_int v0 =(simd_int) *(traceTableRow + i);
                simd_int v1 =(simd_int) *(traceTableRow + i + 1);
                simd_int v2 =(simd_int) *(traceTableRow + i + 2);
                simd_int v3 =(simd_int) *(traceTableRow + i + 3);
//                printf("v0 = %vd\n", v0);
//                printf("v1 = %vd\n", v1);
//                printf("v2 = %vd\n", v2);
//                printf("v3 = %vd\n", v3);

                simd_int v01 = simdi32_pack(v0, v1);        // pack v0, v1 to 16 bits
                simd_int v23 = simdi32_pack(v2, v3);        // pack v2, v3 to 16 bits
                simd_int v0123 = simdi16_pack(v01, v23);    // pack v0...v3 to 8 bits
                size_t pos = (j*segLen8Bit+i/4);
#ifdef AVX2
                v0123 = _mm256_permutevar8x32_epi32(v0123, _mm256_setr_epi32(0,4, 1,5, 2,6, 3,7));
#endif
                simdi_store(((simd_int *)result->trace->trace_table) + pos, v0123);
            }
            if (j+1<s2Len){
//                for (i=0; i<segLen; ++i) {
//                    arr_store((simd_float *) result->trace->trace_table, *(nextTraceTableRow + i), i, segLen, j + 1);
//                }
                for (i=0; i<segLen; i+=4) {
                    simd_int v0 = (simd_int) *(nextTraceTableRow + i);
                    simd_int v1 = (simd_int) *(nextTraceTableRow + i + 1);
                    simd_int v2 = (simd_int) *(nextTraceTableRow + i + 2);
                    simd_int v3 = (simd_int) *(nextTraceTableRow + i + 3);
//                    printf("v0 = %vd\n", v0);
//                    printf("v1 = %vd\n", v1);
//                    printf("v2 = %vd\n", v2);
//                    printf("v3 = %vd\n", v3);

                    simd_int v01 = simdi32_pack(v0, v1);        // pack v0, v1 to 16 bits
                    simd_int v23 = simdi32_pack(v2, v3);        // pack v2, v3 to 16 bits
                    simd_int v0123 = simdi16_pack(v01, v23);    // pack v0...v3 to 8 bits
#ifdef AVX2
                    v0123 = _mm256_permutevar8x32_epi32(v0123, _mm256_setr_epi32(0,4, 1,5, 2,6, 3,7));
#endif
                    size_t pos = (j + 1) * segLen8Bit + i/4;
                    simdi_store(((simd_int *) result->trace->trace_table) + pos,
                                v0123);
                }
            }
        }

//        for (int ti = 0; ti < segLen*segWidth; ti++) {
//            uint32_t * lookup=(uint32_t *)(((simd_float *) result->trace->trace_table)+(1LL*j*segLen));
//            fprintf(stdout, "%3d ", lookup[ti]);
//
//        }
//        fprintf(stdout, "\n");


//        for (int ti = 0; ti < segLen8Bit*segWidth8Bit; ti++) {
//            uint8_t * lookup=(uint8_t *)(((simd_float *) result->trace->sml_trace_table)+(1LL*j*segLen8Bit));
//            fprintf(stdout, "%3d ", lookup[ti]);
//        }
//        fprintf(stdout, "\n");

    }

    /* extract last value from the last column */
    {
        simd_float vH = simdf32_load((float*)(pvHStore + offset));
        for (k=0; k<position; ++k) {
            vH = (simd_float)simdi8_shiftl ((simd_int) vH, 4);
        }
        score = _mm_extract_ps_rpl(vH, VECSIZE_FLOAT-1);
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    return result;
}


/* array index is an ASCII character value from a CIGAR,
   element value is the corresponding integer opcode between 0 and 9 */
const uint8_t cigar_encoded_ops[] = {
        0,         0,         0,         0,
        0,         0,         0,         0,
        0,         0,         0,         0,
        0,         0,         0,         0,
        0,         0,         0,         0,
        0,         0,         0,         0,
        0,         0,         0,         0,
        0,         0,         0,         0,
        0 /*   */, 0 /* ! */, 0 /* " */, 0 /* # */,
        0 /* $ */, 0 /* % */, 0 /* & */, 0 /* ' */,
        0 /* ( */, 0 /* ) */, 0 /* * */, 0 /* + */,
        0 /* , */, 0 /* - */, 0 /* . */, 0 /* / */,
        0 /* 0 */, 0 /* 1 */, 0 /* 2 */, 0 /* 3 */,
        0 /* 4 */, 0 /* 5 */, 0 /* 6 */, 0 /* 7 */,
        0 /* 8 */, 0 /* 9 */, 0 /* : */, 0 /* ; */,
        0 /* < */, 7 /* = */, 0 /* > */, 0 /* ? */,
        0 /* @ */, 0 /* A */, 0 /* B */, 0 /* C */,
        2 /* D */, 0 /* E */, 0 /* F */, 0 /* G */,
        5 /* H */, 1 /* I */, 0 /* J */, 0 /* K */,
        0 /* L */, 0 /* M */, 3 /* N */, 0 /* O */,
        6 /* P */, 0 /* Q */, 0 /* R */, 4 /* S */,
        0 /* T */, 0 /* U */, 0 /* V */, 0 /* W */,
        8 /* X */, 0 /* Y */, 0 /* Z */, 0 /* [ */,
        0 /* \ */, 0 /* ] */, 0 /* ^ */, 0 /* _ */,
        0 /* ` */, 0 /* a */, 0 /* b */, 0 /* c */,
        0 /* d */, 0 /* e */, 0 /* f */, 0 /* g */,
        0 /* h */, 0 /* i */, 0 /* j */, 0 /* k */,
        0 /* l */, 0 /* m */, 0 /* n */, 0 /* o */,
        0 /* p */, 0 /* q */, 0 /* r */, 0 /* s */,
        0 /* t */, 0 /* u */, 0 /* v */, 0 /* w */,
        0 /* x */, 0 /* y */, 0 /* z */, 0 /* { */,
        0 /* | */, 0 /* } */, 0 /* ~ */, 0 /*  */
};

uint32_t AffineNeedlemanWunsch::cigar_encode(uint32_t length, char op_letter)
{
#define BAM_CIGAR_SHIFT 4u
    return (length << BAM_CIGAR_SHIFT) | (cigar_encoded_ops[(int)op_letter]);
#undef BAM_CIGAR_SHIFT
}

uint32_t* AffineNeedlemanWunsch::reverse_uint32_t(const uint32_t *s, size_t length)
{
    uint32_t *r = NULL;
    size_t i = 0;
    size_t j = 0;

    r = reverCigarBuffer;
    for (i=0,j=length-1; i<length; ++i,--j) {
        r[i] = s[j];
    }

    return r;
}



AffineNeedlemanWunsch::cigar_t AffineNeedlemanWunsch::cigar_striped_32 (
        int lena,
        int lenb,
        result_t *result,
        int * j2i)
{

#define UNUSED(expr) do { (void)(expr); } while (0)

#define INC                                                       \
do {                                                              \
    cigar.len += 1;                                              \
    if ((size_t)cigar.len >= size) {                             \
        size = size * 2;                                          \
        cigar.seq = (uint32_t*)realloc(cigar.seq, sizeof(uint32_t)*size);  \
    }                                                             \
} while (0);

#define WRITE(VAL,CHAR)                                         \
do {                                                            \
    INC;                                                        \
    cigar.seq[cigar.len-1] = cigar_encode(VAL,CHAR); \
} while (0)

#define RESET  \
do {           \
    c_mat = 0; \
    c_mis = 0; \
    c_del = 0; \
    c_ins = 0; \
} while (0)

/* internally I accidentally flipped I/D, so rather than go back and
 * rewrite a bunch of code, I fix the problem by just swapping the
 * letters here in the cigar output */
#define WRITE_ANY         \
do {                      \
    if (c_mat) {          \
        WRITE(c_mat,'M'); \
    }                     \
    else if (c_mis) {     \
        WRITE(c_mis,'M'); \
    }                     \
    else if (c_del) {     \
        WRITE(c_del,'I'); \
    }                     \
    else if (c_ins) {     \
        WRITE(c_ins,'D'); \
    }                     \
    RESET;                \
} while (0)

    size_t size = lena+lenb;
    cigar_t cigar;
    uint32_t *cigar_reverse = NULL;
    uint32_t c_mat = 0;
    uint32_t c_mis = 0;
    uint32_t c_del = 0;
    uint32_t c_ins = 0;
    int64_t i = result->end_query;
    int64_t j = result->end_ref;
    int where = DIAG;
    //D *HT = (D*)result->trace->trace_table;
    uint8_t *HT = (uint8_t*)result->trace->trace_table;

    int64_t segWidth = VECSIZE_FLOAT;
    int64_t segLen = (lena + segWidth - 1) / segWidth;
    int64_t segLen8bit = (lena + segWidth*4 - 1) / (segWidth*4);


    cigar.seq = cigarBuffer;
    cigar.len = 0;
    cigar.beg_query = 0;
    cigar.beg_ref = 0;
    /* semi-global alignment includes the end gaps */

    while (i >= 0 || j >= 0) {
        //int64_t loc = j*segLen*segWidth + (i%segLen)*segWidth + (i/segLen);
        int64_t loc = j*segLen8bit*segWidth*4 + (i%segLen)*segWidth + (i/segLen);

        /*assert(i >= 0 && j >= 0);*/
        if (i < 0) {
            if (0 == c_ins) {
                WRITE_ANY;
            }
            while (j >= 0) {
                ++c_ins;
                --j;
            }
            break;
        }
        if (j < 0) {
            if (0 == c_del) {
                WRITE_ANY;
            }
            while (i >= 0) {
                ++c_del;
                --i;
            }
            break;
        }

        if (DIAG == where) {
            if (HT[loc] & DIAG) {
//                char a = case_sensitive ? seqA[i] : toupper(seqA[i]);
//                char b = case_sensitive ? seqB[j] : toupper(seqB[j]);
                j2i[i]=j;
                int matches = true;
//                int matches = (a == b);
//                if (NULL != alphabet_aliases) {
//                    size_t i;
//                    for (i=0; i<aliases_size; i+=1) {
//                        if (alphabet_aliases[i] == a) {
//                            matches |= alphabet_aliases[i+1] == b;
//                        }
//                        else if (alphabet_aliases[i+1] == a) {
//                            matches |= alphabet_aliases[i] == b;
//                        }
//                    }
//                }
                if (matches) {
                    if (0 == c_mat) {
                        WRITE_ANY;
                    }
                    c_mat += 1;
                }
                else {
                    if (0 == c_mis) {
                        WRITE_ANY;
                    }
                    c_mis += 1;
                }
                --i;
                --j;
            }
            else if (HT[loc] & INS) {
                where = INS;
            }
            else if (HT[loc] & DEL) {
                where = DEL;
            }
                /* no bits were set, so this is the zero condition */
            else {
                break;
            }
        }
        else if (INS == where) {
            if (0 == c_ins) {
                WRITE_ANY;
            }
            c_ins += 1;
            --j;
            if (HT[loc] & DIAG_E) {
                where = DIAG;
            }
            else if (HT[loc] & INS_E) {
                where = INS;
            }
            else {
                cigar.beg_ref = -1;
                cigar.beg_query = -1;
                cigar.len = -1;
                return cigar;
            }
        }
        else if (DEL == where) {
            if (0 == c_del) {
                WRITE_ANY;
            }
            c_del += 1;
            --i;
            if (HT[loc] & DIAG_F) {
                where = DIAG;
            }
            else if (HT[loc] & DEL_F) {
                where = DEL;
            }
            else {
                cigar.beg_ref = -1;
                cigar.beg_query = -1;
                cigar.len = -1;
                return cigar;
            }
        }
        else if (ZERO == where) {
            break;
        }
        else {
            return cigar;
        }
    }

    /* in case we missed the last write */
    WRITE_ANY;

    cigar_reverse = reverse_uint32_t(cigar.seq, cigar.len);
    cigar.seq = cigar_reverse;
    cigar.beg_query = i+1;
    cigar.beg_ref = j+1;

    return cigar;
#undef WRITE
#undef WRITE_ANY
#undef INC
#undef RESET
#undef UNUSED
}

