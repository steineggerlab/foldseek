
/* The MIT License
Copyright (c) 2012-1015 Boston College.
Permission is hereby granted, free of charge, to any person obtaining
        a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
        without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
        permit persons to whom the Software is furnished to do so, subject to
the following conditions:
The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
        BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
        ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*
   Written by Michael Farrar, 2006 (alignment), Mengyao Zhao (SSW Library) and Martin Steinegger (change structure add aa composition, profile and AVX2 support).
   Please send bug reports and/or suggestions to martin.steinegger@snu.ac.kr.
*/
#include "Parameters.h"
#include "StructureSmithWaterman.h"
#include "Util.h"
#include "SubstitutionMatrix.h"
#include "Debug.h"
#include "UngappedAlignment.h"
#include "../strucclustutils/EvalueNeuralNet.h"
#include "block_aligner.h"
#include <iostream>

StructureSmithWaterman::StructureSmithWaterman(size_t maxSequenceLength, int aaSize,
                                               bool aaBiasCorrection, float aaBiasCorrectionScale,
                                               SubstitutionMatrix * subAAMat, SubstitutionMatrix * sub3DiMat)
{
    maxSequenceLength += 1;
    this->subMatAA = subAAMat;
    this->subMat3Di = sub3DiMat;
    this->aaBiasCorrection = aaBiasCorrection;
    this->aaBiasCorrectionScale = aaBiasCorrectionScale;
    const int segSize = (maxSequenceLength+7)/8;
    vHStore = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    vHLoad  = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    vE      = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    vHmax   = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile = new s_profile();
    profile->profile_aa_byte = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    profile->profile_aa_word = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    profile->profile_aa_rev_byte = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    profile->profile_aa_rev_word = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    profile->profile_3di_byte = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    profile->profile_3di_word = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    profile->profile_3di_rev_byte = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    profile->profile_3di_rev_word = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    // gap penalties
#ifdef GAP_POS_SCORING
    profile->profile_gDelOpen_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelOpen_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelClose_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelClose_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gIns_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gIns_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelOpen_rev_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelOpen_rev_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelClose_rev_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelClose_rev_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gIns_rev_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gIns_rev_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->gDelOpen = new uint8_t[maxSequenceLength];
    profile->gDelClose = new uint8_t[maxSequenceLength];
    profile->gDelOpen_rev = new uint8_t[maxSequenceLength];
    profile->gDelClose_rev = new uint8_t[maxSequenceLength];
    profile->gIns_rev = new uint8_t[maxSequenceLength];
#endif
    // remaining
    profile->query_aa_rev_sequence = new int8_t[maxSequenceLength];
    profile->query_aa_sequence     = new int8_t[maxSequenceLength];
    profile->query_3di_rev_sequence = new int8_t[maxSequenceLength];
    profile->query_3di_sequence     = new int8_t[maxSequenceLength];
    profile->composition_bias_aa   = new int8_t[maxSequenceLength];
    profile->composition_bias_ss   = new int8_t[maxSequenceLength];
    profile->composition_bias_aa_rev   = new int8_t[maxSequenceLength];
    profile->composition_bias_ss_rev   = new int8_t[maxSequenceLength];
    profile->profile_aa_word_linear = new short*[aaSize];
    profile_aa_word_linear_data = new short[aaSize*maxSequenceLength];
    profile->profile_3di_word_linear = new short*[aaSize];
    profile_3di_word_linear_data = new short[aaSize*maxSequenceLength];
    profile->mat_rev            = new int8_t[std::max(maxSequenceLength, (size_t)aaSize) * aaSize * 2];
    // this is the same alphabet size for both, please change this todo
    profile->mat_aa                = new int8_t[std::max(maxSequenceLength, (size_t)aaSize) * aaSize * 2];
    profile->mat_3di                = new int8_t[std::max(maxSequenceLength, (size_t)aaSize) * aaSize * 2];
    profile->rev_alignment_aa_profile = new int8_t[maxSequenceLength * Sequence::PROFILE_AA_SIZE];
    profile->rev_alignment_3di_profile = new int8_t[maxSequenceLength * Sequence::PROFILE_AA_SIZE];
    profile->alignment_aa_profile = new int8_t[maxSequenceLength * Sequence::PROFILE_AA_SIZE];
    profile->alignment_3di_profile = new int8_t[maxSequenceLength * Sequence::PROFILE_AA_SIZE];
    tmp_composition_bias   = new float[maxSequenceLength];
    /* array to record the largest score of each reference position */
    maxColumn = new uint8_t[maxSequenceLength*sizeof(uint16_t)];
    memset(maxColumn, 0, maxSequenceLength*sizeof(uint16_t));
    memset(profile->query_aa_sequence, 0, maxSequenceLength * sizeof(int8_t));
    memset(profile->query_aa_rev_sequence, 0, maxSequenceLength * sizeof(int8_t));
    memset(profile->query_3di_sequence, 0, maxSequenceLength * sizeof(int8_t));
    memset(profile->query_3di_rev_sequence, 0, maxSequenceLength * sizeof(int8_t));
    memset(profile->mat_rev, 0, maxSequenceLength * aaSize);
    memset(profile->composition_bias_aa, 0, maxSequenceLength * sizeof(int8_t));
    memset(profile->composition_bias_ss, 0, maxSequenceLength * sizeof(int8_t));
    memset(profile->composition_bias_aa_rev, 0, maxSequenceLength * sizeof(int8_t));
    memset(profile->composition_bias_ss_rev, 0, maxSequenceLength * sizeof(int8_t));
    block = block_new_aa_trace_xdrop(maxSequenceLength, maxSequenceLength, 4096);
}

StructureSmithWaterman::~StructureSmithWaterman(){
    free(vHStore);
    free(vHLoad);
    free(vE);
    free(vHmax);
    free(profile->profile_aa_byte);
    free(profile->profile_aa_word);
    free(profile->profile_aa_rev_byte);
    free(profile->profile_aa_rev_word);
    free(profile->profile_3di_byte);
    free(profile->profile_3di_word);
    free(profile->profile_3di_rev_byte);
    free(profile->profile_3di_rev_word);
#ifdef GAP_POS_SCORING
    free(profile->profile_gDelOpen_byte);
    free(profile->profile_gDelOpen_word);
    free(profile->profile_gDelClose_byte);
    free(profile->profile_gDelClose_word);
    free(profile->profile_gIns_byte);
    free(profile->profile_gIns_word);
    free(profile->profile_gDelOpen_rev_byte);
    free(profile->profile_gDelOpen_rev_word);
    free(profile->profile_gDelClose_rev_byte);
    free(profile->profile_gDelClose_rev_word);
    free(profile->profile_gIns_rev_byte);
    free(profile->profile_gIns_rev_word);
    delete[] profile->gDelOpen;
    delete[] profile->gDelClose;
    delete[] profile->gDelOpen_rev;
    delete[] profile->gDelClose_rev;
    delete[] profile->gIns_rev;
#endif
    delete [] profile->query_aa_rev_sequence;
    delete [] profile->query_aa_sequence;
    delete [] profile->query_3di_rev_sequence;
    delete [] profile->query_3di_sequence;
    delete [] profile->composition_bias_aa;
    delete [] profile->composition_bias_ss;
    delete [] profile->composition_bias_aa_rev;
    delete [] profile->composition_bias_ss_rev;
    delete [] profile->profile_aa_word_linear;
    delete [] profile_aa_word_linear_data;
    delete [] profile->profile_3di_word_linear;
    delete [] profile_3di_word_linear_data;
    delete [] profile->mat_rev;
    delete [] profile->mat_aa;
    delete [] profile->mat_3di;
    delete [] profile->rev_alignment_aa_profile;
    delete [] profile->rev_alignment_3di_profile;
    delete [] profile->alignment_aa_profile;
    delete [] profile->alignment_3di_profile;
    delete [] tmp_composition_bias;
    delete [] maxColumn;
    delete profile;
    block_free_aa_trace_xdrop(block);
}


/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch. */
template <typename T, size_t Elements, const unsigned int type>
void StructureSmithWaterman::createQueryProfile(simd_int *profile, const int8_t *query_sequence, const int8_t * composition_bias, const int8_t *mat,
                                                const int32_t query_length, const int32_t aaSize, uint8_t bias,
                                                const int32_t offset, const int32_t entryLength) {

    const int32_t segLen = (query_length+Elements-1)/Elements;
    T* t = (T*)profile;

    /* Generate query profile rearrange query sequence & calculate the weight of match/mismatch */
    for (int32_t nt = 0; LIKELY(nt < aaSize); nt++) {
//		printf("{");
        for (int32_t i = 0; i < segLen; i ++) {
            int32_t  j = i;
//			printf("(");
            for (size_t segNum = 0; LIKELY(segNum < Elements) ; segNum ++) {
                // if will be optmized out by compiler
                if(type == SUBSTITUTIONMATRIX) {     // substitution score for query_seq constrained by nt
                    // query_sequence starts from 1 to n
                    *t++ = ( j >= query_length) ? bias : mat[nt * aaSize + query_sequence[j + offset ]] + composition_bias[j + offset] + bias; // mat[nt][q[j]] mat eq 20*20
//					printf("(%1d, %1d) ", query_sequence[j ], *(t-1));

                } if(type == PROFILE || type == PROFILE_HMM) {
                    // profile starts by 0
                    *t++ = (j >= query_length) ? bias : mat[nt * entryLength + j + offset] + bias;
//					printf("(%1d, %1d) ", j , *(t-1));
                }
                j += segLen;
            }
//			printf(")");
        }
//		printf("}\n");
    }
//	printf("\n");
//	std::flush(std::cout);

}

#ifdef GAP_POS_SCORING
template <typename T, size_t Elements>
void StructureSmithWaterman::createGapProfile(simd_int* profile_gDelOpen, simd_int* profile_gDelClose, simd_int* profile_gIns,
                                              const uint8_t* gDelOpen, const uint8_t* gDelClose, const uint8_t* gIns,
                                              const int32_t query_length, const int32_t offset) {
    const int32_t segLen = (query_length - offset + Elements - 1) / Elements;
    T* delOpen = (T*) profile_gDelOpen;
    T* delClose = (T*) profile_gDelClose;
    T* ins = (T*) profile_gIns;
    for (int32_t i = 0; LIKELY(i < segLen); ++i) {
        int32_t j = i;
        for (size_t segNum = 0; LIKELY(segNum < Elements); ++segNum) {
            *delOpen++ = (j < query_length) ? gDelOpen[j + offset + 1] : 0; // offset + 1 because it calculates F for the next column
            *delClose++ = (j < query_length) ? gDelClose[j + offset + 1] : 0;
            *ins++ = (j < query_length) ? gIns[j + offset] : 0;
            j += segLen;
        }
    }
}
#endif

template <unsigned int profile_type>
StructureSmithWaterman::s_align StructureSmithWaterman::alignScoreEndPos (
        const unsigned char *db_aa_sequence,
        const unsigned char *db_3di_sequence,
        int32_t db_length,
        const uint8_t gap_open,
        const uint8_t gap_extend,
        const int32_t maskLen) {
    int32_t query_length = profile->query_length;
    s_align r;
    r.word = 1;
    r.dbStartPos1 = -1;
    r.qStartPos1 = -1;
    r.cigar = 0;
    r.cigarLen = 0;
    //if (maskLen < 15) {
    //	fprintf(stderr, "When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.\n");
    //}

    std::pair<alignment_end, alignment_end> bests;
    // Find the alignment scores and ending positions
//    if(profile->isProfile){
//        bests = sw_sse2_byte<profile_type>(db_aa_sequence, db_3di_sequence, 0, db_length, query_length,
//                                          gap_open, gap_extend, profile->profile_gDelOpen_byte, profile->profile_gDelClose_byte,
//                                          profile->profile_gIns_byte, profile->profile_aa_byte, profile->profile_3di_byte,UCHAR_MAX, profile->bias, maskLen);
//
//    }else{
//        bests = sw_sse2_byte<SUBSTITUTIONMATRIX>(db_aa_sequence, db_3di_sequence, 0, db_length, query_length,
//                                                 gap_open, gap_extend, NULL, NULL, NULL, profile->profile_aa_byte, profile->profile_3di_byte,UCHAR_MAX, profile->bias, maskLen);
//    }
//    if (bests.first.score == 255) {
        if(profile->isProfile) {
            bests = sw_sse2_word<profile_type>(db_aa_sequence, db_3di_sequence, 0, db_length, query_length,
                                              gap_open, gap_extend,
#ifdef GAP_POS_SCORING
                                              profile->profile_gDelOpen_word, profile->profile_gDelClose_word, profile->profile_gIns_word,
#endif
                                              profile->profile_aa_word,
                                              profile->profile_3di_word, USHRT_MAX, maskLen);

        } else {
            bests = sw_sse2_word<SUBSTITUTIONMATRIX>(db_aa_sequence, db_3di_sequence, 0, db_length, query_length,
                                                     gap_open, gap_extend,
#ifdef GAP_POS_SCORING
                                                     NULL, NULL, NULL,
#endif
                                                     profile->profile_aa_word,
                                                     profile->profile_3di_word, USHRT_MAX, maskLen);
        }
        r.word = 1;
//    } else if (bests.first.score == 255) {
//        fprintf(stderr, "Please set 2 to the score_size parameter of the function ssw_init, otherwise the alignment results will be incorrect.\n");
//        EXIT(EXIT_FAILURE);
//    }

    r.score1 = bests.first.score;
    r.dbEndPos1 = bests.first.ref;
    r.qEndPos1 = bests.first.read;

    if (maskLen >= 15) {
        r.score2 = bests.second.score;
        r.ref_end2 = bests.second.ref;
    } else {
        r.score2 = 0;
        r.ref_end2 = -1;
    }

    // no residue could be aligned
    if (r.dbEndPos1 == -1) {
        return r;
    }
    r.qCov = computeCov(0, r.qEndPos1, query_length);
    r.tCov = computeCov(0, r.dbEndPos1, db_length);

    return r;
}

template
StructureSmithWaterman::s_align StructureSmithWaterman::alignScoreEndPos<StructureSmithWaterman::PROFILE>(const unsigned char*, const unsigned char*, int32_t, const uint8_t, const uint8_t, const int32_t);
template
StructureSmithWaterman::s_align StructureSmithWaterman::alignScoreEndPos<StructureSmithWaterman::PROFILE_HMM>(const unsigned char*, const unsigned char*, int32_t, const uint8_t, const uint8_t, const int32_t);

StructureSmithWaterman::s_align StructureSmithWaterman::alignStartPosBacktraceBlock(
        const unsigned char *db_aa_sequence,
        const unsigned char *db_3di_sequence,
        int32_t db_length,
        const uint8_t gap_open,
        const uint8_t gap_extend,
        std::string & backtrace,
        StructureSmithWaterman::s_align r) {
#define MAX_SIZE 4096 //TODO
    size_t query_len = profile->query_length;
    size_t target_len = db_length;
    Gaps gaps;
    gaps.open   = -gap_open;
    gaps.extend = -gap_extend;
    int32_t target_score = r.score1;

    // note: instead of query_len or target_len, it is possible to use query_aa really large length
    // and reuse data structures to avoid allocations
    PaddedBytes* query_aa = block_new_padded_aa(query_len, MAX_SIZE);
    PaddedBytes* query_3di = block_new_padded_aa(query_len, MAX_SIZE);
    PosBias* query_bias = block_new_pos_bias(query_len, MAX_SIZE);
    PaddedBytes* target_aa = block_new_padded_aa(target_len, MAX_SIZE);
    PaddedBytes* target_3di = block_new_padded_aa(target_len, MAX_SIZE);
    PosBias* target_bias = block_new_pos_bias(target_len, MAX_SIZE);

    int32_t queryStartPos = query_len - (r.qEndPos1 + 1);
    int32_t queryAlnLen = r.qEndPos1 + 1;
    // convert query_aa_rev_sequence to ascii
    std::string query_aa_sequence_str;
    std::string query_3di_sequence_str;
    int16_t * query_bias_arr = new int16_t[queryAlnLen];
    for(int i = 0; i < queryAlnLen; i++){
        query_aa_sequence_str.push_back(subMatAA->num2aa[profile->query_aa_rev_sequence[queryStartPos + i]]);
        query_3di_sequence_str.push_back(subMat3Di->num2aa[profile->query_3di_rev_sequence[queryStartPos + i]]);
        query_bias_arr[i] =  profile->composition_bias_aa_rev[queryStartPos + i] +
                             profile->composition_bias_ss_rev[queryStartPos + i];
    }


    block_set_bytes_padded_aa(query_aa,  (const uint8_t*) query_aa_sequence_str.data(), queryAlnLen, MAX_SIZE);
    block_set_bytes_padded_aa(query_3di, (const uint8_t*) query_3di_sequence_str.data(), queryAlnLen, MAX_SIZE);

    block_set_pos_bias(query_bias, query_bias_arr, queryAlnLen);

    int32_t targetAlnLen = r.dbEndPos1 + 1;
    std::string db_aa_sequence_str;
    std::string db_3di_sequence_str;
    // copy this db_aa_sequence,db_aa_sequence + r.dbEndPos1 + 1 in reverse order to db_aa_sequence_str and mappping to ascii using subMatAA->num2aa
    for(int i = targetAlnLen - 1; i >= 0; i--){
        db_aa_sequence_str.push_back(subMatAA->num2aa[db_aa_sequence[i]]);
        db_3di_sequence_str.push_back(subMat3Di->num2aa[db_3di_sequence[i]]);
    }
    block_set_bytes_padded_aa(target_aa, (const uint8_t*) db_aa_sequence_str.data(), targetAlnLen, MAX_SIZE);
    block_set_bytes_padded_aa(target_3di, (const uint8_t*)db_3di_sequence_str.data(), targetAlnLen, MAX_SIZE);
    int16_t * target_bias_arr = new int16_t[targetAlnLen];
    memset(target_bias_arr, 0, targetAlnLen * sizeof(int16_t));
    block_set_pos_bias(target_bias, target_bias_arr, targetAlnLen);


    AAMatrix* matrix_aa = block_new_simple_aamatrix(1, -1);
    for (int aa1 = 0; aa1 < subMatAA->alphabetSize; aa1++) {
        for (int aa2 = 0; aa2 < subMatAA->alphabetSize; aa2++) {
            // set to actual scores instead of zeros!
            block_set_aamatrix(matrix_aa, subMatAA->num2aa[aa1], subMatAA->num2aa[aa2],
                               subMatAA->subMatrix[aa1][aa2]);
        }
    }

    AAMatrix* matrix_3di = block_new_simple_aamatrix(1, -1);
    for (int aa1 = 0; aa1 < subMat3Di->alphabetSize; aa1++) {

        for (int aa2 = 0; aa2 < subMat3Di->alphabetSize; aa2++) {
            // set to actual scores instead of zeros!
            block_set_aamatrix(matrix_3di, subMat3Di->num2aa[aa1], subMat3Di->num2aa[aa2],
                               subMat3Di->subMatrix[aa1][aa2]);

        }

    }

    Cigar* cigar = block_new_cigar(queryAlnLen, targetAlnLen);

    AlignResult res;
    size_t min_size = 32;
    res.score = -1000000000;
    res.query_idx = -1;
    res.reference_idx = -1;

    // exponential search on min_size until either max_size is reached or target_score is reached
    while (min_size <= MAX_SIZE && res.score < target_score) {
        // allow max block size to grow
        SizeRange range;
        range.min = min_size;
        range.max = MAX_SIZE;
        // estimated x-drop threshold
        int32_t x_drop = -(min_size * gaps.extend + gaps.open);
        block_align_3di_aa_trace_xdrop(block, query_aa, query_3di, query_bias, target_aa, target_3di, target_bias,
                                       matrix_aa, matrix_3di, gaps, range, x_drop);
        res = block_res_aa_trace_xdrop(block);
        min_size *= 2;
    }

    if (res.score != target_score && !(target_score == INT16_MAX && res.score >= target_score)) {
        printf("ERROR: target_score not reached. res.score: %d target_score: %d", res.score, target_score);
        exit(1); // TODO
    }

    block_cigar_aa_trace_xdrop(block, res.query_idx, res.reference_idx, cigar);
//    printf("query_aa: %s\nquery_3di: %s\ntarget_aa: %s\ntarget_3di: %s\nscore: %d\nidx: (%lu, %lu)\n",
//           profile->query_aa_rev_sequence,
//           profile->query_3di_rev_sequence,
//           db_aa_sequence,
//           db_3di_sequence,
//           res.score,
//           res.query_idx,
//           res.reference_idx);


    size_t cigar_len = block_len_cigar(cigar);

    // Note: 'M' signals either query_aa match or mismatch
    uint32_t aaIds = 0;
    size_t queryPos = 0;
    size_t targetPos = 0;
    for (size_t i = 0; i < cigar_len; i++) {
        OpLen o = block_get_cigar(cigar, i);
        //printf("%lu%c", o.len, ops_char[o.op]);
        if(o.op == 1){
            for(size_t j = 0; j < o.len; j++){
                if(query_aa_sequence_str[queryPos + j] == db_aa_sequence_str[targetPos + j]){
                    aaIds++;
                }
            }
            queryPos += o.len;
            targetPos += o.len;
            backtrace.append(o.len,'M');
        }else if(o.op == 4){
            queryPos += o.len;
            backtrace.append(o.len,'I');
        }else if(o.op == 5){
            targetPos += o.len;
            backtrace.append(o.len,'D');
        }
    }
    r.identicalAACnt = aaIds;
    //reverse backtrace
    std::reverse(backtrace.begin(), backtrace.end());
    r.qStartPos1 = (r.qEndPos1 + 1) - queryPos;
    r.dbStartPos1 = (r.dbEndPos1 + 1) - targetPos;

    r.qCov = computeCov(r.qStartPos1, r.qEndPos1, query_len);
    r.tCov = computeCov(r.dbStartPos1, r.dbEndPos1, db_length);

    block_free_cigar(cigar);
    block_free_padded_aa(query_aa);
    block_free_padded_aa(query_3di);
    block_free_pos_bias(query_bias);
    block_free_padded_aa(target_aa);
    block_free_padded_aa(target_3di);
    block_free_pos_bias(target_bias);
    block_free_aamatrix(matrix_3di);
    block_free_aamatrix(matrix_aa);
    delete [] query_bias_arr;
    delete [] target_bias_arr;
    return r;
}


template <unsigned int profile_type>
StructureSmithWaterman::s_align StructureSmithWaterman::alignStartPosBacktrace (
        const unsigned char *db_aa_sequence,
        const unsigned char *db_3di_sequence,
        int32_t db_length,
        const uint8_t gap_open,
        const uint8_t gap_extend,
        const uint8_t alignmentMode,	//  (from high to low) bit 5: return the best alignment beginning position; 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar; 7: if max score >= filters, return cigar; 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
        std::string & backtrace,
        StructureSmithWaterman::s_align r,
        const int covMode, const float covThr,
        const int32_t maskLen) {
    int32_t query_length = profile->query_length;
    int32_t queryOffset = query_length - r.qEndPos1 - 1;
    std::pair<alignment_end, alignment_end> bests_reverse;
    // Find the beginning position of the best alignment.
    if (r.word == 0) {

        if (profile->isProfile) {
            createQueryProfile<int8_t, VECSIZE_INT * 4, profile_type>(profile->profile_aa_rev_byte, profile->query_aa_rev_sequence, NULL, profile->rev_alignment_aa_profile,
                                                                     r.qEndPos1 + 1, profile->alphabetSize, profile->bias, queryOffset, profile->query_length);
            createQueryProfile<int8_t, VECSIZE_INT * 4, profile_type>(profile->profile_3di_rev_byte, profile->query_3di_rev_sequence, NULL, profile->rev_alignment_3di_profile,
                                                                     r.qEndPos1 + 1, profile->alphabetSize, profile->bias, queryOffset, profile->query_length);
#ifdef GAP_POS_SCORING
            createGapProfile<int8_t, VECSIZE_INT * 4>(profile->profile_gDelOpen_rev_byte, profile->profile_gDelClose_rev_byte, profile->profile_gIns_rev_byte,
                                                      profile->gDelOpen_rev, profile->gDelClose_rev, profile->gIns_rev, profile->query_length, queryOffset);
#endif
            bests_reverse = sw_sse2_byte<profile_type>(db_aa_sequence, db_3di_sequence, 1, r.dbEndPos1 + 1, r.qEndPos1 + 1,
                                                      gap_open, gap_extend,
#ifdef GAP_POS_SCORING
                                                      profile->profile_gDelOpen_rev_byte, profile->profile_gDelClose_rev_byte, profile->profile_gIns_rev_byte, 
#endif
                                                      profile->profile_aa_rev_byte, profile->profile_3di_rev_byte,
                                                      r.score1, profile->bias, maskLen);
        }else{
            createQueryProfile<int8_t, VECSIZE_INT * 4, SUBSTITUTIONMATRIX>(profile->profile_aa_rev_byte, profile->query_aa_rev_sequence, profile->composition_bias_aa_rev, profile->mat_aa,
                                                                            r.qEndPos1 + 1, profile->alphabetSize, profile->bias, queryOffset, 0);
            createQueryProfile<int8_t, VECSIZE_INT * 4, SUBSTITUTIONMATRIX>(profile->profile_3di_rev_byte, profile->query_3di_rev_sequence, profile->composition_bias_ss_rev, profile->mat_3di,
                                                                            r.qEndPos1 + 1, profile->alphabetSize, profile->bias, queryOffset, 0);
            bests_reverse = sw_sse2_byte<SUBSTITUTIONMATRIX>(db_aa_sequence, db_3di_sequence, 1, r.dbEndPos1 + 1, r.qEndPos1 + 1,
                                                             gap_open, gap_extend,
#ifdef GAP_POS_SCORING
                                                             NULL, NULL, NULL,
#endif
                                                             profile->profile_aa_rev_byte,profile->profile_3di_rev_byte,
                                                             r.score1, profile->bias, maskLen);
        }
    } else {
        if (profile->isProfile) {
            createQueryProfile<int16_t, VECSIZE_INT * 2, profile_type>(profile->profile_aa_rev_word,profile->query_aa_rev_sequence,NULL,profile->rev_alignment_aa_profile,
                                                                      r.qEndPos1 + 1, profile->alphabetSize, 0,queryOffset, profile->query_length);
            createQueryProfile<int16_t, VECSIZE_INT * 2, profile_type>(profile->profile_3di_rev_word,profile->query_3di_rev_sequence,NULL,profile->rev_alignment_3di_profile,
                                                                      r.qEndPos1 + 1, profile->alphabetSize, 0,queryOffset, profile->query_length);
#ifdef GAP_POS_SCORING
            createGapProfile<int16_t, VECSIZE_INT * 2>(profile->profile_gDelOpen_rev_word,
                                                       profile->profile_gDelClose_rev_word,
                                                       profile->profile_gIns_rev_word, profile->gDelOpen_rev,
                                                       profile->gDelClose_rev, profile->gIns_rev,
                                                       profile->query_length, queryOffset);
#endif
            bests_reverse = sw_sse2_word<profile_type>(db_aa_sequence, db_3di_sequence, 1, r.dbEndPos1 + 1, r.qEndPos1 + 1,
                                                      gap_open, gap_extend,
#ifdef GAP_POS_SCORING
                                                      profile->profile_gDelOpen_rev_word, profile->profile_gDelClose_rev_word, profile->profile_gIns_rev_word,
#endif
                                                      profile->profile_aa_rev_word, profile->profile_3di_rev_word,
                                                      r.score1, maskLen);
        }else{
            createQueryProfile<int16_t, VECSIZE_INT * 2, SUBSTITUTIONMATRIX>(profile->profile_aa_rev_word, profile->query_aa_rev_sequence, profile->composition_bias_aa_rev, profile->mat_aa,
                                                                             r.qEndPos1 + 1, profile->alphabetSize, 0, queryOffset, 0);
            createQueryProfile<int16_t, VECSIZE_INT * 2, SUBSTITUTIONMATRIX>(profile->profile_3di_rev_word, profile->query_3di_rev_sequence, profile->composition_bias_ss_rev, profile->mat_3di,
                                                                             r.qEndPos1 + 1, profile->alphabetSize, 0, queryOffset, 0);
            bests_reverse = sw_sse2_word<SUBSTITUTIONMATRIX>(db_aa_sequence, db_3di_sequence, 1, r.dbEndPos1 + 1, r.qEndPos1 + 1,
                                                             gap_open, gap_extend,
#ifdef GAP_POS_SCORING
                                                             NULL, NULL, NULL,
#endif
                                                             profile->profile_aa_rev_word, profile->profile_3di_rev_word,
                                                             r.score1, maskLen);
        }

    }
    if(bests_reverse.first.score != r.score1){
        fprintf(stderr, "Score of forward/backward SW differ. This should not happen.\n");
        EXIT(EXIT_FAILURE);
    }

    r.dbStartPos1 = bests_reverse.first.ref;
    r.qStartPos1 = r.qEndPos1 - bests_reverse.first.read;

    if (r.dbStartPos1 == -1) {
        fprintf(stderr, "Target start position is -1. This should not happen.\n");
        EXIT(EXIT_FAILURE);
    }

    r.qCov = computeCov(r.qStartPos1, r.qEndPos1, query_length);
    r.tCov = computeCov(r.dbStartPos1, r.dbEndPos1, db_length);
    bool hasLowerCoverage = !(Util::hasCoverage(covThr, covMode, r.qCov, r.tCov));
    // only start and end point are needed
    if (alignmentMode == 1 || hasLowerCoverage) {
        return r;
    }

    // Generate cigar.
    db_length = r.dbEndPos1 - r.dbStartPos1 + 1;
    query_length = r.qEndPos1 - r.qStartPos1 + 1;
    int32_t band_width = abs(db_length - query_length) + 1;

    cigar* path;
    if (profile->isProfile) {
        path = banded_sw<profile_type>(db_aa_sequence + r.dbStartPos1,
                                      db_3di_sequence + r.dbStartPos1,
                                      profile->query_aa_sequence + r.qStartPos1,
                                      profile->query_3di_sequence + r.qStartPos1,
                                      profile->composition_bias_aa + r.qStartPos1,
                                      profile->composition_bias_ss + r.qStartPos1,
                                      db_length, query_length, r.qStartPos1, r.score1,
                                      gap_open, gap_extend,
#ifdef GAP_POS_SCORING
                                      profile->gDelOpen + r.qStartPos1, profile->gDelClose + r.qStartPos1, profile->gIns + r.qStartPos1,
#endif
                                      band_width,profile->alignment_aa_profile,  profile->query_length,
                                      profile->alignment_3di_profile, profile->query_length);
    } else {
        path = banded_sw<SUBSTITUTIONMATRIX>(db_aa_sequence + r.dbStartPos1,
                                             db_3di_sequence + r.dbStartPos1,
                                             profile->query_aa_sequence + r.qStartPos1,
                                             profile->query_3di_sequence + r.qStartPos1,
                                             profile->composition_bias_aa + r.qStartPos1,
                                             profile->composition_bias_ss + r.qStartPos1,
                                             db_length, query_length, r.qStartPos1, r.score1,
                                             gap_open, gap_extend,
#ifdef GAP_POS_SCORING
                                             NULL, NULL, NULL,
#endif
                                             band_width, profile->mat_aa,  profile->alphabetSize,
                                             profile->mat_3di, profile->alphabetSize);
    }
    if (path != NULL) {
        r.cigar = path->seq;
        r.cigarLen = path->length;
    }

    uint32_t aaIds = 0;
    size_t mStateCnt = 0;
    computerBacktrace(profile, db_aa_sequence, r, backtrace, aaIds,  mStateCnt);
    r.identicalAACnt = aaIds;
    if(path != NULL) {
        delete[] path->seq;
        delete path;
    }
    return r;
}

template
StructureSmithWaterman::s_align StructureSmithWaterman::alignStartPosBacktrace<StructureSmithWaterman::PROFILE>(const unsigned char*, const unsigned char*, int32_t, const uint8_t, const uint8_t, const uint8_t, std::string& , StructureSmithWaterman::s_align, const int, const float, const int32_t);
template
StructureSmithWaterman::s_align StructureSmithWaterman::alignStartPosBacktrace<StructureSmithWaterman::PROFILE_HMM>(const unsigned char*, const unsigned char*, int32_t, const uint8_t, const uint8_t, const uint8_t, std::string& , StructureSmithWaterman::s_align, const int, const float, const int32_t);

void StructureSmithWaterman::computerBacktrace(s_profile * query, const unsigned char * db_aa_sequence,
                                               s_align & alignment, std::string & backtrace,
                                               uint32_t & aaIds, size_t & mStatesCnt){
    int32_t targetPos = alignment.dbStartPos1, queryPos = alignment.qStartPos1;
    for (int32_t c = 0; c < alignment.cigarLen; ++c) {
        char letter = StructureSmithWaterman::cigar_int_to_op(alignment.cigar[c]);
        uint32_t length = StructureSmithWaterman::cigar_int_to_len(alignment.cigar[c]);
        backtrace.reserve(length);
        for (uint32_t i = 0; i < length; ++i){
            if (letter == 'M') {
                aaIds += (db_aa_sequence[targetPos] == query->query_aa_sequence[queryPos]);
                ++mStatesCnt;
                ++queryPos;
                ++targetPos;
                backtrace.append("M");
            } else {
                if (letter == 'I') {
                    ++queryPos;
                    backtrace.append("I");
                }
                else{
                    ++targetPos;
                    backtrace.append("D");
                }
            }
        }
    }
}



char StructureSmithWaterman::cigar_int_to_op(uint32_t cigar_int) {
    uint8_t letter_code = cigar_int & 0xfU;
    static const char map[] = {
            'M',
            'I',
            'D',
            'N',
            'S',
            'H',
            'P',
            '=',
            'X',
    };

    if (letter_code >= (sizeof(map)/sizeof(map[0]))) {
        return 'M';
    }

    return map[letter_code];
}

uint32_t StructureSmithWaterman::cigar_int_to_len (uint32_t cigar_int)
{
    uint32_t res = cigar_int >> 4;
    return res;
}

template <const unsigned int type>
std::pair<StructureSmithWaterman::alignment_end, StructureSmithWaterman::alignment_end> StructureSmithWaterman::sw_sse2_byte (const unsigned char* db_aa_sequence,
                                                                                                                              const unsigned char* db_3di_sequence,
                                                                                                                              int8_t ref_dir,	// 0: forward ref; 1: reverse ref
                                                                                                                              int32_t db_length,
                                                                                                                              int32_t query_length,
                                                                                                                              const uint8_t gap_open, /* will be used as - */
                                                                                                                              const uint8_t gap_extend, /* will be used as - */
#ifdef GAP_POS_SCORING
                                                                                                                              const simd_int *gap_open_del,
                                                                                                                              const simd_int *gap_close_del,
                                                                                                                              const simd_int *gap_open_ins,
#endif
                                                                                                                              const simd_int* query_aa_profile_byte,
                                                                                                                              const simd_int* query_3di_profile_byte,
                                                                                                                              uint8_t terminate,	/* the best alignment score: used to terminate
                                                                                                                                                     the matrix calculation when locating the
                                                                                                                                                     alignment beginning point. If this score
                                                                                                                                                     is set to 0, it will not be used */
                                                                                                                              uint8_t bias,  /* Shift 0 point to a positive value. */
                                                                                                                              int32_t maskLen) {
#define max16(m, vm) ((m) = simdi8_hmax((vm)));

    uint8_t max = 0;		                     /* the max alignment score */
    int32_t end_query = query_length - 1;
    int32_t end_db = -1; /* 0_based best alignment ending point; Initialized as isn't aligned -1. */
    const int SIMD_SIZE = VECSIZE_INT * 4;
    int32_t segLen = (query_length + SIMD_SIZE-1) / SIMD_SIZE; /* number of segment */
    /* array to record the largest score of each reference position */
    memset(this->maxColumn, 0, db_length * sizeof(uint8_t));
    uint8_t * maxColumn = (uint8_t *) this->maxColumn;

    /* Define 16 byte 0 vector. */
    simd_int vZero = simdi32_set(0);
    simd_int* pvHStore = vHStore;
    simd_int* pvHLoad = vHLoad;
    simd_int* pvE = vE;
    simd_int* pvHmax = vHmax;
    memset(pvHStore,0,segLen*sizeof(simd_int));
    memset(pvHLoad,0,segLen*sizeof(simd_int));
    memset(pvE,0,segLen*sizeof(simd_int));
    memset(pvHmax,0,segLen*sizeof(simd_int));

    int32_t i, j;
    /* 16 byte insertion begin vector */
#ifdef GAP_POS_SCORING
    simd_int vGapO;
    if (type != PROFILE_HMM) {
        vGapO = simdi8_set(gap_open);
    }
#else
    simd_int vGapO = simdi8_set(gap_open);
#endif

    /* 16 byte insertion extension vector */
    simd_int vGapE = simdi8_set(gap_extend);

    /* 16 byte bias vector */
    simd_int vBias = simdi8_set(2*bias);

    simd_int vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
    simd_int vMaxMark = vZero; /* Trace the highest score till the previous column. */
    simd_int vTemp;
    int32_t edge, begin = 0, end = db_length, step = 1;
    //	int32_t distance = query_length * 2 / 3;
    //	int32_t distance = query_length / 2;
    //	int32_t distance = query_length;

    /* outer loop to process the reference sequence */
    if (ref_dir == 1) {
        begin = db_length - 1;
        end = -1;
        step = -1;
    }
    for (i = begin; LIKELY(i != end); i += step) {
        simd_int e, vF = vZero, vMaxColumn = vZero; /* Initialize F value to 0.
                                                    Any errors to vH values will be corrected in the Lazy_F loop.
                                                    */
        //		max16(maxColumn[i], vMaxColumn);
        //		fprintf(stderr, "middle[%d]: %d\n", i, maxColumn[i]);

        simd_int vH = pvHStore[segLen - 1];
        vH = simdi8_shiftl (vH, 1); /* Shift the 128-bit value in vH left by 1 byte. */
        const simd_int* vPAA = query_aa_profile_byte + db_aa_sequence[i] * segLen; /* Right part of the query_profile_byte */
        const simd_int* vP3Di = query_3di_profile_byte + db_3di_sequence[i] * segLen; /* Right part of the query_profile_byte */
        //	int8_t* t;
        //	int32_t ti;
        //        fprintf(stderr, "i: %d of %d:\t ", i,segLen);
        //for (t = (int8_t*)vP, ti = 0; ti < segLen; ++ti) fprintf(stderr, "%d\t", *t++);
        //fprintf(stderr, "\n");

        /* Swap the 2 H buffers. */
        simd_int* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (j = 0; LIKELY(j < segLen); ++j) {
            simd_int score = simdui8_adds(simdi_load(vPAA + j), simdi_load(vP3Di + j));
//            int8_t* t;
//            int32_t ti;
//            fprintf(stderr, "score: ");
//            for (t = (int8_t*)&score, ti = 0; ti < 16; ++ti) fprintf(stderr, "%d\t", (*t++) - 2*bias);
//            fprintf(stderr, "\n");

            vH = simdui8_adds(vH, score);
            vH = simdui8_subs(vH, vBias); /* vH will be always > 0 */
            //	max16(maxColumn[i], vH);
            //	fprintf(stderr, "H[%d]: %d\n", i, maxColumn[i]);
//            fprintf(stderr, "vh: ");
//            for (t = (int8_t*)&vH, ti = 0; ti < 16; ++ti) fprintf(stderr, "%d\t", *t++);
//            fprintf(stderr, "\n");
            /* Get max from vH, vE and vF. */
            e = simdi_load(pvE + j);
            vH = simdui8_max(vH, e);
#ifdef GAP_POS_SCORING
            if (type == PROFILE_HMM) {
                vH = simdui8_max(vH, simdui8_subs(vF, simdi_load(gap_close_del + j)));
            } else {
#endif
                vH = simdui8_max(vH, vF);
#ifdef GAP_POS_SCORING
            }
#endif
            vMaxColumn = simdui8_max(vMaxColumn, vH);

            //	max16(maxColumn[i], vMaxColumn);
            //	fprintf(stderr, "middle[%d]: %d\n", i, maxColumn[i]);
            //	for (t = (int8_t*)&vMaxColumn, ti = 0; ti < 16; ++ti) fprintf(stderr, "%d\t", *t++);

            /* Save vH values. */
            simdi_store(pvHStore + j, vH);

            /* Update vE value. */
#ifdef GAP_POS_SCORING
            if (type == PROFILE_HMM) {
                // copy vH for update of vF
                vTemp = vH;
                vH = simdui8_subs(vH, simdi_load(gap_open_ins + j)); /* saturation arithmetic, result >= 0 */
            } else {
#endif
                vH = simdui8_subs(vH, vGapO); /* saturation arithmetic, result >= 0 */
#ifdef GAP_POS_SCORING
            }
#endif
            e = simdui8_subs(e, vGapE);
            e = simdui8_max(e, vH);
            simdi_store(pvE + j, e);

            /* Update vF value. */
            vF = simdui8_subs(vF, vGapE);
#ifdef GAP_POS_SCORING
            if (type == PROFILE_HMM) {
                vF = simdui8_max(vF, simdui8_subs(vTemp, simdi_load(gap_open_del + j)));
            } else {
#endif
                vF = simdui8_max(vF, vH);
#ifdef GAP_POS_SCORING
            }
#endif

            /* Load the next vH. */
            vH = simdi_load(pvHLoad + j);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
        /* reset pointers to the start of the saved data */
        j = 0;
        vH = simdi_load (pvHStore + j);

        /*  the computed vF value is for the given column.  since */
        /*  we are at the end, we need to shift the vF value over */
        /*  to the next column. */
        vF = simdi8_shiftl (vF, 1);
#ifdef GAP_POS_SCORING
        if (type == PROFILE_HMM) {
            vTemp = simdui8_subs(vH, simdi_load(gap_open_del + j));
        } else {
#endif
            vTemp = simdui8_subs(vH, vGapO);
#ifdef GAP_POS_SCORING
        }
#endif
        vTemp = simdui8_subs (vF, vTemp);
        vTemp = simdi8_eq (vTemp, vZero);
        uint32_t cmp = simdi8_movemask (vTemp);
        while (cmp != SIMD_MOVEMASK_MAX) {
#ifdef GAP_POS_SCORING
            if (type == PROFILE_HMM) {
                vH = simdui8_max (vH, simdui8_subs(vF, simdi_load(gap_close_del + j)));
                simdi_store(pvE + j, simdui8_max(simdi_load(pvE + j), simdui8_subs(vH, simdi_load(gap_open_ins + j))));
            } else {
#endif
                vH = simdui8_max (vH, vF);
#ifdef GAP_POS_SCORING
            }
#endif
            vMaxColumn = simdui8_max(vMaxColumn, vH);
            simdi_store (pvHStore + j, vH);
            vF = simdui8_subs (vF, vGapE);
            j++;
            if (j >= segLen)
            {
                j = 0;
                vF = simdi8_shiftl (vF, 1);
            }
            vH = simdi_load (pvHStore + j);
#ifdef GAP_POS_SCORING
            if (type == PROFILE_HMM) {
                vTemp = simdui8_subs(vH, simdi_load(gap_open_del + j));
            } else {
#endif
                vTemp = simdui8_subs(vH, vGapO);
#ifdef GAP_POS_SCORING
            }
#endif
            vTemp = simdui8_subs (vF, vTemp);
            vTemp = simdi8_eq (vTemp, vZero);
            cmp  = simdi8_movemask (vTemp);
        }

        vMaxScore = simdui8_max(vMaxScore, vMaxColumn);
        vTemp = simdi8_eq(vMaxMark, vMaxScore);
        cmp = simdi8_movemask(vTemp);
        if (cmp != SIMD_MOVEMASK_MAX) {
            uint8_t temp;
            vMaxMark = vMaxScore;
            max16(temp, vMaxScore);
            vMaxScore = vMaxMark;

            if (LIKELY(temp > max)) {
                max = temp;
                if (max + bias >= 255) break;	//overflow
                end_db = i;

                /* Store the column with the highest alignment score in order to trace the alignment ending position on read. */
                for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
            }
        }

        /* Record the max score of current column. */
        max16(maxColumn[i], vMaxColumn);
        //		fprintf(stderr, "maxColumn[%d]: %d\n", i, maxColumn[i]);
        if (maxColumn[i] == terminate) break;
    }

    /* Trace the alignment ending position on read. */
    uint8_t *t = (uint8_t*)pvHmax;
    int32_t column_len = segLen * SIMD_SIZE;
    for (i = 0; LIKELY(i < column_len); ++i, ++t) {
        int32_t temp;
        if (*t == max) {
            temp = i / SIMD_SIZE + i % SIMD_SIZE * segLen;
            if (temp < end_query) end_query = temp;
        }
    }

    /* Find the most possible 2nd best alignment. */
    alignment_end best0;
    best0.score = max + (2 * bias) >= 255 ? 255 : max;
    best0.ref = end_db;
    best0.read = end_query;

    alignment_end best1;
    best1.score = 0;
    best1.ref = 0;
    best1.read = 0;

    edge = (end_db - maskLen) > 0 ? (end_db - maskLen) : 0;
    for (i = 0; i < edge; i ++) {
        //			fprintf (stderr, "maxColumn[%d]: %d\n", i, maxColumn[i]);
        if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
        }
    }
    edge = (end_db + maskLen) > db_length ? db_length : (end_db + maskLen);
    for (i = edge + 1; i < db_length; i ++) {
        //			fprintf (stderr, "db_length: %d\tmaxColumn[%d]: %d\n", db_length, i, maxColumn[i]);
        if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
        }
    }

    return std::make_pair(best0, best1);
#undef max16
}

template <const unsigned int type>
std::pair<StructureSmithWaterman::alignment_end, StructureSmithWaterman::alignment_end> StructureSmithWaterman::sw_sse2_word (const unsigned char* db_aa_sequence,
                                                                                                                              const unsigned char* db_3di_sequence,
                                                                                                                              int8_t ref_dir,	// 0: forward ref; 1: reverse ref
                                                                                                                              int32_t db_length,
                                                                                                                              int32_t query_length,
                                                                                                                              const uint8_t gap_open, /* will be used as - */
                                                                                                                              const uint8_t gap_extend, /* will be used as - */
#ifdef GAP_POS_SCORING
                                                                                                                              const simd_int *gap_open_del,
                                                                                                                              const simd_int *gap_close_del,
                                                                                                                              const simd_int *gap_open_ins,
#endif
                                                                                                                              const simd_int*query_aa_profile_word,
                                                                                                                              const simd_int*query_3di_profile_word,
                                                                                                                              uint16_t terminate,
                                                                                                                              int32_t maskLen) {

#define max8(m, vm) ((m) = simdi16_hmax((vm)));

    uint16_t max = 0;		                     /* the max alignment score */
    int32_t end_read = query_length - 1;
    int32_t end_ref = 0; /* 1_based best alignment ending point; Initialized as isn't aligned - 0. */
    const unsigned int SIMD_SIZE = VECSIZE_INT * 2;
    int32_t segLen = (query_length + SIMD_SIZE-1) / SIMD_SIZE; /* number of segment */
    /* array to record the alignment read ending position of the largest score of each reference position */
    memset(this->maxColumn, 0, db_length * sizeof(uint16_t));
    uint16_t * maxColumn = (uint16_t *) this->maxColumn;

    /* Define 16 byte 0 vector. */
    simd_int vZero = simdi32_set(0);
    simd_int* pvHStore = vHStore;
    simd_int* pvHLoad = vHLoad;
    simd_int* pvE = vE;
    simd_int* pvHmax = vHmax;
    memset(pvHStore,0,segLen*sizeof(simd_int));
    memset(pvHLoad,0, segLen*sizeof(simd_int));
    memset(pvE,0,     segLen*sizeof(simd_int));
    memset(pvHmax,0,  segLen*sizeof(simd_int));

    int32_t i, j, k;
    /* 16 byte insertion begin vector */
#ifdef GAP_POS_SCORING
    simd_int vGapO;
    if (type != PROFILE_HMM) {
        vGapO = simdi16_set(gap_open);
    }
#else
    simd_int vGapO = simdi16_set(gap_open);
#endif

    /* 16 byte insertion extension vector */
    simd_int vGapE = simdi16_set(gap_extend);

    simd_int vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
    simd_int vMaxMark = vZero; /* Trace the highest score till the previous column. */
    simd_int vTemp;
    int32_t edge, begin = 0, end = db_length, step = 1;

    /* outer loop to process the reference sequence */
    if (ref_dir == 1) {
        begin = db_length - 1;
        end = -1;
        step = -1;
    }
    for (i = begin; LIKELY(i != end); i += step) {
        simd_int e, vF = vZero; /* Initialize F value to 0.
                                Any errors to vH values will be corrected in the Lazy_F loop.
                                */
        simd_int vH = pvHStore[segLen - 1];
        vH = simdi8_shiftl (vH, 2); /* Shift the 128-bit value in vH left by 2 byte. */

        /* Swap the 2 H buffers. */
        simd_int* pv = pvHLoad;

        simd_int vMaxColumn = vZero; /* vMaxColumn is used to record the max values of column i. */

        const simd_int* vPAA = query_aa_profile_word + db_aa_sequence[i] * segLen; /* Right part of the query_profile_byte */
        const simd_int* vP3Di = query_3di_profile_word + db_3di_sequence[i] * segLen; /* Right part of the query_profile_byte */

        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (j = 0; LIKELY(j < segLen); j ++) {
            simd_int score = simdi16_adds(simdi_load(vPAA + j), simdi_load(vP3Di + j));
            vH = simdi16_adds(vH, score);

            /* Get max from vH, vE and vF. */
            e = simdi_load(pvE + j);
            vH = simdi16_max(vH, e);
#ifdef GAP_POS_SCORING
            if (type == PROFILE_HMM) {
                vH = simdi16_max(vH, simdui16_subs(vF, simdi_load(gap_close_del + j)));
            } else {
#endif
                vH = simdi16_max(vH, vF);
#ifdef GAP_POS_SCORING
            }
#endif
            vMaxColumn = simdi16_max(vMaxColumn, vH);

            /* Save vH values. */
            simdi_store(pvHStore + j, vH);

            /* Update vE value. */
#ifdef GAP_POS_SCORING
            if (type == PROFILE_HMM) {
                // copy vH for update of vF
                vTemp = vH;
                vH = simdui16_subs(vH, simdi_load(gap_open_ins + j)); /* saturation arithmetic, result >= 0 */
            } else {
#endif
                vH = simdui16_subs(vH, vGapO); /* saturation arithmetic, result >= 0 */
#ifdef GAP_POS_SCORING
            }
#endif
            e = simdui16_subs(e, vGapE);
            e = simdi16_max(e, vH);
            simdi_store(pvE + j, e);

            /* Update vF value. */
            vF = simdui16_subs(vF, vGapE);
#ifdef GAP_POS_SCORING
            if (type == PROFILE_HMM) {
                vF = simdi16_max(vF, simdui16_subs(vTemp, simdi_load(gap_open_del + j)));
            } else {
#endif
                vF = simdi16_max(vF, vH);
#ifdef GAP_POS_SCORING
            }
#endif

            /* Load the next vH. */
            vH = simdi_load(pvHLoad + j);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
        for (k = 0; LIKELY(k < (int32_t) SIMD_SIZE); ++k) {
            vF = simdi8_shiftl (vF, 2);
            for (j = 0; LIKELY(j < segLen); ++j) {
                vH = simdi_load(pvHStore + j);
#ifdef GAP_POS_SCORING
                if (type == PROFILE_HMM) {
                    vH = simdi16_max(vH, simdui16_subs(vF, simdi_load(gap_close_del + j)));
                    simdi_store(pvE + j, simdi16_max(simdi_load(pvE + j), simdui16_subs(vH, simdi_load(gap_open_ins + j))));
                } else {
#endif
                    vH = simdi16_max(vH, vF);
#ifdef GAP_POS_SCORING
                }
#endif
                vMaxColumn = simdi16_max(vMaxColumn, vH); //newly added line
                simdi_store(pvHStore + j, vH);
#ifdef GAP_POS_SCORING
                if (type == PROFILE_HMM) {
                    vH = simdui16_subs(vH, simdi_load(gap_open_del + j));
                } else {
#endif
                    vH = simdui16_subs(vH, vGapO);
#ifdef GAP_POS_SCORING
                }
#endif
                vF = simdui16_subs(vF, vGapE);
                if (UNLIKELY(! simdi8_movemask(simdi16_gt(vF, vH)))) goto end;
            }
        }

        end:
        vMaxScore = simdi16_max(vMaxScore, vMaxColumn);
        vTemp = simdi16_eq(vMaxMark, vMaxScore);
        uint32_t cmp = simdi8_movemask(vTemp);
        if (cmp != SIMD_MOVEMASK_MAX) {
            uint16_t temp;
            vMaxMark = vMaxScore;
            max8(temp, vMaxScore);
            vMaxScore = vMaxMark;

            if (LIKELY(temp > max)) {
                max = temp;
                end_ref = i;
                for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
            }
        }

        /* Record the max score of current column. */
        max8(maxColumn[i], vMaxColumn);
        if (maxColumn[i] == terminate) break;
    }

    /* Trace the alignment ending position on read. */
    uint16_t *t = (uint16_t*)pvHmax;
    int32_t column_len = segLen * SIMD_SIZE;
    for (i = 0; LIKELY(i < column_len); ++i, ++t) {
        int32_t temp;
        if (*t == max) {
            temp = i / SIMD_SIZE + i % SIMD_SIZE * segLen;
            if (temp < end_read) end_read = temp;
        }
    }

    /* Find the most possible 2nd best alignment. */
    alignment_end best0;
    best0.score = max;
    best0.ref = end_ref;
    best0.read = end_read;

    alignment_end best1;
    best1.score = 0;
    best1.ref = 0;
    best1.read = 0;

    edge = (end_ref - maskLen) > 0 ? (end_ref - maskLen) : 0;
    for (i = 0; i < edge; i ++) {
        if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
        }
    }
    edge = (end_ref + maskLen) > db_length ? db_length : (end_ref + maskLen);
    for (i = edge; i < db_length; i ++) {
        if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
        }
    }

    return std::make_pair(best0, best1);
#undef max8
}

void StructureSmithWaterman::ssw_init(const Sequence* q_aa,
                                      const Sequence* q_3di,
                                      const int8_t* mat_aa,
                                      const int8_t* mat_3di,
                                      const BaseMatrix *m){
    profile->bias = 0;
    const int32_t alphabetSize = m->alphabetSize;
    int32_t compositionBias = 0;
    if (aaBiasCorrection) {
        SubstitutionMatrix::calcLocalAaBiasCorrection(m, q_aa->numSequence, q_aa->L, tmp_composition_bias, 1.0);
        for (int i =0; i < q_aa->L; i++) {
            profile->composition_bias_aa[i] = (int8_t) (tmp_composition_bias[i] < 0.0) ? tmp_composition_bias[i] - 0.5 : tmp_composition_bias[i] + 0.5;
            compositionBias = (compositionBias < profile->composition_bias_aa[i]) ? compositionBias : profile->composition_bias_aa[i];
        }
        SubstitutionMatrix::calcLocalAaBiasCorrection(m, q_3di->numSequence, q_3di->L, tmp_composition_bias, aaBiasCorrectionScale);
        for (int i =0; i < q_aa->L; i++) {
            profile->composition_bias_ss[i] = (int8_t) (tmp_composition_bias[i] < 0.0) ? tmp_composition_bias[i] - 0.5 : tmp_composition_bias[i] + 0.5;
            compositionBias = (compositionBias < profile->composition_bias_ss[i]) ? compositionBias : profile->composition_bias_ss[i];
        }
        compositionBias = std::min(compositionBias, 0);
    } else {
        memset(profile->composition_bias_aa, 0, q_aa->L * sizeof(int8_t));
        memset(profile->composition_bias_ss, 0, q_aa->L * sizeof(int8_t));
    }
    // copy memory to local memory todo: maybe change later
    memcpy(profile->mat_aa, mat_aa, alphabetSize * alphabetSize * sizeof(int8_t));
    memcpy(profile->mat_3di, mat_3di, alphabetSize * alphabetSize * sizeof(int8_t));
    memcpy(profile->query_aa_sequence, q_aa->numSequence, q_aa->L);
    memcpy(profile->query_3di_sequence, q_3di->numSequence, q_3di->L);

    bool isProfile = Parameters::isEqualDbtype(q_3di->getSequenceType(), Parameters::DBTYPE_HMM_PROFILE) &&
                     Parameters::isEqualDbtype(q_aa->getSequenceType(),  Parameters::DBTYPE_HMM_PROFILE);
    if(isProfile){
#ifdef GAP_POS_SCORING
        profile->gIns = q_3di->gIns;
        // insertion penalties are shifted by one position for the reverse direction (2nd to last becomes first)
        std::reverse_copy(q_3di->gIns, q_3di->gIns + q_3di->L - 1, profile->gIns_rev);
        for (int32_t i = 0; i < q_3di->L; ++i) {
            profile->gDelOpen[i] = q_3di->gDel[i] & 0xF;
            profile->gDelClose[i] = q_3di->gDel[i] >> 4;
        }
        profile->gDelClose_rev[0] = 0;
        profile->gDelOpen_rev[0] = 0;
        std::reverse_copy(profile->gDelOpen + 1, profile->gDelOpen + q_3di->L, profile->gDelClose_rev + 1);
        std::reverse_copy(profile->gDelClose + 1, profile->gDelClose + q_3di->L, profile->gDelOpen_rev + 1);
#endif
        memcpy(profile->alignment_aa_profile, q_aa->getAlignmentProfile(), q_aa->L * Sequence::PROFILE_AA_SIZE * sizeof(int8_t));
        memcpy(profile->alignment_3di_profile, q_3di->getAlignmentProfile(), q_3di->L * Sequence::PROFILE_AA_SIZE * sizeof(int8_t));
        // set neutral state 'X' (score=0)
        memset(profile->alignment_aa_profile + ((alphabetSize - 1) * q_aa->L), 0, q_aa->L * sizeof(int8_t));
        memset(profile->alignment_3di_profile + ((alphabetSize - 1) * q_3di->L), 0, q_3di->L * sizeof(int8_t));
    }

    /* Find the bias to use in the substitution matrix */
    int8_t bias = 0;


    if (isProfile) {
        int32_t mat3DiProfileSize = q_3di->L * Sequence::PROFILE_AA_SIZE;
        for (int32_t i = 0; i < mat3DiProfileSize; i++) {
            bias = std::min(profile->alignment_3di_profile[i], bias);
        }
        for (int32_t i = 0; i < mat3DiProfileSize; i++) {
            bias = std::min(profile->alignment_aa_profile[i], bias);
        }
    }else{
        int32_t mat3DiSize = q_3di->subMat->alphabetSize * q_3di->subMat->alphabetSize;
        int8_t bias3Di = 0;
        for (int32_t i = 0; i < mat3DiSize; i++) {
            bias3Di = std::min(bias3Di, mat_3di[i]);
        }
        int32_t matAASize =  q_aa->subMat->alphabetSize * q_aa->subMat->alphabetSize;
        int8_t biasAA = 0;
        for (int32_t i = 0; i < matAASize; i++){
            biasAA = std::min(biasAA, mat_aa[i]);
        }
        bias = bias3Di + biasAA;
    }

    bias = abs(bias) + abs(compositionBias);
    profile->bias = bias;
    profile->isProfile = false;
    if(isProfile){
        profile->isProfile = true;
        createQueryProfile<int8_t, VECSIZE_INT * 4, PROFILE>(profile->profile_aa_byte,profile->query_aa_sequence,NULL,
                                                             profile->alignment_aa_profile, q_aa->L, alphabetSize, profile->bias, 0, q_aa->L);
        createQueryProfile<int8_t, VECSIZE_INT * 4, PROFILE>(profile->profile_3di_byte, profile->query_3di_sequence, NULL,
                                                             profile->alignment_3di_profile, q_3di->L, alphabetSize, profile->bias, 0, q_3di->L);
        createQueryProfile<int16_t, VECSIZE_INT * 2, PROFILE>(profile->profile_aa_word,profile->query_aa_sequence,NULL,
                                                              profile->alignment_aa_profile,q_aa->L, alphabetSize, 0, 0, q_aa->L);
        createQueryProfile<int16_t, VECSIZE_INT * 2, PROFILE>(profile->profile_3di_word,profile->query_3di_sequence,NULL,
                                                              profile->alignment_3di_profile,q_3di->L, alphabetSize, 0, 0, q_3di->L);
#ifdef GAP_POS_SCORING
        createGapProfile<int8_t, VECSIZE_INT * 4>(profile->profile_gDelOpen_byte, profile->profile_gDelClose_byte,
                                                  profile->profile_gIns_byte, profile->gDelOpen, profile->gDelClose, q_3di->gIns, q_3di->L, 0);
        createGapProfile<int16_t, VECSIZE_INT * 2>(profile->profile_gDelOpen_word, profile->profile_gDelClose_word, profile->profile_gIns_word,
                                                   profile->gDelOpen, profile->gDelClose, q_3di->gIns, q_3di->L, 0);
#endif

    }else {
        createQueryProfile<int8_t, VECSIZE_INT * 4, SUBSTITUTIONMATRIX>(profile->profile_aa_byte,
                                                                        profile->query_aa_sequence,
                                                                        profile->composition_bias_aa, profile->mat_aa,
                                                                        q_aa->L, alphabetSize, bias, 0, 0);
        createQueryProfile<int8_t, VECSIZE_INT * 4, SUBSTITUTIONMATRIX>(profile->profile_3di_byte,
                                                                        profile->query_3di_sequence,
                                                                        profile->composition_bias_ss, profile->mat_3di,
                                                                        q_3di->L, alphabetSize, bias, 0, 0);
        createQueryProfile<int16_t, VECSIZE_INT * 2, SUBSTITUTIONMATRIX>(profile->profile_aa_word,
                                                                         profile->query_aa_sequence,
                                                                         profile->composition_bias_aa, profile->mat_aa,
                                                                         q_aa->L, alphabetSize, 0, 0, 0);
        createQueryProfile<int16_t, VECSIZE_INT * 2, SUBSTITUTIONMATRIX>(profile->profile_3di_word,
                                                                         profile->query_3di_sequence,
                                                                         profile->composition_bias_ss, profile->mat_3di,
                                                                         q_3di->L, alphabetSize, 0, 0, 0);
    }
    for(int32_t i = 0; i< alphabetSize; i++) {
        profile->profile_aa_word_linear[i] = &profile_aa_word_linear_data[i*q_aa->L];
        profile->profile_3di_word_linear[i] = &profile_3di_word_linear_data[i*q_3di->L];
        for (int j = 0; j < q_aa->L; j++) {
            profile->profile_aa_word_linear[i][j] = mat_aa[i * alphabetSize + q_aa->numSequence[j]] + profile->composition_bias_ss[j];
            profile->profile_3di_word_linear[i][j] = mat_3di[i * alphabetSize + q_3di->numSequence[j]];
        }
    }

    // create reverse structures
    std::reverse_copy(  profile->query_aa_sequence, profile->query_aa_sequence + q_aa->L, profile->query_aa_rev_sequence);
    std::reverse_copy(  profile->query_3di_sequence, profile->query_3di_sequence + q_3di->L, profile->query_3di_rev_sequence);
    std::reverse_copy(  profile->composition_bias_aa, profile->composition_bias_aa + q_aa->L, profile->composition_bias_aa_rev);
    std::reverse_copy(  profile->composition_bias_ss, profile->composition_bias_ss + q_3di->L, profile->composition_bias_ss_rev);

    profile->query_length = q_aa->L;
    profile->alphabetSize = alphabetSize;

    if (isProfile) {
        for (int32_t i = 0; i < alphabetSize; i++) {
            const int8_t *startToReadAA = &profile->alignment_aa_profile[i * q_aa->L];
            int8_t *startToWriteAA = &profile->rev_alignment_aa_profile[i * q_aa->L];
            std::reverse_copy(startToReadAA, startToReadAA + q_aa->L, startToWriteAA);
            const int8_t *startToRead3Di = &profile->alignment_3di_profile[i * q_3di->L];
            int8_t *startToWrite3Di = &profile->rev_alignment_3di_profile[i * q_3di->L];
            std::reverse_copy(startToRead3Di, startToRead3Di + q_3di->L, startToWrite3Di);
        }
    }
}
template <const unsigned int type>
StructureSmithWaterman::cigar * StructureSmithWaterman::banded_sw(const unsigned char *db_aa_sequence, const unsigned char *db_3di_sequence,
                                                                  const int8_t *query_aa_sequence, const int8_t *query_3di_sequence, const int8_t * compositionBiasAA,
                                                                  const int8_t * compositionBiasSS, int32_t db_length, int32_t query_length, int32_t queryStart, int32_t score, const uint32_t gap_open,
                                                                  const uint32_t gap_extend,
#ifdef GAP_POS_SCORING
                                                                  uint8_t *gDelOpen, uint8_t *gDelClose, uint8_t *gIns,
#endif
                                                                  int32_t band_width, const int8_t *mat_aa, int32_t nAA, const int8_t *mat_3di, int32_t n3Di) {
    /*! @function
     @abstract  Round an integer to the next closest power-2 integer.
     @param  x  integer to be rounded (in place)
     @discussion x will be modified.
     */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

    /* Convert the coordinate in the scoring matrix into the coordinate in one line of the band. */
#define set_u(u, w, i, j) { int x=(i)-(w); x=x>0?x:0; (u)=(j)-x+1; }

    /* Convert the coordinate in the direction matrix into the coordinate in one line of the band. */
#define set_d(u, w, i, j, p) { int x=(i)-(w); x=x>0?x:0; x=(j)-x; (u)=x*3+p; }

    uint32_t *c = (uint32_t*)malloc(16 * sizeof(uint32_t)), *c1;
    int32_t i, j, e, f, temp1, temp2, s = 16, s1 = 8, l, max = 0;
    int64_t s2 = 1024;
    char op, prev_op;
    int64_t width, width_d;
    int32_t *h_b, *e_b, *h_c;
    int8_t *direction, *direction_line;
    cigar* result = new cigar();
    h_b = (int32_t*)malloc(s1 * sizeof(int32_t));
    e_b = (int32_t*)malloc(s1 * sizeof(int32_t));
    h_c = (int32_t*)malloc(s1 * sizeof(int32_t));
    direction = (int8_t*)malloc(s2 * sizeof(int8_t));

    do {
        width = band_width * 2 + 3, width_d = band_width * 2 + 1;
        while (width >= s1) {
            ++s1;
            kroundup32(s1);
            h_b = (int32_t*)realloc(h_b, s1 * sizeof(int32_t));
            e_b = (int32_t*)realloc(e_b, s1 * sizeof(int32_t));
            h_c = (int32_t*)realloc(h_c, s1 * sizeof(int32_t));
        }
        int64_t targetSize = width_d * query_length * 3;
        while (targetSize >= s2) {
            ++s2;
            kroundup32(s2);
            if (s2 < 0) {
                fprintf(stderr, "Alignment score and position are not consensus.\n");
                EXIT(1);
            }
            direction = (int8_t*)realloc(direction, s2 * sizeof(int8_t));
        }
        direction_line = direction;
        for (j = 1; LIKELY(j < width - 1); j ++) h_b[j] = 0;
        for (i = 0; LIKELY(i < query_length); i ++) {
            int32_t beg = 0, end = db_length - 1, u = 0, edge;
            j = i - band_width;	beg = beg > j ? beg : j; // band start
            j = i + band_width; end = end < j ? end : j; // band end
            edge = end + 1 < width - 1 ? end + 1 : width - 1;
            f = h_b[0] = e_b[0] = h_b[edge] = e_b[edge] = h_c[0] = 0;
            int64_t directionOffset = width_d * i * 3;
            direction_line = direction + directionOffset;

            for (j = beg; LIKELY(j <= end); j ++) {
                int32_t b, e1, f1, d, de, df, dh;
                set_u(u, band_width, i, j);	set_u(e, band_width, i - 1, j);
                set_u(b, band_width, i, j - 1); set_u(d, band_width, i - 1, j - 1);
                set_d(de, band_width, i, j, 0);
                set_d(df, band_width, i, j, 1);
                set_d(dh, band_width, i, j, 2);
#ifdef GAP_POS_SCORING
                if (type == PROFILE_HMM) {
                    temp1 = i == 0 ? -gap_open : h_b[e] - gDelOpen[i];
                } else {
#endif
                    temp1 = i == 0 ? -gap_open : h_b[e] - gap_open;
#ifdef GAP_POS_SCORING
                }
#endif
                temp2 = i == 0 ? -gap_extend : e_b[e] - gap_extend;
                e_b[u] = temp1 > temp2 ? temp1 : temp2;
                direction_line[de] = temp1 > temp2 ? 3 : 2;

#ifdef GAP_POS_SCORING
                if (type == PROFILE_HMM) {
                    temp1 = h_c[b] - gIns[i];
                } else {
#endif
                    temp1 = h_c[b] - gap_open;
#ifdef GAP_POS_SCORING
                }
#endif
                temp2 = f - gap_extend;
                f = temp1 > temp2 ? temp1 : temp2;
                direction_line[df] = temp1 > temp2 ? 5 : 4;

#ifdef GAP_POS_SCORING
                if (type == PROFILE_HMM) {
                    e1 = std::max(0, e_b[u] - gDelClose[i + 1]);
                } else {
#endif
                    e1 = e_b[u] > 0 ? e_b[u] : 0;
#ifdef GAP_POS_SCORING
                }
#endif
                f1 = f > 0 ? f : 0;
                temp1 = e1 > f1 ? e1 : f1;
                if(type == SUBSTITUTIONMATRIX){
                    temp2 = h_b[d] + mat_aa[query_aa_sequence[i] * nAA + db_aa_sequence[j]] + compositionBiasAA[i]
                            + mat_3di[query_3di_sequence[i] * n3Di + db_3di_sequence[j]] + compositionBiasSS[i];
                }
                if(type == PROFILE || type == PROFILE_HMM) {
                    temp2 = h_b[d] + mat_aa[db_aa_sequence[j] * nAA + (queryStart + i)]
                            + mat_3di[db_3di_sequence[j] * n3Di + (queryStart + i)];
                }

                h_c[u] = temp1 > temp2 ? temp1 : temp2;

                if (h_c[u] > max) max = h_c[u];

                if (temp1 <= temp2) direction_line[dh] = 1;
                else direction_line[dh] = e1 > f1 ? direction_line[de] : direction_line[df];
            }
            for (j = 1; j <= u; j ++) h_b[j] = h_c[j];
        }
        band_width *= 2;
    } while (LIKELY(max < score));
    band_width /= 2;

    // trace back
    i = query_length - 1;
    j = db_length - 1;
    e = 0;	// Count the number of M, D or I.
    l = 0;	// record length of current cigar
    op = prev_op = 'M';
    temp2 = 2;	// h
    while (LIKELY(i > 0) || LIKELY(j > 0)) {
        set_d(temp1, band_width, i, j, temp2);
        switch (direction_line[temp1]) {
            case 1:
                --i;
                --j;
                temp2 = 2;
                direction_line -= width_d * 3;
                op = 'M';
                break;
            case 2:
                --i;
                temp2 = 0;	// e
                direction_line -= width_d * 3;
                op = 'I';
                break;
            case 3:
                --i;
                temp2 = 2;
                direction_line -= width_d * 3;
                op = 'I';
                break;
            case 4:
                --j;
                temp2 = 1;
                op = 'D';
                break;
            case 5:
                --j;
                temp2 = 2;
                op = 'D';
                break;
            default:
                fprintf(stderr, "Trace back error: %d.\n", direction_line[temp1 - 1]);
                free(direction);
                free(h_c);
                free(e_b);
                free(h_b);
                free(c);
                delete result;
                return 0;
        }
        if (op == prev_op) ++e;
        else {
            ++l;
            while (l >= s) {
                ++s;
                kroundup32(s);
                c = (uint32_t*)realloc(c, s * sizeof(uint32_t));
            }
            c[l - 1] = to_cigar_int(e, prev_op);
            prev_op = op;
            e = 1;
        }
    }
    if (op == 'M') {
        ++l;
        while (l >= s) {
            ++s;
            kroundup32(s);
            c = (uint32_t*)realloc(c, s * sizeof(uint32_t));
        }
        c[l - 1] = to_cigar_int(e + 1, op);
    }else {
        l += 2;
        while (l >= s) {
            ++s;
            kroundup32(s);
            c = (uint32_t*)realloc(c, s * sizeof(uint32_t));
        }
        c[l - 2] = to_cigar_int(e, op);
        c[l - 1] = to_cigar_int(1, 'M');
    }

    // reverse cigar
    c1 = (uint32_t*)new uint32_t[l * sizeof(uint32_t)];
    s = 0;
    e = l - 1;
    while (LIKELY(s <= e)) {
        c1[s] = c[e];
        c1[e] = c[s];
        ++ s;
        -- e;
    }
    result->seq = c1;
    result->length = l;

    free(direction);
    free(h_c);
    free(e_b);
    free(h_b);
    free(c);
    return result;
#undef kroundup32
#undef set_u
#undef set_d
}

uint32_t StructureSmithWaterman::to_cigar_int (uint32_t length, char op_letter)
{
    uint32_t res;
    uint8_t op_code;

    switch (op_letter) {
        case 'M': /* alignment match (can be a sequence match or mismatch */
        default:
            op_code = 0;
            break;
        case 'I': /* insertion to the reference */
            op_code = 1;
            break;
        case 'D': /* deletion from the reference */
            op_code = 2;
            break;
        case 'N': /* skipped region from the reference */
            op_code = 3;
            break;
        case 'S': /* soft clipping (clipped sequences present in SEQ) */
            op_code = 4;
            break;
        case 'H': /* hard clipping (clipped sequences NOT present in SEQ) */
            op_code = 5;
            break;
        case 'P': /* padding (silent deletion from padded reference) */
            op_code = 6;
            break;
        case '=': /* sequence match */
            op_code = 7;
            break;
        case 'X': /* sequence mismatch */
            op_code = 8;
            break;
    }

    res = (length << 4) | op_code;
    return res;
}

void StructureSmithWaterman::printVector(__m128i v){
    for (int i = 0; i < 8; i++)
        printf("%d ", ((short) (sse2_extract_epi16(v, i)) + 32768));
    std::cout << "\n";
}

void StructureSmithWaterman::printVectorUS(__m128i v){
    for (int i = 0; i < 8; i++)
        printf("%d ", (unsigned short) sse2_extract_epi16(v, i));
    std::cout << "\n";
}

unsigned short StructureSmithWaterman::sse2_extract_epi16(__m128i v, int pos) {
    switch(pos){
        case 0: return _mm_extract_epi16(v, 0);
        case 1: return _mm_extract_epi16(v, 1);
        case 2: return _mm_extract_epi16(v, 2);
        case 3: return _mm_extract_epi16(v, 3);
        case 4: return _mm_extract_epi16(v, 4);
        case 5: return _mm_extract_epi16(v, 5);
        case 6: return _mm_extract_epi16(v, 6);
        case 7: return _mm_extract_epi16(v, 7);
    }
    std::cerr << "Fatal error in QueryScore: position in the vector is not in the legal range (pos = " << pos << ")\n";
    EXIT(1);
    // never executed
    return 0;
}

float StructureSmithWaterman::computeCov(unsigned int startPos, unsigned int endPos, unsigned int len) {
    return (std::min(len, std::max(startPos, endPos)) - std::min(startPos, endPos) + 1) / (float) len;
}

template <typename F>
inline F simd_hmax(const F * in, unsigned int n) {
    F current = std::numeric_limits<F>::min();
    do {
        current = std::max(current, *in++);
    } while(--n);

    return current;
}

//template<const int T>
//simd_int StructureSmithWaterman::needlemanWunschScoreVec(const simd_int * subQNNi, const simd_int * target_sub_vec,
//                                 const  simd_int *subQ_dist, const  simd_int *subT_dist,
//                                 const simd_int nwGapPenalty, const simd_int vSubMatBias, const simd_int vDistMatBias){
//
//    simd_int prev_sMat_vec[T+1] __attribute__((aligned(ALIGN_INT)));
//    simd_int sMat_vec[T+1] __attribute__((aligned(ALIGN_INT)));
//    memset(prev_sMat_vec, 0, sizeof(simd_int) * (T+1));
//
//    for(int i = 1; i < T+1; i++){
//        sMat_vec[0]= simdi_setzero();
//        for(int j = 1; j < T+1; j++){
//            // score
//
//            simd_int scoreLookup = UngappedAlignment::Shuffle(target_sub_vec[i-1], subQNNi[j-1]);
//            scoreLookup = simdi_and(scoreLookup, simdi16_set(0x00FF));
//            scoreLookup = simdi16_sub(scoreLookup, vSubMatBias);
//
//            simd_int distLookup = UngappedAlignment::Shuffle(subT_dist[i-1], subQ_dist[j-1]);
//            distLookup = simdi_and(distLookup, simdi16_set(0x00FF));
//            distLookup = simdi16_sub(distLookup, vDistMatBias);
//
//            // add
//            scoreLookup = simdi16_add(scoreLookup, distLookup);
//            sMat_vec[j] = simdi16_max(simdi16_add(prev_sMat_vec[j-1],scoreLookup),
//                                      simdi16_max(simdi16_add(prev_sMat_vec[j],nwGapPenalty),
//                                                  simdi16_add(sMat_vec[j-1],nwGapPenalty)));
//        }
//        std::swap(prev_sMat_vec, sMat_vec);
//    }
//
//    return prev_sMat_vec[T]; /// 4;
//}
