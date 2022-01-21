//
// Created by Charlotte Tumescheit on 2021/12/23.
//

#ifndef FOLDSEEK_STRUCTURESMITHWATERMAN_H
#define FOLDSEEK_STRUCTURESMITHWATERMAN_H

/* $Id: smith_waterman_sse2.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

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
  Written by Michael Farrar, 2006.
  Please send bug reports and/or suggestions to farrar.michael@gmail.com.
*/

#include <climits>
#include "structureto3diseqdist.h"
#include <cstdio>
#include <cstdlib>

#if !defined(__APPLE__) && !defined(__llvm__)
#include <malloc.h>
#endif

#include "simd.h"
#include "BaseMatrix.h"

#include "Sequence.h"
#include "EvalueComputation.h"
#include "EvalueNeuralNet.h"


class StructureSmithWaterman{
public:

    StructureSmithWaterman(size_t maxSequenceLength, int aaSize, bool aaBiasCorrection, float aaBiasCorrectionScale);
    ~StructureSmithWaterman();

    // prints a __m128 vector containing 8 signed shorts
    static void printVector (__m128i v);

    // prints a __m128 vector containing 8 unsigned shorts, added 32768
    static void printVectorUS (__m128i v);

    static unsigned short sse2_extract_epi16(__m128i v, int pos);

    // The dynamic programming matrix entries for the query and database sequences are stored sequentially (the order see the Farrar paper).
    // This function calculates the index within the dynamic programming matrices for the given query and database sequence position.
    static inline int midx (int qpos, int dbpos, int iter){
        return dbpos * (8 * iter) + (qpos % iter) * 8 + (qpos / iter);
    }

    typedef struct {
        short qStartPos;
        short dbStartPos;
        short qEndPos;
        short dbEndPos;
    } aln_t;


    typedef struct {
        uint32_t score1;
        uint32_t score2;
        int32_t dbStartPos1;
        int32_t dbEndPos1;
        int32_t	qStartPos1;
        int32_t qEndPos1;
        int32_t ref_end2;
        float qCov;
        float tCov;
        uint32_t* cigar;
        double evalue;
        int identicalAACnt;
        int32_t cigarLen;
        int word;
    } s_align;

    // @function	ssw alignment.
    /*!	@function	Do Striped Smith-Waterman alignment.

     @param	maskLen	The distance between the optimal and suboptimal alignment ending position >= maskLen. We suggest to use
     readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function will NOT
     return the suboptimal alignment information. Detailed description of maskLen: After locating the optimal
     alignment ending position, the suboptimal alignment score can be heuristically found by checking the second
     largest score in the array that contains the maximal score of each column of the SW matrix. In order to avoid
     picking the scores that belong to the alignments sharing the partial best alignment, SSW C library masks the
     reference loci nearby (mask length = maskLen) the best alignment ending position and locates the second largest
     score from the unmasked elements.

     @return	pointer to the alignment result structure

     @note	Whatever the parameter flag is setted, this function will at least return the optimal and sub-optimal alignment score,
     and the optimal alignment ending positions on target and query sequences. If both bit 6 and 7 of the flag are setted
     while bit 8 is not, the function will return cigar only when both criteria are fulfilled. All returned positions are
     0-based coordinate.
     */
    s_align alignScoreEndPos (
            const unsigned char *db_aa_sequence,
            const unsigned char *db_3di_sequence,
            int32_t db_length,
            const uint8_t gap_open,
            const uint8_t gap_extend,
            const int32_t maskLen);

    s_align alignStartPosBacktrace (
            const unsigned char *db_aa_sequence,
            const unsigned char *db_3di_sequence,
            int32_t db_length,
            const uint8_t gap_open,
            const uint8_t gap_extend,
            const uint8_t alignmentMode,	//  (from high to low) bit 5: return the best alignment beginning position; 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar; 7: if max score >= filters, return cigar; 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
            std::string & backtrace,
            StructureSmithWaterman::s_align r,
            const int covMode, const float covThr,
            const int32_t maskLen);


    /*!	@function	Create the query profile using the query sequence.
     @param	read	pointer to the query sequence; the query sequence needs to be numbers
     @param	readLen	length of the query sequence
     @param	mat	pointer to the substitution matrix; mat needs to be corresponding to the read sequence
     @param	n	the square root of the number of elements in mat (mat has n*n elements)
     @param	score_size	estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 please set 0; if
     your estimated best alignment score >= 255, please set 1; if you don't know, please set 2
     @return	pointer to the query profile structure
     @note	example for parameter read and mat:
     If the query sequence is: ACGTATC, the sequence that read points to can be: 1234142
     Then if the penalty for match is 2 and for mismatch is -2, the substitution matrix of parameter mat will be:
     //A  C  G  T
     2 -2 -2 -2 //A
     -2  2 -2 -2 //C
     -2 -2  2 -2 //G
     -2 -2 -2  2 //T
     mat is the pointer to the array {2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2}
     */
    void ssw_init(const Sequence *q_aa, const Sequence *q_3di, const int8_t *mat_aa, const int8_t *mat_3di, const BaseMatrix *m);


    static char cigar_int_to_op (uint32_t cigar_int);

    static uint32_t cigar_int_to_len (uint32_t cigar_int);


    static float computeCov(unsigned int startPos, unsigned int endPos, unsigned int len);

    static void seq_reverse(int8_t * reverse, const int8_t* seq, int32_t end)	/* end is 0-based alignment ending position */
    {
        int32_t start = 0;
        while (LIKELY(start <= end)) {
            reverse[start] = seq[end];
            reverse[end] = seq[start];
            ++start;
            --end;
        }
    }


private:

    struct s_profile{
        simd_int* profile_aa_byte;	// 0: none
        simd_int* profile_aa_word;	// 0: none
        simd_int* profile_aa_rev_byte;	// 0: none
        simd_int* profile_aa_rev_word;	// 0: none
        int8_t* query_aa_sequence;
        int8_t* query_aa_rev_sequence;
        simd_int* profile_3di_byte;	// 0: none
        simd_int* profile_3di_word;	// 0: none
        simd_int* profile_3di_rev_byte;	// 0: none
        simd_int* profile_3di_rev_word;	// 0: none
        int8_t* query_3di_sequence;
        int8_t* query_3di_rev_sequence;
        int8_t* composition_bias_ss;
        int8_t* composition_bias_aa;
        int8_t* composition_bias_aa_rev;
        int8_t* composition_bias_ss_rev;
        int8_t* mat_aa;
        int8_t* mat_3di;
        // Memory layout of if mat + queryProfile is qL * AA
        //    Query length
        // A  -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
        // C  -1  -4   2   5  -3  -2   0  -3   1  -3  -2   0  -1   2   0   0  -1  -3  -4  -2
        // ...
        // Y -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
        // Memory layout of if mat + sub is AA * AA
        //     A   C    ...                                                                Y
        // A  -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
        // C  -1  -4   2   5  -3  -2   0  -3   1  -3  -2   0  -1   2   0   0  -1  -3  -4  -2
        // ...
        // Y -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
        int8_t* mat_rev; // needed for queryProfile
        int32_t query_length;
        int32_t alphabetSize;
        uint8_t bias;
        short ** profile_aa_word_linear;
        short ** profile_3di_word_linear;
    };
    simd_int* vHStore;
    simd_int* vHLoad;
    simd_int* vE;
    simd_int* vHmax;
    uint8_t * maxColumn;

    typedef struct {
        uint16_t score;
        int32_t ref;	 //0-based position
        int32_t read;    //alignment ending position on read, 0-based
    } alignment_end;


    typedef struct {
        uint32_t* seq;
        int32_t length;
    } cigar;

    /* Striped Smith-Waterman
     Record the highest score of each reference position.
     Return the alignment score and ending position of the best alignment, 2nd best alignment, etc.
     Gap begin and gap extension are different.
     wight_match > 0, all other weights < 0.
     The returned positions are 0-based.
     */
    std::pair<alignment_end, alignment_end> sw_sse2_byte (const unsigned char*db_aa_sequence,
                                                          const unsigned char*db_3di_sequence,
                                                          int8_t ref_dir,	// 0: forward ref; 1: reverse ref
                                                          int32_t db_length,
                                                          int32_t query_length,
                                                          const uint8_t gap_open, /* will be used as - */
                                                          const uint8_t gap_extend, /* will be used as - */
                                                          const simd_int* query_aa_profile_byte,
                                                          const simd_int* query_3di_profile_byte,
                                                          uint8_t terminate,	/* the best alignment score: used to terminate
                                                     the matrix calculation when locating the
                                                     alignment beginning point. If this score
                                                     is set to 0, it will not be used */
                                                          uint8_t bias,  /* Shift 0 point to a positive value. */
                                                          int32_t maskLen);

    std::pair<alignment_end, alignment_end> sw_sse2_word (const unsigned char* db_aa_sequence,
                                                          const unsigned char* db_3di_sequence,
                                                          int8_t ref_dir,	// 0: forward ref; 1: reverse ref
                                                          int32_t db_length,
                                                          int32_t query_length,
                                                          const uint8_t gap_open, /* will be used as - */
                                                          const uint8_t gap_extend, /* will be used as - */
                                                          const simd_int*query_aa_profile_byte,
                                                          const simd_int*query_3di_profile_byte,
                                                          uint16_t terminate,
                                                          int32_t maskLen);

    StructureSmithWaterman::cigar *banded_sw(const unsigned char *db_aa_sequence, const unsigned char *db_3di_sequence,
                                             const int8_t *query_aa_sequence, const int8_t *query_3di_sequence,
                                             const int8_t * compositionBiasAA, const int8_t * compositionBiasSS,
                                             int32_t db_length, int32_t query_length,
                                             int32_t score, const uint32_t gap_open,
                                             const uint32_t gap_extend, int32_t band_width,
                                             const int8_t *mat_aa, const int8_t *mat_3di, int32_t n);

    void computerBacktrace(s_profile * query, const unsigned char * db_sequence,
                           s_align & alignment, std::string & backtrace, uint32_t & aaIds, size_t & mStatesCnt);


//    template<const int T>
//    simd_int needlemanWunschScoreVec(const simd_int * subQNNi, const simd_int * target_sub_vec,
//                                            const  simd_int *subQ_dist, const  simd_int *subT_dist,
//                                            const simd_int nwGapPenalty, const simd_int vSubMatBias, const simd_int vDistMatBias);

    /*!	@function		Produce CIGAR 32-bit unsigned integer from CIGAR operation and CIGAR length
     @param	length		length of CIGAR
     @param	op_letter	CIGAR operation character ('M', 'I', etc)
     @return			32-bit unsigned integer, representing encoded CIGAR operation and length
     */
    inline uint32_t to_cigar_int (uint32_t length, char op_letter);

    s_profile* profile;

    const static unsigned int SUBSTITUTIONMATRIX = 1;
    const static unsigned int PROFILE = 2;

    template <typename T, size_t Elements, const unsigned int type>
    void createQueryProfile(simd_int *profile, const int8_t *query_sequence, const int8_t * composition_bias, const int8_t *mat, const int32_t query_length, const int32_t aaSize, uint8_t bias, const int32_t offset, const int32_t entryLength);

    float *tmp_composition_bias;
    short * profile_aa_word_linear_data;
    short * profile_3di_word_linear_data;
    bool aaBiasCorrection;
    float aaBiasCorrectionScale;
};


const int DIST_MAT_SIZE = Alphabet3diSeqDist::CENTROID_CNT + 1;
//    A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   X
const short distMat[DIST_MAT_SIZE][DIST_MAT_SIZE] = { {14,5, -19, -50, -50, -50, -50,-50, -50 ,-50, -50, -50, -50, -50, -50, -50, -50, -50, -50,-50, 0},
                                                      {5, 11,  2, -11, -15, -18, -23, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, 0},
                                                      {-19,   2,   9,   2,  -7, -13, -16, -23, -50, -24, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                                                      {-50, -11,   2,   7,   2,  -6, -10, -14, -18, -22, -27, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                                                      {-50, -15,  -7,   2,   7,   2,  -5,  -9, -12, -17, -22, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                                                      {-50, -18, -13,  -6,   2,   7,   2,  -4,  -9, -14, -19, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                                                      {-50, -23, -16, -10,  -5,   2,   8,   4,  -3,  -7, -13, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                                                      {-50, -50, -23, -14,  -9,  -4,   4,   9,   1,  -3, -10, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                                                      {-50, -50, -50, -18, -12,  -9,  -3,   1,  11,   2,  -6, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                                                      {-50, -50, -24, -22, -17, -14,  -7,  -3,   2,  10,  -5, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                                                      {-50, -50, -50, -27, -22, -19, -13, -10,  -6,  -5,   6, -50, -50, -50, -50, -50, -50, -50, -50, -50,   0},
                                                      {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50,   8,  -4,  -8, -13, -18, -26, -50, -50, -50,   0},
                                                      {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50,  -4,   9,   0,  -6, -12, -16, -22, -50, -50,   0},
                                                      {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50,  -8,   0,   9,   3,  -5, -12, -17, -23, -50,   0},
                                                      {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -13,  -6,   3,   7,   1,  -8, -12, -19, -50,   0},
                                                      {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -18, -12,  -5,   1,   6,   1,  -6, -15, -50,   0},
                                                      {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -26, -16, -12,  -8,   1,   6,   2, -10, -50,   0},
                                                      {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -22, -17, -12,  -6,   2,   8,   2, -50,   0},
                                                      {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -23, -19, -15, -10,   2,  10,   4,   0},
                                                      {-50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50,   4,  12,   0},
                                                      { 0,   0,   0,   0,   0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,   0,   0,   0}};

#endif //FOLDSEEK_STRUCTURESMITHWATERMAN_H
