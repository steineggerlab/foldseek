//
// Created by Martin Steinegger on 07/07/2022.
//

#include "TMaligner.h"
#include "tmalign/Coordinates.h"
#include <tmalign/TMalign.h>
#include "StructureSmithWaterman.h"
#include "StructureSmithWaterman.h"

TMaligner::TMaligner(unsigned int maxSeqLen, bool tmAlignFast, bool tmScoreOnly)
   : tmAlignFast(tmAlignFast),
     xtm(maxSeqLen), ytm(maxSeqLen), xt(maxSeqLen),
     r1(maxSeqLen), r2(maxSeqLen){
    affineNW = NULL;
    if(tmScoreOnly == false){
        affineNW = new AffineNeedlemanWunsch(maxSeqLen, 20);
    }
    query_x = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float) );
    query_y = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float) );
    query_z = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float) );
    target_x = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float) );
    target_y = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float) );
    target_z = (float*)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float) );
    mem = (float*)mem_align(ALIGN_FLOAT,6*maxSeqLen*4*sizeof(float));
    querySecStruc  = new char[maxSeqLen];
    targetSecStruc = new char[maxSeqLen];
    invmap = new int[maxSeqLen];
}

TMaligner::~TMaligner(){
    if(affineNW != NULL){
        delete affineNW;
    }
    free(query_x);
    free(query_y);
    free(query_z);
    free(target_x);
    free(target_y);
    free(target_z);
    free(mem);
    delete [] querySecStruc;
    delete [] targetSecStruc;
    delete [] invmap;
}

TMaligner::TMscoreResult TMaligner::computeTMscore(float *x, float *y, float *z, unsigned int targetLen,
                                                   int qStartPos, int dbStartPos, const std::string &backtrace,
                                                   int normalizationLen) {
    int qPos = qStartPos;
    int tPos = dbStartPos;
    std::string cigarString = backtrace;
    std::fill(invmap, invmap+queryLen, -1);
    for (size_t btPos = 0; btPos < cigarString.size(); btPos++) {
        if (cigarString[btPos] == 'M') {
            invmap[qPos] = tPos;
            qPos++;
            tPos++;
        }
        else if (cigarString[btPos] == 'I') {
            qPos++;
        }
        else {
            tPos++;
        }
    }

    memcpy(target_x, x, sizeof(float) * targetLen);
    memcpy(target_y, y, sizeof(float) * targetLen);
    memcpy(target_z, z, sizeof(float) * targetLen);
    Coordinates targetCaCords;
    targetCaCords.x = target_x;
    targetCaCords.y = target_y;
    targetCaCords.z = target_z;
    Coordinates queryCaCords;
    queryCaCords.x = query_x;
    queryCaCords.y = query_y;
    queryCaCords.z = query_z;
    float t[3], u[3][3];
    float D0_MIN;

    float rmsd0 = 0.0;
    int L_ali;                // Aligned length in standard_TMscore
    float Lnorm;         //normalization length
    float score_d8,d0,d0_search,dcu0;//for TMscore search
    parameter_set4search(normalizationLen,  normalizationLen, D0_MIN, Lnorm,
                         score_d8, d0, d0_search, dcu0);
    double prevD0_MIN = D0_MIN;// stored for later use
    int prevLnorm = Lnorm;
    double prevd0 = d0;
    double local_d0_search = d0_search;
    double TMalnScore = standard_TMscore(r1, r2, xtm, ytm, xt, targetCaCords, queryCaCords, queryLen, invmap,
                                         L_ali, rmsd0, D0_MIN, Lnorm, d0, score_d8, t, u,  mem);
    D0_MIN = prevD0_MIN;
    Lnorm = prevLnorm;
    d0 = prevd0;
    double TM = detailed_search_standard(r1, r2, xtm, ytm, xt, targetCaCords, queryCaCords, queryLen,
                                         invmap, t, u, 40, local_d0_search, true, Lnorm, score_d8, d0, mem);
    TM = std::max(TM, TMalnScore);
    return TMaligner::TMscoreResult(u, t, TM, rmsd0);
}

void TMaligner::initQuery(float *x, float *y, float *z, char * querySeq, unsigned int queryLen){
    memset(querySecStruc, 0, sizeof(char) * queryLen);
    memcpy(query_x, x, sizeof(float) * queryLen);
    memcpy(query_y, y, sizeof(float) * queryLen);
    memcpy(query_z, z, sizeof(float) * queryLen);
    this->queryLen = queryLen;
    this->querySeq = querySeq;
    Coordinates queryCaCords;
    queryCaCords.x = query_x;
    queryCaCords.y = query_y;
    queryCaCords.z = query_z;
    make_sec(queryCaCords, queryLen, querySecStruc); // secondary structure assignment

}

Matcher::result_t TMaligner::align(unsigned int dbKey, float *x, float *y, float *z, char * targetSeq, unsigned int targetLen, float &TM1){
    backtrace.clear();

    memcpy(target_x, x, sizeof(float) * targetLen);
    memcpy(target_y, y, sizeof(float) * targetLen);
    memcpy(target_z, z, sizeof(float) * targetLen);
    Coordinates targetCaCords;
    targetCaCords.x = target_x;
    targetCaCords.y = target_y;
    targetCaCords.z = target_z;
    Coordinates queryCaCords;
    queryCaCords.x = query_x;
    queryCaCords.y = query_y;
    queryCaCords.z = query_z;

    const bool I_opt = false; // flag for -I, stick to user given alignment
    const bool a_opt = false; // flag for -a, normalized by average length
    const bool u_opt = false; // flag for -u, normalized by user specified length
    const bool d_opt = false; // flag for -d, user specified d0
    double Lnorm_ass = 0.0;
    double  d0_scale = 0.0;

    memset(targetSecStruc, 0, sizeof(char) * targetLen);
    make_sec(targetCaCords, targetLen, targetSecStruc); // secondary structure assignment
    /* entry function for structure alignment */
    float t0[3], u0[3][3];
    float TM2;
    float TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
    float d0_0, TM_0;
    float d0A, d0B, d0u, d0a;
    float d0_out = 5.0;
    float rmsd0 = 0.0;
    float Liden = 0;
    int n_ali = 0;
    int n_ali8 = 0;
    if(queryLen <=5 || targetLen <=5){
        return Matcher::result_t();
    }
    TMalign_main(affineNW,
                 targetCaCords, queryCaCords, targetSeq, querySeq, targetSecStruc, querySecStruc,
                 t0, u0, TM1, TM2, TM3, TM4, TM5,
                 d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                 seqM, seqxA, seqyA,
                 rmsd0, Liden,  n_ali, n_ali8,
                 targetLen, queryLen, Lnorm_ass, d0_scale,
                 I_opt, a_opt, u_opt, d_opt, tmAlignFast, mem, xtm, ytm, xt, r1, r2);
    //std::cout << queryId << "\t" << targetId << "\t" <<  TM_0 << "\t" << TM1 << std::endl;

    //double seqId = (n_ali8 > 0) ? (Liden / (static_cast<double>(n_ali8))) : 0;
    //std::cout << Liden << "\t" << n_ali8 << "\t" << queryLen << std::endl;
    //int rmsdScore = static_cast<int>(rmsd0*1000.0);

    // compute freeshift-backtrace from pairwise global alignment
    // find last match position first
    size_t shiftQ = 0, shiftT = 0, lastM = 0;
    for (size_t i = 0; i < seqxA.length(); ++i) {
        if (seqxA[i] != '-' && seqyA[i] != '-') {
            lastM = i;
        }
    }
    // compute backtrace from first match to last match
    bool hasMatch = false;
    int aaIdCnt = 0;
    for (size_t i = 0; i <= lastM; ++i) {
        if (seqxA[i] != '-' && seqyA[i] != '-') {
            backtrace.append(1, 'M');
            aaIdCnt += (seqxA[i] == seqyA[i]);
            hasMatch = true;
        } else if (seqxA[i] == '-' && seqyA[i] != '-') {
            if (hasMatch == false) {
                shiftQ++;
            } else {
                backtrace.append(1, 'I');
            }
        } else if (seqxA[i] != '-' && seqyA[i] == '-') {
            if (hasMatch == false) {
                shiftT++;
            } else {
                backtrace.append(1, 'D');
            }
        }
    }
    // compute end offsets to fix end positions
    size_t endQ = 0, endT = 0;
    for (size_t i = lastM; i < seqxA.length(); ++i) {
        if (seqxA[i] == '-' && seqyA[i] != '-') {
            endQ++;
        } else if (seqxA[i] != '-' && seqyA[i] == '-') {
            endT++;
        }
    }
    unsigned int alnLength = backtrace.size();
    float seqId = static_cast<float>(aaIdCnt)/static_cast<float>(alnLength);

    float qCov = StructureSmithWaterman::computeCov(shiftQ, queryLen-endQ-1, queryLen);
    float tCov = StructureSmithWaterman::computeCov(shiftT, targetLen-endT-1, targetLen);
    return Matcher::result_t(dbKey, TM_0*100000 , qCov, tCov, seqId, TM2, backtrace.length(), shiftQ, queryLen-endQ-1, queryLen, shiftT, targetLen-endT-1, targetLen, Matcher::compressAlignment(backtrace));
}