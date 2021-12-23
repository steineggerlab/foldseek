
#include <string>
#include <vector>
#include <sstream>
#include <tmalign/TMalign.h>
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "QueryMatcher.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "PareunAlign.h"
#include "simd.h"
#include "structureto3diseqdist.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include <cmath>


//// NN simd version with ungapped alignment to make it faster but sensitiviy suffers :(
//template<const int T>
//simd_int needlemanWunschScoreVec(const simd_int * subQNNi, const simd_int * target_sub_vec,
//                                 const  simd_int *subQ_dist, const  simd_int *subT_dist,
//                                 const simd_int nwGapPenalty, const simd_int vSubMatBias, const simd_int vDistMatBias){
//
//    simd_int prev_sMat_vec[T+1] __attribute__((aligned(ALIGN_INT)));
//    simd_int sMat_vec[T+1] __attribute__((aligned(ALIGN_INT)));
//    memset(prev_sMat_vec, 0, sizeof(simd_int) * (T+1));
//    simd_int value = simdi_setzero();
//    for(int i = 1; i < T+1; i++){
//        sMat_vec[0]= simdi_setzero();
//
//        // score
//        simd_int scoreLookup = UngappedAlignment::Shuffle(target_sub_vec[i-1], subQNNi[i-1]);
//        scoreLookup = simdi_and(scoreLookup, simdi16_set(0x00FF));
//// do this after the inner loop and use simdi8 TODO
//        scoreLookup = simdi16_sub(scoreLookup, vSubMatBias);
////            for(int a = 0; a < 9; a++){
////                cout << ((short *)&scoreLookup)[a] << " ";
////            }
////            cout << endl;
//        // distance
//        simd_int distLookup = UngappedAlignment::Shuffle(subT_dist[i-1], subQ_dist[i-1]);
//        distLookup = simdi_and(distLookup, simdi16_set(0x00FF));
//        distLookup = simdi16_sub(distLookup, vDistMatBias);
//
//        // add
//        scoreLookup = simdi16_add(scoreLookup, distLookup);
//        value = simdi16_add(value,scoreLookup);
//    }
//
//    return value; /// 4;
//}

template<const int T>
simd_int needlemanWunschScoreVec(const simd_int * subQNNi, const simd_int * target_sub_vec,
                                 const  simd_int *subQ_dist, const  simd_int *subT_dist,
                                 const simd_int nwGapPenalty, const simd_int vSubMatBias, const simd_int vDistMatBias){

    simd_int prev_sMat_vec[T+1] __attribute__((aligned(ALIGN_INT)));
    simd_int sMat_vec[T+1] __attribute__((aligned(ALIGN_INT)));
    memset(prev_sMat_vec, 0, sizeof(simd_int) * (T+1));

    for(int i = 1; i < T+1; i++){
        sMat_vec[0]= simdi_setzero();
        for(int j = 1; j < T+1; j++){
            // score
            simd_int scoreLookup = UngappedAlignment::Shuffle(target_sub_vec[i-1], subQNNi[j-1]);
            scoreLookup = simdi_and(scoreLookup, simdi16_set(0x00FF));
            scoreLookup = simdi16_sub(scoreLookup, vSubMatBias);
//            for(int a = 0; a < 9; a++){
//                cout << ((short *)&scoreLookup)[a] << " ";
//            }
//            cout << endl;
            // distance
            simd_int distLookup = UngappedAlignment::Shuffle(subT_dist[i-1], subQ_dist[j-1]);
            distLookup = simdi_and(distLookup, simdi16_set(0x00FF));
            distLookup = simdi16_sub(distLookup, vDistMatBias);

            // add
            scoreLookup = simdi16_add(scoreLookup, distLookup);
            sMat_vec[j] = simdi16_max(simdi16_add(prev_sMat_vec[j-1],scoreLookup),
                                      simdi16_max(simdi16_add(prev_sMat_vec[j],nwGapPenalty),
                                                  simdi16_add(sMat_vec[j-1],nwGapPenalty)));
        }
        std::swap(prev_sMat_vec, sMat_vec);
    }

    return prev_sMat_vec[T]; /// 4;
}


template<const int T>
//Matcher::result_t sw_sse2_word(simdNNAlign alignNNSIMD){
Matcher::result_t sw_sse2_word(const unsigned char *db_sequence, const int8_t ref_dir, const int32_t db_length,
                               const int32_t query_length, const uint8_t gap_open, const uint8_t gap_extend,
                               const simd_int *query_profile_word, const simd_int *query_profile_nn, const simd_int *query_profile_dist, simd_int *vHStore,
                               simd_int *vHLoad, simd_int *vE, simd_int *vHmax, uint8_t *maxColumn, simd_int *nnScoreVec,
                               const int8_t *subMatBiased, const short subMatBias,
                               const int8_t *distMatBiased, const short distMatBias, const char *dbNN,
                               const char *dbDist, const int nwGapPenalty, const uint16_t terminate, short nn_weight) {
#define max8(m, vm) ((m) = simdi16_hmax((vm)));

    uint16_t max = 0;		                     /* the max alignment score */
    int32_t end_read = query_length - 1;
    int32_t end_ref = 0; /* 1_based best alignment ending point; Initialized as isn't aligned - 0. */
    const unsigned int SIMD_SIZE = VECSIZE_INT * 2; // short element size 8 (SEE2) 16 (AVX2)
    int32_t segLen = (query_length + SIMD_SIZE-1) / SIMD_SIZE; /* number of segment */
    /* array to record the alignment read ending position of the largest score of each reference position */

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
//    memset(nnScoreVec, 0, segLen*sizeof(simd_int));
    int32_t i, j, k;

    // for division
    simd_int vb = simdi16_set(32768 / nn_weight);

    /* 16 byte insertion begin vector */
    simd_int vGapO = simdi16_set(gap_open);

    /* 16 byte insertion extension vector */
    simd_int vGapE = simdi16_set(gap_extend);

    /* 16 byte bias vector */
    simd_int vSubMatBias = simdi16_set(subMatBias);
    simd_int vDistMatBias = simdi16_set(distMatBias);
    //simd_int vBias = simdi16_set(-bias);    // set as a negative value for simd use
    simd_int vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
    simd_int vMaxMark = vZero; /* Trace the highest score till the previous column. */
    simd_int vTemp;
    int32_t begin = 0, end = db_length, step = 1;
    simd_int nwGapCPenaltyVec = simdi16_set(-nwGapPenalty);
    //fprintf(stderr, "start alignment of length %d [%d]\n", query_length, segLen * SIMD_SIZE);

    /* outer loop to process the reference sequence */
    if (ref_dir == 1) {
        begin = db_length - 1;
        end = -1;
        step = -1;
    }

    // iterate over db_sequence (target seq)
    for (i = begin; LIKELY(i != end); i += step) {
        simd_int e, vF = vZero, vMaxColumn = vZero; /* Initialize F value to 0.
                                Any errors to vH values will be corrected in the Lazy_F loop.
                                */

        simd_int vH = pvHStore[segLen - 1];
        vH = simdi8_shiftl (vH, 2); /* Shift the 128-bit value in vH left by 2 byte. */
        const simd_int* vP = query_profile_word + db_sequence[i] * segLen; /* Right part of the query_profile_byte, query_profile_word is a pointer */
        simd_int targetSubMat[T];
        simd_int targetDistMat[T];
        for(int x = 0; x < T; x++){
            targetSubMat[x] = simdi_load((simd_int*)&subMatBiased[dbNN[i*T+x] * 32]);
            targetDistMat[x] = simdi_load((simd_int*)&distMatBiased[dbDist[i*T+x] * 32]);
        }

        /* Swap the 2 H buffers. */
        simd_int* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        for (j = 0; LIKELY(j < segLen); j ++) {
            const simd_int* vP_nn = query_profile_nn + j * T;
            const simd_int* subQ_dist = query_profile_dist + j * T;
//             simd_int nn_score = simdi32_set(0);
// todo: add scaling
            simd_int nn_score = needlemanWunschScoreVec<T>(vP_nn, &targetSubMat[0], subQ_dist, &targetDistMat[0], nwGapCPenaltyVec, vSubMatBias, vDistMatBias);

            nn_score = _mm256_mulhrs_epi16(nn_score, vb);

            uint16_t val[8];
            memcpy(val, &nn_score, sizeof(val));

            simdi_store(nnScoreVec+j, nn_score);
        }
        /* inner loop to process the query sequence */
        for (j = 0; LIKELY(j < segLen); j ++) {
            simd_int score = simdi_load(vP + j);
            simd_int nnScore = simdi_load(nnScoreVec + j);
            score = simdi16_adds(score, nnScore);
            // vH = vH + score;
            vH = simdi16_adds(vH, score);

            /* Get max from vH, vE and vF. */
            e = simdi_load(pvE + j);
            vH = simdi16_max(vH, e);
            vH = simdi16_max(vH, vF);


            vMaxColumn = simdi16_max(vMaxColumn, vH);

            /* Save vH values. */
            simdi_store(pvHStore + j, vH);

            /* Update vE value. */
            vH = simdui16_subs(vH, vGapO); /* saturation arithmetic, result >= 0 */
            e = simdui16_subs(e, vGapE);
            e = simdi16_max(e, vH);
            simdi_store(pvE + j, e);

            /* Update vF value. */
            vF = simdui16_subs(vF, vGapE);
            vF = simdi16_max(vF, vH);


            /* Load the next vH. */
            vH = simdi_load(pvHLoad + j);
        }

        /* Lazy_F loop: has been revised to disallow adjacent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
        for (k = 0; LIKELY(k < (int32_t) SIMD_SIZE); ++k) {
            vF = simdi8_shiftl (vF, 2);
            for (j = 0; LIKELY(j < segLen); ++j) {
                vH = simdi_load(pvHStore + j);
                vH = simdi16_max(vH, vF);
                vMaxColumn = simdi16_max(vMaxColumn, vH); //newly added line
                simdi_store(pvHStore + j, vH);
                vH = simdui16_subs(vH, vGapO);
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

    Matcher::result_t result;
    /* Find the most possible 2nd best alignment. */
    result.score = max;

    return result;
#undef max8
}



/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch. */
template <typename T, size_t Elements>
void createQueryProfile(simd_int *profile, const int8_t *query_sequence, const int8_t *mat,
                        const int32_t query_length, const int32_t aaSize, uint8_t bias) {
    const int32_t segmentLen = (query_length + Elements - 1) / Elements;
    T* t = (T*) profile;
    for (int32_t nt = 0; LIKELY(nt < aaSize); nt++) {
        for (int32_t i = 0; i < segmentLen; i++) {
            int32_t  j = i;
            for (size_t segmentNum = 0; LIKELY(segmentNum < Elements) ; segmentNum++) {
                *t = ( j >= query_length) ? bias : mat[nt * aaSize + query_sequence[j]] + bias; // mat[nt][q[j]] mat eq 20*20
                t++;
                j += segmentLen;
            }
        }
    }
}


/* Generate query profile for nearest neighbors. */
template <typename T, size_t Elements, const int NN>
void createNNQueryProfile(simd_int *profile, simd_int *profile_dist, const char *nn_sequence, const char *nn_dist, const int32_t query_length) {
    const int32_t segLen = (query_length + Elements - 1) / Elements;
    T* nn3Di = (T*) profile;
    T* nnDist = (T*) profile_dist;
    size_t pos = 0;
    for (int32_t i = 0; i < segLen; i++) {
        int32_t j = i;
        for (size_t segNum = 0; LIKELY(segNum < Elements); segNum++) {
            // beyond the length of query so pad with neutral consensus values
            for(int x = 0; x < NN; x++){
                *(nn3Di+pos+segNum+Elements*x) = (j >= query_length) ? 20 : nn_sequence[j*NN+x];
                *(nnDist+pos+segNum+Elements*x) = (j >= query_length) ? 20 : nn_dist[j*NN+x];
            }
//            *t++ = (j >= query_length) ? 20 :  consens_sequence[j + (offset - 1)];
            j += segLen;
        }
        pos+=Elements*NN;
    }
    return;
}

// static two dimensional array instead of vectors
template<const int T>
void findNearestNeighbour(char * nn, char * dist, Coordinates & ca,
                          int length, unsigned char * seqInt,
                          std::pair<float, int> * seqDistList, float * seqDistSIMD){
    // (char*)mem_align(ALIGN_INT, par.maxSeqLen * sizeof(char) * 4 );

    // calculate pairwise matrix (euclidean distance)
    std::pair<int,float> pos[T];
    std::pair<int,float> pos2[T];
    char closestStateI;
    int minDistance;
    int seqDistance;

    for(int i = 0; i < length; i++){
        size_t seqDistListPos = 0;
        size_t seqDistListPos2 = 0;

//        //SIMD implementation of calculating distance
        simd_float coords_x_i = simdf32_set(ca.x[i]);
        simd_float coords_y_i = simdf32_set(ca.y[i]);
        simd_float coords_z_i = simdf32_set(ca.z[i]);

        for (int j = 0; j < length; j+=8){
            simd_float loaded_x_j = simdf32_load(&ca.x[j]);
            simd_float loaded_y_j = simdf32_load(&ca.y[j]);
            simd_float loaded_z_j = simdf32_load(&ca.z[j]);

            simd_float dist_simd = simdf32_mul(simdf32_sub(coords_x_i, loaded_x_j), simdf32_sub(coords_x_i, loaded_x_j));
            dist_simd = simdf32_add(dist_simd, simdf32_mul(simdf32_sub(coords_y_i, loaded_y_j), simdf32_sub(coords_y_i, loaded_y_j)));
            dist_simd = simdf32_add(dist_simd, simdf32_mul(simdf32_sub(coords_z_i, loaded_z_j), simdf32_sub(coords_z_i, loaded_z_j)));
            dist_simd = simdf32_sqrt(dist_simd);
            simdf32_store(&seqDistSIMD[j], dist_simd);
        }

        for(int a = 0; a < length; a++){
            if(i != a) {
                seqDistList[seqDistListPos2].first = seqDistSIMD[a];
                seqDistList[seqDistListPos2].second = a;
                seqDistListPos2++;
            }
        }

        // find the T nearest neighbours for each amino acid
        std::partial_sort(seqDistList, seqDistList + T, seqDistList + length - 1);

        for(int m  = 0; m < T; m++){
            pos2[m].first = seqDistList[m].second;
            pos2[m].second = seqDistList[m].first;
        }
        std::sort(pos2, pos2+T);

        for(int n = 0; n < T; n++){
            int neighbour = static_cast<int>(seqInt[pos2[n].first]);
            nn[i*T + n] = neighbour;

            // discretise the distance within the sequence position
            seqDistance = pos2[n].first - i;

            minDistance = INT_MAX;
            closestStateI = Alphabet3diSeqDist::INVALID_STATE;
            for (size_t a = 0; a < Alphabet3diSeqDist::CENTROID_CNT; a++){

                int distToCentroid = abs(Alphabet3diSeqDist::centroids[a] - seqDistance);
                if (distToCentroid < minDistance){
                    closestStateI = a;
                    minDistance = distToCentroid;
                }
            }
            dist[i*T + n] = closestStateI;
        }
    }
}

template<const int T>
short needlemanWunschScore(int subQNNi[T], int subTNNi[T], int subQNNdist[T], int subTNNdist[T], SubstitutionMatrix *subMat, int nwGapPenalty, int m){

    int scoringMatrix[(T+1)*(T+1)]; // = {{0, 0, 0, 0, 0},{0, 0, 0, 0, 0},{0, 0, 0, 0, 0},{0, 0, 0, 0, 0},{0, 0, 0, 0, 0}};
    memset(scoringMatrix, 0, (T+1)*(T+1) * sizeof(int));

    int distMatValue;
    short combinedScore;

    for(int i = 1; i < (T+1); i++){
        for(int j = 1; j < (T+1); j++){

            int score1 = subQNNi[i - 1];
            int score2 = subTNNi[j - 1];
            short scoreTest;

            scoreTest = subMat->subMatrix[score1][score2];

            distMatValue = distMat[subQNNdist[i - 1]][subTNNdist[j - 1]];

            float distMatF = distMatValue/m;
            short distMatS = static_cast<int>(distMatF);

            combinedScore = scoreTest + distMatS;

            // calculate scores - bonus for base pairing and penalty for not
            scoringMatrix[i*(T+1) + j] = std::max(scoringMatrix[(i - 1)*(T+1) + j - 1] + combinedScore, scoringMatrix[(i - 1)*(T+1) + j] - nwGapPenalty);
            scoringMatrix[i*(T+1) + j] = std::max(scoringMatrix[i*(T+1) + j - 1] - nwGapPenalty, scoringMatrix[i*(T+1) + j]);
        }
    }
    short score = scoringMatrix[(T+1)*(T+1)-1]; /// 4;

    return score;
}

template<const int T>
Matcher::result_t alignByNN(char * querynn, unsigned char *querySeqInt, int queryLen, char * queryNNdist, char * targetnn, unsigned char *targetSeqInt, int targetLen, char * targetNNdist, SubstitutionMatrix *subMat, int gapOpen, int gapExtern, int gapNW, int nnWeight, int m){

    Matcher::result_t result;

    // gotoh sw itself

    struct scores{ short H, E, F; };
    uint16_t max_score = 0;

    scores *workspace = new scores[queryLen * 2 + 2];
    scores *curr_sHEF_vec = &workspace[0];
    scores *prev_sHEF_vec = &workspace[queryLen + 1];

    int max = 0; int min = 0;

    int subTNNi[T], subQNNi[T];
    int subTNNdist[T], subQNNdist[T];
    // top row need to be set to a 0 score
    memset(prev_sHEF_vec, 0, sizeof(scores) * (queryLen + 1));
    for (int i = 0; i < targetLen; i++) {

        for(int a = 0; a < T; a++){
            int neighbour = static_cast<int>(targetnn[T*i + a]);
            subTNNi[a] = neighbour;
            subTNNdist[a] = targetNNdist[T*i + a];
        }

        // left outer column need to be set to a 0 score
        prev_sHEF_vec[0].H = 0; prev_sHEF_vec[0].E = 0; prev_sHEF_vec[0].F = 0;
        curr_sHEF_vec[0].H = 0; curr_sHEF_vec[0].E = 0; curr_sHEF_vec[0].E = 0;
        for (int j = 1; j <= queryLen; j++) {
            for(int a = 0; a < T; a++){
                int neighbour = static_cast<int>(querynn[T*(j-1) + a]);
                subQNNi[a] = neighbour;
                subQNNdist[a] = queryNNdist[T*(j-1) + a];
            }

            curr_sHEF_vec[j].E = std::max(curr_sHEF_vec[j - 1].H - gapOpen, curr_sHEF_vec[j - 1].E - gapExtern); // j-1
            curr_sHEF_vec[j].F = std::max(prev_sHEF_vec[j].H - gapOpen, prev_sHEF_vec[j].F - gapExtern); // i-1
            short nnScore = needlemanWunschScore<T>(subTNNi, subQNNi, subTNNdist, subQNNdist, subMat, gapNW, m);

            int subOne = static_cast<int>(targetSeqInt[i]);
            int subTwo = static_cast<int>(querySeqInt[j - 1]);

//            float nnScoreWeighted = nnScore * nnWeight;
//            nnScoreWeighted += (nnScoreWeighted < 0) ? -0.5 : 0.5;
            if(nnScore > max){max = nnScore;}
            if(nnScore < min){min = nnScore;}
            short nnScoreWeighted = nnScore * nnWeight;

            const short tempH = prev_sHEF_vec[j - 1].H + subMat->subMatrix[subOne][subTwo] + nnScoreWeighted;
//            const short tempH = prev_sHEF_vec[j - 1].H + subMat->subMatrix[subOne][subTwo] +  static_cast<short>(nnScoreWeighted); // i - 1, j - 1


            curr_sHEF_vec[j].H = std::max(tempH, curr_sHEF_vec[j].E);
            curr_sHEF_vec[j].H = std::max(curr_sHEF_vec[j].H, curr_sHEF_vec[j].F);
            curr_sHEF_vec[j].H = std::max(curr_sHEF_vec[j].H, static_cast<short>(0));
            max_score = static_cast<uint16_t>(std::max(static_cast<uint16_t>(curr_sHEF_vec[j].H), max_score));
        }
        std::swap(prev_sHEF_vec, curr_sHEF_vec);
    }
    delete [] workspace;
    result.score = max_score;
    return result;
}


int pareunaligner(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, 0, MMseqsParameter::COMMAND_ALIGN);

    Debug(Debug::INFO) << "Sequence database: " << par.db1 << "\n";
    DBReader<unsigned int> qdbr((par.db1+"_ss").c_str(), (par.db1+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    qdbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> qcadbr((par.db1+"_ca").c_str(), (par.db1+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    qcadbr.open(DBReader<unsigned int>::NOSORT);

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);

    DBReader<unsigned int> *tdbr = NULL;
    DBReader<unsigned int> *tcadbr = NULL;

    bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    bool sameDB = false;

    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
        tcadbr = &qcadbr;
    } else {
        tdbr = new DBReader<unsigned int>((par.db2+"_ss").c_str(), (par.db2+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        tdbr->open(DBReader<unsigned int>::NOSORT);
        tcadbr = new DBReader<unsigned int>((par.db2+"_ca").c_str(), (par.db2+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        tcadbr->open(DBReader<unsigned int>::NOSORT);
        if (touch) {
            tdbr->readMmapedDataInMemory();
            tcadbr->readMmapedDataInMemory();
        }
    }

    Debug(Debug::INFO) << "Result database: " << par.db3 << "\n";
    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    Debug(Debug::INFO) << "Output file: " << par.db4 << "\n";
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();

    //temporary output file
    Debug::Progress progress(resultReader.getSize());

    // sub. mat needed for query profile
    int8_t * tinySubMat = (int8_t*) mem_align(ALIGN_INT, subMat.alphabetSize * 32);
    int8_t * tinySubMatBias = (int8_t*) mem_align(ALIGN_INT, subMat.alphabetSize * 32);
    int8_t * tinyDistMat = (int8_t*) mem_align(ALIGN_INT, DIST_MAT_SIZE * 32); // todo: make this a parameter later
    int8_t * tinyDistMatBias = (int8_t*) mem_align(ALIGN_INT, DIST_MAT_SIZE * 32);
    short subMatBias = 0;
    short distMatBias = 0;
    // to find the minimum value of the sub mat to add later so everything is pos
    for (int i = 0; i < subMat.alphabetSize; i++) {
        for (int j = 0; j < subMat.alphabetSize; j++) {
            subMatBias = std::min(subMat.subMatrix[i][j], subMatBias);
        }
    }
    subMatBias = abs(subMatBias);
    for (int i = 0; i < subMat.alphabetSize; i++) {
        memset(tinySubMatBias+i*32, subMatBias, sizeof(char) * 32);
        for (int j = 0; j < subMat.alphabetSize; j++) {
            tinySubMat[i*subMat.alphabetSize + j] = subMat.subMatrix[i][j]; // for farrar profile
            tinySubMatBias[i*32 + j] = subMatBias + subMat.subMatrix[i][j]; // for inner loop of NW
        }
    }

    for(int m = 0; m < DIST_MAT_SIZE; m++){
        for(int n = 0; n < DIST_MAT_SIZE; n++){
            float scaled = distMat[m][n]/par.slope;
            scaled += (scaled < 0) ? -0.5 : 0.5;
            distMatBias = std::min(static_cast<short>(scaled), distMatBias);
        }
    }

    distMatBias = abs(distMatBias);
    for (int i = 0; i < DIST_MAT_SIZE; i++) {
        memset(tinyDistMatBias+i*32, distMatBias, sizeof(char) * 32);
        for (int j = 0; j < DIST_MAT_SIZE; j++) {
            float scaled = distMat[i][j]/par.slope;
            scaled += (scaled < 0) ? -0.5 : 0.5;
            tinyDistMat[i*DIST_MAT_SIZE + j] = static_cast<short>(scaled); // for farrar profile
            tinyDistMatBias[i*32 + j] = distMatBias + static_cast<short>(scaled); // for inner loop of NW
        }
    }


#pragma omp parallel
    {

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        const int segmentSize = (par.maxSeqLen+7)/8;
        float * query_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        simd_int* query_profile_word = (simd_int*)mem_align(ALIGN_INT, subMat.alphabetSize * segmentSize * sizeof(simd_int));
        simd_int* query_profile_nn  = (simd_int*)mem_align(ALIGN_INT, 6 * segmentSize * sizeof(simd_int));
        simd_int* query_profile_dist  = (simd_int*)mem_align(ALIGN_INT, 6 * segmentSize * sizeof(simd_int));

        float * target_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        Sequence qSeq(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat, 0, false, par.compBiasCorrection);
        Sequence tSeq(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat, 0, false, par.compBiasCorrection);

        char * querynn  = (char*)mem_align(ALIGN_INT, par.maxSeqLen * sizeof(char) * 8 ); // 8 to number of neighbours
        char * targetnn = (char*)mem_align(ALIGN_INT, par.maxSeqLen * sizeof(char) * 8 );
        char * querynn_dist  = (char*)mem_align(ALIGN_INT, par.maxSeqLen * sizeof(char) * 8 );
        char * targetnn_dist = (char*)mem_align(ALIGN_INT, par.maxSeqLen * sizeof(char) * 8 );


        simd_int* vHStore = (simd_int*) mem_align(ALIGN_INT, segmentSize * sizeof(simd_int));
        simd_int* vHLoad  = (simd_int*) mem_align(ALIGN_INT, segmentSize * sizeof(simd_int));
        simd_int* vE      = (simd_int*) mem_align(ALIGN_INT, segmentSize * sizeof(simd_int));
        simd_int* vHmax   = (simd_int*) mem_align(ALIGN_INT, segmentSize * sizeof(simd_int));
        simd_int* nnScoreVec   = (simd_int*) mem_align(ALIGN_INT, segmentSize * sizeof(simd_int));
        uint8_t * maxColumn = new uint8_t[par.maxSeqLen*sizeof(uint16_t)];
        std::pair<float, int> * seqDistList = new std::pair<float, int>[par.maxSeqLen];
        float * seqDistSimd = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float));
        std::vector<Matcher::result_t> alignmentResult;
        PareunAlign paruenAlign(par.maxSeqLen, &subMat); // subMat called once, don't need to call it again?
        char buffer[1024+32768];
        std::string resultBuffer;

        // write output file

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            char *data = resultReader.getData(id, thread_idx);
            size_t queryKey = resultReader.getDbKey(id);
            if(*data != '\0') {
                unsigned int queryId = qdbr.getId(queryKey);

                char *querySeq = qdbr.getData(queryId, thread_idx);

                unsigned int querySeqLen = qdbr.getSeqLen(queryId);
                qSeq.mapSequence(id, queryKey, querySeq, querySeqLen);
                // SSE2 = 8 shorts
                // AVX2 = 16 shorts
                float *qdata = (float *) qcadbr.getData(queryId, thread_idx);
                Coordinates queryCaCords;
                memcpy(query_x, qdata, sizeof(float) * qSeq.L);
                memcpy(query_y, &qdata[qSeq.L], sizeof(float) * qSeq.L);
                memcpy(query_z, &qdata[qSeq.L+qSeq.L], sizeof(float) * qSeq.L);

                queryCaCords.x = query_x;
                queryCaCords.y = query_y;
                queryCaCords.z = query_z;
//
                switch (par.numberNN) {
                    case 3:
                        findNearestNeighbour<3>(querynn, querynn_dist, queryCaCords, querySeqLen, qSeq.numSequence, seqDistList, seqDistSimd);
                        createNNQueryProfile<int16_t, VECSIZE_INT * 2, 3>(query_profile_nn, query_profile_dist, (const char*) querynn,(const char*) querynn_dist,  qSeq.L);
                        break;
                    case 4:
                        findNearestNeighbour<4>(querynn, querynn_dist, queryCaCords, querySeqLen, qSeq.numSequence, seqDistList, seqDistSimd);
                        createNNQueryProfile<int16_t, VECSIZE_INT * 2, 4>(query_profile_nn, query_profile_dist, (const char*) querynn,(const char*) querynn_dist,  qSeq.L);
                        break;
                    case 5:
                        findNearestNeighbour<5>(querynn, querynn_dist, queryCaCords, querySeqLen, qSeq.numSequence, seqDistList, seqDistSimd);
                        createNNQueryProfile<int16_t, VECSIZE_INT * 2, 5>(query_profile_nn, query_profile_dist, (const char*) querynn,(const char*) querynn_dist,  qSeq.L);
                        break;
                    case 6:
                        findNearestNeighbour<6>(querynn, querynn_dist, queryCaCords, querySeqLen, qSeq.numSequence, seqDistList, seqDistSimd);
                        createNNQueryProfile<int16_t, VECSIZE_INT * 2, 6>(query_profile_nn, query_profile_dist, (const char*) querynn,(const char*) querynn_dist,  qSeq.L);
                        break;
                    case 7:
                        findNearestNeighbour<7>(querynn, querynn_dist, queryCaCords, querySeqLen, qSeq.numSequence, seqDistList, seqDistSimd);
                        createNNQueryProfile<int16_t, VECSIZE_INT * 2, 7>(query_profile_nn, query_profile_dist, (const char*) querynn,(const char*) querynn_dist,  qSeq.L);
                        break;
                    case 8:
                        findNearestNeighbour<8>(querynn, querynn_dist, queryCaCords, querySeqLen, qSeq.numSequence, seqDistList, seqDistSimd);
                        createNNQueryProfile<int16_t, VECSIZE_INT * 2, 8>(query_profile_nn, query_profile_dist, (const char*) querynn,(const char*) querynn_dist,  qSeq.L);
                        break;
                }

                createQueryProfile<int16_t, VECSIZE_INT * 2>(query_profile_word, (const int8_t*) qSeq.numSequence, tinySubMat, qSeq.L, subMat.alphabetSize, 0);

                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = tdbr->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;

                    char * targetSeq = tdbr->getData(targetId, thread_idx);

                    int targetLen = static_cast<int>(tdbr->getSeqLen(targetId));
                    float * tdata = (float*) tcadbr->getData(targetId, thread_idx);
                    tSeq.mapSequence(targetId, dbKey, targetSeq,tdbr->getSeqLen(targetId));

                    if(Util::canBeCovered(par.covThr, par.covMode, qSeq.L, targetLen)==false){
                        continue;
                    }
                    Coordinates targetCaCords;
                    memcpy(target_x, tdata, sizeof(float) * tSeq.L);
                    memcpy(target_y, &tdata[tSeq.L], sizeof(float) * tSeq.L);
                    memcpy(target_z, &tdata[tSeq.L+tSeq.L], sizeof(float) * tSeq.L);

                    targetCaCords.x = target_x;
                    targetCaCords.y = target_y;
                    targetCaCords.z = target_z;


                    Matcher::result_t res;
                    // switch case
                    switch (par.numberNN) {
                        case 3:
                            findNearestNeighbour<3>(targetnn, targetnn_dist, targetCaCords, tSeq.L, tSeq.numSequence, seqDistList, seqDistSimd);
                            res = sw_sse2_word<3>(tSeq.numSequence, 0, tSeq.L, qSeq.L, par.gapOpen.values.aminoacid(),
                                                  par.gapExtend.values.aminoacid(), query_profile_word, query_profile_nn, query_profile_dist,
                                                  vHStore, vHLoad, vE, vHmax, maxColumn, nnScoreVec, tinySubMatBias, subMatBias, tinyDistMatBias, distMatBias, targetnn, targetnn_dist,
                                                  par.gapNW, SHRT_MAX, par.nnWeight);
                            break;
                        case 4:
                            findNearestNeighbour<4>(targetnn, targetnn_dist, targetCaCords, tSeq.L, tSeq.numSequence, seqDistList, seqDistSimd);
                            res = sw_sse2_word<4>(tSeq.numSequence, 0, tSeq.L, qSeq.L, par.gapOpen.values.aminoacid(),
                                                  par.gapExtend.values.aminoacid(), query_profile_word, query_profile_nn, query_profile_dist,
                                                  vHStore, vHLoad, vE, vHmax, maxColumn, nnScoreVec, tinySubMatBias, subMatBias, tinyDistMatBias, distMatBias, targetnn, targetnn_dist,
                                                  par.gapNW, SHRT_MAX, par.nnWeight);
                            break;
                        case 5:
                            findNearestNeighbour<5>(targetnn, targetnn_dist, targetCaCords, tSeq.L, tSeq.numSequence, seqDistList, seqDistSimd);
                            res = sw_sse2_word<5>(tSeq.numSequence, 0, tSeq.L, qSeq.L, par.gapOpen.values.aminoacid(),
                                                  par.gapExtend.values.aminoacid(), query_profile_word, query_profile_nn, query_profile_dist,
                                                  vHStore, vHLoad, vE, vHmax, maxColumn, nnScoreVec, tinySubMatBias, subMatBias, tinyDistMatBias, distMatBias, targetnn, targetnn_dist,
                                                  par.gapNW, SHRT_MAX, par.nnWeight);
                            break;
                        case 6:
                            findNearestNeighbour<6>(targetnn, targetnn_dist, targetCaCords, tSeq.L, tSeq.numSequence, seqDistList, seqDistSimd);
                            res = sw_sse2_word<6>(tSeq.numSequence, 0, tSeq.L, qSeq.L, par.gapOpen.values.aminoacid(),
                                                  par.gapExtend.values.aminoacid(), query_profile_word, query_profile_nn, query_profile_dist,
                                                  vHStore, vHLoad, vE, vHmax, maxColumn, nnScoreVec, tinySubMatBias, subMatBias, tinyDistMatBias, distMatBias, targetnn, targetnn_dist,
                                                  par.gapNW, SHRT_MAX, par.nnWeight);
                            break;
                        case 7:
                            findNearestNeighbour<7>(targetnn, targetnn_dist, targetCaCords, tSeq.L, tSeq.numSequence, seqDistList, seqDistSimd);
                            res = sw_sse2_word<7>(tSeq.numSequence, 0, tSeq.L, qSeq.L, par.gapOpen.values.aminoacid(),
                                                  par.gapExtend.values.aminoacid(), query_profile_word, query_profile_nn, query_profile_dist,
                                                  vHStore, vHLoad, vE, vHmax, maxColumn, nnScoreVec, tinySubMatBias, subMatBias, tinyDistMatBias, distMatBias, targetnn, targetnn_dist,
                                                  par.gapNW, SHRT_MAX, par.nnWeight);
                            break;
                        case 8:
                            findNearestNeighbour<8>(targetnn, targetnn_dist, targetCaCords, tSeq.L, tSeq.numSequence, seqDistList, seqDistSimd);
                            res = sw_sse2_word<8>(tSeq.numSequence, 0, tSeq.L, qSeq.L, par.gapOpen.values.aminoacid(),
                                                  par.gapExtend.values.aminoacid(), query_profile_word, query_profile_nn, query_profile_dist,
                                                  vHStore, vHLoad, vE, vHmax, maxColumn, nnScoreVec, tinySubMatBias, subMatBias, tinyDistMatBias, distMatBias, targetnn, targetnn_dist,
                                                  par.gapNW, SHRT_MAX, par.nnWeight);
                            break;
                    }

                    unsigned int targetKey = tdbr->getDbKey(targetId);

                    res.dbKey = targetKey;
                    res.eval = 0;
                    //Matcher::result_t res = paruenAlign.align(qSeq, tSeq, &subMat, evaluer);
                    string cigar =  paruenAlign.backtrace2cigar(res.backtrace);

                    if (Alignment::checkCriteria(res, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                        alignmentResult.emplace_back(res);
                        passedNum++;
                        rejected = 0;
                    }

                    else {
                        rejected++;
                    }
                }
            }
            if (alignmentResult.size() > 1) {
                SORT_SERIAL(alignmentResult.begin(), alignmentResult.end(), Matcher::compareHits);
            }
            for (size_t result = 0; result < alignmentResult.size(); result++) {
                size_t len = Matcher::resultToBuffer(buffer, alignmentResult[result], par.addBacktrace);
                resultBuffer.append(buffer, len);
            }
            dbw.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
            resultBuffer.clear();
            alignmentResult.clear();
        }
        free(target_x);
        free(target_y);
        free(target_z);
        free(query_x);
        free(query_y);
        free(query_z);
        free(querynn);
        free(targetnn);
        free(query_profile_word);
        free(vHStore);
        free(vHLoad);
        free(vE);
        free(vHmax);
        free(query_profile_nn);
        free(seqDistSimd);
        free(nnScoreVec);
        free(query_profile_dist);
        free(querynn_dist);
        free(targetnn_dist);
        delete [] maxColumn;
        delete [] seqDistList;
    }
    delete [] tinySubMat;
    delete [] tinySubMatBias;
    delete [] tinyDistMat;
    delete [] tinyDistMatBias;

    dbw.close();
    resultReader.close();
    qdbr.close();
    qcadbr.close();
    if(sameDB == false){
        tdbr->close();
        tcadbr->close();
        delete tdbr;
        delete tcadbr;
    }
    return EXIT_SUCCESS;
}
