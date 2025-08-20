#include "LoLAlign.h"
#include "Fwbw.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "SubstitutionMatrix.h"
#include "Matcher.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Coordinate16.h"
#include "StructureUtil.h"

#include <iostream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <numeric>
#include <cmath>
#include <vector>
#include <fstream>

#ifdef OPENMP
#include <omp.h>
#endif

using namespace std;

bool compareHitsBylolScore(const Matcher::result_t &first, const Matcher::result_t &second) {
    return first.eval > second.eval;
}

lolAlign::lolAlign(unsigned int maxSeqLen, bool computeExactScore)
        : xtm(maxSeqLen), ytm(maxSeqLen), xt(maxSeqLen),
          r1(maxSeqLen), r2(maxSeqLen), computeExactScore(computeExactScore)
{
    query_x = (float *)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    query_y = (float *)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    query_z = (float *)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    target_x = (float *)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    target_y = (float *)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    target_z = (float *)mem_align(ALIGN_FLOAT, maxSeqLen * sizeof(float));
    d_ij = malloc_matrix<float>(maxSeqLen, maxSeqLen);
    d_kl = malloc_matrix<float>(maxSeqLen, maxSeqLen);
    //P = malloc_matrix<float>(maxSeqLen, maxSeqLen);
}



lolAlign::~lolAlign()
{
    free(query_x);
    free(query_y);
    free(query_z);
    free(target_x);
    free(target_y);
    free(target_z);
    free(d_ij);
    free(d_kl);
    free(P);
    free(G);
    free(anchor_query);
    free(anchor_target);
    free(hidden_layer);
    delete[] sa_index;
    delete[] sa_scores;
    delete[] anchor_length;
    delete[] new_anchor_length;
    delete[] gaps;
    free(lol_dist);
    free(lol_seq_dist);
    free(lol_score_vec);
    delete[] final_anchor_query;
    delete[]final_anchor_target;
    free(lol_score_vec_sh);
}


void lolAlign::reallocate_target(size_t targetL){
    free(d_kl);
    d_kl = malloc_matrix<float>(targetL, targetL);
    free(G);
    G = malloc_matrix<float>(queryLen, targetL);
    free(P);
    P = malloc_matrix<float>(queryLen, targetL);
    free(hidden_layer);
    hidden_layer = malloc_matrix<float>(targetL, 3);
    free(anchor_target);
    anchor_target = malloc_matrix<int>(num_sa, targetL);
    free(lol_dist);
    lol_dist = (float *)mem_align(ALIGN_FLOAT, targetL * sizeof(float));
    free(lol_seq_dist);
    lol_seq_dist = (float *)mem_align(ALIGN_FLOAT, targetL * sizeof(float));
    free(lol_score_vec);
    lol_score_vec = (float *)mem_align(ALIGN_FLOAT, targetL * sizeof(float));
    free(lol_score_vec_sh);
    lol_score_vec_sh = (float *)mem_align(ALIGN_FLOAT, targetL * sizeof(float));
    delete[] final_anchor_target;
    final_anchor_target = new int[targetL];
    delete[] final_anchor_query;
    final_anchor_query = new int[queryLen];
}

void lolAlign::calc_gap(int* anchor_query, int* anchor_target, int * gaps,  int queryLen, int targetLen)
{
    gaps[0] = -1;
    int index_q = gaps[1];
    int index_t = gaps[3];
    while(anchor_query[index_q] != 0 || anchor_target[index_t] != 0)
    {
        if (anchor_query[index_q] != 0 && anchor_target[index_t] != 0)
        {
            index_q++;
            index_t++;
        }
        else if (anchor_query[index_q] == 0)
        {
            index_q++;
        }
        else if (anchor_target[index_t] == 0)
        {
            index_t++;
        }
        if (index_q == queryLen || index_t == targetLen)
        {
            gaps[0] = -1;
            return;
        }
    }
    gaps[0] = index_q;
    gaps[2] = index_t;
    while(anchor_query[index_q] == 0 || anchor_target[index_t] == 0){

        if(anchor_query[index_q] == 0 && anchor_target[index_t] == 0){
            index_q++;
            index_t++;
        }
        else if(anchor_query[index_q] == 0){
            index_q++;
        }
        else if(anchor_target[index_t] == 0){
            index_t++;
        }
        if (index_q == queryLen || index_t == targetLen){
            index_q = queryLen;
            index_t = targetLen;
            break;
        }
    }
    gaps[1] = index_q;
    gaps[3] = index_t;
    return;

}
float lolAlign::maxSubArray(float* nums, int numsSize) {

    float currentMax = nums[0];
    float globalMax = nums[0];

    for (int i = 1; i < numsSize; ++i) {
        currentMax = std::max(nums[i], currentMax + nums[i]);
        globalMax = std::max(globalMax, currentMax);
    }

    return globalMax;
}

void lolAlign::index_sort(float* nums, int* index, int numsSize) {
    std::sort(index, index + numsSize, [&nums](int i1, int i2) { return nums[i1] < nums[i2]; });
}




Matcher::result_t lolAlign::align(unsigned int dbKey, float *target_x, float *target_y, float *target_z,
                                  Sequence& tSeqAA, Sequence& tSeq3Di, int targetLen, SubstitutionMatrix &subMatAA, FwBwAligner* fwbwaln, int md)
{
    max_lolmat_idx = 0;
    min_lolmat_idx = queryLen;

    unsigned char *targetNumAA = tSeqAA.numSequence;
    unsigned char *targetNum3Di = tSeq3Di.numSequence;
    const char *targetSeq = tSeqAA.getSeqData();
    lolAlign::computeForwardScoreMatrix(
            targetNumAA,
            targetNum3Di,
            queryLen,
            targetLen,
            subMatAA,
            G);



    int maxIndexX = 0;
    int maxIndexY = 0;
    
    for(int i = 0; i < num_sa; i++){
        sa_index[i] = i;
        sa_scores[i] = 0;
    }

    calc_dist_matrix(target_x, target_y, target_z, targetLen, d_kl, false);

    for(int sa = 0; sa < num_sa; sa++){
        anchor_length[sa] = 0;
        new_anchor_length[sa] = 0;
        for (int i = 0; i < queryLen; i++) {
            anchor_query[sa][i] = 0;
        }
        for (int i = 0; i < targetLen; i++) {
            anchor_target[sa][i] = 0;
        }
    }

    gaps[0] = 0;
    gaps[1] = queryLen;
    gaps[2] = 0;
    gaps[3] = targetLen;

    fwbwaln->resetParams(start_anchor_go, start_anchor_ge, start_anchor_T);
    fwbwaln->initScoreMatrix(G, gaps);
    fwbwaln->runFwBw<false, 0>();
    
    for(int sa = 0; sa < num_sa; sa++){

        /*if(sa % 5 == 0){
            fwbwaln->initScoreMatrix(G, gaps);
            fwbwaln->computeProbabilityMatrix<0>();

        }*/
        maxIndexX = 0;
        maxIndexY = 0;

        float maxScore = 0;

        for (int i = start_anchor_length; i < queryLen - start_anchor_length ; ++i) {
            for (int j = start_anchor_length; j < targetLen - start_anchor_length ; ++j) {
                float alignProb = fwbwaln->getProbability(i, j);
                if(i%2 == 0){
                    if (alignProb >= maxScore) {
                        maxScore = alignProb;
                        maxIndexX = i;
                        maxIndexY = j;
                    }
                }
                else{
                    if (alignProb >= maxScore) {
                        maxScore = alignProb;
                        maxIndexX = i;
                        maxIndexY = j;
                    }
                }

            }
        }



        int start_row = maxIndexX - std::min(maxIndexX, maxIndexY);
        int start_col = maxIndexY - std::min(maxIndexX, maxIndexY);
        int diag_length = std::min(queryLen - start_row, targetLen - start_col);
        for(int i = 0; i < diag_length; i++){
            lol_score_vec[i] = G[start_row + i][start_col + i];
        }
        lol_score_vec[std::min(maxIndexX, maxIndexY)] += 200;
         
        for(int i = -start_anchor_length; i < start_anchor_length; i++){
            for(int j=0; j<diag_length; j++){

                if(d_ij[maxIndexX+i][start_row + j] > 0){

                    lol_dist[j] = std::abs(d_ij[maxIndexX+i][start_row + j] - d_kl[maxIndexY+i][start_col+j]);
                    lol_seq_dist[j] =  std::copysign(1.0f, (maxIndexX+i - start_row + j)) * std::log(1 + std::abs((float)(maxIndexX+i - start_row + j)));

                    //std::cout << "diag_row_dist[" << j<< "]: " << diag_row_dist[j] - diag_col_dist[j] << std::endl;
                }
                else{
                    lol_dist[j] = -1;
                    lol_seq_dist[j] = -1;
                }
            }

            lolscore(lol_dist, lol_seq_dist, lol_score_vec, diag_length, hidden_layer);
        }
        sa_scores[sa] = maxSubArray(lol_score_vec, diag_length);
        align_startAnchors(anchor_query[sa], anchor_target[sa], maxIndexX, maxIndexY, &new_anchor_length[sa], fwbwaln->getProbabiltyMatrix(), G);
        anchor_length[sa] = new_anchor_length[sa];
    }
    
    for (int i = 0; i < queryLen; ++i) {
        std::memset(G[i], 0, targetLen * sizeof(G[0][0]));
    }
    
    index_sort(sa_scores, sa_index, num_sa);

    gaps[0] = 0;
    gaps[1] = 0;
    gaps[2] = 0;
    gaps[3] = 0;

    fwbwaln->resetParams(lol_go, lol_ge, lol_T);
    int sa;
    //bool found_nan = false;


    for (int sa_it = 0; sa_it < SeedNumber; sa_it++){
        sa = sa_index[num_sa - sa_it -1];
        bool add_seq = true;
        for(int iteration = 0; iteration < 1000; iteration++){

            gaps[0] = 0;
            gaps[1] = 0;
            gaps[2] = 0;
            gaps[3] = 0;
            /*if(iteration == 0){
                std::cout << "############  " << sa << " ## " << fwbwaln->temperature << std::endl; 
                for(int te = 0; te< queryLen; te++){
                    if(anchor_query[sa][te] != 0){
                        std::cout << te << std::endl;
                    }
                }
                for(int te = 0; te< targetLen; te++){
                    if(anchor_target[sa][te] != 0){
                        std::cout << te << std::endl;
                    }
                }
            }*/
            while((gaps[1] < queryLen && gaps[3] < targetLen)){

                calc_gap(anchor_query[sa], anchor_target[sa], gaps, queryLen, targetLen);
                

                if (gaps[0] != -1){
                    lolmatrix(anchor_query[sa], anchor_target[sa], new_anchor_length[sa], gaps, d_ij, d_kl, G, queryLen, targetLen, hidden_layer, lol_dist);

                }

                else{
                    break;
                }


            }
            for(int i = 0; i < queryLen; i++){
                if(anchor_query[sa][i] == 2){
                    anchor_query[sa][i] = 1;
                }
            }
            for(int i = 0; i <targetLen; i++){
                if(anchor_target[sa][i] == 2){
                    anchor_target[sa][i] = 1;
                }
            }
            new_anchor_length[0] = 0;

            gaps[0] = min_lolmat_idx;
            gaps[1] = 0;
            gaps[2] = 0;
            gaps[3] = 0;
            float maxP = 0.5;
            //float max_temp = 0.5;

            while((gaps[1] < max_lolmat_idx && gaps[3] < targetLen)){
                calc_gap(anchor_query[sa], anchor_target[sa], gaps, max_lolmat_idx, targetLen);
                if(gaps[0] != -1){
                    fwbwaln->initScoreMatrix(G, gaps);
                    fwbwaln->runFwBw<false, 0>(); 
                    float** fwbwP = fwbwaln->getProbabiltyMatrix();
                    maxP = std::max(maxP, fwbwaln->maxP);
                    if(fwbwaln->maxP == 0){
                        if (fwbwaln->temperature > 30){
                            break;
                        }
                       
                        fwbwaln->temperature += 2;
                        gaps[0] = 0;
                        gaps[1] = 0;
                        gaps[2] = 0;
                        gaps[3] = 0;
                        continue; 
                        
                    }
                    
                    for (int i = 0; i < gaps[1]-gaps[0]; ++i) {
                        // std::copy(fwbwaln->zm[i][0], fwbwaln->zm[i][(gaps[3] - gaps[2])], &P[i + gaps[0]][gaps[2]]);
                        std::copy(fwbwP[i], fwbwP[i] + (gaps[3] - gaps[2]), &P[i + gaps[0]][gaps[2]]);
                    }
                    
                }
                else{
                    break;
                }         
            }
            fwbwaln->temperature = lol_T;


            new_anchor_length[sa] = 0;



            gaps[0] = min_lolmat_idx;
            gaps[1] = 0;
            gaps[2] = 0;
            gaps[3] = 0;

            while((gaps[1] < max_lolmat_idx && gaps[3] < targetLen)){
                calc_gap(anchor_query[sa], anchor_target[sa], gaps, max_lolmat_idx, targetLen);
                if(gaps[0] != -1){
                    for (int i = gaps[0]; i < gaps[1]; i++){
                        for (int j = gaps[2]; j < gaps[3]; j++){
                            if(P[i][j] > maxP -0.1){
                                if(anchor_query[sa][i] == 0 && anchor_target[sa][j] == 0){
                                    anchor_query[sa][i] = 2;
                                    anchor_target[sa][j] = 2;
                                    anchor_length[sa] += 1;
                                    new_anchor_length[sa] += 1;
                                    break;
                                }
                            }
                        }
                    }
                }
                else{
                    break;
                }
            }

            if (new_anchor_length[sa] == 0){
                if(!add_seq && md == 1){
                    add_seq = true;
                    lolAlign::addForwardScoreMatrix(
                        targetNumAA,
                        targetNum3Di,
                        queryLen,
                        targetLen,
                        subMatAA,
                        G);
                }
                else{
                    break;
                }

            }
        }
        //std::cout << "maxP: " << P[queryLen-1][targetLen-1] << std::endl;
        //std::cout << "anchor_query: " << anchor_query[sa][queryLen-1] << std::endl;
        //std::cout << "anchor_target: " << anchor_target[sa][targetLen-1] << std::endl;


        for (int i = 0; i < queryLen; ++i) {
            std::memset(G[i], 0, targetLen * sizeof(G[0][0]));
        }

    }

    float max_lol_score = -10000.0;
    int max_lol_idx = 0;
    

    for (int sa_it = 0; sa_it < SeedNumber; sa_it++){
        sa = sa_index[num_sa - sa_it -1];
        
        int sa_idx = 0;
        for(int i = 0; i < queryLen; i++){
            if(anchor_query[sa][i] != 0){
                final_anchor_query[sa_idx] = i;
                sa_idx++;
               
            }
        }
        sa_idx = 0;
        for(int i = 0; i < targetLen; i++){
            if(anchor_target[sa][i] != 0){
                final_anchor_target[sa_idx] = i;
                sa_idx++;
            }
        }
        computeDi_score(targetNumAA, targetNum3Di, anchor_length[sa], final_anchor_query, final_anchor_target, subMatAA, lol_score_vec);

        for (int i = 0; i < anchor_length[sa]; i++) {
            for (int j = 0;j < anchor_length[sa]; j++) {
                if(d_ij[final_anchor_query[i]][final_anchor_query[j]] > 0.0){
                    lol_dist[j] = std::abs(d_ij[final_anchor_query[i]][final_anchor_query[j]] - d_kl[final_anchor_target[i]][final_anchor_target[j]]);
                    //anchor_dist_target[j] = d_kl[final_anchor_target[i]][final_anchor_target[j]];
                    lol_seq_dist[j] = std::copysign(1.0f, (final_anchor_query[i]-final_anchor_query[j])) * std::log(1 + std::abs((float)(final_anchor_query[i]-final_anchor_query[j])));
                }
                else{
                    lol_dist[j] = -1;
                    lol_seq_dist[j] = -1;
                }
            }
            lolscore(lol_dist, lol_seq_dist, lol_score_vec, anchor_length[sa], hidden_layer);

        }

        float total_lol_score = 0.0;
        for (int i = 0; i < anchor_length[sa]; i++) {
            total_lol_score += lol_score_vec[i];
        }
        if (total_lol_score > max_lol_score){
            max_lol_score = total_lol_score;
            max_lol_idx = sa;
        }
    }
    float seqId = 0.0;

    int sa_idx = 0;
        for(int i = 0; i < queryLen; i++){
            if(anchor_query[max_lol_idx][i] != 0){
                final_anchor_query[sa_idx] = i;
                sa_idx++;
               
            }
        }
        sa_idx = 0;
        for(int i = 0; i < targetLen; i++){
            if(anchor_target[max_lol_idx][i] != 0){
                final_anchor_target[sa_idx] = i;
                sa_idx++;
            }
        }
    
    computeDi_score(targetNumAA, targetNum3Di, anchor_length[max_lol_idx], final_anchor_query, final_anchor_target, subMatAA, lol_score_vec);
    float max_di = 0;
    for (int i = 0; i < anchor_length[max_lol_idx]; i++) {
        max_di += lol_score_vec[i];
    } 

    for(int i = 0; i < queryLen; i++){
        lol_score_vec[i] = 0;
        lol_score_vec_sh[i] = 0;
        
    }
    for (int i = 0; i < anchor_length[max_lol_idx]; i++) {
        for (int j = 0;j < anchor_length[max_lol_idx]; j++) {
            if(d_ij[final_anchor_query[i]][final_anchor_query[j]] > 0.0){
                lol_dist[j] = 0;//std::abs(d_ij[final_anchor_query[i]][final_anchor_query[j]] - d_ij[final_anchor_query[i]][final_anchor_query[j]]);
                //anchor_dist_target[j] = d_kl[final_anchor_target[i]][final_anchor_target[j]];
                lol_seq_dist[j] = std::copysign(1.0f, (final_anchor_query[i]-final_anchor_query[j])) * std::log(1 + std::abs((float)(final_anchor_query[i]-final_anchor_query[j])));
            }
            else{
                lol_dist[j] = -1;
                lol_seq_dist[j] = -1;
            }
        }
        lolscore(lol_dist, lol_seq_dist, lol_score_vec_sh, anchor_length[max_lol_idx], hidden_layer);

    }


    for (int i = 0; i < anchor_length[max_lol_idx]; i++) {
        for (int j = 0;j < anchor_length[max_lol_idx]; j++) {
            if(d_ij[final_anchor_query[i]][final_anchor_query[j]] > 0.0){
                lol_dist[j] = std::abs(d_ij[final_anchor_query[i]][final_anchor_query[j]] - d_kl[final_anchor_target[i]][final_anchor_target[j]]);
                //anchor_dist_target[j] = d_kl[final_anchor_target[i]][final_anchor_target[j]];
                lol_seq_dist[j] = std::copysign(1.0f, (final_anchor_query[i]-final_anchor_query[j])) * std::log(1 + std::abs((float)(final_anchor_query[i]-final_anchor_query[j])));
            }
            else{
                lol_dist[j] = -1;
                lol_seq_dist[j] = -1;
            }
        }
        lolscore(lol_dist, lol_seq_dist, lol_score_vec, anchor_length[max_lol_idx], hidden_layer);

    }

    float norm_lol_sh = 0;
    float temp_nan_check = 0;
    for (int i = 0; i < anchor_length[max_lol_idx]; i++) {
        temp_nan_check = lol_score_vec[i] / (lol_score_vec_sh[i]);
        norm_lol_sh +=  temp_nan_check == temp_nan_check ? temp_nan_check : 0;
    }




    std::string backtrace = "";
    int matches = 0;
    int q_count = 0;
    int t_count = 0;
    



    //build_dataset.append(to_string(total_lol_score));
    //build_dataset.append(to_string(max_di));
    //build_dataset.append(to_string(norm_lol_sh));


    

    /*for(int i = 0 ; i<queryLen; i++){
        std::cout << anchor_query[max_lol_idx][i];
    }
    std::cout << std::endl;
    for(int i = 0 ; i<targetLen; i++){
        std::cout << anchor_target[max_lol_idx][i];
    }
    std::cout << std::endl;
    std::cout << "queryLen: " << queryLen << std::endl;*/
    
    



    while(matches < anchor_length[max_lol_idx]){
        if(anchor_query[max_lol_idx][q_count] != 0 && anchor_target[max_lol_idx][t_count] != 0){
            backtrace.append(1, 'M');
            matches++;
            if(querySeq[q_count] == targetSeq[t_count]){
                seqId += 1.0;
            }
            q_count++;
            t_count++;
        }
        //else if(anchor_query[max_lol_idx][q_count] != 0 && anchor_target[max_lol_idx][t_count] == 0){
        else if(anchor_target[max_lol_idx][t_count] == 0){
            backtrace.append(1, 'D');
            t_count++;
        }
        //else if(anchor_query[max_lol_idx][q_count] == 0 && anchor_target[max_lol_idx][t_count] != 0){
        else if(anchor_query[max_lol_idx][q_count] == 0){
            backtrace.append(1, 'I');
            q_count++;
        }
        //else if(anchor_query[max_lol_idx][q_count] == 0 && anchor_target[max_lol_idx][t_count] == 0){
        //    backtrace.append(1, 'D');
        //    backtrace.append(1, 'I');
        //    q_count++;
        //    t_count++;
        //}
    }
    //std::cout << "backtrace: " << backtrace << " " << backtrace.size() << " " << queryLen << std::endl;



    Matcher::result_t result = Matcher::result_t();
    result.seqId = seqId / (float)anchor_length[max_lol_idx];
    result.qcov = anchor_length[max_lol_idx] / (float)queryLen;
    result.dbcov = anchor_length[max_lol_idx] / (float)targetLen;
    result.score = max_lol_score; 
    if(md == 0){
        result.eval = (((max_lol_score +3 * max_di) * norm_lol_sh/anchor_length[max_lol_idx])/QQ_score)/ std::pow(queryLen * targetLen, 0.25);
        result.score = (((max_lol_score +3 * max_di) * norm_lol_sh/anchor_length[max_lol_idx]))/ std::pow(queryLen * targetLen, 0.25);
    }
    else{
        result.eval = (((max_lol_score +3 * max_di) * norm_lol_sh/anchor_length[max_lol_idx])/QQ_score);
        result.score = (((max_lol_score +3 * max_di) * norm_lol_sh/anchor_length[max_lol_idx]));
    }
    result.dbKey = dbKey;
    result.qStartPos = 0;
    result.dbStartPos = 0;
    result.qEndPos = queryLen - 1;
    result.dbEndPos = targetLen - 1;
    result.qLen = queryLen;
    result.dbLen = targetLen;
    result.alnLength = backtrace.length();
    result.backtrace = Matcher::compressAlignment(backtrace);

    

    // trim backtrace find the first 'M'
    size_t firstM = targetLen;
    int qLeftSkip = 0;
    int tLeftSkip = 0;
    for (size_t i = 0; i < backtrace.size(); i++) {
        if (backtrace[i] == 'M') {
            firstM = i;
            break;
        }
        if (backtrace[i] == 'I') {
            qLeftSkip++;
        } else if (backtrace[i] == 'D') {
            tLeftSkip++;
        }
    }
    /*if(firstM == targetLen){
        Debug(Debug::ERROR) << "No 'M' in backtrace.\n";
        std::cout << backtrace << std::endl;
        EXIT(EXIT_FAILURE);
    }*/

    result.qStartPos  += qLeftSkip;
    result.dbStartPos += tLeftSkip;
    result.qEndPos = result.qStartPos;
    result.dbEndPos = result.dbStartPos;
    for (size_t i = firstM; i < backtrace.size(); i++) {
        switch (backtrace[i]) {
            case 'M': // match consumes 1 base in query + 1 base in target
                result.qEndPos++;
                result.dbEndPos++;
                break;
            case 'I': // insertion consumes 1 base in query
                result.qEndPos++;
                break;
            case 'D': // deletion consumes 1 base in target
                result.dbEndPos++;
                break;
            default:
                // handle unexpected cigar chars if needed
                break;
        }
    }


    result.qEndPos--;
    result.dbEndPos--;
    result.backtrace = Matcher::compressAlignment(backtrace.substr(firstM));
    //result.backtrace = backtrace.substr(firstM);
    result.alnLength = (int)result.backtrace.size();

    return result;

}
void lolAlign::align_startAnchors(int *anchor_query, int *anchor_target, int max_query, int max_target, int *anchor_length, float** fwbwP, float** G){

    for (int i = max_query - this->start_anchor_length; i <= max_query + this->start_anchor_length; ++i)
    {
        int max_ind = max_target + i - max_query;

        anchor_query[i] = 2;
        anchor_target[max_ind] = 2;
        fwbwP[i][max_ind] = 0;
        G[i][max_ind] = 0;
        *anchor_length = *anchor_length + 1;

    }
}

void lolAlign::calc_dist_matrix(float* x, float* y, float* z, size_t len, float** d, bool cutoff){
    const float cutoff_distance = 20.0f;
    const float cutoff_sq = cutoff_distance * cutoff_distance;

    for (size_t i = 0; i < len; ++i){
        d[i][i] = 0.0f;
        for (size_t j = i + 1; j < len; ++j){
            float dx = x[i] - x[j];
            float dy = y[i] - y[j];
            float dz = z[i] - z[j];

            float dist_sq = dx * dx + dy * dy + dz * dz;

            if (cutoff && dist_sq > cutoff_sq){
                d[i][j] = 0.0f;
            }
            else{
                float dist = std::sqrt(dist_sq);
                d[i][j] = dist;
            }
            d[j][i] = d[i][j]; 
        }
    }
}

void lolAlign::computeDi_score(
        unsigned char * targetNumAA,
        unsigned char *targetNum3Di,
        int anchorLen,
        int* final_anchor_query,
        int* final_anchor_target,
        SubstitutionMatrix &subMatAA,
        float *scoreForward)
{


    for (int i = 0; i < anchorLen; ++i)
    {
        int q = final_anchor_query[i];
        int t = final_anchor_target[i];
        scoreForward[i] = static_cast<float>((subMatAA.subMatrix[queryNumAA[q]][targetNumAA[t]] *1.4 ) + (scoringMatrix3Di[queryNum3Di[q]][targetNum3Di[t]] * 2.1));        
    }
}


void lolAlign::initQuery(float *x, float *y, float *z, Sequence& qSeqAA, Sequence& qSeq3Di, int queryLen, SubstitutionMatrix &subMatAA, int maxTLen, int md)
{
    memcpy(query_x, x, sizeof(float) * queryLen);
    memcpy(query_y, y, sizeof(float) * queryLen);
    memcpy(query_z, z, sizeof(float) * queryLen);

    this->queryLen = queryLen;
    this->querySeq = qSeqAA.getSeqData();
    this->query3diSeq = qSeq3Di.getSeqData();
    queryNumAA = qSeqAA.numSequence;
    queryNum3Di = qSeq3Di.numSequence;
    Coordinates queryCaCords;
    queryCaCords.x = query_x;
    queryCaCords.y = query_y;
    queryCaCords.z = query_z;
    //anchor_query = malloc_matrix<int>(num_sa, queryLen);
    G = malloc_matrix<float>(queryLen, maxTLen);
    P = malloc_matrix<float>(queryLen, maxTLen);
    hidden_layer = malloc_matrix<float>(maxTLen, 3);
    anchor_query = malloc_matrix<int>(num_sa, queryLen);
    anchor_target = malloc_matrix<int>(num_sa,maxTLen);
    lol_dist = (float *)mem_align(ALIGN_FLOAT, maxTLen * sizeof(float));
    lol_seq_dist = (float *)mem_align(ALIGN_FLOAT, maxTLen * sizeof(float));
    lol_score_vec = (float *)mem_align(ALIGN_FLOAT, maxTLen * sizeof(float));
    lol_score_vec_sh = (float *)mem_align(ALIGN_FLOAT, maxTLen * sizeof(float));
    final_anchor_query = new int[maxTLen];
    final_anchor_target = new int[maxTLen];
    calc_dist_matrix(query_x, query_y, query_z, queryLen, d_ij, true);
    for(int i = 0; i < queryLen; i++){
        final_anchor_query[i] = i;
    }


    computeDi_score(queryNumAA, queryNum3Di, queryLen, final_anchor_query, final_anchor_query, subMatAA, lol_score_vec);
    float di_score = 0;
    QQ_score = 0;
    for(int i = 0; i < queryLen; i++){
        di_score += lol_score_vec[i];
    }




    for(int i = 0; i < queryLen; i++){
        for(int j=0; j < queryLen; j++){
            if(d_ij[i][j] > 0){

                lol_dist[j] = 0;
                lol_seq_dist[j] =  std::copysign(1.0f, (i-j)) * std::log(1 + std::abs((float)(i-j)));

            }
            else{
                lol_dist[j] = -1;
                lol_seq_dist[j] = -1;
            }
        }
        lolscore(lol_dist, lol_seq_dist, lol_score_vec, queryLen, hidden_layer);
    }



    for(int i = 0; i < queryLen; i++){
        QQ_score += lol_score_vec[i];
    }
    QQ_score = (QQ_score + 3*di_score);
    if(md == 0){
        QQ_score = QQ_score/std::pow(queryLen * queryLen, 0.25);
    }
    
   




    return;
}

void lolAlign::lolmatrix(int *anchor_query, int *anchor_target, int anchor_length, int *gaps, float **d_ij, float **d_kl, float **G, int queryLen, int targetLen, float ** hidden_layer, float * d_dist)
{
    int gap0_start = gaps[0];
    int gap0_end = gaps[1];
    int gap1_start = gaps[2];
    int gap1_end = gaps[3];
    //float *seq_dist = new float[gap1_end - gap1_start];
    float seq_dist = 0.0;
    int anchor_q = 0;
    int anchor_t = 0;

    for (int i = 0; i < anchor_length; i++){
        for(int aq = anchor_q+1; aq < queryLen; aq++){
            if(anchor_query[aq] == 2){
                //anchor_query[aq] = 1;
                anchor_q = aq;
                //std::cout << "anchor_q: " << anchor_q << std::endl;
                break;
            }
        }
        for(int at = anchor_t+1; at < targetLen; at++){
            if(anchor_target[at] == 2){
                //anchor_target[at] = 1;
                anchor_t = at;
                //std::cout << "anchor_t: " << anchor_t << std::endl;
                break;
            }
        }
        

        for (int j = gap0_start; j < gap0_end; j++)
        {
            float dq = d_ij[anchor_q][j];

            if (dq > 0)
            {
                min_lolmat_idx = std::min(min_lolmat_idx, j);
                max_lolmat_idx = std::max(max_lolmat_idx, j+1); 
                seq_dist = std::copysign(1.0f, (anchor_q-j)) * std::log(1 + std::abs((float)(anchor_q-j)));
                for(int l = gap1_start; l < gap1_end; l++)
                {
                    d_dist[l - gap1_start] = std::abs(dq- d_kl[anchor_t][l]);

                    //seq_dist[l - gap1_start] = std::copysign(1.0f, (anchor_q-j)) * std::log(1 + std::abs((float)(anchor_q-j)));
                }
                lolscore(d_dist, seq_dist, G[j], gap1_end - gap1_start, gap1_start, hidden_layer);
            }
        }
    }
    return;
}


void lolAlign::lolscore(float* d_dist, float d_seq, float* score, int length, int start, float** hidden_layer)
{
     // Zero vector for ReLU

    // Process 8 elements at a time
    simd_float seq0 = simdf32_mul(simdf32_set(d_seq), w1_0);
    simd_float seq1 = simdf32_mul(simdf32_set(d_seq), w1_1);
    simd_float seq2 = simdf32_mul(simdf32_set(d_seq), w1_2);
    int i = 0;
    for (; i <= length - VECSIZE_FLOAT; i += VECSIZE_FLOAT) {
        // Load d_dist[i..i+7] into a SIMD register
        simd_float d_dist_vec = simdf32_loadu(&d_dist[i]);

        // Compute hidden_layer[i][k] for k = 0, 1, 2
        simd_float hl_0 = simdf32_add(seq0, simdf32_fmadd(d_dist_vec, w1_d0, b1_0));
        simd_float hl_1 = simdf32_add(seq1, simdf32_fmadd(d_dist_vec, w1_d1, b1_1));
        simd_float hl_2 = simdf32_add(seq2, simdf32_fmadd(d_dist_vec, w1_d2, b1_2));

        // Apply ReLU (max(0, x))
        hl_0 = simdf32_max(hl_0, zero);
        hl_1 = simdf32_max(hl_1, zero);
        hl_2 = simdf32_max(hl_2, zero);

        // Store hidden_layer[i][k] back to memory
        simdf32_storeu(&hidden_layer[i][0], hl_0);
        simdf32_storeu(&hidden_layer[i][1], hl_1);
        simdf32_storeu(&hidden_layer[i][2], hl_2);

        // Compute score[i+start] += hidden_layer[i][k] * w2[k] for k = 0, 1, 2
        simd_float score_vec = simdf32_loadu(&score[i + start]);
        score_vec = simdf32_fmadd(hl_0, w2_0, score_vec);
        score_vec = simdf32_fmadd(hl_1, w2_1, score_vec);
        score_vec = simdf32_fmadd(hl_2, w2_2, score_vec);
        score_vec = simdf32_add(score_vec, b2_vec);

        // Store the updated score back to memory
        simdf32_storeu(&score[i + start], score_vec);
    }

    // Process remaining elements (if length is not a multiple of 8)
    for (; i < length; ++i) {
        for (int k = 0; k < 3; ++k) {
            hidden_layer[i][k] = d_seq * w1[0][k];
            hidden_layer[i][k] += d_dist[i] * w1[1][k];
            hidden_layer[i][k] += b1[k];
            hidden_layer[i][k] = std::max(0.0f, hidden_layer[i][k]);
            score[i + start] += hidden_layer[i][k] * w2[k];
        }
        score[i + start] += b2;
    }
}



void lolAlign::lolscore(float* dist, float* d_seq, float* score, int length, float** hidden_layer)
{

    for (int i = 0; i < length; ++i) {
        if(dist[i] >= 0){
            for (int k = 0; k < 3; ++k) {
                //hidden_layer[i][k] = 0.0;
                hidden_layer[i][k] = d_seq[i] * w1[0][k];
                hidden_layer[i][k] += dist[i] * w1[1][k];
                hidden_layer[i][k] += b1[k];
                hidden_layer[i][k] = std::max(0.0f, hidden_layer[i][k]); // ReLU activation
                score[i] += hidden_layer[i][k] * w2[k];
            }
            score[i] += b2;
        }
    }
}


void lolAlign::computeForwardScoreMatrix(
        unsigned char * targetNumAA,
        unsigned char *targetNum3Di,
        int queryLen,
        int targetLen,
        SubstitutionMatrix &subMatAA,
        float **scoreForward)
{
    for (int i = 0; i < queryLen; ++i)
    {
        for (int j = 0; j < targetLen; ++j)
        {
            scoreForward[i][j] = static_cast<float>((subMatAA.subMatrix[queryNumAA[i]][targetNumAA[j]] *1.4 ) + (scoringMatrix3Di[queryNum3Di[i]][targetNum3Di[j]] * 2.1));
            //scoreForward[i][j] = exp(scoreForward[i][j] / T);
        }
    }
}

void lolAlign::addForwardScoreMatrix(
    unsigned char * targetNumAA,
    unsigned char *targetNum3Di,
    int queryLen,
    int targetLen,
    SubstitutionMatrix &subMatAA,
    float **scoreForward)
{
    for (int i = 0; i < queryLen; ++i)
    {
        for (int j = 0; j < targetLen; ++j)
        {
            scoreForward[i][j] += static_cast<float>((subMatAA.subMatrix[queryNumAA[i]][targetNumAA[j]] *1.4 ) + (scoringMatrix3Di[queryNum3Di[i]][targetNum3Di[j]] * 2.1));
            //scoreForward[i][j] = exp(scoreForward[i][j] / T);
        }
    }
}



int lolalign(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);
    //subMat3Di.subMatrix[2][2] = 0.5;

    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == "blosum62.out")  {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
            std::string matrixName = par.substitutionMatrices[i].name;
            char *serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
            blosum.assign(serializedMatrix);
            free(serializedMatrix);
            break;
        }
    }

    float aaFactor = (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA) ? 1.4 : 0.0;
    SubstitutionMatrix subMatAA(blosum.c_str(), aaFactor, par.scoreBias);

    Debug(Debug::INFO) << "Query database: " << par.db1 << "\n";
    Debug(Debug::INFO) << "Target database: " << par.db2 << "\n";
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader qdbr(par.db1, par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader qcadbr(
        par.db1,
        par.threads,
        IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
        touch ? IndexReader::PRELOAD_INDEX : 0,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
        "_ca"
    );
    IndexReader qdbr3Di(
        par.db1,
        par.threads,
        IndexReader::SEQUENCES,
        touch ? IndexReader::PRELOAD_INDEX : 0,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
        "_ss"
    );

    IndexReader *tdbr = NULL;
    IndexReader *tcadbr = NULL;
    IndexReader *tdbr3Di = NULL;

    bool sameDB = false;
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool alignmentIsExtended = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
        tcadbr = &qcadbr;
        tdbr3Di = &qdbr3Di;
    } else {
        tdbr = new IndexReader(
            par.db2, par.threads,
            alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
            (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0
        );

        tcadbr = new IndexReader(
            par.db2,
            par.threads,
            alignmentIsExtended ? IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB2) : IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
            touch ? IndexReader::PRELOAD_INDEX : 0,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
            alignmentIsExtended ? "_seq_ca" : "_ca"
        );
        
        std::string t3DiDbrName = StructureUtil::getIndexWithSuffix(par.db2, "_ss");
        bool is3DiIdx = Parameters::isEqualDbtype(FileUtil::parseDbType(t3DiDbrName.c_str()), Parameters::DBTYPE_INDEX_DB);
        tdbr3Di = new IndexReader(
            is3DiIdx ? t3DiDbrName : par.db2, par.threads,
            alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
            (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
            alignmentIsExtended ? "_seq_ss" : "_ss"
        );
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbtype = Parameters::DBTYPE_ALIGNMENT_RES;
    if (alignmentIsExtended) {
        dbtype = DBReader<unsigned int>::setExtendedDbtype(dbtype, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbtype);
    dbw.open();

    //fwbw aligner
    size_t blockLen_fwbw = 16;
    float gapOpen_fwbw = -2.0;
    float gapExtend_fwbw = -2.0;
    float temperature_fwbw = 2.0;
    int max_targetLen = static_cast<int>(tdbr->sequenceReader->getMaxSeqLen() + 1);
    max_targetLen = ((max_targetLen + blockLen_fwbw -1) / blockLen_fwbw) * blockLen_fwbw;

    Debug::Progress progress(resultReader.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        
        std::vector<Matcher::result_t> swResults;
        swResults.reserve(300);
        std::string backtrace;
        std::string resultBuffer;
        resultBuffer.reserve(1024 * 1024);
        Coordinate16 qcoords;
        Coordinate16 tcoords;

        Sequence qSeqAA(par.maxSeqLen, qdbr.sequenceReader->getDbtype(), (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence qSeq3Di(par.maxSeqLen, qdbr3Di.sequenceReader->getDbtype(), (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        Sequence tSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence tSeq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        FwBwAligner fwbwaln(gapOpen_fwbw, gapExtend_fwbw, temperature_fwbw, 0, 1, max_targetLen, blockLen_fwbw, 0);
        char buffer[1024 + 32768];

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            size_t queryKey = resultReader.getDbKey(id);
            char *data = resultReader.getData(id, thread_idx);
            if (*data != '\0') {
                unsigned int queryId = qdbr.sequenceReader->getId(queryKey);

                char *querySeq = qdbr.sequenceReader->getData(queryId, thread_idx);
                char *query3diSeq = qdbr3Di.sequenceReader->getData(queryId, thread_idx);
                int queryLen = static_cast<int>(qdbr.sequenceReader->getSeqLen(queryId));
                qSeqAA.mapSequence(id, queryKey, querySeq, queryLen);
                qSeq3Di.mapSequence(id, queryKey, query3diSeq, queryLen);

                char *qcadata = qcadbr.sequenceReader->getData(queryId, thread_idx);
                size_t qCaLength = qcadbr.sequenceReader->getEntryLen(queryId);
                float *qdata = qcoords.read(qcadata, queryLen, qCaLength);

                lolAlign lolaln(std::max(qdbr.sequenceReader->getMaxSeqLen() + 1, tdbr->sequenceReader->getMaxSeqLen() + 1), false);
                lolaln.initQuery(qdata, &qdata[queryLen], &qdata[queryLen + queryLen], qSeqAA, qSeq3Di, queryLen, subMatAA, max_targetLen, par.multiDomain);
                fwbwaln.resizeMatrix<false, 0>(queryLen, max_targetLen);

                if (queryLen <= 10) {
                    lolaln.set_start_anchor_length(0);
                }

                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int)strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = tdbr->sequenceReader->getId(dbKey);

                    char *targetSeq = tdbr->sequenceReader->getData(targetId, thread_idx);
                    int targetLen = static_cast<int>(tdbr->sequenceReader->getSeqLen(targetId));
                    tSeqAA.mapSequence(targetId, dbKey, targetSeq, targetLen);

                    char *target3diSeq = tdbr3Di->sequenceReader->getData(targetId, thread_idx);
                    tSeq3Di.mapSequence(targetId, dbKey, target3diSeq, targetLen);
                    // if (targetLen > max_targetLen) {
                    //     max_targetLen = ((targetLen*2)/16)*16;
                    //     lolaln.reallocate_target(max_targetLen);
                    //     fwbwaln.resizeMatrix(queryLen, max_targetLen);
                    //     std::cout << "Reallocating target to: " << max_targetLen << std::endl;
                    // }
                    // if (Util::canBeCovered(par.covThr, par.covMode, queryLen, targetLen) == false) {
                    //     continue;
                    // }

                    char *tcadata = tcadbr->sequenceReader->getData(targetId, thread_idx);
                    size_t tCaLength = tcadbr->sequenceReader->getEntryLen(targetId);
                    float *tdata = tcoords.read(tcadata, targetLen, tCaLength);

                    Matcher::result_t result;
                    if (targetLen <= 10) {
                        lolaln.set_start_anchor_length(1);
                        if (targetLen < 4) {
                            lolaln.set_start_anchor_length(0);
                        }
                        result = lolaln.align(dbKey, tdata, &tdata[targetLen], &tdata[targetLen + targetLen], tSeqAA, tSeq3Di, targetLen, subMatAA, &fwbwaln, par.multiDomain);
                        if (queryLen > 10) {
                            lolaln.set_start_anchor_length(3);
                        }
                    } else {
                        result = lolaln.align(dbKey, tdata, &tdata[targetLen], &tdata[targetLen + targetLen], tSeqAA, tSeq3Di, targetLen, subMatAA, &fwbwaln, par.multiDomain);
                    }

                    bool hasCov = Util::hasCoverage(par.covThr, par.covMode, 1.0, 1.0);
                    bool hasSeqId = result.seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    //bool hasTMscore = (TMscore >= par.tmScoreThr);

                    if (hasCov && hasSeqId) {
                        swResults.emplace_back(result);
                    }
                }
                
                SORT_SERIAL(swResults.begin(), swResults.end(), compareHitsBylolScore);
                for (size_t i = 0; i < swResults.size(); i++) {
                    size_t len = Matcher::resultToBuffer(buffer, swResults[i], par.addBacktrace, false);
                    resultBuffer.append(buffer, len);
                }
                swResults.clear();
            }

            dbw.writeData(resultBuffer.c_str(), resultBuffer.size(), queryKey, thread_idx);
            resultBuffer.clear();
        }
    }
    dbw.close();
    resultReader.close();
    if (sameDB == false) {
        delete tdbr;
        delete tcadbr;
        delete tdbr3Di;
    }
    return EXIT_SUCCESS;
}
