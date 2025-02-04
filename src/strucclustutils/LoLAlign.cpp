#ifndef LoLAlign
#define LoLAlign
#include <iostream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <numeric>
#include <cmath>
#include <vector>
#include <fstream>
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
#include "simd.h"



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
    free(sa_index);
    free(sa_scores);
    free(anchor_length);
    free(new_anchor_length);
    free(gaps);
    free(lol_dist);
    free(lol_seq_dist);
    free(lol_score_vec);
    free(final_anchor_query);
    free(final_anchor_target);


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
    lol_dist = new float[targetL];
    free(lol_seq_dist);
    lol_seq_dist = new float[targetL];
    free(lol_score_vec);
    lol_score_vec = new float[targetL];  
    free(final_anchor_target);
    final_anchor_target = new int[targetL];
    free(final_anchor_query);
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
    //std::iota(index, index + numsSize, 0);
    std::sort(index, index + numsSize, [&nums](int i1, int i2) { return nums[i1] < nums[i2]; });
}




Matcher::result_t lolAlign::align(unsigned int dbKey, float *target_x, float *target_y, float *target_z,
                                  char * targetSeq, char* target3diSeq, unsigned int targetLen, SubstitutionMatrix &subMatAA, SubstitutionMatrix &subMat3Di, FwBwAligner* fwbwaln)
{


    lolAlign::computeForwardScoreMatrix(
            querySeq,
            query3diSeq,
            targetSeq,
            target3diSeq,
            queryLen,
            targetLen,
            subMatAA,
            subMat3Di,
            start_anchor_T,
            G);



    int maxIndexX = 0;
    int maxIndexY = 0;
    
    for(int i = 0; i < num_sa; i++){
        sa_index[i] = i;
        sa_scores[i] = 0;
    }



    calc_dist_matrix(target_x, target_y, target_z, targetLen, d_kl, false);

    for(int sa = 0; sa < 10; sa++){
        anchor_length[sa] = 0;
        new_anchor_length[sa] = 0;
        for (unsigned int i = 0; i < queryLen; i++) {
            anchor_query[sa][i] = 0;
        }
        for (unsigned int i = 0; i < targetLen; i++) {
            anchor_target[sa][i] = 0;
        }
    }

    gaps[0] = 0;
    gaps[1] = queryLen;
    gaps[2] = 0;
    gaps[3] = targetLen;


    
    fwbwaln->setParams(start_anchor_go, start_anchor_ge, start_anchor_T, 16);
        




    for(int sa = 0; sa < 10; sa++){

        if(sa % 5 == 0){
            fwbwaln->initScoreMatrix(G, targetLen, queryLen, gaps);
            fwbwaln->computeProbabilityMatrix(false);
            //float** fwbwaln_zm = fwbwaln.getZm();
            for (size_t i = 0; i < queryLen; ++i)
            {
                for (size_t j = 0; j < targetLen; ++j)
                {
                    P[i][j] = fwbwaln->zm[i][j];
                }
            }
            //lolAlign::lol_fwbw(G, P, queryLen, targetLen, assignTargetLen, start_anchor_go, start_anchor_ge, start_anchor_T, length, blocks, gaps_start);
        }
        
        


        maxIndexX = 0;
        maxIndexY = 0;

        float maxScore = 0;

        for (int i = this->start_anchor_length; i < static_cast<int>(queryLen) - this->start_anchor_length ; ++i) {
            for (int j = this->start_anchor_length; j < static_cast<int>(targetLen) - this->start_anchor_length ; ++j) {
                if (P[i][j] > maxScore) {
                    maxScore = P[i][j];
                    maxIndexX = i;
                    maxIndexY = j;
                }
            }
        }



        int start_row = maxIndexX - std::min(maxIndexX, maxIndexY);
        int start_col = maxIndexY - std::min(maxIndexX, maxIndexY);
        int diag_length = std::min(queryLen - start_row, targetLen - start_col);
        for(int i = 0; i < diag_length; i++){
            lol_score_vec[i] = G[start_row + i][start_col + i] *start_anchor_T;
        }
        lol_score_vec[std::min(maxIndexX, maxIndexY)] += 100;
         




        for(int i = -start_anchor_length; i < start_anchor_length; i++){
            for(int j=0; j<diag_length; j++){

                if(d_ij[maxIndexX+i][start_row + j] > 0){

                    lol_dist[j] = std::abs(d_ij[maxIndexX+i][start_row + j] - d_kl[maxIndexY+i][start_col+j]);
                    lol_seq_dist[j] =  std::copysign(1.0f, (maxIndexX+i - start_row + j)) * std::log(1 + std::abs((float)(maxIndexX+i - start_row + j)));

                    //std::cout << "diag_row_dist[" << j<< "]: " << diag_row_dist[j] - diag_col_dist[j] << std::endl;
                }
                else{
                    lol_dist[j] = -1;
                    //diag_col_dist[j] = 0;
                    lol_seq_dist[j] = 0;
                }
            }

            lolscore(lol_dist, lol_seq_dist, lol_score_vec, diag_length, hidden_layer);
        }
        sa_scores[sa] = maxSubArray(lol_score_vec, diag_length);
        align_startAnchors(anchor_query[sa], anchor_target[sa], maxIndexX, maxIndexY, &new_anchor_length[sa], P, G);
        anchor_length[sa] = new_anchor_length[sa];


    }
    

    for (size_t i = 0; i < queryLen; ++i)
    {
        for (size_t j = 0; j < targetLen; ++j)
        {
            G[i][j] = 0;
            P[i][j] = 0;
        }
    }
    



    index_sort(sa_scores, sa_index, 10);



    int* gaps = new int[4]{0, 0, 0, 0};
    fwbwaln->setParams(lol_go, lol_ge, lol_T, 16);
    int sa;
    bool found_nan = false;





    for (int sa_it = 0; sa_it < 3; sa_it++){
        sa = sa_index[9 - sa_it];
        //std::cout << "sa: " << sa << std::endl;
        for(int iteration = 0; iteration < 1000; iteration++){
            //std::cout << "iteration: " << iteration << std::endl;

            gaps[0] = 0;
            gaps[1] = 0;
            gaps[2] = 0;
            gaps[3] = 0;
            
            
    
            while((gaps[1] < queryLen && gaps[3] < targetLen)){

                calc_gap(anchor_query[sa], anchor_target[sa], gaps, queryLen, targetLen);

                if (gaps[0] != -1){
                    lolmatrix(anchor_query[sa], anchor_target[sa], new_anchor_length[sa], gaps, d_ij, d_kl, G, queryLen, targetLen, hidden_layer, lol_dist, lol_score_vec);

                }

                else{
                    break;
                }
                //std::cout << "gaps[0]: " << gaps[0] << " gaps[1]: " << gaps[1] << " gaps[2]: " << gaps[2] << " gaps[3]: " << gaps[3] << std::endl;
                //std::cout << "queryLen: " << queryLen << " targetLen: " << targetLen << std::endl;

            }
            for(unsigned int i = 0; i <queryLen; i++){
                if(anchor_query[sa][i] == 2){
                    anchor_query[sa][i] = 1;
                }
            }
            for(unsigned int i = 0; i <targetLen; i++){
                if(anchor_target[sa][i] == 2){
                    anchor_target[sa][i] = 1;
                }
            }
            new_anchor_length[0] = 0;

            gaps[0] = 0;
            gaps[1] = 0;
            gaps[2] = 0;
            gaps[3] = 0;
            float maxP = 0.5;

            while((gaps[1] < queryLen && gaps[3] < targetLen)){
                calc_gap(anchor_query[sa], anchor_target[sa], gaps, queryLen, targetLen);
                if(gaps[0] != -1){
                    fwbwaln->initScoreMatrix(G, gaps[3]-gaps[2], gaps[1]-gaps[0], gaps);
                    fwbwaln->computeProbabilityMatrix(false);
                    //float** fwbwaln_zm = fwbwaln.getZm();
                    for (size_t i = gaps[0]; i < gaps[1]; ++i)
                    {
                        for (size_t j = gaps[2]; j < gaps[3]; ++j)
                        {
                            P[i][j] = fwbwaln->zm[i - gaps[0]][j - gaps[2]];
                        }
                    }



                    //for (size_t i = 0; i < gaps[1] -gaps[0]; ++i) {
                        // Copy entire row segment in one operation
                        //for (size_t j = 0; j < (gaps[3] - gaps[2]) - VECSIZE_FLOAT;  j += VECSIZE_FLOAT)
                        //{
                        //    simd_float vScoreForward = simdf32_loadu(&fwbwaln->zm[i][j]);
                        //    simdf32_store(&P[i+gaps[0]][j+gaps[2]], vScoreForward);
                            //P[i][j] = fwbwaln->zm[i - gaps[0]][j - gaps[2]];
                        //}
                        //for(size_t j = std::max((gaps[3] - gaps[2]) - VECSIZE_FLOAT, 0); j < gaps[3]; j++){
                        //    P[i+gaps[0]][j+gaps[2]] = fwbwaln->zm[i][j];
                        //}
                        //std::copy(&fwbwaln->zm[i][0], &fwbwaln->zm[i][(gaps[3] - gaps[2])], &P[i + gaps[0]][gaps[2]]);
                    //}

                    //lolAlign::lol_fwbw(G, P, queryLen, targetLen, assignTargetLen, start_anchor_go, start_anchor_ge, 2, length, blocks, gaps);

                }
                //std::cout << "gaps[0]: " << gaps[0] << " gaps[1]: " << gaps[1] << " gaps[2]: " << gaps[2] << " gaps[3]: " << gaps[3] << std::endl;
                
                

                if (gaps[0] != -1){
                    for (int i = gaps[0]; i < gaps[1]; i++){
                        for (int j = gaps[2]; j < gaps[3]; j++){
                            if (std::isnan(P[i][j])) {
                                fwbwaln->temperature += 1;
                                gaps[0] = 0;
                                gaps[1] = 0;
                                gaps[2] = 0;
                                gaps[3] = 0;
                                found_nan = true;
                                break;
                            }
                            if (P[i][j] > maxP){
                                maxIndexX = i;
                                maxIndexY = j;
                                maxP = P[i][j];
                            }
                        }
                        if(found_nan){
                            break;
                        }
                    }
                    
                }
                else{
                    break;
                }
                if(found_nan){
                    found_nan = false;
                    break;
                }
                fwbwaln->temperature = lol_T;
                

            }

            new_anchor_length[sa] = 0;



            gaps[0] = 0;
            gaps[1] = 0;
            gaps[2] = 0;
            gaps[3] = 0;

            while((gaps[1] < queryLen && gaps[3] < targetLen)){
                calc_gap(anchor_query[sa], anchor_target[sa], gaps, queryLen, targetLen);
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
                                    //std::cout << "anchor_query[" << i << "]" << "anchor_target[" << j << "] " << P[i][j] << " " << maxP << std::endl;
                                }
                            }
                            //P[i][j] = 0;
                        }
                    }
                }
                else{
                    break;
                }
            }
            if (new_anchor_length[sa] == 0){
                break;

            }
        }

        /*if(true){
        
            std::ofstream outfile("/home/lasse/Desktop/Projects/FB_martin/zmForward.txt");
            if (outfile.is_open())
            {
                for (size_t i = 0; i < queryLen; ++i)
                {
                    for (size_t j = 0; j < targetLen; ++j)
                    {
                        outfile << (G[i][j]) << " ";
                    }
                    outfile << "\n";
                }
                outfile.close();
            }
            else
            {
                std::cerr << "Unable to open file for writing P matrix." << std::endl;
            }

            std::ofstream outfile2("/home/lasse/Desktop/Projects/FB_martin/P_mat.txt");
            if (outfile2.is_open())
            {
                for (size_t i = 0; i < queryLen; ++i)
                {

                    for (size_t j = 0; j < targetLen; ++j){
                        outfile2 <<  P[i][j] << " ";
                        //if(i==j && anchor_query[i] != 0 && anchor_target[j] != 0){
                        //    outfile2 << 1  << " ";
                       // }
                        //else
                        //{
                        //    outfile2 << 0 << " ";
                       // }
                   }
                    outfile2 << "\n";
                }

                outfile2.close();
            }
            else
            {
                std::cerr << "Unable to open file for writing P matrix." << std::endl;
            }
        }*/
    for (size_t i = 0; i < queryLen; ++i)
        {
            for (size_t j = 0; j < targetLen ; ++j)
            {
                G[i][j] = 0;
                P[i][j] = 0;
            }
        }

    }

    float max_lol_score = 0.0;
    int max_lol_idx = 0;
    

    for (int sa_it = 0; sa_it < 3; sa_it++){
        sa = sa_index[9 - sa_it];
        
        int sa_idx = 0;
        for(unsigned int i = 0; i < queryLen; i++){
            if(anchor_query[sa][i] != 0){
                final_anchor_query[sa_idx] = i;
                sa_idx++;
               
            }
        }
        sa_idx = 0;
        for(unsigned int i = 0; i < targetLen; i++){
            if(anchor_target[sa][i] != 0){
                final_anchor_target[sa_idx] = i;
                sa_idx++;
            }
        }

        computeDi_score(querySeq, query3diSeq, targetSeq, target3diSeq, anchor_length[sa], final_anchor_query, final_anchor_target, subMatAA, subMat3Di, lol_score_vec);


        for (int i = 0; i < anchor_length[sa]; i++) {
            lol_score_vec[i] = G[final_anchor_query[i]][final_anchor_target[i]];
        }

        for (int i = 0; i < anchor_length[sa]; i++) {
            for (int j = 0;j < anchor_length[sa]; j++) {
                if(d_ij[final_anchor_query[i]][final_anchor_query[j]] > 0.0){
                    lol_dist[j] = std::abs(d_ij[final_anchor_query[i]][final_anchor_query[j]] - d_kl[final_anchor_target[i]][final_anchor_target[j]]);
                    lol_seq_dist[j] = std::copysign(1.0f, (final_anchor_query[i]-final_anchor_query[j])) * std::log(1 + std::abs((float)(final_anchor_query[i]-final_anchor_query[j])));
                }
                else{
                    lol_dist[j] = -1;
                    //anchor_dist_target[j] = 0;
                    lol_seq_dist[j] = -1;
                }
            }
            lolscore(lol_dist, lol_seq_dist, lol_score_vec, anchor_length[sa], hidden_layer);
            //lol_score[i] += scoreForward[final_anchor_query[i]][final_anchor_target[i]];

        }

        float total_lol_score = 0.0;
        for (int i = 0; i < anchor_length[sa]; i++) {
            total_lol_score += lol_score_vec[i];

        }
        //TODO add normalization to improve performance
        total_lol_score = total_lol_score; // std::sqrt((float)queryLen * (float)targetLen);
        //total_lol_score *= 100;
        
        if (total_lol_score > max_lol_score){
            max_lol_score = total_lol_score;
            max_lol_idx = sa;
        }
    }
    //std::cout << "max_lol_score: " << max_lol_score << std::endl;
    std::string backtrace = "";
    int matches = 0;
    int q_count = 0;
    int t_count = 0;


    while(matches < anchor_length[max_lol_idx]){
        if(anchor_query[max_lol_idx][q_count] != 0 && anchor_target[max_lol_idx][t_count] != 0){
            backtrace.append(1, 'M');
            matches++;
            q_count++;
            t_count++;
        }
        else if(anchor_query[max_lol_idx][q_count] != 0 && anchor_target[max_lol_idx][t_count] == 0){
            backtrace.append(1, 'D');
            t_count++;
        }
        else if(anchor_query[max_lol_idx][q_count] == 0 && anchor_target[max_lol_idx][t_count] != 0){
            backtrace.append(1, 'I');
            q_count++;
        }
        else if(anchor_query[max_lol_idx][q_count] == 0 && anchor_target[max_lol_idx][t_count] == 0){
            backtrace.append(1, 'D');
            backtrace.append(1, 'I');
            q_count++;
            t_count++;
        }
    }

    //free(G);
    //free(P);
    Matcher::result_t result = Matcher::result_t();
    result.score = max_lol_score; 
    result.eval = max_lol_score;
    result.dbKey = dbKey;
    result.qStartPos = 0;
    result.dbStartPos = 0;
    result.qEndPos = queryLen - 1;
    result.dbEndPos = targetLen - 1;
    result.qLen = queryLen;
    result.dbLen = targetLen;
    result.alnLength = backtrace.length();
    result.backtrace = backtrace;

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
    result.backtrace = backtrace.substr(firstM);
    result.alnLength = (int)result.backtrace.size();


    //std::cout << fwbw_time.count() << " " << lol_score_time.count() << std::endl;
    
    return result;

}
void lolAlign::align_startAnchors(int *anchor_query, int *anchor_target, int max_query, int max_target, int *anchor_length, float** P, float** G){

    for (int i = max_query - this->start_anchor_length; i <= max_query + this->start_anchor_length; ++i)
    {
        int max_ind = max_target + i - max_query;

        anchor_query[i] = 2;
        anchor_target[max_ind] = 2;
        P[i][max_ind] = 0;
        G[i][max_ind] = 0;
        *anchor_length = *anchor_length + 1;

    }
}

void lolAlign::calc_dist_matrix(float *x, float *y, float *z, size_t len, float **d, bool cutoff)
{
    for (size_t i = 0; i < len; i++)
    {
        for (size_t j = 0; j < len; j++)
        {
            d[i][j] = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2) + pow(z[i] - z[j], 2));
            if (d[i][j] >= 15.0 && cutoff)
            {
                d[i][j] = 0;
            }

        }
    }
}


void lolAlign::initQuery(float *x, float *y, float *z, char *querySeq, char *query3diSeq, unsigned int queryLen)
{
    memcpy(query_x, x, sizeof(float) * queryLen);
    memcpy(query_y, y, sizeof(float) * queryLen);
    memcpy(query_z, z, sizeof(float) * queryLen);
    this->queryLen = queryLen;
    this->querySeq = querySeq;
    this->query3diSeq = query3diSeq;
    Coordinates queryCaCords;
    queryCaCords.x = query_x;
    queryCaCords.y = query_y;
    queryCaCords.z = query_z;
    anchor_query = malloc_matrix<int>(num_sa, queryLen);
    G = malloc_matrix<float>(queryLen*2, queryLen*2);
    P = malloc_matrix<float>(queryLen*2, queryLen*2);
    hidden_layer = malloc_matrix<float>(queryLen*2, 3);
    anchor_query = malloc_matrix<int>(num_sa, queryLen);
    anchor_target = malloc_matrix<int>(num_sa, queryLen*2);
    lol_dist = new float[queryLen*2];
    lol_seq_dist = new float[queryLen*2];
    lol_score_vec = new float[queryLen*2];
    final_anchor_query = new int[queryLen*2];
    final_anchor_target = new int[queryLen*2];
    calc_dist_matrix(query_x, query_y, query_z, queryLen, d_ij, true);
    return;

}



void lolAlign::lolmatrix(int *anchor_query, int *anchor_target, int anchor_length, int *gaps, float **d_ij, float **d_kl, float **G, int queryLen, int targetLen, float ** hidden_layer, float * d_dist, float * log_score)
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
            //memset(d_dist, 0, sizeof(float) * (gap1_end - gap1_start));
            //memset(seq_dist, 0, sizeof(float) * (gap1_end - gap1_start));
            if (d_ij[anchor_q][j] > 0)
            {

                seq_dist = std::copysign(1.0f, (anchor_q-j)) * std::log(1 + std::abs((float)(anchor_q-j)));
                for(int l = gap1_start; l < gap1_end; l++)
                {
                    d_dist[l - gap1_start] = std::abs(d_ij[anchor_q][j] - d_kl[anchor_t][l]);

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
    for (int i = 0; i < length; ++i) {
        for (int k = 0; k < 3; ++k) {
            //hidden_layer[i][k] = 0.0;
            hidden_layer[i][k] = d_seq * w1[0][k];
            hidden_layer[i][k] += d_dist[i] * w1[1][k];
            hidden_layer[i][k] += b1[k];
            hidden_layer[i][k] = std::max(0.0f, hidden_layer[i][k]); // ReLU activation
            score[i+start] += hidden_layer[i][k] * w2[k];
        }
        score[i+start] += b2;
    }

}



void lolAlign::lolscore(float* dist, float* d_seq, float* score, int length, float** hidden_layer)
{
    //auto start_t = std::chrono::high_resolution_clock::now();


    for (int i = 0; i < length; ++i) {
        if(d_seq[i] > 0){
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
    //auto end = std::chrono::high_resolution_clock::now();
    //lol_score_time += end - start_t;
}


void lolAlign::computeForwardScoreMatrix(
        char *querySeqAA,
        char *querySeq3Di,
        char *targetSeqAA,
        char *targetSeq3Di,
        int queryLen,
        int targetLen,
        SubstitutionMatrix &subMatAA,
        SubstitutionMatrix &subMat3Di,
        float T,
        float **scoreForward)
{
    unsigned char *queryNumAA = seq2num(querySeqAA, subMatAA.aa2num);
    unsigned char *queryNum3Di = seq2num(querySeq3Di, subMat3Di.aa2num);
    unsigned char *targetNumAA = seq2num(targetSeqAA, subMatAA.aa2num);
    unsigned char *targetNum3Di = seq2num(targetSeq3Di, subMat3Di.aa2num);

    for (int i = 0; i < queryLen; ++i)
    {
        for (int j = 0; j < targetLen; ++j)
        {
            scoreForward[i][j] = static_cast<float>((subMatAA.subMatrix[queryNumAA[i]][targetNumAA[j]] *1.4 ) + (subMat3Di.subMatrix[queryNum3Di[i]][targetNum3Di[j]] * 2.1));
            //scoreForward[i][j] = exp(scoreForward[i][j] / T);
        }
    }
}

void lolAlign::computeDi_score(
        char *querySeqAA,
        char *querySeq3Di,
        char *targetSeqAA,
        char *targetSeq3Di,
        int anchorLen,
        int* final_anchor_query,
        int* final_anchor_target,
        SubstitutionMatrix &subMatAA,
        SubstitutionMatrix &subMat3Di,
        float *scoreForward)
{
    unsigned char *queryNumAA = seq2num(querySeqAA, subMatAA.aa2num);
    unsigned char *queryNum3Di = seq2num(querySeq3Di, subMat3Di.aa2num);
    unsigned char *targetNumAA = seq2num(targetSeqAA, subMatAA.aa2num);
    unsigned char *targetNum3Di = seq2num(targetSeq3Di, subMat3Di.aa2num);

    for (int i = 0; i < anchorLen; ++i)
    {
        int q = final_anchor_query[i];
        int t = final_anchor_target[i];
        scoreForward[i] = static_cast<float>((subMatAA.subMatrix[queryNumAA[q]][targetNumAA[t]] *1.4 ) + (subMat3Di.subMatrix[queryNum3Di[q]][targetNum3Di[t]] * 2.1));        
    }
}


int lolalign(int argc, const char **argv, const Command &command)
{
    std::cout << "lolAlign" << std::endl;
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);
    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++)
    {
        if (par.substitutionMatrices[i].name == "blosum62.out")
        {
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
            "_ca");
    DBReader<unsigned int> qdbr3Di((par.db1 + "_ss").c_str(), (par.db1 + "_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    qdbr3Di.open(DBReader<unsigned int>::NOSORT);

    IndexReader *tdbr = NULL;
    IndexReader *tcadbr = NULL;
    DBReader<unsigned int> tdbr3Di((par.db2 + "_ss").c_str(), (par.db2 + "_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    tdbr3Di.open(DBReader<unsigned int>::NOSORT);
    bool sameDB = false;
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool alignmentIsExtended = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;
    if (par.db1.compare(par.db2) == 0)
    {
        sameDB = true;
        tdbr = &qdbr;
        tcadbr = &qcadbr;
        DBReader<unsigned int> tdbr3Di((par.db1 + "_ss").c_str(), (par.db1 + "_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        tdbr3Di.open(DBReader<unsigned int>::NOSORT);
    }
    else
    {
        tdbr = new IndexReader(par.db2, par.threads,
                               alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                               (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        tcadbr = new IndexReader(
                par.db2,
                par.threads,
                alignmentIsExtended ? IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB2) : IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                alignmentIsExtended ? "_seq_ca" : "_ca");
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbtype = Parameters::DBTYPE_ALIGNMENT_RES;
    if (alignmentIsExtended)
    {
        dbtype = DBReader<unsigned int>::setExtendedDbtype(dbtype, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbtype);
    dbw.open();

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

        char buffer[1024 + 32768];
#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++)
        {
            progress.updateProgress();
            char *data = resultReader.getData(id, thread_idx);
            if (*data != '\0')
            {
                lolAlign lolaln(std::max(qdbr.sequenceReader->getMaxSeqLen() + 1, tdbr->sequenceReader->getMaxSeqLen() + 1), false);
                FwBwAligner fwbwaln(16, -2, -2, 1, 1, 1);
                size_t queryKey = resultReader.getDbKey(id);
                unsigned int queryId = qdbr.sequenceReader->getId(queryKey);
                //std::cout << "start Query " << queryId << std::endl;

                char *querySeq = qdbr.sequenceReader->getData(queryId, thread_idx);
                char *query3diSeq = qdbr3Di.getData(queryId, thread_idx);
                int queryLen = static_cast<int>(qdbr.sequenceReader->getSeqLen(queryId));
                int max_targetLen = ((queryLen*2)/16)*16;
                char *qcadata = qcadbr.sequenceReader->getData(queryId, thread_idx);
                size_t qCaLength = qcadbr.sequenceReader->getEntryLen(queryId);
                float *qdata = qcoords.read(qcadata, queryLen, qCaLength);


                lolaln.initQuery(qdata, &qdata[queryLen], &qdata[queryLen + queryLen], querySeq, query3diSeq, queryLen);
                fwbwaln.resizeMatrix(queryLen, max_targetLen);
                

                


                int passedNum = 0;
                int rejected = 0;

                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected)
                {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int)strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = tdbr->sequenceReader->getId(dbKey);
                    char *targetSeq = tdbr->sequenceReader->getData(targetId, thread_idx);

                    int targetLen = static_cast<int>(tdbr->sequenceReader->getSeqLen(targetId));
                    if (targetLen > max_targetLen)
                    {
                        max_targetLen = ((targetLen*2)/16)*16;
                        
                        lolaln.reallocate_target(max_targetLen);
                        fwbwaln.resizeMatrix(queryLen, max_targetLen);
                    }
                    if (Util::canBeCovered(par.covThr, par.covMode, queryLen, targetLen) == false)
                    {
                        continue;
                    }


                    

                    char *tcadata = tcadbr->sequenceReader->getData(targetId, thread_idx);
                    size_t tCaLength = tcadbr->sequenceReader->getEntryLen(targetId);
                    float *tdata = tcoords.read(tcadata, targetLen, tCaLength);
                    char *target3diSeq = tdbr3Di.getData(targetId, thread_idx);
                    
                    if (targetLen < 20)
                    {
                        //std::cout << "Target too short: " << targetLen << std::endl;
                        continue;
                    }

                    Matcher::result_t result = lolaln.align(dbKey, tdata, &tdata[targetLen], &tdata[targetLen + targetLen], targetSeq, target3diSeq, targetLen, subMatAA, subMat3Di, &fwbwaln);



                    bool hasCov = Util::hasCoverage(par.covThr, par.covMode, 1.0, 1.0);
                    bool hasSeqId = result.seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    //bool hasTMscore = (TMscore >= par.tmScoreThr);

                    if(hasCov && hasSeqId) {
                        swResults.emplace_back(result);
                    }
                
                }


                
                SORT_SERIAL(swResults.begin(), swResults.end(), compareHitsBylolScore);
                for(size_t i = 0; i < swResults.size(); i++){
                    size_t len = Matcher::resultToBuffer(buffer, swResults[i], par.addBacktrace, false);
                    resultBuffer.append(buffer, len);
                }


                dbw.writeData(resultBuffer.c_str(), resultBuffer.size(), queryKey, thread_idx);
                resultBuffer.clear();
                swResults.clear();
                //std::cout << "end Query " << queryId << std::endl;

            }
        }

    }
    dbw.close();
    resultReader.close();
    if (sameDB == false)
    {
        delete tdbr;
        delete tcadbr;
    }
    return EXIT_SUCCESS;
}
#endif
