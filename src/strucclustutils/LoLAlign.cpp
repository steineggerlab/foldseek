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



bool compareHitsBylolScore(const Matcher::result_t &first, const Matcher::result_t &second) {
    if (first.eval != second.eval) {
        return first.eval > second.eval;
    }
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
    
    mem = (float *)mem_align(ALIGN_FLOAT, 6 * maxSeqLen * 4 * sizeof(float));
    invmap = new int[maxSeqLen];
}

lolAlign::~lolAlign()
{
    free(query_x);
    free(query_y);
    free(query_z);
    free(target_x);
    free(target_y);
    free(target_z);
    free(mem);
    free(d_ij);
    //free(d_kl);
    delete[] invmap;
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

    for (size_t i = 1; i < numsSize; ++i) {
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
                            char * targetSeq, char* target3diSeq, unsigned int targetLen, SubstitutionMatrix &subMatAA, SubstitutionMatrix &subMat3Di)
{

    int length = 16;
    size_t assignTargetLen = targetLen + (length - targetLen % length) % length;
    int blocks = (int)((assignTargetLen / length));
    float **G = malloc_matrix<float>(queryLen, targetLen + length);
    for (size_t i = 0; i < queryLen; ++i)
    {
        for (size_t j = targetLen; j < targetLen + length; ++j)
        {
            G[i][j] = 1;
        }
    }
    /*for(int i = 0; i<20; i++){
        for(int j = 0; j<20; j++){
            subMatAA.subMatrix[i][j] = (int)subMatAA.subMatrix[i][j] * 1.4;
            subMat3Di.subMatrix[i][j] = (int)subMat3Di.subMatrix[i][j] * 2.1;

        }
    }*/


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
    int* sa_index = new int[num_sa];
    float* sa_scores = new float[num_sa];
    for(int i = 0; i < num_sa; i++){
        sa_index[i] = i;
        sa_scores[i] = 0;
    }
    float ** d_kl = malloc_matrix<float>(targetLen, targetLen);
    calc_dist_matrix(target_x, target_y, target_z, targetLen, d_kl, false);

    
    

    P = malloc_matrix<float>(queryLen, targetLen + length);
    //G = malloc_matrix<float>(queryLen, targetLen + length);
    anchor_target = malloc_matrix<int>(num_sa, targetLen);
    int* anchor_length = new int[num_sa];
    int* new_anchor_length = new int[num_sa];


    
    for(int sa = 0; sa < 10; sa++){
        anchor_length[sa] = 0;
        new_anchor_length[sa] = 0;
        for (int i = 0; i < queryLen; i++) {

            anchor_query[sa][i] = 0;
        }
        for (int i = 0; i < targetLen; i++) {
            anchor_target[sa][i] = 0;
        }
    }


    

    

    int gaps_start[4] = {0, queryLen, 0, targetLen};

    for(int sa = 0; sa < 10; sa++){

        if(sa % 5 == 0){
            lolAlign::lol_fwbw(G, P, queryLen, targetLen, assignTargetLen, start_anchor_go, start_anchor_ge, start_anchor_T, length, blocks, gaps_start);
        }
        maxIndexX = 0;
        maxIndexY = 0;
        
        float maxScore = 0;



        for (int i = this->start_anchor_length; i < queryLen - this->start_anchor_length ; ++i) {
            for (int j = this->start_anchor_length; j < targetLen - this->start_anchor_length ; ++j) {
                
                if (P[i][j] >= maxScore) {
                    maxScore = P[i][j];
                    maxIndexX = i;
                    maxIndexY = j;
                }
            }
        }




        int start_row = maxIndexX - std::min(maxIndexX, maxIndexY);
        int start_col = maxIndexY - std::min(maxIndexX, maxIndexY);
        int diag_length = std::min(queryLen - start_row, targetLen - start_col);
        int* diag_row = new int[diag_length];
        int* diag_col = new int[diag_length];
        float* diag_dist = new float[diag_length];
        //float* diag_col_dist = new float[diag_length];
        float* diag_seq_dist = new float[diag_length];
        float* diag_score = new float[diag_length];
        for(int i = 0; i < diag_length; i++){
            diag_row[i] = start_row + i;
            diag_col[i] = start_col + i;
            diag_score[i] = log(G[start_row + i][start_col + i]) *start_anchor_T;

        }



        for(int i = -start_anchor_length; i < start_anchor_length; i++){
            for(int j=0; j<diag_length; j++){
                if(d_ij[maxIndexX+i][diag_row[j]] > 0){

                    diag_dist[j] = std::abs(d_ij[maxIndexX+i][diag_row[j]] - d_kl[maxIndexY+i][diag_col[j]]);
                    diag_seq_dist[j] = maxIndexX+i -diag_row[j];

                    //std::cout << "diag_row_dist[" << j<< "]: " << diag_row_dist[j] - diag_col_dist[j] << std::endl;
                }
                else{
                    diag_dist[j] = -1;
                    //diag_col_dist[j] = 0;
                    diag_seq_dist[j] = 0;
                }
            }

            lolscore(diag_dist, diag_seq_dist, diag_score, diag_length);

            
        }

        sa_scores[sa] = maxSubArray(diag_score, diag_length); 
        //std::cout << "sa_scores[" << sa << "]: " << sa_scores[sa] << std::endl;



        /*for(int i = 0; i < queryLen; i++){
            if(anchor_query[sa][i] != 0){
                std::cout << "anchor_query[" << i << "] = " << anchor_query[sa][i] << std::endl;
            }
        }

        for(int i = 0; i < targetLen; i++){
            if(anchor_target[sa][i] != 0){
                std::cout << "anchor_target[" << i << "] = " << anchor_target[sa][i] << std::endl;
            }
        }*/
       align_startAnchors(anchor_query[sa], anchor_target[sa], maxIndexX, maxIndexY, &new_anchor_length[sa], P, G);
              


        anchor_length[sa] = new_anchor_length[sa];

    }
    for (size_t i = 0; i < queryLen; ++i)
    {
        for (size_t j = 0; j < targetLen +length; ++j)
        {
            G[i][j] = 0.5;
        }
    }




    for(int i = 0 ; i < queryLen; i++){
        for(int j = 0; j < targetLen; j++){
            P[i][j] = 0;
        }
    }
    index_sort(sa_scores, sa_index, 10);
    
    int* gaps = new int[4]{0, 0, 0, 0};
    int sa; 

    for (int sa_it = 0; sa_it < 3; sa_it++){
        sa = sa_index[9 - sa_it];
        //std::cout << "sa: " << sa << std::endl;
        for(int iteration = 0; iteration < 100; iteration++){


            gaps[0] = 0;
            gaps[1] = 0;
            gaps[2] = 0;
            gaps[3] = 0;
            

            while((gaps[1] < queryLen && gaps[3] < targetLen)){

                calc_gap(anchor_query[sa], anchor_target[sa], gaps, queryLen, targetLen);
                
                if (gaps[0] != -1){
                    lolmatrix(anchor_query[sa], anchor_target[sa], new_anchor_length[sa], gaps, d_ij, d_kl, G, queryLen, targetLen);
                }
                else{
                    break;
                }
                //std::cout << "gaps[0]: " << gaps[0] << " gaps[1]: " << gaps[1] << " gaps[2]: " << gaps[2] << " gaps[3]: " << gaps[3] << std::endl;
                //std::cout << "queryLen: " << queryLen << " targetLen: " << targetLen << std::endl;
                
            }
            for(int i = 0; i <queryLen; i++){
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

            gaps[0] = 0;
            gaps[1] = 0;
            gaps[2] = 0;
            gaps[3] = 0;
            float maxP = 0.5;

            while((gaps[1] < queryLen && gaps[3] < targetLen)){
                calc_gap(anchor_query[sa], anchor_target[sa], gaps, queryLen, targetLen);
                if(gaps[0] != -1){
                    lolAlign::lol_fwbw(G, P, queryLen, targetLen, assignTargetLen, start_anchor_go, start_anchor_ge, 2, length, blocks, gaps);

                }
                //std::cout << "gaps[0]: " << gaps[0] << " gaps[1]: " << gaps[1] << " gaps[2]: " << gaps[2] << " gaps[3]: " << gaps[3] << std::endl;


                if (gaps[0] != -1){
                    for (int i = gaps[0]; i < gaps[1]; i++){
                        for (int j = gaps[2]; j < gaps[3]; j++){
                            if (P[i][j] > maxP){
                                maxIndexX = i;
                                maxIndexY = j;
                                maxP = P[i][j];
                                //std::cout << "maxP: " << maxP <<   std::endl;
                                //std::cout << "maxIndexX: " << maxIndexX << " maxIndexY: " << maxIndexY << std::endl;
                            }
                        }
                    }
                }
                else{
                    break;
                }

            }

            new_anchor_length[sa] = 0;


            
            

            for(int i = 0; i < queryLen; i++){
                for(int j = 0; j < targetLen; j++){
                    if(P[i][j] > maxP -0.1){
                        if(anchor_query[sa][i] == 0 && anchor_target[sa][j] == 0){
                            anchor_query[sa][i] = 2;
                            anchor_target[sa][j] = 2;
                            anchor_length[sa] +=1;
                            new_anchor_length[sa] += 1;
                            //std::cout << "anchor_query[" << i << "]" << "anchor_target[" << j << "] " << P[i][j] << " " << maxP << std::endl;
                        }
                    }
                }
            }
            if (new_anchor_length[sa] == 0){
                break;

            }
        }
        /*if(false){
        
            std::ofstream outfile("/home/lasse/Desktop/Projects/FB_martin/zmForward.txt");
            if (outfile.is_open())
            {
                for (size_t i = 0; i < queryLen; ++i)
                {
                    for (size_t j = 0; j < targetLen; ++j)
                    {
                        outfile << log(G[i][j]) << " ";
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
            for (size_t j = 0; j < targetLen +length; ++j)
            {
                G[i][j] = 0.5;
                P[i][j] = 0;
            }
        }

    }


    float max_lol_score = 0.0;
    int max_lol_idx = 0;
    for (int sa_it = 0; sa_it < 3; sa_it++){
        
        sa = sa_index[9 - sa_it];
    
        float* lol_score = new float[anchor_length[sa]];
        float* seq_distance = new float[anchor_length[sa]];
        float* anchor_dist = new float[anchor_length[sa]];
        //float* anchor_dist_target = new float[anchor_length[sa]];
        for (int i = 0; i < anchor_length[sa]; i++) {
            lol_score[i] = 0.0;
        }
        
        int* final_anchor_query = new int[anchor_length[sa]];
        int* final_anchor_target = new int[anchor_length[sa]];
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
        
        for (int i = 0; i < anchor_length[sa]; i++) {
            for (int j = 0;j < anchor_length[sa]; j++) {
                if(d_ij[final_anchor_query[i]][final_anchor_query[j]] > 0){
                    anchor_dist[j] = std::abs(d_ij[final_anchor_query[i]][final_anchor_query[j]] - d_kl[final_anchor_target[i]][final_anchor_target[j]]);
                    //anchor_dist_target[j] = d_kl[final_anchor_target[i]][final_anchor_target[j]];
                    seq_distance[j] = final_anchor_query[i]-final_anchor_query[j];
                }
                else{
                    anchor_dist[j] = -1;
                    //anchor_dist_target[j] = 0;
                    seq_distance[j] = 0;
                } 
            }
            lolscore(anchor_dist, seq_distance, lol_score, anchor_length[sa]);
            //lol_score[i] += scoreForward[final_anchor_query[i]][final_anchor_target[i]];

        }
        float total_lol_score = 0.0;
        for (int i = 0; i < anchor_length[sa]; i++) {
            total_lol_score += lol_score[i];

        }
        total_lol_score = total_lol_score / sqrt(queryLen * targetLen);
        if (total_lol_score > max_lol_score){
            max_lol_score = total_lol_score;
            max_lol_idx = sa;
        }
    }

    free(G);
    free(P);
    free(d_kl);
    Matcher::result_t result = Matcher::result_t();
    result.score = max_lol_score;
    result.eval = max_lol_score;
    result.dbKey = dbKey;
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
            if (d[i][j] > 15.0 && cutoff)
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
    calc_dist_matrix(query_x, query_y, query_z, queryLen, d_ij, true);
    

}



void lolAlign::lolmatrix(int *anchor_query, int *anchor_target, int anchor_length, int *gaps, float **d_ij, float **d_kl, float **G, int queryLen, int targetLen)
{
    int gap0_start = gaps[0];
    int gap0_end = gaps[1];
    int gap1_start = gaps[2];
    int gap1_end = gaps[3];
    float *d_dist = new float[gap1_end - gap1_start];
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
                float d_ij_row = d_ij[anchor_q][j];
                float * d_kl_row = d_kl[anchor_t];
                seq_dist = std::copysign(1.0f, (anchor_q-j)) * std::log(1 + std::abs((float)(anchor_q-j)));
                for(int l = gap1_start; l < gap1_end; l++)
                {
                    d_dist[l - gap1_start] = std::abs(d_ij_row - d_kl_row[l]);

                    //seq_dist[l - gap1_start] = std::copysign(1.0f, (anchor_q-j)) * std::log(1 + std::abs((float)(anchor_q-j)));
                }
                lolscore(d_dist, seq_dist, G[j], gap1_end - gap1_start, gap1_start);
            }
        }
    }  
}


void lolAlign::lolscore(float* d_dist, float d_seq, float* score, int length, int start)
{
    std::vector<std::vector<float>> x(length, std::vector<float>(2));    

    for (int i = 0; i < length; ++i) {
        x[i][0] = d_seq;
        x[i][1] = d_dist[i];
    }

   std::vector<std::vector<float>> hidden_layer(length, std::vector<float>(3, 0.0f));
    for (int i = 0; i < length; ++i) {
        for (int k = 0; k < 3; ++k) {
            for (int j = 0; j < 2; ++j) {
                hidden_layer[i][k] += x[i][j] * w1[j][k];
            }
            hidden_layer[i][k] += b1[k];
            hidden_layer[i][k] = std::max(0.0f, hidden_layer[i][k]); // ReLU activation
        }
    }
    std::vector<float> log_score(length, 0.0f);
    // Calculate score for each input
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < 3; ++j) {
            log_score[i] += hidden_layer[i][j] * w2[j];
        }
        score[i+start] *= exp(log_score[i] + b2);

        
    }
}



void lolAlign::lolscore(float* dist, float* d_seq, float* score, int length)
{
    std::vector<std::vector<float>> x(length, std::vector<float>(2));

    for (int i = 0; i < length; ++i) {
        x[i][0] = std::copysign(1.0f, d_seq[i]) * std::log(1 + std::abs(d_seq[i]));
        x[i][1] = dist[i];
    }

   std::vector<std::vector<float>> hidden_layer(length, std::vector<float>(3, 0.0f));
    for (int i = 0; i < length; ++i) {
        for (int k = 0; k < 3; ++k) {
            for (int j = 0; j < 2; ++j) {
                hidden_layer[i][k] += x[i][j] * w1[j][k];
            }
            hidden_layer[i][k] += b1[k];
            hidden_layer[i][k] = std::max(0.0f, hidden_layer[i][k]); // ReLU activation
        }
    }
    // Calculate score for each input
    for (int i = 0; i < length; ++i) {
        if(dist[i] >=0){
            for (int j = 0; j < 3; ++j) {
                score[i] += hidden_layer[i][j] * w2[j];
            }
            score[i] += b2;
        }
        
    }
    
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

    for (size_t i = 0; i < queryLen; ++i)
    {
        for (size_t j = 0; j < targetLen; ++j)
        {
            scoreForward[i][j] = static_cast<float>(subMatAA.subMatrix[queryNumAA[i]][targetNumAA[j]] + subMat3Di.subMatrix[queryNum3Di[i]][targetNum3Di[j]]);
            scoreForward[i][j] = exp(scoreForward[i][j] / T);
        }
    }
}

void lolAlign::lol_fwbw(float **scoreForward_s, float **P, size_t queryLen_s, size_t targetLen_s, size_t assignTargetLen, float go, float ge, float T, int length, int blocks, int* gaps)
{

    int queryLen = gaps[1] - gaps[0]; 
    int targetLen = gaps[3] - gaps[2];
    assignTargetLen = targetLen + (length - targetLen % length) % length;
    blocks = (int)((assignTargetLen / length));
    float **scoreForward = allocateMemory(queryLen, assignTargetLen);
    float **scoreBackward = allocateMemory(queryLen, assignTargetLen);
    float **zmForward = allocateMemory(queryLen, assignTargetLen);
    float **zmBackward = allocateMemory(queryLen, assignTargetLen);
    float **zmBlock = allocateMemory(queryLen + 1, length + 1);
    float **zmaxForward = allocateMemory(blocks, queryLen);
    float **zmaxBackward = allocateMemory(blocks, queryLen);
    float *zeBlock = new float[length + 1];
    float *zfBlock = new float[length + 1];

    for (size_t i = 0; i < queryLen; ++i)
    {
        for (size_t j = 0; j < targetLen; ++j)
        {
            scoreForward[i][j] = scoreForward_s[i+gaps[0]][j+gaps[2]];
            scoreBackward[i][j] = scoreForward_s[(queryLen - 1 - i)+gaps[0]][targetLen - 1 - j+gaps[2]];
        }
    }


    float *zInit[3];
    zInit[0] = new float[queryLen];
    zInit[1] = new float[queryLen];
    zInit[2] = new float[queryLen];

    for (unsigned int i = 0; i < queryLen; ++i)
    {
        zInit[0][i] = -DBL_MAX;
        zInit[1][i] = -DBL_MAX;
        zInit[2][i] = -DBL_MAX;
    }

    for (size_t b = 0; b < blocks; ++b)
    {

        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, static_cast<size_t>(targetLen)) - start;

        forwardBackwardSaveBlockMaxLocal(scoreForward, zInit, T, go, ge, queryLen, start, end, memcpy_cols, targetLen,
                                         zmForward, zmaxForward[b], zmBlock, zeBlock, zfBlock, gaps);
    }


    ///////////////////////////////////Backward////////////////////////////////////////

    for (unsigned int i = 0; i < queryLen; ++i)
    {
        zInit[0][i] = -DBL_MAX;
        zInit[1][i] = -DBL_MAX;
        zInit[2][i] = -DBL_MAX;
    }

    for (size_t b = 0; b < blocks; ++b)
    {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, static_cast<size_t>(targetLen)) - start;

        forwardBackwardSaveBlockMaxLocal(scoreBackward, zInit, T, go, ge, queryLen, start, end, memcpy_cols, targetLen,
                                         zmBackward, zmaxBackward[b], zmBlock, zeBlock, zfBlock, gaps);
    }

    ///////////////////////////////////Rescale////////////////////////////////////////
    // Rescale the values by the maximum in the log space for each block
    // This turns the matrix into log space
  

    rescaleBlocks(zmForward, zmaxForward, queryLen, length, blocks, targetLen);
    rescaleBlocks(zmBackward, zmaxBackward, queryLen, length, blocks, targetLen);


    /*std::ofstream outfile("/home/lasse/Desktop/Projects/FB_martin/zmForward.txt");
    if (outfile.is_open())
    {
        for (size_t i = 0; i < queryLen; ++i)
        {
            for (size_t j = 0; j < assignTargetLen; ++j)
            {
                outfile << zmForward[i][j] << " ";
            }
            outfile << "\n";
        }
        outfile.close();
    }
    else
    {
        std::cerr << "Unable to open file for writing P matrix." << std::endl;
    }*/

    float max_zm = -std::numeric_limits<float>::max();
    for (size_t i = 0; i < queryLen; ++i)
    {
        for (size_t j = 0; j < targetLen; ++j)
        {
            max_zm = std::max(max_zm, zmForward[i][j]);
        }
    }

    float sum_exp = 0.0;
    for (size_t i = 0; i < queryLen; ++i)
    {
        for (size_t j = 0; j < targetLen; ++j)
        {
            sum_exp += exp(zmForward[i][j] - max_zm);
        }
    }
    float logsumexp_zm = max_zm + log(sum_exp);

    for (size_t i = 0; i < queryLen; ++i)
    {
        for (size_t j = 0; j < targetLen; ++j)
        {
            P[i+gaps[0]][j+gaps[2]] = exp(
                zmForward[i][j] + zmBackward[queryLen - 1 - i][targetLen - 1 - j] - log(scoreForward_s[i+gaps[0]][j+gaps[2]]) - logsumexp_zm);
        }
    }

    free(scoreBackward);
    free(zmForward);
    free(zmBackward);
    free(zmBlock);
    free(zmaxForward);
    free(zmaxBackward);
    free(zeBlock);
    free(zfBlock);
    return;
}

void lolAlign::forwardBackwardSaveBlockMaxLocal(float** S, float** z_init,
                                                float T, float go, float ge,
                                                size_t rows, size_t start, size_t end, size_t memcpy_cols, size_t targetlen,
                                                float **zm, float *zmax, float **zmBlock, float *zeBlock, float *zfBlock, int* gaps)
{
    float exp_go = exp(go / T);
    float exp_ge = exp(ge / T);

    memset(zeBlock, 0, (end - start + 1) * sizeof(float));
    memset(zfBlock, 0, (end - start + 1) * sizeof(float));

    std::vector<float> ze_first(rows + 1, 0);
    std::vector<float> zf_first(rows + 1, 0);

    // Init blocks
    memset(zmBlock[0], 0, (end - start + 1) * sizeof(float));

    for (size_t i = 0; i < rows; ++i)
    {
        zmBlock[i + 1][0] = z_init[0][i];
        ze_first[i + 1] = z_init[1][i];
        zf_first[i + 1] = z_init[2][i];
    }

    size_t cols = memcpy_cols;

    float current_max = 0;

    for (size_t i = 1; i <= rows; ++i)
    {
        if (i != 1)
        {
            zmBlock[i - 1][0] = exp(zmBlock[i - 1][0]);
            ze_first[i - 1] = exp(ze_first[i - 1]);
            zf_first[i - 1] = exp(zf_first[i - 1]);
            // Debug(Debug::INFO) << zmBlock[i - 1][0] << '\t';
        }
        const float expMax = exp(-current_max);
        // #pragma omp simd
        for (size_t j = 1; j <= cols; ++j)
        {
            // std::cout << rows << " " << cols << "\n" << std::endl;
            if (j == 1)
            {
                float tmp = (zmBlock[i - 1][j - 1] + ze_first[i - 1] + zf_first[i - 1] + expMax);
                zmBlock[i][j] = tmp * S[(i - 1)][start + j - 1];
            }
            else
            {
                float tmp = (zmBlock[i - 1][j - 1] + zeBlock[j - 1] + zfBlock[j - 1] + expMax);
                zmBlock[i][j] = tmp * S[i - 1][start + j - 1];
            }
        }

#pragma omp simd
        for (size_t j = 1; j <= cols; ++j)
        {
            if (j == 1)
            {
                zeBlock[j] = exp(zmBlock[i][j - 1]) * exp_go + exp(ze_first[i]) * exp_ge;
            }
            else
            {
                zeBlock[j] = zmBlock[i][j - 1] * exp_go + zeBlock[j - 1] * exp_ge;
            }
        }
#pragma omp simd
        for (size_t j = 1; j <= cols; ++j)
        {
            zfBlock[j] = zmBlock[i - 1][j] * exp_go + zfBlock[j] * exp_ge;
        }

        float z_temp_m = *std::max_element(zmBlock[i] + 1, zmBlock[i] + cols + 1);
        float z_temp_e = *std::max_element(zeBlock, zeBlock + cols + 1);
        float z_temp_f = *std::max_element(zfBlock + 1, zfBlock + cols + 1);
        float z_temp = std::max(z_temp_m, std::max(z_temp_e, z_temp_f));
        zmax[i - 1] = log(z_temp);
        current_max += zmax[i - 1];
#pragma omp simd
        for (size_t j = 1; j <= cols; ++j)
        {
            zmBlock[i][j] /= z_temp;
            zeBlock[j] /= z_temp;
            zfBlock[j] /= z_temp;
        }

        zmBlock[i][0] -= zmax[i - 1];
        ze_first[i] -= zmax[i - 1];
        zf_first[i] -= current_max;
        if (i < rows)
        {
            zmBlock[i + 1][0] -= current_max;
            ze_first[i + 1] -= current_max;
            z_init[0][i - 1] = log(zmBlock[i][cols]) + current_max;
            z_init[1][i - 1] = log(zeBlock[cols]) + current_max;
            z_init[2][i - 1] = log(zfBlock[cols]) + current_max;
        }
    }

    std::vector<float> rescale(rows);
    std::partial_sum(zmax, zmax + rows, rescale.begin());

    for (size_t i = 0; i < rows; ++i)
    {
        memcpy(zm[i] + start, zmBlock[i + 1] + 1, memcpy_cols * sizeof(float));
    }
}

void lolAlign::rescaleBlocks(float **matrix, float **scale, size_t rows, size_t length, size_t blocks, size_t targetLen)
{
    for (size_t b = 0; b < blocks; ++b)
    {
        size_t start = b * length;
        size_t end = std::min((b + 1) * length, targetLen);
        std::vector<float> cumsum(rows);
        std::partial_sum(scale[b], scale[b] + rows, cumsum.begin());

        for (size_t i = 0; i < rows; ++i)
        {
            for (size_t j = start; j < end; ++j)
            {
                //matrix[i][j] = log(matrix[i][j] + std::numeric_limits<float>::min()) + cumsum[i];
                matrix[i][j] = log(matrix[i][j]) + cumsum[i];
            }
        }
    }
}

float **lolAlign::allocateMemory(size_t rows, size_t cols)
{
    // Allocate memory for an array of pointers to rows
    float **array = (float **)malloc(rows * sizeof(float *));
    if (array == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Allocate memory for each row
    for (int i = 0; i < rows; i++)
    {
        array[i] = (float *)malloc(cols * sizeof(float));
        if (array[i] == NULL)
        {
            fprintf(stderr, "Memory allocation failed\n");
            exit(1);
        }
    }

    return array;
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
        lolAlign lolaln(std::max(qdbr.sequenceReader->getMaxSeqLen() + 1, tdbr->sequenceReader->getMaxSeqLen() + 1), false);

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
                size_t queryKey = resultReader.getDbKey(id);

                unsigned int queryId = qdbr.sequenceReader->getId(queryKey);
                char *querySeq = qdbr.sequenceReader->getData(queryId, thread_idx);
                char *query3diSeq = qdbr3Di.getData(queryId, thread_idx);
                int queryLen = static_cast<int>(qdbr.sequenceReader->getSeqLen(queryId));
                char *qcadata = qcadbr.sequenceReader->getData(queryId, thread_idx);
                size_t qCaLength = qcadbr.sequenceReader->getEntryLen(queryId);
                float *qdata = qcoords.read(qcadata, queryLen, qCaLength);
                lolaln.initQuery(qdata, &qdata[queryLen], &qdata[queryLen + queryLen], querySeq, query3diSeq, queryLen);

                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected)
                {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int)strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = tdbr->sequenceReader->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB)) ? true : false;
                    if (isIdentity == true)
                    {
                        backtrace.append(SSTR(queryLen));
                        backtrace.append(1, 'M');
                        Matcher::result_t result(dbKey, 0, 1.0, 1.0, 1.0, 1.0, std::max(queryLen, queryLen), 0, queryLen - 1, queryLen, 0, queryLen - 1, queryLen, backtrace);
                        size_t len = Matcher::resultToBuffer(buffer, result, par.addBacktrace, false);
                        resultBuffer.append(buffer, len);
                        backtrace.clear();
                        continue;
                    }
                    char *targetSeq = tdbr->sequenceReader->getData(targetId, thread_idx);
                    int targetLen = static_cast<int>(tdbr->sequenceReader->getSeqLen(targetId));
                    if (Util::canBeCovered(par.covThr, par.covMode, queryLen, targetLen) == false)
                    {
                        continue;
                    }

                    char *tcadata = tcadbr->sequenceReader->getData(targetId, thread_idx);
                    size_t tCaLength = tcadbr->sequenceReader->getEntryLen(targetId);
                    float *tdata = tcoords.read(tcadata, targetLen, tCaLength);
                    char *target3diSeq = tdbr3Di.getData(targetId, thread_idx);
                    Matcher::result_t result = lolaln.align(dbKey, tdata, &tdata[targetLen], &tdata[targetLen + targetLen], targetSeq, target3diSeq, targetLen, subMatAA, subMat3Di);
                    swResults.emplace_back(result);
                    SORT_SERIAL(swResults.begin(), swResults.end(), compareHitsBylolScore);
                    for(size_t i = 0; i < swResults.size(); i++){
                        size_t len = Matcher::resultToBuffer(buffer, swResults[i], par.addBacktrace, false);
                        resultBuffer.append(buffer, len);
                    }

                    dbw.writeData(resultBuffer.c_str(), resultBuffer.size(), queryKey, thread_idx);
                    resultBuffer.clear();
                    swResults.clear();


                    /*std::ofstream outfile("/home/lasse/Desktop/Projects/FB_martin/P_mat.txt");
                    if (outfile.is_open()) {
                        for (size_t i = 0; i < queryLen; ++i) {
                            for (size_t j = 0; j < assignTargetLen; ++j) {
                                outfile << P[i][j] << " ";
                            }
                            outfile << "\n";
                        }
                        outfile.close();
                    } else {
                        std::cerr << "Unable to open file for writing P matrix." << std::endl;
                    }*/
                }
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
