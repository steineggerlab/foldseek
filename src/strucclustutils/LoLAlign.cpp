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





lolAlign::lolAlign(unsigned int maxSeqLen, bool computeExactScore)
   : xtm(maxSeqLen), ytm(maxSeqLen), xt(maxSeqLen),
     r1(maxSeqLen), r2(maxSeqLen), computeExactScore(computeExactScore){
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

lolAlign::~lolAlign(){
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


void lolAlign::initQuery(float *x, float *y, float *z, char * querySeq, char* query3diSeq, unsigned int queryLen){
    memset(querySecStruc, 0, sizeof(char) * queryLen);
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
}



float lolAlign::alignment_lolScore(std::vector<float> d_ij, std::vector<float> d_kl, std::vector<float> anchorpoints, size_t anchor_length){
    std::vector<float> score;

    float score_sum = 0.0;
    for(size_t an = 0; an < anchor_length; an++){
        int i = anchorpoints[an];
        int k = anchorpoints[an + anchor_length];

        for(size_t an2 = 0; an2 <  anchor_length; an++){
            int j = anchorpoints[an2];
            int l = anchorpoints[an2 + anchor_length];
            lolscore(d_ij, d_kl, d_kl, score);
            score_sum += score[0];
        }
    }
    return score_sum;
}




void lolAlign::lolmatrix(int *anchor, int anchor_length, int *gaps, float **d_ij, float **d_kl, float **G){
    for (int i = 0; i < anchor_length; i++) {
        int anchor_q = anchor[i];
        int anchor_t = anchor[i + anchor_length];
        for(int j = gaps[0]; j<= gaps[1]; j++){
            for(int l = gaps[2]; l<= gaps[3]; l++){    
                G[j][l] = d_ij[anchor_q][j] + d_kl[anchor_q][l];
                
            }
        }
    }
}

void lolAlign::lolscore(std::vector<float> d_ij, std::vector<float> d_kl, std::vector<float> d_seq, std::vector<float> score) {
    size_t n = d_ij.size();
    std::vector<float> seq_distance_feat(n);
    std::vector<std::vector<float>> x(n, std::vector<float>(2));

    // Compute seq_distance_feat and x
    for (size_t i = 0; i < n; ++i) {
        seq_distance_feat[i] = std::copysign(1.0, d_seq[i]) * std::log(1.0 + std::abs(d_seq[i]));
        x[i][0] = seq_distance_feat[i];
        x[i][1] = std::abs(d_ij[i] - d_kl[i]);
    }

    std::vector<std::vector<float>> hidden_layer(n, std::vector<float>(3, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (int j = 0; j < 3; ++j) {
            hidden_layer[i][j] = b1[j];
            for (int k = 0; k < 2; ++k) {
                hidden_layer[i][j] += x[i][k] * w1[j][k];
            }
            // ReLU activation
            hidden_layer[i][j] = (0.0 <  hidden_layer[i][j]) ? hidden_layer[i][j] : 0.0;
        }
    }
    for (size_t i = 0; i < n; ++i) {
        score[i] = b2;
        for (int j = 0; j < 3; ++j) {
            score[i] += hidden_layer[i][j] * w2[j];
        }
    }
    return;
}


void lolAlign::computeForwardScoreMatrix(
    char* querySeqAA,
    char* querySeq3Di,
    char* targetSeqAA,
    char* targetSeq3Di,
    int queryLen,
    int targetLen,
    SubstitutionMatrix &subMatAA,
    SubstitutionMatrix &subMat3Di,
    float T,
    float** scoreForward
) {
    unsigned char* queryNumAA = seq2num(querySeqAA, subMatAA.aa2num);
    unsigned char* queryNum3Di = seq2num(querySeq3Di, subMat3Di.aa2num);
    unsigned char* targetNumAA = seq2num(targetSeqAA, subMatAA.aa2num);
    unsigned char* targetNum3Di = seq2num(targetSeq3Di, subMat3Di.aa2num);
    
    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            scoreForward[i][j] = static_cast<float>(subMatAA.subMatrix[queryNumAA[i]][targetNumAA[j]] + subMat3Di.subMatrix[queryNum3Di[i]][targetNum3Di[j]]);
            scoreForward[i][j] = exp(scoreForward[i][j] / T);
        }
    }
    
}

void lolAlign::lol_fwbw(float** scoreForward, float** P, size_t queryLen, size_t targetLen, size_t assignTargetLen, float go, float ge, float T, int length, int blocks) {
    
    float** scoreBackward = allocateMemory(queryLen, assignTargetLen);
    float** zmForward = allocateMemory(queryLen, assignTargetLen);
    float** zmBackward = allocateMemory(queryLen, assignTargetLen);
    float** zmBlock = allocateMemory(queryLen + 1, length + 1);
    float** zmaxForward = allocateMemory(blocks, queryLen);
    float** zmaxBackward = allocateMemory(blocks, queryLen);
    float* zeBlock = new float[length + 1];
    float* zfBlock = new float[length + 1];
    

    for(size_t i = 0; i < queryLen; ++i){
        for(size_t j = 0; j < assignTargetLen; ++j){
            scoreBackward[i][j] = scoreForward[(queryLen - 1 - i)] [targetLen - 1 - j];
        }
    }

    

    float* zInit[3];
    zInit[0] = new float[queryLen];
    zInit[1] = new float[queryLen];
    zInit[2] = new float[queryLen];


    for (unsigned int i=0; i < queryLen; ++i){
        zInit[0][i] = -DBL_MAX;
        zInit[1][i] = -DBL_MAX;
        zInit[2][i] = -DBL_MAX;
    }
 
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, static_cast<size_t>(targetLen)) - start;

        forwardBackwardSaveBlockMaxLocal(scoreForward, zInit, T, go, ge, queryLen, start, end, memcpy_cols, targetLen,
                                         zmForward, zmaxForward[b], zmBlock ,zeBlock, zfBlock);
        
    }
    

    ///////////////////////////////////Backward////////////////////////////////////////

    for (unsigned int i=0; i < queryLen; ++i){
        zInit[0][i] = -DBL_MAX;
        zInit[1][i] = -DBL_MAX;
        zInit[2][i] = -DBL_MAX;
    }

    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, static_cast<size_t>(targetLen)) - start;
        
        forwardBackwardSaveBlockMaxLocal(scoreBackward, zInit, T, go, ge, queryLen, start, end, memcpy_cols, targetLen,
                                         zmBackward, zmaxBackward[b], zmBlock, zeBlock, zfBlock);
    }

    ///////////////////////////////////Rescale////////////////////////////////////////
    // Rescale the values by the maximum in the log space for each block
    // This turns the matrix into log space

    
    rescaleBlocks(zmForward, zmaxForward, queryLen, length, blocks, targetLen);
    rescaleBlocks(zmBackward, zmaxBackward, queryLen, length, blocks, targetLen);


    /*std::ofstream outfile("/home/lasse/Desktop/Projects/FB_martin/zmForward.txt");
    if (outfile.is_open()) {
        for (size_t i = 0; i < queryLen; ++i) {
            for (size_t j = 0; j < assignTargetLen; ++j) {
                outfile << zmBackward[i][j] << " ";
            }
            outfile << "\n";
        }
        outfile.close();
    } else {
        std::cerr << "Unable to open file for writing P matrix." << std::endl;
    }*/
    

    float max_zm = -std::numeric_limits<float>::max();
    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            max_zm = std::max(max_zm, zmForward[i][j]);
        }
    }

    float sum_exp = 0.0;
    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            sum_exp += exp(zmForward[i][j] - max_zm);
        }
    }
    float logsumexp_zm = max_zm + log(sum_exp);



    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            P[i][j] = exp(
                zmForward[i][j]
                + zmBackward[queryLen - 1 - i][targetLen - 1 - j]
                - log(scoreForward[i][j]) 
                - logsumexp_zm
            );
            
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
                                                   float** zm, float* zmax, float** zmBlock, float* zeBlock, float* zfBlock) {
    float exp_go = exp(go / T);
    float exp_ge = exp(ge / T);

    memset(zeBlock, 0, (end - start + 1) * sizeof(float)); 
    memset(zfBlock, 0, (end - start + 1) * sizeof(float)); 
 
    std::vector<float> ze_first(rows+1, 0);
    std::vector<float> zf_first(rows+1, 0);
    
    //Init blocks
    memset(zmBlock[0], 0, (end - start + 1) * sizeof(float));

    for (size_t i = 0; i < rows; ++i) {
        zmBlock[i+1][0] = z_init[0][i];
        ze_first[i+1] = z_init[1][i];
        zf_first[i+1] = z_init[2][i];
    }


    size_t cols = memcpy_cols;


    float current_max = 0;

    for (size_t i = 1; i <= rows; ++i) {
        if (i != 1) {
            zmBlock[i - 1][0] = exp(zmBlock[i - 1][0]);
            ze_first[i - 1] = exp(ze_first[i - 1]);
            zf_first[i - 1] = exp(zf_first[i - 1]);
            // Debug(Debug::INFO) << zmBlock[i - 1][0] << '\t';
        }
        const float expMax = exp(-current_max);
        //#pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
        //std::cout << rows << " " << cols << "\n" << std::endl;
            if (j == 1){
                float tmp = (zmBlock[i - 1][j - 1] + ze_first[i-1] + zf_first[i - 1] + expMax);
                zmBlock[i][j] = tmp * S[(i - 1)][start + j - 1];
            }
            else{
                float tmp = (zmBlock[i - 1][j - 1] + zeBlock[j - 1] + zfBlock[j - 1] + expMax);
                zmBlock[i][j] = tmp * S[i - 1][start + j - 1];
            }
        }
        

        

        #pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
            if (j == 1) {
                zeBlock[j] = exp(zmBlock[i][j - 1]) * exp_go + exp(ze_first[i]) * exp_ge;
            } else {
                zeBlock[j] = zmBlock[i][j - 1] * exp_go + zeBlock[j - 1] * exp_ge;
            }

        }
        #pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
                zfBlock[j] = zmBlock[i - 1][j] * exp_go + zfBlock[j] * exp_ge;
        }

        float z_temp = *std::max_element(zmBlock[i] + 1, zmBlock[i] + cols + 1);
        zmax[i-1] = log(z_temp);
        current_max += zmax[i-1];
        #pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
            zmBlock[i][j] /= z_temp;
            zeBlock[j] /= z_temp;
            zfBlock[j] /= z_temp;
        }
        
        
        zmBlock[i][0] -= zmax[i-1];
        ze_first[i] -= zmax[i-1];
        zf_first[i] -= current_max;
        if (i < rows) {
            zmBlock[i+1][0] -= current_max;
            ze_first[i+1] -= current_max;
            z_init[0][i-1] = log(zmBlock[i][cols]) + current_max;
            z_init[1][i-1] = log(zeBlock[cols]) + current_max;
            z_init[2][i-1] = log(zfBlock[cols]) + current_max;

        }

    }

    std::vector<float> rescale(rows);
    std::partial_sum(zmax, zmax + rows, rescale.begin());


    for (size_t i = 0; i < rows; ++i) {
        memcpy(zm[i] + start, zmBlock[i+1]+1, memcpy_cols * sizeof(float));
    }
}

void lolAlign::rescaleBlocks(float **matrix, float **scale, size_t rows, size_t length, size_t blocks, size_t targetLen){
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = std::min((b + 1) * length, targetLen);
        std::vector<float> cumsum(rows);
        std::partial_sum(scale[b], scale[b] + rows, cumsum.begin());

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = start; j < end; ++j) {
                matrix[i][j] = log(matrix[i][j]) + cumsum[i];
            }
        }
    }
}


float** lolAlign::allocateMemory(size_t rows, size_t cols) {
    // Allocate memory for an array of pointers to rows
    float** array = (float**)malloc(rows * sizeof(float*));
    if (array == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Allocate memory for each row
    for (int i = 0; i < rows; i++) {
        array[i] = (float*)malloc(cols * sizeof(float));
        if (array[i] == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(1);
        }
    }

    return array;
}









int lolalign(int argc, const char **argv, const Command &command) {
    std::cout << "lolAlign" << std::endl;    
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);
    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == "blosum62.out") {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
            std::string matrixName = par.substitutionMatrices[i].name;
            char * serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
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
    DBReader<unsigned int> qdbr3Di((par.db1 + "_ss").c_str(), (par.db1 + "_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    qdbr3Di.open(DBReader<unsigned int>::NOSORT);


    IndexReader *tdbr = NULL;
    IndexReader *tcadbr = NULL;
    DBReader<unsigned int> tdbr3Di((par.db2 + "_ss").c_str(), (par.db2 + "_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    tdbr3Di.open(DBReader<unsigned int>::NOSORT);
    bool sameDB = false;
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool alignmentIsExtended = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
        tcadbr = &qcadbr;
        DBReader<unsigned int> tdbr3Di((par.db1 + "_ss").c_str(), (par.db1 + "_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        tdbr3Di.open(DBReader<unsigned int>::NOSORT);

    } else {
        tdbr = new IndexReader(par.db2, par.threads,
                        alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                        (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        tcadbr = new IndexReader(
                par.db2,
                par.threads,
                alignmentIsExtended ?  IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB2) :
                IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                alignmentIsExtended ? "_seq_ca" : "_ca"
        );
        
        
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    

    int dbtype =  Parameters::DBTYPE_ALIGNMENT_RES;
    if(alignmentIsExtended){
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
        lolAlign lolaln(std::max(qdbr.sequenceReader->getMaxSeqLen() + 1,tdbr->sequenceReader->getMaxSeqLen() + 1), false);

        std::vector<Matcher::result_t> swResults;
        swResults.reserve(300);
        std::string backtrace;
        std::string resultBuffer;
        resultBuffer.reserve(1024*1024);
        Coordinate16 qcoords;
        Coordinate16 tcoords;

        char buffer[1024+32768];
#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            char *data = resultReader.getData(id, thread_idx);
            if(*data != '\0') {
                size_t queryKey = resultReader.getDbKey(id);

                unsigned int queryId = qdbr.sequenceReader->getId(queryKey);
                char *querySeq = qdbr.sequenceReader->getData(queryId, thread_idx);
                char *query3diSeq = qdbr3Di.getData(queryId, thread_idx);
                int queryLen = static_cast<int>(qdbr.sequenceReader->getSeqLen(queryId));
                char *qcadata = qcadbr.sequenceReader->getData(queryId, thread_idx);
                size_t qCaLength = qcadbr.sequenceReader->getEntryLen(queryId);
                float* qdata = qcoords.read(qcadata, queryLen, qCaLength);
                lolaln.initQuery(qdata, &qdata[queryLen], &qdata[queryLen+queryLen], querySeq, query3diSeq, queryLen);

                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = tdbr->sequenceReader->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;
                    if(isIdentity == true){
                        backtrace.append(SSTR(queryLen));
                        backtrace.append(1, 'M');
                        Matcher::result_t result(dbKey, 0 , 1.0, 1.0, 1.0, 1.0, std::max(queryLen,queryLen), 0, queryLen-1, queryLen, 0, queryLen-1, queryLen, backtrace);
                        size_t len = Matcher::resultToBuffer(buffer, result, par.addBacktrace, false);
                        resultBuffer.append(buffer, len);
                        backtrace.clear();
                        continue;
                    }
                    char * targetSeq = tdbr->sequenceReader->getData(targetId, thread_idx);
                    int targetLen = static_cast<int>(tdbr->sequenceReader->getSeqLen(targetId));
                    if(Util::canBeCovered(par.covThr, par.covMode, queryLen, targetLen)==false){
                        continue;
                    }

                    char *tcadata = tcadbr->sequenceReader->getData(targetId, thread_idx);
                    size_t tCaLength = tcadbr->sequenceReader->getEntryLen(targetId);
                    float* tdata = tcoords.read(tcadata, targetLen, tCaLength);
                    char* target3diSeq = tdbr3Di.getData(targetId, thread_idx);
                    int length = 16;
                    size_t assignTargetLen = targetLen + (length - targetLen % length) % length;
                    int blocks = (int)((assignTargetLen/length));
                    float ** scoreForward = malloc_matrix<float>(queryLen, assignTargetLen);
                    for (size_t i = 0; i < queryLen; ++i) {
                        for (size_t j = targetLen; j < assignTargetLen; ++j) {
                            scoreForward[i][j] = 0;
                        }
                    }
                    
                    lolaln.computeForwardScoreMatrix(
                        querySeq,
                        query3diSeq,
                        targetSeq,
                        target3diSeq,
                        queryLen,
                        targetLen,
                        subMatAA,
                        subMat3Di,
                        10.0,
                        scoreForward
                    );
                    
                    float** P = malloc_matrix<float>(queryLen, assignTargetLen);

                    lolaln.lol_fwbw(scoreForward, P, queryLen, targetLen, assignTargetLen, -10.0, -1.0, 2.0, length, blocks);

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
    if(sameDB == false){
        delete tdbr;
        delete tcadbr;
    }
    return EXIT_SUCCESS;
}
#endif
