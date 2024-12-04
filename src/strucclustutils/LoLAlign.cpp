#ifndef LoLAlign
#define LoLAlign
#include <iostream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <numeric>
#include <cmath>
#include <vector>
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





lolAlign::lolAlign(unsigned int maxSeqLen, bool tmAlignFast, bool tmScoreOnly, bool computeExactScore)
   : tmAlignFast(tmAlignFast),
     xtm(maxSeqLen), ytm(maxSeqLen), xt(maxSeqLen),
     r1(maxSeqLen), r2(maxSeqLen), computeExactScore(computeExactScore){
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

lolAlign::~lolAlign(){
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



int lolalign(int argc, const char **argv, const Command &command) {
    std::cout << "lolAlign" << std::endl;    
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

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
    bool sameDB = false;
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool alignmentIsExtended = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
        tcadbr = &qcadbr;
        DBReader<unsigned int> qdbr3Di((par.db1 + "_ss").c_str(), (par.db1 + "_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        qdbr3Di.open(DBReader<unsigned int>::NOSORT);

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
        DBReader<unsigned int> qdbr3Di((par.db2 + "_ss").c_str(), (par.db2 + "_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        qdbr3Di.open(DBReader<unsigned int>::NOSORT);
        
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
