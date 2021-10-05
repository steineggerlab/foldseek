
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

#ifdef OPENMP
#include <omp.h>
#endif

#include <cmath>

// static two dimensional array instead of vectors
void findNearestNeighbour(char * nn, Coordinates & ca, int length, unsigned char * seqInt){
    // (char*)mem_align(ALIGN_INT, par.maxSeqLen * sizeof(char) * 4 );

    // calculate pairwise matrix (euclidean distance)
    std::vector<vector<float>> nnDist;
    std::vector<vector<int>> nnPos;
    std::vector<int> nnInt;
    float seqDist;

    for(int i = 0; i < length; i++){

        std::vector<std::pair<float, int>> seqDistList;
        for(int j = 0; j < length; j++){
            if(i != j) {
                // calculate distance
                seqDist = sqrt(
                        ((ca.x[i] - ca.x[j]) * (ca.x[i] - ca.x[j])) + ((ca.y[i] - ca.y[j]) * (ca.y[i] - ca.y[j])) +
                        ((ca.z[i] - ca.z[j]) * (ca.z[i] - ca.z[j])));
                seqDistList.emplace_back(make_pair(seqDist, j));
            }
        }
        // find the four nearest neighbours for each amino acid
        std::sort(seqDistList.begin(), seqDistList.end());
        int pos[4];

        for(int m  = 0; m < 4; m++){
            pos[m]  = seqDistList[m].second;
        }
        std::sort(begin(pos), end(pos));
        for(int n = 0; n < 4; n++){
            int neighbour = static_cast<int>(seqInt[pos[n]]);
            nn[i*4 + n] = neighbour;
        }
        seqDistList.clear();
    }
}

short needlemanWunschScore(int subQNNi[4], int subTNNi[4], SubstitutionMatrix *subMat, int nwGapPenalty){

    int scoringMatrix[5][5] = {{0, 0, 0, 0, 0},{0, 0, 0, 0, 0},{0, 0, 0, 0, 0},{0, 0, 0, 0, 0},{0, 0, 0, 0, 0}};

    for(int i = 1; i < 5; i++){
        for(int j = 1; j < 5; j++){

            int score1 = subQNNi[i - 1];
            int score2 = subTNNi[j - 1];
            short scoreTest;

//            if(score1 < 0 or score1 > 15){cout << " NO " << score1;}
//            if(score2 < 0 or score2 > 15){cout << " NO " << score2;}

            scoreTest = subMat->subMatrix[score1][score2];

            // calculate scores - bonus for base pairing and penalty for not
            scoringMatrix[i][j] = std::max(scoringMatrix[i-1][j-1] + scoreTest, scoringMatrix[i-1][j] - nwGapPenalty);
            scoringMatrix[i][j] = std::max(scoringMatrix[i][j-1] - nwGapPenalty, scoringMatrix[i][j]);

        }
    }

    // for testing: print matrix
//    cout << "Scoring Matrix" << endl;
//    for (int a = 0; a < 5; a++) {
//        for(int b = 0; b < 5; b++){
//            cout << scoringMatrix[a][b] << " ";
//        }
//        cout << endl;
//    }
    short score = scoringMatrix[4][4]; /// 4;

    return score;
}

Matcher::result_t alignByNN(char * querynn, unsigned char *querySeqInt, int queryLen, char * targetnn, unsigned char *targetSeqInt, int targetLen, SubstitutionMatrix *subMat, int gapOpen, int gapExtern, int gapNW, float nnWeight){

    Matcher::result_t result;
    // gotoh sw itself

        struct scores{ short H, E, F; };
        uint16_t max_score = 0;

        scores *workspace = new scores[queryLen * 2 + 2];
        scores *curr_sHEF_vec = &workspace[0];
        scores *prev_sHEF_vec = &workspace[queryLen + 1];

        int subTNNi[4], subQNNi[4];
        // top row need to be set to a 0 score
        memset(prev_sHEF_vec, 0, sizeof(scores) * (queryLen + 1));
        for (int i = 0; i < targetLen; i++) {

            for(int a = 0; a < 4; a++){
                int neighbour = static_cast<int>(targetnn[4*i + a]);
                subTNNi[a] =  neighbour;
            }

            // left outer column need to be set to a 0 score
            prev_sHEF_vec[0].H = 0; prev_sHEF_vec[0].E = 0; prev_sHEF_vec[0].F = 0;
            curr_sHEF_vec[0].H = 0; curr_sHEF_vec[0].E = 0; curr_sHEF_vec[0].E = 0;
            for (int j = 1; j <= queryLen; j++) {
                for(int a = 0; a < 4; a++){
                    int neighbour = static_cast<int>(querynn[4*(j-1) + a]);
                    subQNNi[a] = neighbour;
                }

                curr_sHEF_vec[j].E = std::max(curr_sHEF_vec[j - 1].H - gapOpen, curr_sHEF_vec[j - 1].E - gapExtern); // j-1
                curr_sHEF_vec[j].F = std::max(prev_sHEF_vec[j].H - gapOpen, prev_sHEF_vec[j].F - gapExtern); // i-1
                short bla = needlemanWunschScore(subTNNi, subQNNi, subMat, gapNW);
                bla = nnWeight * bla;
                int subOne = static_cast<int>(targetSeqInt[i]);
                int subTwo = static_cast<int>(querySeqInt[j - 1]);
//                if(subOne < 0 or subOne > 15){cout << " NO " << subOne;}
//                if(subTwo < 0 or subTwo > 15){cout << " NO " << subTwo;}
                int test_shit = 0;
                if(subOne == subTwo){test_shit = 2;}else{test_shit = -1;}
                const short tempH = prev_sHEF_vec[j - 1].H + subMat->subMatrix[subOne][subTwo] + bla; // i - 1, j - 1
//                const short tempH = prev_sHEF_vec[j - 1].H + subMat->subMatrix[subOne][subTwo]  + bla; // i - 1, j - 1

                curr_sHEF_vec[j].H = std::max(tempH, curr_sHEF_vec[j].E);
                curr_sHEF_vec[j].H = std::max(curr_sHEF_vec[j].H, curr_sHEF_vec[j].F);
                curr_sHEF_vec[j].H = std::max(curr_sHEF_vec[j].H, static_cast<short>(0));
                max_score = static_cast<uint16_t>(std::max(static_cast<uint16_t>(curr_sHEF_vec[j].H), max_score));
            }
            std::swap(prev_sHEF_vec, curr_sHEF_vec);
        }
        delete [] workspace;

//        cout << max_score << " ";
//    result.backtrace = optAlnResult.cigar_string;
    result.score = max_score;
//    result.qStartPos = optAlnResult.query_start;
//    result.qEndPos = optAlnResult.query_end;
//    result.dbEndPos = optAlnResult.target_end;
//    result.dbStartPos = optAlnResult.query_start;
    //result.qCov = SmithWaterman::computeCov(result.qStartPos1, result.qEndPos1, querySeqObj->L);
//    result.qcov = SmithWaterman::computeCov(result.qStartPos, result.qEndPos, qSeq.L);
    //result.tCov = SmithWaterman::computeCov(result.dbStartPos1, result.dbEndPos1, targetSeqObj->L);
//    result.dbcov = SmithWaterman::computeCov(result.dbStartPos, result.dbEndPos, tSeq.L);
//    result.eval = evaluer.computeEvalue(result.score, queryLen);

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

    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, par.scoreBias);

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

#pragma omp parallel
    {

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        float * query_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
//        float *mem = (float*)mem_align(ALIGN_FLOAT,6*par.maxSeqLen*4*sizeof(float));
        Sequence qSeq(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat, 0, false, par.compBiasCorrection);
        Sequence tSeq(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat, 0, false, par.compBiasCorrection);
        int * ires = new int[par.maxSeqLen];
        char * querynn  = (char*)mem_align(ALIGN_INT, par.maxSeqLen * sizeof(char) * 4 );
        char * targetnn = (char*)mem_align(ALIGN_INT, par.maxSeqLen * sizeof(char) * 4 );

        std::vector<Matcher::result_t> alignmentResult;
        PareunAlign paruenAlign(par.maxSeqLen, &subMat); // subMat called once, don't need to call it again?
        char buffer[1024+32768];
        std::string resultBuffer;
        //int gapOpen = 15; int gapExtend = 3; // 3di gap optimization
//        EvalueComputation evaluer(tdbr->getAminoAcidDBSize(), &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
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

                float *qdata = (float *) qcadbr.getData(queryId, thread_idx);
                Coordinates queryCaCords;
                memcpy(query_x, qdata, sizeof(float) * qSeq.L);
                memcpy(query_y, &qdata[qSeq.L], sizeof(float) * qSeq.L);
                memcpy(query_z, &qdata[qSeq.L+qSeq.L], sizeof(float) * qSeq.L);

                queryCaCords.x = query_x;
                queryCaCords.y = query_y;
                queryCaCords.z = query_z;

                findNearestNeighbour(querynn, queryCaCords, querySeqLen, qSeq.numSequence);

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

                    findNearestNeighbour(targetnn, targetCaCords, tSeq.L, tSeq.numSequence);
                    Matcher::result_t res = alignByNN(querynn, qSeq.numSequence, qSeq.L, targetnn, tSeq.numSequence, tSeq.L, &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids, par.gapNW, par.nnWeight);
                    unsigned int targetKey = tdbr->getDbKey(targetId);
                    res.dbKey = targetKey;
                    //Matcher::result_t res = paruenAlign.align(qSeq, tSeq, &subMat, evaluer);
                    string cigar =  paruenAlign.backtrace2cigar(res.backtrace);


                    // float t[3], u[3][3];
                    //double TMalnScore = get_score4pareun(r1, r2,  xtm, ytm, queryCaCords, targetCaCords, ires,
                    //                                     queryLen, t, u, mem);

                    //writing temp output file
                    //tempfile << queryId << "\t" << targetId << "\t" << res.score << "\t" << cigar << "\n";

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
        delete [] ires;
        free(querynn);
        free(targetnn);
    }

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
