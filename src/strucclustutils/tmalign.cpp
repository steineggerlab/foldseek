
#include <string>
#include <vector>
#include <tmalign/TMalign.h>
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "QueryMatcher.h"
#include "LocalParameters.h"
#include "Matcher.h"

#ifdef OPENMP
#include <omp.h>
#endif

// need for sorting the results
bool compareHitsByTMScore(const Matcher::result_t &first, const Matcher::result_t &second) {
        if (first.eval != second.eval) {
            return first.eval > second.eval;
        }
        if (first.score != second.score) {
            return first.score > second.score;
        }
        if (first.dbLen != second.dbLen) {
            return first.dbLen < second.dbLen;
        }
        return first.dbKey < second.dbKey;
}


int tmalign(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);


    //Debug(Debug::INFO) << "Sequence database: " << par.db1 << "\n";
    DBReader<unsigned int> qdbr((par.db1+"_ss").c_str(), (par.db1+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    qdbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> qcadbr((par.db1+"_ca").c_str(), (par.db1+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    qcadbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *tdbr = NULL;
    DBReader<unsigned int> *tcadbr = NULL;

    bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
        tcadbr = &qcadbr;
    } else {
        tdbr = new DBReader<unsigned int>((par.db2).c_str(), (par.db2 + ".index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        tdbr->open(DBReader<unsigned int>::NOSORT);
        tcadbr = new DBReader<unsigned int>((par.db2+"_ca").c_str(), (par.db2+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        tcadbr->open(DBReader<unsigned int>::NOSORT);
        if (touch) {
            tdbr->readMmapedDataInMemory();
            tcadbr->readMmapedDataInMemory();
        }
    }

    //Debug(Debug::INFO) << "Result database: " << par.db3 << "\n";
    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    //Debug(Debug::INFO) << "Output file: " << par.db4 << "\n";
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();

    const bool I_opt = false; // flag for -I, stick to user given alignment
    const bool a_opt = false; // flag for -a, normalized by average length
    const bool u_opt = false; // flag for -u, normalized by user specified length
    const bool d_opt = false; // flag for -d, user specified d0
    const bool fast_opt = false; // flags for -fast, fTM-align algorithm
    double Lnorm_ass = 0.0;
    double  d0_scale = 0.0;


    // parsing code
    std::string atom_opt = "auto";// use C alpha atom for protein and C3' for RNA
    std::string suffix_opt = "";    // set -suffix to empty
    std::string dir_opt = "";    // set -dir to empty
    std::string dir1_opt = "";    // set -dir1 to empty
    std::string dir2_opt = "";    // set -dir2 to empty
    Debug::Progress progress(resultReader.getSize());

#pragma omp parallel
    {

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string seqM, seqxA, seqyA;// for output alignment
        // ya[0...ylen-1][0..2], in general,
        // ya is regarded as native structure
        // --> superpose xa onto ya
        char * querySecStruc  = new char[par.maxSeqLen];
        char * targetSecStruc = new char[par.maxSeqLen];
        float * query_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );

        float * target_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float *mem = (float*)mem_align(ALIGN_FLOAT,6*par.maxSeqLen*4*sizeof(float));
        AffineNeedlemanWunsch affineNW(std::max(qdbr.getMaxSeqLen() + 1,tdbr->getMaxSeqLen() + 1), 20);
        std::vector<Matcher::result_t> swResults;
        char buffer[1024+32768];
        std::string resultBuffer;
#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            char *data = resultReader.getData(id, thread_idx);
            if(*data != '\0') {
                size_t queryKey = resultReader.getDbKey(id);
                unsigned int queryId = qdbr.getId(queryKey);
                char *querySeq = qdbr.getData(queryId, thread_idx);
                int queryLen = static_cast<int>(qdbr.getSeqLen(queryId));
                float *qdata = (float *) qcadbr.getData(queryId, thread_idx);
                memset(querySecStruc, 0, sizeof(int) * queryLen);
                Coordinates queryCaCords;
                memcpy(query_x, qdata, sizeof(float) * queryLen);
                memcpy(query_y, &qdata[queryLen], sizeof(float) * queryLen);
                memcpy(query_z, &qdata[queryLen+queryLen], sizeof(float) * queryLen);
                queryCaCords.x = query_x;
                queryCaCords.y = query_y;
                queryCaCords.z = query_z;
                make_sec(queryCaCords, queryLen, querySecStruc); // secondary structure assignment

                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = tdbr->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;
                    if(isIdentity == true){
                        std::string backtrace = "";
                        int anat= (queryLen%4) ? (queryLen/4)*4+4 : queryLen;
                        for(int i=queryLen;i<anat;i++){
                            query_x[i]=0.0f;
                            query_x[i]=0.0f;
                            query_x[i]=0.0f;
                        }

//                        float tmscore= tmscore_cpu_soa_sse2(queryLen,
//                                                            query_x, query_y, query_z,
//                                                            query_x, query_y, query_z,
//                                                            0,0,0);
                        Matcher::result_t result(dbKey, 0 , 1.0, 1.0, 1.0, 1.0, std::max(queryLen,queryLen), 0, queryLen-1, queryLen, 0, queryLen-1, queryLen, backtrace);
                        size_t len = Matcher::resultToBuffer(buffer, result, true, false);
                        resultBuffer.append(buffer, len);
                        continue;
                    }
                    char * targetSeq = tdbr->getData(targetId, thread_idx);
                    int targetLen = static_cast<int>(tdbr->getSeqLen(targetId));
                    float * tdata = (float*) tcadbr->getData(targetId, thread_idx);
                    if(Util::canBeCovered(par.covThr, par.covMode, queryLen, targetLen)==false){
                        continue;
                    }
                    memset(targetSecStruc, 0, sizeof(int)*targetLen);
                    Coordinates targetCaCords;
                    memcpy(target_x, tdata, sizeof(float) * targetLen);
                    memcpy(target_y, &tdata[targetLen], sizeof(float) * targetLen);
                    memcpy(target_z, &tdata[targetLen+targetLen], sizeof(float) * targetLen);

                    targetCaCords.x = target_x;
                    targetCaCords.y = target_y;
                    targetCaCords.z = target_z;
                    make_sec(targetCaCords, targetLen, targetSecStruc); // secondary structure assignment
                    /* entry function for structure alignment */
                    float t0[3], u0[3][3];
                    float TM1, TM2;
                    float TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                    float d0_0, TM_0;
                    float d0A, d0B, d0u, d0a;
                    float d0_out = 5.0;
                    float rmsd0 = 0.0;
                    float Liden = 0;
                    int n_ali = 0;
                    int n_ali8 = 0;
                    if(queryLen <=5 || targetLen <=5){
                        continue;
                    }
                    TMalign_main(&affineNW,
                            targetCaCords, queryCaCords, targetSeq, querySeq, targetSecStruc, querySecStruc,
                            t0, u0, TM1, TM2, TM3, TM4, TM5,
                            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                            seqM, seqxA, seqyA,
                            rmsd0, Liden,  n_ali, n_ali8,
                            targetLen, queryLen, Lnorm_ass, d0_scale,
                            I_opt, a_opt, u_opt, d_opt, fast_opt, mem);
                    //std::cout << queryId << "\t" << targetId << "\t" <<  TM_0 << "\t" << TM1 << std::endl;

                    double seqId = Liden/(static_cast<double>(n_ali8));
                    //int rmsdScore = static_cast<int>(rmsd0*1000.0);
                    std::string backtrace = "";
                    Matcher::result_t result(dbKey, static_cast<int>(TM_0*100) , 1.0, 1.0, seqId, TM_0, std::max(queryLen,targetLen), 0, queryLen-1, queryLen, 0, targetLen-1, targetLen, backtrace);


                    bool hasCov = Util::hasCoverage(par.covThr, par.covMode, 1.0, 1.0);
                    bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    bool hasTMscore = (TM1 >= par.tmScoreThr);
                    if(hasCov && hasSeqId  && hasTMscore){
                        swResults.emplace_back(result);
                    }
                }
                SORT_SERIAL(swResults.begin(), swResults.end(), compareHitsByTMScore);

                for(size_t i = 0; i < swResults.size(); i++){
                    size_t len = Matcher::resultToBuffer(buffer, swResults[i], true, false);
                    resultBuffer.append(buffer, len);
                }

                dbw.writeData(resultBuffer.c_str(), resultBuffer.size(), queryKey, thread_idx);
                resultBuffer.clear();
                swResults.clear();
            }
        }

        delete [] querySecStruc;
        delete [] targetSecStruc;
        free(query_x);
        free(query_y);
        free(query_z);
        free(target_x);
        free(target_y);
        free(target_z);
        free(mem);
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


