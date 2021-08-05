//
// Created by Stephanie Kim on 6/9/21.
//

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


int pareunaln(int argc, const char **argv, const Command& command) {
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
    ofstream tempfile;
    tempfile.open (par.db4+".temp.out");
    Debug::Progress progress(resultReader.getSize());

#pragma omp parallel
    {

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        unsigned int queryDbKey;
        float * query_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float *mem = (float*)mem_align(ALIGN_FLOAT,6*par.maxSeqLen*4*sizeof(float));
        Sequence qSeq(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat, 0, false, par.compBiasCorrection);
        Sequence tSeq(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat, 0, false, par.compBiasCorrection);

        std::vector<Matcher::result_t> alignmentResult;
        PareunAlign paruenAlign(par.maxSeqLen, &subMat); // subMat called once, don't need to call it again?
        char buffer[1024+32768];
        std::string resultBuffer;
        //int gapOpen = 15; int gapExtend = 3; // 3di gap optimization
        EvalueComputation evaluer(tdbr->getAminoAcidDBSize(), &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
        // write output file

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            char *data = resultReader.getData(id, thread_idx);
            if(*data != '\0') {
                size_t queryKey = resultReader.getDbKey(id);
                unsigned int queryId = qdbr.getId(queryKey);
                char *querySeq = qdbr.getData(queryId, thread_idx);
                unsigned int querySeqLen = qdbr.getSeqLen(id);
                qSeq.mapSequence(id, queryKey, querySeq, querySeqLen);

                int queryLen = static_cast<int>(qdbr.getSeqLen(queryId));
                float *qdata = (float *) qcadbr.getData(queryId, thread_idx);
                Coordinates queryCaCords;
                memcpy(query_x, qdata, sizeof(float) * queryLen);
                memcpy(query_y, &qdata[queryLen], sizeof(float) * queryLen);
                memcpy(query_z, &qdata[queryLen+queryLen], sizeof(float) * queryLen);

                queryCaCords.x = query_x;
                queryCaCords.y = query_y;
                queryCaCords.z = query_z;

                size_t passedNum = 0;
                unsigned int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    char dbKeyBuffer[255 + 1];
                    const char* words[10];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = tdbr->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;
                    if(isIdentity == true){
                        std::string backtrace = "";
                        float rmsdbla;
                        float bR[3][3];
                        float bt[3];
                        int anat= (queryLen%4) ? (queryLen/4)*4+4 : queryLen;
                        for(int i=queryLen;i<anat;i++){
                            query_x[i]=0.0f;
                            query_x[i]=0.0f;
                            query_x[i]=0.0f;
                        }

                        Matcher::result_t result(dbKey, 0 , 1.0, 1.0, 1.0, 1.0, std::max(queryLen,queryLen), 0, queryLen-1, queryLen, 0, queryLen-1, queryLen, backtrace);
                        size_t len = Matcher::resultToBuffer(buffer, result, true, false);
                        resultBuffer.append(buffer, len);
                        continue;
                    }
                    char * targetSeq = tdbr->getData(targetId, thread_idx);
                    unsigned int targetSeqLen = tdbr->getSeqLen(targetId);

                    int targetLen = static_cast<int>(tdbr->getSeqLen(targetId));
                    float * tdata = (float*) tcadbr->getData(targetId, thread_idx);
                    tSeq.mapSequence(targetId, dbKey, targetSeq, targetSeqLen);

                    if(Util::canBeCovered(par.covThr, par.covMode, queryLen, targetLen)==false){
                        continue;
                    }
                    Coordinates targetCaCords;
                    memcpy(target_x, tdata, sizeof(float) * targetLen);
                    memcpy(target_y, &tdata[targetLen], sizeof(float) * targetLen);
                    memcpy(target_z, &tdata[targetLen+targetLen], sizeof(float) * targetLen);

                    targetCaCords.x = target_x;
                    targetCaCords.y = target_y;
                    targetCaCords.z = target_z;

                    Matcher::result_t res = paruenAlign.align(qSeq, tSeq, &subMat, &queryCaCords, &targetCaCords, evaluer);
                    string cigar =  paruenAlign.backtrace2cigar(res.backtrace);

                    if (isIdentity) {
                        // set coverage and seqid of identity
                        res.qcov = 1.0f;
                        res.dbcov = 1.0f;
                        res.seqId = 1.0f;
                    }

                    bool hasCov = Util::hasCoverage(par.covThr, par.covMode, 1.0, 1.0);
                    bool hasSeqId = res.seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    //bool hasTMscore = (TM1 >= par.tmScoreThr);

                    //Add TM align score
                    Coordinates xtm(queryLen); Coordinates ytm(targetLen);
                    Coordinates r1(queryLen); Coordinates r2(targetLen);
                    float t[3], u[3][3];
                    float d0=0.5;

                    // Matching residue index collection
                    vector<int> alignedResidueIndex = paruenAlign.AlignedResidueIndex(res);

                    double TMalnScore = get_score4pareun(r1, r2, xtm, ytm, queryCaCords, targetCaCords, alignedResidueIndex,
                                                         d0,t, u, mem);

                    //writing temp output file
                    tempfile << queryId << "\t" << targetId << "\t" << res.score << "\t" << TMalnScore << "\t" << cigar << "\n";
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
            dbw.writeData(resultBuffer.c_str(), resultBuffer.length(), queryDbKey, thread_idx);
            resultBuffer.clear();
            alignmentResult.clear();
        }

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

int tmalign(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, 0, MMseqsParameter::COMMAND_ALIGN);


    Debug(Debug::INFO) << "Sequence database: " << par.db1 << "\n";
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

    const bool i_opt = false; // flag for -i, with user given initial alignment
    const bool I_opt = false; // flag for -I, stick to user given alignment
    const bool a_opt = false; // flag for -a, normalized by average length
    const bool u_opt = false; // flag for -u, normalized by user specified length
    const bool d_opt = false; // flag for -d, user specified d0
    const bool fast_opt = true; // flags for -fast, fTM-align algorithm
    double Lnorm_ass = 0.0;
    double  d0_scale = 0.0;


    // parsing code
    string atom_opt = "auto";// use C alpha atom for protein and C3' for RNA
    string suffix_opt = "";    // set -suffix to empty
    string dir_opt = "";    // set -dir to empty
    string dir1_opt = "";    // set -dir1 to empty
    string dir2_opt = "";    // set -dir2 to empty
    int ter_opt = 1;     // TER, END, or different chainID
    int split_opt = 2;     // split chain
    int infmt1_opt = 0;     // PDB format for chain_1
    int infmt2_opt = 0;
    int byresi_opt = 0;     // set -byresi to 0
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
        AffineNeedlemanWunsch affineNW(par.maxSeqLen, 3);

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

                size_t passedNum = 0;
                unsigned int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    char dbKeyBuffer[255 + 1];
                    const char* words[10];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = tdbr->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;
                    if(isIdentity == true){
                        std::string backtrace = "";
                        float rmsdbla;
                        float bR[3][3];
                        float bt[3];
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
                    //make_sec(targetCaCords, targetLen, targetSecStruc); // secondary structure assignment
                    /* entry function for structure alignment */
                    float t0[3], u0[3][3];
                    float TM1, TM2;
                    float TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                    float d0_0, TM_0;
                    float d0A, d0B, d0u, d0a;
                    float d0_out = 5.0;
                    float rmsd0 = 0.0;
                    int L_ali;                // Aligned length in standard_TMscore
                    float Liden = 0;
                    float TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                    int n_ali = 0;
                    int n_ali8 = 0;

                    TMalign_main(&affineNW,
                            queryCaCords, targetCaCords, querySeq, targetSeq, querySeq, targetSeq,
                            t0, u0, TM1, TM2, TM3, TM4, TM5,
                            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                            seqM, seqxA, seqyA,
                            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                            queryLen, targetLen, Lnorm_ass, d0_scale,
                            i_opt, I_opt, a_opt, u_opt, d_opt, fast_opt, mem);
                    std::cout << queryId << "\t" << targetId << "\t" <<  TM_0 << "\t" << TM1 << std::endl;

                    double seqId = Liden/(static_cast<double>(n_ali8));
                    int rmsdScore = static_cast<int>(rmsd0*1000.0);
                    std::string backtrace = "";
                    Matcher::result_t result(dbKey, static_cast<int>(TM_0*100) , 1.0, 1.0, seqId, TM1, std::max(queryLen,targetLen), 0, queryLen-1, queryLen, 0, targetLen-1, targetLen, backtrace);


                    bool hasCov = Util::hasCoverage(par.covThr, par.covMode, 1.0, 1.0);
                    bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    bool hasTMscore = (TM1 >= par.tmScoreThr);
                    if(hasCov && hasSeqId  && hasTMscore){
                        size_t len = Matcher::resultToBuffer(buffer, result, true, false);
                        resultBuffer.append(buffer, len);
                    }
                }
                dbw.writeData(resultBuffer.c_str(), resultBuffer.size(), queryKey, thread_idx);
                resultBuffer.clear();
            }
        }
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

int aln2TMscore(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    // never allow deletions
    par.allowDeletion = false;

    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qdbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> qcadbr((par.db1+"_ca").c_str(), (par.db1+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    qcadbr.open(DBReader<unsigned int>::NOSORT);

    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        qdbr.readMmapedDataInMemory();
    }

    bool sameDB = false;
    DBReader<unsigned int> *tdbr = NULL;
    DBReader<unsigned int> *tcadbr = NULL;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
        tcadbr = &qcadbr;
    } else {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tdbr->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            tdbr->readMmapedDataInMemory();
        }
    }

    DBReader<unsigned int> alndbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alndbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, tdbr->getDbtype());
    dbw.open();
    Debug::Progress progress(alndbr.getSize());

    const char newline = '\n';
    const char tab = '\t';
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<Matcher::result_t> results;
        results.reserve(300);

        unsigned int queryDbKey;
        float * query_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float *mem = (float*)mem_align(ALIGN_FLOAT,6*par.maxSeqLen*4*sizeof(float));

#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < alndbr.getSize(); i++) {
            progress.updateProgress();

            unsigned int queryKey = alndbr.getDbKey(i);
            char *qSeq = NULL;
            if (par.extractMode == Parameters::EXTRACT_QUERY) {
                qSeq = qdbr.getDataByDBKey(queryKey, thread_idx);
            }

            char *data = alndbr.getData(i, thread_idx);
            Matcher::readAlignmentResults(results, data, true);
            dbw.writeStart(thread_idx);
            for (size_t j = 0; j < results.size(); j++) {
                Matcher::result_t& res = results[j];
                size_t length = 0;
                char* seq = NULL;
                if (qSeq) {
                    seq = qSeq + res.qStartPos;
                    length = res.qEndPos - res.qStartPos+1;
                } else if (par.extractMode == Parameters::EXTRACT_TARGET) {
                    seq = tdbr->getDataByDBKey(res.dbKey, thread_idx) + res.dbStartPos;
                    length = res.dbEndPos - res.dbStartPos+1;
                } else {
                    Debug(Debug::ERROR) << "Missing extraction type!\n";
                    EXIT(EXIT_FAILURE);
                }
                unsigned int queryId = qdbr.getId(queryKey);
                int queryLen = static_cast<int>(qdbr.getSeqLen(queryId));
                float *qdata = (float *) qcadbr.getData(queryId, thread_idx);
                Coordinates queryCaCords;
                memcpy(query_x, qdata, sizeof(float) * queryLen);
                memcpy(query_y, &qdata[queryLen], sizeof(float) * queryLen);
                memcpy(query_z, &qdata[queryLen+queryLen], sizeof(float) * queryLen);

                queryCaCords.x = query_x;
                queryCaCords.y = query_y;
                queryCaCords.z = query_z;

                Coordinates targetCaCords;
                char dbKeyBuffer[255 + 1];
                const char* words[10];
                Util::parseKey(data, dbKeyBuffer);
                data = Util::skipLine(data);
                const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                unsigned int targetId = tdbr->getId(dbKey);
                int targetLen = static_cast<int>(tdbr->getSeqLen(targetId));
                float * tdata = (float*) tcadbr->getData(targetId, thread_idx);
                memcpy(target_x, tdata, sizeof(float) * targetLen);
                memcpy(target_y, &tdata[targetLen], sizeof(float) * targetLen);
                memcpy(target_z, &tdata[targetLen+targetLen], sizeof(float) * targetLen);

                targetCaCords.x = target_x;
                targetCaCords.y = target_y;
                targetCaCords.z = target_z;



                // Matching residue index collection
                int iressize = res.qStartPos + res.backtrace.size();
                vector<int> ires(iressize, -1);
                int qPos = res.qStartPos;
                int tPos = res.dbStartPos;
                string cigarString = res.backtrace;
                for (size_t btPos = 0; btPos < cigarString.size(); btPos++) {
                    if (cigarString[btPos] == 'M') {
                        ires[qPos] = tPos;
                        qPos++;
                        tPos++;
                    }
                    else if (cigarString[btPos] == 'I') {
                        tPos++;
                    }
                   else {
                        qPos++;
                    }
                }
                //Add TM align score
                Coordinates xtm(queryLen); Coordinates ytm(targetLen);
                Coordinates r1(queryLen); Coordinates r2(targetLen);
                float t[3], u[3][3];
                double TMalnScore = get_score4pareun(r1, r2, xtm, ytm, queryCaCords, targetCaCords, ires,
                                                     ires.size(), t, u, mem);
                stringstream ss;
                ss << TMalnScore;
                const char* strTMscore = ss.str().c_str();
                stringstream ss2;
                ss2 << targetId;
                const char* strtargetID = ss2.str().c_str();

                dbw.writeAdd(strtargetID, 8, thread_idx);
                dbw.writeAdd(&tab, 1, thread_idx);
                dbw.writeAdd(strTMscore, 4, thread_idx);
                dbw.writeAdd(&newline, 1, thread_idx);

            }
            dbw.writeEnd(queryKey, thread_idx);
            results.clear();
        }
    }
    dbw.close();

    if (par.extractMode == Parameters::EXTRACT_QUERY) {
        DBReader<unsigned int>::softlinkDb(par.db1, par.db4, DBFiles::SEQUENCE_ANCILLARY);
    } else {
        DBReader<unsigned int>::softlinkDb(par.db2, par.db4, DBFiles::SEQUENCE_ANCILLARY);
    }

    alndbr.close();
    qdbr.close();
    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }

    return EXIT_SUCCESS;
}