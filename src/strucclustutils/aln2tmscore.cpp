//
// Created by Martin Steinegger on 11/14/20.
//
#include "Command.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "DBWriter.h"
#include "Matcher.h"
#include "structureto3di.h"
#include "SubstitutionMatrix.h"
#include "GemmiWrapper.h"
#include "tmalign/TMalign.h"
#include <iostream>
#include <dirent.h>

#ifdef OPENMP
#include <omp.h>
#endif

int aln2tmscore(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    // never allow deletions
    par.allowDeletion = false;

    DBReader<unsigned int> qdbr((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(),
                                par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qdbr.open(DBReader<unsigned int>::NOSORT);

    bool sameDB = false;
    DBReader<unsigned int> *tdbr = NULL;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
    } else {
        tdbr = new DBReader<unsigned int>((par.db2 + "_ca").c_str(), (par.db2 + "_ca.index").c_str(),
                                          par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tdbr->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            tdbr->readMmapedDataInMemory();
        }
    }

    DBReader<unsigned int> alndbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alndbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_TMSCORE);
    dbw.open();
    Debug::Progress progress(alndbr.getSize());


#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<Matcher::result_t> results;
        results.reserve(300);

        float * query_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * query_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_x = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_y = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float * target_z = (float*)mem_align(ALIGN_FLOAT, par.maxSeqLen * sizeof(float) );
        float *mem = (float*)mem_align(ALIGN_FLOAT,6*par.maxSeqLen*4*sizeof(float));
        int * invmap = new int[par.maxSeqLen];
        Coordinates targetCaCords;
        Coordinates queryCaCords;
        std::string resultsStr;
        Coordinates xtm(par.maxSeqLen); Coordinates ytm(par.maxSeqLen);
        Coordinates xt(par.maxSeqLen);
        Coordinates r1(par.maxSeqLen); Coordinates r2(par.maxSeqLen);
#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < alndbr.getSize(); i++) {
            progress.updateProgress();
            unsigned int queryKey = alndbr.getDbKey(i);
            char *data = alndbr.getData(i, thread_idx);
            Matcher::readAlignmentResults(results, data, false);

            unsigned int queryId = qdbr.getId(queryKey);
            int queryLen = static_cast<int>((qdbr.getEntryLen(queryId)-1)/(3*sizeof(float)));
            float *qdata = (float *) qdbr.getData(queryId, thread_idx);
            memcpy(query_x, qdata, sizeof(float) * queryLen);
            memcpy(query_y, &qdata[queryLen], sizeof(float) * queryLen);
            memcpy(query_z, &qdata[queryLen+queryLen], sizeof(float) * queryLen);
            queryCaCords.x = query_x;
            queryCaCords.y = query_y;
            queryCaCords.z = query_z;

            for (size_t j = 0; j < results.size(); j++) {
                Matcher::result_t& res = results[j];
                if (res.backtrace.empty()) {
                    Debug(Debug::ERROR) << "Backtrace cigar is missing in the alignment result. Please recompute the alignment with the -a flag.\n"
                                           "Command: foldseek structurealign " << par.db1 << " " << par.db2 << " " << par.db3 << " " << "alnNew -a\n";
                    EXIT(EXIT_FAILURE);
                }
                const unsigned int dbKey = res.dbKey;
                unsigned int targetId = tdbr->getId(dbKey);
                int targetLen = static_cast<int>((tdbr->getEntryLen(targetId)-1)/(3*sizeof(float)));
                float * tdata = (float*) tdbr->getData(targetId, thread_idx);
                memcpy(target_x, tdata, sizeof(float) * targetLen);
                memcpy(target_y, &tdata[targetLen], sizeof(float) * targetLen);
                memcpy(target_z, &tdata[targetLen+targetLen], sizeof(float) * targetLen);
                if(queryKey == 585 && dbKey == 3493){
//                    std::cout << "bal" << std::endl;
                }
                targetCaCords.x = target_x;
                targetCaCords.y = target_y;
                targetCaCords.z = target_z;

                // Matching residue index collection
                int qPos = res.qStartPos;
                int tPos = res.dbStartPos;
                std::string cigarString = res.backtrace;
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
                //Add TM align score

                float t[3], u[3][3];

                float D0_MIN;

                float rmsd0 = 0.0;
                int L_ali;                // Aligned length in standard_TMscore
                float Lnorm;         //normalization length
                float score_d8,d0,d0_search,dcu0;//for TMscore search
                parameter_set4search(targetLen,  queryLen, D0_MIN, Lnorm,
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
                                              invmap, t, u, 40, 8, local_d0_search, true, Lnorm, score_d8, d0, mem);
                TM = std::max(TM, TMalnScore);
                //std::cout << TMalnScore << std::endl;
                resultsStr.append(SSTR(dbKey));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(TM));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(t[0]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(t[1]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(t[2]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(u[0][0]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(u[0][1]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(u[0][2]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(u[1][0]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(u[1][1]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(u[1][2]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(u[2][0]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(u[2][1]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(u[2][2]));
                resultsStr.push_back('\n');
            }
            dbw.writeData(resultsStr.c_str(), resultsStr.size(), queryKey, thread_idx);

            results.clear();
            resultsStr.clear();
        }
        free(query_x);
        free(query_y);
        free(query_z);
        free(target_x);
        free(target_y);
        free(target_z);
        free(mem);
        delete [] invmap;
    }
    dbw.close();


    alndbr.close();
    qdbr.close();
    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }

    return EXIT_SUCCESS;
}