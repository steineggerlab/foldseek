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

    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qdbr.open(DBReader<unsigned int>::NOSORT);



    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        qdbr.readMmapedDataInMemory();
    }

    bool sameDB = false;
    DBReader<unsigned int> *tdbr = NULL;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
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
        int * invmap = new int[par.maxSeqLen];
        std::string resultsStr;
#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < alndbr.getSize(); i++) {
            progress.updateProgress();

            unsigned int queryKey = alndbr.getDbKey(i);
            char *qSeq = NULL;
            if (par.extractMode == Parameters::EXTRACT_QUERY) {
                qSeq = qdbr.getDataByDBKey(queryKey, thread_idx);
            }

            char *data = alndbr.getData(i, thread_idx);
            Matcher::readAlignmentResults(results, data, false);
            unsigned int queryId = qdbr.getId(queryKey);
            int queryLen = static_cast<int>((qdbr.getEntryLen(queryId)-1)/(3*sizeof(float)));
            float *qdata = (float *) qdbr.getData(queryId, thread_idx);
            Coordinates queryCaCords;
            memcpy(query_x, qdata, sizeof(float) * queryLen);
            memcpy(query_y, &qdata[queryLen], sizeof(float) * queryLen);
            memcpy(query_z, &qdata[queryLen+queryLen], sizeof(float) * queryLen);
            queryCaCords.x = query_x;
            queryCaCords.y = query_y;
            queryCaCords.z = query_z;

            for (size_t j = 0; j < results.size(); j++) {
                Matcher::result_t& res = results[j];
                Coordinates targetCaCords;
                char dbKeyBuffer[255 + 1];
                const char* words[10];
                Util::parseKey(data, dbKeyBuffer);
                data = Util::skipLine(data);
                const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                unsigned int targetId = tdbr->getId(dbKey);
                int targetLen = static_cast<int>((tdbr->getEntryLen(targetId)-1)/(3*sizeof(float)));
                float * tdata = (float*) tdbr->getData(targetId, thread_idx);
                memcpy(target_x, tdata, sizeof(float) * targetLen);
                memcpy(target_y, &tdata[targetLen], sizeof(float) * targetLen);
                memcpy(target_z, &tdata[targetLen+targetLen], sizeof(float) * targetLen);

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
                double TMalnScore = get_score4pareun(r1, r2,  xtm, ytm, queryCaCords, targetCaCords, invmap,
                                                     queryLen, t, u, mem);
                resultsStr.append(SSTR(dbKey));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(TMalnScore));
                resultsStr.push_back('\n');
            }
            dbw.writeData(resultsStr.c_str(), resultsStr.size(), queryKey, thread_idx);

            results.clear();
            resultsStr.clear();
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