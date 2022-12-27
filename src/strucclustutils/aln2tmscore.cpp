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
#include "TMaligner.h"
#include "Coordinate16.h"
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
        std::string resultsStr;
        TMaligner tmaln(std::max(qdbr.getMaxSeqLen() + 1,tdbr->getMaxSeqLen() + 1), false);
        Coordinate16 qcoords;
        Coordinate16 tcoords;

#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < alndbr.getSize(); i++) {
            progress.updateProgress();
            unsigned int queryKey = alndbr.getDbKey(i);
            char *data = alndbr.getData(i, thread_idx);
            Matcher::readAlignmentResults(results, data, false);

            unsigned int queryId = qdbr.getId(queryKey);
            int queryLen = static_cast<int>((qdbr.getEntryLen(queryId)-1)/(3*sizeof(float)));
            char *qcadata = qdbr.getData(queryId, thread_idx);
            size_t qCaLength = qdbr.getEntryLen(queryId);
            float* qdata = qcoords.read(qcadata, queryLen, qCaLength);
            tmaln.initQuery(qdata, &qdata[queryLen], &qdata[queryLen+queryLen], NULL, queryLen);

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
                char *tcadata = tdbr->getData(targetId, thread_idx);
                size_t tCaLength = tdbr->getEntryLen(targetId);
                float* tdata = tcoords.read(tcadata, targetLen, tCaLength);

                // Matching residue index collection
                TMaligner::TMscoreResult tmres = tmaln.computeTMscore(tdata, &tdata[targetLen], &tdata[targetLen + targetLen], targetLen, res.qStartPos,
                                     res.dbStartPos, res.backtrace);

                //std::cout << TMalnScore << std::endl;
                resultsStr.append(SSTR(dbKey));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.tmscore));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.t[0]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.t[1]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.t[2]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.u[0][0]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.u[0][1]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.u[0][2]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.u[1][0]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.u[1][1]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.u[1][2]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.u[2][0]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.u[2][1]));
                resultsStr.push_back(' ');
                resultsStr.append(SSTR(tmres.u[2][2]));
                resultsStr.push_back('\n');
            }
            dbw.writeData(resultsStr.c_str(), resultsStr.size(), queryKey, thread_idx);

            results.clear();
            resultsStr.clear();
        }
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