//
// Created by Martin Steinegger on 11/14/20.
//
#include "Command.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "DBWriter.h"
#include "Matcher.h"
#include "tmalign/TMalign.h"
#include "TMaligner.h"
#include "Coordinate16.h"

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


    DBReader<unsigned int> qSeqReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    qSeqReader.open(DBReader<unsigned int>::NOSORT);

    bool sameDB = false;
    DBReader<unsigned int> *tdbr = NULL;
    DBReader<unsigned int> *tSeqReader = NULL;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
        tSeqReader = &qSeqReader;
    } else {
        tdbr = new DBReader<unsigned int>((par.db2 + "_ca").c_str(), (par.db2 + "_ca.index").c_str(),
                                          par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tdbr->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            tdbr->readMmapedDataInMemory();
        }

        tSeqReader = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
        tSeqReader->open(DBReader<unsigned int>::NOSORT);
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
        std::string resultsStr;
        resultsStr.reserve(10 * 1024);

        TMaligner tmaln(std::max(qdbr.getMaxSeqLen() + 1,tdbr->getMaxSeqLen() + 1), false, true);
        Coordinate16 qcoords;
        Coordinate16 tcoords;

#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < alndbr.getSize(); i++) {
            progress.updateProgress();
            unsigned int queryKey = alndbr.getDbKey(i);

            unsigned int qSeqId = qSeqReader.getId(queryKey);
            int queryLen = qSeqReader.getSeqLen(qSeqId);
            unsigned int queryId = qdbr.getId(queryKey);
            char *qcadata = qdbr.getData(queryId, thread_idx);
            size_t qCaLength = qdbr.getEntryLen(queryId);
            float* qdata = qcoords.read(qcadata, queryLen, qCaLength);

            tmaln.initQuery(qdata, &qdata[queryLen], &qdata[queryLen+queryLen], NULL, queryLen);

            char *data = alndbr.getData(i, thread_idx);
            while (*data != '\0') {
                Matcher::result_t res = Matcher::parseAlignmentRecord(data, false);
                data = Util::skipLine(data);

                if (res.backtrace.empty()) {
                    Debug(Debug::ERROR) << "Backtrace cigar is missing in the alignment result. Please recompute the alignment with the -a flag.\n"
                                           "Command: foldseek structurealign " << par.db1 << " " << par.db2 << " " << par.db3 << " " << "alnNew -a\n";
                    EXIT(EXIT_FAILURE);
                }
                const unsigned int dbKey = res.dbKey;
                unsigned int tSeqId = tSeqReader->getId(dbKey);
                int targetLen = tSeqReader->getSeqLen(tSeqId);
                unsigned int targetId = tdbr->getId(dbKey);
                char *tcadata = tdbr->getData(targetId, thread_idx);
                size_t tCaLength = tdbr->getEntryLen(targetId);
                float* tdata = tcoords.read(tcadata, targetLen, tCaLength);

                // Matching residue index collection
                TMaligner::TMscoreResult tmres = tmaln.computeTMscore(tdata, &tdata[targetLen], &tdata[targetLen + targetLen], targetLen,
                                                                      res.qStartPos, res.dbStartPos, res.backtrace,
                                                                      std::min(static_cast<unsigned int>(res.backtrace.length()), std::min(res.dbLen, res.qLen)));

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

            resultsStr.clear();
        }
    }
    dbw.close();

    alndbr.close();
    qSeqReader.close();
    qdbr.close();
    if (sameDB == false) {
        tSeqReader->close();
        delete tSeqReader;
        tdbr->close();
        delete tdbr;
    }

    return EXIT_SUCCESS;
}