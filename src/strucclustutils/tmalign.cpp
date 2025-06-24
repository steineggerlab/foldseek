#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "Debug.h"
#include "Util.h"
#include "QueryMatcher.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "StructureUtil.h"
#include "StructureSmithWaterman.h"
#include "TMaligner.h"
#include "Coordinate16.h"

#ifdef OPENMP
#include <omp.h>
#endif

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

    Debug(Debug::INFO) << "Query database: " << par.db1 << "\n";
    Debug(Debug::INFO) << "Target database: " << par.db2 << "\n";
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader qdbr(par.db1, par.threads, IndexReader::SEQUENCES,
                     touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader qcadbr(
            par.db1,
            par.threads,
            IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
            touch ? IndexReader::PRELOAD_INDEX : 0,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
            "_ca"
    );

    IndexReader *tdbr = NULL;
    IndexReader *tcadbr = NULL;
    bool sameDB = false;
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool alignmentIsExtended = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
        tcadbr = &qcadbr;
    } else {
        tdbr = new IndexReader(
                par.db2, par.threads,
                alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        tcadbr = new IndexReader(
                par.db2,
                par.threads,
                alignmentIsExtended
                ? IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB2)
                : IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                alignmentIsExtended ? "_seq_ca" : "_ca"
        );
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(),
                                        par.threads,
                                        DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbtype = Parameters::DBTYPE_ALIGNMENT_RES;
    if (alignmentIsExtended) {
        dbtype = DBReader<unsigned int>::setExtendedDbtype(dbtype,Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), 1, par.compressed, dbtype);
    dbw.open();

    Debug::Progress progress(resultReader.getSize());

    std::vector<TMaligner *> tmaligner;
    std::vector<Coordinate16 *> tcoords;
    tmaligner.resize(par.threads);
    tcoords.resize(par.threads);

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        tmaligner[thread_idx] = new TMaligner(std::max(qdbr.sequenceReader->getMaxSeqLen() + 1,
                                                       tdbr->sequenceReader->getMaxSeqLen() + 1),
                                              par.tmAlignFast, false, false);
        tcoords[thread_idx] = new Coordinate16();
    }


    std::vector<Matcher::result_t> swResults;
    std::vector<Matcher::result_t> finalHits;
    std::vector<unsigned int> dbKeys;
    std::string resultBuffer;

    for (size_t id = 0; id < resultReader.getSize(); id++) {
        progress.updateProgress();
        swResults.clear();
        finalHits.clear();
        dbKeys.clear();
        char *data = resultReader.getData(id,0);
        if (*data == '\0') {
            continue;
        }

        size_t queryKey = resultReader.getDbKey(id);
        unsigned int queryId = qdbr.sequenceReader->getId(queryKey);
        char *querySeq = qdbr.sequenceReader->getData(queryId, 0);
        int queryLen = static_cast<int>(qdbr.sequenceReader->getSeqLen(queryId));

        char *qcadata = qcadbr.sequenceReader->getData(queryId, 0);
        size_t qCaLength = qcadbr.sequenceReader->getEntryLen(queryId);

        Coordinate16 qcoords;
        float* qdata = qcoords.read(qcadata, queryLen, qCaLength);

        while (*data != '\0') {
            char dbKeyBuffer[256];
            Util::parseKey(data, dbKeyBuffer);
            const unsigned int dbKey = static_cast<unsigned int>(strtoul(dbKeyBuffer, NULL, 10));
            dbKeys.push_back(dbKey);
            data = Util::skipLine(data);
        }

        swResults.resize(dbKeys.size());
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            tmaligner[thread_idx]->initQuery(qdata, &qdata[queryLen], &qdata[queryLen + queryLen],
                                             querySeq, queryLen);
        }
        int passedNum = 0;
        int rejected = 0;
        size_t chunkSize = (par.maxAccept == INT_MAX &&
                            par.maxRejected == INT_MAX  ) ? dbKeys.size() : par.threads;
        for (size_t chunkStart = 0; chunkStart < dbKeys.size(); chunkStart += chunkSize) {
            if (passedNum >= par.maxAccept || rejected >= par.maxRejected) {
                break;
            }
            size_t chunkEnd = std::min(chunkStart + chunkSize, dbKeys.size());
#pragma omp parallel
            {
                unsigned int thread_idx = 0;
#ifdef OPENMP
                thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
                std::string backtrace;

#pragma omp for schedule(dynamic, 1)
                for (size_t i = chunkStart; i < chunkEnd; i++) {
                    Matcher::result_t tmpResult;
                    unsigned int targetId = tdbr->sequenceReader->getId(dbKeys[i]);
                    bool isIdentity = ((queryId == targetId)
                                       && (par.includeIdentity || sameDB));
                    if (isIdentity) {
                        backtrace.clear();
                        backtrace.append(SSTR(queryLen));
                        backtrace.append(1, 'M');
                        tmpResult = Matcher::result_t(dbKeys[i], 100, 1.0, 1.0, 1.0,
                                                      1.0, std::max(queryLen, queryLen), 0, queryLen - 1,
                                                      queryLen, 0, queryLen - 1, queryLen, backtrace);
                    } else {
                        tmpResult.dbKey = dbKeys[i];
                        char *targetSeq = tdbr->sequenceReader->getData(targetId, thread_idx);
                        int targetLen = static_cast<int>(tdbr->sequenceReader->getSeqLen(targetId));
                        if (!Util::canBeCovered(par.covThr, par.covMode, queryLen, targetLen)) {
                            tmpResult.eval = -1.0f; // this should avoid that the hit is added
                            tmpResult.score = -1.0f;
                            tmpResult.seqId = -1.0f;
                            tmpResult.qcov = 0.0f;
                            tmpResult.dbcov = 0.0f;
                        } else {
                            char *tcadata = tcadbr->sequenceReader->getData(targetId, thread_idx);
                            size_t tCaLength = tcadbr->sequenceReader->getEntryLen(targetId);

                            float *tdata = tcoords[thread_idx]->read(tcadata, targetLen, tCaLength);

                            float TMscore;
                            tmpResult = tmaligner[thread_idx]->align(dbKeys[i],
                                                                     tdata, &tdata[targetLen],
                                                                     &tdata[targetLen + targetLen],
                                                                     targetSeq, targetLen, TMscore);
                            // TM-align could not align
                            if (TMscore == std::numeric_limits<float>::min()) {
                                tmpResult.eval = -1.0f; // this should avoid that the hit is added
                                tmpResult.score = -1.0f;
                                tmpResult.seqId = -1.0f;
                                tmpResult.qcov = 0.0f;
                                tmpResult.dbcov = 0.0f;
                            } else {
                                float qTM = (float) tmpResult.score / 100000.0f;
                                float tTM = tmpResult.eval;
                                switch (par.tmAlignHitOrder) {
                                    case LocalParameters::TMALIGN_HIT_ORDER_AVG:
                                        tmpResult.eval = (qTM + tTM) / 2.0f;
                                        break;
                                    case LocalParameters::TMALIGN_HIT_ORDER_QUERY:
                                        tmpResult.eval = qTM;
                                        break;
                                    case LocalParameters::TMALIGN_HIT_ORDER_TARGET:
                                        tmpResult.eval = tTM;
                                        break;
                                    case LocalParameters::TMALIGN_HIT_ORDER_MIN:
                                        tmpResult.eval = std::min(qTM, tTM);
                                        break;
                                    case LocalParameters::TMALIGN_HIT_ORDER_MAX:
                                        tmpResult.eval = std::max(qTM, tTM);
                                        break;
                                }
                                tmpResult.score = static_cast<int>(qTM * 100.0f); // e.g. scaled
                            }
                        }
                    }

                    swResults[i] = tmpResult;
                }
            } // end parallel

            for (size_t i = chunkStart; i < chunkEnd; i++) {
                if (passedNum >= par.maxAccept || rejected >= par.maxRejected) {
                    break;
                }
                const Matcher::result_t &r = swResults[i];

                bool hasCov    = Util::hasCoverage(par.covThr, par.covMode, r.qcov, r.dbcov);
                bool hasSeqId  = (r.seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon()));
                bool hasTMscore= (r.eval >= par.tmScoreThr);

                if (hasCov && hasSeqId && hasTMscore) {
                    finalHits.push_back(r);
                    passedNum++;
                    rejected = 0;
                } else {
                    rejected++;
                }
            }
        } // end chunk

        SORT_SERIAL(finalHits.begin(), finalHits.end(), compareHitsByTMScore);
        resultBuffer.clear();

        char buffer[32768];
        for (size_t i = 0; i < finalHits.size(); i++) {
            size_t len = Matcher::resultToBuffer(buffer, finalHits[i], par.addBacktrace, false);
            resultBuffer.append(buffer, len);
        }

        dbw.writeData(resultBuffer.c_str(), resultBuffer.size(), queryKey, 0);
    }

    dbw.close();
    resultReader.close();

    for (int i = 0; i < par.threads; i++) {
        delete tmaligner[i];
        delete tcoords[i];
    }

    if (sameDB == false) {
        delete tdbr;
        delete tcadbr;
    }

    return EXIT_SUCCESS;
}
