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
#include "LoLAlign.h"
#include "Coordinate16.h"

#ifdef OPENMP
#include <omp.h>
#endif

bool compareHitsByScore(const Matcher::result_t &first, const Matcher::result_t &second) {
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

int runStructureAligner(int argc, const char **argv, const Command& command, bool runLoLAlign) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    SubstitutionMatrix *subMat3Di = nullptr;
    SubstitutionMatrix *subMatAA  = nullptr;
    if (runLoLAlign) {
        subMat3Di = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);

        std::string blosum;
        for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
            if (par.substitutionMatrices[i].name == "blosum62.out")  {
                std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
                std::string matrixName = par.substitutionMatrices[i].name;
                char *serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
                blosum.assign(serializedMatrix);
                free(serializedMatrix);
                break;
            }
        }
        float aaFactor = (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA) ? 1.4 : 0.0;
        subMatAA = new SubstitutionMatrix(blosum.c_str(), aaFactor, par.scoreBias);
    }

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

    IndexReader *qdbr3Di = NULL;

    if (runLoLAlign) {
        qdbr3Di = new IndexReader(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    }

    IndexReader *tdbr = NULL;
    IndexReader *tcadbr = NULL;
    IndexReader *tdbr3Di = NULL;

    bool sameDB = false;
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool alignmentIsExtended = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
        tcadbr = &qcadbr;
        if (runLoLAlign) {
            tdbr3Di = qdbr3Di;
        }
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
        if (runLoLAlign) {
            tdbr3Di = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        }
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

    std::vector<Coordinate16 *> tcoords;

    std::vector<LoLAlign *> lolaligner;
    std::vector<FwBwAligner *> fwbwAligner;
    std::vector<Sequence *> tSeqAAs;
    std::vector<Sequence *> tSeq3Dis;

    std::vector<TMaligner *> tmaligner;
    tcoords.resize(par.threads);

    int max_targetLen = static_cast<int>(tdbr->sequenceReader->getMaxSeqLen() + 1);
    max_targetLen = ((max_targetLen + par.blocklen -1) / par.blocklen) * par.blocklen;
        
    if (runLoLAlign){
        lolaligner.resize(par.threads);
        fwbwAligner.resize(par.threads);
        tSeqAAs.resize(par.threads);
        tSeq3Dis.resize(par.threads);
    } else {
        tmaligner.resize(par.threads);
    }
    

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        if (runLoLAlign) {
            lolaligner[thread_idx] = new LoLAlign(std::max(qdbr.sequenceReader->getMaxSeqLen() + 1,
                                                       tdbr->sequenceReader->getMaxSeqLen() + 1), false);
            fwbwAligner[thread_idx] = new FwBwAligner(-par.fwbwGapopen, -par.fwbwGapextend, par.temperature, 0, qdbr.sequenceReader->getMaxSeqLen() + 1, tdbr->sequenceReader->getMaxSeqLen() + 1, par.blocklen, 0);                
            tSeqAAs[thread_idx] = new Sequence(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMatAA, 0, false, par.compBiasCorrection);
            tSeq3Dis[thread_idx] = new Sequence(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat3Di, 0, false, par.compBiasCorrection);
           
        } else {
            tmaligner[thread_idx] = new TMaligner(std::max(qdbr.sequenceReader->getMaxSeqLen() + 1,
                                                        tdbr->sequenceReader->getMaxSeqLen() + 1),
                                                par.tmAlignFast, false, false);
        }
        tcoords[thread_idx] = new Coordinate16();
    }


    std::vector<Matcher::result_t> swResults;
    std::vector<Matcher::result_t> finalHits;
    std::vector<unsigned int> dbKeys;
    std::string resultBuffer;

    Sequence * qSeqAA = nullptr;
    Sequence * qSeq3Di = nullptr;

    if (runLoLAlign) {
        qSeqAA = new Sequence(par.maxSeqLen, qdbr.sequenceReader->getDbtype(), (const BaseMatrix *) subMatAA, 0, false, par.compBiasCorrection);
        qSeq3Di = new Sequence(par.maxSeqLen, qdbr3Di->getDbtype(), (const BaseMatrix *) subMat3Di, 0, false, par.compBiasCorrection);
    }

    for (size_t id = 0; id < resultReader.getSize(); id++) {
        progress.updateProgress();
        swResults.clear();
        finalHits.clear();
        dbKeys.clear();
        resultBuffer.clear();
        size_t queryKey = resultReader.getDbKey(id);
        char *data = resultReader.getData(id,0);
        if (*data == '\0') {
            dbw.writeData(resultBuffer.c_str(), resultBuffer.size(), queryKey, 0);
            continue;
        }

        unsigned int queryId = qdbr.sequenceReader->getId(queryKey);
        char *querySeq = qdbr.sequenceReader->getData(queryId, 0);
        int queryLen = static_cast<int>(qdbr.sequenceReader->getSeqLen(queryId));

        if (runLoLAlign) {
            char *query3diSeq = qdbr3Di->sequenceReader->getData(queryId, 0);
            qSeqAA->mapSequence(id, queryKey, querySeq, queryLen);
            qSeq3Di->mapSequence(id, queryKey, query3diSeq, queryLen);
        }
    
        char *qcadata = qcadbr.sequenceReader->getData(queryId, 0);
        size_t qCaLength = qcadbr.sequenceReader->getEntryLen(queryId);

        Coordinate16 qcoords; // gyuri // gg
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
            if (runLoLAlign) {
                lolaligner[thread_idx]->initQuery(qdata, &qdata[queryLen], &qdata[queryLen + queryLen], *qSeqAA, *qSeq3Di, queryLen, *subMatAA, max_targetLen, par.multiDomain);
            } else {
                tmaligner[thread_idx]->initQuery(qdata, &qdata[queryLen], &qdata[queryLen + queryLen],
                                             querySeq, queryLen);
            }
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
                        if (!Util::canBeCovered(par.covThr, par.covMode, queryLen, targetLen) && runLoLAlign == false) {
                            tmpResult.eval = -1.0f; // this should avoid that the hit is added
                            tmpResult.score = -1.0f;
                            tmpResult.seqId = -1.0f;
                            tmpResult.qcov = 0.0f;
                            tmpResult.dbcov = 0.0f;
                        } else {
                            if (runLoLAlign) {
                                tSeqAAs[thread_idx]->mapSequence(targetId, dbKeys[i], targetSeq, targetLen);

                                char *target3diSeq = tdbr3Di->sequenceReader->getData(targetId, thread_idx);
                                tSeq3Dis[thread_idx]->mapSequence(targetId, dbKeys[i], target3diSeq, targetLen);
                            }
                            char *tcadata = tcadbr->sequenceReader->getData(targetId, thread_idx);
                            size_t tCaLength = tcadbr->sequenceReader->getEntryLen(targetId);

                            float *tdata = tcoords[thread_idx]->read(tcadata, targetLen, tCaLength);

                            float TMscore;


                            if (runLoLAlign) {
                                if (targetLen <= 10) {
                                    lolaligner[thread_idx]->setStartAnchorLength(1);
                                    if (targetLen < 4) {
                                        lolaligner[thread_idx]->setStartAnchorLength(0);
                                    }
                                    tmpResult = lolaligner[thread_idx]->align(dbKeys[i], tdata, &tdata[targetLen], &tdata[targetLen + targetLen], *tSeqAAs[thread_idx], *tSeq3Dis[thread_idx], targetLen, *subMatAA, fwbwAligner[thread_idx], par.multiDomain);
                                    if (queryLen > 10) {
                                        lolaligner[thread_idx]->setStartAnchorLength(3);
                                    }
                                } else {
                                    tmpResult = lolaligner[thread_idx]->align(dbKeys[i], tdata, &tdata[targetLen], &tdata[targetLen + targetLen], *tSeqAAs[thread_idx], *tSeq3Dis[thread_idx], targetLen, *subMatAA, fwbwAligner[thread_idx], par.multiDomain);
                                }

                            } else {
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
		        bool hasTMscore = (r.eval >= par.tmScoreThr);
                if(runLoLAlign){
			        hasTMscore = true;
		        }
                if (hasCov && hasSeqId && hasTMscore) {
                    finalHits.push_back(r);
                    passedNum++;
                    rejected = 0;
                } else {
                    rejected++;
                }
            }
        } // end chunk
        SORT_SERIAL(finalHits.begin(), finalHits.end(), compareHitsByScore);
        
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
        if (runLoLAlign) {
            delete lolaligner[i];
            delete fwbwAligner[i];
            delete tSeqAAs[i];
            delete tSeq3Dis[i];
        } else {
            delete tmaligner[i];
        }
        
        delete tcoords[i];
    }

    if (sameDB == false) {
        delete tdbr;
        delete tcadbr;
        if (runLoLAlign) {
            delete tdbr3Di;
        }
    } 

    if (runLoLAlign) {
        delete subMat3Di;
        delete subMatAA;

        delete qdbr3Di;

        delete qSeqAA;
        delete qSeq3Di;
    }

    return EXIT_SUCCESS;
}


int tmalign(int argc, const char **argv, const Command &command) {
    return runStructureAligner(argc, argv, command, false);
}

int lolalign(int argc, const char **argv, const Command &command) {
    return runStructureAligner(argc, argv, command, true);
}
