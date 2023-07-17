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

    Debug(Debug::INFO) << "Query database: " << par.db1 << "\n";
    Debug(Debug::INFO) << "Target database: " << par.db2 << "\n";
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader qdbr(par.db1, par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader qcadbr(
            par.db1,
            par.threads,
            IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY),
            touch ? IndexReader::PRELOAD_INDEX : 0,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
            "_ca"
    );

    IndexReader *tdbr = NULL;
    IndexReader *tcadbr = NULL;
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
        tcadbr = &qcadbr;
    } else {
        tdbr = new IndexReader(par.db2, par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        tcadbr = new IndexReader(
                par.db2,
                par.threads,
                IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY),
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                "_ca"
        );
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();

    Debug::Progress progress(resultReader.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        TMaligner tmaln(std::max(qdbr.sequenceReader->getMaxSeqLen() + 1,tdbr->sequenceReader->getMaxSeqLen() + 1), par.tmAlignFast, false);
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
                int queryLen = static_cast<int>(qdbr.sequenceReader->getSeqLen(queryId));
                char *qcadata = qcadbr.sequenceReader->getData(queryId, thread_idx);
                size_t qCaLength = qcadbr.sequenceReader->getEntryLen(queryId);
                float* qdata = qcoords.read(qcadata, queryLen, qCaLength);
                tmaln.initQuery(qdata, &qdata[queryLen], &qdata[queryLen+queryLen], querySeq, queryLen);

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

                    // align here
                    float TMscore;
                    Matcher::result_t result = tmaln.align(dbKey, tdata, &tdata[targetLen], &tdata[targetLen+targetLen], targetSeq, targetLen, TMscore);
                    float qTM = (static_cast<float>(result.score) / 100000);
                    float tTM = result.eval;
                    switch(par.tmAlignHitOrder){
                        case LocalParameters::TMALIGN_HIT_ORDER_AVG:
                            result.eval = (qTM + tTM) / 2.0;
                            break;
                        case LocalParameters::TMALIGN_HIT_ORDER_QUERY:
                            result.eval = qTM;
                            break;
                        case LocalParameters::TMALIGN_HIT_ORDER_TARGET:
                            result.eval = tTM;
                            break;
                        case LocalParameters::TMALIGN_HIT_ORDER_MIN:
                            result.eval = std::min(qTM, tTM);
                            break;
                        case LocalParameters::TMALIGN_HIT_ORDER_MAX:
                            result.eval = std::max(qTM, tTM);
                            break;
                    }
                    result.score = static_cast<int>(qTM * 100);
                    bool hasCov = Util::hasCoverage(par.covThr, par.covMode, 1.0, 1.0);
                    bool hasSeqId = result.seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    bool hasTMscore = (TMscore >= par.tmScoreThr);
                    if(hasCov && hasSeqId  && hasTMscore){
                        swResults.emplace_back(result);
                        passedNum++;
                        rejected = 0;
                    }else{
                        rejected++;
                    }
                }
                SORT_SERIAL(swResults.begin(), swResults.end(), compareHitsByTMScore);

                for(size_t i = 0; i < swResults.size(); i++){
                    size_t len = Matcher::resultToBuffer(buffer, swResults[i], par.addBacktrace, false);
                    resultBuffer.append(buffer, len);
                }

                dbw.writeData(resultBuffer.c_str(), resultBuffer.size(), queryKey, thread_idx);
                resultBuffer.clear();
                swResults.clear();
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
