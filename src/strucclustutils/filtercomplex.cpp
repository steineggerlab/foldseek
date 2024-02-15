#include "DBWriter.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "TranslateNucl.h"
#include "MemoryMapped.h"
#include "createcomplexreport.h"
#include "LDDT.h"
#include "CalcProbTP.h"
#include <map>

#ifdef OPENMP
#include <omp.h>
#endif
/*
bool checkFilterCriteria(Matcher::result_t &ressum, double seqIdThr, int alnLenThr, int covMode, float covThr) {
    const bool seqIdOK = (ressum.seqId >= seqIdThr);
    const bool covOK = Util::hasCoverage(covThr, covMode, ressum.qcov, ressum.dbcov);
    const bool alnLenOK = Util::hasAlignmentLength(alnLenThr, ressum.alnLength);
    //const bool tmOK = (re.)
    if
        // general accaptance criteria
        (
          seqIdOK    &&
          covOK      &&
          alnLenOK  
        ) {
        return true;
    } else {
        return false;
    }
}
*/
bool checkFilterCriteria(float qcov, float dbcov, int covMode, float covThr) {
    const bool covOK = Util::hasCoverage(covThr, covMode, qcov, dbcov);
    if (
          covOK
        ) {
        return true;
    } else {
        return false;
    }
}
int filtercomplex(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    int dbaccessMode = (DBReader<unsigned int>::USE_INDEX);
    std::map<unsigned int, unsigned int> qKeyToSet;
    std::map<unsigned int, unsigned int> tKeyToSet;
    IndexReader qDbr(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader tDbr(par.db2, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif
    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), 1, shouldCompress, dbType);
    resultWriter.open();
    const bool isDb = par.dbOut;
    std::string qLookupFile = par.db1 + ".lookup";
    std::string tLookupFile = par.db2 + ".lookup";
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));
    Matcher::result_t res;
    std::map<unsigned int, unsigned int> qChainKeyToComplexIdMap, tChainKeyToComplexIdMap;
    std::map<unsigned int, std::vector<unsigned int>> qComplexIdToChainKeyMap, tComplexIdToChainKeyMap;
    std::vector<unsigned int> qComplexIdVec, tComplexIdVec;
    getKeyToIdMapIdToKeysMapIdVec(qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeyMap, qComplexIdVec);
    getKeyToIdMapIdToKeysMapIdVec(tLookupFile, tChainKeyToComplexIdMap, tComplexIdToChainKeyMap, tComplexIdVec);
    qChainKeyToComplexIdMap.clear();
    Debug::Progress progress(qComplexIdVec.size());
    std::vector<ScoreComplexResult> complexResults;
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Matcher::result_t res;
        std::vector<ScoreComplexResult> localComplexResults;
        std::map<unsigned int, unsigned int> tComplexLength;
#pragma omp for schedule(dynamic, 10) nowait
        for (size_t tComplexIdx = 0; tComplexIdx < tComplexIdVec.size(); tComplexIdx++) {
            unsigned int tComplexId = tComplexIdVec[tComplexIdx];
            std::vector<unsigned int> &tChainKeys = tComplexIdToChainKeyMap[tComplexId];
            unsigned int tSeqLen = 0;
            for (size_t tChainIdx = 0; tChainIdx < tChainKeys.size(); tChainIdx++ ) {
                unsigned int tChainKey = tChainKeys[tChainIdx];
                unsigned int tChainDbKey = alnDbr.getId(tChainKey);
                if (tChainDbKey == NOT_AVAILABLE_CHAIN_KEY) {
                    continue;
                }
                const char *entry[255];
                size_t columns = Util::getWordsOfLine(tDbr.sequenceReader->getDataByDBKey(tChainKey, thread_idx), entry, 255);
                unsigned int curSeqLen = Util::fast_atoi<unsigned int>(entry[2]);
                tSeqLen += curSeqLen;
            }
            tComplexLength.emplace(tComplexId, tSeqLen);
        }
        tComplexIdToChainKeyMap.clear();

        for (size_t queryComplexIdx = 0; queryComplexIdx < qComplexIdVec.size(); queryComplexIdx++) {
            progress.updateProgress();
            std::vector<unsigned int> assIdVec;
            unsigned int qComplexId = qComplexIdVec[queryComplexIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap[qComplexId];
            std::map<unsigned int, unsigned int> covSum;
            unsigned int qSeqLen = 0;
            std::map<unsigned int, unsigned int> assIdTodbKey;
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ) {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainDbKey = alnDbr.getId(qChainKey);
                if (qChainDbKey == NOT_AVAILABLE_CHAIN_KEY) {
                    continue;
                }
                const char *entry[255];
                size_t columns = Util::getWordsOfLine(qDbr.sequenceReader->getDataByDBKey(qChainKey, thread_idx), entry, 255);
                unsigned int curSeqLen = Util::fast_atoi<unsigned int>(entry[2]);
                qSeqLen += curSeqLen;
                char *data = alnDbr.getData(qChainDbKey, thread_idx);
                while (*data != '\0') {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                    if (!retComplex.isValid){
                        Debug(Debug::ERROR) << "No scorecomplex result provided";
                        EXIT(EXIT_FAILURE);
                    }
                    data = Util::skipLine(data);
                    unsigned int assId = retComplex.assId;
//                    unsigned int compAlnIdx = std::find(assIdVec.begin(), assIdVec.end(), assId) - assIdVec.begin();
                    if (covSum.find(assId) == covSum.end()) {
                        covSum[assId] = (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                        assIdTodbKey.emplace(assId, res.dbKey);
                        }
                    else{
                    covSum[assId] += (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                    }
                } // while end
            }
            // assID : alnLen
            std::string result;
            for (const auto& pair : covSum){
                float qcov = static_cast<float>(pair.second) / static_cast<float>(qSeqLen);
                float dbcov = static_cast<float>(pair.second) / static_cast<float>(tComplexLength[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]]);
                if (checkFilterCriteria(qcov, dbcov, par.covMode, par.covThr)){
                    result += std::to_string(qComplexId)+ "\t"  + std::to_string(tChainKeyToComplexIdMap[assIdTodbKey[pair.first]])+ "\n" ;
                }
            }
            // CHECK CRIETERIA HERE
            // WRITE RESULT here
            // WRITE TARGET COMPLEX IDS FOR COMPLEX THAT FULLFIL THE CRITERIA
            resultWriter.writeData(result.c_str(), result.length(), qComplexId, 0, isDb, isDb);
        } // for end
    } // MP end
    resultWriter.close(true);
    if (isDb == false) {
        FileUtil::remove(par.db4Index.c_str());
    }
    alnDbr.close();
    return EXIT_SUCCESS;
}