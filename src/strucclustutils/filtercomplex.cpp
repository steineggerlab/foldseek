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

unsigned int getQueryResidueLength( IndexReader& qDbr, std::vector<unsigned int> &qChainKeys) {
        unsigned int qResidueLen = 0;
        for (auto qChainKey: qChainKeys) {
            size_t id = qDbr.sequenceReader->getId(qChainKey);
            // Not accessible
            if (id == NOT_AVAILABLE_CHAIN_KEY)
                return 0;
            qResidueLen += qDbr.sequenceReader->getSeqLen(id);
        }
        return qResidueLen;
}

unsigned int getTargetResidueLength( IndexReader *qDbr, std::vector<unsigned int> &qChainKeys) {
        unsigned int qResidueLen = 0;
        for (auto qChainKey: qChainKeys) {
            size_t id = qDbr->sequenceReader->getId(qChainKey);
            // Not accessible
            if (id == NOT_AVAILABLE_CHAIN_KEY)
                return 0;
            qResidueLen += qDbr->sequenceReader->getSeqLen(id);
        }
        return qResidueLen;
}

std::vector<unsigned int> selecHighestCoverage( std::map<unsigned int , std::map<unsigned int, unsigned int>> &covMap){
        std::vector<unsigned int> assIdvec;
        for (auto pair : covMap){
            assIdvec.push_back(pair.second.rbegin()->second);
        }
        return assIdvec;
}


static void getlookupInfo(
        const std::string &file,
        std::map<unsigned int, std::string> &complexIdtoName,
        std::map<unsigned int, unsigned int> &chainKeyToComplexIdLookup,
        std::map<unsigned int, std::vector<unsigned int>> &complexIdToChainKeysLookup,
        std::vector<unsigned int> &complexIdVec
) {
    if (file.length() == 0) {
        return;
    }
    MemoryMapped lookupDB(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char *data = (char *) lookupDB.getData();
    char *end = data + lookupDB.mappedSize();
    const char *entry[255];
    int prevComplexId =  -1;
    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        auto chainKey = Util::fast_atoi<int>(entry[0]);
        std::string chainName(entry[1], (entry[2] - entry[1]) - 1);
        auto complexId = Util::fast_atoi<int>(entry[2]);
        chainKeyToComplexIdLookup.emplace(chainKey, complexId);
        
        size_t lastUnderscoreIndex = chainName.find_last_of('_');
        std::string complexName = chainName.substr(0, lastUnderscoreIndex);

        if (complexId != prevComplexId) {
            complexIdToChainKeysLookup.emplace(complexId, std::vector<unsigned int>());
            complexIdVec.emplace_back(complexId);
            complexIdtoName.emplace(complexId, complexName);
            prevComplexId = complexId;
        }
        complexIdToChainKeysLookup.at(complexId).emplace_back(chainKey);
        data = Util::skipLine(data);
    }
    lookupDB.close();
}


int filtercomplex(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    int dbaccessMode = (DBReader<unsigned int>::USE_INDEX);
    std::map<unsigned int, unsigned int> qKeyToSet;
    std::map<unsigned int, unsigned int> tKeyToSet;
    char buffer[32];

    IndexReader* qDbr;
    qDbr = new IndexReader(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader* tDbr;
    if (sameDB) {
        tDbr = qDbr;
    }
    else{
        tDbr = new IndexReader(par.db2, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    }
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    size_t localThreads = 1;

#ifdef OPENMP
    //localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif
    const bool shouldCompress = (par.compressed == true);
    const int db4Type = Parameters::DBTYPE_CLUSTER_RES;
    
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), 1, shouldCompress, db4Type);
    resultWriter.open();

    const int db5Type = Parameters::DBTYPE_GENERIC_DB;
    DBWriter resultWrite5(par.db5.c_str(), par.db5Index.c_str(), 1, shouldCompress, db5Type);
    resultWrite5.open();

    std::string qLookupFile = par.db1 + ".lookup";
    std::string tLookupFile = par.db2 + ".lookup";
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));
    Matcher::result_t res;
    std::map<unsigned int, unsigned int> qChainKeyToComplexIdMap, tChainKeyToComplexIdMap;
    std::map<unsigned int, std::vector<unsigned int>> qComplexIdToChainKeyMap, tComplexIdToChainKeyMap;
    std::map<unsigned int, std::string> qcomplexIdToName, tcomplexIdToName;
    std::vector<unsigned int> qComplexIdVec, tComplexIdVec;
    getlookupInfo(qLookupFile, qcomplexIdToName,qChainKeyToComplexIdMap, qComplexIdToChainKeyMap, qComplexIdVec);
    getlookupInfo(tLookupFile, tcomplexIdToName, tChainKeyToComplexIdMap, tComplexIdToChainKeyMap, tComplexIdVec);
    qChainKeyToComplexIdMap.clear();
    Debug::Progress progress(qComplexIdVec.size());
    std::vector<ScoreComplexResult> complexResults;
    std::map<unsigned int, unsigned int> tComplexLength;
    std::map<unsigned int, unsigned int> qComplexLength;

#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
    Matcher::result_t res;
    std::vector<ScoreComplexResult> localComplexResults;
#pragma omp for schedule(dynamic, 10) nowait
        for (size_t tComplexIdx = 0; tComplexIdx < tComplexIdVec.size(); tComplexIdx++) {
            unsigned int tComplexId = tComplexIdVec[tComplexIdx];
            std::vector<unsigned int> &tChainKeys = tComplexIdToChainKeyMap[tComplexId];
            if (tChainKeys.empty()) {
                continue;
            }
            unsigned int reslen = getTargetResidueLength(tDbr, tChainKeys);
            tComplexLength[tComplexId] =reslen;
        }
        for (size_t qComplexIdx = 0; qComplexIdx < qComplexIdVec.size(); qComplexIdx++) {
            unsigned int qComplexId = qComplexIdVec[qComplexIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap[qComplexId];
            if (qChainKeys.empty()) {
                continue;
            }
            unsigned int reslen = getTargetResidueLength(qDbr, qChainKeys);
            qComplexLength[qComplexId] = reslen;
        }
        
        for (size_t queryComplexIdx = 0; queryComplexIdx < qComplexIdVec.size(); queryComplexIdx++) {
            std::map<unsigned int, unsigned int> qcovSum, tcovSum;
            std::map<unsigned int, double> qtmScores, ttmScores;
            unsigned int qComplexId = qComplexIdVec[queryComplexIdx];
            std::map<unsigned int, unsigned int> assIdTodbKey;
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap[qComplexId];
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ) {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainDbKey = alnDbr.getId(qChainKey);
                if (qChainDbKey == NOT_AVAILABLE_CHAIN_KEY) {
                    continue;
                }
                char *data = alnDbr.getData(qChainDbKey, thread_idx);
                while (*data) {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                    if (!retComplex.isValid){
                        Debug(Debug::ERROR) << "No scorecomplex result provided";
                        EXIT(EXIT_FAILURE);
                    }
                    data = Util::skipLine(data);
                    unsigned int assId = retComplex.assId;
                    if (qcovSum.find(assId) == qcovSum.end()) {
                        qcovSum[assId] = (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                        assIdTodbKey.emplace(assId, res.dbKey);
                        qtmScores.emplace(assId, retComplex.qTmScore);
                        }
                    else{
                        qcovSum[assId] += (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                    }
                    if (tcovSum.find(assId) == tcovSum.end()) {
                        tcovSum[assId] = (std::max(res.dbStartPos, res.dbEndPos) - std::min(res.dbStartPos, res.dbEndPos) + 1);
                        assIdTodbKey.emplace(assId, res.dbKey);
                        ttmScores.emplace(assId, retComplex.tTmScore);
                        }
                    else{
                        tcovSum[assId] += (std::max(res.dbStartPos, res.dbEndPos) - std::min(res.dbStartPos, res.dbEndPos) + 1);
                    }
                }
            }
            std::string result;
            std::string result5;
            std::vector<unsigned int> keysToDelete;
            for (const auto& pair : qcovSum){
                float qcov = static_cast<float>(pair.second) / static_cast<float>(qComplexLength[qComplexId]);
                float dbcov = static_cast<float>(tcovSum[pair.first]) / static_cast<float>(tComplexLength[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]]);                
                if (!checkFilterCriteria(qcov, dbcov, par.covMode, par.covThr)){
                    keysToDelete.push_back(pair.first);
                }         
            }
            for (const auto& key : keysToDelete) {
                qcovSum.erase(key);
                tcovSum.erase(key);
            }

            std::map<unsigned int , std::map<unsigned int, unsigned int>> qcompIdToassIdToalnSum, tcompIdToassIdToalnSum, avgcompIdToassIdToalnSum;
            for (const auto& pair : qcovSum){
                if (qcompIdToassIdToalnSum.find(tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]) == qcompIdToassIdToalnSum.end()){
                    qcompIdToassIdToalnSum[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]] = {{pair.second, pair.first}};
                    tcompIdToassIdToalnSum[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]] ={{ tcovSum[pair.first], pair.first}};
                    avgcompIdToassIdToalnSum[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]] = {{(pair.second+tcovSum[pair.first])/2, pair.first}};
                }
                else{
                    qcompIdToassIdToalnSum[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]][pair.second] = pair.first;
                    tcompIdToassIdToalnSum[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]][tcovSum[pair.first]] = pair.first;    
                    avgcompIdToassIdToalnSum[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]][(pair.second + tcovSum[pair.first])/2] = pair.first ;
                }
            }
            std::vector<unsigned int> selectedAssIDs;
            switch (par.covMode) {
                case Parameters::COV_MODE_BIDIRECTIONAL:
                    selectedAssIDs = selecHighestCoverage(avgcompIdToassIdToalnSum);
                    break;
                case Parameters::COV_MODE_TARGET:
                    selectedAssIDs = selecHighestCoverage(tcompIdToassIdToalnSum);
                    break;
                case Parameters::COV_MODE_QUERY:
                    selectedAssIDs = selecHighestCoverage(qcompIdToassIdToalnSum);
                    break;
                case Parameters::COV_MODE_LENGTH_QUERY :
                    break;
                case Parameters::COV_MODE_LENGTH_TARGET :
                    break;
                case Parameters::COV_MODE_LENGTH_SHORTER :
                    break;
            }
            for (const auto& pair : qcovSum){
                if (std::find(selectedAssIDs.begin(), selectedAssIDs.end(), pair.first) != selectedAssIDs.end()){
                    char *outpos = Itoa::u32toa_sse2(tChainKeyToComplexIdMap[assIdTodbKey[pair.first]], buffer);
                    result.append(buffer, (outpos - buffer - 1));
                    result.push_back('\n');
                    
                }
                if (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL) {
                            result5.append(std::to_string(pair.first) + "\t" +qcomplexIdToName[qComplexId] + "\t" + tcomplexIdToName[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]] + "\t" + std::to_string(pair.second/static_cast<float>(qComplexLength[qComplexId])) + "\t" + std::to_string(tcovSum[pair.first]/ static_cast<float>(tComplexLength[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]])) +"\t"+ std::to_string(qtmScores[pair.first])+"\t"+ std::to_string(ttmScores[pair.first])+ "\n");
                }
                else if (par.covMode == Parameters::COV_MODE_TARGET){
                    result5.append(std::to_string(pair.first) + "\t" +qcomplexIdToName[qComplexId] + "\t" + tcomplexIdToName[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]] + "\t" + std::to_string(tcovSum[pair.first]/ static_cast<float>(tComplexLength[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]])) +"\t"+ std::to_string(ttmScores[pair.first])+ "\n");
                
                }
                else if (par.covMode == Parameters::COV_MODE_QUERY) {
                    result5.append(std::to_string(pair.first) + "\t" +qcomplexIdToName[qComplexId] + "\t" + tcomplexIdToName[tChainKeyToComplexIdMap[assIdTodbKey[pair.first]]] + "\t" + std::to_string(pair.second/static_cast<float>(qComplexLength[qComplexId])) +"\t"+ std::to_string(qtmScores[pair.first])+"\n");
                }
            }

            resultWriter.writeData(result.c_str(), result.length(), qComplexId);
            resultWrite5.writeData(result5.c_str(), result5.length(), 0);
        }
    }
    resultWriter.close(true);
    resultWrite5.close(true);
    alnDbr.close();
    delete qDbr;
    if (sameDB == false) {
        delete tDbr;
    }
    return EXIT_SUCCESS;
}