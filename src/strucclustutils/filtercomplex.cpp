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



unsigned int adjustAlnLen(unsigned int qcov, unsigned int tcov, int covMode) {
    switch (covMode) {
        case Parameters::COV_MODE_BIDIRECTIONAL:
            return (qcov+tcov)/2;
        case Parameters::COV_MODE_TARGET:
            return qcov;
        case Parameters::COV_MODE_QUERY:
            return tcov;
        case Parameters::COV_MODE_LENGTH_QUERY :
        case Parameters::COV_MODE_LENGTH_TARGET :
        case Parameters::COV_MODE_LENGTH_SHORTER :
            return 0;
        default:
            return 0;
    }
}

static bool hasTM(float TMThr, int covMode, double qTM, double tTM){
    switch (covMode) {
        case Parameters::COV_MODE_BIDIRECTIONAL:
            return ((qTM>= TMThr) && (tTM >= TMThr));
        case Parameters::COV_MODE_TARGET:
            return (tTM >= TMThr);
        case Parameters::COV_MODE_QUERY:
            return (qTM >= TMThr);
        case Parameters::COV_MODE_LENGTH_QUERY :
        case Parameters::COV_MODE_LENGTH_TARGET :
        case Parameters::COV_MODE_LENGTH_SHORTER :
            return true;
        default:
            return true;
    }
}

struct ComplexFilterCriteria {
    ComplexFilterCriteria() {}
    ComplexFilterCriteria(unsigned int dbKey, unsigned int qTotalAlnLen, unsigned int tTotalAlnLen, double qTM, double tTM, std::vector<double> qChainTmScores, std::vector<double> tChainTmScores) :
                        dbKey(dbKey), qTotalAlnLen(qTotalAlnLen), tTotalAlnLen(tTotalAlnLen), qTM(qTM), tTM(tTM), qChainTmScores(qChainTmScores), tChainTmScores(tChainTmScores) {}

    bool satisfyFilterCriteria(int covMode, float covThr, float TMThr) {
        const bool covOK = Util::hasCoverage(covThr, covMode, qCov, tCov);
        const bool TMOK = hasTM(TMThr, covMode, qTM, tTM);
        return (covOK && TMOK);
    }

    unsigned int dbKey;
    unsigned int qTotalAlnLen;
    unsigned int tTotalAlnLen;
    float qCov;
    float tCov;
    double qTM;
    double tTM;
    std::vector<double> qChainTmScores; //TODO
    std::vector<double> tChainTmScores; //TODO
    // std::vector<unsigned int> alignedQChainKeys; //TODO
    // std::vector<unsigned int> alignedTChainKeys; //TODO
};

//TODO
// std::vector<unsigned double> computeChainTmScore() {
//     std::vector<unsigned double> tmScores;
//     return tmScores;
// }
unsigned int getComplexResidueLength( IndexReader *qDbr, std::vector<unsigned int> &qChainKeys) {
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
    std::map<unsigned int, unsigned int> tComplexLength;
    std::map<unsigned int, unsigned int> qComplexLength;

#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
    Matcher::result_t res;
#pragma omp for schedule(dynamic, 10) nowait

        for (size_t tComplexIdx = 0; tComplexIdx < tComplexIdVec.size(); tComplexIdx++) {
            unsigned int tComplexId = tComplexIdVec[tComplexIdx];
            std::vector<unsigned int> &tChainKeys = tComplexIdToChainKeyMap[tComplexId];
            if (tChainKeys.empty()) {
                continue;
            }
            unsigned int reslen = getComplexResidueLength(tDbr, tChainKeys);
            tComplexLength[tComplexId] =reslen;
        }
        for (size_t qComplexIdx = 0; qComplexIdx < qComplexIdVec.size(); qComplexIdx++) {
            unsigned int qComplexId = qComplexIdVec[qComplexIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap[qComplexId];
            if (qChainKeys.empty()) {
                continue;
            }
            unsigned int reslen = getComplexResidueLength(qDbr, qChainKeys);
            qComplexLength[qComplexId] = reslen;
        }
        
        for (size_t queryComplexIdx = 0; queryComplexIdx < qComplexIdVec.size(); queryComplexIdx++) {
            std::map<unsigned int, ComplexFilterCriteria> localComplexCriteria;
            unsigned int qComplexId = qComplexIdVec[queryComplexIdx];
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
                
                    if (localComplexCriteria.find(assId) == localComplexCriteria.end()) {
                        localComplexCriteria[assId] = ComplexFilterCriteria();
                        localComplexCriteria[assId].dbKey = res.dbKey;
                        localComplexCriteria[assId].qTotalAlnLen = (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                        localComplexCriteria[assId].tTotalAlnLen = (std::max(res.dbStartPos, res.dbEndPos) - std::min(res.dbStartPos, res.dbEndPos) + 1);
                        localComplexCriteria[assId].qTM = retComplex.qTmScore;
                        localComplexCriteria[assId].tTM = retComplex.tTmScore;
                        // localComplexCriteria[assId].qChainTmScores = computeChainTmScore(); //TODO
                        // localComplexCriteria[assId].tChainTmScores = computeChainTmScore(); //TODO
                    } else {
                        localComplexCriteria[assId].qTotalAlnLen += (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                        localComplexCriteria[assId].tTotalAlnLen += (std::max(res.dbStartPos, res.dbEndPos) - std::min(res.dbStartPos, res.dbEndPos) + 1);
                    }
                }
            }
            std::string result;
            std::string result5;
            std::vector<unsigned int> assIdsToDelete;

            for (auto& assId_res : localComplexCriteria){
                unsigned int tComplexId = tChainKeyToComplexIdMap[assId_res.second.dbKey];
                assId_res.second.qCov = static_cast<float>(assId_res.second.qTotalAlnLen) / static_cast<float>(qComplexLength[qComplexId]);
                assId_res.second.tCov = static_cast<float>(assId_res.second.tTotalAlnLen) / static_cast<float>(tComplexLength[tComplexId]);
                if (!assId_res.second.satisfyFilterCriteria(par.covMode, par.covThr, par.filtTmThr)){
                    assIdsToDelete.push_back(assId_res.first);
                }
            }
            for (const auto& key : assIdsToDelete) {
                localComplexCriteria.erase(key);
            }
            
            std::map<unsigned int, std::vector<unsigned int>> cmplIdToBestAssId; // cmplId : [assId, alnSum]
            for (const auto& assId_res : localComplexCriteria){
                unsigned int tComplexId = tChainKeyToComplexIdMap[assId_res.second.dbKey];
                unsigned int alnlen = adjustAlnLen(assId_res.second.qTotalAlnLen, assId_res.second.tTotalAlnLen, par.covMode);
                if (cmplIdToBestAssId.find(tComplexId) == cmplIdToBestAssId.end()){
                    cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                }
                else{
                    if (alnlen > cmplIdToBestAssId[tComplexId][1]){
                        cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                    }
                }
            }

            std::vector<unsigned int> selectedAssIDs;
            for (const auto& pair : cmplIdToBestAssId){
                selectedAssIDs.push_back(pair.second[0]);
            }

            for (const auto& assId_res : localComplexCriteria){
                unsigned int tComplexId = tChainKeyToComplexIdMap[assId_res.second.dbKey];
                if (std::find(selectedAssIDs.begin(), selectedAssIDs.end(), assId_res.first) != selectedAssIDs.end()){
                    char *outpos = Itoa::u32toa_sse2(tComplexId, buffer);
                    result.append(buffer, (outpos - buffer - 1));
                    result.push_back('\n');
                    result5.append(std::to_string(assId_res.first) + "\t" +qcomplexIdToName[qComplexId] + "\t" + tcomplexIdToName[tComplexId] + "\t" + std::to_string(assId_res.second.qCov) + "\t" + std::to_string(assId_res.second.tCov) + "\t"+ std::to_string(assId_res.second.qTM)+"\t"+ std::to_string(assId_res.second.tTM)+ "\n");
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
