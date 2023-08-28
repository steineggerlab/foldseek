#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "TranslateNucl.h"
#include "MemoryMapped.h"
#include "Coordinate16.h"
#include "createcomplexreport.h"
#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>
#include "LDDT.h"
#include "CalcProbTP.h"
#include <map>
#include <map>

#ifdef OPENMP
#include <omp.h>
#endif

typedef std::pair<std::string, std::string> compNameChainName_t;

void getComplexNameChainName(std::string &chainName, compNameChainName_t &compAndChainName) {
    size_t pos = chainName.rfind('_');
    std::string comp = chainName.substr(0, pos);
    std::string chain = chainName.substr(pos + 1);
    compAndChainName = {comp, chain};
}

void getResult (std::vector<ComplexResult> &complexResVec, std::vector<std::string> &qChainVector, std::vector<std::string> &tChainVector,  float qTMScore, float tTMScore, std::string t, std::string u, int assId) {
    char buffer[1024];
    std::string result;
    std::string qComplexName;
    std::string tComplexName;
    std::string qChainString;
    std::string tChainString;
    compNameChainName_t compAndChainName;
    getComplexNameChainName(qChainVector[0], compAndChainName);
    qComplexName = compAndChainName.first;
    qChainString = compAndChainName.second;
    getComplexNameChainName(tChainVector[0], compAndChainName);
    tComplexName = compAndChainName.first;
    tChainString = compAndChainName.second;
    for (size_t qChainId = 1; qChainId < qChainVector.size(); qChainId++) {
        getComplexNameChainName(qChainVector[qChainId], compAndChainName);
        qChainString += ',' + compAndChainName.second;
    }
    for (size_t tChainId = 1; tChainId < tChainVector.size(); tChainId++) {
        getComplexNameChainName(tChainVector[tChainId], compAndChainName);
        tChainString += ',' + compAndChainName.second;
    }
    int count = snprintf(buffer,sizeof(buffer),"%s\t%s\t%s\t%s\t%1.5f\t%1.5f\t%s\t%s\t%d\n", qComplexName.c_str(), tComplexName.c_str(), qChainString.c_str(), tChainString.c_str(), qTMScore, tTMScore, t.c_str(), u.c_str(), assId);
    result.append(buffer, count);
    complexResVec.emplace_back(ComplexResult(assId, result));
}

struct ComplexAlignment {
    ComplexAlignment(){};
    ComplexAlignment(std::string qChain, std::string tChain, double qTMscore, double tTMscore, std::string t, std::string u, unsigned int assId) : qTMScore(qTMscore), tTMScore(tTMscore), t(t), u(u), assId(assId){
        qChainVector = {qChain};
        tChainVector = {tChain};
    };
    std::vector<std::string> qChainVector;
    std::vector<std::string> tChainVector;
    double qTMScore;
    double tTMScore;
    std::string t;
    std::string u;
    unsigned int assId;
};

int createcomplexreport(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    int dbaccessMode = (DBReader<unsigned int>::USE_INDEX);
    std::map<unsigned int, unsigned int> qKeyToSet;
    std::map<unsigned int, unsigned int> tKeyToSet;
    IndexReader qDbr(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader qDbrHeader(par.db1, par.threads, IndexReader::SRC_HEADERS , (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    IndexReader *tDbrHeader;
    if(sameDB){
        tDbrHeader = &qDbrHeader;
    } else{
        tDbrHeader = new IndexReader(par.db2, par.threads, IndexReader::SRC_HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif

    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), localThreads, shouldCompress, dbType);
    resultWriter.open();
    const bool isDb = par.dbOut;
    std::string qLookupFile = par.db1 + ".lookup";
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));

    Matcher::result_t res;
    std::map<unsigned int, unsigned int> qChainKeyToComplexIdMap;
    std::map<unsigned int, std::vector<unsigned int>> qComplexIdToChainKeyMap;
    std::vector<unsigned int> qComplexIdVec;
    getKeyToIdMapIdToKeysMapIdVec(qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeyMap, qComplexIdVec);
    qChainKeyToComplexIdMap.clear();
    Debug::Progress progress(qComplexIdVec.size());

#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Matcher::result_t res;
        std::vector<ComplexResult> complexResVec;

#pragma omp  for schedule(dynamic, 10)
        for (size_t queryIdx = 0; queryIdx < qComplexIdVec.size(); queryIdx++) {
            progress.updateProgress();
            std::vector<unsigned int> assIdVec;
            std::vector<ComplexAlignment> compAlns;
            unsigned int qComplexId = qComplexIdVec[queryIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap[qComplexId];
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ) {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int queryKey = alnDbr.getId(qChainKey);
                if (queryKey == NOT_AVAILABLE_CHAIN_KEY)
                    continue;
                size_t qHeaderId = qDbrHeader.sequenceReader->getId(queryKey);
                const char *qHeader = qDbrHeader.sequenceReader->getData(qHeaderId, thread_idx);
                compNameChainName_t qCompAndChainName;
                std::string queryId = Util::parseFastaHeader(qHeader);
                getComplexNameChainName(queryId, qCompAndChainName);
                char *data = alnDbr.getData(queryKey, thread_idx);
                while (*data != '\0') {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                    if (retComplex.isValid == false){
                        Debug(Debug::ERROR) << "No scorecomplex result provided";
                        EXIT(EXIT_FAILURE);
                    }
                    data = Util::skipLine(data);
                    size_t tHeaderId = tDbrHeader->sequenceReader->getId(res.dbKey);
                    const char *tHeader = tDbrHeader->sequenceReader->getData(tHeaderId, thread_idx);
                    std::string targetId = Util::parseFastaHeader(tHeader);
                    unsigned int assId = retComplex.assId;
                    unsigned int compAlnIdx = find(assIdVec.begin(), assIdVec.end(), assId) - assIdVec.begin();
                    if (compAlnIdx == compAlns.size()) {
                        assIdVec.emplace_back(assId);
                        compAlns.emplace_back(ComplexAlignment(queryId, targetId, retComplex.qTmScore, retComplex.tTmScore, retComplex.tString, retComplex.uString, assId));
                    } else {
                        compAlns[compAlnIdx].qChainVector.emplace_back(queryId);
                        compAlns[compAlnIdx].tChainVector.emplace_back(targetId);
                    }
                } // while end
            }
            for (size_t compAlnIdx = 0; compAlnIdx < compAlns.size(); compAlnIdx++) {
                ComplexAlignment &compAln = compAlns[compAlnIdx];
                getResult(complexResVec,compAln.qChainVector,compAln.tChainVector,compAln.qTMScore,compAln.tTMScore,compAln.t,compAln.u,compAln.assId);
            }
        } // for end
        SORT_SERIAL(complexResVec.begin(), complexResVec.end(), compareComplexResult);
        for (size_t i=0; i < complexResVec.size(); i++) {
            resultWriter.writeData(complexResVec[i].result.c_str(), complexResVec[i].result.length(), 0, localThreads - 1, false, false);
        }
    } // MP end
    resultWriter.close(true);
    if (isDb == false) {
        FileUtil::remove(par.db4Index.c_str());
    }
    alnDbr.close();
    if (sameDB == false) {
        delete tDbrHeader;
    }
    return EXIT_SUCCESS;
}
