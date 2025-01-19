#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "TranslateNucl.h"
#include "MemoryMapped.h"
#include "MultimerUtil.h"
#include "LDDT.h"
#include "CalcProbTP.h"
#include <map>
#ifdef OPENMP
#include <omp.h>
#endif

void getComplexNameChainName(const chainName_t &chainName, compNameChainName_t &compAndChainName) {
    size_t pos = chainName.rfind('_');
    std::string comp = chainName.substr(0, pos);
    std::string chain = chainName.substr(pos + 1);
    compAndChainName = {comp, chain};
}

void getScoreComplexResults(
        std::vector<ScoreComplexResult> &scoreComplexResults,
        const std::vector<std::string> &qChainVector,
        const std::vector<std::string> &tChainVector,
        double qTMScore,
        double tTMScore,
        const std::string &u,
        const std::string &t,
        unsigned int qComplexId,
        unsigned int assId
) {
    char buffer[1024];
    resultToWrite_t resultToWrite;
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
    int count = snprintf(buffer,sizeof(buffer),"%s\t%s\t%s\t%s\t%1.5f\t%1.5f\t%s\t%s\t%d\n", qComplexName.c_str(), tComplexName.c_str(), qChainString.c_str(), tChainString.c_str(), qTMScore, tTMScore, u.c_str(), t.c_str(), assId);
    resultToWrite.append(buffer, count);
    scoreComplexResults.emplace_back(qComplexId, assId, resultToWrite);
}

struct ComplexAlignment {
    ComplexAlignment(){};
    ComplexAlignment(chainName_t &qChainName, chainName_t &tChainName, double qTmScore, double tTmScore, std::string &u, std::string &t,  unsigned int assId) : qTMScore(qTmScore), tTMScore(tTmScore), u(u), t(t), assId(assId){
        qChainNames = {qChainName};
        tChainNames = {tChainName};
    };
    std::vector<chainName_t> qChainNames;
    std::vector<chainName_t> tChainNames;
    double qTMScore;
    double tTMScore;
    std::string u;
    std::string t;
    unsigned int assId;
};

int createmultimerreport(int argc, const char **argv, const Command &command) {
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
    if (sameDB) {
        tDbrHeader = &qDbrHeader;
    } else {
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
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), 1, shouldCompress, dbType);
    resultWriter.open();
    const bool isDb = par.dbOut;
    std::string qLookupFile = par.db1 + ".lookup";
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));

    std::map<unsigned int, unsigned int> qChainKeyToComplexIdMap;
    std::map<unsigned int, std::vector<unsigned int>> qComplexIdToChainKeyMap;
    std::vector<unsigned int> qComplexIdVec;
    getKeyToIdMapIdToKeysMapIdVec_index(qDbr, qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeyMap, qComplexIdVec);
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
#pragma omp for schedule(dynamic, 10) nowait
        for (size_t queryComplexIdx = 0; queryComplexIdx < qComplexIdVec.size(); queryComplexIdx++) {
            progress.updateProgress();
            std::vector<unsigned int> assIdVec;
            std::vector<ComplexAlignment> compAlns;
            //
            unsigned int qComplexId = qComplexIdVec[queryComplexIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap[qComplexId];
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ) {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainDbKey = alnDbr.getId(qChainKey);
                if (qChainDbKey == NOT_AVAILABLE_CHAIN_KEY) {
                    continue;
                }
                size_t qHeaderId = qDbrHeader.sequenceReader->getId(qChainKey);
                const char *qHeader = qDbrHeader.sequenceReader->getData(qHeaderId, thread_idx);
                compNameChainName_t qCompAndChainName;
                chainName_t queryChainName = Util::parseFastaHeader(qHeader);
                getComplexNameChainName(queryChainName, qCompAndChainName);
                char *data = alnDbr.getData(qChainDbKey, thread_idx);
                while (*data != '\0') {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                    if (!retComplex.isValid){
                        Debug(Debug::ERROR) << "No scorecomplex result provided";
                        EXIT(EXIT_FAILURE);
                    }
                    data = Util::skipLine(data);
                    size_t tHeaderId = tDbrHeader->sequenceReader->getId(res.dbKey);
                    const char *tHeader = tDbrHeader->sequenceReader->getData(tHeaderId, thread_idx);
                    chainName_t targetChainName = Util::parseFastaHeader(tHeader);
                    unsigned int assId = retComplex.assId;
                    unsigned int compAlnIdx = std::find(assIdVec.begin(), assIdVec.end(), assId) - assIdVec.begin();
                    if (compAlnIdx == compAlns.size()) {
                        assIdVec.emplace_back(assId);
                        compAlns.emplace_back(queryChainName, targetChainName, retComplex.qTmScore, retComplex.tTmScore, retComplex.uString, retComplex.tString, assId);
                    } else {
                        compAlns[compAlnIdx].qChainNames.emplace_back(queryChainName);
                        compAlns[compAlnIdx].tChainNames.emplace_back(targetChainName);
                    }
                } // while end
            }
            for (size_t compAlnIdx = 0; compAlnIdx < compAlns.size(); compAlnIdx++) {
                const ComplexAlignment &aln = compAlns[compAlnIdx];
                getScoreComplexResults(localComplexResults, aln.qChainNames, aln.tChainNames, aln.qTMScore, aln.tTMScore, aln.u, aln.t, qComplexId, aln.assId);
            }
        } // for end
#pragma omp critical
        {
            complexResults.insert(complexResults.end(), localComplexResults.begin(), localComplexResults.end());
        }
    } // MP end
    SORT_PARALLEL(complexResults.begin(), complexResults.end(), compareComplexResultByQuery);
    for (size_t complexResIdx = 0; complexResIdx < complexResults.size(); complexResIdx++) {
        const ScoreComplexResult& cRes = complexResults[complexResIdx];
        const resultToWrite_t& data = cRes.resultToWrite;
        resultWriter.writeData(data.c_str(), data.length(), cRes.assId, 0, isDb, isDb);
    }
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