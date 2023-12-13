#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "MemoryMapped.h"
#include "createcomplexreport.h"

#ifdef OPENMP
#include <omp.h>
#endif

typedef std::pair<unsigned int, unsigned int> ChainKeyPair_t; // queryChain, dbChain

bool compareChainKeyPair_t(const ChainKeyPair_t &first, const ChainKeyPair_t &second) {
    if (first.first < second.first)
        return true;
    if (first.first > second.first)
        return false;
    if (first.second < second.second)
        return true;
    if (first.second > second.second)
        return false;
    return false;
}

class ComplexExpander {
public:
    ComplexExpander(DBReader<unsigned int> &alnDbr, chainKeyToComplexId_t &dbChainKeyToComplexIdMap, complexIdToChainKeys_t &dbComplexIdToChainKeysMap, unsigned int thread_idx)
    : alnDbr(alnDbr), thread_idx(thread_idx), dbChainKeyToComplexIdMap(dbChainKeyToComplexIdMap), dbComplexIdToChainKeysMap(dbComplexIdToChainKeysMap) {}

    void getQueryDbKeyPairs(std::vector<unsigned int> &qChainKeys, std::vector<ChainKeyPair_t> & queryDbKeyPairs) {
        for (auto qChainKey: qChainKeys) {
            unsigned int qKey = alnDbr.getId(qChainKey);
            if (qKey == NOT_AVAILABLE_CHAIN_KEY)
                continue;
            char *data = alnDbr.getData(qKey, thread_idx);
            while (*data != '\0') {
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const auto dbChainKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                const unsigned int dbComplexId = dbChainKeyToComplexIdMap.at(dbChainKey);
                if (std::find(dbFoundComplexIndexes.begin(), dbFoundComplexIndexes.end(), dbComplexId) == dbFoundComplexIndexes.end())
                    dbFoundComplexIndexes.emplace_back(dbComplexId);
                data = Util::skipLine(data);
            }
        }
        if (dbFoundComplexIndexes.empty())
            return;
        for (auto dbComplexId: dbFoundComplexIndexes) {
            auto &dbChainKeys = dbComplexIdToChainKeysMap.at(dbComplexId);
            for (auto qChainKey: qChainKeys) {
                for (auto dbChainKey: dbChainKeys) {
                    queryDbKeyPairs.emplace_back(qChainKey, dbChainKey);
                }
            }
        }
        dbFoundComplexIndexes.clear();
    }

private:
    DBReader<unsigned int> &alnDbr;
    unsigned int thread_idx;
    std::vector<unsigned int> dbFoundComplexIndexes;
    chainKeyToComplexId_t &dbChainKeyToComplexIdMap;
    complexIdToChainKeys_t &dbComplexIdToChainKeysMap;
};

int expandcomplex(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    std::string qLookupFile = par.db1 + ".lookup";
    std::string dbLookupFile = par.db2 + ".lookup";
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    // TEMP
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_PREFILTER_RES);
    resultWriter.open();
    std::vector<unsigned int> qComplexIndices;
    std::vector<unsigned int> dbComplexIndices;
    std::vector<ChainKeyPair_t>  chainKeyPairs;
    chainKeyToComplexId_t qChainKeyToComplexIdMap;
    chainKeyToComplexId_t dbChainKeyToComplexIdMap;
    complexIdToChainKeys_t dbComplexIdToChainKeysMap;
    complexIdToChainKeys_t qComplexIdToChainKeysMap;
    getKeyToIdMapIdToKeysMapIdVec(qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIndices);
    getKeyToIdMapIdToKeysMapIdVec(dbLookupFile, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, dbComplexIndices);
    dbComplexIndices.clear();
    qChainKeyToComplexIdMap.clear();
    Debug::Progress progress(qComplexIndices.size());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        resultToWrite_t result;
        ComplexExpander complexExpander(alnDbr, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, thread_idx);
#pragma omp for schedule(dynamic, 1)
        // for each q complex
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
            unsigned int qComplexId = qComplexIndices[qCompIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
            complexExpander.getQueryDbKeyPairs(qChainKeys, chainKeyPairs);
            SORT_SERIAL(chainKeyPairs.begin(), chainKeyPairs.end(), compareChainKeyPair_t);
            unsigned int qPrevChainKey = chainKeyPairs[0].first;
            for (size_t chainKeyPairIdx=0; chainKeyPairIdx<chainKeyPairs.size(); chainKeyPairIdx++) {
                unsigned int qCurrChainKey = chainKeyPairs[chainKeyPairIdx].first;
                if (qCurrChainKey != qPrevChainKey) {
                    resultWriter.writeData(result.c_str(),result.length(),qPrevChainKey,thread_idx);
                    result.clear();
                    qPrevChainKey = chainKeyPairs[chainKeyPairIdx].first;
                }
                result.append(SSTR(chainKeyPairs[chainKeyPairIdx].second));
                result.push_back('\n');
            }
            progress.updateProgress();
        }
    }
    qComplexIndices.clear();
    dbChainKeyToComplexIdMap.clear();
    qComplexIdToChainKeysMap.clear();
    dbComplexIdToChainKeysMap.clear();
    alnDbr.close();
    resultWriter.close(true);
    return EXIT_SUCCESS;
}