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
typedef std::pair<unsigned int, unsigned int> Chain_t; // complexId, chainKey
typedef std::pair<Chain_t, Chain_t> ChainAln_t; // queryChain, dbChain
typedef std::vector<ChainAln_t>  SearchResult_t; // dbComplexId, vector of ChainAln_t

bool compareChainAln_tByDbComplexId(const ChainAln_t &first, const ChainAln_t &second) {
    if (first.second.first < second.second.first)
        return true;
    if (first.second.first > second.second.first)
        return false;
    return false;
}

bool compareChainAln_tByQueryChainKeyDbChainKey(const ChainAln_t &first, const ChainAln_t &second) {
    if (first.first.second < second.first.second)
        return true;
    if (first.first.second > second.first.second)
        return false;
    if (first.second.second < second.second.second)
        return true;
    if (first.second.second > second.second.second)
        return false;
    return false;
}

class ComplexExpander {
public:
    ComplexExpander(DBReader<unsigned int> &alnDbr, chainKeyToComplexId_t &dbChainKeyToComplexIdLookup, unsigned int thread_idx)
    : alnDbr(alnDbr), thread_idx(thread_idx), dbChainKeyToComplexIdLookup(dbChainKeyToComplexIdLookup){}

    void getSearchResults(unsigned int qComplexId, std::vector<unsigned int> &qChainKeys, std::vector<SearchResult_t> &searchResults) {
        for (auto qChainKey: qChainKeys) {
            unsigned int qKey = alnDbr.getId(qChainKey);
            if (qKey == NOT_AVAILABLE_CHAIN_KEY)
                continue;
            char *data = alnDbr.getData(qKey, thread_idx);
            qChain = Chain_t(qComplexId, qChainKey);
            while (*data != '\0') {
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const auto dbChainKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                const unsigned int dbComplexId = dbChainKeyToComplexIdLookup.at(dbChainKey);
                currAln =  ChainAln_t(qChain, Chain_t(dbComplexId, dbChainKey));
                currAlns.emplace_back(currAln);
                data = Util::skipLine(data);
            }
        }
        if (currAlns.empty())
            return;
        SORT_SERIAL(currAlns.begin(), currAlns.end(), compareChainAln_tByDbComplexId);
        unsigned int currDbComplexId = currAlns[0].second.first;
        searchResult = {};
        for (auto &aln: currAlns) {
            if (aln.second.first != currDbComplexId) {
                if (!searchResult.empty())
                    searchResults.emplace_back(searchResult);
                currDbComplexId = aln.second.first;
                searchResult = {};
            }
            searchResult.emplace_back(aln);
        }
        if (!searchResult.empty())
            searchResults.emplace_back(searchResult);
        currAlns.clear();
        searchResult.clear();
    }

private:
    DBReader<unsigned int> &alnDbr;
    unsigned int thread_idx;
    Chain_t qChain;
    ChainAln_t currAln;
    std::vector<ChainAln_t> currAlns;
    SearchResult_t searchResult;
    chainKeyToComplexId_t &dbChainKeyToComplexIdLookup;
};

int expandcomplex(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    std::string qLookupFile = par.db1 + ".lookup";
    std::string dbLookupFile = par.db2 + ".lookup";
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();
    std::vector<unsigned int> qComplexIndices;
    std::vector<unsigned int> dbComplexIndices;
    chainKeyToComplexId_t qChainKeyToComplexIdMap;
    chainKeyToComplexId_t dbChainKeyToComplexIdMap;
    complexIdToChainKeys_t dbComplexIdToChainKeysMap;
    complexIdToChainKeys_t qComplexIdToChainKeysMap;
    getKeyToIdMapIdToKeysMapIdVec(qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIndices);
    getKeyToIdMapIdToKeysMapIdVec(dbLookupFile, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, dbComplexIndices);
    dbComplexIndices.clear();
    qChainKeyToComplexIdMap.clear();
    dbComplexIdToChainKeysMap.clear();
    Debug::Progress progress(qComplexIndices.size());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<SearchResult_t> searchResults;
        resultToWrite_t result;
        ComplexExpander complexExpander(alnDbr, dbChainKeyToComplexIdMap, thread_idx);
#pragma omp for schedule(dynamic, 1)
        // for each q complex
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
            unsigned int qComplexId = qComplexIndices[qCompIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
            complexExpander.getSearchResults(qComplexId, qChainKeys, searchResults);
            for (unsigned int searchResultIdx = 0; searchResultIdx < searchResults.size(); searchResultIdx++) {
                auto &currentSearchResult = searchResults[searchResultIdx];
                SORT_SERIAL(currentSearchResult.begin(), currentSearchResult.end(),compareChainAln_tByQueryChainKeyDbChainKey);
                unsigned int qPrevChainKey = currentSearchResult[0].first.second;
                for (size_t alnIdx=0; alnIdx<currentSearchResult.size(); alnIdx++) {
                    auto &aln = currentSearchResult[alnIdx];
                    if (aln.first.second != qPrevChainKey) {
                        resultWriter.writeData(result.c_str(),result.length(),qPrevChainKey,thread_idx);
                        result.clear();
                        qPrevChainKey = aln.first.second;
                    }
                    result.append(SSTR(aln.second.second));
                    result.push_back('\n');
                }
                resultWriter.writeData(result.c_str(),result.length(),qPrevChainKey,thread_idx);
                result.clear();
            }
            searchResults.clear();
            progress.updateProgress();
        }
    }
    qComplexIndices.clear();
    dbChainKeyToComplexIdMap.clear();
    qComplexIdToChainKeysMap.clear();
    alnDbr.close();
    resultWriter.close(true);
    return EXIT_SUCCESS;
}