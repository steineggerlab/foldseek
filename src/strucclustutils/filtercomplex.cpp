//
// Created by Martin Steinegger on 2/6/24.
//
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
#include "createcomplexreport.h"
#include "LDDT.h"
#include "CalcProbTP.h"
#include <map>
#ifdef OPENMP
#include <omp.h>
#endif

int filtercomplex(int argc, const char **argv, const Command &command) {
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
#pragma omp for schedule(dynamic, 10) nowait
        for (size_t queryComplexIdx = 0; queryComplexIdx < qComplexIdVec.size(); queryComplexIdx++) {
            progress.updateProgress();
            std::vector<unsigned int> assIdVec;
//            std::vector<ComplexAlignment> compAlns;
            unsigned int qComplexId = qComplexIdVec[queryComplexIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap[qComplexId];
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ) {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainDbKey = alnDbr.getId(qChainKey);
                if (qChainDbKey == NOT_AVAILABLE_CHAIN_KEY) {
                    continue;
                }

//                getComplexNameChainName(queryChainName, qCompAndChainName);
                char *data = alnDbr.getData(qChainDbKey, thread_idx);
                while (*data != '\0') {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                    if (!retComplex.isValid){
                        Debug(Debug::ERROR) << "No scorecomplex result provided";
                        EXIT(EXIT_FAILURE);
                    }
                    data = Util::skipLine(data);
                    unsigned int assId = retComplex.assId;
                    unsigned int compAlnIdx = std::find(assIdVec.begin(), assIdVec.end(), assId) - assIdVec.begin();
//                    if (compAlnIdx == compAlns.size()) {
//                        assIdVec.emplace_back(assId);
//                        compAlns.emplace_back(queryChainName, targetChainName, retComplex.qTmScore, retComplex.tTmScore, retComplex.uString, retComplex.tString, assId);
//                    } else {
//                        compAlns[compAlnIdx].qChainNames.emplace_back(queryChainName);
//                        compAlns[compAlnIdx].tChainNames.emplace_back(targetChainName);
//                    }
                } // while end
            }
            // CHECK CRIETERIA HERE
            
            // WRITE RESULT here
            std::string result;\
            // WRITE TARGET COMPLEX IDS FOR COMPLEX THAT FULLFIL THE CRITERIA
            result.push_back("1\n");
            resultWriter.writeData(result.c_str(), result.length(), qComplexId, 0, isDb, isDb);

//            for (size_t compAlnIdx = 0; compAlnIdx < compAlns.size(); compAlnIdx++) {
//                const ComplexAlignment &aln = compAlns[compAlnIdx];
//                getScoreComplexResults(localComplexResults, aln.qChainNames, aln.tChainNames, aln.qTMScore, aln.tTMScore, aln.u, aln.t, aln.assId);
//            }
        } // for end

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