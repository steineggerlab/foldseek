#include "DBWriter.h"
#include "DBReader.h"
#include "IndexReader.h"
#include "createcomplexreport.h"
#include "LocalParameters.h"
#include "Coordinate16.h"

// #include "Util.h"
// #include "Matcher.h"
// #include "Debug.h"
// #include "FileUtil.h"
// #include "MemoryMapped.h"
// #include "tmalign/basic_fun.h"
// #include "LDDT.h"
// #include "CalcProbTP.h"
// #include <map>

#ifdef OPENMP
#include <omp.h>
#endif

struct Complex {
    int complexId;
    std::string complexName;

    unsigned int nChain;
    std::vector<unsigned int> chainKeys;

    unsigned int complexLength;

    // Coordinate16 Coords;
}

static void getlookupInfo(
        const std::string &file,
        std::vector<Complex> &complexes,
        std::map<unsigned int, unsigned int> &chainKeyToComplexIdLookup,
        std::map<unsigned int, unsigned int> &complexIdtoIdx
        // std::map<unsigned int, std::vector<unsigned int>> &complexIdToChainKeysLookup,
        // std::map<unsigned int, std::string> &complexIdtoName,
        // std::vector<unsigned int> &complexIdVec
) {
    if (file.length() == 0) {
        return;
    }
    MemoryMapped lookupDB(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char *data = (char *) lookupDB.getData();
    char *end = data + lookupDB.mappedSize();
    const char *entry[255];
    int prevComplexId =  -1;
    int nComplex = 0;
    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        auto chainKey = Util::fast_atoi<int>(entry[0]);
        auto complexId = Util::fast_atoi<int>(entry[2]);
        chainKeyToComplexIdLookup.emplace(chainKey, complexId);

        if (complexId != prevComplexId) {
            std::string chainName(entry[1], (entry[2] - entry[1]) - 1);
            size_t lastUnderscoreIndex = chainName.find_last_of('_');
            std::string complexName = chainName.substr(0, lastUnderscoreIndex);

            Complex complex;
            complex.complexId = complexId;
            complex.nChain = 1;
            complex.complexName = complexName;
            complex.chainKeys.emplace_back(chainKey);
            complexes.emplace_back(complex);
            complexIdtoIdx.emplace(complexId, nComplex);

            prevComplexId = complexId;
            nComplex++;
        }
        else {
            complexes.back().nChain++;
            complexes.back().chainKeys.emplace_back(chainKey);
        }

        data = Util::skipLine(data);
    }
    lookupDB.close();
}

static void sumComplexLength (DBReader<unsigned int> &structDbr, std::vector<Complex> &complexes) {
    // Fill in the complex length
    for (size_t complexIdx = 0; complexIdx < complexes.size(); complexIdx++) {
        Complex &cmpl = complexes[complexIdx];
        if (cmpl.chainKeys.size() == 0) {
            continue;
        }
        unsigned int cmplId = cmpl.complexId;
        unsigned int cmplLen = 0;
        for (size_t chainIdx = 0; chainIdx < cmpl.chainKeys.size(); chainIdx++) {
            unsigned int chainKey = cmpl.chainKeys[chainIdx];
            structDbr.get(chainKey);
            cmplLen += structDbr.getSequenceLength();
        }
        cmpl.complexLength = cmplLen;
    }

}

int filtercomplex(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    int dbaccessMode = (DBReader<unsigned int>::USE_INDEX);

    IndexReader* qDbr = NULL;
    qDbr = new IndexReader(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    DBReader<unsigned int> qStructDbr((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(), 
                                par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qStructDbr.open(DBReader<unsigned int>::NOSORT);

    IndexReader* tDbr = NULL;
    DBReader<unsigned int> *tStructDbr = NULL;
    if (sameDB) {
        tDbr = qDbr;
        tStructDbr = &qStructDbr;
    }
    else{
        tDbr = new IndexReader(par.db2, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
        tStructDbr = new DBReader<unsigned int>((par.db2 + "_ca").c_str(), (par.db2 + "_ca.index").c_str(),
                                           par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tStructDbr->open(DBReader<unsigned int>::NOSORT);
    }

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX| DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    size_t localThreads = 1;

#ifdef OPENMP
localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif

    const bool shouldCompress = (par.compressed == true);
    const int db4Type = Parameters::DBTYPE_CLUSTER_RES;
    
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), 1, shouldCompress, db4Type);
    resultWriter.open();

    //TODO: remove resultWrite5 when done
    const int db5Type = Parameters::DBTYPE_GENERIC_DB;
    DBWriter resultWrite5(par.db5.c_str(), par.db5Index.c_str(), 1, shouldCompress, db5Type);
    resultWrite5.open();

    std::string qLookupFile = par.db1 + ".lookup";
    std::string tLookupFile = par.db2 + ".lookup";

    chainKeyToComplexId_t qChainKeyToComplexIdMap, tChainKeyToComplexIdMap;
    // complexIdToChainKeys_t qComplexIdToChainKeyMap, tComplexIdToChainKeyMap;
    // std::map<unsigned int, std::string> qcomplexIdToName, tcomplexIdToName;
    // std::vector<unsigned int> qComplexIdVec, tComplexIdVec;
    std::vector<Complex> qComplexes, tComplexes;
    std::map<unsigned int, unsigned int> qComplexIdtoIdx, tComplexIdtoIdx;
    // getlookupInfo(qLookupFile, qcomplexIdToName,qChainKeyToComplexIdMap, qComplexIdToChainKeyMap, qComplexIdVec);
    // std::map<unsigned int, unsigned int> qComplexLength, tComplexLength;
    // std::map<unsigned int, std::string> qComplexIdResult;
    getlookupInfo(qLookupFile, qComplexes, qChainKeyToComplexIdMap, qComplexIdtoIdx);

    // Fill in the complex length
    sumComplexLength(qStructDbr, qComplexes);

    if (sameDB) {
        tComplexes = qComplexes;
        tChainKeyToComplexIdMap = qChainKeyToComplexIdMap;
        tComplexIdtoIdx = qComplexIdtoIdx;
    }
    else {
        getlookupInfo(tLookupFile, tComplexes, tChainKeyToComplexIdMap, tComplexIdtoIdx);
        sumComplexLength(tStructDbr, tComplexes);
    }
    
    return EXIT_SUCCESS;
}
