#ifndef FOLDSEEK_CREATECOMPLEXREPORT_H
#define FOLDSEEK_CREATECOMPLEXREPORT_H
#include "Matcher.h"
#include "MemoryMapped.h"

const unsigned int NOT_AVAILABLE_CHAIN_KEY = 4294967295;
const double MAX_ASSIGNED_CHAIN_RATIO = 1.0;
const double TOO_SMALL_MEAN = 1.0;
const double TOO_SMALL_CV = 0.1;
const double FILTERED_OUT = 0.0;
const unsigned int UNCLUSTERED = 0;
const unsigned int CLUSTERED = 1;
const unsigned int MIN_PTS = 2;
const unsigned int SINGLE_CHAINED_COMPLEX = 1;
const float BIT_SCORE_MARGIN = 0.9;
//const float CLUSTERING_STEPS = 100.0;
//const float DEF_DIST = -1.0;
const float DEF_BIT_SCORE = -1.0;
const int UNINITIALIZED = 0;
const float LEARNING_RATE = 0.1;
const float DEFAULT_EPS = 0.1;
const unsigned int FINISH_CLUSTERING = 2;
typedef std::pair<std::string, std::string> compNameChainName_t;
typedef std::map<unsigned int, unsigned int> chainKeyToComplexId_t;
typedef std::map<unsigned int, std::vector<unsigned int>> complexIdToChainKeys_t;
typedef std::vector<unsigned int> cluster_t;
typedef std::map<std::pair<unsigned int, unsigned int>, float> distMap_t;
typedef std::string resultToWrite_t;
typedef std::string chainName_t;
typedef std::pair<unsigned int, resultToWrite_t> resultToWriteWithKey_t;

struct ScoreComplexResult {
    ScoreComplexResult() {}
    ScoreComplexResult(unsigned int assId, resultToWrite_t &resultToWrite) : assId(assId), resultToWrite(resultToWrite) {}
    unsigned int assId;
    resultToWrite_t resultToWrite;
};

static bool compareComplexResult(const ScoreComplexResult &first, const ScoreComplexResult &second) {
    if (first.assId < second.assId) {
        return true;
    }
    if (first.assId > second.assId) {
        return false;
    }
    return false;
}

struct ComplexDataHandler {
    ComplexDataHandler(bool isValid): assId(UINT_MAX), qTmScore(0.0f), tTmScore(0.0f), isValid(isValid) {}
    ComplexDataHandler(unsigned int assId, double qTmScore, double tTmScore, std::string &uString, std::string &tString, bool isValid)
            : assId(assId), qTmScore(qTmScore), tTmScore(tTmScore), uString(uString), tString(tString), isValid(isValid) {}
    unsigned int assId;
    double qTmScore;
    double tTmScore;
    std::string uString;
    std::string tString;
    bool isValid;
};


static void getKeyToIdMapIdToKeysMapIdVec(
        const std::string &file,
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
        auto complexId = Util::fast_atoi<int>(entry[2]);
        chainKeyToComplexIdLookup.emplace(chainKey, complexId);
        if (complexId != prevComplexId) {
            complexIdToChainKeysLookup.emplace(complexId, std::vector<unsigned int>());
            complexIdVec.emplace_back(complexId);
            prevComplexId = complexId;
        }
        complexIdToChainKeysLookup.at(complexId).emplace_back(chainKey);
        data = Util::skipLine(data);
    }
    lookupDB.close();
}

static ComplexDataHandler parseScoreComplexResult(const char *data, Matcher::result_t &res) {
    const char *entry[255];
    size_t columns = Util::getWordsOfLine(data, entry, 255);
    if (columns!=16)
        return {false};
    char key[255];
    ptrdiff_t keySize =  (entry[1] - data);
    strncpy(key, data, keySize);
    key[keySize] = '\0';
    auto dbKey = Util::fast_atoi<unsigned int>(key);
    int score = Util::fast_atoi<int>(entry[1]);
    float seqId = strtof(entry[2],NULL);
    double eval = strtod(entry[3],NULL);
    int qStartPos =  Util::fast_atoi<int>(entry[4]);
    int qEndPos = Util::fast_atoi<int>(entry[5]);
    int qLen = Util::fast_atoi<int>(entry[6]);
    int dbStartPos = Util::fast_atoi<int>(entry[7]);
    int dbEndPos = Util::fast_atoi<int>(entry[8]);
    int dbLen = Util::fast_atoi<int>(entry[9]);
    auto backtrace = std::string(entry[10], entry[11] - entry[10]);
    float qCov = SmithWaterman::computeCov(qStartPos==-1 ? 0 : qStartPos, qEndPos, qLen);
    float dbCov = SmithWaterman::computeCov(dbStartPos==-1 ? 0 : dbStartPos, dbEndPos, dbLen);
    size_t alnLength = Matcher::computeAlnLength(qStartPos==-1 ? 0 : qStartPos, qEndPos, dbStartPos==-1 ? 0 : dbStartPos, dbEndPos);
    double qTmScore = std::stod(entry[11]);
    double tTmScore = std::stod(entry[12]);
    std::string uString = std::string(entry[13], entry[14] - entry[13]-1);
    std::string tString = std::string(entry[14], entry[15] - entry[14]-1);
    auto assId = Util::fast_atoi<unsigned int>(entry[15]);
    res = Matcher::result_t(dbKey, score, qCov, dbCov, seqId, eval, alnLength, qStartPos, qEndPos, qLen, dbStartPos, dbEndPos, dbLen, -1, -1, -1, -1, backtrace);
    return {assId, qTmScore, tTmScore, uString, tString, true};
}

#endif //FOLDSEEK_CREATECOMPLEXREPORT_H