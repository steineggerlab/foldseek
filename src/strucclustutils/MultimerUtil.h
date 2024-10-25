#ifndef FOLDSEEK_MULTIMERUTIL_H
#define FOLDSEEK_MULTIMERUTIL_H
#include "Matcher.h"
#include "MemoryMapped.h"
#include "TMaligner.h"

const unsigned int NOT_AVAILABLE_CHAIN_KEY = 4294967295;
const float MAX_ASSIGNED_CHAIN_RATIO = 1.0;
const double TOO_SMALL_MEAN = 1.0;
const double TOO_SMALL_CV = 0.1;
const double FILTERED_OUT = 0.0;
const unsigned int INITIALIZED_LABEL = 0;
const unsigned int MIN_PTS = 2;
const float LEARNING_RATE = 0.1;
const float TM_SCORE_MARGIN = 0.7;
const unsigned int MULTIPLE_CHAINED_COMPLEX = 2;
const unsigned int SIZE_OF_SUPERPOSITION_VECTOR = 12;
const int SKIP_MONOMERS = 1;
typedef std::pair<std::string, std::string> compNameChainName_t;
typedef std::map<unsigned int, unsigned int> chainKeyToComplexId_t;
typedef std::map<unsigned int, std::vector<unsigned int>> complexIdToChainKeys_t;
typedef std::vector<unsigned int> cluster_t;
typedef std::string resultToWrite_t;
typedef std::string chainName_t;
typedef std::pair<unsigned int, resultToWrite_t> resultToWriteWithKey_t;

struct Chain {
    Chain() {}
    Chain(unsigned int complexId, unsigned int chainKey) : complexId(complexId), chainKey(chainKey) {}
    unsigned int complexId;
    unsigned int chainKey;
    std::vector<float> caVecX;
    std::vector<float> caVecY;
    std::vector<float> caVecZ;
};

struct ChainToChainAln {
    ChainToChainAln() {}
    ChainToChainAln(Chain &queryChain, Chain &targetChain, float *qCaData, float *dbCaData, Matcher::result_t &alnResult, TMaligner::TMscoreResult &tmResult) : qChain(queryChain), dbChain(targetChain), tmScore((float)tmResult.tmscore) {
        alnLength = alnResult.alnLength;
        matches = 0;
        unsigned int qPos = alnResult.qStartPos;
        unsigned int dbPos = alnResult.dbStartPos;
        unsigned int qXPos = 0;
        unsigned int qYPos = alnResult.qLen;
        unsigned int qZPos = alnResult.qLen * 2;
        unsigned int dbXPos = 0;
        unsigned int dbYPos = alnResult.dbLen;
        unsigned int dbZPos = alnResult.dbLen * 2;
        for (char cigar : alnResult.backtrace) {
            switch (cigar) {
                case 'M':
                    matches++;
                    qChain.caVecX.emplace_back(qCaData[qXPos + qPos]);
                    qChain.caVecY.emplace_back(qCaData[qYPos + qPos]);
                    qChain.caVecZ.emplace_back(qCaData[qZPos + qPos++]);
                    dbChain.caVecX.emplace_back(dbCaData[dbXPos + dbPos]);
                    dbChain.caVecY.emplace_back(dbCaData[dbYPos + dbPos]);
                    dbChain.caVecZ.emplace_back(dbCaData[dbZPos + dbPos++]);
                    break;
                case 'I':
                    qPos++;
                    break;
                case 'D':
                    dbPos++;
                    break;
            }
        }
        char buffer[4096];
        label = INITIALIZED_LABEL;
        superposition[0] = tmResult.u[0][0];
        superposition[1] = tmResult.u[0][1];
        superposition[2] = tmResult.u[0][2];
        superposition[3] = tmResult.u[1][0];
        superposition[4] = tmResult.u[1][1];
        superposition[5] = tmResult.u[1][2];
        superposition[6] = tmResult.u[2][0];
        superposition[7] = tmResult.u[2][1];
        superposition[8] = tmResult.u[2][2];
        superposition[9] = tmResult.t[0];
        superposition[10] = tmResult.t[1];
        superposition[11] = tmResult.t[2];
        size_t len = Matcher::resultToBuffer(buffer, alnResult, true, true, false);
        resultToWrite.append(buffer, len-1);
    }

    Chain qChain;
    Chain dbChain;
    unsigned int matches;
    unsigned int alnLength;
    resultToWrite_t resultToWrite;
    double superposition[SIZE_OF_SUPERPOSITION_VECTOR];
    unsigned int label;
    float tmScore;

    float getDistance(const ChainToChainAln &o) const {
        float dist = 0;
        for (size_t i=0; i<SIZE_OF_SUPERPOSITION_VECTOR; i++) {
            dist += std::pow(superposition[i] - o.superposition[i], 2);
        }
        return std::sqrt(dist);
    }

    void free() {
        qChain.caVecX.clear();
        qChain.caVecY.clear();
        qChain.caVecZ.clear();
        dbChain.caVecX.clear();
        dbChain.caVecY.clear();
        dbChain.caVecZ.clear();
        resultToWrite.clear();
    }
};

struct ScoreComplexResult {
    ScoreComplexResult() {}
    ScoreComplexResult(unsigned int queryComplexId, unsigned int assId, resultToWrite_t &resultToWrite) : queryComplexId(queryComplexId), assId(assId), resultToWrite(resultToWrite) {}
    unsigned int queryComplexId;
    unsigned int assId;
    resultToWrite_t resultToWrite;
};

//static bool compareComplexResult(const ScoreComplexResult &first, const ScoreComplexResult &second) {
//    if (first.assId < second.assId) {
//        return true;
//    }
//    if (first.assId > second.assId) {
//        return false;
//    }
//    return false;
//}

static bool compareComplexResultByQuery(const ScoreComplexResult &first, const ScoreComplexResult &second) {
    if (first.queryComplexId < second.queryComplexId) {
        return true;
    }
    if (first.queryComplexId > second.queryComplexId) {
        return false;
    }
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

static char* fastfloatToBuffer(float value, char* buffer) {
    if (value < 0) {
        value *= -1;
        *(buffer) = '-';
        buffer++;
    }
    int value1 = (int)(value);
    buffer = Itoa::i32toa_sse2(value1, buffer);
    *(buffer) = '.';
    buffer++;

    double value2 = value - value1;
    if (value2 < 0.1){
        *(buffer) = '0';
        buffer++;
    }
    if (value2 < 0.01){
        *(buffer) = '0';
        buffer++;
    }
    buffer = Itoa::i32toa_sse2((int)(value2 * 1000), buffer);
    return buffer;
}

#endif //FOLDSEEK_MULTIMERUTIL_H
