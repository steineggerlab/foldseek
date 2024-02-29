#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "StructureUtil.h"
#include "TMaligner.h"
#include "Coordinate16.h"
#include "createcomplexreport.h"
#include "set"

#ifdef OPENMP
#include <omp.h>
#endif

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

    float getDistance(const ChainToChainAln &o) {
        float dist = 0;
        for (size_t i=0; i<SIZE_OF_SUPERPOSITION_VECTOR; i++) {
            dist += std::pow(superposition[i] - o.superposition[i], 2);
        }
        dist = std::sqrt(dist);
        return dist;
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

// carrying chainToChainAlignments from the same query and target complex
struct SearchResult {
    SearchResult() {}
    SearchResult(std::vector<unsigned int> &chainKeys) : qChainKeys(chainKeys), alnVec({}) {}
    SearchResult(std::vector<unsigned int> &chainKeys, unsigned int qResidueLen) : qChainKeys(chainKeys), qResidueLen(qResidueLen), alnVec({}) {}
    std::vector<unsigned int> qChainKeys;
    std::vector<unsigned int> dbChainKeys;
    unsigned int qResidueLen;
    unsigned int dbResidueLen;
    std::vector<ChainToChainAln> alnVec;

    void resetDbComplex(std::vector<unsigned int> &chainKeys, unsigned int residueLen) {
        dbChainKeys = chainKeys;
        dbResidueLen = residueLen;
    }

    void standardize() {
        if (dbResidueLen == 0)
            alnVec.clear();

        if (alnVec.empty())
            return;

        double length = alnVec.size();
        double mean;
        double var;
        double sd;
        double cv;
        for (size_t i = 0; i < SIZE_OF_SUPERPOSITION_VECTOR; i++) {
            mean = 0.0;
            var = 0.0;
            for (auto &aln: alnVec) {
                mean += aln.superposition[i] / length;
            }
            for (auto &aln: alnVec) {
                var += std::pow(aln.superposition[i] - mean, 2) / length;
            }
            sd = std::sqrt(var);
            cv = (abs(mean) > TOO_SMALL_MEAN) ? sd / std::abs(mean) : sd;
            for (auto &aln: alnVec) {
                aln.superposition[i] = cv < TOO_SMALL_CV ? FILTERED_OUT : (aln.superposition[i] - mean) / sd;
            }
        }
    }
};

// compute complex tm score
// carrying final output lines
struct Assignment {
    Assignment() {}
    Assignment(unsigned int qLength, unsigned int dbLength): qResidueLength(qLength), dbResidueLength(dbLength), matches(0) {}
    unsigned int qResidueLength;
    unsigned int dbResidueLength;
    unsigned int matches;
    std::vector<float> qCaXVec;
    std::vector<float> qCaYVec;
    std::vector<float> qCaZVec;
    std::vector<float> dbCaXVec;
    std::vector<float> dbCaYVec;
    std::vector<float> dbCaZVec;
    double qTmScore;
    double dbTmScore;
    std::string tString;
    std::string uString;
    std::string backtrace;
    std::string assignmentInfo;
    std::vector<resultToWriteWithKey_t> resultToWriteLines;
    TMaligner::TMscoreResult tmResult;

    void appendChainToChainAln(ChainToChainAln &aln) {
        matches += aln.matches;
        qCaXVec.insert(qCaXVec.end(), aln.qChain.caVecX.begin(), aln.qChain.caVecX.end());
        qCaYVec.insert(qCaYVec.end(), aln.qChain.caVecY.begin(), aln.qChain.caVecY.end());
        qCaZVec.insert(qCaZVec.end(), aln.qChain.caVecZ.begin(), aln.qChain.caVecZ.end());
        dbCaXVec.insert(dbCaXVec.end(), aln.dbChain.caVecX.begin(), aln.dbChain.caVecX.end());
        dbCaYVec.insert(dbCaYVec.end(), aln.dbChain.caVecY.begin(), aln.dbChain.caVecY.end());
        dbCaZVec.insert(dbCaZVec.end(), aln.dbChain.caVecZ.begin(), aln.dbChain.caVecZ.end());
        resultToWriteLines.emplace_back(aln.qChain.chainKey, aln.resultToWrite);
    }

    void reset() {
        matches = 0;
        qCaXVec.clear();
        qCaYVec.clear();
        qCaZVec.clear();
        dbCaXVec.clear();
        dbCaYVec.clear();
        dbCaZVec.clear();
        resultToWriteLines.clear();
        uString.clear();
        tString.clear();
    }

    void getTmScore(TMaligner &tmAligner) {
        backtrace = std::string(matches, 'M');
        unsigned int normLen = std::min(qResidueLength, dbResidueLength);
        tmAligner.initQuery(&qCaXVec[0], &qCaYVec[0], &qCaZVec[0], NULL, matches);
        tmResult = tmAligner.computeTMscore(&dbCaXVec[0], &dbCaYVec[0], &dbCaZVec[0],matches,0,0, backtrace,normLen);
        qTmScore = tmResult.tmscore * normLen / qResidueLength;
        dbTmScore = tmResult.tmscore * normLen / dbResidueLength;
        qCaXVec.clear();
        qCaYVec.clear();
        qCaZVec.clear();
        dbCaXVec.clear();
        dbCaYVec.clear();
        dbCaZVec.clear();
        backtrace.clear();
    }

    void updateResultToWriteLines() {
        char sep = ',';
        char tab = '\t';
        tString.append(std::to_string(tmResult.t[0]) + sep + std::to_string(tmResult.t[1]) + sep + std::to_string(tmResult.t[2]));
        uString.append(std::to_string(tmResult.u[0][0]) + sep + std::to_string(tmResult.u[0][1]) + sep + std::to_string(tmResult.u[0][2]) + sep);
        uString.append(std::to_string(tmResult.u[1][0]) + sep + std::to_string(tmResult.u[1][1]) + sep + std::to_string(tmResult.u[1][2]) + sep);
        uString.append(std::to_string(tmResult.u[2][0]) + sep + std::to_string(tmResult.u[2][1]) + sep + std::to_string(tmResult.u[2][2]));
        assignmentInfo = tab + std::to_string(qTmScore) + tab + std::to_string(dbTmScore) + tab + uString + tab + tString;

        for (auto &resultToWrite: resultToWriteLines) {
            resultToWrite.second.append(assignmentInfo);
        }

        uString.clear();
        tString.clear();
        assignmentInfo.clear();
    }
};

struct NeighborsWithDist {
    NeighborsWithDist(unsigned int neighbor, float dist) : neighbor(neighbor), dist(dist) {}
    unsigned int neighbor;
    float dist;
};

bool compareChainToChainAlnByDbComplexId(const ChainToChainAln &first, const ChainToChainAln &second) {
    if (first.dbChain.complexId < second.dbChain.complexId)
        return true;
    if (first.dbChain.complexId > second.dbChain.complexId)
        return false;
    return false;
}

bool compareAssignment(const Assignment &first, const Assignment &second) {
    if (first.qTmScore > second.qTmScore)
        return true;
    if (first.qTmScore < second.qTmScore)
        return false;
    if (first.dbTmScore > second.dbTmScore)
        return true;
    if (first.dbTmScore < second.dbTmScore)
        return false;
    return false;
}

bool compareNeighborWithDist(const NeighborsWithDist &first, const NeighborsWithDist &second) {
    if (first.dist < second.dist)
        return true;
    if (first.dist > second.dist)
        return false;
    return false;
}

class DBSCANCluster {
public:
    DBSCANCluster(SearchResult &searchResult, std::set<cluster_t> &finalClusters, double minCov) : searchResult(searchResult), finalClusters(finalClusters) {
        cLabel = 0;
        minClusterSize = (unsigned int) ((double) searchResult.qChainKeys.size() * minCov);
        idealClusterSize = std::min(searchResult.qChainKeys.size(), searchResult.dbChainKeys.size());
        prevMaxClusterSize = 0;
        maxDist = 0;
        eps = DEFAULT_EPS;
        learningRate = LEARNING_RATE;
    }

    bool getAlnClusters() {
        // rbh filter
        filterAlnsByRBH();
        fillDistMap();
        // To skip DBSCAN clustering when alignments are few enough.
        if (searchResult.alnVec.size() <= idealClusterSize)
            return checkClusteringNecessity();

        return runDBSCAN();
    }

private:
    SearchResult &searchResult;
    float eps;
    float maxDist;
    float learningRate;
    unsigned int cLabel;
    unsigned int prevMaxClusterSize;
    unsigned int maxClusterSize;
    unsigned int idealClusterSize;
    unsigned int minClusterSize;
    std::vector<unsigned int> neighbors;
    std::vector<unsigned int> neighborsOfCurrNeighbor;
    std::vector<NeighborsWithDist> neighborsWithDist;
    std::set<unsigned int> qFoundChainKeys;
    std::set<unsigned int> dbFoundChainKeys;
    distMap_t distMap;
    std::vector<cluster_t> currClusters;
    std::set<cluster_t> &finalClusters;
    std::map<unsigned int, float> qBestTmScore;
    std::map<unsigned int, float> dbBestTmScore;

    bool runDBSCAN() {
        initializeAlnLabels();
        if (eps >= maxDist)
            return finishDBSCAN();

        for (size_t centerAlnIdx=0; centerAlnIdx < searchResult.alnVec.size(); centerAlnIdx++) {
            ChainToChainAln &centerAln = searchResult.alnVec[centerAlnIdx];
            if (centerAln.label != 0)
                continue;

            getNeighbors(centerAlnIdx, neighbors);
            if (neighbors.size() < MIN_PTS)
                continue;

            centerAln.label = ++cLabel;
            size_t neighborIdx = 0;
            while (neighborIdx < neighbors.size()) {
                unsigned int neighborAlnIdx = neighbors[neighborIdx++];
                if (centerAlnIdx == neighborAlnIdx)
                    continue;

                ChainToChainAln &neighborAln = searchResult.alnVec[neighborAlnIdx];
                neighborAln.label = cLabel;
                getNeighbors(neighborAlnIdx, neighborsOfCurrNeighbor);
                if (neighborsOfCurrNeighbor.size() < MIN_PTS)
                    continue;

                for (auto neighbor : neighborsOfCurrNeighbor) {
                    if (std::find(neighbors.begin(), neighbors.end(), neighbor) == neighbors.end())
                        neighbors.emplace_back(neighbor);
                }
            }
            if (neighbors.size() > idealClusterSize || checkChainRedundancy())
                getNearestNeighbors(centerAlnIdx);

            // too small cluster
            if (neighbors.size() < maxClusterSize)
                continue;

            // new Biggest cluster
            if (neighbors.size() > maxClusterSize) {
                maxClusterSize = neighbors.size();
                currClusters.clear();
            }
            SORT_SERIAL(neighbors.begin(), neighbors.end());
            currClusters.emplace_back(neighbors);
        }

        if (!finalClusters.empty() && currClusters.empty())
            return finishDBSCAN();

        if (maxClusterSize < prevMaxClusterSize)
            return finishDBSCAN();

        if (maxClusterSize > prevMaxClusterSize) {
            finalClusters.clear();
            prevMaxClusterSize = maxClusterSize;
        }
        finalClusters.insert(currClusters.begin(), currClusters.end());
        eps += learningRate;
        return runDBSCAN();
    }

    void fillDistMap() {
        float dist;
        distMap.clear();
        for (size_t i=0; i < searchResult.alnVec.size(); i++) {
            ChainToChainAln &prevAln = searchResult.alnVec[i];
            for (size_t j = i+1; j < searchResult.alnVec.size(); j++) {
                ChainToChainAln &currAln = searchResult.alnVec[j];
                dist = prevAln.getDistance(currAln);
                maxDist = std::max(maxDist, dist);
                distMap.insert({{i,j}, dist});
            }
        }
    }

    void getNeighbors(size_t centerIdx, std::vector<unsigned int> &neighborVec) {
        neighborVec.clear();
        neighborVec.emplace_back(centerIdx);
        for (size_t neighborIdx = 0; neighborIdx < searchResult.alnVec.size(); neighborIdx++) {

            if (neighborIdx == centerIdx)
                continue;

            if (distMap[{std::min(centerIdx, neighborIdx), std::max(centerIdx, neighborIdx)}] >= eps)
                continue;

            neighborVec.emplace_back(neighborIdx);
        }
    }

    void initializeAlnLabels() {
        for (auto &aln : searchResult.alnVec) {
            aln.label = INITIALIZED_LABEL;
        }
        cLabel = INITIALIZED_LABEL;
        maxClusterSize = 0;
        currClusters.clear();
    }

    bool checkChainRedundancy() {
        qFoundChainKeys.clear();
        dbFoundChainKeys.clear();

        for (auto neighborIdx : neighbors) {
            if (!qFoundChainKeys.insert(searchResult.alnVec[neighborIdx].qChain.chainKey).second)
                return true;

            if (!dbFoundChainKeys.insert(searchResult.alnVec[neighborIdx].dbChain.chainKey).second)
                return true;
        }
        return false;
    }

    unsigned int checkClusteringNecessity() {
        if (searchResult.alnVec.empty())
            return false;
        for (size_t alnIdx=0; alnIdx<searchResult.alnVec.size(); alnIdx++) {
            neighbors.emplace_back(alnIdx);
        }
        if (checkChainRedundancy()) {
            neighbors.clear();
            if (searchResult.alnVec.size() < MULTIPLE_CHAINED_COMPLEX)
                return finishDBSCAN();

            return runDBSCAN();
        }
        prevMaxClusterSize = neighbors.size();
        finalClusters.insert(neighbors);
        return finishDBSCAN();
    }

    bool finishDBSCAN() {
        initializeAlnLabels();
        neighbors.clear();
        neighborsOfCurrNeighbor.clear();
        neighborsWithDist.clear();
        qBestTmScore.clear();
        dbBestTmScore.clear();
        qFoundChainKeys.clear();
        dbFoundChainKeys.clear();
        distMap.clear();
        bool isAlnFound = !finalClusters.empty();
        bool isBigEnoughAln = prevMaxClusterSize >= minClusterSize;
        return isAlnFound && isBigEnoughAln;
    }

    void filterAlnsByRBH() {
        float tmScore;
        unsigned int qKey;
        unsigned int dbKey;
        qBestTmScore.clear();
        dbBestTmScore.clear();
        qFoundChainKeys.clear();
        dbFoundChainKeys.clear();
        for (auto qChainKey: searchResult.qChainKeys) {
            qBestTmScore.insert({qChainKey, DEF_TM_SCORE});
        }
        for (auto dbChainKey: searchResult.dbChainKeys) {
            dbBestTmScore.insert({dbChainKey, DEF_TM_SCORE});
        }
        for (auto &aln: searchResult.alnVec) {
            qKey = aln.qChain.chainKey;
            dbKey = aln.dbChain.chainKey;
            tmScore = aln.tmScore;
            qBestTmScore[qKey] = qBestTmScore[qKey] < UNINITIALIZED ? tmScore : std::max(tmScore, qBestTmScore[qKey]);
            dbBestTmScore[dbKey] = dbBestTmScore[dbKey] < UNINITIALIZED ? tmScore : std::max(tmScore, dbBestTmScore[dbKey]);
        }
        size_t alnIdx = 0;
        while (alnIdx < searchResult.alnVec.size()) {
            qKey = searchResult.alnVec[alnIdx].qChain.chainKey;
            dbKey = searchResult.alnVec[alnIdx].dbChain.chainKey;
            tmScore = searchResult.alnVec[alnIdx].tmScore;
            if (tmScore < std::max(qBestTmScore[qKey], dbBestTmScore[dbKey]) * TM_SCORE_MARGIN) {
                searchResult.alnVec.erase(searchResult.alnVec.begin() + alnIdx);
                continue;
            }
            qFoundChainKeys.insert(qKey);
            dbFoundChainKeys.insert(dbKey);
            alnIdx ++;
        }

        if (std::min(qFoundChainKeys.size(), dbFoundChainKeys.size()) < minClusterSize)
            searchResult.alnVec.clear();
    }

    void getNearestNeighbors(unsigned int centerIdx) {
        qFoundChainKeys.clear();
        dbFoundChainKeys.clear();
        neighborsWithDist.clear();
        neighborsWithDist.emplace_back(centerIdx, 0.0);
        for (auto neighborIdx: neighbors) {
            if (neighborIdx == centerIdx)
                continue;
            neighborsWithDist.emplace_back(neighborIdx, distMap[{std::min(centerIdx, neighborIdx), std::max(centerIdx, neighborIdx)}]);
        }
        SORT_SERIAL(neighborsWithDist.begin(), neighborsWithDist.end(), compareNeighborWithDist);
        neighbors.clear();
        for (auto neighborWithDist : neighborsWithDist) {
            if (!qFoundChainKeys.insert(searchResult.alnVec[neighborWithDist.neighbor].qChain.chainKey).second)
                break;
            if (!dbFoundChainKeys.insert(searchResult.alnVec[neighborWithDist.neighbor].dbChain.chainKey).second)
                break;
            neighbors.emplace_back(neighborWithDist.neighbor);
        }
    }
};

class ComplexScorer {
public:
    ComplexScorer(IndexReader *qDbr3Di, IndexReader *tDbr3Di, DBReader<unsigned int> &alnDbr, IndexReader *qCaDbr, IndexReader *tCaDbr, unsigned int thread_idx, double minAssignedChainsRatio) : alnDbr(alnDbr), qCaDbr(qCaDbr), tCaDbr(tCaDbr), thread_idx(thread_idx), minAssignedChainsRatio(minAssignedChainsRatio) {
        maxChainLen = std::max(qDbr3Di->sequenceReader->getMaxSeqLen()+1, tDbr3Di->sequenceReader->getMaxSeqLen()+1);
        q3diDbr = qDbr3Di;
        t3diDbr = tDbr3Di;
        maxResLen = maxChainLen * 2;
        tmAligner = new TMaligner(maxResLen, false, true, false);
    }

    void getSearchResults(unsigned int qComplexId, std::vector<unsigned int> &qChainKeys, chainKeyToComplexId_t &dbChainKeyToComplexIdLookup, complexIdToChainKeys_t &dbComplexIdToChainKeysLookup, std::vector<SearchResult> &searchResults) {
        hasBacktrace = false;
        unsigned int qResLen = getQueryResidueLength(qChainKeys);
        if (qResLen == 0) return;
        paredSearchResult = SearchResult(qChainKeys, qResLen);
        // for each chain from the query Complex
        for (auto qChainKey: qChainKeys) {
            unsigned int qKey = alnDbr.getId(qChainKey);
            if (qKey == NOT_AVAILABLE_CHAIN_KEY) continue;
            char *data = alnDbr.getData(qKey, thread_idx);
            if (*data == '\0') continue;
            qAlnResult = Matcher::parseAlignmentRecord(data);
            size_t qDbId = qCaDbr->sequenceReader->getId(qChainKey);
            char *qCaData = qCaDbr->sequenceReader->getData(qDbId, thread_idx);
            size_t qCaLength = qCaDbr->sequenceReader->getEntryLen(qDbId);
            unsigned int &qLen = qAlnResult.qLen;
            float *queryCaData = qCoords.read(qCaData, qAlnResult.qLen, qCaLength);
            qChain = Chain(qComplexId, qChainKey);
            tmAligner->initQuery(queryCaData, &queryCaData[qLen], &queryCaData[qLen * 2], NULL, qLen);
            // for each alignment from the query chain
            while (*data != '\0') {
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const auto dbChainKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                const unsigned int dbComplexId = dbChainKeyToComplexIdLookup.at(dbChainKey);
                dbAlnResult = Matcher::parseAlignmentRecord(data);
                data = Util::skipLine(data);
                if (dbAlnResult.backtrace.empty()) continue;
                hasBacktrace = true;
                size_t tCaId = tCaDbr->sequenceReader->getId(dbChainKey);
                char *tCaData = tCaDbr->sequenceReader->getData(tCaId, thread_idx);
                size_t tCaLength = tCaDbr->sequenceReader->getEntryLen(tCaId);
                unsigned int & dbLen = dbAlnResult.dbLen;
                float *targetCaData = tCoords.read(tCaData, dbLen, tCaLength);
                dbChain = Chain(dbComplexId, dbChainKey);
                tmResult = tmAligner->computeTMscore(targetCaData,&targetCaData[dbLen],&targetCaData[dbLen * 2],dbLen,dbAlnResult.qStartPos,dbAlnResult.dbStartPos,Matcher::uncompressAlignment(dbAlnResult.backtrace),dbAlnResult.qLen);
                currAln =  ChainToChainAln(qChain, dbChain, queryCaData, targetCaData, dbAlnResult, tmResult);
                currAlns.emplace_back(currAln);
                currAln.free();
            } // while end
        } // for end
        // When alignments have no backtrace
        if (!hasBacktrace) {
            Debug(Debug::ERROR) << "Backtraces are required. Please run search with '-a' option.\n";
            EXIT(EXIT_FAILURE);
        }
        if (currAlns.empty())
            return;
        SORT_SERIAL(currAlns.begin(), currAlns.end(), compareChainToChainAlnByDbComplexId);
        unsigned int currDbComplexId = currAlns[0].dbChain.complexId;
        std::vector<unsigned int> currDbChainKeys = dbComplexIdToChainKeysLookup.at(currDbComplexId);
        unsigned int currDbResLen = getDbResidueLength(currDbChainKeys);
        paredSearchResult.resetDbComplex(currDbChainKeys, currDbResLen);
        for (auto &aln: currAlns) {
            if (aln.dbChain.complexId == currDbComplexId) {
                paredSearchResult.alnVec.emplace_back(aln);
                continue;
            }
            paredSearchResult.standardize();
            if (!paredSearchResult.alnVec.empty())
                searchResults.emplace_back(paredSearchResult);

            paredSearchResult.alnVec.clear();
            currDbComplexId = aln.dbChain.complexId;
            currDbChainKeys = dbComplexIdToChainKeysLookup.at(currDbComplexId);
            currDbResLen = getDbResidueLength(currDbChainKeys);
            paredSearchResult.resetDbComplex(currDbChainKeys, currDbResLen);
            paredSearchResult.alnVec.emplace_back(aln);
        }
        currAlns.clear();
        paredSearchResult.standardize();
        if (!paredSearchResult.alnVec.empty() && currDbChainKeys.size() >= MULTIPLE_CHAINED_COMPLEX)
            searchResults.emplace_back(paredSearchResult);

        paredSearchResult.alnVec.clear();
    }

    void getAssignments(SearchResult &searchResult, std::vector<Assignment> &assignments) {
        if (maxResLen < maxChainLen * std::min(searchResult.qChainKeys.size(),  searchResult.dbChainKeys.size())) {
            delete tmAligner;
            maxResLen = std::max(searchResult.qChainKeys.size(), searchResult.dbChainKeys.size()) * maxChainLen;
            tmAligner = new TMaligner(maxResLen, false, true, false);
        }
        finalClusters.clear();
        DBSCANCluster dbscanCluster = DBSCANCluster(searchResult, finalClusters, minAssignedChainsRatio);
        if (!dbscanCluster.getAlnClusters()) {
            finalClusters.clear();
            return;
        }
        assignment = Assignment(searchResult.qResidueLen, searchResult.dbResidueLen);
        for (auto &cluster: finalClusters) {
            for (auto alnIdx: cluster) {
                assignment.appendChainToChainAln(searchResult.alnVec[alnIdx]);
            }
            assignment.getTmScore(*tmAligner);
            assignment.updateResultToWriteLines();
            assignments.emplace_back(assignment);
            assignment.reset();
        }
        finalClusters.clear();
    }

    void free() {
        delete tmAligner;
    }

private:
    TMaligner *tmAligner;
    TMaligner::TMscoreResult tmResult;
    Matcher::result_t qAlnResult;
    Matcher::result_t dbAlnResult;
    unsigned int maxChainLen;
    DBReader<unsigned int> &alnDbr;
    IndexReader *qCaDbr;
    IndexReader *tCaDbr;
    IndexReader *q3diDbr;
    IndexReader *t3diDbr;
    Coordinate16 qCoords;
    Coordinate16 tCoords;
    unsigned int thread_idx;
    double minAssignedChainsRatio;
    unsigned int maxResLen;
    Chain qChain;
    Chain dbChain;
    ChainToChainAln currAln;
    std::vector<ChainToChainAln> currAlns;
    Assignment assignment;
    SearchResult paredSearchResult;
    std::set<cluster_t> finalClusters;
    bool hasBacktrace;

    unsigned int getQueryResidueLength(std::vector<unsigned int> &qChainKeys) {
        unsigned int qResidueLen = 0;
        size_t qDbId;
        for (auto qChainKey: qChainKeys) {
            qDbId = q3diDbr->sequenceReader->getId(qChainKey);
            // Not accessible
            if (qDbId == NOT_AVAILABLE_CHAIN_KEY)
                return 0;

            qResidueLen += q3diDbr->sequenceReader->getSeqLen(qDbId);
        }
        return qResidueLen;
    }

    unsigned int getDbResidueLength(std::vector<unsigned int> &dbChainKeys) {
        unsigned int dbResidueLen = 0;
        size_t tDbId;
        for (auto dbChainKey: dbChainKeys) {
            tDbId = t3diDbr->sequenceReader->getId(dbChainKey);
            // Not accessible
            if (tDbId == NOT_AVAILABLE_CHAIN_KEY)
                return 0;

            dbResidueLen += t3diDbr->sequenceReader->getSeqLen(tDbId);
        }
        return dbResidueLen;
    }
};

int scorecomplex(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(alnDbr.getDbtype());
    int dbType = Parameters::DBTYPE_ALIGNMENT_RES;
    bool needSrc = false;
    if (extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC) {
        needSrc = true;
        dbType = DBReader<unsigned int>::setExtendedDbtype(dbType, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbType);
    resultWriter.open();

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    std::string t3DiDbrName =  StructureUtil::getIndexWithSuffix(par.db2, "_ss");
    bool is3DiIdx = Parameters::isEqualDbtype(FileUtil::parseDbType(t3DiDbrName.c_str()), Parameters::DBTYPE_INDEX_DB);
    IndexReader t3DiDbr(
            is3DiIdx ? t3DiDbrName : par.db2,
            par.threads,
            needSrc ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
            touch ? IndexReader::PRELOAD_INDEX : 0,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
            needSrc ? "_seq_ss" : "_ss"
    );
    IndexReader tCaDbr(
            par.db2,
            par.threads,
            needSrc
            ? IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB2)
            : IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
            touch ? IndexReader::PRELOAD_INDEX : 0,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
            needSrc ? "_seq_ca" : "_ca"
    );
    IndexReader* q3DiDbr = NULL;
    IndexReader* qCaDbr = NULL;
    bool sameDB = false;
    if (par.db1 == par.db2) {
        sameDB = true;
        q3DiDbr = &t3DiDbr;
        qCaDbr = &tCaDbr;
    } else {
        q3DiDbr = new IndexReader(
                StructureUtil::getIndexWithSuffix(par.db1, "_ss"),
                par.threads, IndexReader::SEQUENCES,
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA
        );
        qCaDbr = new IndexReader(
                par.db1,
                par.threads,
                IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                "_ca"
        );
    }

    double minAssignedChainsRatio = par.minAssignedChainsThreshold > MAX_ASSIGNED_CHAIN_RATIO ? MAX_ASSIGNED_CHAIN_RATIO: par.minAssignedChainsThreshold;

    std::vector<unsigned int> qComplexIndices;
    std::vector<unsigned int> dbComplexIndices;
    chainKeyToComplexId_t qChainKeyToComplexIdMap;
    chainKeyToComplexId_t dbChainKeyToComplexIdMap;
    complexIdToChainKeys_t dbComplexIdToChainKeysMap;
    complexIdToChainKeys_t qComplexIdToChainKeysMap;
    std::string qLookupFile = par.db1 + ".lookup";
    std::string dbLookupFile = par.db2 + ".lookup";
    getKeyToIdMapIdToKeysMapIdVec(qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIndices);
    getKeyToIdMapIdToKeysMapIdVec(dbLookupFile, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, dbComplexIndices);
    qChainKeyToComplexIdMap.clear();
    dbComplexIndices.clear();
    Debug::Progress progress(qComplexIndices.size());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
        char buffer[4096];
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<SearchResult> searchResults;
        std::vector<Assignment> assignments;
        std::vector<resultToWrite_t> resultToWriteLines;
        ComplexScorer complexScorer(q3DiDbr, &t3DiDbr, alnDbr, qCaDbr, &tCaDbr, thread_idx, minAssignedChainsRatio);
#pragma omp for schedule(dynamic, 1)
        // for each q complex
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
            unsigned int qComplexId = qComplexIndices[qCompIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
            if (qChainKeys.size() < MULTIPLE_CHAINED_COMPLEX)
                continue;
            complexScorer.getSearchResults(qComplexId, qChainKeys, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, searchResults);
            // for each db complex
            for (size_t dbId = 0; dbId < searchResults.size(); dbId++) {
                complexScorer.getAssignments(searchResults[dbId], assignments);
            }
            SORT_SERIAL(assignments.begin(), assignments.end(), compareAssignment);
            // for each query chain key
            for (size_t qChainKeyIdx = 0; qChainKeyIdx < qChainKeys.size(); qChainKeyIdx++) {
                resultToWriteLines.emplace_back("");
            }
            // for each assignment
            for (unsigned int assignmentId = 0; assignmentId < assignments.size(); assignmentId++){
                Assignment &assignment = assignments[assignmentId];
                // for each output line from this assignment
                for (size_t resultToWriteIdx = 0; resultToWriteIdx < assignment.resultToWriteLines.size(); resultToWriteIdx++) {
                    unsigned int &qKey = assignment.resultToWriteLines[resultToWriteIdx].first;
                    resultToWrite_t &resultToWrite = assignment.resultToWriteLines[resultToWriteIdx].second;
                    snprintf(buffer, sizeof(buffer), "%s\t%d\n", resultToWrite.c_str(), assignmentId);
                    unsigned int currIdx = find(qChainKeys.begin(), qChainKeys.end(), qKey) - qChainKeys.begin();
                    resultToWriteLines[currIdx].append(buffer);
                }
            }
            for (size_t qChainKeyIdx = 0; qChainKeyIdx < qChainKeys.size(); qChainKeyIdx++) {
                resultToWrite_t &resultToWrite = resultToWriteLines[qChainKeyIdx];
                unsigned int & qKey = qChainKeys[qChainKeyIdx];
                resultWriter.writeData(resultToWrite.c_str(),resultToWrite.length(),qKey,thread_idx);
            }
            assignments.clear();
            resultToWriteLines.clear();
            searchResults.clear();
            progress.updateProgress();
        }
        complexScorer.free();
    }

    qComplexIndices.clear();
    dbChainKeyToComplexIdMap.clear();
    dbComplexIdToChainKeysMap.clear();
    qComplexIdToChainKeysMap.clear();
    alnDbr.close();
    if (!sameDB) {
        delete q3DiDbr;
        delete qCaDbr;
    }
    resultWriter.close(true);
    return EXIT_SUCCESS;
}
