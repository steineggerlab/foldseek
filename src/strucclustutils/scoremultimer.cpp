#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "StructureUtil.h"
#include "Coordinate16.h"
#include "MultimerUtil.h"
#include "set"
#include "unordered_set"
#include "tmalign/basic_fun.h"
#include "FileUtil.h"
#include "LDDT.h"
#include <map>

#ifdef OPENMP
#include <omp.h>
#endif
#define INTERFACE_THRESHOLD 8

// carrying chainToChainAlignments from the same query and target complex
struct SearchResult {
    SearchResult() {}
    SearchResult(std::vector<unsigned int> &qChainKeys, unsigned int qResidueLen, std::vector<unsigned int> &dbChainKeys, unsigned int dbResidueLen, const std::vector<ChainToChainAln> &alnVec) : qChainKeys(qChainKeys), dbChainKeys(dbChainKeys), qResidueLen(qResidueLen), dbResidueLen(dbResidueLen), alnVec(alnVec) {}

    std::vector<unsigned int> qChainKeys;
    std::vector<unsigned int> dbChainKeys;
    unsigned int qResidueLen;
    unsigned int dbResidueLen;
    std::vector<ChainToChainAln> alnVec;

    void resetDbComplex(std::vector<unsigned int> &chainKeys, unsigned int residueLen) {
        dbChainKeys = chainKeys;
        dbResidueLen = residueLen;
    }

    void standardize(int MonomerIncludeMode) {
        if (dbResidueLen == 0)
            alnVec.clear();

        if (MonomerIncludeMode == SKIP_MONOMERS && dbChainKeys.size() < MULTIPLE_CHAINED_COMPLEX)
            alnVec.clear();

        if (alnVec.empty())
            return;

        double length = alnVec.size();
        for (size_t i = 0; i < SIZE_OF_SUPERPOSITION_VECTOR; i++) {
            double mean = 0.0;
            double var = 0.0;
            for (auto &aln : alnVec) {
                mean += aln.superposition[i] / length;
            }
            for (auto &aln : alnVec) {
                var += std::pow(aln.superposition[i] - mean, 2) / length;
            }
            double sd = std::sqrt(var);
            double cv = (abs(mean) > TOO_SMALL_MEAN) ? sd / std::abs(mean) : sd;
            for (auto &aln : alnVec) {
                aln.superposition[i] = cv < TOO_SMALL_CV ? FILTERED_OUT : (aln.superposition[i] - mean) / sd;
            }
        }
    }

    void clear() {
        alnVec.clear();
        qChainKeys.clear();
        dbChainKeys.clear();
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
    std::vector<unsigned int> matchLenVec; 
    double qTmScore;
    double dbTmScore;
    std::string tString;
    std::string uString;
    std::string backtrace;
    resultToWrite_t assignmentResult;
    TMaligner::TMscoreResult tmResult;
    std::vector<resultToWriteWithKey_t> chainToChainResults;
    resultToWrite_t filterResult;
    unsigned int assignmentId;

    void appendChainToChainAln(ChainToChainAln &aln) {
        matches += aln.matches;
        qCaXVec.insert(qCaXVec.end(), aln.qChain.caVecX.begin(), aln.qChain.caVecX.end());
        qCaYVec.insert(qCaYVec.end(), aln.qChain.caVecY.begin(), aln.qChain.caVecY.end());
        qCaZVec.insert(qCaZVec.end(), aln.qChain.caVecZ.begin(), aln.qChain.caVecZ.end());
        dbCaXVec.insert(dbCaXVec.end(), aln.dbChain.caVecX.begin(), aln.dbChain.caVecX.end());
        dbCaYVec.insert(dbCaYVec.end(), aln.dbChain.caVecY.begin(), aln.dbChain.caVecY.end());
        dbCaZVec.insert(dbCaZVec.end(), aln.dbChain.caVecZ.begin(), aln.dbChain.caVecZ.end());
        matchLenVec.push_back(aln.matches);
        chainToChainResults.emplace_back(aln.qChain.chainKey, aln.resultToWrite);
    }

    void reset() {
        matches = 0;
        qCaXVec.clear();
        qCaYVec.clear();
        qCaZVec.clear();
        dbCaXVec.clear();
        dbCaYVec.clear();
        dbCaZVec.clear();
        uString.clear();
        tString.clear();
        chainToChainResults.clear();
        assignmentResult.clear();
        filterResult.clear();
    }

    bool getTmScore(TMaligner &tmAligner) {
        // for safety
        if (matches == 0) {
            return false;
        }
        backtrace = std::string(matches, 'M');
        unsigned int normLen = std::min(qResidueLength, dbResidueLength);
        tmAligner.initQuery(&qCaXVec[0], &qCaYVec[0], &qCaZVec[0], NULL, matches);
        tmResult = tmAligner.computeTMscore(&dbCaXVec[0], &dbCaYVec[0], &dbCaZVec[0],matches,0,0, backtrace,normLen);
        qTmScore = tmResult.tmscore * normLen / qResidueLength;
        dbTmScore = tmResult.tmscore * normLen / dbResidueLength;
        return true;
    }

    void updateResultToWriteLines() {
        char sep = ',';
        char tab = '\t';
        tString.append(std::to_string(tmResult.t[0]) + sep + std::to_string(tmResult.t[1]) + sep + std::to_string(tmResult.t[2]));
        uString.append(std::to_string(tmResult.u[0][0]) + sep + std::to_string(tmResult.u[0][1]) + sep + std::to_string(tmResult.u[0][2]) + sep);
        uString.append(std::to_string(tmResult.u[1][0]) + sep + std::to_string(tmResult.u[1][1]) + sep + std::to_string(tmResult.u[1][2]) + sep);
        uString.append(std::to_string(tmResult.u[2][0]) + sep + std::to_string(tmResult.u[2][1]) + sep + std::to_string(tmResult.u[2][2]));
        assignmentResult = tab + std::to_string(qTmScore) + tab + std::to_string(dbTmScore) + tab + uString + tab + tString;
        uString.clear();
        tString.clear();
    }

    void getChainToChainResult(const unsigned int qChainKey, std::string &chainToChainResult) const {
        for (auto & chainResultWithKey: chainToChainResults) {
            if (chainResultWithKey.first != qChainKey) {
                continue;
            }
            chainToChainResult = chainResultWithKey.second;
            return;
        }
    }
};

struct NeighborsWithDist {
    NeighborsWithDist(unsigned int neighbor, float dist) : neighbor(neighbor), dist(dist) {}
    unsigned int neighbor;
    float dist;
};

bool compareChainToChainAlnByDbComplexId(const ChainToChainAln &first, const ChainToChainAln &second) {
    if (first.dbChain.complexId < second.dbChain.complexId) {
        return true;
    }
    if (first.dbChain.complexId > second.dbChain.complexId) {
        return false;
    }

    if (first.qChain.chainKey < second.qChain.chainKey) {
        return true;
    }
    if (first.qChain.chainKey > second.qChain.chainKey) {
        return false;
    }

    if (first.dbChain.chainKey < second.dbChain.chainKey) {
        return true;
    }
    if (first.dbChain.chainKey > second.dbChain.chainKey) {
        return false;
    }
    return false;
}

bool compareChainToChainAln(const ChainToChainAln &first, const ChainToChainAln &second) {
    if (first.qChain.chainKey < second.qChain.chainKey) {
        return true;
    }
    if (first.qChain.chainKey > second.qChain.chainKey) {
        return false;
    }

    if (first.dbChain.chainKey < second.dbChain.chainKey) {
        return true;
    }
    if (first.dbChain.chainKey > second.dbChain.chainKey) {
        return false;
    }
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

void getResult(std::string &result, std::string &currentResult, const Assignment &assignment) {
    currentResult.append("\t" + assignment.assignmentResult + "\t" + std::to_string(assignment.assignmentId));
    if (!assignment.filterResult.empty()) {
        currentResult.append("\t" + assignment.filterResult);
    }
    currentResult.append("\n");
    result.append(currentResult);
}

class DBSCANCluster {
public:
    DBSCANCluster(SearchResult &searchResult, std::set<cluster_t> &finalClusters, float minCov) : searchResult(searchResult), finalClusters(finalClusters) {
        cLabel = 0;
        minimumClusterSize = std::ceil((float) searchResult.qChainKeys.size() * minCov);
        maximumClusterSize = std::min(searchResult.qChainKeys.size(), searchResult.dbChainKeys.size());
        maximumClusterNum = searchResult.alnVec.size() / maximumClusterSize;
        prevMaxClusterSize = 0;
        maxDist = FLT_MIN;
        minDist = FLT_MAX;
        learningRate = LEARNING_RATE;
    }

    bool getAlnClusters() {
        // if Query or Target is a Monomer Complex.
        if (std::min(searchResult.qChainKeys.size(), searchResult.dbChainKeys.size()) < MULTIPLE_CHAINED_COMPLEX)
            return earlyStopForMonomers();

        // rbh filter
        filterAlnsByRBH();
        fillDistMatrix();
        // To skip DBSCAN clustering when alignments are few enough.
        if (searchResult.alnVec.size() <= maximumClusterSize)
            return checkClusteringNecessity();

        return runDBSCAN();
    }

private:
    SearchResult &searchResult;
    float eps;
    float maxDist;
    float minDist;
    float learningRate;
    unsigned int cLabel;
    unsigned int maximumClusterNum;
    unsigned int prevMaxClusterSize;
    unsigned int currMaxClusterSize;
    unsigned int maximumClusterSize;
    unsigned int minimumClusterSize;
    std::vector<unsigned int> neighbors;
    std::vector<unsigned int> neighborsOfCurrNeighbor;
    std::unordered_set<unsigned int> foundNeighbors;
    std::vector<NeighborsWithDist> neighborsWithDist;
    std::unordered_set<unsigned int> qFoundChainKeys;
    std::unordered_set<unsigned int> dbFoundChainKeys;
    std::vector<float> distMatrix;
    std::vector<cluster_t> currClusters;
    std::set<cluster_t> &finalClusters;
    std::map<unsigned int, float> qBestTmScore;
    std::map<unsigned int, float> dbBestTmScore;

    bool earlyStopForMonomers() {
        if (minimumClusterSize >= MULTIPLE_CHAINED_COMPLEX)
            return finishDBSCAN();

        getSingleChainedCluster();
        return finishDBSCAN();
    }

    void getSingleChainedCluster() {
        finalClusters.clear();
        for (unsigned int alnIdx = 0; alnIdx < searchResult.alnVec.size(); alnIdx++ ) {
            neighbors = {alnIdx};
            finalClusters.insert(neighbors);
        }
    }

    bool runDBSCAN() {
        unsigned int neighborIdx;
        unsigned int neighborAlnIdx;
        while (eps < maxDist) {
            initializeAlnLabels();
            for (size_t centerAlnIdx = 0; centerAlnIdx < searchResult.alnVec.size(); centerAlnIdx++) {
                ChainToChainAln &centerAln = searchResult.alnVec[centerAlnIdx];
                if (centerAln.label != 0)
                    continue;

                getNeighbors(centerAlnIdx, neighbors);
                if (neighbors.size() < MIN_PTS)
                    continue;

                centerAln.label = ++cLabel;
                foundNeighbors.clear();
                foundNeighbors.insert(neighbors.begin(), neighbors.end());
                neighborIdx = 0;
                while (neighborIdx < neighbors.size()) {
                    neighborAlnIdx = neighbors[neighborIdx++];
                    if (centerAlnIdx == neighborAlnIdx)
                        continue;

                    ChainToChainAln &neighborAln = searchResult.alnVec[neighborAlnIdx];
                    neighborAln.label = cLabel;
                    getNeighbors(neighborAlnIdx, neighborsOfCurrNeighbor);
                    if (neighborsOfCurrNeighbor.size() < MIN_PTS)
                        continue;

                    for (auto neighbor : neighborsOfCurrNeighbor) {
                        if (foundNeighbors.insert(neighbor).second)
                            neighbors.emplace_back(neighbor);
                    }
                }
                if (neighbors.size() > maximumClusterSize || checkChainRedundancy())
                    getNearestNeighbors(centerAlnIdx);

                // too small cluster
                if (neighbors.size() < currMaxClusterSize)
                    continue;

                // new Biggest cluster
                if (neighbors.size() > currMaxClusterSize) {
                    currMaxClusterSize = neighbors.size();
                    currClusters.clear();
                }
                SORT_SERIAL(neighbors.begin(), neighbors.end());
                currClusters.emplace_back(neighbors);
            }

            if (!finalClusters.empty() && currClusters.empty())
                return finishDBSCAN();

            if (currMaxClusterSize < prevMaxClusterSize)
                return finishDBSCAN();

            if (currMaxClusterSize > prevMaxClusterSize) {
                finalClusters.clear();
                prevMaxClusterSize = currMaxClusterSize;
            }

            if (currMaxClusterSize >= minimumClusterSize)
                finalClusters.insert(currClusters.begin(), currClusters.end());

            if (currMaxClusterSize == maximumClusterSize && finalClusters.size() == maximumClusterNum)
                return finishDBSCAN();

            eps += learningRate;
        }

        if (minimumClusterSize < MULTIPLE_CHAINED_COMPLEX && prevMaxClusterSize < MULTIPLE_CHAINED_COMPLEX)
            getSingleChainedCluster();

        return finishDBSCAN();
    }

    size_t getDistMatrixIndex(size_t i, size_t j) const {
            if (i > j) std::swap(i, j); // Ensure i <= j for symmetry
            size_t n = searchResult.alnVec.size();
            return (2 * n *i - i - i * i) / 2 + j - i - 1;
    }

    void fillDistMatrix() {
        size_t size = searchResult.alnVec.size();
        float dist;
        distMatrix.resize(size * (size - 1) / 2, 0.0f);
        for (size_t i = 0; i < searchResult.alnVec.size(); i++) {
            const ChainToChainAln &prevAln = searchResult.alnVec[i];
            for (size_t j = i + 1; j < searchResult.alnVec.size(); j++) {
                const ChainToChainAln &currAln = searchResult.alnVec[j];
                dist = prevAln.getDistance(currAln);
                maxDist = std::max(maxDist, dist);
                minDist = std::min(minDist, dist);
                distMatrix[getDistMatrixIndex(i, j)] = dist;
            }
        }
        eps = minDist;
    }

    void getNeighbors(size_t centerIdx, std::vector<unsigned int> &neighborVec) const {
        neighborVec.clear();
        neighborVec.emplace_back(centerIdx);
        for (size_t neighborIdx = 0; neighborIdx < searchResult.alnVec.size(); neighborIdx++) {

            if (neighborIdx == centerIdx)
                continue;

            if (distMatrix[getDistMatrixIndex(centerIdx, neighborIdx)] >= eps)
                continue;

            neighborVec.emplace_back(neighborIdx);
        }
    }

    void initializeAlnLabels() {
        for (auto &aln : searchResult.alnVec) {
            aln.label = INITIALIZED_LABEL;
        }
        cLabel = INITIALIZED_LABEL;
        currMaxClusterSize = 0;
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

    bool checkClusteringNecessity() {
        // Too few alns => do nothing and finish it
        if (searchResult.alnVec.size() < minimumClusterSize)
            return finishDBSCAN();
        // All alns as a cluster
        for (size_t alnIdx=0; alnIdx<searchResult.alnVec.size(); alnIdx++) {
            neighbors.emplace_back(alnIdx);
        }
        // Redundant chains => DBSCAN clustering
        if (checkChainRedundancy()) {
            neighbors.clear();
            return runDBSCAN();
        }
        // Already good => finish it without clustering
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
        distMatrix.clear();
        return !finalClusters.empty();
    }

    void filterAlnsByRBH() {
        qBestTmScore.clear();
        dbBestTmScore.clear();
        qFoundChainKeys.clear();
        dbFoundChainKeys.clear();
        for (auto qChainKey: searchResult.qChainKeys) {
            qBestTmScore.insert({qChainKey, FLT_MIN});
        }
        for (auto dbChainKey: searchResult.dbChainKeys) {
            dbBestTmScore.insert({dbChainKey, FLT_MIN});
        }
        for (auto &aln: searchResult.alnVec) {
            unsigned int qKey = aln.qChain.chainKey;
            unsigned int dbKey = aln.dbChain.chainKey;
            float tmScore = aln.tmScore;
            qBestTmScore[qKey] = std::max(tmScore, qBestTmScore[qKey]);
            dbBestTmScore[dbKey] = std::max(tmScore, dbBestTmScore[dbKey]);
        }
        size_t alnIdx = 0;
        while (alnIdx < searchResult.alnVec.size()) {
            unsigned int qKey = searchResult.alnVec[alnIdx].qChain.chainKey;
            unsigned int dbKey = searchResult.alnVec[alnIdx].dbChain.chainKey;
            float tmScore = searchResult.alnVec[alnIdx].tmScore;
            if (tmScore < std::max(qBestTmScore[qKey], dbBestTmScore[dbKey]) * TM_SCORE_MARGIN) {
                searchResult.alnVec.erase(searchResult.alnVec.begin() + alnIdx);
                continue;
            }
            qFoundChainKeys.insert(qKey);
            dbFoundChainKeys.insert(dbKey);
            alnIdx++;
        }

        if (std::min(qFoundChainKeys.size(), dbFoundChainKeys.size()) < minimumClusterSize)
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
            neighborsWithDist.emplace_back(neighborIdx, distMatrix[getDistMatrixIndex(centerIdx, neighborIdx)]);
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
    ComplexScorer(IndexReader *q3diDbr, IndexReader *t3diDbr, DBReader<unsigned int> &alnDbr, DBReader<unsigned int> *qCaDbr, DBReader<unsigned int> *tCaDbr, unsigned int thread_idx, float minAssignedChainsRatio, int monomerIncludeMode)
        : q3diDbr(q3diDbr), t3diDbr(t3diDbr), alnDbr(alnDbr), qCaDbr(qCaDbr), tCaDbr(tCaDbr), thread_idx(thread_idx), minAssignedChainsRatio(minAssignedChainsRatio), monomerIncludeMode(monomerIncludeMode)  {
        maxChainLen = std::max(q3diDbr->sequenceReader->getMaxSeqLen()+1, t3diDbr->sequenceReader->getMaxSeqLen()+1);
        maxResLen = maxChainLen * 2;
        tmAligner = new TMaligner(maxResLen, false, true, false);
    }
    ~ComplexScorer() {
        delete tmAligner;
    }

    void getSearchResultLinesMap(std::vector<unsigned int> &qChainKeys, alignmentLinesMap_t &alignmentLinesMap) {
        qResLen = getQueryResidueLength(qChainKeys);
        if (qResLen == 0) {
            return;
        }

        for (auto qChainKey: qChainKeys) {
            unsigned int qDbKey = alnDbr.getId(qChainKey);
            if (qDbKey == NOT_AVAILABLE_CHAIN_KEY) {
                continue;
            }
            char *data = alnDbr.getData(qDbKey, thread_idx);
            size_t dataSize = alnDbr.getDataSize();
            if (*data == '\0') {
                continue;
            }
            char dbKeyBuffer[256];
            char lineBuffer[1024];
            while (*data != '\0') {
                Util::parseKey(data, dbKeyBuffer);
                const auto dbChainKey = static_cast<unsigned int>(strtoul(dbKeyBuffer, NULL, 10));
                Util::getLine(data, dataSize, lineBuffer, 1024);
                alignmentLinesMap.insert({{qChainKey, dbChainKey}, static_cast<std::string>(lineBuffer)});
                data = Util::skipLine(data);
            } // while end
        } // for end
    }

    void getSearchResultByDbComplex(unsigned int qComplexId, unsigned int dbComplexId, std::vector<unsigned int> &qChainKeys, std::vector<unsigned int> &dbChainKeys, alignmentLinesMap_t &alignmentLinesMap, SearchResult &searchResult) {
        bool hasBacktrace = false;
        unsigned int dbResLen = getDbResidueLength(dbChainKeys);
        if (dbResLen == 0) {
            removeCurrAlnLines(qChainKeys, dbChainKeys, alignmentLinesMap);
            return;
        }

        Coordinate16 qCoords;
        Coordinate16 tCoords;
        currAlns.clear();
        Chain qChain;
        Chain dbChain;
        // for each chain from the query Complex
        for (auto qChainKey: qChainKeys) {
            Matcher::result_t qAlnResult;
            if (!getQueryAlnResult(qChainKey, dbChainKeys, alignmentLinesMap, qAlnResult)) {
                return;
            }
            size_t qDbId = qCaDbr->getId(qChainKey);
            char *qCaData = qCaDbr->getData(qDbId, thread_idx);
            size_t qCaLength = qCaDbr->getEntryLen(qDbId);
            unsigned int qLen = qAlnResult.qLen;
            float *queryCaData = qCoords.read(qCaData, qAlnResult.qLen, qCaLength);
            qChain = Chain(qComplexId, qChainKey);
            tmAligner->initQuery(queryCaData, &queryCaData[qLen], &queryCaData[qLen * 2], NULL, qLen);
            // for each alignment from the query chain
            for (auto dbChainKey: dbChainKeys) {
                std::string &data = alignmentLinesMap[{qChainKey, dbChainKey}];
                if (data.empty()) {
                    alignmentLinesMap.erase({qChainKey, dbChainKey});
                    continue;
                }
                Matcher::result_t dbAlnResult = Matcher::parseAlignmentRecord(data.c_str());
                if (dbAlnResult.backtrace.empty()) {
                    alignmentLinesMap.erase({qChainKey, dbChainKey});
                    continue;
                }
                hasBacktrace = true;
                size_t tCaId = tCaDbr->getId(dbChainKey);
                char *tCaData = tCaDbr->getData(tCaId, thread_idx);
                size_t tCaLength = tCaDbr->getEntryLen(tCaId);
                unsigned int dbLen = dbAlnResult.dbLen;
                float *targetCaData = tCoords.read(tCaData, dbLen, tCaLength);
                dbChain = Chain(dbComplexId, dbChainKey);
                TMaligner::TMscoreResult tmResult = tmAligner->computeTMscore(targetCaData, &targetCaData[dbLen], &targetCaData[dbLen * 2], dbLen, dbAlnResult.qStartPos, dbAlnResult.dbStartPos, Matcher::uncompressAlignment(dbAlnResult.backtrace), dbAlnResult.qLen);
                currAlns.emplace_back(qChain, dbChain, queryCaData, targetCaData, dbAlnResult, tmResult);
                alignmentLinesMap.erase({qChainKey, dbChainKey});
            }
        } // for end

        // When no alignment is found.
        if (currAlns.empty()) {
            return;
        }

        // When alignments have no backtrace
        if (!hasBacktrace) {
            Debug(Debug::ERROR) << "Backtraces are required. Please run search with '-a' option.\n";
            EXIT(EXIT_FAILURE);
        }

        SORT_SERIAL(currAlns.begin(), currAlns.end(), compareChainToChainAln);
        searchResult = SearchResult(qChainKeys, qResLen, dbChainKeys, dbResLen, currAlns);
        searchResult.standardize(monomerIncludeMode);
    }

    void getAssignments(SearchResult &searchResult, std::vector<Assignment> &assignments) {
        if (maxResLen < maxChainLen * std::min(searchResult.qChainKeys.size(),  searchResult.dbChainKeys.size())) {
            delete tmAligner;
            maxResLen = std::max(searchResult.qChainKeys.size(), searchResult.dbChainKeys.size()) * maxChainLen;
            tmAligner = new TMaligner(maxResLen, false, true, false);
        }
        finalClusters.clear();
        DBSCANCluster dbscanCluster(searchResult, finalClusters, minAssignedChainsRatio);
        if (!dbscanCluster.getAlnClusters()) {
            finalClusters.clear();
            return;
        }
        for (auto &cluster : finalClusters) {
            Assignment assignment(searchResult.qResidueLen, searchResult.dbResidueLen);
            for (auto alnIdx : cluster) {
                assignment.appendChainToChainAln(searchResult.alnVec[alnIdx]);
            }
            // matches==0 ? skip
            if (!assignment.getTmScore(*tmAligner)) {
                continue;
            }
            assignment.updateResultToWriteLines();
            assignments.emplace_back(std::move(assignment));
        }
        finalClusters.clear();
    }

private:
    IndexReader *q3diDbr;
    IndexReader *t3diDbr;
    DBReader<unsigned int> &alnDbr;
    DBReader<unsigned int> *qCaDbr;
    DBReader<unsigned int> *tCaDbr;
    const unsigned int thread_idx;
    const float minAssignedChainsRatio;
    const int monomerIncludeMode;
    TMaligner *tmAligner;
    unsigned int maxChainLen;
    unsigned int maxResLen;
    std::vector<ChainToChainAln> currAlns;
    std::set<cluster_t> finalClusters;
    unsigned int qResLen;

    unsigned int getQueryResidueLength(std::vector<unsigned int> &qChainKeys) const {
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

    unsigned int getDbResidueLength(const std::vector<unsigned int> &dbChainKeys) const {
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

    static void removeCurrAlnLines(const std::vector<unsigned int> &qChainKeys, const std::vector<unsigned int> &dbChainKeys,  alignmentLinesMap_t &alignmentLinesMap) {
        for (auto qChainKey: qChainKeys) {
            for (auto dbChainKey: dbChainKeys) {
                alignmentLinesMap.erase({qChainKey, dbChainKey});
            }
        }
    }

    static bool getQueryAlnResult(unsigned int qChainKey, const std::vector<unsigned int> &dbChainKeys,  alignmentLinesMap_t &alignmentLinesMap, Matcher::result_t &qAlnResult) {
        for (auto dbChainKey: dbChainKeys) {
            if (alignmentLinesMap[{qChainKey, dbChainKey}].empty()) {
                continue;
            }
            qAlnResult = Matcher::parseAlignmentRecord(alignmentLinesMap[{qChainKey, dbChainKey}].c_str());
            return true;
        }
        return false;
    }
};

class ComplexFilter {
public:
    ComplexFilter(
        const std::vector<unsigned int> &qChainKeys,
        IndexReader *qDbr,
        DBReader<unsigned int> *qStructDbr,
        LocalParameters &par,
        unsigned int thread_idx
    )
        : qChainKeys_(qChainKeys),
          qDbr_(qDbr),
          qStructDbr_(qStructDbr),
          par_(par),
          thread_idx(thread_idx)
    {
        qInterfaceVec_.resize(qChainKeys_.size());
    }

    void computeInterfaceRegion()
    {
        float threshold = INTERFACE_THRESHOLD;
        float d2 = threshold * threshold;
        Coordinate16 coords, coords2;

        for (size_t chainIdx = 0; chainIdx < qChainKeys_.size(); chainIdx++) {
            unsigned int chainKey = qChainKeys_[chainIdx];
            qChainKeyTochainIdx_[chainKey] = chainIdx;

            unsigned int chainDbId = qDbr_->sequenceReader->getId(chainKey);
            char *cadata = qStructDbr_->getData(chainDbId, thread_idx);
            size_t caLength = qStructDbr_->getEntryLen(chainDbId);
            size_t chainLen = qDbr_->sequenceReader->getSeqLen(chainDbId);
            float *chainData = coords.read(cadata, chainLen, caLength);

            for (size_t chainIdx2 = 0; chainIdx2 < qChainKeys_.size(); chainIdx2++) {
                if (chainIdx == chainIdx2) continue;

                unsigned int chainKey2 = qChainKeys_[chainIdx2];
                unsigned int chainDbId2 = qDbr_->sequenceReader->getId(chainKey2);
                char *cadata2 = qStructDbr_->getData(chainDbId2, thread_idx);
                size_t caLength2 = qStructDbr_->getEntryLen(chainDbId2);
                size_t chainLen2 = qDbr_->sequenceReader->getSeqLen(chainDbId2);
                float *chainData2 = coords2.read(cadata2, chainLen2, caLength2);

                for (size_t res1 = 0; res1 < chainLen; res1++) {
                    for (size_t res2 = 0; res2 < chainLen2; res2++) {
                        float dist = BasicFunction::dist(
                            chainData[res1], chainData[chainLen + res1], chainData[2 * chainLen + res1],
                            chainData2[res2], chainData2[chainLen2 + res2], chainData2[2 * chainLen2 + res2]
                        );
                        if (dist < d2) {
                            qInterfaceVec_[chainIdx].push_back(res1);
                            break;
                        }
                    }
                }
            }
        }
    }
    void filterAssignment(
        unsigned int assignmentId,
        Assignment &assignment,
        std::map<unsigned int, std::pair<Assignment, unsigned int>> &tCompBestAssignment,
        chainKeyToComplexId_t &dbChainKeyToComplexIdMap,
        complexIdToChainKeys_t &dbComplexIdToChainKeysMap
    ) 
    {
        unsigned int tComplexId;
        unsigned int qalnlen = 0;
        unsigned int talnlen = 0;
        unsigned int adjustAlnLen;
        assignment.assignmentId = assignmentId;
        assignment.filterResult.clear();
        size_t alnChainNum = assignment.chainToChainResults.size();

        std::vector<int> qChainLengths(alnChainNum);
        std::vector<int> dbChainLengths(alnChainNum);
        std::vector<std::string> backtraceVec(alnChainNum);
        std::vector<int> qStartPosVec(alnChainNum);
        std::vector<int> dbStartPosVec(alnChainNum);
        std::vector<float> qTMscores(alnChainNum);
        std::vector<float> dbTMscores(alnChainNum);
        std::vector<size_t> qidxToqChainKey(alnChainNum);
        std::string uString, tString;

        if (alnChainNum < (size_t)par_.minAlignedChains)
            return;

        unsigned int dbChainNum = 0;

        for (size_t i = 0; i < alnChainNum; i++) {

            const unsigned int &qKey = assignment.chainToChainResults[i].first;
            std::string resultToWrite = assignment.chainToChainResults[i].second+"\t"+assignment.assignmentResult+"\t"+std::to_string(assignment.assignmentId);
            const char* data = resultToWrite.c_str();

            Matcher::result_t res;
            ComplexDataHandler retComplex = parseScoreComplexResult(data, res);

            unsigned int tChainKey = res.dbKey;
            tComplexId = dbChainKeyToComplexIdMap.at(tChainKey);

            qalnlen += abs(res.qEndPos - res.qStartPos) + 1;
            talnlen += abs(res.dbEndPos - res.dbStartPos) + 1;

            qChainLengths[i] = res.qLen;
            dbChainLengths[i] = res.dbLen;
            backtraceVec[i] = res.backtrace;
            qStartPosVec[i] = res.qStartPos;
            dbStartPosVec[i] = res.dbStartPos;

            uString = retComplex.uString;
            tString = retComplex.tString;

            dbChainNum = dbComplexIdToChainKeysMap.at(tComplexId).size();
            qidxToqChainKey[i] = qKey;
        }

        // check if multimer tm matches the threshold parameter
        if (par_.covMode == Parameters::COV_MODE_BIDIRECTIONAL && (assignment.qTmScore < par_.filtMultTmThr || assignment.dbTmScore < par_.filtMultTmThr)) {
            return;
        } else if (par_.covMode == Parameters::COV_MODE_TARGET && assignment.dbTmScore < par_.filtMultTmThr) {
            return;
        } else if (par_.covMode == Parameters::COV_MODE_QUERY && assignment.qTmScore < par_.filtMultTmThr) {
            return;
        }

        // check if multimer coverage matches the threshold parameter
        float qCov = static_cast<float>(qalnlen)/static_cast<float>(assignment.qResidueLength);
        float tCov = static_cast<float>(talnlen)/static_cast<float>(assignment.dbResidueLength);
        if (par_.covMode == Parameters::COV_MODE_BIDIRECTIONAL) {
            adjustAlnLen = (qCov+tCov)/2;
        } else if (par_.covMode == Parameters::COV_MODE_TARGET) {
            adjustAlnLen = tCov;
        } else if (par_.covMode == Parameters::COV_MODE_QUERY) {
            adjustAlnLen = qCov;
        }
        if(! Util::hasCoverage(par_.covThr, par_.covMode, qCov, tCov)){
            return;
        }
        // chain by chain tmscore
        unsigned int matchLen = 0;
        for (unsigned int qidx = 0; qidx < alnChainNum; qidx++) {
            unsigned int qmatchLen = assignment.matchLenVec[qidx];
            Coordinates tmt(qmatchLen), tchain(qmatchLen);
            for (unsigned int residx = matchLen; residx< qmatchLen + matchLen; residx++) {
                tchain.x[residx - matchLen] = assignment.dbCaXVec[residx];
                tchain.y[residx - matchLen] = assignment.dbCaYVec[residx];
                tchain.z[residx - matchLen] = assignment.dbCaZVec[residx];
            }
            // fill t and u
            float t[3];
            float u[3][3];
            std::string tmp;
            int ti = 0;
            const int tlen = static_cast<int>(tString.size());
            for (int k=0; k<tlen; k++) {
                if (k ==tlen-1) {
                    t[ti] = std::stof(tmp);
                } else if (tString[k] == ',') {
                    t[ti] = std::stof(tmp);
                    tmp.clear();
                    ti++;
                } else {
                    tmp.push_back(tString[k]);
                }
            }
            std::string tmp2;
            int ui = 0;
            int uj = 0;
            const int ulen = static_cast<int>(uString.size());
            for (int k=0; k < ulen; k++) {
                if (k==ulen-1) {
                    u[ui][uj] = std::stof(tmp2);
                } else if (uString[k] == ',') {
                    u[ui][uj] = std::stof(tmp2);
                    tmp2.clear();
                    uj++;
                } else {
                    tmp2.push_back(uString[k]);
                }
                if (uj == 3) {
                    ui++;
                    uj = 0;
                }
            }
            // based on t and u, calculate chain tm
            BasicFunction::do_rotation(tchain, tmt, qmatchLen, t, u);
            float d0 = 1.24*(cbrt(dbChainLengths[qidx]-15)) -1.8;
            float d02 = d0*d0;
            float tmScore = 0;
            for (unsigned int ci=0; ci<qmatchLen; ci++) {
                float xa_x = assignment.qCaXVec[matchLen + ci];
                float xa_y = assignment.qCaYVec[matchLen + ci];
                float xa_z = assignment.qCaZVec[matchLen + ci];
                float ya_x = tmt.x[ci];
                float ya_y = tmt.y[ci];
                float ya_z = tmt.z[ci];
                float di = BasicFunction::dist(xa_x, xa_y, xa_z, ya_x, ya_y, ya_z);
                float oneDividedDist = 1/(1+di/d02);
                tmScore += oneDividedDist;
            }
            float qTmScore = tmScore / qChainLengths[qidx];
            float dbTmScore = tmScore / dbChainLengths[qidx];
            qTMscores[qidx] = qTmScore;
            dbTMscores[qidx] = dbTmScore;
            matchLen+=qmatchLen;
        }
        // check chain-tm-threshold 
        // if cov-mode 0, every chain should be aligned, every tm scores should be higher than the threshold
        int chainpassNum = 0;
        if (par_.covMode == Parameters::COV_MODE_BIDIRECTIONAL) {
            if (dbChainNum != qChainKeys_.size() || dbChainNum != alnChainNum) {
                return;
            }
            for (size_t qidx = 0; qidx < alnChainNum; qidx++) {
                if(qTMscores[qidx] < par_.filtChainTmThr) {
                    return;
                }
                if(dbTMscores[qidx] < par_.filtChainTmThr) {
                    return;
                }
            }
        // if cov-mode 1 or 2, min-aligned-chains should have tm higher than the threshold.
        } else if (par_.covMode == Parameters::COV_MODE_TARGET) {
            for (size_t qidx = 0; qidx < alnChainNum; qidx++) {
                if (dbTMscores[qidx]>= par_.filtChainTmThr) {
                chainpassNum++; 
                }
            }
            if(chainpassNum < par_.minAlignedChains) {
                return;
            }
        } else if (par_.covMode == Parameters::COV_MODE_QUERY) {
            for (size_t qidx = 0; qidx < alnChainNum; qidx++) {
                if (qTMscores[qidx]>= par_.filtChainTmThr) {
                chainpassNum++; 
                }
            }
            if(chainpassNum < par_.minAlignedChains) {
                return;
            }
        }
        
        // interface-lddt, only check if aligned chains > 1
        // if the interface lddt parameter is set and aligned chain num ==1, then this assignment doesn't pass
        // if the interface lddt parameter isn't set, interfacelddt is printed out as 0 if aligned chain num ==1
        // if aligned chain num > 1, always calculate
        float interfaceLddt = 0;
        if (alnChainNum == 1 && par_.filtInterfaceLddtThr > 0) {
            return;
        } else if (alnChainNum > 1) {
            std::vector<int> qcaToResidue(assignment.qResidueLength, -1);
            std::vector<int> dbcaToResidue(assignment.dbResidueLength, -1);
            unsigned int numIncrease = 0;
            unsigned int alnChainSum = 0;
            for (size_t qidx = 0; qidx < alnChainNum; qidx++) {
                std::string backtraceString = Matcher::uncompressAlignment(backtraceVec[qidx]);
                unsigned int qStart = qStartPosVec[qidx];
                unsigned int dbStart = dbStartPosVec[qidx];
                unsigned int qPos = 0;
                unsigned int dbPos = 0;
                for (size_t btPos = 0; btPos < backtraceString.size(); btPos++) {
                    if (backtraceString[btPos] == 'M') {
                        qcaToResidue[alnChainSum + qStart + qPos] = numIncrease;
                        dbcaToResidue[alnChainSum + dbStart + dbPos] = numIncrease++;
                        qPos++;
                        dbPos++;
                    } else if (backtraceString[btPos] == 'I') {
                        qPos++;
                    } else if (backtraceString[btPos] == 'D') {
                        dbPos++;
                    }
                }
                alnChainSum += assignment.matchLenVec[qidx];
            }
            std::vector<float> qIntVecX, qIntVecY, qIntVecZ, dbIntVecX, dbIntVecY, dbIntVecZ;
            alnChainSum = 0;
            unsigned int wholeIntLen = 0;
            for (size_t qidx = 0; qidx < alnChainNum; qidx++) {
                std::string backtraceString = Matcher::uncompressAlignment(backtraceVec[qidx]);
                std::vector<int> qchainIsMatch(qChainLengths[qidx], -1);
                std::vector<bool> qchainIsVisited(qChainLengths[qidx], 0);
                unsigned int qStart = qStartPosVec[qidx];
                unsigned int dbStart = dbStartPosVec[qidx];
                unsigned int qPos = 0;
                unsigned int dbPos = 0;
                for (size_t btPos = 0; btPos < backtraceString.size(); btPos++) {
                    if (backtraceString[btPos] == 'M') {
                        qchainIsMatch[qStart + qPos] = dbStart + dbPos;
                        qPos++;
                        dbPos++;
                    } else if (backtraceString[btPos] == 'I') {
                        qPos++;
                    } else if (backtraceString[btPos] == 'D') {
                        dbPos++;
                    }
                }
                std::vector<unsigned int>& qChainInterfaceVec =  qInterfaceVec_[qChainKeyTochainIdx_[qidxToqChainKey[qidx]]];
                wholeIntLen += qChainInterfaceVec.size();
                for (size_t resIdIdx = 0; resIdIdx< qChainInterfaceVec.size(); resIdIdx++) {
                    if(qchainIsMatch[qChainInterfaceVec[resIdIdx]] > -1 && qchainIsVisited[qChainInterfaceVec[resIdIdx]] == 0) {
                        qIntVecX.push_back(assignment.qCaXVec[qcaToResidue[alnChainSum + qChainInterfaceVec[resIdIdx]]]);
                        qIntVecY.push_back(assignment.qCaYVec[qcaToResidue[alnChainSum + qChainInterfaceVec[resIdIdx]]]);
                        qIntVecZ.push_back(assignment.qCaZVec[qcaToResidue[alnChainSum + qChainInterfaceVec[resIdIdx]]]);
                        dbIntVecX.push_back(assignment.dbCaXVec[dbcaToResidue[alnChainSum + qchainIsMatch[qChainInterfaceVec[resIdIdx]]]]);
                        dbIntVecY.push_back(assignment.dbCaYVec[dbcaToResidue[alnChainSum + qchainIsMatch[qChainInterfaceVec[resIdIdx]]]]);
                        dbIntVecZ.push_back(assignment.dbCaZVec[dbcaToResidue[alnChainSum + qchainIsMatch[qChainInterfaceVec[resIdIdx]]]]);
                        qchainIsVisited[qChainInterfaceVec[resIdIdx]] = 1;
                    }
                }
                alnChainSum += assignment.matchLenVec[qidx];
            }
            unsigned int intAlnLen = dbIntVecX.size();
            if (intAlnLen > 0) {
                std::string alnbt(intAlnLen, 'M');
                LDDTCalculator lddtcalculator(intAlnLen+1, intAlnLen+1);
                lddtcalculator.initQuery(intAlnLen, &qIntVecX[0], &qIntVecY[0], &qIntVecZ[0]);
                LDDTCalculator::LDDTScoreResult lddtres = lddtcalculator.computeLDDTScore(intAlnLen, 0, 0, alnbt, &dbIntVecX[0], &dbIntVecY[0], &dbIntVecZ[0]);
                interfaceLddt = lddtres.avgLddtScore * lddtres.scoreLength / wholeIntLen;

                if(interfaceLddt < par_.filtInterfaceLddtThr) {
                    return;
                }
            } else if (par_.filtInterfaceLddtThr > 0) {
                return;
            }
        }
        // write down qcov, tcov, qchaintms, tchaintms, interface lddt if everything passed
        std::string result;
        result.append(SSTR(qCov));
        result.append("\t");
        result.append(SSTR(tCov));
        result.append("\t");
        for (unsigned int i = 0; i < qTMscores.size(); i++) {
            result.append(SSTR(qTMscores[i]));
            if(i < qTMscores.size() - 1) {
                result.append(",");
            }
        }
        result.append("\t");
        for (unsigned int i = 0; i < dbTMscores.size(); i++) {
            result.append(SSTR(dbTMscores[i]));
            if(i < dbTMscores.size() - 1) {
                result.append(",");
            }
        }
        result.append("\t");
        result.append(SSTR(interfaceLddt));
        assignment.filterResult = result;

        auto it = tCompBestAssignment.find(tComplexId);
        if (it == tCompBestAssignment.end() || adjustAlnLen > it->second.second)
            tCompBestAssignment[tComplexId] = {assignment, adjustAlnLen};
    }

private:
    const std::vector<unsigned int> &qChainKeys_;
    IndexReader *qDbr_;
    DBReader<unsigned int> *qStructDbr_;
    LocalParameters &par_;

    std::vector<std::vector<unsigned int>> qInterfaceVec_;
    std::map<unsigned int, unsigned int> qChainKeyTochainIdx_;
    unsigned int thread_idx;
};


static void getlookupInfo(
        IndexReader* dbr,
        const std::string &file,
        std::map<unsigned int, unsigned int> &chainKeyToComplexIdLookup,
        std::map<unsigned int, std::vector<unsigned int>> &complexIdToChainKeysLookup,
        std::vector<unsigned int> &complexIdVec,
        std::map<unsigned int, std::string> &chainKeyToChainNameMap
) {
    if (file.length() == 0) {
        return;
    }
    MemoryMapped lookupDB(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char *data = (char *) lookupDB.getData();
    char *end = data + lookupDB.mappedSize();
    const char *entry[255];
    unsigned int maxSet = 0;
    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        unsigned int complexId = Util::fast_atoi<int>(entry[2]);
        if (complexId > maxSet) {
            maxSet = complexId;
        }
        data = Util::skipLine(data);
    }

    data = (char *) lookupDB.getData();
    end = data + lookupDB.mappedSize();
    std::vector<bool> isVistedSet(maxSet + 1, false);

    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        unsigned int chainKey = Util::fast_atoi<int>(entry[0]);
        unsigned int chainDbId = dbr->sequenceReader->getId(chainKey);
        if (chainDbId != NOT_AVAILABLE_CHAIN_KEY) {
            size_t complexId = Util::fast_atoi<int>(entry[2]);
            chainKeyToComplexIdLookup.emplace(chainKey, complexId);
            std::string chainName(entry[1], (entry[2] - entry[1]) - 1);
            size_t lastUnderscoreIndex = chainName.find_last_of('_');
            std::string complexName = chainName.substr(0, lastUnderscoreIndex);
            chainName = chainName.substr(lastUnderscoreIndex + 1, chainName.size()); // 7soy_1.pdb_A -> A
            chainKeyToChainNameMap.emplace(chainKey, chainName);
            if (isVistedSet[complexId] == 0){
                complexIdToChainKeysLookup.emplace(complexId, std::vector<unsigned int>());
                complexIdVec.emplace_back(complexId);
                isVistedSet[complexId] = 1;
            }
            complexIdToChainKeysLookup.at(complexId).emplace_back(chainKey);
        }
        data = Util::skipLine(data);
    }
    lookupDB.close();
}

int scoremultimer(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbType = alnDbr.getDbtype();
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(dbType);
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
    IndexReader* t3DiDbr = new IndexReader(
        is3DiIdx ? t3DiDbrName : par.db2,
        par.threads,
        needSrc ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
        touch ? IndexReader::PRELOAD_INDEX : 0,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
        needSrc ? "_seq_ss" : "_ss"
    );
    DBReader<unsigned int>* tCaDbr = new DBReader<unsigned int>(
        needSrc? (par.db2 + "_seq_ca").c_str() : (par.db2 + "_ca").c_str(), 
        needSrc? (par.db2 + "_seq_ca.index").c_str() : (par.db2 + "_ca.index").c_str(),
        par.threads,
        DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA
    );
    tCaDbr->open(DBReader<unsigned int>::NOSORT);

    IndexReader* q3DiDbr = NULL;
    DBReader<unsigned int> *qCaDbr = NULL;
    bool sameDB = false;
    if (par.db1 == par.db2) {
        sameDB = true;
        q3DiDbr = t3DiDbr;
        qCaDbr = tCaDbr;
    } else {
        q3DiDbr = new IndexReader(
                StructureUtil::getIndexWithSuffix(par.db1, "_ss"),
                par.threads, IndexReader::SEQUENCES,
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA
        );
        qCaDbr = new DBReader<unsigned int>((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(),
        par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        qCaDbr->open(DBReader<unsigned int>::NOSORT);
    }

    const float minAssignedChainsRatio = par.minAssignedChainsThreshold > MAX_ASSIGNED_CHAIN_RATIO ? MAX_ASSIGNED_CHAIN_RATIO : par.minAssignedChainsThreshold;
    int monomerIncludeMode = par.monomerIncludeMode;
    std::vector<unsigned int> qComplexIndices;
    std::vector<unsigned int> dbComplexIndices;
    chainKeyToComplexId_t qChainKeyToComplexIdMap, dbChainKeyToComplexIdMap;
    complexIdToChainKeys_t dbComplexIdToChainKeysMap, qComplexIdToChainKeysMap;
    chainKeyToChainName_t qChainKeyToChainNameMap, dbChainKeyToChainNameMap;
    std::string qLookupFile = par.db1 + ".lookup";
    std::string dbLookupFile = par.db2 + ".lookup";
    getlookupInfo(q3DiDbr, qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIndices, qChainKeyToChainNameMap);
    if (sameDB) {
        dbChainKeyToComplexIdMap = qChainKeyToComplexIdMap;
        dbComplexIdToChainKeysMap = qComplexIdToChainKeysMap;
        dbComplexIndices = qComplexIndices;
        dbChainKeyToChainNameMap = qChainKeyToChainNameMap;
    } else {
        getlookupInfo(t3DiDbr, dbLookupFile, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, dbComplexIndices, dbChainKeyToChainNameMap);
    }
    // seems not used
    // qChainKeyToChainNameMap.clear();
    // dbChainKeyToChainNameMap.clear();
    Debug::Progress progress(qComplexIndices.size());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
        Coordinate16 qcoords;
        Coordinate16 tcoords;
        Matcher::result_t res;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        resultToWrite_t resultToWrite;
        resultToWrite_t currentResultToWrite;
        alignmentLinesMap_t alignmentLinesMap;
        SearchResult searchResult;
        std::vector<Assignment> assignments;
        std::map<unsigned int, std::pair<Assignment, unsigned int>> tCompBestAssignment;
        ComplexScorer complexScorer(q3DiDbr, t3DiDbr, alnDbr, qCaDbr, tCaDbr, thread_idx, minAssignedChainsRatio, monomerIncludeMode);
#pragma omp for schedule(dynamic, 1)
        // for each q complex
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
            unsigned int qComplexId = qComplexIndices[qCompIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
            if (monomerIncludeMode == SKIP_MONOMERS && qChainKeys.size() < MULTIPLE_CHAINED_COMPLEX) {
                progress.updateProgress();
                continue;
            }
            // read the search file only once
            complexScorer.getSearchResultLinesMap(qChainKeys, alignmentLinesMap);
            if (alignmentLinesMap.empty()) {
                continue;
            }
            // for each db complex
            for (size_t dbId = 0; dbId < dbComplexIndices.size(); dbId++) {
                unsigned int dbComplexId = dbComplexIndices[dbId];
                std::vector<unsigned int> &dbChainKeys = dbComplexIdToChainKeysMap.at(dbComplexId);
                complexScorer.getSearchResultByDbComplex(qComplexId, dbComplexId, qChainKeys, dbChainKeys, alignmentLinesMap, searchResult);
                if (searchResult.alnVec.empty()) {
                    continue;
                }
                complexScorer.getAssignments(searchResult, assignments);
                searchResult.clear();
            }
            SORT_SERIAL(assignments.begin(), assignments.end(), compareAssignment);
            // for each query chain key
            ComplexFilter filter(qChainKeys, q3DiDbr, qCaDbr, par, thread_idx);
            filter.computeInterfaceRegion();
            // for each assignment, filter
            for (unsigned int assignmentId = 0; assignmentId < assignments.size(); assignmentId++){
                filter.filterAssignment(assignmentId, assignments[assignmentId], tCompBestAssignment, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap);
            }
            for (size_t qChainKeyIdx = 0; qChainKeyIdx < qChainKeys.size(); qChainKeyIdx++) {
                resultToWrite.clear();
                unsigned int & qKey = qChainKeys[qChainKeyIdx];
                for(const auto &pair : tCompBestAssignment) {
                    const Assignment &assignment = pair.second.first;
                    currentResultToWrite.clear();
                    assignment.getChainToChainResult(qKey, currentResultToWrite);
                    if (currentResultToWrite.empty()) {
                        continue;
                    }
                    getResult(resultToWrite, currentResultToWrite, assignment);
                }
                //TODO: not filtered
                // for (auto &assignment: assignments) {
                //     assignment.getChainToChainResult(qKey, currentResultToWrite);
                //     if (currentResultToWrite.empty()) {
                //         continue;
                //     }
                //     getResult(resultToWrite, currentResultToWrite, assignment);
                // }
                resultWriter.writeData(resultToWrite.c_str(), resultToWrite.length(), qKey, thread_idx);
            }
            alignmentLinesMap.clear();
            assignments.clear();
            currentResultToWrite.clear();
            resultToWrite.clear();
            tCompBestAssignment.clear();
            progress.updateProgress();
        }
    }
    alnDbr.close();
    delete t3DiDbr;
    tCaDbr->close();
    delete tCaDbr;
    if (!sameDB) {
        delete q3DiDbr;
        qCaDbr->close();
        delete qCaDbr;
    }
    resultWriter.close(false);
    return EXIT_SUCCESS;
}
