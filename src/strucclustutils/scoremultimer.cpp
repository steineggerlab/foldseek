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

    void standardize(int MonomerIncludeMode) {
        if (dbResidueLen == 0)
            alnVec.clear();

        if (MonomerIncludeMode == SKIP_MONOMERS && dbChainKeys.size() < MULTIPLE_CHAINED_COMPLEX)
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

    void getNeighbors(size_t centerIdx, std::vector<unsigned int> &neighborVec) {
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
        float tmScore;
        unsigned int qKey;
        unsigned int dbKey;
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
            qKey = aln.qChain.chainKey;
            dbKey = aln.dbChain.chainKey;
            tmScore = aln.tmScore;
            qBestTmScore[qKey] = std::max(tmScore, qBestTmScore[qKey]);
            dbBestTmScore[dbKey] = std::max(tmScore, dbBestTmScore[dbKey]);
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
    ComplexScorer(IndexReader *qDbr3Di, IndexReader *tDbr3Di, DBReader<unsigned int> &alnDbr, DBReader<unsigned int> *qCaDbr, DBReader<unsigned int> *tCaDbr, unsigned int thread_idx, float minAssignedChainsRatio, int monomerIncludeMode) : alnDbr(alnDbr), qCaDbr(qCaDbr), tCaDbr(tCaDbr), thread_idx(thread_idx), minAssignedChainsRatio(minAssignedChainsRatio), monomerIncludeMode(monomerIncludeMode)  {
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
            size_t qDbId = qCaDbr->getId(qChainKey);
            char *qCaData = qCaDbr->getData(qDbId, thread_idx);
            size_t qCaLength = qCaDbr->getEntryLen(qDbId);
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
                size_t tCaId = tCaDbr->getId(dbChainKey);
                char *tCaData = tCaDbr->getData(tCaId, thread_idx);
                size_t tCaLength = tCaDbr->getEntryLen(tCaId);
                unsigned int & dbLen = dbAlnResult.dbLen;
                float *targetCaData = tCoords.read(tCaData, dbLen, tCaLength);
                dbChain = Chain(dbComplexId, dbChainKey);
                tmResult = tmAligner->computeTMscore(targetCaData,&targetCaData[dbLen],&targetCaData[dbLen * 2],dbLen,dbAlnResult.qStartPos,dbAlnResult.dbStartPos,Matcher::uncompressAlignment(dbAlnResult.backtrace),dbAlnResult.qLen);
                currAln =  ChainToChainAln(qChain, dbChain, queryCaData, targetCaData, dbAlnResult, tmResult);
                currAlns.emplace_back(currAln);
                currAln.free();
            } // while end
        } // for end
        if (currAlns.empty())
            return;
        // When alignments have no backtrace
        if (!hasBacktrace) {
            Debug(Debug::ERROR) << "Backtraces are required. Please run search with '-a' option.\n";
            EXIT(EXIT_FAILURE);
        }
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
            paredSearchResult.standardize(monomerIncludeMode);
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
        paredSearchResult.standardize(monomerIncludeMode);
        if (!paredSearchResult.alnVec.empty())
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
        DBSCANCluster dbscanCluster(searchResult, finalClusters, minAssignedChainsRatio);
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
    DBReader<unsigned int> *qCaDbr;
    DBReader<unsigned int> *tCaDbr;
    IndexReader *q3diDbr;
    IndexReader *t3diDbr;
    Coordinate16 qCoords;
    Coordinate16 tCoords;
    unsigned int thread_idx;
    float minAssignedChainsRatio;
    unsigned int maxResLen;
    Chain qChain;
    Chain dbChain;
    ChainToChainAln currAln;
    std::vector<ChainToChainAln> currAlns;
    Assignment assignment;
    SearchResult paredSearchResult;
    std::set<cluster_t> finalClusters;
    bool hasBacktrace;
    int monomerIncludeMode;

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

unsigned int cigarToAlignedLength(const std::string &cigar) {
    std::string backtrace = Matcher::uncompressAlignment(cigar);
    unsigned int alni = 0;
    for (size_t btPos = 0; btPos < backtrace.size(); btPos++) {
        if (backtrace[btPos] == 'M') {
            alni++;
        }
    }
    return alni;
}

unsigned int getInterfaceLength(std::vector<unsigned int> &qChainKeys, IndexReader *qDbr, DBReader<unsigned int> *qStructDbr, unsigned int thread_idx, float threshold = INTERFACE_THRESHOLD) {
    float d2 = threshold * threshold;
    Coordinate16 coords, coords2;
    std::set<chainToResidue> local_interface = std::set<chainToResidue>();
    for (size_t chainIdx = 0; chainIdx < qChainKeys.size(); chainIdx++) {
        unsigned int chainKey = qChainKeys[chainIdx];
        unsigned int chainDbId = qDbr->sequenceReader->getId(chainKey);
        char *cadata = qStructDbr->getData(chainDbId, thread_idx);
        size_t caLength = qStructDbr->getEntryLen(chainDbId);
        size_t chainLen = qDbr->sequenceReader->getSeqLen(chainDbId);
        float* chainData = coords.read(cadata, chainLen, caLength);
        
        for (size_t chainIdx2 = chainIdx+1; chainIdx2 < qChainKeys.size(); chainIdx2++) {
            unsigned int chainKey2 = qChainKeys[chainIdx2];
            unsigned int chainDbId2 = qDbr->sequenceReader->getId(chainKey2);
            char *cadata2 = qStructDbr->getData(chainDbId2, thread_idx);
            size_t caLength2 = qStructDbr->getEntryLen(chainDbId2);
            size_t chainLen2 = qDbr->sequenceReader->getSeqLen(chainDbId2);
            float* chainData2 = coords2.read(cadata2, chainLen2, caLength2);
            for (size_t chainResIdx=0; chainResIdx < chainLen; chainResIdx++) {
                bool isInterface = false;
                for (size_t chainResIdx2=0; chainResIdx2 < chainLen2; chainResIdx2++) {
                    float dist = BasicFunction::dist(chainData[chainResIdx], chainData[chainLen + chainResIdx], chainData[2*chainLen + chainResIdx],
                                                    chainData2[chainResIdx2], chainData2[chainLen2 + chainResIdx2], chainData2[2*chainLen2 + chainResIdx2]);
                    if (dist < d2) {
                        isInterface = true;
                        local_interface.insert({chainKey2, chainResIdx2});
                    }
                }
                if (isInterface) {
                    local_interface.insert({chainKey, chainResIdx});
                }
            }
        }
    }
    return local_interface.size();
}

struct Complex {
    int complexId;
    unsigned int nChain;
    unsigned int complexLength;
    std::string complexName;
    // std::vector<unsigned int> chainLengths;
    std::vector<unsigned int> chainKeys;
    
    Complex() : complexId(0), nChain(0), complexLength(0), complexName("") {}
    ~Complex() {
        chainKeys.clear();
        // chainLengths.clear();
    }
};

typedef Coordinates AlignedCoordinate;

unsigned int adjustAlnLen(unsigned int qcov, unsigned int tcov, int covMode) {
    switch (covMode) {
        case Parameters::COV_MODE_BIDIRECTIONAL:
            return (qcov+tcov)/2;
        case Parameters::COV_MODE_TARGET:
            return qcov;
        case Parameters::COV_MODE_QUERY:
            return tcov;
        default:
            return 0;
    }
}

struct chainAlignment {
    unsigned int qKey;
    unsigned int tKey;
    unsigned int qLen;
    unsigned int tLen;
    unsigned int alnLen;
    unsigned int qStartPos;
    unsigned int tStartPos;
    std::string cigar;
    chainAlignment() : qKey(0), tKey(0), qLen(0), tLen(0), alnLen(0), qStartPos(0), tStartPos(0), cigar("") {}
    chainAlignment(unsigned int qKey, unsigned int tKey, unsigned int qLen, unsigned int tLen, unsigned int alnLen, unsigned int qStartPos, unsigned int tStartPos, const std::string &cigar) : 
        qKey(qKey), tKey(tKey), qLen(qLen), tLen(tLen), alnLen(alnLen), qStartPos(qStartPos), tStartPos(tStartPos), cigar(cigar) {}
    ~chainAlignment() {}
};

class ComplexFilterCriteria {
public:
    unsigned int targetComplexId;

    // per complex
    unsigned int qTotalAlnLen;
    unsigned int tTotalAlnLen;
    float qCov;
    float tCov;
    float interfaceLddt;
    float qTm;
    float tTm;
    float avgTm;
    float t[3];
    float u[3][3];

    // per chain : criteria for chainTmThr & lddtThr
    std::vector<float> qAlnChainTms;
    std::vector<float> tAlnChainTms;
    std::vector<chainAlignment> alignedChains;

    ComplexFilterCriteria() {}
    ComplexFilterCriteria(
        unsigned int targetComplexId, float qTm, float tTm, float tstring[3], float ustring[3][3]
    ) :
        targetComplexId(targetComplexId), qTotalAlnLen(0), tTotalAlnLen(0),
        qCov(0), tCov(0), interfaceLddt(0), qTm(qTm), tTm(tTm), avgTm(0)
    {
        std::copy(tstring, tstring + 3, t);
        for (int i = 0; i < 3; i++) {
            std::copy(ustring[i], ustring[i] + 3, u[i]);
        }
    }

    ~ComplexFilterCriteria() {
        qAlnChainTms.clear();
        tAlnChainTms.clear();
        alignedChains.clear();
    }

    bool hasTm(float TmThr, int covMode) {
        switch (covMode) {
            case Parameters::COV_MODE_BIDIRECTIONAL:
                return ((qTm>= TmThr) && (tTm >= TmThr));
            case Parameters::COV_MODE_TARGET:   
                return (tTm >= TmThr);
            case Parameters::COV_MODE_QUERY:
                return (qTm >= TmThr);
            default:
                return true;
        }
    }

    bool hasChainTm(float chainTmThr, int covMode, int minAlignedChains, unsigned int qChainNum, unsigned int tChainNum) {
        int num = 0;
        switch (covMode) {
            case Parameters::COV_MODE_BIDIRECTIONAL:
                if (qAlnChainTms.size()<std::min(qChainNum, tChainNum)) {
                    return false;
                }
                for (size_t i = 0; i < qAlnChainTms.size(); i++) {
                    if (qAlnChainTms[i] < chainTmThr || tAlnChainTms[i] < chainTmThr) {
                        return false;
                    }
                }
                return true;
                break;
            case Parameters::COV_MODE_TARGET:
                for (size_t i = 0; i < tAlnChainTms.size(); i++) {
                    if (tAlnChainTms[i] >= chainTmThr) {
                        num++;
                    }
                }
                if(num >= minAlignedChains) {
                    return true;
                } else {
                    return false;
                }
                break;
            case Parameters::COV_MODE_QUERY:
                for (size_t i = 0; i < qAlnChainTms.size(); i++) {
                    if (qAlnChainTms[i] >= chainTmThr) {
                        num++;
                    }
                }
                if(num >= minAlignedChains) {
                    return true;
                } else {
                    return false;
                }
                break;
            default:
                return false;
        }
    }

    bool hasChainNum(int covMode, unsigned int qChainNum, unsigned int tChainNum) {
        switch (covMode) {
            case Parameters::COV_MODE_BIDIRECTIONAL:
                if (qChainNum != tChainNum) {
                    return false;
                }
                break;
            default:
                return true;
        }
        return true;
    }

    bool hasInterfaceLDDT(float iLddtThr) {
        return(interfaceLddt >= iLddtThr);
    }

    bool hasAlnChainNum(unsigned int minAlignedChains) {
        if (minAlignedChains <= alignedChains.size()) {
            return true;
        }
        return false;
    }

    bool satisfy_first(int covMode, float covThr, float TmThr, int minAlignedChains, unsigned int qChainNum, unsigned int tChainNum) {
        const bool covOK = covThr ? Util::hasCoverage(covThr, covMode, qCov, tCov) : true;
        const bool TmOK = TmThr ? hasTm(TmThr, covMode) : true;
        const bool chainNumOK = hasChainNum(covMode, qChainNum, tChainNum);
        const bool alnChainNumOK = hasAlnChainNum(minAlignedChains);
        return (covOK && TmOK && alnChainNumOK && chainNumOK); 
    }

    bool satisfy_second(int covMode, float chainTmThr, float iLddtThr, int minAlignedChains, unsigned int qChainNum, unsigned int tChainNum) {
        const bool chainTmOK = chainTmThr ? hasChainTm(chainTmThr, covMode, minAlignedChains, qChainNum, tChainNum) : true; 
        const bool lddtOK = iLddtThr ? hasInterfaceLDDT(iLddtThr) : true; 
        return (chainTmOK && lddtOK); 
    }

    void updateAln(unsigned int qAlnLen, unsigned int tAlnLen) {
        qTotalAlnLen += qAlnLen;
        tTotalAlnLen += tAlnLen;
    }

    void computeChainTmScore(AlignedCoordinate &qchain, AlignedCoordinate &tchain, unsigned int totalAlnLen) {
        AlignedCoordinate tmt(totalAlnLen);
        BasicFunction::do_rotation(tchain, tmt, totalAlnLen, t, u);

        unsigned int chainOffset = 0;
        for (unsigned int i=0; i<alignedChains.size(); i++) {
            chainAlignment &chainaln = alignedChains[i];
            unsigned int qLen = chainaln.qLen;
            unsigned int tLen = chainaln.tLen;
            unsigned int alnLen = chainaln.alnLen;
    
            float d0 = 1.24*(cbrt(tLen-15)) -1.8;
            float d02 = d0*d0;

            float tmScore = 0;
            for (unsigned int ci=chainOffset; ci<chainOffset+alnLen; ci++) {
                float xa_x = qchain.x[ci];
                float xa_y = qchain.y[ci];
                float xa_z = qchain.z[ci];
                float ya_x = tmt.x[ci];
                float ya_y = tmt.y[ci];
                float ya_z = tmt.z[ci];
                float di = BasicFunction::dist(xa_x, xa_y, xa_z, ya_x, ya_y, ya_z);
                float oneDividedDist = 1/(1+di/d02);
                tmScore += oneDividedDist;
            }

            float qtmscore = tmScore / qLen;
            float ttmscore = tmScore / tLen;
            updateChainTmScore(qtmscore, ttmscore);
            chainOffset += alnLen;
        }
    }

    void updateChainTmScore(float qChainTm, float tChainTm) {
        qAlnChainTms.push_back(qChainTm);
        tAlnChainTms.push_back(tChainTm);
    }

    void fillComplexAlignment(chainAlignment &alnchain, unsigned int &chainOffset, float *qdata, float *tdata, 
        AlignedCoordinate &qAlnCoords, AlignedCoordinate &tAlnCoords) {
        int mi = chainOffset;
        int qi = alnchain.qStartPos;
        int ti = alnchain.tStartPos;
        int qLen = alnchain.qLen;
        int tLen = alnchain.tLen;
        std::string backtrace = Matcher::uncompressAlignment(alnchain.cigar);
                
        for (size_t btPos = 0; btPos < backtrace.size(); btPos++) {
            if (backtrace[btPos] == 'M') {
                qAlnCoords.x[mi] = qdata[qi];
                qAlnCoords.y[mi] = qdata[qLen + qi];
                qAlnCoords.z[mi] = qdata[2*qLen + qi];
                tAlnCoords.x[mi] = tdata[ti];
                tAlnCoords.y[mi] = tdata[tLen + ti];
                tAlnCoords.z[mi] = tdata[2*tLen + ti];
                qi++;
                ti++;
                mi++;
            }
            else if (backtrace[btPos] == 'I') {
                qi++;
            }
            else {
                ti++;
            }
        }
        chainOffset = mi;
    }

    void calcCov(unsigned int qLen, unsigned int tLen) {
        qCov = static_cast<float>(qTotalAlnLen) / static_cast<float>(qLen);
        tCov = static_cast<float>(tTotalAlnLen) / static_cast<float>(tLen);
    }

    void computeInterfaceLddt(AlignedCoordinate &qAlnCoords, AlignedCoordinate &tAlnCoords, unsigned int interfaceLength, float threshold = INTERFACE_THRESHOLD) {
        if (alignedChains.size() == 1) { // No interface if only one chain aligned
            interfaceLddt = 1;
            return;
        }
        std::vector<unsigned int> chainOffsets(alignedChains.size(), 0);
        unsigned int acc = 0;
        for (size_t i = 0; i < alignedChains.size(); i++) {
            chainOffsets[i] = acc;
            acc += alignedChains[i].alnLen;
        }
        
        float t2 = threshold * threshold;

        std::set<unsigned int> interfacePos;    
        unsigned int intAlnLen = 0;

        // Find and save interface Coordinates
        for (size_t chainIdx = 0; chainIdx < chainOffsets.size(); chainIdx++) {
            unsigned int c1_start = chainOffsets[chainIdx];
            unsigned int c1_end = c1_start + alignedChains[chainIdx].alnLen;
            for (size_t resIdx1 = c1_start; resIdx1 < c1_end; resIdx1++) {
                bool isInterface = false;
                for (size_t resIdx2 = c1_end; resIdx2 < acc; resIdx2++) { // Rest of the chainss
                    float dist = BasicFunction::dist(qAlnCoords.x[resIdx1], qAlnCoords.y[resIdx1], qAlnCoords.z[resIdx1],
                                                    qAlnCoords.x[resIdx2], qAlnCoords.y[resIdx2], qAlnCoords.z[resIdx2]);
                    if (dist < t2) {
                        isInterface = true;
                        if (interfacePos.find(resIdx2) == interfacePos.end()) {
                            interfacePos.insert(resIdx2);
                            intAlnLen++;
                        }
                    }
                }
                if (isInterface && interfacePos.find(resIdx1) == interfacePos.end()) {
                    interfacePos.insert(resIdx1);
                    intAlnLen++;
                }
            }
        }

        if (intAlnLen == 0) {
            return;
        }

        AlignedCoordinate qInterface(intAlnLen);
        AlignedCoordinate tInterface(intAlnLen);
        size_t idx = 0;
        //     // if (qInterfacePos[chainIdx].size() >= 4) { // TODO: Is it important? then change interfacePos into vector. But it can cause (intLen > idx) + downstream errors in lddt calculation
        for (size_t resIdx: interfacePos) {
            qInterface.x[idx] = qAlnCoords.x[resIdx];
            qInterface.y[idx] = qAlnCoords.y[resIdx];
            qInterface.z[idx] = qAlnCoords.z[resIdx];
            tInterface.x[idx] = tAlnCoords.x[resIdx];
            tInterface.y[idx] = tAlnCoords.y[resIdx];
            tInterface.z[idx] = tAlnCoords.z[resIdx];
            idx++;
        }
            // }    

        std::string bt(intAlnLen, 'M');
        LDDTCalculator *lddtcalculator = NULL;
        lddtcalculator = new LDDTCalculator(intAlnLen+1, intAlnLen+1);
        lddtcalculator->initQuery(intAlnLen, &qInterface.x[0], &qInterface.y[0], &qInterface.z[0]);
        LDDTCalculator::LDDTScoreResult lddtres = lddtcalculator->computeLDDTScore(intAlnLen, 0, 0, bt, &tInterface.x[0], &tInterface.y[0], &tInterface.z[0]);
        interfaceLddt = lddtres.avgLddtScore * lddtres.scoreLength / interfaceLength;
        delete lddtcalculator;
    }
};


std::string filterToBuffer(ComplexFilterCriteria cmplfiltcrit , float filtinterfacelddt){
    std::string result;
    result.append(SSTR(cmplfiltcrit.qCov));
    result.append("\t");
    result.append(SSTR(cmplfiltcrit.tCov));
    result.append("\t");

    for (unsigned int i = 0; i < cmplfiltcrit.qAlnChainTms.size(); i++) {
        result.append(SSTR(cmplfiltcrit.qAlnChainTms[i]));
        if(i < cmplfiltcrit.qAlnChainTms.size() - 1) {
            result.append(",");
        }
    }
    if (cmplfiltcrit.qAlnChainTms.size() == 0) {
        result.append(".");
    }
    result.append("\t");
    
    for (unsigned int i = 0; i < cmplfiltcrit.tAlnChainTms.size(); i++) {
        result.append(SSTR(cmplfiltcrit.tAlnChainTms[i]));
        if(i < cmplfiltcrit.tAlnChainTms.size() - 1) {
            result.append(",");
        }
    }
    if (cmplfiltcrit.tAlnChainTms.size() == 0) {
        result.append(".");
    }
    result.append("\t");

    if(filtinterfacelddt == 0){
        result.append(".");
    } else {
        result.append(SSTR(cmplfiltcrit.interfaceLddt));
    }
    result.append("\t");
    return result;
}

void fillUArr(const std::string &uString, float (&u)[3][3]) {
    std::string tmp;
    int i = 0;
    int j=0;
    const int ulen = static_cast<int>(uString.size());
    for (int k=0; k < ulen; k++) {
        if (k==ulen-1) {
            u[i][j] = std::stof(tmp);
        } else if (uString[k] == ',') {
            u[i][j] = std::stof(tmp);
            tmp.clear();
            j++;
        } else {
            tmp.push_back(uString[k]);
        }
        if (j == 3) {
            i++;
            j = 0;
        }
    }
}

void fillTArr(const std::string &tString, float (&t)[3]) {
    std::string tmp;
    int i = 0;
    const int tlen = static_cast<int>(tString.size());
    for (int k=0; k<tlen; k++) {
        if (k ==tlen-1) {
            t[i] = std::stof(tmp);
        } else if (tString[k] == ',') {
            t[i] = std::stof(tmp);
            tmp.clear();
            i++;
        } else {
            tmp.push_back(tString[k]);
        }
    }
}

void getComplexResidueLength( IndexReader *Dbr, std::vector<Complex> &complexes) {
    for (size_t complexIdx = 0; complexIdx < complexes.size(); complexIdx++) {
        Complex *complex = &complexes[complexIdx];
        std::vector<unsigned int> &chainKeys = complex->chainKeys;
        if (chainKeys.empty()) {
            continue;
        }
        unsigned int cmpllen = 0;
        for (auto chainKey: chainKeys) {
            size_t id = Dbr->sequenceReader->getId(chainKey);
            if (id == NOT_AVAILABLE_CHAIN_KEY) {
                break;
            }
            unsigned int reslen = Dbr->sequenceReader->getSeqLen(id);
            // complex->chainLengths.push_back(reslen);
            cmpllen += reslen;
        }
        complex->complexLength = cmpllen;
    }
}

static void getlookupInfo(
        IndexReader* dbr,
        const std::string &file,
        std::map<unsigned int, unsigned int> &chainKeyToComplexIdLookup,
        std::vector<Complex> &complexes,
        std::map<unsigned int, unsigned int> &complexIdtoIdx,
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

    int nComplex = 0;
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
                Complex complex;
                complex.complexId = complexId;
                complex.complexName = complexName;
                complexIdtoIdx.emplace(complexId, nComplex);
                complexes.emplace_back(complex);
                isVistedSet[complexId] = 1;
                nComplex++;
            }
            complexIdToChainKeysLookup.at(complexId).emplace_back(chainKey);
            complexes[complexIdtoIdx.at(complexId)].chainKeys.emplace_back(chainKey);
            complexes[complexIdtoIdx.at(complexId)].nChain++;
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
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(alnDbr.getDbtype());
    int dbType = Parameters::DBTYPE_ALIGNMENT_RES;
    // int dbinfoType = Parameters::DBTYPE_CLUSTER_RES;
    bool needSrc = false;
    if (extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC) {
        needSrc = true;
        dbType = DBReader<unsigned int>::setExtendedDbtype(dbType, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbType);
    resultWriter.open();
    // DBWriter resultinfoWriter((par.db4 + "info").c_str(), (par.db4 + "info.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbinfoType);
    // resultinfoWriter.open();

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    std::string t3DiDbrName =  StructureUtil::getIndexWithSuffix(par.db2, "_ss");
    bool is3DiIdx = Parameters::isEqualDbtype(FileUtil::parseDbType(t3DiDbrName.c_str()), Parameters::DBTYPE_INDEX_DB);
    IndexReader* t3DiDbr = NULL;
    DBReader<unsigned int> *tCaDbr = NULL;
    t3DiDbr = new IndexReader(
            is3DiIdx ? t3DiDbrName : par.db2,
            par.threads,
            needSrc ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
            touch ? IndexReader::PRELOAD_INDEX : 0,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
            needSrc ? "_seq_ss" : "_ss"
    );
    tCaDbr = new DBReader<unsigned int>(
        needSrc? (par.db2 + "_seq_ca").c_str() : (par.db2 + "_ca").c_str(), 
        needSrc? (par.db2 + "_seq_ca.index").c_str() : (par.db2 + "_ca.index").c_str(),
        par.threads,
        DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);

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

    float minAssignedChainsRatio = par.minAssignedChainsThreshold > MAX_ASSIGNED_CHAIN_RATIO ? MAX_ASSIGNED_CHAIN_RATIO: par.minAssignedChainsThreshold;
    int monomerIncludeMode = par.monomerIncludeMode;
    std::vector<Complex> qComplexes, dbComplexes;
    std::map<unsigned int, unsigned int> qComplexIdToIdx, dbComplexIdToIdx;
    std::vector<unsigned int> qComplexIndices;
    std::vector<unsigned int> dbComplexIndices;
    chainKeyToComplexId_t qChainKeyToComplexIdMap, dbChainKeyToComplexIdMap;
    complexIdToChainKeys_t dbComplexIdToChainKeysMap, qComplexIdToChainKeysMap;
    chainKeyToChainName_t qChainKeyToChainNameMap, dbChainKeyToChainNameMap;
    std::string qLookupFile = par.db1 + ".lookup";
    std::string dbLookupFile = par.db2 + ".lookup";
    getlookupInfo(q3DiDbr, qLookupFile, qChainKeyToComplexIdMap, qComplexes, qComplexIdToIdx, qComplexIdToChainKeysMap, qComplexIndices, qChainKeyToChainNameMap);
    getComplexResidueLength(q3DiDbr, qComplexes);
    if (sameDB) {
        dbChainKeyToComplexIdMap = qChainKeyToComplexIdMap;
        dbComplexes = qComplexes;
        dbComplexIdToIdx = qComplexIdToIdx;
        dbComplexIdToChainKeysMap = qComplexIdToChainKeysMap;
        dbComplexIndices = qComplexIndices;
        dbChainKeyToChainNameMap = qChainKeyToChainNameMap;
    } else {
        getlookupInfo(t3DiDbr, dbLookupFile, dbChainKeyToComplexIdMap, dbComplexes, dbComplexIdToIdx, dbComplexIdToChainKeysMap, dbComplexIndices, dbChainKeyToChainNameMap);
        getComplexResidueLength(t3DiDbr, dbComplexes);
    }
    Debug::Progress progress(qComplexIndices.size());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
        char buffer[4096];
        Coordinate16 qcoords;
        Coordinate16 tcoords;
        Matcher::result_t res;
        std::map<unsigned int, ComplexFilterCriteria> localComplexMap;
        std::map<unsigned int, std::vector<unsigned int>> cmplIdToBestAssId;
        std::vector<unsigned int> selectedAssIDs;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<SearchResult> searchResults;
        std::vector<Assignment> assignments;
        std::vector<std::vector<resultToWrite_t>> resultToWriteLines;
        std::vector<resultToWrite_t> resultToWriteLinesFinal;
        ComplexScorer complexScorer(q3DiDbr, t3DiDbr, alnDbr, qCaDbr, tCaDbr, thread_idx, minAssignedChainsRatio, monomerIncludeMode);
#pragma omp for schedule(dynamic, 1)
        // for each q complex
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
            Complex qComplex = qComplexes[qCompIdx];
            unsigned int qComplexId = qComplexIndices[qCompIdx];
            std::map<std::vector<unsigned int>, unsigned int> qalnchain2intlen;
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
            if (monomerIncludeMode == SKIP_MONOMERS && qChainKeys.size() < MULTIPLE_CHAINED_COMPLEX)
                continue;
            complexScorer.getSearchResults(qComplexId, qChainKeys, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, searchResults);
            // for each db complex
            for (size_t dbId = 0; dbId < searchResults.size(); dbId++) {
                complexScorer.getAssignments(searchResults[dbId], assignments);
            }
            SORT_SERIAL(assignments.begin(), assignments.end(), compareAssignment);
            // for each query chain key
            resultToWriteLines.resize(qChainKeys.size());
            resultToWriteLinesFinal.resize(qChainKeys.size());
            // for each assignment
            for (unsigned int assignmentId = 0; assignmentId < assignments.size(); assignmentId++){
                Assignment &assignment = assignments[assignmentId];

                // for each output line from this assignment
                for (size_t resultToWriteIdx = 0; resultToWriteIdx < assignment.resultToWriteLines.size(); resultToWriteIdx++) {
                    unsigned int &qKey = assignment.resultToWriteLines[resultToWriteIdx].first;
                    resultToWrite_t &resultToWrite = assignment.resultToWriteLines[resultToWriteIdx].second;
                    snprintf(buffer, sizeof(buffer), "%s\t%d\n", resultToWrite.c_str(), assignmentId);
                    unsigned int currIdx = find(qChainKeys.begin(), qChainKeys.end(), qKey) - qChainKeys.begin();
                    resultToWriteLines[currIdx].emplace_back(buffer);
                }
            }
            // Writing and reading resultToWriteLines could be revised further. Not efficient now.
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++) {
                std::vector<resultToWrite_t> &resultToWrites = resultToWriteLines[qChainIdx];
                for (size_t resultIdx = 0; resultIdx < resultToWrites.size(); resultIdx++) {
                    resultToWrite_t& resultToWrite = resultToWrites[resultIdx];

                    unsigned int & qChainKey = qChainKeys[qChainIdx];
                    const char* data = resultToWrite.c_str();
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                    unsigned int assId = retComplex.assId;
                    unsigned int tChainKey = res.dbKey;
                    unsigned int tComplexId = dbChainKeyToComplexIdMap.at(tChainKey);
                    unsigned int dbComplexIdx = dbComplexIdToIdx.at(tComplexId);
                    std::vector<unsigned int> tChainKeys = dbComplexes[dbComplexIdx].chainKeys;
                    float u[3][3];
                    float t[3];
                    fillUArr(retComplex.uString, u);
                    fillTArr(retComplex.tString, t);
                    unsigned int qalnlen = (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                    unsigned int talnlen = (std::max(res.dbStartPos, res.dbEndPos) - std::min(res.dbStartPos, res.dbEndPos) + 1);

                    if (localComplexMap.find(assId) == localComplexMap.end()) {
                        ComplexFilterCriteria cmplfiltcrit(tComplexId, retComplex.qTmScore, retComplex.tTmScore, t, u);
                        localComplexMap[assId] = cmplfiltcrit;
                    }
                    ComplexFilterCriteria &cmplfiltcrit = localComplexMap.at(assId);
                    cmplfiltcrit.updateAln(qalnlen, talnlen);
    
                    unsigned int matchLen = cigarToAlignedLength(res.backtrace);
                    chainAlignment chainaln = chainAlignment(qChainKey, tChainKey, res.qLen, res.dbLen, matchLen, res.qStartPos, res.dbStartPos, res.backtrace);
                    cmplfiltcrit.alignedChains.push_back(chainaln);
                }
            }

            for (auto& assId_res : localComplexMap) {
                ComplexFilterCriteria &cmplfiltcrit = assId_res.second;
                unsigned int tComplexId  = cmplfiltcrit.targetComplexId;                
                unsigned int dbComplexIdx = dbComplexIdToIdx.at(tComplexId);
                Complex  &tComplex = dbComplexes[dbComplexIdx];
                cmplfiltcrit.calcCov(qComplex.complexLength, tComplex.complexLength);
                if (!(cmplfiltcrit.satisfy_first(par.covMode, par.covThr, par.tmScoreThr, par.minAlignedChains, qComplex.nChain, tComplex.nChain))) {
                    continue;
                }
                if (par.filtChainTmThr || par.filtInterfaceLddtThr) { // TODO: Recover
                    // Fill aligned coords
                    unsigned int totalAlnLen = 0;
                    for (size_t i = 0; i < cmplfiltcrit.alignedChains.size(); i++) {
                        totalAlnLen += cmplfiltcrit.alignedChains[i].alnLen;
                    }

                    AlignedCoordinate qAlnCoords = AlignedCoordinate(totalAlnLen);
                    AlignedCoordinate tAlnCoords = AlignedCoordinate(totalAlnLen);
                    Coordinate16 qcoords, tcoords;
                    unsigned int chainOffset = 0;
                    
                    for (size_t chainIdx = 0; chainIdx < cmplfiltcrit.alignedChains.size(); chainIdx++) {
                        // Bring Coordinates from cadb
                        chainAlignment &alnchain = cmplfiltcrit.alignedChains[chainIdx];
                        unsigned int qChainKey = alnchain.qKey;
                        unsigned int qChainDbId = q3DiDbr->sequenceReader->getId(qChainKey);
                        char *qcadata = qCaDbr->getData(qChainDbId, thread_idx);
                        size_t qCaLength = qCaDbr->getEntryLen(qChainDbId);
                        size_t qChainLen = q3DiDbr->sequenceReader->getSeqLen(qChainDbId);
                        float* qdata = qcoords.read(qcadata, qChainLen, qCaLength);
                        
                        unsigned int tChainKey = alnchain.tKey;
                        unsigned int tChainDbId = t3DiDbr->sequenceReader->getId(tChainKey);
                        size_t tCaLength = tCaDbr->getEntryLen(tChainDbId);
                        size_t tChainLen = t3DiDbr->sequenceReader->getSeqLen(tChainDbId);
                        char *tcadata = tCaDbr->getData(tChainDbId, thread_idx);
                        float* tdata = tcoords.read(tcadata, tChainLen, tCaLength);

                        // Save each chain into Alignedcoords
                        cmplfiltcrit.fillComplexAlignment(alnchain, chainOffset, qdata, tdata, qAlnCoords, tAlnCoords);
                    }

                    if (par.filtChainTmThr > 0.0) { // TODO: Recover
                        cmplfiltcrit.computeChainTmScore(qAlnCoords, tAlnCoords, totalAlnLen);
                    }

                    if (par.filtInterfaceLddtThr > 0.0) { // TODO: Recover
                        std::vector<unsigned int> qAlnChainKeys(cmplfiltcrit.alignedChains.size());
                        for (size_t i = 0; i < cmplfiltcrit.alignedChains.size(); i++) {
                            qAlnChainKeys[i] = cmplfiltcrit.alignedChains[i].qKey;
                        }
                        sort(qAlnChainKeys.begin(), qAlnChainKeys.end());
                        if (qalnchain2intlen.find(qAlnChainKeys) == qalnchain2intlen.end()) {
                            unsigned int interfaceLength = getInterfaceLength(qAlnChainKeys, q3DiDbr, qCaDbr, thread_idx);
                            qalnchain2intlen[qAlnChainKeys] = interfaceLength;
                        }
                        unsigned int interfaceLength = qalnchain2intlen.at(qAlnChainKeys);

                        cmplfiltcrit.computeInterfaceLddt(qAlnCoords, tAlnCoords, interfaceLength);
                    }

                    if (!(cmplfiltcrit.satisfy_second(par.covMode, par.filtChainTmThr, par.filtInterfaceLddtThr, par.minAlignedChains, qComplex.nChain, tComplex.nChain))) {
                        continue;
                    }
                }

                unsigned int alnlen = adjustAlnLen(cmplfiltcrit.qTotalAlnLen, cmplfiltcrit.tTotalAlnLen, par.covMode);
                
                if (cmplIdToBestAssId.find(tComplexId) == cmplIdToBestAssId.end()) {
                    cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                } else {
                    if (alnlen > cmplIdToBestAssId.at(tComplexId)[1]) {
                        cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                    }
                }
            }

            for (const auto& pair : cmplIdToBestAssId) {
                selectedAssIDs.push_back(pair.second[0]);
            }

            if (selectedAssIDs.size() == 0 && sameDB) {
                float t[3];
                float u[3][3];
                for (int i=0; i < 3; i++) {
                    t[i] = 0.0;
                }
                for (int i=0; i < 3; i++) {
                    for (int j=0; j < 3; j++) {
                        u[i][j] = 0.0;
                    }
                }
                ComplexFilterCriteria cmplfiltcrit(qComplexId, 1.0, 1.0, t, u);
                cmplfiltcrit.qCov = 1.0;
                cmplfiltcrit.tCov = 1.0;
                cmplfiltcrit.interfaceLddt = 1.0;

                selectedAssIDs.push_back(0);
                localComplexMap.insert({0, cmplfiltcrit});
            }
            for (unsigned int assIdidx = 0; assIdidx < selectedAssIDs.size(); assIdidx++) {
                unsigned int assId = selectedAssIDs.at(assIdidx);
                ComplexFilterCriteria &cmplfiltcrit = localComplexMap.at(assId);
                
                std::string result1 = filterToBuffer(cmplfiltcrit, par.filtInterfaceLddtThr);
                Assignment &assignment = assignments[assId];
                for (size_t resultToWriteIdx = 0; resultToWriteIdx < assignment.resultToWriteLines.size(); resultToWriteIdx++) {
                    unsigned int &qKey = assignment.resultToWriteLines[resultToWriteIdx].first;
                    resultToWrite_t &resultToWrite = assignment.resultToWriteLines[resultToWriteIdx].second;
                    
                    snprintf(buffer, sizeof(buffer), "%s\t%s\t%d\n", resultToWrite.c_str(), result1.c_str(), assId);
                    unsigned int currIdx = find(qChainKeys.begin(), qChainKeys.end(), qKey) - qChainKeys.begin();
                    resultToWriteLinesFinal[currIdx].append(buffer);
                }
            }
            for (size_t qChainKeyIdx = 0; qChainKeyIdx < qChainKeys.size(); qChainKeyIdx++) {
                resultToWrite_t &resultToWrite = resultToWriteLinesFinal[qChainKeyIdx];
                unsigned int & qKey = qChainKeys[qChainKeyIdx];
                resultWriter.writeData(resultToWrite.c_str(),resultToWrite.length(),qKey,thread_idx);
            }
            resultToWriteLinesFinal.clear();
            localComplexMap.clear();
            cmplIdToBestAssId.clear();
            selectedAssIDs.clear();

            assignments.clear();
            resultToWriteLines.clear();
            searchResults.clear();
            progress.updateProgress();
        }
        complexScorer.free();
    }

    qComplexIndices.clear();
    dbComplexIndices.clear();
    qChainKeyToComplexIdMap.clear();
    dbChainKeyToComplexIdMap.clear();
    dbComplexIdToChainKeysMap.clear();
    qComplexIdToChainKeysMap.clear();
    alnDbr.close();
    delete t3DiDbr;
    if (!sameDB) {
        delete q3DiDbr;
        delete qCaDbr;
    }
    qComplexes.clear();
    dbComplexes.clear();
    resultWriter.close(false);
    return EXIT_SUCCESS;
}
