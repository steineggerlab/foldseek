#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "structureto3diseqdist.h"
#include "StructureUtil.h"
#include "TMaligner.h"
#include "Coordinate16.h"
#include "MemoryMapped.h"
#include "createcomplexreport.h"

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
    ChainToChainAln(Chain &queryChain, Chain &targetChain, float *qCaData, float *dbCaData, Matcher::result_t &alnResult, TMaligner::TMscoreResult &tmResult) : qChain(queryChain), dbChain(targetChain) {
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
        char buffer[1024];
        label = 0;
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
    double superposition[12];
    unsigned int label;

    double getDistance(const ChainToChainAln &o) {
        double dist = 0;
        for (size_t i=0; i<12; i++) {
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

    void filterAlnVec(double minAlignedQueryChainRatio) {
        if (dbResidueLen == 0)
            alnVec.clear();
        if (alnVec.empty())
            return;
        ChainToChainAln &firstAln = alnVec[0];
        unsigned int minAlignedQueryChainNum = (unsigned int)(qChainKeys.size() * minAlignedQueryChainRatio);
        unsigned int qChainCount = 1;
        unsigned int qPrevChainKey = firstAln.qChain.chainKey;
        for (auto &aln : alnVec) {
            if (aln.qChain.chainKey == qPrevChainKey)
                continue;
            qPrevChainKey = aln.qChain.chainKey;
            if (++qChainCount >= minAlignedQueryChainNum) return;
        }
        alnVec.clear();
    }

    void standardize() {
        if (alnVec.empty())
            return;

        double length = (double) alnVec.size();
        double mean[12] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double var[12] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double sd[12] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double cv[12] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        for (auto &aln: alnVec) {
            for (size_t i = 0; i < 12; i++ ) {
                mean[i] += aln.superposition[i] / length;
            }
        }

        for (auto &aln: alnVec) {
            for (size_t i = 0; i < 12; i++ ) {
                var[i] += std::pow(aln.superposition[i] - mean[i], 2) / length;
            }
        }

        for (size_t i=0; i<12; i++) {
            sd[i] = std::sqrt(var[i]);
            cv[i] = (abs(mean[i]) > TOO_SMALL_MEAN) ? sd[i]/std::abs(mean[i]) : sd[i];
        }

        for (auto &aln: alnVec) {
            for (size_t i = 0; i < 12; i++ ) {
                aln.superposition[i] = cv[i] < TOO_SMALL_CV ? FILTERED_OUT : (aln.superposition[i] - mean[i]) / sd[i];
            }
        }
        return;
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
        resultToWriteLines.emplace_back(resultToWriteWithKey_t(aln.qChain.chainKey, aln.resultToWrite));
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

bool compareChainToChainAlnByDbComplexId(const ChainToChainAln &first, const ChainToChainAln &second) {
    if (first.dbChain.complexId < second.dbChain.complexId)
        return true;
    if (first.dbChain.complexId > second.dbChain.complexId)
        return false;
    return false;
}

bool compareChainToChainAlnByClusterLabel(const ChainToChainAln &first, const ChainToChainAln &second) {
    if (first.label < second.label)
        return true;
    if (first.label > second.label)
        return false;
    if (first.qChain.chainKey < second.qChain.chainKey)
        return true;
    if (first.qChain.chainKey > second.qChain.chainKey)
        return false;
    if (first.dbChain.chainKey < second.dbChain.chainKey)
        return true;
    if (first.dbChain.chainKey > second.dbChain.chainKey)
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

class DBSCANCluster {
public:
    DBSCANCluster(SearchResult &searchResult, double minCov) : searchResult(searchResult) {
        recursiveNum = 0;
        cLabel = 0;
        minClusterSize = (unsigned int) ((double)searchResult.qChainKeys.size() * minCov);
        eps = DEFAULT_EPS;
        idealClusterSize = std::min(searchResult.qChainKeys.size(), searchResult.dbChainKeys.size());
        bestClusters.clear();
        prevMaxClusterSize = 0;
        fillDistMap();
    }

    unsigned int clusterAlns() {
        initializeAlnLabels();
        if (++recursiveNum > MAX_RECURSIVE_NUM) return UNCLUSTERED;
        for (size_t centerAlnIdx=0; centerAlnIdx < searchResult.alnVec.size(); centerAlnIdx++) {
            ChainToChainAln &centerAln = searchResult.alnVec[centerAlnIdx];
            if (centerAln.label != 0) continue;
            getNeighbors(centerAlnIdx, neighbors);
            if (neighbors.size() < MIN_PTS) continue;
            centerAln.label = ++cLabel;
            unsigned int i = 0;
            while (i < neighbors.size()) {
                unsigned int neighborAlnIdx = neighbors[i++];
                if (centerAlnIdx == neighborAlnIdx) continue;
                ChainToChainAln &neighborAln = searchResult.alnVec[neighborAlnIdx];
                neighborAln.label = cLabel;
                getNeighbors(neighborAlnIdx, neighborsOfCurrNeighbor);
                if (neighborsOfCurrNeighbor.size() < MIN_PTS) continue;
                for (auto neighbor : neighborsOfCurrNeighbor) {
                    if (std::find(neighbors.begin(), neighbors.end(), neighbor) == neighbors.end())
                        neighbors.emplace_back(neighbor);
                }
            }
            // too big cluster
            if (neighbors.size() > idealClusterSize )
                continue;
            // redundant chains
            if (checkChainRedundancy())
                continue;
            // too small cluster
            if (neighbors.size() < maxClusterSize)
                continue;
            else if (neighbors.size() == maxClusterSize)
                currClusters.emplace_back(neighbors);
            // new Biggest cluster
            else if (neighbors.size() > maxClusterSize) {
                maxClusterSize = neighbors.size();
                currClusters.clear();
                currClusters.emplace_back(neighbors);
            }
        }

        if (maxClusterSize == idealClusterSize) {
            bestClusters = currClusters;
            return finishDBSCAN();
        } else if (maxClusterSize < prevMaxClusterSize)
            return finishDBSCAN();
        else if (maxClusterSize == prevMaxClusterSize && currClusters.size() < bestClusters.size())
            return finishDBSCAN();
        bestClusters = currClusters;
        prevMaxClusterSize = maxClusterSize;
        eps += LEARNING_RATE;
        return clusterAlns();
    }

private:
    const double LEARNING_RATE = 0.1;
    const double DEFAULT_EPS = 0.1;
    SearchResult &searchResult;
    double eps;
    unsigned int cLabel;
    unsigned int prevMaxClusterSize;
    unsigned int maxClusterSize;
    unsigned int recursiveNum;
    unsigned int idealClusterSize;
    unsigned int minClusterSize;
    std::vector<unsigned int> neighbors;
    std::vector<unsigned int> neighborsOfCurrNeighbor;
    std::vector<unsigned int> qFoundChainKeys;
    std::vector<unsigned int> dbFoundChainKeys;
    distMap_t distMap;
    cluster_t currClusters;
    cluster_t bestClusters;

    void fillDistMap() {
        distMap.clear();
//        float x = -1;
        for (size_t i=0; i < searchResult.alnVec.size(); i++) {
            ChainToChainAln &prevAln = searchResult.alnVec[i];
            for (size_t j = i+1; j < searchResult.alnVec.size(); j++) {
                ChainToChainAln &currAln = searchResult.alnVec[j];
//                if (prevAln.qChain.chainKey == currAln.qChain.chainKey || prevAln.dbChain.chainKey == currAln.dbChain.chainKey)
//                    continue;
                double dist = prevAln.getDistance(currAln);
                distMap.insert({{i,j}, dist});
//                x = x < 0 || x > dist ? dist : x;
            }
        }
//        eps = x;
    }

    void getNeighbors(unsigned int centerIdx, std::vector<unsigned int> &neighborVec) {
        neighborVec.clear();
        neighborVec.emplace_back(centerIdx);
        for (size_t neighborIdx = 0; neighborIdx < searchResult.alnVec.size(); neighborIdx++) {
            double dist = centerIdx < neighborIdx ? distMap[{centerIdx, neighborIdx}] : distMap[{neighborIdx, centerIdx}];
            if (neighborIdx == centerIdx || dist >= eps)
                continue;
            neighborVec.emplace_back(neighborIdx);
        }
    }

    void initializeAlnLabels() {
        for (auto &aln : searchResult.alnVec) {
            aln.label = UNCLUSTERED;
        }
        cLabel = UNCLUSTERED;
        maxClusterSize = 0;
        currClusters.clear();
    }

    bool checkChainRedundancy() {
        qFoundChainKeys.clear();
        dbFoundChainKeys.clear();
        for (auto neighborIdx : neighbors) {
            unsigned int qChainKey = searchResult.alnVec[neighborIdx].qChain.chainKey;
            unsigned int dbChainKey = searchResult.alnVec[neighborIdx].dbChain.chainKey;
            if (std::find(qFoundChainKeys.begin(), qFoundChainKeys.end(), qChainKey) != qFoundChainKeys.end())
                return true;
            if (std::find(dbFoundChainKeys.begin(), dbFoundChainKeys.end(), dbChainKey) != dbFoundChainKeys.end())
                return true;
            qFoundChainKeys.emplace_back(qChainKey);
            dbFoundChainKeys.emplace_back(dbChainKey);
        }
        return false;
    }

    unsigned int finishDBSCAN() {
        initializeAlnLabels();
        if (prevMaxClusterSize < minClusterSize) return UNCLUSTERED;
        cLabel = CLUSTERED;
        for (auto &cluster: bestClusters) {
            for (auto alnIdx: cluster) {
                searchResult.alnVec[alnIdx].label = cLabel;
            }
            cLabel++;
        }
        SORT_SERIAL(searchResult.alnVec.begin(), searchResult.alnVec.end(), compareChainToChainAlnByClusterLabel);
        return CLUSTERED;
    }

//    void free() {
//        neighbors.clear();
//        neighborsOfCurrNeighbor.clear();
//        qFoundChainKeys.clear();
//        dbFoundChainKeys.clear();
//        distMap.clear();
//        currClusters.clear();
//        bestClusters.clear();
//    }
};

class ComplexScorer {
public:
    ComplexScorer(IndexReader *qDbr3Di, IndexReader *tDbr3Di, DBReader<unsigned int> &alnDbr, IndexReader *qCaDbr, IndexReader *tCaDbr, unsigned int thread_idx, double minAssignedChainsRatio) : alnDbr(alnDbr), qCaDbr(qCaDbr), tCaDbr(tCaDbr), thread_idx(thread_idx), minAssignedChainsRatio(minAssignedChainsRatio) {
        maxChainLen = std::max(qDbr3Di->sequenceReader->getMaxSeqLen()+1, tDbr3Di->sequenceReader->getMaxSeqLen()+1);
        q3diDbr = qDbr3Di;
        t3diDbr = tDbr3Di;
        maxResLen = maxChainLen * 2;
        tmAligner = new TMaligner(maxResLen, false, true);
    }

    void getSearchResults(unsigned int qComplexId, std::vector<unsigned int> &qChainKeys, chainKeyToComplexId_t &dbChainKeyToComplexIdLookup, complexIdToChainKeys_t &dbComplexIdToChainKeysLookup, std::vector<SearchResult> &searchResults) {
        searchResults.clear();
        unsigned int qResLen = getQueryResidueLength(qChainKeys);
        if (qResLen == 0) return;
        SearchResult paredSearchResult = SearchResult(qChainKeys, qResLen);
        for (auto qChainKey: qChainKeys) {
            unsigned int qKey = alnDbr.getId(qChainKey);
            if (qKey == NOT_AVAILABLE_CHAIN_KEY) continue;
            char *data = alnDbr.getData(qKey, thread_idx);
            if (*data == '\0') continue;
            qAlnResult = Matcher::parseAlignmentRecord(data);
            size_t qId = qCaDbr->sequenceReader->getId(qChainKey);
            char *qCaData = qCaDbr->sequenceReader->getData(qId, thread_idx);
            size_t qCaLength = qCaDbr->sequenceReader->getEntryLen(qId);
            unsigned int &qLen = qAlnResult.qLen;
            float *queryCaData = qCoords.read(qCaData, qAlnResult.qLen, qCaLength);
            qChain = Chain(qComplexId, qChainKey);
            tmAligner->initQuery(queryCaData, &queryCaData[qLen], &queryCaData[qLen * 2], NULL, qLen);
            while (*data != '\0') {
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const auto dbChainKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                const unsigned int dbComplexId = dbChainKeyToComplexIdLookup.at(dbChainKey);
                dbAlnResult = Matcher::parseAlignmentRecord(data);
                if (dbAlnResult.backtrace.empty()) {
                    Debug(Debug::ERROR) << "Backtraces are required. Please run search with '-a' option.\n";
                    EXIT(EXIT_FAILURE);
                }
                size_t tCaId = tCaDbr->sequenceReader->getId(dbChainKey);
                char *tCaData = tCaDbr->sequenceReader->getData(tCaId, thread_idx);
                size_t tCaLength = tCaDbr->sequenceReader->getEntryLen(tCaId);
                unsigned int & dbLen = dbAlnResult.dbLen;
                float *targetCaData = tCoords.read(tCaData, dbLen, tCaLength);
                dbChain = Chain(dbComplexId, dbChainKey);
                tmResult = tmAligner->computeTMscore(targetCaData,&targetCaData[dbLen],&targetCaData[dbLen * 2],dbLen,dbAlnResult.qStartPos,dbAlnResult.dbStartPos,Matcher::uncompressAlignment(dbAlnResult.backtrace),dbAlnResult.alnLength);
                alnFromCurrentQuery =  ChainToChainAln(qChain, dbChain, queryCaData, targetCaData, dbAlnResult, tmResult);
                alnsFromCurrentQuery.emplace_back(alnFromCurrentQuery);
                alnFromCurrentQuery.free();
                data = Util::skipLine(data);
            }
        }
        if (alnsFromCurrentQuery.empty()) return;
        SORT_SERIAL(alnsFromCurrentQuery.begin(), alnsFromCurrentQuery.end(), compareChainToChainAlnByDbComplexId);
        unsigned int currDbComplexId = alnsFromCurrentQuery[0].dbChain.complexId;
        std::vector<unsigned int> &currDbChainKeys = dbComplexIdToChainKeysLookup.at(currDbComplexId);
        unsigned int currDbResLen = getDbResidueLength(currDbChainKeys);
        paredSearchResult.resetDbComplex(currDbChainKeys, currDbResLen);
        for (auto &aln: alnsFromCurrentQuery) {
            if (aln.dbChain.complexId == currDbComplexId) {
                paredSearchResult.alnVec.emplace_back(aln);
                continue;
            }
            paredSearchResult.filterAlnVec(minAssignedChainsRatio);
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
        alnsFromCurrentQuery.clear();
        paredSearchResult.filterAlnVec(minAssignedChainsRatio);
        paredSearchResult.standardize();
        if (!paredSearchResult.alnVec.empty())
            searchResults.emplace_back(paredSearchResult);
        paredSearchResult.alnVec.clear();
    }

    void getAssignments(SearchResult &searchResult, std::vector<Assignment> &assignments) {
        if (maxResLen < maxChainLen * std::min(searchResult.qChainKeys.size(),  searchResult.dbChainKeys.size())) {
            delete tmAligner;
            maxResLen = std::max(searchResult.qChainKeys.size(), searchResult.dbChainKeys.size()) * maxChainLen;
            tmAligner = new TMaligner(maxResLen, false, true);
        }
        unsigned int currLabel;
        DBSCANCluster dbscanCluster = DBSCANCluster(searchResult,  minAssignedChainsRatio);
        currLabel = dbscanCluster.clusterAlns();
        if (currLabel == UNCLUSTERED) return;
        assignment = Assignment(searchResult.qResidueLen, searchResult.dbResidueLen);
        for (auto &currAln: searchResult.alnVec) {
            if (currAln.label == UNCLUSTERED) continue;
            if (currAln.label != currLabel) {
                assignment.getTmScore(*tmAligner);
                assignment.updateResultToWriteLines();
                assignments.emplace_back(assignment);
                assignment.reset();
                currLabel = currAln.label;
            }
            assignment.appendChainToChainAln(currAln);
        }
        assignment.getTmScore(*tmAligner);
        assignment.updateResultToWriteLines();
        assignments.emplace_back(assignment);
        assignment.reset();
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
    ChainToChainAln alnFromCurrentQuery;
    std::vector<ChainToChainAln> alnsFromCurrentQuery;
    Assignment assignment;

    unsigned int getQueryResidueLength(std::vector<unsigned int> &qChainKeys) {
        unsigned int qResidueLen = 0;
        for (auto qChainKey: qChainKeys) {
            size_t id = q3diDbr->sequenceReader->getId(qChainKey);
            // Not accessible
            if (id == NOT_AVAILABLE_CHAIN_KEY) return 0;
            qResidueLen += q3diDbr->sequenceReader->getSeqLen(id);
        }
        return qResidueLen;
    }

    unsigned int getDbResidueLength(std::vector<unsigned int> &dbChainKeys) {
        unsigned int dbResidueLen = 0;
        for (auto dbChainKey: dbChainKeys) {
            size_t id = t3diDbr->sequenceReader->getId(dbChainKey);
            // Not accessible
            if (id == NOT_AVAILABLE_CHAIN_KEY) return 0;
            dbResidueLen += t3diDbr->sequenceReader->getSeqLen(id);
        }
        return dbResidueLen;
    }
};

int scorecomplex(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader q3DiDbr(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader *t3DiDbr = NULL;
    auto *qCaDbr = new IndexReader(par.db1, par.threads, IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1), touch ? IndexReader::PRELOAD_INDEX : 0, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA, "_ca" );
    IndexReader *tCaDbr = NULL;
    bool sameDB = false;
    if (par.db1 == par.db2) {
        sameDB = true;
        t3DiDbr = &q3DiDbr;
        tCaDbr = qCaDbr;
    } else {
        t3DiDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        tCaDbr = new IndexReader(par.db2, par.threads, IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1), touch ? IndexReader::PRELOAD_INDEX : 0, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA, "_ca");
    }

    std::string qLookupFile = par.db1 + ".lookup";
    std::string dbLookupFile = par.db2 + ".lookup";

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();
    double minAssignedChainsRatio = par.minAssignedChainsThreshold > MAX_ASSIGNED_CHAIN_RATIO ? MAX_ASSIGNED_CHAIN_RATIO: par.minAssignedChainsThreshold;

    std::vector<unsigned int> qComplexIndices;
    std::vector<unsigned int> dbComplexIndices;
    chainKeyToComplexId_t qChainKeyToComplexIdMap;
    chainKeyToComplexId_t dbChainKeyToComplexIdMap;
    complexIdToChainKeys_t dbComplexIdToChainKeysMap;
    complexIdToChainKeys_t qComplexIdToChainKeysMap;
    getKeyToIdMapIdToKeysMapIdVec(qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIndices);
    getKeyToIdMapIdToKeysMapIdVec(dbLookupFile, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, dbComplexIndices);
    qChainKeyToComplexIdMap.clear();
    dbComplexIndices.clear();
    Debug::Progress progress(qComplexIndices.size());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
        char buffer[1024];
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<SearchResult> searchResults;
        std::vector<Assignment> assignments;
        std::vector<resultToWrite_t> resultToWriteLines;
        ComplexScorer complexScorer(&q3DiDbr, t3DiDbr, alnDbr, qCaDbr, tCaDbr, thread_idx, minAssignedChainsRatio);
#pragma omp for schedule(dynamic, 1)
        // for each q complex
        for (size_t qId = 0; qId < qComplexIndices.size(); qId++) {
            progress.updateProgress();
            unsigned int qComplexId = qComplexIndices[qId];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
            complexScorer.getSearchResults(qComplexId, qChainKeys, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, searchResults);
            // for each db complex
            for (size_t dbId = 0; dbId < searchResults.size(); dbId++) {
                complexScorer.getAssignments(searchResults[dbId], assignments);
            }
            SORT_SERIAL(assignments.begin(), assignments.end(), compareAssignment);
            // for each query chain key
            for (unsigned int qChainKeyIdx = 0; qChainKeyIdx < qChainKeys.size(); qChainKeyIdx++) {
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
        }
        complexScorer.free();
    }
    qComplexIndices.clear();
    dbChainKeyToComplexIdMap.clear();
    dbComplexIdToChainKeysMap.clear();
    qComplexIdToChainKeysMap.clear();
    alnDbr.close();
    delete qCaDbr;
    if (!sameDB) {
        delete t3DiDbr;
        delete tCaDbr;
    }
    resultWriter.close(true);
    return EXIT_SUCCESS;
}