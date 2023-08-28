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

    void setCaData(std::vector<float> &caX, std::vector<float> &caY, std::vector<float> &caZ) {
        caVecX = caX;
        caVecY = caY;
        caVecZ = caZ;
    }
};

struct FeatureVector{
    FeatureVector() {}
    FeatureVector(float u[3][3], const float t[3]) {
        for (size_t i=0; i<3; i++) {
            for (size_t j=0; j<3; j++) {
                features[3*i+j] =  u[i][j];
            }
            features[9+i] = t[i];
        }
        label = 0;
    }
    double features[12];
    int label;

    double getDistance(const FeatureVector &o) {
        double dist = 0;
        for (size_t i=0; i<12; i++) {
            dist += std::pow(features[i] - o.features[i], 2);
        }
        dist = std::sqrt(dist);
        return dist;
    }
};

struct ChainToChainAln {
    ChainToChainAln() {}
    ChainToChainAln(
            Chain &queryChain,
            Chain &targetChain,
            float *queryCaData,
            float *targetCaData,
            Matcher::result_t &alnResult,
            TMaligner::TMscoreResult &tmResult
            ) : qChain(queryChain), dbChain(targetChain) {
        std::vector<float> qCaXVec;
        std::vector<float> qCaYVec;
        std::vector<float> qCaZVec;
        std::vector<float> dbCaXVec;
        std::vector<float> dbCaYVec;
        std::vector<float> dbCaZVec;
        matches = 0;
        unsigned int qPos = alnResult.qStartPos;
        unsigned int dbPos = alnResult.dbStartPos;
        unsigned int qXPos = 0;
        unsigned int qYPos = alnResult.qLen;
        unsigned int qZPos = alnResult.qLen*2;
        unsigned int dbXPos = 0;
        unsigned int dbYPos = alnResult.dbLen;
        unsigned int dbZPos = alnResult.dbLen*2;
        for (char cigar : alnResult.backtrace) {
            switch (cigar) {
                case 'M':
                    matches ++;
                    qCaXVec.emplace_back(queryCaData[qXPos + qPos]);
                    qCaYVec.emplace_back(queryCaData[qYPos + qPos]);
                    qCaZVec.emplace_back(queryCaData[qZPos + qPos++]);
                    dbCaXVec.emplace_back(targetCaData[dbXPos + dbPos]);
                    dbCaYVec.emplace_back(targetCaData[dbYPos + dbPos]);
                    dbCaZVec.emplace_back(targetCaData[dbZPos + dbPos++]);
                    break;
                case 'I':
                    qPos++;
                    break;
                case 'D':
                    dbPos++;
                    break;
//                default:
//                    Debug(Debug::WARNING) << "backtrace ???" << "\n";
//                    break;
            }
        }
        qChain.setCaData(qCaXVec, qCaYVec, qCaZVec);
        dbChain.setCaData(dbCaXVec, dbCaYVec, dbCaZVec);
        alnLength = alnResult.alnLength;
        featVec = FeatureVector(tmResult.u, tmResult.t);
        char buffer[1024];
        size_t len = Matcher::resultToBuffer(buffer, alnResult, true, true, false);
        result.append(buffer, len-1);
    }

    Chain qChain;
    Chain dbChain;
    unsigned int matches;
    unsigned int alnLength;
    FeatureVector featVec;
    std::string result;
};

struct IndexPairWithDist{
    IndexPairWithDist(unsigned int index1, unsigned int index2, double distance) : indexPair({index1, index2}), distance(distance) {}
    std::pair<unsigned int, unsigned int> indexPair;
    double distance;
};

struct Complex {
    Complex() {}
    Complex(std::vector<unsigned int> &chainKeys) : qChainKeys(chainKeys), alnVec({}) {}
    Complex(std::vector<unsigned int> &chainKeys, unsigned int qResidueLen) : qChainKeys(chainKeys), qResidueLen(qResidueLen), alnVec({}) {}
    std::vector<unsigned int> qChainKeys;
    std::vector<unsigned int> dbChainKeys;
    unsigned int qResidueLen;
    unsigned int dbResidueLen;
    std::vector<ChainToChainAln> alnVec;

    void resetDbComplex(std::vector<unsigned int> &chainKeys, unsigned int residueLen) {
        dbChainKeys = chainKeys;
        dbResidueLen = residueLen;
    }

    void filterAlnVec(float compatibleCheckRatio) {
        if (alnVec.empty())
            return;
        std::vector<ChainToChainAln> newAlnVec;
        std::vector<ChainToChainAln> currDbComplexAlnVec;
        ChainToChainAln firstAln = alnVec[0];
        unsigned int dbPrevComplexId = firstAln.dbChain.complexId;
        unsigned int qChainCount = 1;
        unsigned int qPrevChainKey = firstAln.qChain.chainKey;

        for (auto &aln : alnVec) {
            if (aln.dbChain.complexId != dbPrevComplexId) {
                if (qChainCount >= (unsigned int)(qChainKeys.size() * compatibleCheckRatio))
                    newAlnVec.insert(newAlnVec.end(), currDbComplexAlnVec.begin(), currDbComplexAlnVec.end());
                currDbComplexAlnVec.clear();
                dbPrevComplexId = aln.dbChain.complexId;
                qChainCount = 1;
                qPrevChainKey = aln.qChain.chainKey;
                currDbComplexAlnVec.emplace_back(aln);
            } else if (aln.qChain.chainKey != qPrevChainKey) {
                qChainCount ++;
                qPrevChainKey = aln.qChain.chainKey;
                currDbComplexAlnVec.emplace_back(aln);
            } else {
                currDbComplexAlnVec.emplace_back(aln);
            }
        }
        if (qChainCount >= (unsigned int)(qChainKeys.size() * compatibleCheckRatio))
            newAlnVec.insert(newAlnVec.end(), currDbComplexAlnVec.begin(), currDbComplexAlnVec.end());
        alnVec = newAlnVec;
    }

    void normalize() {
        unsigned int length = alnVec.size();
        double mean[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
        double var[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
        double sd[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
        double cv[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
        for (size_t i=0; i < length; i++) {
            for (size_t j=0; j<12; j++) {
                mean[j] += alnVec[i].featVec.features[j] / (double)length;
            }
        }
        for (size_t i=0; i < length; i++) {
            for (size_t j=0; j<12; j++) {
                double value = alnVec[i].featVec.features[j];
                value -= mean[j];
                value *= value;
                var[j] += value/(double)length;
            }
        }
        for (size_t j=0; j<12; j++) {
            sd[j] = std::sqrt(var[j]);
            cv[j] = abs(mean[j]) > 1.0 ? sd[j]/std::abs(mean[j]) : sd[j];
        }
        for (size_t i=0; i < length; i++) {
            for (size_t j=0; j<12; j++) {
                alnVec[i].featVec.features[j] = (cv[j] < 0.1) ? 0.0 : (alnVec[i].featVec.features[j] - mean[j]) / sd[j];
            }
        }
    }
};

struct ComplexToComplexAln {
    ComplexToComplexAln() {}
    ComplexToComplexAln(unsigned int qLength, unsigned int dbLength): qResidueLength(qLength), dbResidueLength(dbLength), matches(0) {}
    std::vector<unsigned int> qChainKeys;
    std::vector<unsigned int> dbChainKeys;
    unsigned int qResidueLength;
    unsigned int dbResidueLength;
    unsigned matches;
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
    std::vector<std::string> alnResVec;

    void appendChainToChainAln(ChainToChainAln &aln) {
        qChainKeys.emplace_back(aln.qChain.chainKey);
        dbChainKeys.emplace_back(aln.dbChain.chainKey);
        matches += aln.matches;
        qCaXVec.insert(qCaXVec.end(), aln.qChain.caVecX.begin(), aln.qChain.caVecX.end());
        qCaYVec.insert(qCaYVec.end(), aln.qChain.caVecY.begin(), aln.qChain.caVecY.end());
        qCaZVec.insert(qCaZVec.end(), aln.qChain.caVecZ.begin(), aln.qChain.caVecZ.end());
        dbCaXVec.insert(dbCaXVec.end(), aln.dbChain.caVecX.begin(), aln.dbChain.caVecX.end());
        dbCaYVec.insert(dbCaYVec.end(), aln.dbChain.caVecY.begin(), aln.dbChain.caVecY.end());
        dbCaZVec.insert(dbCaZVec.end(), aln.dbChain.caVecZ.begin(), aln.dbChain.caVecZ.end());
        alnResVec.emplace_back(aln.result);
    }
    void reset() {
        qChainKeys.clear();
        dbChainKeys.clear();
        matches = 0;
        qCaXVec.clear();
        qCaYVec.clear();
        qCaZVec.clear();
        dbCaXVec.clear();
        dbCaYVec.clear();
        dbCaZVec.clear();
        alnResVec.clear();
        uString.clear();
        tString.clear();
    }

    void getTmScore(TMaligner &tmAligner) {
        std::string backtrace(matches, 'M');
        unsigned int normLen = std::min(qResidueLength, dbResidueLength);
        tmAligner.initQuery(&qCaXVec[0], &qCaYVec[0], &qCaZVec[0], NULL, matches);
        TMaligner::TMscoreResult tmResult= tmAligner.computeTMscore(&dbCaXVec[0], &dbCaYVec[0], &dbCaZVec[0],
                                                                    matches,0,0, backtrace,normLen);
        qTmScore = tmResult.tmscore * normLen / qResidueLength;
        dbTmScore = tmResult.tmscore * normLen / dbResidueLength;
        std::string sep;
        for (size_t i=0; i<3; i++) {
            sep = i==2 ? "" : ",";
            tString.append(std::to_string(tmResult.t[i]) + sep);
            for (size_t j=0; j<3; j++) {
                sep = i==2 && j==2 ? "" : ",";
                uString.append(std::to_string(tmResult.u[i][j]) + sep);
            }
        }
        backtrace.clear();
        sep.clear();

    }

};

bool compareChainToChainAlnByDbComplexId(const ChainToChainAln &first, const ChainToChainAln &second) {
    if (first.dbChain.complexId < second.dbChain.complexId)
        return true;
    if (first.dbChain.complexId > second.dbChain.complexId)
        return false;
    return false;
}

bool compareChainToChainAlnByComplexIdAndChainKey(const ChainToChainAln &first, const ChainToChainAln &second) {
    if (first.qChain.complexId < second.qChain.complexId)
        return true;
    if (first.qChain.complexId > second.qChain.complexId)
        return false;
    if (first.dbChain.complexId < second.dbChain.complexId)
        return true;
    if (first.dbChain.complexId > second.dbChain.complexId)
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

bool compareChainToChainAlnByClusterLabel(const ChainToChainAln &first, const ChainToChainAln &second) {
    if (first.featVec.label < second.featVec.label)
        return true;
    if (first.featVec.label > second.featVec.label)
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

bool compareIndexPairWithDistByDist(const IndexPairWithDist &first, const IndexPairWithDist &second) {
    if (first.distance < second.distance) {
        return true;
    }
    if (first.distance > second.distance) {
        return false;
    }
    return false;
}

bool compareComplexToComplexAln(const ComplexToComplexAln &first, const ComplexToComplexAln &second) {
    if (first.qTmScore > second.qTmScore) {
        return true;
    }
    if (first.qTmScore < second.qTmScore) {
        return false;
    }
    if (first.dbTmScore > second.dbTmScore) {
        return true;
    }
    if (first.dbTmScore < second.dbTmScore) {
        return false;
    }
    return false;
}

class DBSCANCluster {
public:
    DBSCANCluster(unsigned int minSize, float defEps) {
        recursiveNum = 0;
        minClusterSize = minSize;
        defaultEps = defEps;
    }
    void  clusterAlns(Complex &qComplex, float eps, unsigned int clusterSize) {
        initializeAlnLabeling(qComplex);
        if (++recursiveNum > MAX_RECURSIVE_NUM)
            return;
        if (clusterSize < minClusterSize)
            return;
        if (clusterSize==1)
            return labelAlns(qComplex);
        if (clusterSize==2)
            return getAlnPairs(qComplex);
        int cLabel = 0;
        clearClusterVectors();
        for (size_t i=0; i<qComplex.alnVec.size(); i++) {
            ChainToChainAln &centerAln = qComplex.alnVec[i];
            if (centerAln.featVec.label != 0) continue;
            std::vector<unsigned int> neighbors = getNeighbors(qComplex, i, eps);
            if (neighbors.size() < MIN_PTS) {
                centerAln.featVec.label = -1;
                continue;
            }
            centerAln.featVec.label = ++cLabel;
            std::vector<unsigned int> neighborsOfNeighbors = neighbors;
            unsigned int j = 0;
            neighbors.clear();
            while (j < neighborsOfNeighbors.size()) {
                unsigned int q = neighborsOfNeighbors[j++];
                if (i==q) continue;
                ChainToChainAln &neighborAln = qComplex.alnVec[q];
                neighborAln.featVec.label = (neighborAln.featVec.label == 0) || (neighborAln.featVec.label == -1) ? cLabel : neighborAln.featVec.label;
                std::vector<unsigned int> neighborsOfCurrNeighbor = getNeighbors(qComplex, q, eps);
                if (neighborsOfCurrNeighbor.size() >= MIN_PTS) {
                    for (unsigned int &neighbor : neighborsOfCurrNeighbor) {
                        if (std::find(neighborsOfNeighbors.begin(), neighborsOfNeighbors.end(), neighbor) == neighborsOfNeighbors.end())
                            neighborsOfNeighbors.emplace_back(neighbor);
                    }
                }
            }
            std::vector<unsigned int> qFoundChainKeys;
            std::vector<unsigned int> dbFoundChainKeys;
            bool isDefectiveCluster = false;
            for (auto neighborIdx : neighborsOfNeighbors) {
                ChainToChainAln currAln = qComplex.alnVec[neighborIdx];
                unsigned int qChainKey = currAln.qChain.chainKey;
                unsigned int dbChainKey = currAln.dbChain.chainKey;
                bool isNewQChainKey = std::find(qFoundChainKeys.begin(), qFoundChainKeys.end(), qChainKey) == qFoundChainKeys.end();
                bool isNewDbChainKey = std::find(dbFoundChainKeys.begin(), dbFoundChainKeys.end(), dbChainKey)==dbFoundChainKeys.end();
                if (isNewQChainKey && isNewDbChainKey) {
                    qFoundChainKeys.emplace_back(qChainKey);
                    dbFoundChainKeys.emplace_back(dbChainKey);
                } else {
                    isDefectiveCluster = true;
                    break;
                }
            }
            if (isDefectiveCluster)
                defectiveClusters.emplace_back(cLabel);
            else if (neighborsOfNeighbors.size() > clusterSize)
                bigClusters.emplace_back(cLabel);
            else if (neighborsOfNeighbors.size() < clusterSize)
                smallClusters.emplace_back(cLabel);
            else
                validClusters.emplace_back(cLabel);
        }

        if (!validClusters.empty()) {
            keepValidClustersOnly(qComplex);
            SORT_SERIAL(qComplex.alnVec.begin(), qComplex.alnVec.end(), compareChainToChainAlnByClusterLabel);
            return;
        }
        else if (cLabel==0 || !smallClusters.empty()) return clusterAlns(qComplex, eps * (1 + LEARNING_RATE), clusterSize);
        else if (!bigClusters.empty()) return clusterAlns(qComplex, eps * (1 - LEARNING_RATE), clusterSize);
        else if (!defectiveClusters.empty()) return clusterAlns(qComplex, defaultEps, clusterSize - 1);
    }

private:
    const unsigned int MAX_RECURSIVE_NUM = 1000;
    const float LEARNING_RATE = 0.05;
    const unsigned int MIN_PTS = 2;
    float defaultEps;
    unsigned int recursiveNum;
    unsigned int minClusterSize;
    std::vector<unsigned int> validClusters;
    std::vector<unsigned int> smallClusters;
    std::vector<unsigned int> bigClusters;
    std::vector<unsigned int> defectiveClusters;

    static std::vector<unsigned int> getNeighbors(Complex &qComplex, unsigned int centerIdx, float eps) {
        ChainToChainAln centerAln = qComplex.alnVec[centerIdx];
        std::vector<unsigned int> neighbors;
        neighbors.emplace_back(centerIdx);
        for (size_t neighborIdx=0; neighborIdx < qComplex.alnVec.size(); neighborIdx++) {
            if (neighborIdx == centerIdx) continue;
            ChainToChainAln neighborAln = qComplex.alnVec[neighborIdx];
            double dist = centerAln.featVec.getDistance(neighborAln.featVec);
            if (dist<eps)
                neighbors.emplace_back(neighborIdx);
        }
        return neighbors;
    }

    static void labelAlns(Complex &qComplex) {
        int cLabel = 0;
        for (auto & aln : qComplex.alnVec) {
            aln.featVec.label = ++cLabel;
        }
    }

    static void getAlnPairs(Complex &qComplex) {
        int cLabel = 0;
        std::vector<IndexPairWithDist> IndexPairs;
        for (size_t i=0; i<qComplex.alnVec.size(); i++) {
            ChainToChainAln prevAln = qComplex.alnVec[i];
            for (size_t j = i+1; j < qComplex.alnVec.size(); j++) {
                ChainToChainAln currAln = qComplex.alnVec[j];
                if (qComplex.alnVec[i].qChain.chainKey==qComplex.alnVec[j].qChain.chainKey || qComplex.alnVec[i].dbChain.chainKey==qComplex.alnVec[j].dbChain.chainKey) continue;
                double dist = prevAln.featVec.getDistance(currAln.featVec);
                IndexPairs.emplace_back(IndexPairWithDist(i, j, dist));
            }
        }
        SORT_SERIAL(IndexPairs.begin(), IndexPairs.end(), compareIndexPairWithDistByDist);
        if (IndexPairs.empty())
            return labelAlns(qComplex);
        for (auto &IndexPair : IndexPairs) {
            unsigned int alnIdx1 = IndexPair.indexPair.first;
            unsigned int alnIdx2 = IndexPair.indexPair.second;
            if (qComplex.alnVec[alnIdx1].featVec.label > 0 || qComplex.alnVec[alnIdx2].featVec.label > 0) continue;
            qComplex.alnVec[alnIdx1].featVec.label = ++cLabel;
            qComplex.alnVec[alnIdx2].featVec.label = cLabel;
        }
        SORT_SERIAL(qComplex.alnVec.begin(), qComplex.alnVec.end(), compareChainToChainAlnByClusterLabel);
    }

    static void initializeAlnLabeling(Complex &qComplex) {
        for (auto &aln : qComplex.alnVec) {
            aln.featVec.label=0;
        }
    }

    void clearClusterVectors() {
        validClusters.clear();
        smallClusters.clear();
        bigClusters.clear();
        defectiveClusters.clear();
    }

    void keepValidClustersOnly(Complex &qComplex) {
        for (auto &aln : qComplex.alnVec) {
            aln.featVec.label = std::find(validClusters.begin(), validClusters.end(), aln.featVec.label) == validClusters.end() ? -1 : aln.featVec.label;
        }
    }
};

class ComplexScorer {
public:
    ComplexScorer(
            IndexReader *qDbr3Di, IndexReader *tDbr3Di, DBReader<unsigned int> &alnDbr,
            IndexReader *qCaDbr, IndexReader *tCaDbr, unsigned int thread_idx, float minAssignedChainsRatio
            ) : alnDbr(alnDbr), qCaDbr(qCaDbr), tCaDbr(tCaDbr), thread_idx(thread_idx),
        minAssignedChainsRatio(minAssignedChainsRatio) {
        maxChainLen = std::max(qDbr3Di->sequenceReader->getMaxSeqLen() + 1,
                               tDbr3Di->sequenceReader->getMaxSeqLen() + 1);
        q3diDbr = qDbr3Di;
        t3diDbr = tDbr3Di;
        maxResLen = maxChainLen * 2;
        tmAligner = new TMaligner(maxResLen, false, true);
    }

    std::vector<Complex> getQComplexes(
            unsigned int qComplexId, std::vector<unsigned int> &qChainKeys,
            std::map<unsigned int, unsigned int> &dbChainKeyToComplexIdLookup,
            std::map<unsigned int, std::vector<unsigned int>> &dbComplexIdToChainKeysLookup
            ) {
        std::vector<Complex> qComplexes;
        std::vector<ChainToChainAln> parsedAlns;
        unsigned int qResLen = getQueryResidueLength(qChainKeys);
        if (qResLen == 0) {
            return qComplexes;
        }
        Complex qComplex = Complex(qChainKeys, qResLen);
        for (auto qChainKey: qChainKeys) {
            unsigned int qKey = alnDbr.getId(qChainKey);
            if (qKey == NOT_AVAILABLE_CHAIN_KEY)
                continue;
            char *data = alnDbr.getData(qKey, thread_idx);
            if (*data == '\0') continue;
            Matcher::result_t qAlnResult = Matcher::parseAlignmentRecord(data);
            size_t qId = qCaDbr->sequenceReader->getId(qChainKey);
            char *qCaData = qCaDbr->sequenceReader->getData(qId, thread_idx);
            size_t qCaLength = qCaDbr->sequenceReader->getEntryLen(qId);
            float *queryCaData = qCoords.read(qCaData, qAlnResult.qLen, qCaLength);
            Chain qChain = Chain(qComplexId, qChainKey);
            tmAligner->initQuery(queryCaData, &queryCaData[qAlnResult.qLen], &queryCaData[qAlnResult.qLen * 2], NULL,
                                 qAlnResult.qLen);
            while (*data != '\0') {
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const auto dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                const unsigned int dbComplexId = dbChainKeyToComplexIdLookup.at(dbKey);
                Matcher::result_t alnResult = Matcher::parseAlignmentRecord(data);
                size_t tCaId = tCaDbr->sequenceReader->getId(dbKey);
                char *tCaData = tCaDbr->sequenceReader->getData(tCaId, thread_idx);
                size_t tCaLength = tCaDbr->sequenceReader->getEntryLen(tCaId);
                float *targetCaData = tCoords.read(tCaData, alnResult.dbLen, tCaLength);
                Chain dbChain = Chain(dbComplexId, dbKey);
                TMaligner::TMscoreResult tmResult = tmAligner->computeTMscore(targetCaData,
                                                                              &targetCaData[alnResult.dbLen],
                                                                              &targetCaData[alnResult.dbLen +
                                                                                            alnResult.dbLen],
                                                                              alnResult.dbLen, alnResult.qStartPos,
                                                                              alnResult.dbStartPos,
                                                                              Matcher::uncompressAlignment(
                                                                                      alnResult.backtrace),
                                                                              alnResult.alnLength);
                ChainToChainAln chainAln(qChain, dbChain, queryCaData, targetCaData, alnResult, tmResult);
                parsedAlns.emplace_back(chainAln);

                data = Util::skipLine(data);
            }
        }
        if (parsedAlns.empty())
            return qComplexes;
        SORT_SERIAL(parsedAlns.begin(), parsedAlns.end(), compareChainToChainAlnByDbComplexId);
        unsigned int currDbComplexId = parsedAlns[0].dbChain.complexId;
        std::vector<unsigned int> currDbChainKeys = dbComplexIdToChainKeysLookup.at(currDbComplexId);
        unsigned int currDbResLen = getDbResidueLength(currDbChainKeys);
        qComplex.resetDbComplex(currDbChainKeys, currDbResLen);
        for (auto &aln: parsedAlns) {
            if (aln.dbChain.complexId != currDbComplexId) {
                qComplex.filterAlnVec(1.0);
                if (!qComplex.alnVec.empty() && currDbResLen > 0) {
                    qComplex.normalize();
                    qComplexes.emplace_back(qComplex);
                }
                qComplex.alnVec.clear();
                currDbComplexId = aln.dbChain.complexId;
                currDbChainKeys = dbComplexIdToChainKeysLookup.at(currDbComplexId);
                currDbResLen = getDbResidueLength(currDbChainKeys);
                qComplex.resetDbComplex(currDbChainKeys, currDbResLen);
            }
            qComplex.alnVec.emplace_back(aln);
        }
        qComplex.filterAlnVec(1.0);
        if (!qComplex.alnVec.empty() && currDbResLen > 0) {
            qComplex.normalize();
            SORT_SERIAL(qComplex.alnVec.begin(), qComplex.alnVec.end(), compareChainToChainAlnByComplexIdAndChainKey);
            qComplexes.emplace_back(qComplex);
        }
        return qComplexes;
    }

    void getComplexAlns(Complex &qComplex, std::vector<ComplexToComplexAln> &results) {
        if (maxResLen < maxChainLen * qComplex.qChainKeys.size() ||
            maxResLen < maxChainLen * qComplex.dbChainKeys.size()) {
            delete tmAligner;
            maxResLen = std::max(maxChainLen * qComplex.qChainKeys.size(), maxChainLen * qComplex.dbChainKeys.size());
            tmAligner = new TMaligner(maxResLen, false, true);
        }
        DBSCANCluster dbscanCluster = DBSCANCluster(
                (unsigned int) std::ceil(qComplex.qChainKeys.size() * minAssignedChainsRatio), 0.5);
        dbscanCluster.clusterAlns(qComplex, 0.5, std::min(qComplex.qChainKeys.size(), qComplex.dbChainKeys.size()));
        int currLabel = 0;
        ComplexToComplexAln complexAln = ComplexToComplexAln(qComplex.qResidueLen, qComplex.dbResidueLen);
        for (auto &currAln: qComplex.alnVec) {
            if (currAln.featVec.label == 0 || currAln.featVec.label == -1) continue;
            if (currAln.featVec.label == currLabel) {
                complexAln.appendChainToChainAln(currAln);
            } else {
                if (currLabel > 0) {
                    complexAln.getTmScore(*tmAligner);
                    results.emplace_back(complexAln);
                }
                complexAln.reset();
                complexAln.appendChainToChainAln(currAln);
                currLabel = currAln.featVec.label;
            }
        }
        if (currLabel > 0) {
            complexAln.getTmScore(*tmAligner);
            results.emplace_back(complexAln);
        }
    }

    void free() {
        delete tmAligner;
    }

private:
    TMaligner *tmAligner;
    unsigned int maxChainLen;
    DBReader<unsigned int> &alnDbr;
    IndexReader *qCaDbr;
    IndexReader *tCaDbr;
    IndexReader *q3diDbr;
    IndexReader *t3diDbr;
    Coordinate16 qCoords;
    Coordinate16 tCoords;
    unsigned int thread_idx;
    float minAssignedChainsRatio;
    unsigned int maxResLen;

    unsigned int getQueryResidueLength(std::vector<unsigned int> &qChainKeys) {
        unsigned int qResidueLen = 0;
        for (auto qChainKey: qChainKeys) {
            size_t id = q3diDbr->sequenceReader->getId(qChainKey);
            if (id == NOT_AVAILABLE_CHAIN_KEY) {
                // Not accessible
                return 0;
            }
            qResidueLen += q3diDbr->sequenceReader->getSeqLen(id);
        }
        return qResidueLen;
    }

    unsigned int getDbResidueLength(std::vector<unsigned int> &dbChainKeys) {
        unsigned int dbResidueLen = 0;
        for (auto dbChainKey: dbChainKeys) {
            size_t id = t3diDbr->sequenceReader->getId(dbChainKey);
            if (id == NOT_AVAILABLE_CHAIN_KEY) {
                // Not accessible
                return 0;
            }
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
    float minAssignedChainsRatio = par.minAssignedChainsThreshold > 1.0 ? 1.0 : par.minAssignedChainsThreshold;

    std::vector<unsigned int> qComplexIdVec;
    std::vector<unsigned int> dbComplexIdVec;
    std::map<unsigned int, unsigned int> qChainKeyToComplexIdMap;
    std::map<unsigned int, unsigned int> dbChainKeyToComplexIdMap;
    std::map<unsigned int, std::vector<unsigned int>> dbComplexIdToChainKeysMap;
    std::map<unsigned int, std::vector<unsigned int>> qComplexIdToChainKeysMap;
    getKeyToIdMapIdToKeysMapIdVec(qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIdVec);
    getKeyToIdMapIdToKeysMapIdVec(dbLookupFile, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, dbComplexIdVec);
    qChainKeyToComplexIdMap.clear();
    dbComplexIdVec.clear();
    Debug::Progress progress(qComplexIdVec.size());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
        char buffer[1024];
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<ComplexToComplexAln> assignments;
        std::vector<std::string> resultToWriteLines;
        ComplexScorer complexScorer(&q3DiDbr, t3DiDbr, alnDbr, qCaDbr, tCaDbr, thread_idx, minAssignedChainsRatio);
#pragma omp for schedule(dynamic, 1)
        // for each q complex
        for (size_t queryIdx = 0; queryIdx < qComplexIdVec.size(); queryIdx++) {
            progress.updateProgress();
            unsigned int qComplexId = qComplexIdVec[queryIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
            std::vector<Complex> qComplexes = complexScorer.getQComplexes(qComplexId, qChainKeys, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap);
            // for each db complex
            for (size_t qComplexIdx = 0; qComplexIdx < qComplexes.size(); qComplexIdx++) {
                complexScorer.getComplexAlns(qComplexes[qComplexIdx], assignments);
            }
            SORT_SERIAL(assignments.begin(), assignments.end(), compareComplexToComplexAln);
            // for each assignment
            for (unsigned int qChainKeyIdx = 0; qChainKeyIdx < qChainKeys.size(); qChainKeyIdx++) {
                resultToWriteLines.emplace_back("");
            }
            for (unsigned int assignmentId=0; assignmentId < assignments.size(); assignmentId++){
                ComplexToComplexAln &assignment = assignments[assignmentId];
                // for each output line from this assignment
                for (size_t alnIdx=0; alnIdx < assignment.alnResVec.size(); alnIdx++) {
                    unsigned int qKey = assignment.qChainKeys[alnIdx];
                    snprintf(buffer,sizeof(buffer), "%s\t%1.5f\t%1.5f\t%s\t%s\t%d\n", assignment.alnResVec[alnIdx].c_str(), assignment.qTmScore, assignment.dbTmScore, assignment.tString.c_str(), assignment.uString.c_str(), assignmentId);
                    unsigned int currIdx = find(qChainKeys.begin(), qChainKeys.end(), qKey) - qChainKeys.begin();
                    resultToWriteLines[currIdx].append(buffer);
                }
            }
            for (size_t qChainKeyIdx=0; qChainKeyIdx < qChainKeys.size(); qChainKeyIdx++) {
                resultWriter.writeData(resultToWriteLines[qChainKeyIdx].c_str(), resultToWriteLines[qChainKeyIdx].length(), qChainKeys[qChainKeyIdx], thread_idx);
            }
            assignments.clear();
            resultToWriteLines.clear();
        }
        complexScorer.free();
    }
    qComplexIdVec.clear();
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
