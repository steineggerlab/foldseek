#include "DBWriter.h"
#include "Util.h"
#include "LocalParameters.h"
#include "DBReader.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "Coordinate16.h"
#include "tmalign/basic_fun.h"
#include "MultimerUtil.h"
#include "LDDT.h"
#include <map>
#include <set>

#ifdef OPENMP
#include <omp.h>
#endif


void findIntereface(std::vector<size_t> & resIdx1, float squareThreshold, float* qdata, float* tdata, unsigned int qChainLen, unsigned int tChainLen) {
    bool noSameRes = true;
    size_t sameRes = 0;
    
    for (unsigned int qRes = 0; qRes < qChainLen; qRes ++) {
        for (unsigned int tRes = 0; tRes < tChainLen; tRes ++) {
            float distance = MathUtil::squareDist(qdata[qRes], qdata[qChainLen + qRes], qdata[qChainLen * 2 + qRes], tdata[tRes], tdata[tChainLen + tRes], tdata[tChainLen * 2 + tRes]);
            if (distance < 0.01) {
                noSameRes = false;
                sameRes++;
                break;
            }
            if (distance < squareThreshold) {
                resIdx1.push_back(qRes);
                if (noSameRes) {
                    break;
                }
            } 
        }
    }
    if (sameRes / chainLen > 0.9) {
        resIdx1.clear();
    }
}

//one diemr db as an input, one interface db as an output
int createinterfacedb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    IndexReader* qDbr = NULL;
    qDbr = new IndexReader(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    DBReader<unsigned int> qStructDbr((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(), 
                                par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qStructDbr.open(DBReader<unsigned int>::NOSORT);

    std::string qLookupFile = par.db1 + ".lookup";
    std::vector<unsigned int> qComplexIndices;
    chainKeyToComplexId_t qChainKeyToComplexIdMap;
    complexIdToChainKeys_t qComplexIdToChainKeysMap;
    getKeyToIdMapIdToKeysMapIdVec(qDbr, qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIndices);
    qChainKeyToComplexIdMap.clear();

    const bool shouldCompress = (par.compressed == true);
    
    DBWriter ssdbw((outputName+"_ss").c_str(), (outputName+"_ss.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    ssdbw.open();
    DBWriter cadbw((outputName+"_ca").c_str(), (outputName+"_ca.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
    cadbw.open();
    DBWriter aadbw((outputName).c_str(), (outputName+".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    aadbw.open();
    std::vector<Vec3> ca;
    for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
        unsigned int qComplexId = qComplexIndices[qCompIdx];
        std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
        if(qChainKeys.size() != 2) {
            continue;
        }
        //TODO: for each dimer, store only interface c-alphas
        unsigned int qChainIdx = 0;
        unsigned int qChainKey = qChainKeys[qChainIdx];
        unsigned int qChainDbId = qDbr->sequenceReader->getId(qChainKey);
        char *qcadata = qStructDbr.getData(qChainDbId, thread_idx);
        size_t qCaLength = qStructDbr.getEntryLen(qChainDbId);
        size_t qChainLen = qDbr->sequenceReader->getSeqLen(qChainDbId);
        float* qdata = qcoords.read(qcadata, qChainLen, qCaLength);

        unsigned int tChainIdx = 1;
        unsigned int tChainKey = qChainKeys[tChainIdx];
        unsigned int tChainDbId = qDbr->sequenceReader->getId(tChainKey);
        char *tcadata = qStructDbr.getData(tChainDbId, thread_idx);
        size_t tCaLength = qStructDbr.getEntryLen(tChainDbId);
        size_t tChainLen = qDbr->sequenceReader->getSeqLen(tChainDbId);
        float* tdata = qcoords.read(tcadata, tChainLen, tCaLength);
        
        float distanceThreshold = 10;
        const float squareThreshold = distanceThreshold * distanceThreshold;

        std::vector<size_t> resIdx1, resIdx2;
        findIntereface(resIdx1, squareThreshold, qdata, tdata, qChainLen, tChainLen);
        findIntereface(resIdx2, squareThreshold, tdata, qdata, tChainLen, qChainLen);
        if (resIdx1.size() >= 4 && resIdx2.size() >= 4) {
            for (size_t i = 0; i < resIdx1.size(); i++) {
                Vec3 caAtom = {qdata[i], qdata[i + qChainLen], qdata[i + qChainLen * 2]};
                ca.push_back(caAtom);
            }
            //TODO: 2nd chain's interface
            
            // for (size_t i = 0; i < resIdx2.size(); i++) {
            //     Vec3 caAtom = {tdata[i], tdata[i + tChainLen], tdata[i + tChainLen * 2]};
            //     ca.push_back(caAtom);
            // }
        }

        ca.clear();
    }


    //TODO: multithreading
    
    //TODO: rebuildBackbone, and then write them in aa, ca, ss db
    //TODO: header, mapping, lookup, source: do not change
    //TODO: set parameters in LocalParameters.cpp
}