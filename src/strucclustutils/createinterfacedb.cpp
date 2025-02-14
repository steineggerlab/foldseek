#include "DBWriter.h"
#include "Util.h"
#include "LocalParameters.h"
#include "DBReader.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "Coordinate16.h"
#include "tmalign/basic_fun.h"
#include "MultimerUtil.h"
#include "structureto3di.h"
#include "SubstitutionMatrix.h"
#include "GemmiWrapper.h"
#include "PulchraWrapper.h"
#include "LDDT.h"
#include <map>
#include "itoa.h"
#include "MathUtil.h"
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
    if (sameRes / qChainLen > 0.9) {
        resIdx1.clear();
    }
}


//one dimer db as an input, one interface db as an output
int createinterfacedb(int argc, const char **argv, const Command &command) {
    int thread_idx = 0;
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    int dbaccessMode = (DBReader<unsigned int>::USE_INDEX);
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

    // const bool shouldCompress = (par.compressed == true);
    
    DBWriter ssdbw((par.db2 + "_ss").c_str(), (par.db2 + "_ss.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    ssdbw.open();
    DBWriter cadbw((par.db2 + "_ca").c_str(), (par.db2 + "_ca.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
    cadbw.open();
    DBWriter aadbw((par.db2).c_str(), (par.db2 + ".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    aadbw.open();
    std::vector<Vec3> ca, n, c, cb;
    std::vector<char> ami;
    std::vector<char> alphabet3di1, alphabet3di2;
    StructureTo3Di structureTo3Di;
    SubstitutionMatrix mat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    PulchraWrapper pulchra;
    Coordinate16 qcoords;
    Coordinate16 tcoords;
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
            for (size_t i = 0; i < resIdx2.size(); i++) {
                Vec3 caAtom = {tdata[i], tdata[i + tChainLen], tdata[i + tChainLen * 2]};
                ca.push_back(caAtom);
                //TODO: should I store amino acid too?
            }
            //TODO: is it right?
            pulchra.rebuildBackbone(&ca[0],
                                    &n[0],
                                    &c[0],
                                    &ami[0],
                                    resIdx1.size() + resIdx2.size());

            char *states = structureTo3Di.structure2states(ca.data(),
                                                                n.data(),
                                                                c.data(),
                                                                cb.data(),
                                                                resIdx1.size() + resIdx2.size());
            for (size_t pos = 0; pos < resIdx1.size(); pos++) {
                alphabet3di1.push_back(mat.num2aa[static_cast<int>(states[pos])]);
            }
            for (size_t pos = 0; pos < resIdx2.size(); pos++) {
                alphabet3di2.push_back(mat.num2aa[static_cast<int>(states[resIdx1.size() + pos])]);
            }
            ssdbw.writeData(alphabet3di1.data(), alphabet3di1.size(), qChainKey, thread_idx);
            ssdbw.writeData(alphabet3di2.data(), alphabet3di2.size(), tChainKey, thread_idx);
        
        }
        alphabet3di1.clear();
        alphabet3di2.clear();
        ca.clear();
        n.clear();
        c.clear();
        ami.clear();
    }


    return EXIT_SUCCESS;
    //TODO: multithreading
    //TODO: set parameters in LocalParameters.cpp
}