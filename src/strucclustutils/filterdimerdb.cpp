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
#include "itoa.h"

#ifdef OPENMP
#include <omp.h>
#endif

//one dimer db as an input, one dimer db with 
int filterdimerdb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    DBReader<unsigned int> qDbr((par.db1).c_str(), (par.db1 + ".index").c_str(), 
                                par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qDbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> qStructDbr((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(), 
                                par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qStructDbr.open(DBReader<unsigned int>::NOSORT);

    std::string qLookupFile = par.db1 + ".lookup";
    std::vector<unsigned int> qComplexIndices;
    chainKeyToComplexId_t qChainKeyToComplexIdMap;
    complexIdToChainKeys_t qComplexIdToChainKeysMap;
    getKeyToIdMapIdToKeysMapIdVec(qDbr, qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIndices);
    qChainKeyToComplexIdMap.clear();


    const int db2Type = Parameters::DBTYPE_GENERIC_DB;
    DBWriter resultWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, 0, db2Type);
    resultWriter.open();

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
    #ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
    #endif
        std::vector<int8_t> camol1, camol2;
        std::vector<float> ca1, ca2;
        Coordinate16 qcoords;
        Coordinate16 tcoords;
#pragma omp for schedule(static)
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
            unsigned int qComplexId = qComplexIndices[qCompIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
            if (qChainKeys.size() != 2) {
                continue;
            }
            unsigned int qChainIdx = 0;
            unsigned int qChainKey = qChainKeys[qChainIdx];
            unsigned int qChainDbId = qDbr.getId(qChainKey);
            char *qcadata = qStructDbr.getData(qChainDbId, thread_idx);
            size_t qCaLength = qStructDbr.getEntryLen(qChainDbId);
            size_t qChainLen = qDbr.getSeqLen(qChainDbId);
            float* qdata = qcoords.read(qcadata, qChainLen, qCaLength);

            unsigned int tChainIdx = 1;
            unsigned int tChainKey = qChainKeys[tChainIdx];
            unsigned int tChainDbId = qDbr.getId(tChainKey);
            char *tcadata = qStructDbr.getData(tChainDbId, thread_idx);
            size_t tCaLength = qStructDbr.getEntryLen(tChainDbId);
            size_t tChainLen = qDbr.getSeqLen(tChainDbId);
            float* tdata = qcoords.read(tcadata, tChainLen, tCaLength);
            
            float distanceThreshold = 10;
            std::vector<size_t> resIdx1, resIdx2;
            const float squareThreshold = distanceThreshold * distanceThreshold;
            findInterface(resIdx1, squareThreshold, qdata, tdata, qChainLen, tChainLen);
            findInterface(resIdx2, squareThreshold, tdata, qdata, tChainLen, qChainLen);
            if (resIdx1.size() >= 4 && resIdx2.size() >= 4) {  
                // std::cout<<qChainKey<<"\t"<<tChainKey<<std::endl;
                resultWriter.writeData("0\n", 2, qChainKey, thread_idx);
                resultWriter.writeData("0\n", 2, tChainKey, thread_idx);
            }
            camol1.clear();
            camol2.clear();
            ca1.clear();
            ca2.clear();
        }
    }
    resultWriter.close(true);
    qStructDbr.close();
    qDbr.close();

    return EXIT_SUCCESS;
}
