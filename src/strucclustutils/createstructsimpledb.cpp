#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "MultimerUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif


void writeMultimerDb(std::vector<unsigned int> &chainKeys, DBReader<unsigned int> &aadb, bool isCompressed, unsigned int complexId, DBWriter &aaWriter) {
    if (chainKeys.size() == 1) {
        unsigned int chainDbId = aadb.getId(chainKeys[0]);
        char *aadata = aadb.getData(chainDbId, 0);
        if (isCompressed) {
            aaWriter.writeData(aadata, *(reinterpret_cast<unsigned int *>(aadata)) + sizeof(unsigned int) + 1,  complexId, 0, false, false);
        } else {
            aaWriter.writeData(aadata, aadb.getEntryLen(chainDbId), complexId, 0, true, false);
        }
        aaWriter.writeIndexEntry(complexId, aaWriter.getStart(0), aadb.getEntryLen(chainDbId), 0);
    } else {
        int totalLength = 0;
        int wholeLength = 0;
        for (size_t chainKey = 0; chainKey < chainKeys.size(); chainKey++) {
            unsigned int chainDbId = aadb.getId(chainKeys[chainKey]);
            size_t entryLength = std::max(aadb.getEntryLen(chainDbId), static_cast<size_t>(1));
            if (chainKey == chainKeys.size() -1) {
                wholeLength += entryLength - 1;
            } else {
                wholeLength += entryLength - 2;
            }
            
        }
        char nullByte = '\0';
        char* result;
        result = new char[wholeLength];
        for (size_t chainKey = 0; chainKey < chainKeys.size(); chainKey++) {
            unsigned int chainDbId = aadb.getId(chainKeys[chainKey]);
            char *aadata = aadb.getData(chainDbId, 0);
            size_t entryLength = std::max(aadb.getEntryLen(chainDbId), static_cast<size_t>(1));
            if (chainKey == chainKeys.size() -1) {
                strncpy(result + totalLength, aadata, entryLength - 1);
                totalLength += entryLength - 1;
            } else {
                strncpy(result + totalLength, aadata, entryLength - 2);
                totalLength += entryLength - 2;
            }
        }
        if (isCompressed) {
            aaWriter.writeData(result, totalLength + 1, complexId, 0, true, false);
        } else {
            aaWriter.writeData(result, totalLength, complexId, 0, false, false);
            aaWriter.writeAdd(&nullByte, sizeof(char), 0);
        }
        aaWriter.writeIndexEntry(complexId, aaWriter.getStart(0), totalLength + 1, 0);  
    }
}

int createstructsimpledb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    DBReader<unsigned int> aadb(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    aadb.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBReader<unsigned int> ssdb((par.db1 + "_ss").c_str(), (par.db1 + "_ss.index").c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    ssdb.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBReader<unsigned int> cadb((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    cadb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int db1type = Parameters::DBTYPE_AMINO_ACIDS;
    int db2type = Parameters::DBTYPE_AMINO_ACIDS;
    int db3type = LocalParameters::DBTYPE_CA_ALPHA;

    DBWriter aaWriter(par.db2.c_str(), par.db2Index.c_str(), 1, par.compressed, db1type);
    aaWriter.open();
    DBWriter ssWriter((par.db2 + "_ss").c_str(), (par.db2 + "_ss.index").c_str(), 1, par.compressed, db2type);
    ssWriter.open();
    DBWriter caWriter((par.db2 + "_ca").c_str(), (par.db2 + "_ca.index").c_str(), 1, par.compressed, db3type);
    caWriter.open();

    chainKeyToComplexId_t chainKeyToComplexIdMap;
    complexIdToChainKeys_t complexIdToChainKeysMap;
    std::vector<unsigned int> complexIndices;
    std::string lookupFile = par.db1 + ".lookup";
    getKeyToIdMapIdToKeysMapIdVec(aadb, lookupFile, chainKeyToComplexIdMap, complexIdToChainKeysMap, complexIndices);
    chainKeyToComplexIdMap.clear();

    const bool isCompressed = aadb.isCompressed();
    Debug::Progress progress(complexIndices.size());


    for (size_t compIdx = 0; compIdx < complexIndices.size(); compIdx++) {
        progress.updateProgress();
        unsigned int complexId = complexIndices[compIdx];
        std::vector<unsigned int> &chainKeys = complexIdToChainKeysMap.at(complexId);
        writeMultimerDb(chainKeys, aadb, isCompressed, complexId, aaWriter);
        writeMultimerDb(chainKeys, ssdb, isCompressed, complexId, ssWriter);
        writeMultimerDb(chainKeys, cadb, isCompressed, complexId, caWriter);
    }   
    complexIdToChainKeysMap.clear();
    complexIndices.clear();
    aadb.close();
    aaWriter.close(true);
    ssdb.close();
    ssWriter.close(true);
    cadb.close();
    caWriter.close(true);

    return EXIT_SUCCESS;
}