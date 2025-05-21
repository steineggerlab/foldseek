#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "MultimerUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

int createstructsimpledb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    DBReader<unsigned int> aadb(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    aadb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbtype = par.dbType;
    // db1Type = Parameters::DBTYPE_AMINO_ACIDS;
    // db3Type = LocalParameters::DBTYPE_CA_ALPHA;

    DBWriter aaWriter(par.db2.c_str(), par.db2Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbtype);
    aaWriter.open();

    chainKeyToComplexId_t chainKeyToComplexIdMap;
    complexIdToChainKeys_t complexIdToChainKeysMap;
    std::vector<unsigned int> complexIndices;
    std::string lookupFile = par.db1 + ".lookup";
    getKeyToIdMapIdToKeysMapIdVec(aadb, lookupFile, chainKeyToComplexIdMap, complexIdToChainKeysMap, complexIndices);
    chainKeyToComplexIdMap.clear();

    char newLine = '\n';
    char nullByte = '\0';
    char* result;
    const bool isCompressed = aadb.isCompressed();
    //TODO: multithreading
    unsigned int thread_idx = 0;

    for (size_t compIdx = 0; compIdx < complexIndices.size(); compIdx++) {
        unsigned int complexId = complexIndices[compIdx];
        std::vector<unsigned int> &chainKeys = complexIdToChainKeysMap.at(complexId);
        if (chainKeys.size() == 1) {
            //TODO: write directly
            unsigned int chainDbId = aadb.getId(chainKeys[0]);
            char *aadata = aadb.getData(chainDbId, thread_idx);
            if (isCompressed) {
                aaWriter.writeData(aadata, *(reinterpret_cast<unsigned int *>(aadata)) + sizeof(unsigned int) + 1,  complexId, thread_idx, false, false);
            } else {
                aaWriter.writeData(aadata, aadb.getEntryLen(chainDbId), complexId, thread_idx, true, false);
            }
            aaWriter.writeIndexEntry(complexId, aaWriter.getStart(0), aadb.getEntryLen(chainDbId), thread_idx);
        } else {
            int totalLength = 0;
            for (size_t chainKey = 0; chainKey < chainKeys.size(); chainKey++) {
                unsigned int chainDbId = aadb.getId(chainKeys[chainKey]);
                char *aadata = aadb.getData(chainDbId, thread_idx);
                size_t entryLength = std::max(aadb.getEntryLen(chainDbId), static_cast<size_t>(1));
                result = new char[entryLength];
                strncpy(result + totalLength, aadata, entryLength - 1);
                totalLength += entryLength - 1;
            }
            result[totalLength] = newLine;
            if (isCompressed) {
                aaWriter.writeData(result, totalLength + 1, complexId, thread_idx, true, false);
            } else {
                aaWriter.writeData(result, totalLength, complexId, thread_idx, false, false);
                aaWriter.writeAdd(&newLine, sizeof(char), thread_idx);
                aaWriter.writeAdd(&nullByte, sizeof(char), thread_idx);
            }
            delete [] result;
            result = nullptr;
            aaWriter.writeIndexEntry(complexId, aaWriter.getStart(0), totalLength + 2, thread_idx);  
        }
            
    }

    return EXIT_SUCCESS;
}