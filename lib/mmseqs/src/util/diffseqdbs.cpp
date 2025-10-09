#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <utility>
#include <climits>
#include "Parameters.h"

#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "FastSort.h"

#ifdef OPENMP
#include <omp.h>
#endif


struct compareSecondEntry {
    bool
    operator()(const std::pair<std::string, unsigned int> &lhs, const std::pair<std::string, unsigned int> &rhs) const {
        return (lhs.second < rhs.second);
    }
};

struct compareFirstEntry {
    bool
    operator()(const std::pair<std::string, unsigned int> &lhs, const std::pair<std::string, unsigned int> &rhs) const {
        return (lhs.first < rhs.first) ||
               (lhs.first == rhs.first && lhs.second < rhs.second);
    }
};

struct compareKeyToFirstEntry {
    bool operator()(const std::pair<std::string, unsigned int> &lhs, const std::string &rhs) const {
        return  (lhs.first < rhs);
    }

    bool operator()(const std::string &lhs, const std::pair<std::string, unsigned int> &rhs) const {
        return  (lhs < rhs.first);
    }
};

int diffseqdbs(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> oldReader(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    oldReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> newReader(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    newReader.open(DBReader<unsigned int>::NOSORT);

    std::ofstream removedSeqDBWriter, keptSeqDBWriter, newSeqDBWriter;
    removedSeqDBWriter.open(par.db3);
    keptSeqDBWriter.open(par.db4);
    newSeqDBWriter.open(par.db5);

    // Fill up the hash tables for the old and new DB
    size_t indexSizeOld = oldReader.getSize();
    // key pairs contain (headerID, key) where key is the DB key corresponding to the header
    std::pair<std::string, unsigned int> *keysOld
            = new std::pair<std::string, unsigned int>[indexSizeOld];
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < indexSizeOld; ++id) {
            if (par.useSequenceId) {
                keysOld[id] = std::make_pair(
                        Util::parseFastaHeader(oldReader.getData(id, thread_idx)),
                        oldReader.getDbKey(id)
                );
            } else {
                keysOld[id] = std::make_pair(
                        Util::removeWhiteSpace(oldReader.getData(id, thread_idx)),
                        oldReader.getDbKey(id)
                );
            }
        }
    }

    size_t indexSizeNew = newReader.getSize();
    std::pair<std::string, unsigned int> *keysNew
            = new std::pair<std::string, unsigned int>[indexSizeNew];

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < indexSizeNew; ++id) {
            if (par.useSequenceId) {
                keysNew[id] = std::make_pair(
                        Util::parseFastaHeader(newReader.getData(id,thread_idx)),
                        newReader.getDbKey(id));
            } else {
                keysNew[id] = std::make_pair(
                        Util::removeWhiteSpace(newReader.getData(id, thread_idx)),
                        newReader.getDbKey(id));
            }
        }
    }

    //sort by header for binary search
    SORT_PARALLEL(keysNew, keysNew + indexSizeNew, compareFirstEntry());
    // remove duplicates in new DB by setting the dbkey to UINT_MAX
    for(size_t i = 0; i + 1  < indexSizeNew; ++i) {
        if(keysNew[i].first == keysNew[i+1].first) {
           keysNew[i+1].second = UINT_MAX;
        }
    }
    // default initialized with false
    bool* checkedNew = new bool[indexSizeNew]();
    // doesn't need to be initialized
    size_t *mappedIds = new size_t[indexSizeNew];
    bool* deletedIds = new bool[indexSizeOld]();

    // copy the orignal dbKey from keysOld to originalOldKeys
    unsigned int* originalOldKeys = new unsigned int[indexSizeOld]();
    for (size_t i = 0; i < indexSizeOld; ++i) {
        originalOldKeys[i] = keysOld[i].second;
        keysOld[i].second = i;
    }

    // sorting should be the same as with orignal dbKeys since they are monotonically increasing
    SORT_PARALLEL(keysOld, keysOld + indexSizeOld, compareFirstEntry());
    for (size_t i = 0; i + 1 < indexSizeOld; ++i) {
        if(keysOld[i].first == keysOld[i+1].first) {
            deletedIds[keysOld[i+1].second] = true;
        }
    }
    for (size_t i = 0; i  < indexSizeOld; ++i) {
        keysOld[i].second = originalOldKeys[keysOld[i].second];
    }
    delete [] originalOldKeys;
    // restore original order
    SORT_PARALLEL(keysOld, keysOld + indexSizeOld, compareSecondEntry());

#pragma omp parallel for schedule(dynamic, 10)
    for (size_t id = 0; id < indexSizeOld; ++id) {
        if (deletedIds[id]) {
            continue;
        }
        const std::string &keyToSearch = keysOld[id].first;
        std::pair<std::string, unsigned int> *mappedKey
                = std::lower_bound(keysNew, keysNew + indexSizeNew, keyToSearch, compareKeyToFirstEntry());

        if (mappedKey != (keysNew + indexSizeNew) && keyToSearch.compare(mappedKey->first) == 0) {
            // Found
            size_t indexInNewDB = (mappedKey - keysNew);
            checkedNew[indexInNewDB] = true;
            mappedIds[indexInNewDB] = id;
        } else {
            // Not found
            deletedIds[id] = true;
        }
    }

    for (size_t id = 0; id < indexSizeNew; ++id) {
        if (keysNew[id].second == UINT_MAX) {
            continue;
        }
        if (checkedNew[id]) {
            keptSeqDBWriter << keysOld[mappedIds[id]].second << "\t" << keysNew[id].second << std::endl;
        } else {
            newSeqDBWriter << keysNew[id].second << std::endl;
        }
    }

    for (size_t i = 0; i < indexSizeOld; ++i) {
        if(deletedIds[i]) {
            removedSeqDBWriter << keysOld[i].second << std::endl;
        }
    }
    removedSeqDBWriter.close();


    newSeqDBWriter.close();
    keptSeqDBWriter.close();

    delete[] deletedIds;
    delete[] mappedIds;
    delete[] checkedNew;
    delete[] keysNew;
    delete[] keysOld;
    newReader.close();
    oldReader.close();

    return EXIT_SUCCESS;
}

