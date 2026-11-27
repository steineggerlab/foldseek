#include <cstring>
#include <cstdio>
#include <fstream>
#include <set>

#include "LocalParameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Coordinate16.h"
#include "MultimerUtil.h"
#include "GemmiWrapper.h"
#include "PdbWriter.h"

#ifdef OPENMP
#include <omp.h>
#endif

class KeyIterator {
public:
    virtual ~KeyIterator() {}
    virtual size_t getSize() const = 0;
    virtual std::pair<const unsigned int*, size_t> getDbKeys(size_t index) = 0;
};

class DbKeyIterator : public KeyIterator {
public:
    DbKeyIterator(DBReader<unsigned int>& db) : db(db) {}
    ~DbKeyIterator() {}

    size_t getSize() const override {
        return db.getSize();
    }

    std::pair<const unsigned int*, size_t> getDbKeys(size_t index) override {
        key = db.getDbKey(index);
        return std::make_pair(&key, 1ull);
    }

private:
    DBReader<unsigned int>& db;
    unsigned int key;
};

class MapIterator : public KeyIterator {
public:
    MapIterator(const complexIdToChainKeys_t& map, const std::vector<unsigned int>& complexIndices)
        : map(map), complexIndices(complexIndices) {}
    ~MapIterator() {}

    size_t getSize() const override {
        return complexIndices.size();
    }

    std::pair<const unsigned int*, size_t> getDbKeys(size_t index) override {
        unsigned int key = complexIndices[index];
        const std::vector<unsigned int>& currentKeys = map.at(key);
        return std::make_pair(currentKeys.data(), currentKeys.size());
    }

private:
    const complexIdToChainKeys_t& map;
    const std::vector<unsigned int>& complexIndices;
};

void writeTitle(FILE* handle, const char* headerData, size_t headerLen) {
    int remainingHeader = headerLen;
    fprintf(handle, "TITLE     %.*s\n",  std::min(70, (int)remainingHeader), headerData);
    remainingHeader -= 70;
    int continuation = 2;
    while (remainingHeader > 0) {
        fprintf(handle, "TITLE  % 3d%.*s\n", continuation, std::min(70, (int)remainingHeader), headerData + (headerLen - remainingHeader));
        remainingHeader -= 70;
        continuation++;
    }
}

void writeModelRemark(FILE* handle, int model, const char* headerData, size_t headerLen) {
    int remainingHeader = headerLen;
    fprintf(handle, "REMARK 999 MODEL %4d %.*s\n", model, std::min(58, remainingHeader), headerData);
    remainingHeader -= 58;
    while (remainingHeader > 0) {
        fprintf(handle, "REMARK 999 %.*s\n", std::min(69, remainingHeader), headerData + (headerLen - remainingHeader));
        remainingHeader -= 69;
    }
}

std::string getChainName(const std::string& name) {
    size_t pos = name.find_last_of('_');
    // entries without a chain suffix become chain A
    if (pos == std::string::npos || pos + 1 == name.size()) {
        return "A";
    }
    return name.substr(pos + 1);
}

int convert2pdb(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    int outputMode = par.pdbOutputMode;
    const bool cifOutput = par.pdbOutputFormat == LocalParameters::PDB_OUTPUT_FORMAT_MMCIF;
    const char* extension = cifOutput ? ".cif" : ".pdb";
    int localThreads = par.threads;
    if (outputMode == LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
        localThreads = 1;
    } else {
        localThreads = par.threads;
    }

    int mode = DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX;
    if (outputMode != LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
        mode |= DBReader<unsigned int>::USE_LOOKUP;
    }
    DBReader<unsigned int> db(par.db1.c_str(), par.db1Index.c_str(), localThreads, mode);
    db.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> db_header(par.hdr1.c_str(), par.hdr1Index.c_str(), localThreads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    db_header.open(DBReader<unsigned int>::NOSORT);

    std::string dbCa = par.db1 + "_ca";
    std::string dbCaIndex = par.db1 + "_ca.index";
    DBReader<unsigned int> db_ca(dbCa.c_str(), dbCaIndex.c_str(), localThreads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    db_ca.open(DBReader<unsigned int>::NOSORT);

    std::string lookupFile = par.db1 + ".lookup";
    chainKeyToComplexId_t chainKeyToComplexIdMap;
    complexIdToChainKeys_t complexIdToChainKeysMap;
    std::vector<unsigned int> complexIndices;
    if (outputMode == LocalParameters::PDB_OUTPUT_MODE_COMPLEX) {
        getKeyToIdMapIdToKeysMapIdVec(db, lookupFile, chainKeyToComplexIdMap, complexIdToChainKeysMap, complexIndices);
    }

    Debug(Debug::INFO) << "Start writing file to " << par.db2 << "\n";
    FILE* handle = NULL;
    GemmiWriter multiModelWriter;
    std::ofstream multiModelStream;
    if (outputMode == LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
        if (cifOutput) {
            // all entries share one data block, gemmi serializes it after the loop
            multiModelWriter.reset(Util::remove_extension(FileUtil::baseName(par.db2)), "");
            multiModelStream.open(par.db2.c_str());
            if (multiModelStream.is_open() == false) {
                perror(par.db2.c_str());
                EXIT(EXIT_FAILURE);
            }
        } else {
            handle = fopen(par.db2.c_str(), "w");
            if (handle == NULL) {
                perror(par.db2.c_str());
                EXIT(EXIT_FAILURE);
            }
            if (db.getSize() > 9999) {
                Debug(Debug::WARNING) << "PDB format cannot contain more than 9999 models, use --pdb-output-format 1 to write mmCIF instead\n";
            }
            for (size_t i = 0; i < db.getSize(); ++i) {
                unsigned int headerId = db_header.getId(db.getDbKey(i));
                const char* headerData = db_header.getData(headerId, 0);
                const size_t headerLen = db_header.getEntryLen(headerId) - 2;
                writeModelRemark(handle, (int)(i + 1), headerData, headerLen);
            }
        }
    } else {
        if (mkdir(par.db2.c_str(), 0755) != 0) {
            if (errno != EEXIST) {
                perror(par.db2.c_str());
                EXIT(EXIT_FAILURE);
            }
        }
    }

    size_t overflowEntries = 0;

#pragma omp parallel num_threads(localThreads) reduction(+:overflowEntries)
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        FILE* threadHandle = NULL;
        if (outputMode == LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
            threadHandle = handle;
        }
        // the multi model output shares one structure, every other mode writes one per file
        GemmiWriter threadWriter;
        GemmiWriter* writer = (outputMode == LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL)
                                  ? &multiModelWriter : &threadWriter;
        std::ofstream cifStream;
        Coordinate16 coords;
        std::vector<std::string> chainIds;

        KeyIterator* keyIterator;
        if (outputMode == LocalParameters::PDB_OUTPUT_MODE_COMPLEX) {
            keyIterator = new MapIterator(complexIdToChainKeysMap, complexIndices);
        } else {
            keyIterator = new DbKeyIterator(db);
        }
        const size_t size = keyIterator->getSize();
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < size; ++i) {
            std::pair<const unsigned int*, size_t> keys = keyIterator->getDbKeys(i);
            int serial = 0;
            bool overflow = false;
            chainIds.clear();
            if (outputMode != LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
                unsigned int key = keys.first[0];
                // std::string name = db.getLookupEntryName(key);
                size_t lookupKey = db.getLookupIdByKey(key);
                if (lookupKey == DB_ENTRY_NOT_FOUND) {
                    continue;
                }
                std::string name = db.getLookupEntryName(lookupKey);
                if (outputMode == LocalParameters::PDB_OUTPUT_MODE_COMPLEX) {
                    std::string chain = name.substr(name.find_last_of('_') + 1);
                    if (chain.size() == 0) {
                        // keep the last underscore, if there is no chain
                        name = name.substr(0, name.find_last_of('_') + 1);
                    } else {
                        name = name.substr(0, name.find_last_of('_'));
                    }
                }
                std::string filename = par.db2 + "/" + name + extension;

                unsigned int headerId = db_header.getId(key);
                const char* headerData = db_header.getData(headerId, thread_idx);
                const size_t headerLen = db_header.getEntryLen(headerId) - 2;
                if (cifOutput) {
                    cifStream.open(filename.c_str());
                    if (cifStream.is_open() == false) {
                        perror(filename.c_str());
                        EXIT(EXIT_FAILURE);
                    }
                    writer->reset(name, std::string(headerData, headerLen));
                    writer->addModel(1, "");
                } else {
                    threadHandle = fopen(filename.c_str(), "w");
                    if (threadHandle == NULL) {
                        perror(filename.c_str());
                        EXIT(EXIT_FAILURE);
                    }
                    writeTitle(threadHandle, headerData, headerLen);

                    // PDB has one column for the chain, keep names that fit unique and give the rest unused ones.
                    // Past 62 chains use two characters in columns 21-22 like gemmi.
                    const char* idChars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
                    const size_t idWidth = keys.second > 62 ? 2 : 1;
                    const size_t idCount = idWidth == 1 ? 62 : 62 * 62;
                    overflow |= idWidth > 1;
                    chainIds.assign(keys.second, std::string());
                    std::set<std::string> used;
                    for (size_t j = 0; j < keys.second; ++j) {
                        std::string chain = getChainName(db.getLookupEntryName(db.getLookupIdByKey(keys.first[j])));
                        if (chain.size() <= idWidth && used.insert(chain).second) {
                            chainIds[j] = chain;
                        }
                    }
                    size_t nextId = 0;
                    for (size_t j = 0; j < keys.second; ++j) {
                        if (chainIds[j].empty() == false) {
                            continue;
                        }
                        std::string id;
                        while (nextId < idCount) {
                            id = idWidth == 1 ? std::string(1, idChars[nextId]) : std::string(1, idChars[nextId / 62]) + idChars[nextId % 62];
                            nextId++;
                            if (used.insert(id).second) {
                                chainIds[j] = id;
                                break;
                            }
                        }
                    }
                }
            }

            std::string chainName = "A";
            for (size_t j = 0; j < keys.second; ++j) {
                unsigned int key = keys.first[j];

                unsigned int seqId = db.getId(key);
                const char* seqData = db.getData(seqId, thread_idx);
                const size_t seqLen = std::max(db.getEntryLen(seqId), (size_t)2) - 2;

                unsigned int caId = db_ca.getId(key);
                const char* caData = db_ca.getData(caId, thread_idx);
                const size_t caLen = db_ca.getEntryLen(caId);

                float* ca = coords.read(caData, seqLen, caLen);

                if (outputMode == LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
                    unsigned int headerId = db_header.getId(key);
                    const char* headerData = db_header.getData(headerId, thread_idx);
                    const size_t headerLen = db_header.getEntryLen(headerId) - 2;
                    if (cifOutput) {
                        // mmCIF has no MODEL record, the model number carries the db key
                        writer->addModel((int)key, std::string(headerData, headerLen));
                    } else {
                        fprintf(threadHandle, "MODEL     %4d\n", (int)(i + 1));
                    }
                } else {
                    chainName = getChainName(db.getLookupEntryName(db.getLookupIdByKey(key)));
                }
                if (cifOutput) {
                    // mmCIF keeps the full chain name, the PDB record has room for one character
                    writer->addChain(chainName);
                }
                const std::string chainId = (chainIds.empty() || chainIds[j].empty()) ? chainName.substr(0, 1) : chainIds[j];
                char line[82];
                const char* aa3 = NULL;
                for (size_t j = 0; j < seqLen; ++j) {
                    aa3 = threeLetterResidue(seqData[j]);
                    if (cifOutput) {
                        writer->addCalpha(aa3, (int)(j + 1), ca[j], ca[j + (1 * seqLen)], ca[j + (2 * seqLen)]);
                    } else {
                        writeCalphaAtom(line, ++serial, aa3, chainId.c_str(), (int)(j + 1), ca[j], ca[j + (1 * seqLen)], ca[j + (2 * seqLen)], overflow);
                        fputs(line, threadHandle);
                    }
                }
                if (cifOutput == false && seqLen > 0) {
                    writeTer(line, ++serial, aa3, chainId.c_str(), (int)seqLen);
                    fputs(line, threadHandle);
                }
                if (outputMode == LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL && cifOutput == false) {
                    fprintf(threadHandle, "ENDMDL\n");
                }
            }
            overflowEntries += overflow;
            if (outputMode != LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
                if (cifOutput) {
                    const bool written = writer->writeCif(cifStream);
                    cifStream.close();
                    if (written == false || cifStream.fail()) {
                        Debug(Debug::ERROR) << "Cannot write file " << par.db2 << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                } else {
                    fprintf(threadHandle, "END   \n");
                    fclose(threadHandle);
                }
            }
        }
        delete keyIterator;
    }
    if (overflowEntries > 0) {
        Debug(Debug::WARNING) << overflowEntries << " entries do not fit the PDB columns and were written with hybrid-36 numbers, fewer decimals or two character chain IDs, use --pdb-output-format 1 to write mmCIF instead\n";
    }

    if (outputMode == LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
        if (cifOutput) {
            const bool written = multiModelWriter.writeCif(multiModelStream);
            multiModelStream.close();
            if (written == false || multiModelStream.fail()) {
                Debug(Debug::ERROR) << "Cannot write file " << par.db2 << "\n";
                EXIT(EXIT_FAILURE);
            }
        } else if (fprintf(handle, "END   \n") < 0 || fclose(handle) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << par.db2 << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    db_ca.close();
    db_header.close();
    db.close();

    return EXIT_SUCCESS;
}
