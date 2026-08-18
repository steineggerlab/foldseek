#include <cstring>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "LocalParameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Coordinate16.h"
#include "MultimerUtil.h"

// gemmi throws, so this single file is compiled with -fexceptions
// (see set_source_files_properties in src/CMakeLists.txt)
#define GEMMI_WRITE_IMPLEMENTATION
#include "to_mmcif.hpp"
#include "to_cif.hpp"
#include "polyheur.hpp"

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

// gemmi prints coordinates with the shortest round trip representation of the
// double, so a float 3.8 would end up as 3.79999995. The source files and the
// PDB output carry three decimals, keep the same precision here.
static double roundToPdbPrecision(float value) {
    return std::round((double)value * 1000.0) / 1000.0;
}

// The defaults would emit a bogus 1x1x1 unit cell and an empty _symmetry for a
// coordinate only structure, everything else gemmi writes is real information.
static gemmi::MmcifOutputGroups calphaOutputGroups() {
    gemmi::MmcifOutputGroups groups(true);
    groups.cell = false;
    groups.symmetry = false;
    groups.group_pdb = true;
    return groups;
}

static void addCalpha(gemmi::Structure& st, const char* residueName, int seqId, float x, float y, float z) {
    gemmi::Residue residue;
    residue.name = residueName;
    residue.seqid = gemmi::SeqId(seqId, ' ');
    residue.label_seq = seqId;
    residue.het_flag = 'A';
    residue.entity_type = gemmi::EntityType::Polymer;

    gemmi::Atom atom;
    atom.name = "CA";
    atom.element = gemmi::El::C;
    atom.occ = 1.0f;
    atom.b_iso = 0.0f;
    atom.pos = gemmi::Position(roundToPdbPrecision(x), roundToPdbPrecision(y), roundToPdbPrecision(z));
    residue.atoms.push_back(atom);

    st.models.back().chains.back().residues.push_back(residue);
}

// modelTitles keeps the source entry per model, mmCIF has no per model title item
static bool writeStructureCif(gemmi::Structure& st,
                              const std::vector<std::pair<std::string, std::string>>& modelTitles,
                              const std::string& path) {
    try {
        gemmi::setup_entities(st);
        gemmi::cif::Document doc = gemmi::make_mmcif_document(st, calphaOutputGroups());
        if (modelTitles.empty() == false && doc.blocks.empty() == false) {
            gemmi::cif::Loop& loop = doc.blocks[0].init_loop("_foldseek_model.",
                                                             {"pdbx_PDB_model_num", "title"});
            for (size_t i = 0; i < modelTitles.size(); ++i) {
                loop.add_row({gemmi::cif::quote(modelTitles[i].first),
                              gemmi::cif::quote(modelTitles[i].second)});
            }
        }
        std::ofstream os(path.c_str());
        if (os.is_open() == false) {
            return false;
        }
        gemmi::cif::write_cif_to_stream(os, doc, gemmi::cif::Style::Pdbx);
        os.close();
        return os.fail() == false;
    } catch (const std::exception&) {
        return false;
    }
}

int convert2pdb(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    int outputMode = par.pdbOutputMode;
    const bool cifOutput = par.outputFormat == LocalParameters::STRUCTURE_OUTPUT_FORMAT_MMCIF;
    const char* extension = cifOutput ? ".cif" : ".pdb";
    // one file per chain, the PDB output keeps its existing grouping in every mode
    const bool splitPerChain = cifOutput && outputMode == LocalParameters::PDB_OUTPUT_MODE_SINGLECHAIN;
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
    if (outputMode == LocalParameters::PDB_OUTPUT_MODE_COMPLEX || LocalParameters::PDB_OUTPUT_MODE_SINGLECHAIN) {
        getKeyToIdMapIdToKeysMapIdVec(db, lookupFile, chainKeyToComplexIdMap, complexIdToChainKeysMap, complexIndices);
    }

    Debug(Debug::INFO) << "Start writing file to " << par.db2 << "\n";
    FILE* handle = NULL;
    gemmi::Structure multiModelSt;
    std::vector<std::pair<std::string, std::string>> multiModelTitles;
    if (outputMode == LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
        if (cifOutput) {
            // all entries share one data block, gemmi serializes it after the loop
            multiModelSt.name = Util::remove_extension(FileUtil::baseName(par.db2));
            multiModelSt.info["_entry.id"] = multiModelSt.name;
        } else {
            handle = fopen(par.db2.c_str(), "w");
            if (handle == NULL) {
                perror(par.db2.c_str());
                EXIT(EXIT_FAILURE);
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

    const char* threeLetterLookup[26] = { "ALA", "ASX", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "XLE", "LYS", "LEU", "MET", "ASN", "PYL", "PRO", "GLN", "ARG", "SER", "THR", "SEC", "VAL", "TRP", "XAA", "TYR", "GLX" };

#pragma omp parallel num_threads(localThreads)
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
        gemmi::Structure threadSt;
        std::vector<std::pair<std::string, std::string>> threadTitles;
        const bool sharedStructure = outputMode == LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL;
        gemmi::Structure& st = sharedStructure ? multiModelSt : threadSt;
        std::vector<std::pair<std::string, std::string>>& titles = sharedStructure ? multiModelTitles : threadTitles;
        std::string cifPath;
        Coordinate16 coords;

        KeyIterator* keyIterator;
        if (splitPerChain) {
            keyIterator = new DbKeyIterator(db);
        } else if (outputMode == LocalParameters::PDB_OUTPUT_MODE_COMPLEX || LocalParameters::PDB_OUTPUT_MODE_SINGLECHAIN) {
            keyIterator = new MapIterator(complexIdToChainKeysMap, complexIndices);
        } else {
            keyIterator = new DbKeyIterator(db);
        }
        const size_t size = keyIterator->getSize();
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < size; ++i) {
            std::pair<const unsigned int*, size_t> keys = keyIterator->getDbKeys(i);
            if (outputMode != LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
                unsigned int key = keys.first[0];
                // std::string name = db.getLookupEntryName(key);
                size_t lookupKey = db.getLookupIdByKey(key);
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
                cifPath = par.db2 + "/" + name + extension;

                unsigned int headerId = db_header.getId(key);
                const char* headerData = db_header.getData(headerId, thread_idx);
                const size_t headerLen = db_header.getEntryLen(headerId) - 2;
                if (cifOutput) {
                    st = gemmi::Structure();
                    titles.clear();
                    st.name = name;
                    st.info["_entry.id"] = name;
                    st.info["_struct.title"] = std::string(headerData, headerLen);
                    st.models.emplace_back("1");
                } else {
                    threadHandle = fopen(cifPath.c_str(), "w");
                    if (threadHandle == NULL) {
                        perror(cifPath.c_str());
                        EXIT(EXIT_FAILURE);
                    }
                    writeTitle(threadHandle, headerData, headerLen);
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
                        // mmCIF has no MODEL record, the model number carries the db key
                        std::string modelName = std::to_string(key);
                        st.models.emplace_back(modelName);
                        titles.emplace_back(modelName, std::string(headerData, headerLen));
                    } else {
                        fprintf(threadHandle, "MODEL % 8d\n", key);
                        writeTitle(threadHandle, headerData, headerLen);
                    }
                } else {
                    // std::string name = db.getLookupEntryName(key);
                    size_t lookupKey = db.getLookupIdByKey(key);
                    std::string name = db.getLookupEntryName(lookupKey);
                    chainName = name.substr(name.find_last_of('_') + 1);
                    if (chainName.size() == 0) {
                        chainName = "A";
                    }
                }
                if (cifOutput) {
                    // mmCIF keeps the full chain name, the PDB record has room for one character
                    st.models.back().chains.emplace_back(chainName);
                }
                for (size_t j = 0; j < seqLen; ++j) {
                    // make AA upper case
                    char aa = seqData[j] & ~0x20;
                    if (aa == '\0') {
                        aa = 'X';
                    }
                    const char* aa3 = threeLetterLookup[(int)(aa - 'A')];
                    if (cifOutput) {
                        addCalpha(st, aa3, (int)(j + 1), ca[j], ca[j + (1 * seqLen)], ca[j + (2 * seqLen)]);
                    } else {
                        fprintf(threadHandle, "ATOM  %5d  CA  %s %c%4d    %8.3f%8.3f%8.3f\n", (int)(j + 1), aa3, chainName[0], int(j + 1), ca[j], ca[j + (1 * seqLen)], ca[j + (2 * seqLen)]);
                    }
                }
                if (outputMode == LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL && cifOutput == false) {
                    fprintf(threadHandle, "ENDMDL\n");
                }
            }
            if (outputMode != LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
                if (cifOutput) {
                    if (writeStructureCif(st, titles, cifPath) == false) {
                        Debug(Debug::ERROR) << "Cannot write file " << cifPath << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                } else {
                    fclose(threadHandle);
                }
            }
        }
        delete keyIterator;
    }

    if (outputMode == LocalParameters::PDB_OUTPUT_MODE_MULTIMODEL) {
        if (cifOutput) {
            if (writeStructureCif(multiModelSt, multiModelTitles, par.db2) == false) {
                Debug(Debug::ERROR) << "Cannot write file " << par.db2 << "\n";
                EXIT(EXIT_FAILURE);
            }
        } else if (fclose(handle) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << par.db2 << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    db_ca.close();
    db_header.close();
    db.close();

    return EXIT_SUCCESS;
}
