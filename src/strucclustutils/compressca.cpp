#include "LocalParameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Coordinate16.h"

#ifdef OPENMP
#include <omp.h>
#endif

int compressca(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.overrideParameterDescription(par.PARAM_COORD_STORE_MODE, "Coordinate storage mode: \n1: C-alpha as float\n2: C-alpha as difference (uint16_t)\n3: Plain text list of floats", "^[1-3]{1}$", 0);
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string caDbData = par.db1 + "_ca";
    std::string caDbIndex = par.db1 + "_ca.index";
    DBReader<unsigned int> caDb(caDbData.c_str(), caDbIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    caDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    const bool isPlainTextCA = caDb.getDbtype() == LocalParameters::DBTYPE_GENERIC_DB;

    DBReader<unsigned int> seqDb(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    seqDb.open(DBReader<unsigned int>::NOSORT);

    int dbtype = LocalParameters::DBTYPE_CA_ALPHA;
    if (par.coordStoreMode == LocalParameters::COORD_STORE_MODE_CA_PLAIN_TEXT) {
        dbtype = LocalParameters::DBTYPE_GENERIC_DB;
    }
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, dbtype);
    writer.open();

    Debug::Progress progress(caDb.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Coordinate16 coords;
        std::vector<int8_t> camol;

        std::string plain;
        if (par.coordStoreMode == LocalParameters::COORD_STORE_MODE_CA_PLAIN_TEXT) {
            plain.reserve(10 * 1024);
        }
        std::vector<float> floatCoords;

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < caDb.getSize(); i++) {
            progress.updateProgress();
            unsigned int key = caDb.getDbKey(i);
            char *data = caDb.getData(i, thread_idx);
            size_t length = caDb.getEntryLen(i) -1;
            unsigned int seqId = seqDb.getId(key);
            size_t chainLen = seqDb.getSeqLen(seqId);

            if (isPlainTextCA) {
                floatCoords.clear();
                char *next;
                char *current = data;
                while (current < data + length) {
                    floatCoords.push_back(strtof(current, &next));
                    if (next == current) {
                        break;
                    }
                    current = next + 1;
                }

                data = reinterpret_cast<char*>(floatCoords.data());
                length = floatCoords.size() * sizeof(float);
            }

            const size_t uncompressedSize = chainLen * (3 * sizeof(float));
            if (par.coordStoreMode == LocalParameters::COORD_STORE_MODE_CA_DIFF) {
                if (length >= uncompressedSize) {
                    // coords in COORD_STORE_MODE_CA_FLOAT, so convert to diff
                    camol.resize((chainLen - 1) * 3 * sizeof(int16_t) + 3 * sizeof(float));
                    int16_t* camolf16 = reinterpret_cast<int16_t*>(camol.data());
                    // check if any of the coordinates is too large to be stored as int16_t
                    if (!Coordinate16::convertToDiff16(chainLen, (float*)(data) + 0 * chainLen, camolf16, 1)
                    && !Coordinate16::convertToDiff16(chainLen, (float*)(data) + 1 * chainLen, camolf16 + 1 * (chainLen + 1), 1)
                    && !Coordinate16::convertToDiff16(chainLen, (float*)(data) + 2 * chainLen, camolf16 + 2 * (chainLen + 1), 1)) {
                        writer.writeData((const char*)camol.data(), (chainLen - 1) * 3 * sizeof(uint16_t) + 3 * sizeof(float) + 1 * sizeof(uint8_t), key, thread_idx);
                        continue;
                    }
                }
                // if we didn't get into the if above, the coordinates were already in diff mode or conversion failed
                writer.writeData(data, length, key, thread_idx);
            } else if (par.coordStoreMode == LocalParameters::COORD_STORE_MODE_CA_FLOAT) {
                if (length >= uncompressedSize) {
                    // already in COORD_STORE_MODE_CA_FLOAT, so directly write the data
                    writer.writeData(data, length, key, thread_idx);
                } else {
                    // need to convert the coordinates to COORD_STORE_MODE_CA_FLOAT
                    float* uncompressed = coords.read(data, chainLen, length);
                    writer.writeData((const char*)uncompressed, uncompressedSize, key, thread_idx);
                }
            } else if (par.coordStoreMode == LocalParameters::COORD_STORE_MODE_CA_PLAIN_TEXT) {
                float* uncompressed = coords.read(data, chainLen, length);
                plain.append(SSTR(uncompressed[0]));
                for (size_t i = 1; i < (chainLen * 3); i++) {
                    plain.append(",");
                    plain.append(SSTR(uncompressed[i]));
                }
                plain.append("\n");
                writer.writeData(plain.c_str(), plain.size(), key, thread_idx);
                plain.clear();
            } else {
                Debug(Debug::ERROR) << "Unknown storage mode: " << par.coordStoreMode << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }

    writer.close(true);
    seqDb.close();
    caDb.close();

    return EXIT_SUCCESS;
}
