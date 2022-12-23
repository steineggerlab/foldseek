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
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string caDbData = par.db1 + "_ca";
    std::string caDbIndex = par.db1 + "_ca.index";
    DBReader<unsigned int> caDb(caDbData.c_str(), caDbIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    caDb.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> seqDb(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    seqDb.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
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

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < caDb.getSize(); i++) {
            progress.updateProgress();
            unsigned int key = caDb.getDbKey(i);
            char *data = caDb.getData(i, thread_idx);
            size_t length = caDb.getEntryLen(i) -1;
            size_t chainLen = seqDb.getSeqLen(key);

            if (length >= (chainLen * (3 * sizeof(float)))) {
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
            writer.writeData(data, length, key, thread_idx);
        }
    }

    writer.close(true);
    seqDb.close();
    caDb.close();

    return EXIT_SUCCESS;
}
