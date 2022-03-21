#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Alignment.h"
#include "structureto3diseqdist.h"
#include "StructureSmithWaterman.h"
#include "StructureUtil.h"
#include <tuple>

#ifdef OPENMP
#include <omp.h>
#endif

#define	EXIT_FAILURE	1
#define	EXIT_SUCCESS	0


int generatetree(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    std::vector<Matcher::result_t> results;
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    size_t dbFrom = 0;
    size_t dbSize = reader.getSize();
    unsigned int thread_idx = 0;
    std::vector<Matcher::result_t> result;

    int score_list[dbSize][dbSize];
    int queryKey  = 0;

    for(int n = 0; n < dbSize; n++){
        for(int m = 0; m < dbSize; m++){
            score_list[n][m] = 0;
        }
    }

    for (size_t i = dbFrom; i < dbFrom + dbSize; ++i) {
        char *data = reader.getData(i, thread_idx);
        Matcher::readAlignmentResults(result, data, 0);
        queryKey = reader.getDbKey(i);
        for (int j = 0; j < result.size(); j++) {
            score_list[queryKey][result[j].dbKey] = result[j].score;
        }
        result.clear();
    }

    for(int n = 0; n < dbSize; n++){
        for(int m = 0; m < dbSize; m++){
            std::cout << score_list[n][m] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Hello, worldyMcWorldface!" << std::endl;
    return EXIT_SUCCESS;
}
