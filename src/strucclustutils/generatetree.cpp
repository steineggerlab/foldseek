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

    for (size_t i = dbFrom; i < dbFrom + dbSize; ++i) {
        char *data = reader.getData(i, thread_idx);

        bool readFirst = false;
//        writer.writeStart(thread_idx);
        while (*data != '\0') {
            Matcher::result_t domain = Matcher::parseAlignmentRecord(data, true);
            data = Util::skipLine(data);
            std::cout << domain.dbKey << ": " << domain.score << std::endl;
            }
        }
//    char *data = alndbr.getData(i, thread_idx);
//    Matcher::readAlignmentResults(results, data, false);

//    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
//    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    std::cout << "Hello, worldyMcWorldface!" << std::endl;
    return EXIT_SUCCESS;
}
