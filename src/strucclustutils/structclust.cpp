#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "structclust.sh.h"
#include "DBReader.h"

int structclust(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    CommandCaller cmd;

    DBReader<unsigned int> *reader;
    reader = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);

    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(reader->getDbtype());
    if (extended & Parameters::DBTYPE_EXTENDED_SET) {
        cmd.addVariable("NEEDSET", "1");
    } else {
        cmd.addVariable("NEEDSET", "0");
    }


    cmd.addVariable("RESULT",par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("ALN", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("INPUT", par.filenames.back().c_str());
    par.filenames.pop_back();

    cmd.addVariable("CLUST_PAR", par.createParameterString(par.clust).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = par.db2 + ".sh";
    FileUtil::writeFile(program, structclust_sh, structclust_sh_len);
    cmd.execProgram(FileUtil::getRealPathFromSymLink(program).c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
