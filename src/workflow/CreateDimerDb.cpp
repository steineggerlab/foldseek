#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "createdimerdb.sh.h"

void setCreateDimerDefaults(Parameters *p) {
    p->removeTmpFiles = true;
}
void setCreateDimerMustPassAlong(Parameters *p) {
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
}

int createdimerdb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setCreateDimerDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setCreateDimerMustPassAlong(&par);

    CommandCaller cmd;
    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("OUT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("IN", par.filenames.back().c_str());
    par.filenames.pop_back();

    cmd.addVariable("FILTERDIMERDB_PAR", par.createParameterString(par.filterdimerdb).c_str());
    cmd.addVariable("CREATEDIMERDB_PAR", par.createParameterString(par.createdimerdb).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

    std::string program = tmpDir + "/createdimerdb.sh";
    FileUtil::writeFile(program, createdimerdb_sh, createdimerdb_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}