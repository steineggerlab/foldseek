#include <cassert>

#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "DBReader.h"
#include "Debug.h"

#include "easyinterfacecluster.sh.h"

void setEasyInterfaceClusterDefaults(Parameters *p) {
}

void setEasyInterfaceClusterMustPassAlong(Parameters *p) {
}

int easyinterfacecluster(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();

    setEasyInterfaceClusterDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasyInterfaceClusterMustPassAlong(&par);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }

    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULT", par.filenames.back().c_str());
    par.filenames.pop_back();
    std::string input = par.filenames.back().c_str();
    cmd.addVariable("INPUT", input.c_str());
    par.filenames.pop_back();

    int dbtype = FileUtil::parseDbType((input+"_ss").c_str());
    int isInterfacedb = (DBReader<unsigned int>::getExtendedDbtype(dbtype) & LocalParameters::DBTYPE_EXTENDED_INTERFACE);
    int padded = (DBReader<unsigned int>::getExtendedDbtype(dbtype) & Parameters::DBTYPE_EXTENDED_GPU);
    cmd.addVariable("ISINTERFACEDB", isInterfacedb ? "TRUE" : NULL);
    cmd.addVariable("NOTPADDED", padded ? NULL : "TRUE");

    cmd.addVariable("GPU", par.gpu ? "TRUE" : NULL);
    cmd.addVariable("MAKEPADDEDSEQDB_PAR", par.createParameterString(par.makepaddeddb).c_str());
    cmd.addVariable("EASYMULTIMERCLUSTER_PAR", par.createParameterString(par.easymultimerclusterworkflow,true).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

    std::string program = tmpDir + "/easyinterfacecluster.sh";
    FileUtil::writeFile(program, easyinterfacecluster_sh, easyinterfacecluster_sh_len);

    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
