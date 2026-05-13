#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "interfacecluster.sh.h"

void setInterfaceClusterDefaults(LocalParameters *p) {
    p->removeTmpFiles = true;
    p->distanceThreshold = 10;
    p->minResidueNum = 4;
    p->filtInterfaceLddtThr = 0.0;
    p->filtMultTmThr = 0.4;
    p->filtChainTmThr = 0.0;
}

void mustseInterfaceCluster(LocalParameters *p) {
    p->clusteringSetMode = 1;
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
}

int interfacecluster(int argc, const char **argv, const Command &command) {

    LocalParameters &par = LocalParameters::getLocalInstance();

    setInterfaceClusterDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    mustseInterfaceCluster(&par);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    CommandCaller cmd;

    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("OUT", par.filenames.back().c_str());
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
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("MULTIMERCLUSTER_PAR", par.createParameterString(par.multimerclusterworkflow, true).c_str()); 

    std::string program = tmpDir + "/interfacecluster.sh";
    FileUtil::writeFile(program, interfacecluster_sh, interfacecluster_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}