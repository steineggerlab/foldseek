#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "interfacesearch.sh.h"

void setInterfaceSearchDefaults(LocalParameters *p) {
    p->removeTmpFiles = true;
    p->distanceThreshold = 10;
    p->minResidueNum = 4;
    //TODO: set paramters
    // p->filtInterfaceLddtThr = 0.0;
    // p->filtMultTmThr = 0.4;
    // p->filtChainTmThr = 0.0;
}

void mustseInterfaceSearch(LocalParameters *p) {
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
}

int interfacesearch(int argc, const char **argv, const Command &command) {

    LocalParameters &par = LocalParameters::getLocalInstance();

    setInterfaceSearchDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    mustseInterfaceSearch(&par);

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
    std::string input_target = par.filenames.back().c_str();
    cmd.addVariable("INPUT_TARGET", input_target.c_str());
    par.filenames.pop_back();
    std::string input_query = par.filenames.back().c_str();
    cmd.addVariable("INPUT_QUERY", input_query.c_str());
    par.filenames.pop_back();

    int dbtype_target = FileUtil::parseDbType((input_target+"_ss").c_str());
    int isInterfacedb_target = (DBReader<unsigned int>::getExtendedDbtype(dbtype_target) & LocalParameters::DBTYPE_EXTENDED_INTERFACE);
    int padded_target = (DBReader<unsigned int>::getExtendedDbtype(dbtype_target) & Parameters::DBTYPE_EXTENDED_GPU);
    cmd.addVariable("ISINTERFACEDB_TARGET", isInterfacedb_target ? "TRUE" : NULL);
    cmd.addVariable("NOTPADDED_TARGET", padded_target ? NULL : "TRUE");

    int dbtype_query = FileUtil::parseDbType((input_query+"_ss").c_str());
    int isInterfacedb_query = (DBReader<unsigned int>::getExtendedDbtype(dbtype_query) & LocalParameters::DBTYPE_EXTENDED_INTERFACE);
    cmd.addVariable("ISINTERFACEDB_QUERY", isInterfacedb_query ? "TRUE" : NULL);
    int padded_query = (DBReader<unsigned int>::getExtendedDbtype(dbtype_query) & Parameters::DBTYPE_EXTENDED_GPU);
    cmd.addVariable("NOTPADDED_QUERY", padded_query ? NULL : "TRUE");

    cmd.addVariable("GPU", par.gpu ? "TRUE" : NULL);
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("MULTIMERSEARCH_PAR", par.createParameterString(par.multimersearchworkflow, true).c_str()); 

    std::string program = tmpDir + "/interfacesearch.sh";
    FileUtil::writeFile(program, interfacesearch_sh, interfacesearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}