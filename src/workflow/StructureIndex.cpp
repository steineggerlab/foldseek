#include "LocalParameters.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"

#include "structureindex.sh.h"

extern void setStructureSearchWorkflowDefaults(LocalParameters *p);
extern void setStructureSearchMustPassAlong(LocalParameters *p);

int structureindex(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();

    setStructureSearchWorkflowDefaults(&par);
    par.PARAM_COV_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_C.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MIN_SEQ_ID.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_SEQS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SPLIT.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.splitsequence.size(); i++) {
        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.extractorfs.size(); i++) {
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.splitsequence.size(); i++) {
        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++) {
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);
    setStructureSearchMustPassAlong(&par);

    std::string tmpDir = par.db2;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.createindex));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    std::vector<MMseqsParameter*> createIndexWithoutIndexSubset;
    for (size_t i = 0; i < par.createindex.size(); i++) {
        if (par.createindex[i]->uniqid != par.PARAM_INDEX_SUBSET.uniqid) {
            createIndexWithoutIndexSubset.push_back(par.createindex[i]);
        }
    }

    const bool excludeKmers = (par.indexExclude & LocalParameters::INDEX_EXCLUDE_KMER_INDEX) != 0;
    const bool excludeCa = (par.indexExclude & LocalParameters::INDEX_EXCLUDE_CA) != 0;

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("CREATEINDEX_PAR", par.createParameterString(createIndexWithoutIndexSubset, true).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("INDEX_DB_CA_KEY", SSTR(LocalParameters::INDEX_DB_CA_KEY).c_str());
    cmd.addVariable("SS_SUBSET_MODE", SSTR(excludeKmers ? 3 : 1).c_str());
    cmd.addVariable("INCLUDE_CA", excludeCa == false ? "TRUE" : NULL);

    std::string program(tmpDir + "/structureindex.sh");
    FileUtil::writeFile(program, structureindex_sh, structureindex_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // should never get here
    return EXIT_FAILURE;
}
