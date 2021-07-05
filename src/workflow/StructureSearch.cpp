#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "PrefilteringIndexReader.h"
#include "strucclust.sh.h"

void setStructureSearchWorkflowDefaults(LocalParameters *p) {
    p->sensitivity = 5.7;
    p->removeTmpFiles = true;
    p->writeLookup = false;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
}


int structuresearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.createdb.size(); i++){
        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.splitsequence.size(); i++) {
        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.result2profile.size(); i++){
        par.result2profile[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setStructureSearchWorkflowDefaults(&par);

    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);


    bool needBacktrace = false;
    bool needTaxonomy = false;
    bool needTaxonomyMapping = false;
    bool needLookup = false;

    {
        bool needSequenceDB = false;
        bool needFullHeaders = false;
        bool needSource = false;
        Parameters::getOutputFormat(par.formatAlignmentMode, par.outfmt, needSequenceDB, needBacktrace, needFullHeaders,
                                    needLookup, needSource, needTaxonomyMapping, needTaxonomy);
    }

    if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM || par.greedyBestHits) {
        needBacktrace = true;
    }
    if (needBacktrace) {
        Debug(Debug::INFO) << "Alignment backtraces will be computed, since they were requested by output format.\n";
        par.addBacktrace = true;
        par.PARAM_ADD_BACKTRACE.wasSet = true;
    }
    if(needLookup){
        par.writeLookup = true;
    }

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULTS", par.filenames.back().c_str());
    par.filenames.pop_back();
    std::string target = par.filenames.back().c_str();
    cmd.addVariable("TARGET", target.c_str());
    par.filenames.pop_back();

    if (needTaxonomy || needTaxonomyMapping) {
        std::vector<std::string> missingFiles = Parameters::findMissingTaxDbFiles(target);
        if (missingFiles.empty() == false) {
            Parameters::printTaxDbError(target, missingFiles);
            EXIT(EXIT_FAILURE);
        }
    }


    const bool isIndex = PrefilteringIndexReader::searchForIndex(target).empty() == false;
    cmd.addVariable("INDEXEXT", isIndex ? ".idx" : NULL);
    cmd.addVariable("SEARCH_MODULE", "search");
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow, true).c_str());


    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("GREEDY_BEST_HITS", par.greedyBestHits ? "TRUE" : NULL);
    cmd.addVariable("LEAVE_INPUT", par.dbOut ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    cmd.addVariable("CREATEDB_QUERY_PAR", par.createParameterString(par.createdb).c_str());
    par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("CONVERT_PAR", par.createParameterString(par.convertalignments).c_str());
    cmd.addVariable("SUMMARIZE_PAR", par.createParameterString(par.summarizeresult).c_str());

    std::string program = tmpDir + "/structuresearch.sh";
    FileUtil::writeFile(program, easystructuresearch_sh, easystructuresearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
