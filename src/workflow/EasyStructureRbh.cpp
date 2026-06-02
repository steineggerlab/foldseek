#include <cassert>
#include "LinsearchIndexReader.h"
#include "PrefilteringIndexReader.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "LocalParameters.h"
namespace structureRbh{
#include "easyrbh.sh.h"
}
#include "easystructurerbh.sh.h"
#include "structty.h"


int structureeasyrbh(int argc, const char **argv, const Command &command) {
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
    for (size_t i = 0; i < par.result2structprofile.size(); i++){
        par.result2structprofile[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.removeTmpFiles = true;
    par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    par.writeLookup = false;
    par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_SOFT;
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    par.PARAM_REMOVE_TMP_FILES.wasSet = true;
    par.PARAM_ALIGNMENT_MODE.wasSet = true;

    bool needBacktrace = false;
    bool needTaxonomy = false;
    bool needTaxonomyMapping = false;
    bool needLookup = false;

    {
        bool needSequenceDB = false;
        bool need3DiDB = false;
        bool needFullHeaders = false;
        bool needSource = false;
        bool needQCA = false;
        bool needTCA = false;
        bool needTMalign = false;
        bool needLDDT = false;
        LocalParameters::getOutputFormat(
            par.formatAlignmentMode, par.outfmt, needSequenceDB, need3DiDB, needBacktrace, needFullHeaders,
            needLookup, needSource, needTaxonomyMapping, needTaxonomy, needQCA, needTCA, needTMalign, needLDDT
        );
    }

    if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM || par.greedyBestHits) {
        needBacktrace = true;
    }
    if (needBacktrace) {
        Debug(Debug::INFO) << "Alignment backtraces will be computed, since they were requested by output format.\n";
        par.addBacktrace = true;
        par.PARAM_ADD_BACKTRACE.wasSet = true;
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
    std::string resultsPath = par.filenames.back();
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

    cmd.addVariable("QUERY", par.filenames.back().c_str());
    {
        std::vector<MMseqsParameter*> searchParams = par.removeParameter(par.structuresearchworkflow, par.PARAM_VIEW_RESULTS);
        cmd.addVariable("SEARCH_PAR", par.createParameterString(searchParams, true).c_str());
    }
    // When viewing results, defer tmp cleanup so StrucTTY can read the tmp DBs (D10).
    // Actual cleanup is re-invoked after launch in Step 3.
    cmd.addVariable("REMOVE_TMP", (par.removeTmpFiles && !par.viewResults) ? "TRUE" : NULL);
    cmd.addVariable("LEAVE_INPUT", par.dbOut ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    cmd.addVariable("CREATEDB_QUERY_PAR", par.createParameterString(par.structurecreatedb).c_str());
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.structurecreatedb).c_str());
    cmd.addVariable("CONVERT_PAR", par.createParameterString(par.convertalignments).c_str());

    // Write helper scripts to tmpDir
    std::string easyrbhProgram = tmpDir + "/easyrbh.sh";
    FileUtil::writeFile(easyrbhProgram, structureRbh::easyrbh_sh, structureRbh::easyrbh_sh_len);

    std::string program = tmpDir + "/easystructurerbh.sh";
    FileUtil::writeFile(program, easystructurerbh_sh, easystructurerbh_sh_len);
    std::string argString = program;
    for (const auto& s : par.filenames) { argString += " "; argString += s; }
    if (std::system(argString.c_str()) != EXIT_SUCCESS) { EXIT(EXIT_FAILURE); }
    if (par.viewResults) {
        structty::RunOptions opts;
        const std::string queryInput = par.filenames[0];
        const bool queryIsDb  = FileUtil::fileExists((queryInput + ".dbtype").c_str());
        const bool targetIsDb = FileUtil::fileExists((target + ".dbtype").c_str());
        // Fallback plaintext query load (Step 4 will switch to query DB read).
        opts.input_files.push_back(queryInput);
        opts.foldseek_query_db = queryIsDb  ? queryInput : (tmpDir + "/query");
        opts.foldseek_db       = targetIsDb ? target     : (tmpDir + "/target");
        opts.foldseek_file     = resultsPath;
        structty::run(opts);
        // D10: cleanup was deferred past launch (Step 2); re-invoke the workflow
        // script in cleanup-only mode now that the viewer has closed.
        if (par.removeTmpFiles) {
            cmd.addVariable("CLEANUP_ONLY", "TRUE");
            if (std::system(argString.c_str()) != EXIT_SUCCESS) { EXIT(EXIT_FAILURE); }
        }
    }
    return EXIT_SUCCESS;
}

