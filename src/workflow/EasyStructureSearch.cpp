#include <dirent.h>
#include <fstream>
#include <cassert>
#include <LocalParameters.h>
#include "LinsearchIndexReader.h"
#include "PrefilteringIndexReader.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "DBReader.h"
#include "Parameters.h"
#include "easystructuresearch.sh.h"
#include "structty.h"

// 17-column layout that StrucTTY's FoldseekParser parses as fmt 17
// (see lib/structty/src/structure/FoldseekParser.hpp). Column order matters:
// a different order makes lddt/qtmscore swap places and qaln/taln be read as numbers.
static const char VIEWER_OUTFMT[] =
    "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
    "evalue,bits,lddt,qtmscore,ttmscore,qaln,taln";

// StrucTTY color modes accepted by structty::RunOptions::mode
// (lib/structty/src/structure/Parameters.cpp). mmseqs does not regex-check
// string parameters, so validate here instead of relying on the regex.
static void validateStructtyMode(const std::string &mode) {
    if (mode == "protein" || mode == "chain" || mode == "rainbow" || mode == "plddt"
        || mode == "interface" || mode == "conservation" || mode == "aligned") {
        return;
    }
    Debug(Debug::ERROR) << "Invalid --structty-mode: " << mode << "\n"
                        << "Choose one of: protein, chain, rainbow, plddt, interface, "
                        << "conservation, aligned\n";
    EXIT(EXIT_FAILURE);
}

// --view-structty needs 3D coordinates, so reject inputs that cannot supply any
// before createdb runs -- otherwise the user waits out a whole search first.
// Keep this identical in EasyStructureSearch, StructureSearch, EasyMultimerSearch
// and MultimerSearch (the four workflows that support the viewer).
static bool isSequenceFasta(const std::string &path) {
    std::ifstream ifs(path.c_str());
    if (ifs.is_open() == false) {
        return false;
    }
    std::string line;
    while (std::getline(ifs, line)) {
        const size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            continue;
        }
        return line[first] == '>';
    }
    return false;
}

// A directory is judged by one representative file; scanning every entry would be
// wasteful on large inputs and one FASTA is enough to know the input is sequences.
static std::string firstFileInDirectory(const std::string &dir) {
    DIR *handle = opendir(dir.c_str());
    if (handle == NULL) {
        return "";
    }
    std::string found;
    while (dirent *entry = readdir(handle)) {
        const std::string name(entry->d_name);
        if (name == "." || name == "..") {
            continue;
        }
        const std::string path = dir + "/" + name;
        if (FileUtil::fileExists(path.c_str()) && FileUtil::directoryExists(path.c_str()) == false) {
            found = path;
            break;
        }
    }
    closedir(handle);
    return found;
}

static void validateStructtyInputs(const std::vector<std::string> &inputs) {
    for (size_t i = 0; i < inputs.size(); i++) {
        const std::string &input = inputs[i];
        if (FileUtil::fileExists((input + ".dbtype").c_str())) {
            if (FileUtil::fileExists((input + "_ca.dbtype").c_str()) == false) {
                Debug(Debug::ERROR) << "--view-structty needs C-alpha coordinates, but the database "
                                    << input << " has none (" << input << "_ca is missing).\n"
                                    << "It was built from sequences, or with --index-exclude 2. "
                                    << "Rebuild it from PDB/mmCIF structures.\n";
                EXIT(EXIT_FAILURE);
            }
            continue;
        }
        const std::string probe = FileUtil::directoryExists(input.c_str())
                                  ? firstFileInDirectory(input) : input;
        if (probe.empty() == false && isSequenceFasta(probe)) {
            Debug(Debug::ERROR) << "--view-structty cannot render " << input
                                << ": it is a sequence FASTA and carries no 3D coordinates.\n"
                                << "Pass PDB/mmCIF structures instead. createdb --prostt5-model "
                                << "predicts 3Di from sequence, but writes no C-alpha coordinates.\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

void setEasyStructureSearchDefaults(Parameters *p) {
    // TODO: 7-mer sensitivity is not optimized yet
    p->kmerSize = 6;
    p->sensitivity = 9.5;
    p->maxResListLen = 1000;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->removeTmpFiles = true;
    p->reportMode = 2;
}
void setEasyStructureSearchMustPassAlong(Parameters *p) {
    p->PARAM_K.wasSet = true;
    p->PARAM_S.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
}

int easystructuresearch(int argc, const char **argv, const Command &command) {
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

    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setEasyStructureSearchDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasyStructureSearchMustPassAlong(&par);
    validateStructtyMode(par.structtyMode);
    if (par.viewResults && par.filenames.size() > 2) {
        // <query...> <target> <results> <tmpDir>
        validateStructtyInputs(std::vector<std::string>(par.filenames.begin(),
                                                        par.filenames.end() - 2));
    }
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
        if (par.reportMode != 2) {
            needTaxonomy = true;
            needTaxonomyMapping = true;
        }
    }

    if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM ||
        par.formatAlignmentMode == LocalParameters::FORMAT_ALIGNMENT_PDB_SUPERPOSED  ||
        par.greedyBestHits) {
        needBacktrace = true;
    }
    if (needBacktrace) {
        Debug(Debug::INFO) << "Alignment backtraces will be computed, since they were requested by output format.\n";
        par.addBacktrace = true;
        par.PARAM_ADD_BACKTRACE.wasSet = true;
    }
    // The StrucTTY viewer m8 carries qaln/taln/lddt, which all need alignment
    // backtraces. needBacktrace above is derived from the *user's* --format-output
    // and knows nothing about the viewer format, so force backtraces on here.
    if (par.viewResults && par.addBacktrace == false) {
        Debug(Debug::INFO) << "Alignment backtraces will be computed for the StrucTTY viewer.\n";
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
    cmd.addVariable("LEAVE_INPUT", par.dbOut ? "TRUE" : NULL);
    par.filenames.pop_back();

    if (needTaxonomy || needTaxonomyMapping) {
        std::vector<std::string> missingFiles = Parameters::findMissingTaxDbFiles(target);
        if (missingFiles.empty() == false) {
            Parameters::printTaxDbError(target, missingFiles);
            EXIT(EXIT_FAILURE);
        }
    }

    int dbtype = FileUtil::parseDbType((target+"_ss").c_str());
    int padded = (DBReader<unsigned int>::getExtendedDbtype(dbtype) & Parameters::DBTYPE_EXTENDED_GPU);
    cmd.addVariable("NOTPADDED", padded ? NULL : "TRUE");

    const bool isIndex = PrefilteringIndexReader::searchForIndex(target).empty() == false;
    cmd.addVariable("INDEXEXT", isIndex ? ".idx" : NULL);

    cmd.addVariable("CREATELININDEX_PAR", NULL);
    {
        std::vector<MMseqsParameter*> searchParams = par.removeParameter(par.structuresearchworkflow, par.PARAM_VIEW_RESULTS);
        searchParams = par.removeParameter(searchParams, par.PARAM_STRUCTTY_MODE);
        searchParams = par.removeParameter(searchParams, par.PARAM_STRUCTTY_SS);
        cmd.addVariable("SEARCH_PAR", par.createParameterString(searchParams, true).c_str());
    }
    cmd.addVariable("LNDB_PAR", par.createParameterString(par.verbandcompression, true).c_str());

    // When viewing results, defer tmp cleanup so StrucTTY can read the tmp DBs (D10).
    // Actual cleanup is re-invoked after launch in Step 3.
    cmd.addVariable("REMOVE_TMP", (par.removeTmpFiles && !par.viewResults) ? "TRUE" : NULL);
    cmd.addVariable("GREEDY_BEST_HITS", par.greedyBestHits ? "TRUE" : NULL);
    cmd.addVariable("GPU", par.gpu ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("MAKEPADDEDSEQDB_PAR", par.createParameterString(par.makepaddeddb).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    cmd.addVariable("CREATEDB_QUERY_PAR", par.createParameterString(par.structurecreatedb).c_str());
    par.prostt5Model = "";
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.structurecreatedb).c_str());
    cmd.addVariable("CONVERT_PAR", par.createParameterString(par.convertalignments).c_str());
    // Viewer m8: same convertalis parameters, but with the fixed viewer layout.
    // Swap par.outfmt instead of appending --format-output, because convertalis
    // rejects a duplicated flag ("Duplicate parameter --format-output").
    std::string viewerConvertPar;
    {
        const std::string userOutfmt = par.outfmt;
        const bool userOutfmtWasSet = par.PARAM_FORMAT_OUTPUT.wasSet;
        par.outfmt = VIEWER_OUTFMT;
        par.PARAM_FORMAT_OUTPUT.wasSet = true;
        viewerConvertPar = par.createParameterString(par.convertalignments);
        par.outfmt = userOutfmt;
        par.PARAM_FORMAT_OUTPUT.wasSet = userOutfmtWasSet;
    }
    cmd.addVariable("VIEW_RESULTS", par.viewResults ? "TRUE" : NULL);
    cmd.addVariable("VIEWER_CONVERT_PAR", par.viewResults ? viewerConvertPar.c_str() : NULL);
    cmd.addVariable("SUMMARIZE_PAR", par.createParameterString(par.summarizeresult).c_str());
    
    cmd.addVariable("TAXONOMY", needTaxonomy && needTaxonomyMapping && par.reportMode != 2 ? "TRUE" : NULL);
    cmd.addVariable("TAXONOMYREPORT_PAR", par.createParameterString(par.taxonomyreport).c_str());
    std::string program = tmpDir + "/easystructuresearch.sh";
    FileUtil::writeFile(program, easystructuresearch_sh, easystructuresearch_sh_len);
    std::string argString = program;
    for (const auto& s : par.filenames) { argString += " "; argString += s; }
    if (std::system(argString.c_str()) != EXIT_SUCCESS) { EXIT(EXIT_FAILURE); }
    if (par.viewResults) {
        structty::RunOptions opts;
        const std::string queryInput = par.filenames[0];
        const bool queryIsDb  = FileUtil::fileExists((queryInput + ".dbtype").c_str());
        const bool targetIsDb = FileUtil::fileExists((target + ".dbtype").c_str());
        // The viewer probes each path and picks its scene from what it finds, so the
        // query DB goes in as the positional query and the target DB as -fst.
        opts.input_files.push_back(queryIsDb ? queryInput : (tmpDir + "/query"));
        opts.foldseek_target   = targetIsDb ? target : (tmpDir + "/target");
        opts.foldseek_result   = tmpDir + "/viewer_results.m8";
        opts.mode              = par.structtyMode;
        opts.show_structure    = par.structtyShowStructure;
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


