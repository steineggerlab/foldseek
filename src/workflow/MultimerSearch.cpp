#include <dirent.h>
#include <fstream>
#include <cassert>

#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#include "multimersearch.sh.h"
#include "structty.h"

// StrucTTY color modes accepted by structty::RunOptions::mode
// (lib/structty/src/structure/Parameters.cpp). mmseqs does not regex-check
// string parameters, so validate here instead of relying on the regex.
// Keep this identical in the other workflows that support the viewer.
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

int multimersearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, false, Parameters::PARSE_VARIADIC, 0);
    validateStructtyMode(par.structtyMode);
    if (par.viewResults) {
        if (par.multimerReportMode == 0) {
            Debug(Debug::ERROR) << "--view-structty needs the multimer report, but "
                                << "--multimer-report-mode is 0.\n"
                                << "Drop --view-structty, or leave the report mode at 1.\n";
            EXIT(EXIT_FAILURE);
        }
        if (par.filenames.size() > 2) {
            // <queryDB> <targetDB> <alignmentDB> <tmpDir>
            validateStructtyInputs(std::vector<std::string>(par.filenames.begin(),
                                                            par.filenames.begin() + 2));
        }
    }
    if(par.PARAM_FORMAT_OUTPUT.wasSet == false){
        par.outfmt = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid";
    }
    par.addBacktrace = true;
    par.PARAM_ADD_BACKTRACE.wasSet = true;
    par.printParameters(command.cmd, argc, argv, *command.params);

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

    if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM ||
        par.formatAlignmentMode == LocalParameters::FORMAT_ALIGNMENT_PDB_SUPERPOSED  ||
        par.greedyBestHits) {
        needBacktrace = true;
    }
    if (needLookup) {
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
    double eval =  par.evalThr;
    if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_TMALIGN){
        par.evalThr = par.eValueThrExpandMultimer;
        cmd.addVariable("MULTIMER_ALIGNMENT_ALGO", "tmalign");
//        cmd.addVariable("MULTIMER_ALIGN_PREF_PAR", par.createParameterString(par.structurealign).c_str());
        cmd.addVariable("MULTIMER_ALIGN_PAR", par.createParameterString(par.tmalign).c_str());
    }else if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA || par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI){
        par.evalThr = par.eValueThrExpandMultimer;
        cmd.addVariable("MULTIMER_ALIGNMENT_ALGO", "structurealign");
        cmd.addVariable("MULTIMER_ALIGN_PAR", par.createParameterString(par.structurealign).c_str());
    }
    par.evalThr = eval;
    switch(par.prefMode){
        case LocalParameters::PREF_MODE_KMER:
            cmd.addVariable("PREFMODE", "KMER");
            break;
        case LocalParameters::PREF_MODE_UNGAPPED:
            cmd.addVariable("PREFMODE", "UNGAPPED");
            break;
        case LocalParameters::PREF_MODE_EXHAUSTIVE:
            cmd.addVariable("PREFMODE", "EXHAUSTIVE");
            break;
    }
    if(par.exhaustiveSearch){
        cmd.addVariable("PREFMODE", "EXHAUSTIVE");
    }
    cmd.addVariable("NO_REPORT", par.multimerReportMode == 0 ? "TRUE" : NULL);
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("OUTPUT", par.filenames.back().c_str());
    par.filenames.pop_back();
    const std::string targetDb = par.filenames.back();
    cmd.addVariable("TARGETDB", targetDb.c_str());
    par.filenames.pop_back();
    const std::string queryDb = par.filenames.back();
    cmd.addVariable("QUERYDB", queryDb.c_str());
    cmd.addVariable("LEAVE_INPUT", par.dbOut ? "TRUE" : NULL);
    par.filenames.pop_back();

    // initial search speed up!
    par.addBacktrace = par.exhaustiveSearch;
    par.alignmentType = par.exhaustiveSearch ? par.alignmentType : LocalParameters::ALIGNMENT_TYPE_3DI_AA;

    {
        // The inner `search` call must not inherit the viewer flags, or it launches a
        // second StrucTTY of its own and waits for input before this workflow finishes.
        std::vector<MMseqsParameter*> searchParams = par.removeParameter(par.structuresearchworkflow, par.PARAM_VIEW_RESULTS);
        searchParams = par.removeParameter(searchParams, par.PARAM_STRUCTTY_MODE);
        searchParams = par.removeParameter(searchParams, par.PARAM_STRUCTTY_SS);
        cmd.addVariable("SEARCH_PAR", par.createParameterString(searchParams, true).c_str());
    }
    cmd.addVariable("SCOREMULTIMER_PAR", par.createParameterString(par.scoremultimer).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("REPORT_PAR", par.createParameterString(par.createmultimerreport).c_str());
    cmd.addVariable("VIEW_RESULTS", par.viewResults ? "TRUE" : NULL);
    // Defer tmp cleanup past the viewer launch: it reads the query/target complex DBs
    // and the viewer report from tmp. Cleanup runs via the CLEANUP_ONLY re-invocation.
    cmd.addVariable("REMOVE_TMP", (par.removeTmpFiles && !par.viewResults) ? "TRUE" : NULL);
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
//    cmd.addVariable("EXP_MULTIMER_PAR", ("-e " + std::to_string(par.eValueThrExpandMultimer)).c_str());
    std::string program = tmpDir + "/multimersearch.sh";
    FileUtil::writeFile(program, multimersearch_sh, multimersearch_sh_len);
    // The other viewer-capable workflows run the script with std::system so control
    // comes back and the viewer can be launched; execProgram would replace this
    // process and never return.
    std::string argString = program;
    for (const auto& s : par.filenames) { argString += " "; argString += s; }
    if (std::system(argString.c_str()) != EXIT_SUCCESS) { EXIT(EXIT_FAILURE); }
    if (par.viewResults) {
        structty::RunOptions opts;
        // query/target are complex DBs (chain-level entries); the 14-column viewer
        // report puts StrucTTY on the multimer path.
        opts.input_files.push_back(queryDb);
        opts.foldseek_target   = targetDb;
        opts.foldseek_result   = tmpDir + "/viewer_report";
        opts.mode              = par.structtyMode;
        opts.show_structure    = par.structtyShowStructure;
        structty::run(opts);
        if (par.removeTmpFiles) {
            cmd.addVariable("CLEANUP_ONLY", "TRUE");
            if (std::system(argString.c_str()) != EXIT_SUCCESS) { EXIT(EXIT_FAILURE); }
        }
    }
    return EXIT_SUCCESS;
}
