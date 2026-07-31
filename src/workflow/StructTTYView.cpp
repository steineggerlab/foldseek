#include <dirent.h>
#include <fstream>
#include "structty.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"

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

// The viewer needs 3D coordinates, so reject inputs that cannot supply any before
// anything else happens. Keep this identical in EasyStructureSearch, StructureSearch,
// EasyMultimerSearch and MultimerSearch.
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
                Debug(Debug::ERROR) << "The StrucTTY viewer needs C-alpha coordinates, but the database "
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
            Debug(Debug::ERROR) << "The StrucTTY viewer cannot render " << input
                                << ": it is a sequence FASTA and carries no 3D coordinates.\n"
                                << "Pass PDB/mmCIF structures instead. createdb --prostt5-model "
                                << "predicts 3Di from sequence, but writes no C-alpha coordinates.\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

int structtyview(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    validateStructtyMode(par.structtyMode);

    // --foldseek-target and --foldseek-result are a pair: a target source without a
    // result has nothing to look up, and a result without one leaves the hit
    // structures unreachable. Mirrors StrucTTY's own -fst/-fsr rule.
    if (par.foldseekTarget.empty() != par.foldseekResult.empty()) {
        Debug(Debug::ERROR) << "--foldseek-target and --foldseek-result must be given together.\n"
                            << "Use --foldseek-target auto to download the hit structures.\n";
        EXIT(EXIT_FAILURE);
    }

    std::vector<std::string> toValidate = par.filenames;
    // "auto" is a reserved value, not a path.
    if (par.foldseekTarget.empty() == false && par.foldseekTarget != "auto") {
        toValidate.push_back(par.foldseekTarget);
    }
    validateStructtyInputs(toValidate);

    structty::RunOptions opts;
    opts.input_files     = par.filenames;        // query: structures or a Foldseek DB
    opts.foldseek_target = par.foldseekTarget;
    opts.foldseek_result = par.foldseekResult;
    opts.mode            = par.structtyMode;
    opts.show_structure  = par.structtyShowStructure;

    structty::run(opts);
    return EXIT_SUCCESS;
}
