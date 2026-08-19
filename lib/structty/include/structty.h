#pragma once
#include <string>
#include <vector>

namespace structty {

struct RunOptions {
    std::vector<std::string> input_files;
    std::string mode        = "protein";
    bool show_structure     = false;
    bool no_panel           = false;
    bool benchmark          = false;

    std::string chains_file;
    std::string msa_file;
    std::string foldmason_file;

    // Foldseek DB, structure directory, structure file or "auto" to download
    std::string foldseek_target;
    // Foldseek m8 (12/17/21/29 columns) or multimer report (14 columns)
    std::string foldseek_result;
};

// Launch the interactive viewer. Blocks until the user presses Q.
// Returns false when the inputs cannot be rendered; the reason is printed.
bool run(const RunOptions& opts);

} // namespace structty
