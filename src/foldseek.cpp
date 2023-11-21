#include <cstddef>
#include "Command.h"
#include "LocalParameters.h"

const char* binary_name = "foldseek";
const char* tool_name = "foldseek";
const char* tool_introduction = "Foldseek enables fast and sensitive comparisons of large structure sets. It reaches sensitivities similar to state-of-the-art structural aligners while being at least 20,000 times faster.\n\nPlease cite:\nvan Kempen, M., Kim, S.S., Tumescheit, C., Mirdita, M., Lee, J., Gilchrist, C.L.M., Söding, J., and Steinegger, M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology, doi:10.1038/s41587-023-01773-0 (2023)";
const char* main_author = "Michel van Kempen, Stephanie Kim, Charlotte Tumescheit, Milot Mirdita, Jeongjae Lee, Cameron L. M. Gilchrist, Johannes Söding, Martin Steinegger";
const char* show_extended_help = "1";
const char* show_bash_info = NULL;
const char* index_version_compatible = "fs1";
bool hide_base_commands = true;
bool hide_base_downloads = true;

extern std::vector<Command> baseCommands;
extern std::vector<Command> foldseekCommands;
void init() {
    registerCommands(&baseCommands);
    registerCommands(&foldseekCommands);
}
void (*initCommands)(void) = init;
void initParameterSingleton() { new LocalParameters; }