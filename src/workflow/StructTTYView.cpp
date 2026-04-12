#include <cstdlib>
#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
namespace structtyViewerStandalone {
#include "structty_viewer.sh.h"
}

int structtyview(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string queryFile = par.filenames[0];
    std::string resultM8 = par.filenames[1];
    std::string targetDb;
    if (par.filenames.size() > 2) {
        targetDb = par.filenames[2];
    }

    // Write structty_viewer.sh to a temporary location
    std::string structtyViewer = "/tmp/structty_viewer_" + SSTR(getpid()) + ".sh";
    FileUtil::writeFile(structtyViewer, structtyViewerStandalone::structty_viewer_sh, structtyViewerStandalone::structty_viewer_sh_len);

    std::string viewerCmd = "sh \"" + structtyViewer + "\""
        + " \"" + par.structtyPath + "\""
        + " \"" + queryFile + "\""
        + " \"" + resultM8 + "\""
        + " \"" + targetDb + "\"";

    int ret = system(viewerCmd.c_str());

    // Clean up
    std::remove(structtyViewer.c_str());

    return ret == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
