#include "LocalParameters.h"

#ifdef HAVE_STRUCTTY
#include "structty.h"

static const char *STRUCTTY_MODES[] = {
        "protein", "chain", "rainbow", "plddt", "interface", "conservation", "align",
        "align-fs", "align-near"
};

int structtyview(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    structty::RunOptions opts;
    opts.input_files.push_back(par.db1);
    opts.foldseek_target = par.db2;
    opts.foldseek_result = par.db3;
    opts.mode           = STRUCTTY_MODES[par.structtyMode];
    opts.show_structure = par.structtyShowStructure;

    return structty::run(opts) ? EXIT_SUCCESS : EXIT_FAILURE;
}
#else
#include "Debug.h"

int structtyview(int, const char **, const Command &) {
    Debug(Debug::ERROR) << "This foldseek was built without the StrucTTY viewer. "
                        << "Rebuild with -DENABLE_STRUCTTY=1.\n";
    return EXIT_FAILURE;
}
#endif
