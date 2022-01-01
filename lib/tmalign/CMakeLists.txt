add_library(tmalign
        affineneedlemanwunsch.cpp
        affineneedlemanwunsch.h
        basic_fun.h
        TMalign.h
        TMalign.cpp
        Kabsch.h
        Coordinates.h
        )


mmseqs_setup_derived_target(tmalign)

