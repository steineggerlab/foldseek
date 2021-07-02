add_library(tmalign
        affineneedlemanwunsch.cpp
        affineneedlemanwunsch.h
        basic_fun.h
        TMalign.h
        TMalign.cpp
        NW.h
        NW.cpp
        Kabsch.h
        se.h
        xyz_sfetch.cpp
        pstream.h
        param_set.h
        Coordinates.h
        )


mmseqs_setup_derived_target(tmalign)

