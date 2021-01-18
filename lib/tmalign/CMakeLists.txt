add_library(tmalign
        affineneedlemanwunsch.cpp
        affineneedlemanwunsch.h
        TMalign.h
        NW.h
        Kabsch.h
        se.h
        xyz_sfetch.cpp
        pstream.h
        pdb2fasta.cpp
        pdb2xyz.cpp
        param_set.h
        Coordinates.h
        )


mmseqs_setup_derived_target(tmalign)

