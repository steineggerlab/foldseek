add_library(tmalign
        TMalign.h
        TMalign.cpp
        NW.h
        Kabsch.h
        se.cpp
        se.h
        xyz_sfetch.cpp
        pstream.h
        pdb2fasta.cpp
        pdb2xyz.cpp
        param_set.h)
mmseqs_setup_derived_target(tmalign)

