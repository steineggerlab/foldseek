add_library(dssp
        structure.h
        structure.cpp
        utils.cpp
        utils.h
        primitives3d.cpp
        primitives3d.h
        quaternion.h
        quaternion_utils.h
        mas.cpp
        mas.h
        iocif.cpp
        iocif.h)
mmseqs_setup_derived_target(dssp)

