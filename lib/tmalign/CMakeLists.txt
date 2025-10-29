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
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-ffp-contract=off" HAS_FFP_CONTRACT_OFF)
if (HAS_FFP_CONTRACT_OFF)
  # automatic FMA causes different values on different architectures
  append_target_property(tmalign COMPILE_FLAGS -w -ffp-contract=off)
else()
  append_target_property(tmalign COMPILE_FLAGS -w)
endif()
append_target_property(tmalign LINK_FLAGS -w)
