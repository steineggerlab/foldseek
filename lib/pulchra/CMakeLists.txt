
add_library(pulchra
        pulchra.cpp
        pulchra.h
        pulchra_common.h
        nco_data.h
        )

mmseqs_setup_derived_target(pulchra)
