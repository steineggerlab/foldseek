set(strucclustutils_source_files
        strucclustutils/tmalign.cpp
        strucclustutils/aln2tmscore.cpp
        strucclustutils/structcreatedb.cpp
        strucclustutils/structurealign.cpp
        strucclustutils/samplemulambda.cpp
        strucclustutils/structureconvertalis.cpp
        strucclustutils/structureto3didescriptor.cpp
        strucclustutils/EvalueNeuralNet.cpp
        strucclustutils/EvalueNeuralNet.h
        strucclustutils/PulchraWrapper.cpp
        strucclustutils/PulchraWrapper.h
        strucclustutils/structurerescorediagonal.cpp
        strucclustutils/convert2pdb.cpp
        strucclustutils/compressca.cpp
        strucclustutils/scoremultimer.cpp
        strucclustutils/filtermultimer.cpp
        strucclustutils/createmultimerreport.cpp
        strucclustutils/MultimerUtil.h
        strucclustutils/expandmultimer.cpp
        strucclustutils/makepaddeddb.cpp
        strucclustutils/result2structprofile.cpp
        strucclustutils/createstructsubdb.cpp
        strucclustutils/createinterfacedb.cpp
        strucclustutils/createinterfacedbs.cpp
        PARENT_SCOPE
        )

add_library(gemmiwrapper
    GemmiWrapper.cpp
    GemmiWrapper.h
    )

get_target_property(COMPILE_TMP mmseqs-framework COMPILE_FLAGS)
get_target_property(LINK_TMP mmseqs-framework LINK_FLAGS)
get_target_property(DEF_TMP mmseqs-framework COMPILE_DEFINITIONS)
get_target_property(INCL_TMP mmseqs-framework INCLUDE_DIRECTORIES)
set_target_properties(gemmiwrapper PROPERTIES COMPILE_FLAGS "${COMPILE_TMP} -w" LINK_FLAGS "${LINK_TMP} -w")
set_property(TARGET gemmiwrapper APPEND PROPERTY COMPILE_DEFINITIONS ${DEF_TMP})
set_property(TARGET gemmiwrapper APPEND PROPERTY INCLUDE_DIRECTORIES ${INCL_TMP})
restore_exceptions(gemmiwrapper)
target_link_libraries(gemmiwrapper foldcomp)
if(HAVE_ZLIB)
    target_link_libraries(gemmiwrapper z)
endif()
target_include_directories(gemmiwrapper PRIVATE ../../lib/gemmi)
target_include_directories(gemmiwrapper PRIVATE ../../lib/foldcomp/src)

set(prostt5_source_files
    strucclustutils/ProstT5.h
    strucclustutils/ProstT5.cpp
    PARENT_SCOPE)