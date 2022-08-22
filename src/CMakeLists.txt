include_directories(commons)

add_subdirectory(strucclustutils)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(workflow)
add_subdirectory(test)

add_library(foldseek-framework
        ${commons_source_files}
        ${strucclustutils_source_files}
        ${workflow_source_files})

mmseqs_setup_derived_target(foldseek-framework)

add_executable(foldseek
        LocalCommandDeclarations.h
        foldseek.cpp)

mmseqs_setup_derived_target(foldseek)
target_link_libraries(foldseek version gemmiwrapper 3di pulchra kerasify tmalign foldseek-framework)

add_dependencies(foldseek local-generated)

install(TARGETS foldseek DESTINATION bin)
