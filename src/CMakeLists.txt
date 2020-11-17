include_directories(commons)
add_subdirectory(strucclustutils)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(strucclust
        ${commons_source_files}
        ${strucclustutils_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        strucclust.cpp
        )

mmseqs_setup_derived_target(strucclust)
target_link_libraries(strucclust version dssp pdbtools)

add_dependencies(strucclust local-generated)

install(TARGETS strucclust DESTINATION bin)
