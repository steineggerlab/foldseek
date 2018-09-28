include_directories(commons)
add_subdirectory(strucclust)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(plass
        ${commons_source_files}
        ${strucclust_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        plass.cpp
        )

mmseqs_setup_derived_target(strucclust)
target_link_libraries(strucclust version)

add_dependencies(plass local-generated)

install(TARGETS plass DESTINATION bin)
