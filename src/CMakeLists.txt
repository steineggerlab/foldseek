set(HAVE_GCS 0 CACHE BOOL "Have Google Cloud Storage SDK")

include_directories(commons)
add_subdirectory(strucclustutils)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(workflow)

add_library(foldseek-framework
        ${commons_source_files}
        ${strucclustutils_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        FoldseekBase.cpp
)
mmseqs_setup_derived_target(foldseek-framework)
target_link_libraries(foldseek-framework version gemmiwrapper 3di pulchra kerasify tmalign block-aligner-c)
add_dependencies(foldseek-framework local-generated)

if(HAVE_GCS)
        find_package(google_cloud_cpp_storage REQUIRED)
        target_link_libraries(foldseek-framework google-cloud-cpp::storage)
endif()

if (NOT FOLDSEEK_FRAMEWORK_ONLY)
        add_executable(foldseek foldseek.cpp)
        mmseqs_setup_derived_target(foldseek foldseek-framework)
        install(TARGETS foldseek DESTINATION bin)
endif()