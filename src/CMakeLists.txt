set(HAVE_GCS 0 CACHE BOOL "Have Google Cloud Storage SDK")
set(HAVE_F16C 0 CACHE BOOL "Have CPU with F16C")

include_directories(commons)
add_subdirectory(strucclustutils)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(foldseek
        ${commons_source_files}
        ${strucclustutils_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        foldseek.cpp)

mmseqs_setup_derived_target(foldseek)
target_link_libraries(foldseek version gemmiwrapper 3di pulchra kerasify tmalign foldcomp)

if(HAVE_F16C)
        include(AppendTargetProperty)
        append_target_property(foldseek COMPILE_FLAGS "-mf16c")
        append_target_property(foldseek LINK_FLAGS "-mf16c")
endif()

if(HAVE_GCS)
        find_package(google_cloud_cpp_storage REQUIRED)
        target_link_libraries(foldseek google-cloud-cpp::storage)
endif()

add_dependencies(foldseek local-generated)

install(TARGETS foldseek DESTINATION bin)
