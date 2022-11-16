# Author: Milot Mirdita (milot@mirdita.de), Hyunbin Kim (khb7840@gmail.com)
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(foldcomp C CXX)

option(BUILD_LIBRARY "Build library" OFF)
option(BUILD_PYTHON "Build python support" OFF)
option(BUILD_WERROR "Build with warnings as errors" OFF)
option(BUILD_IWYU "Build with include-what-you-use" OFF)
option(GCS_SUPPORT "Enable Google Cloud Storage support" OFF)

include_directories(src)
add_subdirectory(src)

if(HAVE_SANITIZER)
    include(FindUBSan)
    include(FindASan)
    include(FindMSan)
    include(FindTSan)
endif()

if(BUILD_LIBRARY)
    add_library(foldcomp
        ${foldcomp_header_files}
        ${foldcomp_source_files})
elseif(BUILD_PYTHON)
        find_package(PythonExtensions REQUIRED)
        add_library(foldcomp MODULE
            ${foldcomp_header_files}
            ${foldcomp_source_files}
            foldcomp/foldcomp.cxx)
        python_extension_module(foldcomp)
        install(TARGETS foldcomp LIBRARY DESTINATION foldcomp)
else()
    include_directories(lib)

    add_executable(foldcomp
        ${foldcomp_header_files}
        ${foldcomp_source_files}
        ${foldcomp_exe_header_files}
        ${foldcomp_exe_source_files})

    target_compile_definitions(foldcomp PUBLIC FOLDCOMP_EXECUTABLE)
    # For debugging
    # target_compile_definitions(foldcomp PUBLIC _GLIBCXX_DEBUG=1 _LIBCPP_DEBUG=1)

    find_package(OpenMP REQUIRED)
    if(OPENMP_CXX_FOUND)
        if((CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
            target_link_libraries(foldcomp PUBLIC OpenMP::OpenMP_CXX)
        else()
            target_link_libraries(foldcomp PRIVATE "${OpenMP_CXX_FLAGS}")
            target_compile_options(foldcomp PRIVATE "${OpenMP_CXX_FLAGS}")
            target_compile_definitions(foldcomp PUBLIC OPENMP)
        endif()
    endif()

    find_package(ZLIB REQUIRED)
    target_link_libraries(foldcomp PUBLIC ZLIB::ZLIB)

    include_directories(lib/gemmi)

    include_directories(lib/microtar)
    add_subdirectory(lib/microtar)

    target_link_libraries(foldcomp PUBLIC microtar)

    if(GCS_SUPPORT)
        find_package(google_cloud_cpp_storage REQUIRED)
        target_link_libraries(foldcomp PUBLIC google-cloud-cpp::storage)
        target_compile_definitions(foldcomp PUBLIC HAVE_GCS)
    endif()
endif()

# always enable math defines like M_PI
target_compile_definitions(foldcomp PUBLIC _USE_MATH_DEFINES=1)

if(MSVC)
    target_compile_options(foldcomp PRIVATE /W4)
else()
    target_compile_options(foldcomp PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(BUILD_WERROR)
    if(MSVC)
        target_compile_options(foldcomp PRIVATE /WX)
    else()
        target_compile_options(foldcomp PRIVATE -Werror)
    endif()
endif()

if(BUILD_IWYU)
    find_program(iwyu_path NAMES include-what-you-use iwyu REQUIRED)
    set_property(TARGET foldcomp PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path} -Xiwyu --mapping_file=${CMAKE_CURRENT_LIST_DIR}/cmake/include-what-you-use.imp)
endif()

enable_testing()

