cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(foldseek C CXX)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs/cmake")

set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
include(MMseqsSetupDerivedTarget)
add_subdirectory(lib/mmseqs EXCLUDE_FROM_ALL)

add_subdirectory(lib/dssp)
add_subdirectory(lib/3di)
add_subdirectory(lib/tmalign)
add_subdirectory(data)
include_directories(lib)
include_directories(lib/3di)
include_directories(lib/gemmi)

add_subdirectory(src)
