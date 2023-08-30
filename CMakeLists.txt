cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
# keep old policy for setting implicit link libraries
# zlib causes issues in static builds otherwise
cmake_policy(SET CMP0060 OLD)
project(foldseek C CXX)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs/cmake")

set(IGNORE_RUST_VERSION 0 CACHE BOOL "Ignore Rust version check")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
include(MMseqsSetupDerivedTarget)
add_subdirectory(lib/mmseqs EXCLUDE_FROM_ALL)

add_subdirectory(lib/3di)
add_subdirectory(lib/pulchra)
add_subdirectory(lib/kerasify)
add_subdirectory(lib/tmalign)
set(BUILD_LIBRARY 1 CACHE INTERNAL "" FORCE)
add_subdirectory(lib/foldcomp EXCLUDE_FROM_ALL)
add_subdirectory(data)
include_directories(lib)
include_directories(lib/3di)
include_directories(lib/pulchra)
include_directories(lib/kerasify)
include_directories(lib/gemmi)

add_subdirectory(lib/corrosion)
# don't try to link to gcc_s, its not needed anyway
list(REMOVE_ITEM Rust_CARGO_TARGET_LINK_NATIVE_LIBS "gcc_s")

if(Rust_VERSION VERSION_LESS "1.68.0" AND NOT IGNORE_RUST_VERSION)
    message(FATAL_ERROR "Rust >= 1.68.0 is required")
endif()

set(RUST_FEATURE "")
if(HAVE_AVX2)
    set(RUST_FEATURE simd_avx2)
elseif(HAVE_SSE4_1 OR HAVE_SSE2)
    set(RUST_FEATURE simd_sse2)
elseif(HAVE_ARM8)
    set(RUST_FEATURE simd_neon)
endif()

if (NATIVE_ARCH AND (RUST_FEATURE STREQUAL ""))
    if(EMSCRIPTEN)
        set(RUST_FEATURE simd_wasm)
    elseif(ARM)
        set(RUST_FEATURE simd_neon)
    elseif(X86 OR X64)
        if (HAVE_AVX2_EXTENSIONS)
            set(RUST_FEATURE simd_avx2)
        elseif (HAVE_SSE2_EXTENSIONS)
            set(RUST_FEATURE simd_sse2)
        endif()
    endif()
endif()

if(RUST_FEATURE STREQUAL "")
    message(FATAL_ERROR "No SIMD support for this architecture")
endif()
message("-- Rust Feature: ${RUST_FEATURE}")
corrosion_import_crate(
    MANIFEST_PATH lib/block-aligner/c/Cargo.toml
    FEATURES "${RUST_FEATURE}"
    CRATE_TYPES staticlib
)
include_directories(lib/block-aligner/c)
if(EMSCRIPTEN)
    corrosion_add_target_local_rustflags(block-aligner-c "-Clink-args=--no-entry -sRELOCATABLE=1")
endif()


add_subdirectory(src)
