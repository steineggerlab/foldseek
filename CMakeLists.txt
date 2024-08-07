cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
# keep old policy for setting implicit link libraries
# zlib causes issues in static builds otherwise
cmake_policy(SET CMP0060 OLD)
project(foldseek C CXX)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs/cmake")

set(ENABLE_PROSTT5 1 CACHE BOOL "Enable ProstT5")
set(ENABLE_CUDA 0 CACHE BOOL "Enable CUDA")
set(IGNORE_RUST_VERSION 0 CACHE BOOL "Ignore Rust version check")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
include(MMseqsSetupDerivedTarget)
add_subdirectory(lib/mmseqs EXCLUDE_FROM_ALL)

set(FOLDSEEK_FRAMEWORK_ONLY 0 CACHE BOOL "Framework mode (don't create foldseek executable)")
if (FOLDSEEK_FRAMEWORK_ONLY)
    set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
endif()
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

set(ENV{CARGO_NET_OFFLINE} true)
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
    corrosion_add_target_local_rustflags(block_aligner_c "-Clink-args=--no-entry -sRELOCATABLE=1")
endif()

set(CANDLE_FEATURE "")
if(ENABLE_CUDA)
    if(NOT DEFINED CUDAToolkit_ROOT)
        message(FATAL_ERROR "please set -DCUDAToolkit_ROOT=path-to-cuda")
    endif()
    set(CANDLE_FEATURE cuda)
endif()

if(ENABLE_PROSTT5)
    corrosion_import_crate(
        MANIFEST_PATH lib/prostt5/c/Cargo.toml
        CRATE_TYPES staticlib
        FEATURES "${CANDLE_FEATURE}"
    )
    if(ENABLE_CUDA)
        corrosion_set_env_vars(cprostt5 CUDA_ROOT=${CUDAToolkit_ROOT})
    endif()
    include_directories(lib/prostt5/c)
    if(EMSCRIPTEN)
        corrosion_add_target_local_rustflags(cprostt5 "-Clink-args=--no-entry -sRELOCATABLE=1")
    endif()
endif()

add_subdirectory(src)
