cmake_minimum_required(VERSION 3.7)
project(ctranslate2 C CXX)

option(WITH_MKL "Compile with Intel MKL backend" ON)
option(WITH_DNNL "Compile with DNNL backend" OFF)
option(WITH_ACCELERATE "Compile with Accelerate backend" OFF)
option(WITH_OPENBLAS "Compile with OpenBLAS backend" OFF)
option(WITH_RUY "Compile with ruy" OFF)
option(WITH_CUDA "Compile with CUDA backend" OFF)
option(WITH_CUDNN "Compile with cuDNN backend" OFF)
option(WITH_FLASHATTENTION "Compile with flash-attention" ON)
option(WITH_SPDLOG "Compile with spdlog" ON)
# option(WITH_LLAMAFILE_SGEMM "Compile with llammafile sgemm backend" ON)
option(CUDA_DYNAMIC_LOADING "Dynamically load CUDA libraries at runtime" OFF)
option(ENABLE_CPU_DISPATCH "Compile CPU kernels for multiple ISA and dispatch at runtime" ON)
option(ENABLE_PROFILING "Compile with profiling support" OFF)
option(BUILD_CLI "Compile the clients" ON)
option(BUILD_TESTS "Compile the tests" OFF)
option(WITH_TENSOR_PARALLEL "Compile with NCCL and MPI backend" OFF)
option(BUILD_PROSTT5 "Compile with ProstT5" OFF)
option(ENABLE_WALL "Compile with -Wall" ON)

if(WITH_CUDA)
  enable_language(CUDA)
  if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 17)
  endif()
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

if(ENABLE_PROFILING)
  message(STATUS "Enable profiling support")
  add_definitions(-DCT2_ENABLE_PROFILING)
endif()

set(PRIVATE_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/mock
  )

set(FLASH_ATTENTION_SOURCES
  src/layers/flash_attention.cc
  src/ops/flash_attention.cc
  src/ops/flash_attention_cpu.cc
  )
set(SOURCES
  src/allocator.cc
  src/batch_reader.cc
  src/buffered_translation_wrapper.cc
  src/cpu/allocator.cc
  src/cpu/backend.cc
  src/cpu/cpu_info.cc
  src/cpu/cpu_isa.cc
  src/cpu/kernels.cc
  src/cpu/parallel.cc
  src/cpu/primitives.cc
  src/decoding.cc
  src/decoding_utils.cc
  src/devices.cc
  src/dtw.cc
  src/encoder.cc
  src/env.cc
  src/filesystem.cc
  src/generator.cc
  src/layers/attention_layer.cc
  src/layers/attention.cc
  src/layers/common.cc
  src/layers/decoder.cc
  src/layers/prostt5cnn.cc
  src/layers/transformer.cc
  src/layers/wav2vec2.cc
  src/layers/whisper.cc
  src/logging.cc
  src/models/language_model.cc
  src/models/model.cc
  src/models/model_factory.cc
  src/models/model_reader.cc
  src/models/sequence_to_sequence.cc
  src/models/transformer.cc
  src/models/wav2vec2.cc
  src/models/whisper.cc
  src/ops/activation.cc
  src/ops/add.cc
  src/ops/alibi_add.cc
  src/ops/alibi_add_cpu.cc
  src/ops/bias_add.cc
  src/ops/bias_add_cpu.cc
  src/ops/concat.cc
  src/ops/concat_split_slide_cpu.cc
  src/ops/conv1d.cc
  src/ops/conv1d_cpu.cc
  src/ops/conv2d.cc
  src/ops/conv2d_cpu.cc
  src/ops/cos.cc
  src/ops/dequantize.cc
  src/ops/dequantize_cpu.cc
  src/ops/gather.cc
  src/ops/gather_cpu.cc
  src/ops/gelu.cc
  src/ops/gemm.cc
  src/ops/gumbel_max.cc
  src/ops/gumbel_max_cpu.cc
  src/ops/layer_norm.cc
  src/ops/layer_norm_cpu.cc
  src/ops/log.cc
  src/ops/matmul.cc
  src/ops/mean.cc
  src/ops/mean_cpu.cc
  src/ops/median_filter.cc
  src/ops/min_max.cc
  src/ops/mul.cc
  src/ops/multinomial.cc
  src/ops/multinomial_cpu.cc
  src/ops/quantize.cc
  src/ops/quantize_cpu.cc
  src/ops/relu.cc
  src/ops/rms_norm.cc
  src/ops/rms_norm_cpu.cc
  src/ops/rotary.cc
  src/ops/rotary_cpu.cc
  src/ops/sin.cc
  src/ops/softmax.cc
  src/ops/softmax_cpu.cc
  src/ops/split.cc
  src/ops/slide.cc
  src/ops/sub.cc
  src/ops/swish.cc
  src/ops/tanh.cc
  src/ops/tile.cc
  src/ops/tile_cpu.cc
  src/ops/topk.cc
  src/ops/topk_cpu.cc
  src/ops/topp_mask.cc
  src/ops/topp_mask_cpu.cc
  src/ops/transpose.cc
  src/ops/nccl_ops.cc
  src/ops/nccl_ops_cpu.cc
  src/padder.cc
  src/profiler.cc
  src/random.cc
  src/sampling.cc
  src/scoring.cc
  src/storage_view.cc
  src/thread_pool.cc
  src/translator.cc
  src/types.cc
  src/utils.cc
  src/vocabulary.cc
  src/vocabulary_map.cc
)
if(WITH_FLASHATTENTION)
  list(APPEND SOURCES ${FLASH_ATTENTION_SOURCES})
endif()

if(DEFINED ARM)
  add_definitions(-DCT2_ARM64_BUILD)
  set(CT2_BUILD_ARCH "arm64")
elseif(DEFINED X86 OR DEFINED X64)
  set(CT2_BUILD_ARCH "x86_64")
endif()

# macro(ct2_compile_kernels_for_isa isa flag)
#   configure_file(
#     src/cpu/kernels.cc
#     ${CMAKE_CURRENT_BINARY_DIR}/kernels_${isa}.cc
#     COPYONLY)
#   set_source_files_properties(
#     ${CMAKE_CURRENT_BINARY_DIR}/kernels_${isa}.cc
#     PROPERTIES COMPILE_FLAGS ${flag})
#   list(APPEND SOURCES ${CMAKE_CURRENT_BINARY_DIR}/kernels_${isa}.cc)
# endmacro()

# if(ENABLE_CPU_DISPATCH)
#   message(STATUS "Compiling for multiple CPU ISA and enabling runtime dispatch")
#   add_definitions(-DCT2_WITH_CPU_DISPATCH)
#   if(CT2_BUILD_ARCH STREQUAL "x86_64")
#     if(WIN32)
#       ct2_compile_kernels_for_isa(avx "/arch:AVX")
#       ct2_compile_kernels_for_isa(avx2 "/arch:AVX2")
#       ct2_compile_kernels_for_isa(avx512 "/arch:AVX512")
#     else()
#       ct2_compile_kernels_for_isa(avx "-mavx")
#       ct2_compile_kernels_for_isa(avx2 "-mavx2 -mfma")
#       ct2_compile_kernels_for_isa(avx512 "-mavx512f -mavx512cd -mavx512vl -mavx512bw -mavx512dq")
#     endif()
#   elseif(CT2_BUILD_ARCH STREQUAL "arm64")
#     ct2_compile_kernels_for_isa(neon "-DUSE_NEON")
#   endif()
# endif()

# if(DEFINED OpenMP_CXX_LIBRARIES)
#   list(APPEND LIBRARIES ${OpenMP_CXX_LIBRARIES})
# endif()

# if(WITH_MKL)
#   find_path(MKL_ROOT include/mkl.h DOC "Path to MKL root directory" PATHS
#     $ENV{MKLROOT}
#     ${INTEL_ROOT}/mkl
#     ${INTEL_ROOT}/oneAPI/mkl/latest
#     ${INTEL_ROOT}/oneapi/mkl/latest
#     )

#   # Find MKL includes.
#   find_path(MKL_INCLUDE_DIR NAMES mkl.h HINTS ${MKL_ROOT}/include/)
#   if(MKL_INCLUDE_DIR)
#     message(STATUS "Found MKL include directory: ${MKL_INCLUDE_DIR}")
#   else()
#     message(FATAL_ERROR "MKL include directory not found")
#   endif()

#   # Find MKL libraries.
#   find_library(MKL_CORE_LIBRARY NAMES mkl_core HINTS ${MKL_ROOT}/lib ${MKL_ROOT}/lib/intel64)
#   if(MKL_CORE_LIBRARY)
#     get_filename_component(MKL_LIBRARY_DIR ${MKL_CORE_LIBRARY} DIRECTORY)
#     message(STATUS "Found MKL library directory: ${MKL_LIBRARY_DIR}")
#   else()
#     message(FATAL_ERROR "MKL library directory not found")
#   endif()

#   add_definitions(-DCT2_WITH_MKL -DMKL_ILP64)
#   set(MKL_LIBRARIES
#     ${MKL_LIBRARY_DIR}/libmkl_core.a
#     ${MKL_LIBRARY_DIR}/libmkl_intel_ilp64.a
#   )

#   if(DEFINED OpenMP_CXX_LIBRARIES)
#     list(APPEND MKL_LIBRARIES ${MKL_LIBRARY_DIR}/libmkl_gnu_thread.a)
#   else()
#     list(APPEND MKL_LIBRARIES ${MKL_LIBRARY_DIR}/libmkl_sequential.a)
#   endif()
#   list(APPEND PRIVATE_INCLUDE_DIRECTORIES ${MKL_INCLUDE_DIR})
#   if(WIN32 OR APPLE)
#     list(APPEND LIBRARIES ${MKL_LIBRARIES})
#   else()
#     list(APPEND LIBRARIES -Wl,--start-group ${MKL_LIBRARIES} -Wl,--end-group)
#   endif()
# endif()

# if(WITH_DNNL)
#   set(ONEAPI_DNNL_PATH ${INTEL_ROOT}/oneapi/dnnl/latest)
#   if(OPENMP_RUNTIME STREQUAL "INTEL")
#     set(ONEAPI_DNNL_PATH ${ONEAPI_DNNL_PATH}/cpu_iomp)
#   else()
#     set(ONEAPI_DNNL_PATH ${ONEAPI_DNNL_PATH}/cpu_gomp)
#   endif()

#   find_path(DNNL_INCLUDE_DIR NAMES dnnl.h PATHS ${ONEAPI_DNNL_PATH}/include)
#   if(DNNL_INCLUDE_DIR)
#     message(STATUS "Found DNNL include directory: ${DNNL_INCLUDE_DIR}")
#   else()
#     message(FATAL_ERROR "DNNL include directory not found")
#   endif()

#   find_library(DNNL_LIBRARY NAMES dnnl PATHS ${ONEAPI_DNNL_PATH}/lib)
#   if(DNNL_LIBRARY)
#     message(STATUS "Found DNNL library: ${DNNL_LIBRARY}")
#   else()
#     message(FATAL_ERROR "DNNL library not found")
#   endif()

#   add_definitions(-DCT2_WITH_DNNL)
#   list(APPEND PRIVATE_INCLUDE_DIRECTORIES ${DNNL_INCLUDE_DIR})
#   list(APPEND LIBRARIES ${DNNL_LIBRARY})
# endif()

# if (WITH_ACCELERATE)
#   set(BLA_VENDOR Apple)
#   find_package(BLAS REQUIRED)
#   add_definitions(-DCT2_WITH_ACCELERATE)
#   list(APPEND LIBRARIES ${BLAS_LIBRARIES})
# endif()

if (WITH_OPENBLAS)
  find_path(OPENBLAS_INCLUDE_DIR NAMES cblas.h)
  if(OPENBLAS_INCLUDE_DIR)
    message(STATUS "Found OpenBLAS include directory: ${OPENBLAS_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR "OpenBLAS include directory not found")
  endif()

  find_library(OPENBLAS_LIBRARY NAMES openblas)
  if(OPENBLAS_LIBRARY)
    message(STATUS "Found OpenBLAS library: ${OPENBLAS_LIBRARY}")
  else()
    message(FATAL_ERROR "OpenBLAS library not found")
  endif()

  add_definitions(-DCT2_WITH_OPENBLAS)
  list(APPEND PRIVATE_INCLUDE_DIRECTORIES ${OPENBLAS_INCLUDE_DIR})
  list(APPEND LIBRARIES ${OPENBLAS_LIBRARY})
endif()

if (WITH_RUY)
  add_definitions(-DCT2_WITH_RUY)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  # set(CPUINFO_LIBRARY_TYPE static CACHE STRING "cpuinfo library type")
  add_subdirectory(third_party/ruy EXCLUDE_FROM_ALL)
  unset(CMAKE_POSITION_INDEPENDENT_CODE)
  list(APPEND LIBRARIES ruy)
endif()

# if (WITH_LLAMAFILE_SGEMM)
#   add_definitions(-DCT2_WITH_LLAMAFILE_SGEMM)
#   # list(APPEND PRIVATE_INCLUDE_DIRECTORIES third_party/llamafile)
#   list(APPEND SOURCES third_party/llamafile/sgemm.cpp)
# endif()

if(WITH_CUDA)
  add_definitions(-DCT2_WITH_CUDA)
  list(APPEND CUDA_NVCC_FLAGS "-Xfatbin=-compress-all")

  # We should ensure that the Thrust include directories appear before
  # -I/usr/local/cuda/include for both GCC and NVCC, so that the headers
  # are coming from the submodule and not the system.
  set(THRUST_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/thrust/dependencies/cub
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/thrust
    )
  add_definitions(-DTHRUST_IGNORE_CUB_VERSION_CHECK=1)
  list(APPEND PRIVATE_INCLUDE_DIRECTORIES ${THRUST_INCLUDE_DIRS})

  list(APPEND SOURCES
    src/cuda/allocator.cc
    src/cuda/primitives.cu
    src/cuda/random.cu
    src/cuda/utils.cc
    src/ops/alibi_add_gpu.cu
    src/ops/bias_add_gpu.cu
    src/ops/concat_split_slide_gpu.cu
    src/ops/conv1d_gpu.cu
    src/ops/conv2d_gpu.cu
    src/ops/dequantize_gpu.cu
    src/ops/gather_gpu.cu
    src/ops/gumbel_max_gpu.cu
    src/ops/layer_norm_gpu.cu
    src/ops/mean_gpu.cu
    src/ops/multinomial_gpu.cu
    src/ops/rms_norm_gpu.cu
    src/ops/rotary_gpu.cu
    src/ops/softmax_gpu.cu
    src/ops/tile_gpu.cu
    src/ops/topk_gpu.cu
    src/ops/topp_mask_gpu.cu
    src/ops/quantize_gpu.cu
    src/ops/nccl_ops_gpu.cu)

  if(WITH_FLASHATTENTION)
    list(APPEND CUDA_NVCC_FLAGS "--expt-relaxed-constexpr")
    list(APPEND CUDA_NVCC_FLAGS "--expt-extended-lambda")
    set(CUTLASS_INCLUDE_DIRS
      ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cutlass/include
    )
    list(APPEND PRIVATE_INCLUDE_DIRECTORIES ${CUTLASS_INCLUDE_DIRS})
    add_definitions(-DCT2_WITH_FLASHATTENTION)

    list(APPEND SOURCES 
      src/ops/flash_attention_gpu.cu
      src/ops/flash-attention/flash_fwd_hdim32_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim32_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim64_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim64_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim96_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim96_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim128_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim128_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim160_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim160_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim192_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim192_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim224_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim224_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim256_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_hdim256_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim32_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim32_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim64_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim64_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim96_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim96_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim128_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim128_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim160_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim160_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim192_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim192_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim224_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim224_fp16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim256_bf16_sm80.cu
      src/ops/flash-attention/flash_fwd_split_hdim256_fp16_sm80.cu
    )
  endif()
endif()

add_library(${PROJECT_NAME} STATIC ${SOURCES})
target_compile_options(
  ${PROJECT_NAME} PRIVATE
  ${MMSEQS_ARCH}
  $<$<COMPILE_LANGUAGE:CUDA>:
      ${CUDA_NVCC_FLAGS}
      >
)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
target_link_libraries(${PROJECT_NAME} PRIVATE ${LIBRARIES})
target_include_directories(${PROJECT_NAME} BEFORE
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>
  PRIVATE ${PRIVATE_INCLUDE_DIRECTORIES}
)

add_subdirectory(prostt5)