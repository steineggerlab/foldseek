#include "ctranslate2/ops/conv2d.h"
#include "ctranslate2/ops/gemm.h"

#include "type_dispatch.h"
#include "cuda/helpers.h"
#include "cuda/utils.h"

// #include <iostream>

namespace ctranslate2 {
  namespace ops {

    template <typename T>
    __global__ void add_bias_kernel(const dim_t size, const dim_t channels, const dim_t output_height, const dim_t output_width, const T* bias, T* output) {
      const dim_t index = blockIdx.x * blockDim.x + threadIdx.x;
      if (index < size) {
        output[index] += bias[(index / (output_height * output_width)) % channels];
      }
    }

    template <Device D, typename T>
    void Conv2D::compute(const StorageView& input,
                         const StorageView& weight,
                         const StorageView* bias,
                         StorageView& output,
                         const StorageView* qscale) const {
      if (qscale) {
        throw std::runtime_error("Quantization is not supported in this Conv2D implementation");
      }

      // Input dimensions: (batch_size, in_channels, input_height, input_width)
      const dim_t batch_size = input.dim(0);
      const dim_t in_channels = input.dim(1);
      const dim_t input_height = input.dim(2);
      const dim_t input_width = input.dim(3);

      // Weight dimensions: (out_channels, in_channels, kernel_height, kernel_width)
      const dim_t out_channels = weight.dim(0);
      const dim_t kernel_height = weight.dim(2);
      const dim_t kernel_width = weight.dim(3);

      // Calculate output dimensions
      const dim_t output_height = (input_height + 2 * _padding_h - (_dilation * (kernel_height - 1) + 1)) / _stride_h + 1;
      const dim_t output_width = (input_width + 2 * _padding_w - (_dilation * (kernel_width - 1) + 1)) / _stride_w + 1;

      // Resize output tensor
      output.resize({batch_size, out_channels, output_height, output_width});

      StorageView im2col_output({batch_size, output_height * output_width, in_channels * kernel_height * kernel_width}, 0.0f, weight.device());
      im2col_transposed_gpu(input, im2col_output, kernel_height, kernel_width);
      // std::cout << "im2col_output: " << im2col_output.to(Device::CPU) << ", " << output_height << ", " << output_width << std::endl;

      StorageView weight_view(weight.dtype(), weight.device());
      weight_view.view(const_cast<void*>(weight.buffer()), {weight.dim(0), in_channels * kernel_height * kernel_width});
      // std::cout << "weight.dtype: " << (int)(weight.dtype()) << " weight.device: " << (int)(weight.device()) << std::endl;
      // std::cout << "weight_view: " << weight_view.to(Device::CPU) << ", " << output_height << ", " << output_width << std::endl;

      const dim_t m = out_channels;
      const dim_t n = output_height * output_width;
      const dim_t k = in_channels * kernel_height * kernel_width;
      const dim_t strideb = k * output_height * output_width;
      const dim_t stridec = out_channels * output_height * output_width;
      auto* b = im2col_output.data<float>();
      auto* c = output.data<float>();

      const Gemm gemm(1.0, 0.0, false, true);
      const auto device = im2col_output.device();
      for (dim_t i = 0; i < batch_size; ++i) {
        float* b_i = b + (i * strideb);
        float* c_i = c + (i * stridec);
        StorageView bb({n, k}, b_i, device); // transposed
        StorageView cc({m, n}, c_i, device);

        // std::cout << "bb: " << bb.to(Device::CPU) << std::endl;
        // std::cout << "cc: " << cc.to(Device::CPU) << std::endl;
        gemm(weight_view, bb, cc);
        // std::cout << "cc: " << cc.to(Device::CPU) << std::endl;

      }
      // std::cout << "gemm: " << output.to(Device::CPU) << std::endl;

      // Add bias
      if (bias) {
        const float* bias_data = bias->data<float>();
        float* output_data = output.data<float>();

        const dim_t size = output.size();
        const dim_t channels = bias->size();
        // std::cout << "size: " << size << " channels: " << channels << " out_channels: " << out_channels << std::endl;

        const dim3 block_dim(256);
        const dim3 grid_dim((size + block_dim.x - 1) / block_dim.x);

        add_bias_kernel<<<grid_dim, block_dim, 0, cuda::get_cuda_stream()>>>(
          size,
          out_channels,
          output_height,
          output_width,
          cuda::device_cast(bias_data),
          cuda::device_cast(output_data)
        );
        cudaDeviceSynchronize();
        // std::cout << "bias: " << bias->to(Device::CPU) << std::endl;
        // std::cout << "output: " << output.to(Device::CPU) << std::endl;
      }
    }

    // im2col kernel adapted from candle:
    // https://github.com/huggingface/candle/blob/6f0b807ffd553fed27325a2a118b0e30bb6d9cbd/candle-kernels/src/conv.cu
    //  Apache-2.0 OR MIT license
    template <typename T>
    __device__ void im2col(
        const size_t dst_numel,
        const size_t h_out,
        const size_t w_out,
        const size_t h_k,
        const size_t w_k,
        const size_t stride_h,
        const size_t stride_w,
        const size_t padding_h,
        const size_t padding_w,
        const size_t dilation,
        const int64_t b_size,
        const int64_t c_in,
        const int64_t h_in,
        const int64_t w_in,
        const T *src,
        T *dst
    ) {
        const size_t dst_i = blockIdx.x * blockDim.x + threadIdx.x;
        // dst: (b_size, h_out, w_out, c_in, h_k, w_k)
        // src: (b_size, c_in, h_in, w_in)
        if (dst_i >= dst_numel) {
            return;
        }
        const size_t src_s0 = c_in * h_in * w_in;
        const size_t src_s1 = h_in * w_in;
        const size_t src_s2 = w_in;
        const size_t src_s3 = 1;

        const size_t dst_s4 = w_k;
        const size_t dst_s3 = h_k * dst_s4;
        const size_t dst_s2 = c_in * dst_s3;
        const size_t dst_s1 = w_out * dst_s2;
        const size_t dst_s0 = h_out * dst_s1;

        size_t tmp_dst_i = dst_i;
        const size_t b_idx = tmp_dst_i / dst_s0;
        tmp_dst_i -= b_idx * dst_s0;
        const size_t h_idx = tmp_dst_i / dst_s1;
        tmp_dst_i -= h_idx * dst_s1;
        const size_t w_idx = tmp_dst_i / dst_s2;
        tmp_dst_i -= w_idx * dst_s2;
        const size_t c_idx = tmp_dst_i / dst_s3;
        tmp_dst_i -= c_idx * dst_s3;
        const size_t h_k_idx = tmp_dst_i / dst_s4;
        tmp_dst_i -= h_k_idx * dst_s4;
        const size_t w_k_idx = tmp_dst_i;

        size_t src_h_idx = h_idx * stride_h + h_k_idx * dilation;
        size_t src_w_idx = w_idx * stride_w + w_k_idx * dilation;
        
        if (src_h_idx < padding_h || src_h_idx >= h_in + padding_h) {
            dst[dst_i] = static_cast<T>(0);
        }
        else if (src_w_idx < padding_w || src_w_idx >= w_in + padding_w) {
            dst[dst_i] = static_cast<T>(0);
        }
        else {
            src_h_idx -= padding_h;
            src_w_idx -= padding_w;
            const size_t src_i =
            b_idx * src_s0
            + c_idx * src_s1
            + src_h_idx * src_s2
            + src_w_idx * src_s3;
            dst[dst_i] = src[src_i];
        }
    }

    template <typename T>
    __global__ void im2col_kernel(
        const size_t dst_numel,
        const size_t h_out,
        const size_t w_out,
        const size_t h_k,
        const size_t w_k,
        const size_t stride_h,
        const size_t stride_w,
        const size_t padding_h,
        const size_t padding_w,
        const size_t dilation,
        const int64_t info0,
        const int64_t info1,
        const int64_t info2,
        const int64_t info3,
        const T *src,
        T *dst) {
        im2col<T>(
            dst_numel,
            h_out,
            w_out,
            h_k,
            w_k,
            stride_h,
            stride_w,
            padding_h,
            padding_w,
            dilation,
            info0,
            info1,
            info2,
            info3,
            src,
            dst);
    }

    void Conv2D::im2col_transposed_gpu(const StorageView& input, StorageView& output, const dim_t kernel_height, const dim_t kernel_width) const {
      const size_t dst_numel = output.size();
      const size_t h_out = (input.dim(2) + 2 * _padding_h - (_dilation * (kernel_height - 1) + 1)) / _stride_h + 1;
      const size_t w_out = (input.dim(3) + 2 * _padding_w - (_dilation * (kernel_width - 1) + 1)) / _stride_w + 1;
      const int64_t *info = input.shape().data();
      const float* src = input.data<float>();
      float* dst = output.data<float>();

      const dim3 block_dim(256);
      const dim3 grid_dim((dst_numel + block_dim.x - 1) / block_dim.x);

      im2col_kernel<<<grid_dim, block_dim, 0, cuda::get_cuda_stream()>>>(
          dst_numel,
          h_out,
          w_out,
          kernel_height,
          kernel_width,
          _stride_h,
          _stride_w,
          _padding_h,
          _padding_w,
          _dilation,
          info[0],
          info[1],
          info[2],
          info[3],
          src,
          dst);
      cudaError_t err = cudaGetLastError();
      if (err != cudaSuccess) {
          std::cerr << "CUDA error after kernel launch: " << cudaGetErrorString(err) << std::endl;
          abort();
      }

      cudaDeviceSynchronize();
      err = cudaGetLastError();
      if (err != cudaSuccess) {
          std::cerr << "CUDA error after device synchronization: " << cudaGetErrorString(err) << std::endl;
          abort();
      }
    }

#define DECLARE_IMPL(T)                                                 \
    template void                                                       \
    Conv2D::compute<Device::CUDA, T>(const StorageView& input,          \
                                     const StorageView& weight,         \
                                     const StorageView* bias,           \
                                     StorageView& output,               \
                                     const StorageView* qscale) const;

    DECLARE_IMPL(float)
    DECLARE_IMPL(float16_t)
    DECLARE_IMPL(bfloat16_t)

  }
}
