#include "ctranslate2/ops/conv2d.h"

#include "ctranslate2/ops/gemm.h"
#include "cpu/parallel.h"
#include "ctranslate2/ops/quantize.h"
#include "ctranslate2/ops/dequantize.h"

namespace ctranslate2 {
  namespace ops {

    template<>
    void
    Conv2D::compute<Device::CPU, float>(const StorageView &input, const StorageView &weight, const StorageView *bias,
                                        StorageView &output, const StorageView *qscale) const {
      if (_dilation != 1)
        throw std::runtime_error("Dilation is not supported in this Conv2D implementation");

      const dim_t batch_size = input.dim(0);
      const dim_t in_channels = input.dim(1);
      const dim_t input_height = input.dim(2);
      const dim_t input_width = input.dim(3);
      const dim_t out_channels = weight.dim(0);
      const dim_t kernel_height = weight.dim(2);
      const dim_t kernel_width = weight.dim(3);

      // Correct output height and width calculation
      const dim_t output_height = (input_height + 2 * _padding_h - (_dilation * (kernel_height - 1) + 1)) / _stride_h + 1;
      const dim_t output_width = (input_width + 2 * _padding_w - (_dilation * (kernel_width - 1) + 1)) / _stride_w + 1;

      output.resize({batch_size, out_channels, output_height, output_width});

      compute_with_gemm(input, weight, output, qscale);
      
      // Add bias
      if (bias) {
        const auto a = bias->data<float>();
        const auto b = output.data<float>();
        cpu::parallel_for(0, batch_size * out_channels, 1, [&](dim_t begin, dim_t end) {
          for (dim_t i = begin; i < end; ++i) {
            const auto a_i = a[i % out_channels];
            const auto b_i = b + i * output_height * output_width;
            primitives<>::add(a_i, b_i, b_i, output_height * output_width);
          }
        });
      }
    }

    void Conv2D::compute_with_gemm(const StorageView &input, const StorageView &weight, StorageView &output,
                                   const StorageView *qscale) const {
      const dim_t batch_size = input.dim(0);
      const dim_t in_channels = input.dim(1);
      const dim_t out_channels = weight.dim(0);
      const dim_t kernel_height = weight.dim(2);
      const dim_t kernel_width = weight.dim(3);
      const dim_t output_height = output.dim(2);
      const dim_t output_width = output.dim(3);

      StorageView im2col_output({batch_size, output_height * output_width, in_channels * kernel_height * kernel_width}, 0.0f, weight.device());
      im2col_transposed(input, im2col_output, kernel_height, kernel_width);
      StorageView weight_view(weight.dtype(), weight.device());
      weight_view.view(const_cast<void*>(weight.buffer()), {weight.dim(0), in_channels * kernel_height * kernel_width});

      const dim_t m = out_channels;
      const dim_t n = output_height * output_width;
      const dim_t k = in_channels * kernel_height * kernel_width;
      const dim_t strideb = k * output_height * output_width;
      const dim_t stridec = out_channels * output_height * output_width;
      auto* b = im2col_output.data<float>();
      auto* c = output.data<float>();
      const Gemm gemm(1.0, 0.0, false, true);
      const Quantize quantize_op(Quantize::ScaleType::PER_LAYER,
                                 /*shift_to_uint8=*/false,
                                 /*round_before_cast=*/true);
      const Dequantize dequantize_op;
      const auto device = im2col_output.device();
      cpu::parallel_for(0, batch_size, 1, [&](dim_t begin, dim_t end) {
        StorageView qinput(weight.dtype(), device);
        StorageView qinput_scale(device);
        if (qscale)
          qinput_scale.to(qscale->dtype());
        StorageView qoutput(DataType::INT32, device);
        for (dim_t i = begin; i < end; ++i) {
          float* b_i = b + (i * strideb);
          float* c_i = c + (i * stridec);
          StorageView bb({n, k}, b_i); // transposed
          StorageView cc({m, n}, c_i);

          if (qscale) {
            quantize_op(bb, qinput, qinput_scale);
            gemm(weight_view, qinput, qoutput);
            dequantize_op(qoutput,
                          *qscale,
                          qinput_scale,
                          /*trans_a=*/false,
                          /*trans_b=*/true,
                          cc);
          } else {
            gemm(weight_view, bb, cc);
          }
        }
      });
    }

    void Conv2D::im2col_transposed(const StorageView& input, StorageView& output, const dim_t kernel_height, const dim_t kernel_width) const {
      const dim_t batch_size = input.dim(0);
      const dim_t in_channels = input.dim(1);
      const dim_t input_height = input.dim(2);
      const dim_t input_width = input.dim(3);
      auto* out = output.data<float>();
      const auto* in = input.data<float>();
      dim_t out_offset = 0;
      const auto in_batch_stride = in_channels * input_height * input_width;
      for (dim_t batch_offset = 0; batch_offset < batch_size * in_batch_stride; batch_offset += in_batch_stride) {
        for (int i = -_padding_h; i <= input_height - kernel_height + _padding_h; i += _stride_h) {
          for (int j = -_padding_w; j <= input_width - kernel_width + _padding_w; j += _stride_w) {
            for (dim_t c = batch_offset; c < batch_offset + in_channels * input_height * input_width; c += input_height * input_width) {
              for (int kh = 0; kh < kernel_height; ++kh) {
                for (int kw = 0; kw < kernel_width; ++kw) {
                  auto window_i = kh + i;
                  auto window_j = kw + j;
                  if (window_i >= 0 && window_i < input_height && window_j >= 0 && window_j < input_width) {
                    out[out_offset] = in[window_i * input_width + window_j + c];
                  }
                  out_offset += 1;
                }
              }
            }
          }
        }
      }
    }

  }
}
