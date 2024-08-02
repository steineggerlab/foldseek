#include "ctranslate2/ops/conv2d.h"

#include "dispatch.h"

namespace ctranslate2 {
  namespace ops {

    Conv2D::Conv2D(dim_t stride_h, dim_t stride_w, dim_t padding_h, dim_t padding_w, dim_t dilation)
      : _stride_h(stride_h)
      , _stride_w(stride_w)
      , _padding_h(padding_h)
      , _padding_w(padding_w)
      , _dilation(dilation)
    {
    }

    void Conv2D::operator()(const StorageView& input,
                            const StorageView& weight,
                            const StorageView& bias,
                            StorageView& output,
                            const StorageView* qscale) const {
      operator()(input, weight, &bias, output, qscale);
    }

    void Conv2D::operator()(const StorageView& input,
                            const StorageView& weight,
                            StorageView& output,
                            const StorageView* qscale) const {
      operator()(input, weight, nullptr, output, qscale);
    }

    void Conv2D::operator()(const StorageView& input,
                            const StorageView& weight,
                            const StorageView* bias,
                            StorageView& output,
                            const StorageView* qscale) const {
      PROFILE("Conv2D");
      const dim_t batch_size = input.dim(0);
      const dim_t input_height = input.dim(2);
      const dim_t input_width = input.dim(3);
      const dim_t out_channels = weight.dim(0);
      const dim_t kernel_height = weight.dim(2);
      const dim_t kernel_width = weight.dim(3);
      const dim_t output_height = 
        // (input_height + (2 * _padding_h) - (_dilation * (kernel_height - 1) + 1)) / _stride_h + 1;
        (input_height + 2 * _padding_h - _dilation * (kernel_height - 1) - 1) / _stride_h + 1;
      const dim_t output_width = 
        // (input_width + (2 * _padding_w) - (_dilation * (kernel_width - 1) + 1)) / _stride_w + 1;
        (input_width + 2 * _padding_w - _dilation * (kernel_width - 1) - 1) / _stride_w + 1;

      output.resize({batch_size, out_channels, output_height, output_width});

      DEVICE_AND_FLOAT_DISPATCH("Conv2D", input.device(), input.dtype(),
                                (compute<D, T>(input, weight, bias, output, qscale)));
    }

  }
}
