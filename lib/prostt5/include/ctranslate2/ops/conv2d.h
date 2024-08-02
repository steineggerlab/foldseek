#pragma once

#include "activation.h"
#include "op.h"

namespace ctranslate2 {
  namespace ops {

    class Conv2D : public Op {
    public:
      Conv2D(dim_t stride_h = 1, dim_t stride_w = 1, dim_t padding_h = 0, dim_t padding_w = 0, dim_t dilation = 1);

      void operator()(const StorageView& input,
                      const StorageView& weight,
                      const StorageView& bias,
                      StorageView& output,
                      const StorageView* qscale = nullptr) const;

      void operator()(const StorageView& input,
                      const StorageView& weight,
                      StorageView& output,
                      const StorageView* qscale = nullptr) const;

    private:
      dim_t _stride_h;
      dim_t _stride_w;
      dim_t _padding_h;
      dim_t _padding_w;
      dim_t _dilation;

      void operator()(const StorageView& input,
                      const StorageView& weight,
                      const StorageView* bias,
                      StorageView& output,
                      const StorageView* qscale) const;

      template <Device D, typename T>
      void compute(const StorageView& input,
                   const StorageView& weight,
                   const StorageView* bias,
                   StorageView& output,
                   const StorageView* qscale = nullptr) const;

      void compute_with_gemm(const StorageView& input, const StorageView& weight, StorageView& output,
                             const StorageView* qscale) const;

      void im2col_transposed(const StorageView& input, StorageView& output, dim_t kernel_height, dim_t kernel_width) const;
#ifdef CT2_WITH_CUDA
      void im2col_transposed_gpu(const StorageView& input, StorageView& output, dim_t kernel_height, dim_t kernel_width) const;
#endif
    };

  }
}
