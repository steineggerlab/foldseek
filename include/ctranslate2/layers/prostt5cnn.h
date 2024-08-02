#pragma once

#include "ctranslate2/layers/common.h"
#include "ctranslate2/ops/conv2d.h"
#include "ctranslate2/ops/relu.h"
#include "ctranslate2/ops/slide.h"
#include "ctranslate2/ops/concat.h"
#include "ctranslate2/ops/transpose.h"

namespace ctranslate2 {
  namespace layers {
    class ProstT5CNN : public Layer {
    public:
      ProstT5CNN(const models::Model& model) :
        ProstT5CNN(
          model.get_variable("cnv1/weight"),
          model.get_variable("cnv1/bias"),
          model.get_variable("cnv2/weight"),
          model.get_variable("cnv2/bias")
        ) {}
      ProstT5CNN(const StorageView& conv0, const StorageView& bias0,
                 const StorageView& conv1, const StorageView& bias1);
      void operator()(const StorageView& input, StorageView& output) const;
      DataType output_type() const override;
      dim_t output_size() const override;

    private:
      ops::Concat _concat_op;
      ops::Transpose _transpose_op;
      ops::Conv2D _conv_op;
      ops::ReLU _relu_op;
      const StorageView& _conv0;
      const StorageView& _bias0;
      const StorageView& _conv1;
      const StorageView& _bias1;
    };

  }
}
