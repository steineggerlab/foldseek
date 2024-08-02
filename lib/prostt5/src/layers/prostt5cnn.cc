#include "ctranslate2/layers/prostt5cnn.h"

namespace ctranslate2 {
  namespace layers {

    ProstT5CNN::ProstT5CNN(const StorageView& conv0, const StorageView& bias0,
                           const StorageView& conv1, const StorageView& bias1)
      : _concat_op(1),
        _transpose_op(std::vector<dim_t>{0, 2, 1}),
        _conv_op(1, 1, 3, 0, 1), // kernel_size = 1x1, stride = 1, padding = 3x0, dilation = 1
        _relu_op(),
        _conv0(conv0),
        _bias0(bias0),
        _conv1(conv1),
        _bias1(bias1) {}

    void ProstT5CNN::operator()(const StorageView& input, StorageView& output) const {
      const DataType dtype = input.dtype();
      const Device device = input.device();
      const dim_t batch_size = input.dim(0);
      const dim_t seq_len = input.dim(1);
      const dim_t hidden_dim = input.dim(2);

      // Adjust _slide_op size based on seq_len
      ops::Slide _slide_op(1, 1, seq_len - 2, false);

      // Step 1: Slide operation
      StorageView sliced_ys(device, dtype);
      _slide_op(input, sliced_ys);

      // Step 2: Padding
      StorageView pad_right(std::vector<dim_t>{1, 1, hidden_dim}, dtype, device);
      pad_right.zero();

      StorageView padded(device, dtype);
      _concat_op({&sliced_ys, &pad_right}, padded);

      // Step 3: Transpose and expand dims
      StorageView transposed(device, dtype);
      _transpose_op(padded, transposed);
      transposed.expand_dims(3);

      // Step 4: First convolution
      StorageView conv0_out(device, dtype);
      _conv_op(transposed, _conv0, _bias0, conv0_out);

      // Step 5: ReLU activation
      StorageView relu0_out(device, dtype);
      _relu_op(conv0_out, relu0_out);

      // Step 6: Second convolution
      StorageView conv1_out(device, dtype);
      _conv_op(relu0_out, _conv1, _bias1, conv1_out);

      // Step 7: Squeeze the last dimension
      conv1_out.squeeze(-1);

      // Set output
      output = std::move(conv1_out);
    }

    DataType ProstT5CNN::output_type() const {
      return _conv0.dtype();
    }

    dim_t ProstT5CNN::output_size() const {
      return _conv1.dim(0);
    }

   }
}