/*
 * Copyright (c) 2016 Robert W. Rose
 *
 * MIT License, see LICENSE file.
 */

#include "keras_model.h"

#include <cmath>
#include <istream>
#include <limits>
#include <stdio.h>
#include <utility>
#include <sstream>

#include <simde/simde-common.h>

bool ReadUnsignedInt(std::istream* file, unsigned int* i) {
    KASSERT(file, "Invalid file stream");
    KASSERT(i, "Invalid pointer");

    file->read((char*)i, sizeof(unsigned int));
#if SIMDE_ENDIAN_ORDER == SIMDE_ENDIAN_BIG
    *i = __builtin_bswap32(*i);
#endif
    KASSERT(file->gcount() == sizeof(unsigned int), "Expected unsigned int");

    return true;
}

bool ReadFloat(std::istream* file, float* f) {
    KASSERT(file, "Invalid file stream");
    KASSERT(f, "Invalid pointer");

    file->read((char*)f, sizeof(float));
#if SIMDE_ENDIAN_ORDER == SIMDE_ENDIAN_BIG
    unsigned int fi;
    memcpy(&fi, &f, sizeof(unsigned int));
    fi = __builtin_bswap32(fi);
    memcpy(f, &fi, sizeof(float));
#endif
    KASSERT(file->gcount() == sizeof(float), "Expected float");

    return true;
}

bool ReadFloats(std::istream* file, float* f, size_t n) {
    KASSERT(file, "Invalid file stream");
    KASSERT(f, "Invalid pointer");

    file->read((char*)f, sizeof(float) * n);
#if SIMDE_ENDIAN_ORDER == SIMDE_ENDIAN_BIG
    for (size_t i = 0; i < n; ++i) {
       unsigned int fi;
       memcpy(&fi, &f[i], sizeof(unsigned int));
       fi = __builtin_bswap32(fi);
       memcpy(f + i, &fi, sizeof(float));
    }
#endif
    KASSERT(((unsigned int)file->gcount()) == sizeof(float) * n,
            "Expected floats");

    return true;
}

bool KerasLayerActivation::LoadLayer(std::istream* file) {
    KASSERT(file, "Invalid file stream");

    unsigned int activation = 0;
    KASSERT(ReadUnsignedInt(file, &activation),
            "Failed to read activation type");

    switch (activation) {
    case kLinear:
        activation_type_ = kLinear;
        break;
    case kRelu:
        activation_type_ = kRelu;
        break;
    case kSoftPlus:
        activation_type_ = kSoftPlus;
        break;
    case kHardSigmoid:
        activation_type_ = kHardSigmoid;
        break;
    case kSigmoid:
        activation_type_ = kSigmoid;
        break;
    case kTanh:
        activation_type_ = kTanh;
        break;
    default:
        KASSERT(false, "Unsupported activation type %d", activation);
    }

    return true;
}

bool KerasLayerActivation::Apply(Tensor* in, Tensor* out) {
    KASSERT(in, "Invalid input");
    KASSERT(out, "Invalid output");

    *out = *in;

    switch (activation_type_) {
    case kLinear:
        break;
    case kRelu:
        for (size_t i = 0; i < out->data_.size(); i++) {
            if (out->data_[i] < 0.0) {
                out->data_[i] = 0.0;
            }
        }
        break;
    case kSoftPlus:
        for (size_t i = 0; i < out->data_.size(); i++) {
            out->data_[i] = std::log(1.0 + std::exp(out->data_[i]));
        }
        break;
    case kHardSigmoid:
        for (size_t i = 0; i < out->data_.size(); i++) {
            float x = (out->data_[i] * 0.2) + 0.5;

            if (x <= 0) {
                out->data_[i] = 0.0;
            } else if (x >= 1) {
                out->data_[i] = 1.0;
            } else {
                out->data_[i] = x;
            }
        }
        break;
    case kSigmoid:
        for (size_t i = 0; i < out->data_.size(); i++) {
            float& x = out->data_[i];

            if (x >= 0) {
                out->data_[i] = 1.0 / (1.0 + std::exp(-x));
            } else {
                float z = std::exp(x);
                out->data_[i] = z / (1.0 + z);
            }
        }
        break;
    case kTanh:
        for (size_t i = 0; i < out->data_.size(); i++) {
            out->data_[i] = std::tanh(out->data_[i]);
        }
        break;
    default:
        break;
    }

    return true;
}

bool KerasLayerDense::LoadLayer(std::istream* file) {
    KASSERT(file, "Invalid file stream");

    unsigned int weights_rows = 0;
    KASSERT(ReadUnsignedInt(file, &weights_rows), "Expected weight rows");
    KASSERT(weights_rows > 0, "Invalid weights # rows");

    unsigned int weights_cols = 0;
    KASSERT(ReadUnsignedInt(file, &weights_cols), "Expected weight cols");
    KASSERT(weights_cols > 0, "Invalid weights shape");

    unsigned int biases_shape = 0;
    KASSERT(ReadUnsignedInt(file, &biases_shape), "Expected biases shape");
    KASSERT(biases_shape > 0, "Invalid biases shape");

    weights_.Resize(weights_rows, weights_cols);
    KASSERT(
        ReadFloats(file, weights_.data_.data(), weights_rows * weights_cols),
        "Expected weights");

    biases_.Resize(biases_shape);
    KASSERT(ReadFloats(file, biases_.data_.data(), biases_shape),
            "Expected biases");

    KASSERT(activation_.LoadLayer(file), "Failed to load activation");

    return true;
}

bool KerasLayerDense::Apply(Tensor* in, Tensor* out) {
    KASSERT(in, "Invalid input");
    KASSERT(out, "Invalid output");
    KASSERT(in->dims_.size() <= 2, "Invalid input dimensions");

    if (in->dims_.size() == 2) {
        KASSERT(in->dims_[1] == weights_.dims_[0], "Dimension mismatch %d %d",
                in->dims_[1], weights_.dims_[0]);
    }

    Tensor tmp(weights_.dims_[1]);

    for (int i = 0; i < weights_.dims_[0]; i++) {
        for (int j = 0; j < weights_.dims_[1]; j++) {
            tmp(j) += (*in)(i)*weights_(i, j);
        }
    }

    for (int i = 0; i < biases_.dims_[0]; i++) {
        tmp(i) += biases_(i);
    }

    KASSERT(activation_.Apply(&tmp, out), "Failed to apply activation");

    return true;
}

bool KerasLayerConvolution2d::LoadLayer(std::istream* file) {
    KASSERT(file, "Invalid file stream");

    unsigned int weights_i = 0;
    KASSERT(ReadUnsignedInt(file, &weights_i), "Expected weights_i");
    KASSERT(weights_i > 0, "Invalid weights # i");

    unsigned int weights_j = 0;
    KASSERT(ReadUnsignedInt(file, &weights_j), "Expected weights_j");
    KASSERT(weights_j > 0, "Invalid weights # j");

    unsigned int weights_k = 0;
    KASSERT(ReadUnsignedInt(file, &weights_k), "Expected weights_k");
    KASSERT(weights_k > 0, "Invalid weights # k");

    unsigned int weights_l = 0;
    KASSERT(ReadUnsignedInt(file, &weights_l), "Expected weights_l");
    KASSERT(weights_l > 0, "Invalid weights # l");

    unsigned int biases_shape = 0;
    KASSERT(ReadUnsignedInt(file, &biases_shape), "Expected biases shape");
    KASSERT(biases_shape > 0, "Invalid biases shape");

    weights_.Resize(weights_i, weights_j, weights_k, weights_l);
    KASSERT(ReadFloats(file, weights_.data_.data(),
                       weights_i * weights_j * weights_k * weights_l),
            "Expected weights");

    biases_.Resize(biases_shape);
    KASSERT(ReadFloats(file, biases_.data_.data(), biases_shape),
            "Expected biases");

    KASSERT(activation_.LoadLayer(file), "Failed to load activation");

    return true;
}

bool KerasLayerConvolution2d::Apply(Tensor* in, Tensor* out) {
    KASSERT(in, "Invalid input");
    KASSERT(out, "Invalid output");

    KASSERT(in->dims_[0] == weights_.dims_[1],
            "Input 'depth' doesn't match kernel 'depth'");

    int st_nj = (weights_.dims_[2] - 1) / 2;
    int st_pj = (weights_.dims_[2]) / 2;
    int st_nk = (weights_.dims_[3] - 1) / 2;
    int st_pk = (weights_.dims_[3]) / 2;

    Tensor tmp(weights_.dims_[0], in->dims_[1] - st_nj - st_pj,
               in->dims_[2] - st_nk - st_pk);

    // Iterate over each kernel.
    for (int i = 0; i < weights_.dims_[0]; i++) {
        // Iterate over each 'depth'.
        for (int j = 0; j < weights_.dims_[1]; j++) {
            // 2D convolution in x and y (k and l in Tensor dimensions).
            for (int tj = st_nj; tj < in->dims_[1] - st_pj; tj++) {
                for (int tk = st_nk; tk < in->dims_[2] - st_pk; tk++) {
                    // Iterate over kernel.
                    for (int k = 0; k < weights_.dims_[2]; k++) {
                        for (int l = 0; l < weights_.dims_[3]; l++) {
                            const float& weight = weights_(i, j, k, l);
                            const float& value =
                                (*in)(j, tj - st_nj + k, tk - st_nk + l);

                            tmp(i, tj - st_nj, tk - st_nk) += weight * value;
                        }
                    }
                }
            }
        }

        // Apply kernel bias to all points in output.
        for (int j = 0; j < tmp.dims_[1]; j++) {
            for (int k = 0; k < tmp.dims_[2]; k++) {
                tmp(i, j, k) += biases_(i);
            }
        }
    }

    KASSERT(activation_.Apply(&tmp, out), "Failed to apply activation");

    return true;
}

bool KerasLayerFlatten::LoadLayer(std::istream* file) {
    KASSERT(file, "Invalid file stream");
    return true;
}

bool KerasLayerFlatten::Apply(Tensor* in, Tensor* out) {
    KASSERT(in, "Invalid input");
    KASSERT(out, "Invalid output");

    *out = *in;
    out->Flatten();

    return true;
}

bool KerasLayerElu::LoadLayer(std::istream* file) {
    KASSERT(file, "Invalid file stream");

    KASSERT(ReadFloat(file, &alpha_), "Failed to read alpha");

    return true;
}

bool KerasLayerElu::Apply(Tensor* in, Tensor* out) {
    KASSERT(in, "Invalid input");
    KASSERT(out, "Invalid output");

    *out = *in;

    for (size_t i = 0; i < out->data_.size(); i++) {
        if (out->data_[i] < 0.0) {
            out->data_[i] = alpha_ * (exp(out->data_[i]) - 1.0);
        }
    }

    return true;
}

bool KerasLayerMaxPooling2d::LoadLayer(std::istream* file) {
    KASSERT(file, "Invalid file stream");

    KASSERT(ReadUnsignedInt(file, &pool_size_j_), "Expected pool size j");
    KASSERT(ReadUnsignedInt(file, &pool_size_k_), "Expected pool size k");

    return true;
}

bool KerasLayerMaxPooling2d::Apply(Tensor* in, Tensor* out) {
    KASSERT(in, "Invalid input");
    KASSERT(out, "Invalid output");

    KASSERT(in->dims_.size() == 3, "Input must have 3 dimensions");

    Tensor tmp(in->dims_[0], in->dims_[1] / pool_size_j_,
               in->dims_[2] / pool_size_k_);

    for (int i = 0; i < tmp.dims_[0]; i++) {
        for (int j = 0; j < tmp.dims_[1]; j++) {
            const int tj = j * pool_size_j_;

            for (int k = 0; k < tmp.dims_[2]; k++) {
                const int tk = k * pool_size_k_;

                // Find maximum value over patch starting at tj, tk.
                float max_val = -std::numeric_limits<float>::infinity();

                for (unsigned int pj = 0; pj < pool_size_j_; pj++) {
                    for (unsigned int pk = 0; pk < pool_size_k_; pk++) {
                        const float& pool_val = (*in)(i, tj + pj, tk + pk);
                        if (pool_val > max_val) {
                            max_val = pool_val;
                        }
                    }
                }

                tmp(i, j, k) = max_val;
            }
        }
    }

    *out = tmp;

    return true;
}

bool KerasLayerLSTM::LoadLayer(std::istream* file) {
    KASSERT(file, "Invalid file stream");

    unsigned int wi_rows = 0;
    KASSERT(ReadUnsignedInt(file, &wi_rows), "Expected Wi rows");
    KASSERT(wi_rows > 0, "Invalid Wi # rows");

    unsigned int wi_cols = 0;
    KASSERT(ReadUnsignedInt(file, &wi_cols), "Expected Wi cols");
    KASSERT(wi_cols > 0, "Invalid Wi shape");

    unsigned int ui_rows = 0;
    KASSERT(ReadUnsignedInt(file, &ui_rows), "Expected Ui rows");
    KASSERT(ui_rows > 0, "Invalid Ui # rows");

    unsigned int ui_cols = 0;
    KASSERT(ReadUnsignedInt(file, &ui_cols), "Expected Ui cols");
    KASSERT(ui_cols > 0, "Invalid Ui shape");

    unsigned int bi_shape = 0;
    KASSERT(ReadUnsignedInt(file, &bi_shape), "Expected bi shape");
    KASSERT(bi_shape > 0, "Invalid bi shape");

    unsigned int wf_rows = 0;
    KASSERT(ReadUnsignedInt(file, &wf_rows), "Expected Wf rows");
    KASSERT(wf_rows > 0, "Invalid Wf # rows");

    unsigned int wf_cols = 0;
    KASSERT(ReadUnsignedInt(file, &wf_cols), "Expected Wf cols");
    KASSERT(wf_cols > 0, "Invalid Wf shape");

    unsigned int uf_rows = 0;
    KASSERT(ReadUnsignedInt(file, &uf_rows), "Expected Uf rows");
    KASSERT(uf_rows > 0, "Invalid Uf # rows");

    unsigned int uf_cols = 0;
    KASSERT(ReadUnsignedInt(file, &uf_cols), "Expected Uf cols");
    KASSERT(uf_cols > 0, "Invalid Uf shape");

    unsigned int bf_shape = 0;
    KASSERT(ReadUnsignedInt(file, &bf_shape), "Expected bf shape");
    KASSERT(bf_shape > 0, "Invalid bf shape");

    unsigned int wc_rows = 0;
    KASSERT(ReadUnsignedInt(file, &wc_rows), "Expected Wc rows");
    KASSERT(wc_rows > 0, "Invalid Wc # rows");

    unsigned int wc_cols = 0;
    KASSERT(ReadUnsignedInt(file, &wc_cols), "Expected Wc cols");
    KASSERT(wc_cols > 0, "Invalid Wc shape");

    unsigned int uc_rows = 0;
    KASSERT(ReadUnsignedInt(file, &uc_rows), "Expected Uc rows");
    KASSERT(uc_rows > 0, "Invalid Uc # rows");

    unsigned int uc_cols = 0;
    KASSERT(ReadUnsignedInt(file, &uc_cols), "Expected Uc cols");
    KASSERT(uc_cols > 0, "Invalid Uc shape");

    unsigned int bc_shape = 0;
    KASSERT(ReadUnsignedInt(file, &bc_shape), "Expected bc shape");
    KASSERT(bc_shape > 0, "Invalid bc shape");

    unsigned int wo_rows = 0;
    KASSERT(ReadUnsignedInt(file, &wo_rows), "Expected Wo rows");
    KASSERT(wo_rows > 0, "Invalid Wo # rows");

    unsigned int wo_cols = 0;
    KASSERT(ReadUnsignedInt(file, &wo_cols), "Expected Wo cols");
    KASSERT(wo_cols > 0, "Invalid Wo shape");

    unsigned int uo_rows = 0;
    KASSERT(ReadUnsignedInt(file, &uo_rows), "Expected Uo rows");
    KASSERT(uo_rows > 0, "Invalid Uo # rows");

    unsigned int uo_cols = 0;
    KASSERT(ReadUnsignedInt(file, &uo_cols), "Expected Uo cols");
    KASSERT(uo_cols > 0, "Invalid Uo shape");

    unsigned int bo_shape = 0;
    KASSERT(ReadUnsignedInt(file, &bo_shape), "Expected bo shape");
    KASSERT(bo_shape > 0, "Invalid bo shape");

    // Load Input Weights and Biases
    Wi_.Resize(wi_rows, wi_cols);
    KASSERT(ReadFloats(file, Wi_.data_.data(), wi_rows * wi_cols),
            "Expected Wi weights");

    Ui_.Resize(ui_rows, ui_cols);
    KASSERT(ReadFloats(file, Ui_.data_.data(), ui_rows * ui_cols),
            "Expected Ui weights");

    bi_.Resize(1, bi_shape);
    KASSERT(ReadFloats(file, bi_.data_.data(), bi_shape), "Expected bi biases");

    // Load Forget Weights and Biases
    Wf_.Resize(wf_rows, wf_cols);
    KASSERT(ReadFloats(file, Wf_.data_.data(), wf_rows * wf_cols),
            "Expected Wf weights");

    Uf_.Resize(uf_rows, uf_cols);
    KASSERT(ReadFloats(file, Uf_.data_.data(), uf_rows * uf_cols),
            "Expected Uf weights");

    bf_.Resize(1, bf_shape);
    KASSERT(ReadFloats(file, bf_.data_.data(), bf_shape), "Expected bf biases");

    // Load State Weights and Biases
    Wc_.Resize(wc_rows, wc_cols);
    KASSERT(ReadFloats(file, Wc_.data_.data(), wc_rows * wc_cols),
            "Expected Wc weights");

    Uc_.Resize(uc_rows, uc_cols);
    KASSERT(ReadFloats(file, Uc_.data_.data(), uc_rows * uc_cols),
            "Expected Uc weights");

    bc_.Resize(1, bc_shape);
    KASSERT(ReadFloats(file, bc_.data_.data(), bc_shape), "Expected bc biases");

    // Load Output Weights and Biases
    Wo_.Resize(wo_rows, wo_cols);
    KASSERT(ReadFloats(file, Wo_.data_.data(), wo_rows * wo_cols),
            "Expected Wo weights");

    Uo_.Resize(uo_rows, uo_cols);
    KASSERT(ReadFloats(file, Uo_.data_.data(), uo_rows * uo_cols),
            "Expected Uo weights");

    bo_.Resize(1, bo_shape);
    KASSERT(ReadFloats(file, bo_.data_.data(), bo_shape), "Expected bo biases");

    KASSERT(innerActivation_.LoadLayer(file),
            "Failed to load inner activation");
    KASSERT(activation_.LoadLayer(file), "Failed to load activation");

    unsigned int return_sequences = 0;
    KASSERT(ReadUnsignedInt(file, &return_sequences),
            "Expected return_sequences param");
    return_sequences_ = (bool)return_sequences;

    return true;
}

bool KerasLayerLSTM::Apply(Tensor* in, Tensor* out) {
    // Assume bo always keeps the output shape and we will always receive one
    // single sample.
    int outputDim = bo_.dims_[1];
    Tensor ht_1 = Tensor(1, outputDim);
    Tensor ct_1 = Tensor(1, outputDim);

    ht_1.Fill(0.0f);
    ct_1.Fill(0.0f);

    int steps = in->dims_[0];

    Tensor outputs, lastOutput;

    if (return_sequences_) {
        outputs.dims_ = {steps, outputDim};
        outputs.data_.reserve(steps * outputDim);
    }

    for (int s = 0; s < steps; s++) {
        Tensor x = in->Select(s);

        KASSERT(Step(&x, &lastOutput, &ht_1, &ct_1), "Failed to execute step");

        if (return_sequences_) {
            outputs.data_.insert(outputs.data_.end(), lastOutput.data_.begin(),
                                 lastOutput.data_.end());
        }
    }

    if (return_sequences_) {
        *out = outputs;
    } else {
        *out = lastOutput;
    }

    return true;
}

bool KerasLayerEmbedding::LoadLayer(std::istream* file) {
    KASSERT(file, "Invalid file stream");

    unsigned int weights_rows = 0;
    KASSERT(ReadUnsignedInt(file, &weights_rows), "Expected weight rows");
    KASSERT(weights_rows > 0, "Invalid weights # rows");

    unsigned int weights_cols = 0;
    KASSERT(ReadUnsignedInt(file, &weights_cols), "Expected weight cols");
    KASSERT(weights_cols > 0, "Invalid weights shape");

    weights_.Resize(weights_rows, weights_cols);
    KASSERT(
        ReadFloats(file, weights_.data_.data(), weights_rows * weights_cols),
        "Expected weights");

    return true;
}

bool KerasLayerEmbedding::Apply(Tensor* in, Tensor* out) {
    int output_rows = in->dims_[1];
    int output_cols = weights_.dims_[1];
    out->dims_ = {output_rows, output_cols};
    out->data_.reserve(output_rows * output_cols);

    std::for_each(in->data_.begin(), in->data_.end(), [=](float i) {
        std::vector<float>::const_iterator first =
            this->weights_.data_.begin() + (i * output_cols);
        std::vector<float>::const_iterator last =
            this->weights_.data_.begin() + (i + 1) * output_cols;

        out->data_.insert(out->data_.end(), first, last);
    });

    return true;
}

bool KerasLayerLSTM::Step(Tensor* x, Tensor* out, Tensor* ht_1, Tensor* ct_1) {
    Tensor xi = x->Dot(Wi_) + bi_;
    Tensor xf = x->Dot(Wf_) + bf_;
    Tensor xc = x->Dot(Wc_) + bc_;
    Tensor xo = x->Dot(Wo_) + bo_;

    Tensor i_ = xi + ht_1->Dot(Ui_);
    Tensor f_ = xf + ht_1->Dot(Uf_);
    Tensor c_ = xc + ht_1->Dot(Uc_);
    Tensor o_ = xo + ht_1->Dot(Uo_);

    Tensor i, f, cc, o;

    KASSERT(innerActivation_.Apply(&i_, &i),
            "Failed to apply inner activation on i");
    KASSERT(innerActivation_.Apply(&f_, &f),
            "Failed to apply inner activation on f");
    KASSERT(activation_.Apply(&c_, &cc), "Failed to apply activation on c_");
    KASSERT(innerActivation_.Apply(&o_, &o),
            "Failed to apply inner activation on o");

    *ct_1 = f.Multiply(*ct_1) + i.Multiply(cc);

    KASSERT(activation_.Apply(ct_1, &cc), "Failed to apply activation on c");
    *out = *ht_1 = o.Multiply(cc);

    return true;
}

bool KerasModel::LoadModel(const std::string& data) {
    std::stringstream file(data);

    unsigned int num_layers = 0;
    KASSERT(ReadUnsignedInt(&file, &num_layers), "Expected number of layers");

    for (unsigned int i = 0; i < num_layers; i++) {
        unsigned int layer_type = 0;
        KASSERT(ReadUnsignedInt(&file, &layer_type), "Expected layer type");

        KerasLayer* layer = NULL;

        switch (layer_type) {
        case kDense:
            layer = new KerasLayerDense();
            break;
        case kConvolution2d:
            layer = new KerasLayerConvolution2d();
            break;
        case kFlatten:
            layer = new KerasLayerFlatten();
            break;
        case kElu:
            layer = new KerasLayerElu();
            break;
        case kActivation:
            layer = new KerasLayerActivation();
            break;
        case kMaxPooling2D:
            layer = new KerasLayerMaxPooling2d();
            break;
        case kLSTM:
            layer = new KerasLayerLSTM();
            break;
        case kEmbedding:
            layer = new KerasLayerEmbedding();
            break;
        default:
            break;
        }

        KASSERT(layer, "Unknown layer type %d", layer_type);

        bool result = layer->LoadLayer(&file);
        if (!result) {
            printf("Failed to load layer %d", i);
            delete layer;
            return false;
        }

        layers_.push_back(layer);
    }

    return true;
}

bool KerasModel::Apply(Tensor* in, Tensor* out) {
    Tensor temp_in, temp_out;

    for (unsigned int i = 0; i < layers_.size(); i++) {
        if (i == 0) {
            temp_in = *in;
        }

        KASSERT(layers_[i]->Apply(&temp_in, &temp_out),
                "Failed to apply layer %d", i);

        temp_in = temp_out;
    }

    *out = temp_out;

    return true;
}
