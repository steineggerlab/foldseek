/*
 * Copyright (c) 2016 Robert W. Rose
 *
 * MIT License, see LICENSE file.
 */

#ifndef KERAS_MODEL_H_
#define KERAS_MODEL_H_

#include <algorithm>
#include <chrono>
#include <math.h>
#include <numeric>
#include <string>
#include <vector>

#define KASSERT(x, ...)                                                        \
    if (!(x)) {                                                                \
        printf("KASSERT: %s(%d): ", __FILE__, __LINE__);                       \
        printf(__VA_ARGS__);                                                   \
        printf("\n");                                                          \
        return false;                                                          \
    }

#define KASSERT_EQ(x, y, eps)                                                  \
    if (fabs(x - y) > eps) {                                                   \
        printf("KASSERT: Expected %f, got %f\n", y, x);                        \
        return false;                                                          \
    }

#ifdef DEBUG
#define KDEBUG(x, ...)                                                         \
    if (!(x)) {                                                                \
        printf("%s(%d): ", __FILE__, __LINE__);                                \
        printf(__VA_ARGS__);                                                   \
        printf("\n");                                                          \
        exit(-1);                                                              \
    }
#else
#define KDEBUG(x, ...) ;
#endif

class Tensor {
  public:
    Tensor() {}

    Tensor(int i) { Resize(i); }

    Tensor(int i, int j) { Resize(i, j); }

    Tensor(int i, int j, int k) { Resize(i, j, k); }

    Tensor(int i, int j, int k, int l) { Resize(i, j, k, l); }

    void Resize(int i) {
        dims_ = {i};
        data_.resize(i);
    }

    void Resize(int i, int j) {
        dims_ = {i, j};
        data_.resize(i * j);
    }

    void Resize(int i, int j, int k) {
        dims_ = {i, j, k};
        data_.resize(i * j * k);
    }

    void Resize(int i, int j, int k, int l) {
        dims_ = {i, j, k, l};
        data_.resize(i * j * k * l);
    }

    inline void Flatten() {
        KDEBUG(dims_.size() > 0, "Invalid tensor");

        int elements = dims_[0];
        for (unsigned int i = 1; i < dims_.size(); i++) {
            elements *= dims_[i];
        }
        dims_ = {elements};
    }

    inline float& operator()(int i) {
        KDEBUG(dims_.size() == 1, "Invalid indexing for tensor");
        KDEBUG(i < dims_[0] && i >= 0, "Invalid i: %d (max %d)", i, dims_[0]);

        return data_[i];
    }

    inline float& operator()(int i, int j) {
        KDEBUG(dims_.size() == 2, "Invalid indexing for tensor");
        KDEBUG(i < dims_[0] && i >= 0, "Invalid i: %d (max %d)", i, dims_[0]);
        KDEBUG(j < dims_[1] && j >= 0, "Invalid j: %d (max %d)", j, dims_[1]);

        return data_[dims_[1] * i + j];
    }

    inline float operator()(int i, int j) const {
        KDEBUG(dims_.size() == 2, "Invalid indexing for tensor");
        KDEBUG(i < dims_[0] && i >= 0, "Invalid i: %d (max %d)", i, dims_[0]);
        KDEBUG(j < dims_[1] && j >= 0, "Invalid j: %d (max %d)", j, dims_[1]);

        return data_[dims_[1] * i + j];
    }

    inline float& operator()(int i, int j, int k) {
        KDEBUG(dims_.size() == 3, "Invalid indexing for tensor");
        KDEBUG(i < dims_[0] && i >= 0, "Invalid i: %d (max %d)", i, dims_[0]);
        KDEBUG(j < dims_[1] && j >= 0, "Invalid j: %d (max %d)", j, dims_[1]);
        KDEBUG(k < dims_[2] && k >= 0, "Invalid k: %d (max %d)", k, dims_[2]);

        return data_[dims_[2] * (dims_[1] * i + j) + k];
    }

    inline float& operator()(int i, int j, int k, int l) {
        KDEBUG(dims_.size() == 4, "Invalid indexing for tensor");
        KDEBUG(i < dims_[0] && i >= 0, "Invalid i: %d (max %d)", i, dims_[0]);
        KDEBUG(j < dims_[1] && j >= 0, "Invalid j: %d (max %d)", j, dims_[1]);
        KDEBUG(k < dims_[2] && k >= 0, "Invalid k: %d (max %d)", k, dims_[2]);
        KDEBUG(l < dims_[3] && l >= 0, "Invalid l: %d (max %d)", l, dims_[3]);

        return data_[dims_[3] * (dims_[2] * (dims_[1] * i + j) + k) + l];
    }

    inline void Fill(float value) {
        std::fill(data_.begin(), data_.end(), value);
    }

    Tensor Unpack(int row) const {
        KASSERT(dims_.size() >= 2, "Invalid tensor");
        std::vector<int> pack_dims =
            std::vector<int>(dims_.begin() + 1, dims_.end());
        int pack_size = std::accumulate(pack_dims.begin(), pack_dims.end(), 0);

        std::vector<float>::const_iterator first =
            data_.begin() + (row * pack_size);
        std::vector<float>::const_iterator last =
            data_.begin() + (row + 1) * pack_size;

        Tensor x = Tensor();
        x.dims_ = pack_dims;
        x.data_ = std::vector<float>(first, last);

        return x;
    }

    Tensor Select(int row) const {
        Tensor x = Unpack(row);
        x.dims_.insert(x.dims_.begin(), 1);

        return x;
    }

    Tensor operator+(const Tensor& other) {
        KASSERT(dims_ == other.dims_,
                "Cannot add tensors with different dimensions");

        Tensor result;
        result.dims_ = dims_;
        result.data_.reserve(data_.size());

        std::transform(data_.begin(), data_.end(), other.data_.begin(),
                       std::back_inserter(result.data_),
                       [](float x, float y) { return x + y; });

        return result;
    }

    Tensor Multiply(const Tensor& other) {
        KASSERT(dims_ == other.dims_,
                "Cannot multiply elements with different dimensions");

        Tensor result;
        result.dims_ = dims_;
        result.data_.reserve(data_.size());

        std::transform(data_.begin(), data_.end(), other.data_.begin(),
                       std::back_inserter(result.data_),
                       [](float x, float y) { return x * y; });

        return result;
    }

    Tensor Dot(const Tensor& other) {
        KDEBUG(dims_.size() == 2, "Invalid tensor dimensions");
        KDEBUG(other.dims_.size() == 2, "Invalid tensor dimensions");
        KASSERT(dims_[1] == other.dims_[0],
                "Cannot multiply with different inner dimensions");

        Tensor tmp(dims_[0], other.dims_[1]);

        for (int i = 0; i < dims_[0]; i++) {
            for (int j = 0; j < other.dims_[1]; j++) {
                for (int k = 0; k < dims_[1]; k++) {
                    tmp(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }

        return tmp;
    }

    void Print() {
        if (dims_.size() == 1) {
            printf("[ ");
            for (int i = 0; i < dims_[0]; i++) {
                printf("%f ", (*this)(i));
            }
            printf("]\n");
        } else if (dims_.size() == 2) {
            printf("[\n");
            for (int i = 0; i < dims_[0]; i++) {
                printf(" [ ");
                for (int j = 0; j < dims_[1]; j++) {
                    printf("%f ", (*this)(i, j));
                }
                printf("]\n");
            }
            printf("]\n");
        } else if (dims_.size() == 3) {
            printf("[\n");
            for (int i = 0; i < dims_[0]; i++) {
                printf(" [\n");
                for (int j = 0; j < dims_[1]; j++) {
                    printf("  [ ");
                    for (int k = 0; k < dims_[2]; k++) {
                        printf("%f ", (*this)(i, j, k));
                    }
                    printf("  ]\n");
                }
                printf(" ]\n");
            }
            printf("]\n");
        } else if (dims_.size() == 4) {
            printf("[\n");
            for (int i = 0; i < dims_[0]; i++) {
                printf(" [\n");
                for (int j = 0; j < dims_[1]; j++) {
                    printf("  [\n");
                    for (int k = 0; k < dims_[2]; k++) {
                        printf("   [");
                        for (int l = 0; l < dims_[3]; l++) {
                            printf("%f ", (*this)(i, j, k, l));
                        }
                        printf("]\n");
                    }
                    printf("  ]\n");
                }
                printf(" ]\n");
            }
            printf("]\n");
        }
    }

    void PrintShape() {
        printf("(");
        for (unsigned int i = 0; i < dims_.size(); i++) {
            printf("%d ", dims_[i]);
        }
        printf(")\n");
    }

    std::vector<int> dims_;
    std::vector<float> data_;
};

class KerasLayer {
  public:
    KerasLayer() {}

    virtual ~KerasLayer() {}

    virtual bool LoadLayer(std::istream* file) = 0;

    virtual bool Apply(Tensor* in, Tensor* out) = 0;
};

class KerasLayerActivation : public KerasLayer {
  public:
    enum ActivationType {
        kLinear = 1,
        kRelu = 2,
        kSoftPlus = 3,
        kSigmoid = 4,
        kTanh = 5,
        kHardSigmoid = 6
    };

    KerasLayerActivation() : activation_type_(ActivationType::kLinear) {}

    virtual ~KerasLayerActivation() {}

    virtual bool LoadLayer(std::istream* file);

    virtual bool Apply(Tensor* in, Tensor* out);

  private:
    ActivationType activation_type_;
};

class KerasLayerDense : public KerasLayer {
  public:
    KerasLayerDense() {}

    virtual ~KerasLayerDense() {}

    virtual bool LoadLayer(std::istream* file);

    virtual bool Apply(Tensor* in, Tensor* out);

  private:
    Tensor weights_;
    Tensor biases_;

    KerasLayerActivation activation_;
};

class KerasLayerConvolution2d : public KerasLayer {
  public:
    KerasLayerConvolution2d() {}

    virtual ~KerasLayerConvolution2d() {}

    virtual bool LoadLayer(std::istream* file);

    virtual bool Apply(Tensor* in, Tensor* out);

  private:
    Tensor weights_;
    Tensor biases_;

    KerasLayerActivation activation_;
};

class KerasLayerFlatten : public KerasLayer {
  public:
    KerasLayerFlatten() {}

    virtual ~KerasLayerFlatten() {}

    virtual bool LoadLayer(std::istream* file);

    virtual bool Apply(Tensor* in, Tensor* out);

  private:
};

class KerasLayerElu : public KerasLayer {
  public:
    KerasLayerElu() : alpha_(1.0f) {}

    virtual ~KerasLayerElu() {}

    virtual bool LoadLayer(std::istream* file);

    virtual bool Apply(Tensor* in, Tensor* out);

  private:
    float alpha_;
};

class KerasLayerMaxPooling2d : public KerasLayer {
  public:
    KerasLayerMaxPooling2d() : pool_size_j_(0), pool_size_k_(0) {}

    virtual ~KerasLayerMaxPooling2d() {}

    virtual bool LoadLayer(std::istream* file);

    virtual bool Apply(Tensor* in, Tensor* out);

  private:
    unsigned int pool_size_j_;
    unsigned int pool_size_k_;
};

class KerasLayerLSTM : public KerasLayer {
  public:
    KerasLayerLSTM() : return_sequences_(false) {}

    virtual ~KerasLayerLSTM() {}

    virtual bool LoadLayer(std::istream* file);

    virtual bool Apply(Tensor* in, Tensor* out);

  private:
    bool Step(Tensor* x, Tensor* out, Tensor* ht_1, Tensor* ct_1);

    Tensor Wi_;
    Tensor Ui_;
    Tensor bi_;
    Tensor Wf_;
    Tensor Uf_;
    Tensor bf_;
    Tensor Wc_;
    Tensor Uc_;
    Tensor bc_;
    Tensor Wo_;
    Tensor Uo_;
    Tensor bo_;

    KerasLayerActivation innerActivation_;
    KerasLayerActivation activation_;
    bool return_sequences_;
};

class KerasLayerEmbedding : public KerasLayer {
  public:
    KerasLayerEmbedding() {}

    virtual ~KerasLayerEmbedding() {}

    virtual bool LoadLayer(std::istream* file);

    virtual bool Apply(Tensor* in, Tensor* out);

  private:
    Tensor weights_;
};

class KerasModel {
  public:
    enum LayerType {
        kDense = 1,
        kConvolution2d = 2,
        kFlatten = 3,
        kElu = 4,
        kActivation = 5,
        kMaxPooling2D = 6,
        kLSTM = 7,
        kEmbedding = 8
    };

    KerasModel() {}

    virtual ~KerasModel() {
        for (unsigned int i = 0; i < layers_.size(); i++) {
            delete layers_[i];
        }
    }

    virtual bool LoadModel(const std::string& filename);

    virtual bool Apply(Tensor* in, Tensor* out);

  private:
    std::vector<KerasLayer*> layers_;
};

class KerasTimer {
  public:
    KerasTimer() {}

    void Start() { start_ = std::chrono::high_resolution_clock::now(); }

    double Stop() {
        std::chrono::time_point<std::chrono::high_resolution_clock> now =
            std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> diff = now - start_;

        return diff.count();
    }

  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

#endif // KERAS_MODEL_H_
