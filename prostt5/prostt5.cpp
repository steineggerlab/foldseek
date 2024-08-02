#include <vector>
#include <cstring>
#include <cstdint>
#include <string>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <ctranslate2/types.h>
#include <ctranslate2/random.h>
#include <ctranslate2/devices.h>
#include <ctranslate2/encoder.h>
#include <ctranslate2/layers/prostt5cnn.h>

#include "prostt5.h"

std::vector<std::vector<uint8_t>> compute_argmax_per_batch(const ctranslate2::StorageView& tensor) {
  if (tensor.rank() != 3) {
    throw std::runtime_error("Input tensor must have 3 dimensions (batch_size, channels, seq_length).");
  }

  const ctranslate2::dim_t batch_size = tensor.dim(0);
  const ctranslate2::dim_t channels = tensor.dim(1);
  const ctranslate2::dim_t seq_length = tensor.dim(2);

  std::vector<float> flattened_tensor = tensor.to_vector<float>();
  std::vector<std::vector<uint8_t>> argmax_indices(batch_size, std::vector<uint8_t>(seq_length));
  for (ctranslate2::dim_t b = 0; b < batch_size; ++b) {
    for (ctranslate2::dim_t s = 0; s < seq_length; ++s) {
      float max_value = -std::numeric_limits<float>::infinity();
      uint8_t max_index = 0;
      for (ctranslate2::dim_t c = 0; c < channels; ++c) {
        float value = flattened_tensor[b * channels * seq_length + c * seq_length + s];
        if (value > max_value) {
          max_value = value;
          max_index = c;
        }
      }
      argmax_indices[b][s] = max_index;
    }
  }
  return argmax_indices;
}

char number_to_char(uint32_t n) {
    switch (n) {
        case 0:  return 'A';
        case 1:  return 'C';
        case 2:  return 'D';
        case 3:  return 'E';
        case 4:  return 'F';
        case 5:  return 'G';
        case 6:  return 'H';
        case 7:  return 'I';
        case 8:  return 'K';
        case 9:  return 'L';
        case 10: return 'M';
        case 11: return 'N';
        case 12: return 'P';
        case 13: return 'Q';
        case 14: return 'R';
        case 15: return 'S';
        case 16: return 'T';
        case 17: return 'V';
        case 18: return 'W';
        case 19: return 'Y';
        default: return 'X';
    }
}

struct ProstT5::Impl {
    ctranslate2::Device device;
    ctranslate2::ComputeType compute_type;
    std::unique_ptr<ctranslate2::Encoder> enc;
    const ctranslate2::models::LanguageModel* lm;

    Impl(
        const std::string& model_dir,
        const std::string& device_str,
        const std::vector<int>& device_indices,
        size_t inter_threads,
        size_t intra_threads,
        long max_queued_batches,
        int cpu_core_offset,
        const std::string& cpu_compute_type,
        const std::string& cuda_compute_type,
        unsigned int seed
    ) {
        if (model_dir.empty()) {
            throw std::invalid_argument("Model directory is required to run translation");
        }
        if (seed != 0) {
            ctranslate2::set_random_seed(seed);
        }

        device = ctranslate2::str_to_device(device_str);
        if (device == ctranslate2::Device::CPU) {
            compute_type = ctranslate2::str_to_compute_type(cpu_compute_type);
        } else if (device == ctranslate2::Device::CUDA) {
            compute_type = ctranslate2::str_to_compute_type(cuda_compute_type);
        }

        ctranslate2::ReplicaPoolConfig pool_config;
        pool_config.num_threads_per_replica = intra_threads;
        pool_config.max_queued_batches = max_queued_batches;
        pool_config.cpu_core_offset = cpu_core_offset;

        ctranslate2::models::ModelLoader model_loader(model_dir);
        model_loader.device = device;
        model_loader.device_indices = device_indices;
        model_loader.compute_type = compute_type;
        model_loader.num_replicas_per_device = inter_threads;

        enc = std::make_unique<ctranslate2::Encoder>(model_loader, pool_config);
        const auto* model = enc->get_first_replica().model().get();
        lm = reinterpret_cast<const ctranslate2::models::LanguageModel*>(model);
    }
};

ProstT5::ProstT5(
    const std::string& model_dir,
    const std::string& device_str,
    const std::vector<int>& device_indices,
    size_t inter_threads,
    size_t intra_threads,
    long max_queued_batches,
    int cpu_core_offset,
    const std::string& cpu_compute_type,
    const std::string& cuda_compute_type,
    unsigned int seed
) : pimpl(std::make_unique<Impl>(
    model_dir,
    device_str,
    device_indices,
    inter_threads,
    intra_threads,
    max_queued_batches,
    cpu_core_offset,
    cpu_compute_type,
    cuda_compute_type,
    seed)) {}

ProstT5::~ProstT5() = default;

std::vector<std::string> ProstT5::predict(const std::vector<std::string>& input_seqs) {
    std::vector<std::string> predictions;
    const ctranslate2::Vocabulary& vocab = pimpl->lm->get_vocabulary();
    std::vector<std::vector<size_t>> input;
    for (const auto& input_seq : input_seqs) {
        std::vector<size_t> seq;
        seq.emplace_back(149);
        for (const char& c : input_seq) {
            std::string curr(1, c);
            size_t id = vocab.to_id(curr, true);
            seq.emplace_back(id);
        }
        seq.emplace_back(1);
        input.emplace_back(seq);
    }

    auto forward = pimpl->enc->forward_batch_async(input);
    auto ys = forward.get().last_hidden_state;

    std::vector<std::vector<uint8_t>> argmax_results = compute_argmax_per_batch(ys);
    for (size_t i = 0; i < argmax_results.size(); ++i) {
        std::string pred;
        size_t seq_length = input_seqs[i].length();
        for (size_t j = 0; j < argmax_results[i].size() && j < seq_length; ++j) {
            pred.append(1, number_to_char(argmax_results[i][j]));
        }
        predictions.push_back(pred);
    }
    return predictions;
}
