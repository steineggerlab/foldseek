#include "ProstT5.h"

#include "llama.h"

#include <limits>
#include <vector>

static char number_to_char(unsigned int n) {
    switch(n) {
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
        default: return 'X'; // Default case for numbers not in the list
    }
}

static int encode(llama_context * ctx, std::vector<llama_token> & enc_input, std::string & result) {
    const struct llama_model * model = llama_get_model(ctx);

    if (llama_encode(ctx, llama_batch_get_one(enc_input.data(), enc_input.size())) < 0) {
        // LOG_ERR("%s : failed to encode\n", __func__);
        return 1;
    }

    // LOG_INF("%s: n_tokens = %zu, n_seq = %d\n", __func__, enc_input.size(), 1);
    float* embeddings = llama_get_embeddings(ctx);
    if (embeddings == nullptr) {
        return 1;
    }
    int * arg_max_idx = new int[enc_input.size()];
    float * arg_max = new float[enc_input.size()];
    std::fill(arg_max, arg_max + enc_input.size(), std::numeric_limits<float>::lowest());
    int seq_len = enc_input.size() - 1;
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < seq_len; ++j) {
            if(embeddings[i*seq_len + j] > arg_max[j]){
               arg_max_idx[j] = i;
               arg_max[j] = embeddings[i*seq_len + j];
            }
        }
    }
    for (int i = 0; i < seq_len - 1; ++i) {
        result.push_back(number_to_char(arg_max_idx[i]));
    }
    delete[] arg_max_idx;
    delete[] arg_max;
    return 0;
}

static std::vector<std::string> string_split(const std::string & input, char separator) {
    std::vector<std::string> parts;
    size_t begin_pos = 0;
    size_t separator_pos = input.find(separator);
    while (separator_pos != std::string::npos) {
        std::string part = input.substr(begin_pos, separator_pos - begin_pos);
        parts.emplace_back(part);
        begin_pos = separator_pos + 1;
        separator_pos = input.find(separator, begin_pos);
    }
    parts.emplace_back(input.substr(begin_pos, separator_pos - begin_pos));
    return parts;
}

static std::vector<ggml_backend_dev_t> parse_device_list(const std::string & value) {
    std::vector<ggml_backend_dev_t> devices;
    auto dev_names = string_split(value, ',');
    if (dev_names.empty()) {
        devices.push_back(nullptr);
        return devices;
    }
    if (dev_names.size() == 1 && dev_names[0] == "none") {
        devices.push_back(nullptr);
    } else {
        for (const auto & device : dev_names) {
            auto * dev = ggml_backend_dev_by_name(device.c_str());
            if (!dev || ggml_backend_dev_type(dev) != GGML_BACKEND_DEVICE_TYPE_GPU) {
                continue;
            }
            devices.push_back(dev);
        }
        devices.push_back(nullptr);
    }
    return devices;
}

// struct lora_adapter_info {
//     std::string path;
//     float scale;
// };

// struct lora_adapter_container : lora_adapter_info {
//     struct llama_lora_adapter* adapter;
// };

// struct init_result {
//     std::vector<lora_adapter_container> lora_adapters;
// };

LlamaInitGuard::LlamaInitGuard(bool verbose) {
    if (!verbose) {
        llama_log_set([](ggml_log_level, const char *, void *) {}, nullptr);
    }
    llama_backend_init();
    llama_numa_init(GGML_NUMA_STRATEGY_DISABLED);
}

LlamaInitGuard::~LlamaInitGuard() {
    llama_backend_free();
}

ProstT5Model::ProstT5Model(const std::string& model_file, std::string& device) {
    auto mparams = llama_model_default_params();
    std::vector<ggml_backend_dev_t> devices = parse_device_list(device);
    if (!devices.empty()) {
        mparams.devices = devices.data();
    }

    int gpus = 0;
    for (const auto& dev : devices) {
        if (!dev) {
            continue;
        }
        gpus += ggml_backend_dev_type(dev) == GGML_BACKEND_DEVICE_TYPE_GPU;
    }
    if (gpus > 0) {
        mparams.n_gpu_layers = 24;
    } else {
        mparams.n_gpu_layers = 0;
    }
    mparams.use_mmap        = true;
    model = llama_load_model_from_file(model_file.c_str(), mparams);

    // for (auto & la : params.lora_adapters) {
    //     lora_adapter_container loaded_la;
    //     loaded_la.path = la.path;
    //     loaded_la.scale = la.scale;
    //     loaded_la.adapter = llama_lora_adapter_init(model, la.path.c_str());
    //     if (loaded_la.adapter == nullptr) {
    //         llama_free_model(model);
    //         return;
    //     }
    //     lora_adapters.push_back(loaded_la); // copy to list of loaded adapters
    // }
}

ProstT5Model::~ProstT5Model() {
    llama_free_model(model);
}

ProstT5::ProstT5(ProstT5Model& model, int threads) : model(model) {
    auto cparams = llama_context_default_params();
    cparams.n_threads = threads;
    cparams.n_threads_batch = threads;
    cparams.n_ubatch = 2048;
    cparams.n_batch = 2048;
    cparams.n_ctx = 2048;
    cparams.embeddings = true;
    cparams.attention_type = LLAMA_ATTENTION_TYPE_NON_CAUSAL;

    ctx = llama_new_context_with_model(model.model, cparams);
    // batch = llama_batch_init(4096, 0, 1);
    // if (!params.lora_init_without_apply) {
    //     llama_lora_adapter_clear(lctx);
    //     for (auto & la : iparams.lora_adapters) {
    //         if (la.scale != 0.0f) {
    //             llama_lora_adapter_set(lctx, la.adapter, la.scale);
    //         }
    //     }
    // }
};

ProstT5::~ProstT5() {
    llama_free(ctx);
}

std::string ProstT5::predict(const std::string& aa) {
    std::string result;
    std::vector<llama_token> embd_inp;
    embd_inp.reserve(aa.length() + 2);
    embd_inp.emplace_back(llama_token_get_token(model.model, "<AA2fold>"));
    llama_token unk_aa = llama_token_get_token(model.model, "▁X");
    for (size_t i = 0; i < aa.length(); ++i) {
        std::string current_char("▁");
        current_char.append(1, toupper(aa[i]));
        llama_token token = llama_token_get_token(model.model, current_char.c_str());
        if (token == LLAMA_TOKEN_NULL) {
            embd_inp.emplace_back(unk_aa);
        } else {
            embd_inp.emplace_back(token);
        }
    }
    embd_inp.emplace_back(llama_token_get_token(model.model, "</s>"));
    encode(ctx, embd_inp, result);
    return result;
}

std::vector<std::string> ProstT5::getDevices() {
    std::vector<std::string> devices;
    for (size_t i = 0; i < ggml_backend_dev_count(); ++i) {
        ggml_backend_dev_t dev = ggml_backend_dev_get(i);
        std::string name = ggml_backend_dev_name(dev);
        std::string description = ggml_backend_dev_description(dev);
        // ignore Metal in CI
        if (name == "Metal" && description.find("Paravirtual") != std::string::npos) {
            continue;
        }
        devices.push_back(name);
    }
    return devices;
}

void ProstT5::perf() {
    llama_perf_context_print(ctx);
}
