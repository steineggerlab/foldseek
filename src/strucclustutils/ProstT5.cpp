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

    // clear previous kv_cache values (irrelevant for embeddings)
    llama_kv_cache_clear(ctx);
    llama_set_embeddings(ctx, true);
    // run model
    if (llama_model_has_encoder(model) && !llama_model_has_decoder(model)) {
        if (llama_encode(ctx, llama_batch_get_one(enc_input.data(), enc_input.size())) < 0) {
            // LOG_ERR("%s : failed to encode\n", __func__);
            return 1;
        }
    } else { 
        // LOG_ERR("%s : no encoder\n", __func__);
        return 1;
    }
    // Log the embeddings (assuming n_embd is the embedding size per token)
    // LOG_INF("%s: n_tokens = %zu, n_seq = %d\n", __func__, enc_input.size(), 1);
    float* embeddings = llama_get_embeddings(ctx);
    if (embeddings == nullptr) {
        // LOG_ERR("%s : failed to retrieve embeddings\n", __func__);
        return 1;
    }
    int * arg_max_idx = new int[enc_input.size()];
    float * arg_max = new float[enc_input.size()];
    std::fill(arg_max, arg_max + enc_input.size(), std::numeric_limits<float>::min());
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
    // printf("\n");
    delete [] arg_max_idx;
    delete [] arg_max;
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
        // throw std::invalid_argument("no devices specified");
        return {};
    }
    if (dev_names.size() == 1 && dev_names[0] == "none") {
        devices.push_back(nullptr);
    } else {
        for (const auto & device : dev_names) {
            auto * dev = ggml_backend_dev_by_name(device.c_str());
            if (!dev || ggml_backend_dev_type(dev) != GGML_BACKEND_DEVICE_TYPE_GPU) {
                // throw std::invalid_argument(string_format("invalid device: %s", device.c_str()));
                return {};
            }
            devices.push_back(dev);
        }
        devices.push_back(nullptr);
    }
    return devices;
}

struct lora_adapter_info {
    std::string path;
    float scale;
};

struct lora_adapter_container : lora_adapter_info {
    struct llama_lora_adapter * adapter;
};

struct init_result {
    struct llama_model   * model   = nullptr;
    struct llama_context * context = nullptr;
    std::vector<lora_adapter_container> lora_adapters;
};

struct cpu_params {
    int      n_threads                   = -1;
    bool     cpumask[GGML_MAX_N_THREADS] = {false}; // CPU affinity mask.
    bool     mask_valid                  = false;   // Default: any CPU
    enum ggml_sched_priority  priority   = GGML_SCHED_PRIO_NORMAL;  // Scheduling prio : (0 - normal, 1 - medium, 2 - high, 3 - realtime)
    bool     strict_cpu                  = false;   // Use strict CPU placement
    uint32_t poll                        = 50;      // Polling (busywait) level (0 - no polling, 100 - mostly polling)
};

struct common_params {
    int32_t n_ctx                 =  4096; // context size
    int32_t n_batch               =  2048; // logical batch size for prompt processing (must be >=32 to use BLAS)
    int32_t n_ubatch              =   512; // physical batch size for prompt processing (must be >=32 to use BLAS)
    int32_t n_parallel            =     1; // number of parallel sequences to decode
    // float   rope_freq_base        =  0.0f; // RoPE base frequency
    // float   rope_freq_scale       =  0.0f; // RoPE frequency scaling factor
    // float   yarn_ext_factor       = -1.0f; // YaRN extrapolation mix factor
    // float   yarn_attn_factor      =  1.0f; // YaRN magnitude scaling factor
    // float   yarn_beta_fast        = 32.0f; // YaRN low correction dim
    // float   yarn_beta_slow        =  1.0f; // YaRN high correction dim
    // int32_t yarn_orig_ctx         =     0; // YaRN original context length
    float   defrag_thold          =  0.1f; // KV cache defragmentation threshold

    // // offload params
    std::vector<ggml_backend_dev_t> devices; // devices to use for offloading

    int32_t n_gpu_layers      = -1;  // number of layers to store in VRAM (-1 - use default)
    int32_t main_gpu          = 0;   // the GPU that is used for scratch and small tensors
    float   tensor_split[128] = {0}; // how split tensors should be distributed across GPUs

    enum llama_split_mode split_mode = LLAMA_SPLIT_MODE_LAYER; // how to split the model across GPUs

    struct cpu_params cpuparams;
    struct cpu_params cpuparams_batch;

    ggml_backend_sched_eval_callback cb_eval = nullptr;
    void * cb_eval_user_data                 = nullptr;

    ggml_numa_strategy numa = GGML_NUMA_STRATEGY_DISABLED;

    enum llama_rope_scaling_type rope_scaling_type = LLAMA_ROPE_SCALING_TYPE_UNSPECIFIED;
    enum llama_pooling_type      pooling_type      = LLAMA_POOLING_TYPE_UNSPECIFIED; // pooling type for embeddings
    enum llama_attention_type    attention_type    = LLAMA_ATTENTION_TYPE_UNSPECIFIED; // attention type for embeddings

    std::string model                = ""; // model path                                                    // NOLINT
    std::string rpc_servers          = ""; // comma separated list of RPC servers                           // NOLINT

    bool lora_init_without_apply = false; // only load lora to memory, but do not apply it to ctx (user can manually apply lora later using llama_lora_adapter_apply)
    std::vector<lora_adapter_info> lora_adapters; // lora adapter path with user defined scale

    bool flash_attn        = false; // flash attention
    bool no_perf           = false; // disable performance metrics
    bool logits_all        = false; // return logits for all tokens in the batch
    bool use_mmap          = true;  // use mmap for faster loads
    bool use_mlock         = false; // use mlock to keep model in memory
    bool no_kv_offload     = false; // disable KV offloading
    bool warmup            = true;  // warmup run
    bool check_tensors     = false; // validate tensor data

    ggml_type cache_type_k = GGML_TYPE_F16; // KV cache data type for the K
    ggml_type cache_type_v = GGML_TYPE_F16; // KV cache data type for the V

    bool embedding         = true; // get only sentence embedding
};

static struct init_result init_from_params(common_params & params) {
    init_result iparams;
    auto mparams = llama_model_default_params();

    if (!params.devices.empty()) {
        mparams.devices = params.devices.data();
    }
    if (params.n_gpu_layers != -1) {
        mparams.n_gpu_layers = params.n_gpu_layers;
    }
    mparams.rpc_servers     = params.rpc_servers.c_str();
    mparams.main_gpu        = params.main_gpu;
    mparams.split_mode      = params.split_mode;
    mparams.tensor_split    = params.tensor_split;
    mparams.use_mmap        = params.use_mmap;
    mparams.use_mlock       = params.use_mlock;
    mparams.check_tensors   = params.check_tensors;
    mparams.n_gpu_layers = 24;
    mparams.kv_overrides = NULL;

    llama_model * model = nullptr;

    model = llama_load_model_from_file(params.model.c_str(), mparams);

    if (model == NULL) {
        // LOG_ERR("%s: failed to load model '%s'\n", __func__, params.model.c_str());
        return iparams;
    }

    auto cparams = llama_context_default_params();

    cparams.n_ctx             = params.n_ctx;
    cparams.n_seq_max         = params.n_parallel;
    cparams.n_batch           = params.n_batch;
    cparams.n_ubatch          = params.n_ubatch;
    cparams.n_threads         = params.cpuparams.n_threads;
    cparams.n_threads_batch   = params.cpuparams_batch.n_threads == -1 ?
                                params.cpuparams.n_threads : params.cpuparams_batch.n_threads;
    cparams.logits_all        = params.logits_all;
    cparams.embeddings        = params.embedding;
    // cparams.rope_scaling_type = params.rope_scaling_type;
    // cparams.rope_freq_base    = params.rope_freq_base;
    // cparams.rope_freq_scale   = params.rope_freq_scale;
    // cparams.yarn_ext_factor   = params.yarn_ext_factor;
    // cparams.yarn_attn_factor  = params.yarn_attn_factor;
    // cparams.yarn_beta_fast    = params.yarn_beta_fast;
    // cparams.yarn_beta_slow    = params.yarn_beta_slow;
    // cparams.yarn_orig_ctx     = params.yarn_orig_ctx;
    cparams.pooling_type      = params.pooling_type;
    cparams.attention_type    = params.attention_type;
    cparams.defrag_thold      = params.defrag_thold;
    cparams.cb_eval           = params.cb_eval;
    cparams.cb_eval_user_data = params.cb_eval_user_data;
    cparams.offload_kqv       = !params.no_kv_offload;
    cparams.flash_attn        = params.flash_attn;
    cparams.no_perf           = params.no_perf;

    cparams.type_k = params.cache_type_k;
    cparams.type_v = params.cache_type_v;

    llama_context * lctx = llama_new_context_with_model(model, cparams);
    if (lctx == NULL) {
        // LOG_ERR("%s: failed to create context with model '%s'\n", __func__, params.model.c_str());
        llama_free_model(model);
        return iparams;
    }


    // load and optionally apply lora adapters
    for (auto & la : params.lora_adapters) {
        lora_adapter_container loaded_la;
        loaded_la.path = la.path;
        loaded_la.scale = la.scale;
        loaded_la.adapter = llama_lora_adapter_init(model, la.path.c_str());
        if (loaded_la.adapter == nullptr) {
            // LOG_ERR("%s: failed to apply lora adapter '%s'\n", __func__, la.path.c_str());
            llama_free(lctx);
            llama_free_model(model);
            return iparams;
        }
        iparams.lora_adapters.push_back(loaded_la); // copy to list of loaded adapters
    }
    if (!params.lora_init_without_apply) {
        llama_lora_adapter_clear(lctx);
        for (auto & la : iparams.lora_adapters) {
            if (la.scale != 0.0f) {
                llama_lora_adapter_set(lctx, la.adapter, la.scale);
            }
        }
    }

    iparams.model   = model;
    iparams.context = lctx;

    return iparams;
}

struct llama_model;
struct llama_context;

ProstT5::ProstT5(const std::string& model_file, std::string & device) {
    llama_log_set([](ggml_log_level, const char *, void *) {}, NULL);

    ggml_backend_load_all();

    common_params params;
    params.n_ubatch = params.n_batch;
    params.warmup = false;
    params.model = model_file;
    params.cpuparams.n_threads = 1;
    params.use_mmap = true;
    params.devices = parse_device_list(device);


    llama_backend_init();
    llama_numa_init(params.numa);

    // load the model
    init_result llama_init = init_from_params(params);

    model = llama_init.model;
    ctx = llama_init.context;

};

ProstT5::~ProstT5() {
    // clean up
    llama_free(ctx);
    llama_free_model(model);
    llama_backend_free();
}

std::string ProstT5::predict(const std::string& aa) {
    std::string result;
    std::vector<llama_token> embd_inp;

    embd_inp.reserve(aa.length() + 2);
    embd_inp.emplace_back(llama_token_get_token(model, "<AA2fold>"));
    llama_token unk_aa = llama_token_get_token(model, "▁X");
    for (size_t i = 0; i < aa.length(); ++i) {
        std::string current_char("▁");
        current_char.append(1, toupper(aa[i]));
        llama_token token = llama_token_get_token(model, current_char.c_str());
        if (token == LLAMA_TOKEN_NULL) {
            embd_inp.emplace_back(unk_aa);
        } else {
            embd_inp.emplace_back(token);
        }
    }
    embd_inp.emplace_back(llama_token_get_token(model, "</s>"));


    encode(ctx, embd_inp, result);
    return result;
}

std::vector<std::string > ProstT5::getDevices() {
    std::vector<std::string> devices;
    for (size_t i = 0; i < ggml_backend_dev_count(); ++i) {
        ggml_backend_dev_t dev = ggml_backend_dev_get(i);
        std::string name = ggml_backend_dev_name(dev);
        devices.push_back(name);
    }
    return devices;
}

void ProstT5::perf() {
    llama_perf_context_print(ctx);
}
