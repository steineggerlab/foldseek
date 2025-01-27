#ifndef PROST_T5_H
#define PROST_T5_H

#include <string>
#include <vector>

struct llama_model;
struct llama_context;

class LlamaInitGuard {
public:
    explicit LlamaInitGuard(bool verbose = false);
    ~LlamaInitGuard();

    LlamaInitGuard(const LlamaInitGuard&) = delete;
    LlamaInitGuard& operator=(const LlamaInitGuard&) = delete;
};

class ProstT5Model {
public:
    ProstT5Model(const std::string& model_file, std::string& device);
    ~ProstT5Model();

    llama_model* model;
};

class ProstT5 {
public:
    ProstT5(ProstT5Model& model, int threads);
    ~ProstT5();

    static std::vector<std::string> getDevices();
    
    std::string predict(const std::string& aa);
    void perf();

    ProstT5Model& model;
    llama_context* ctx;
};


#endif // PROST_T5_H