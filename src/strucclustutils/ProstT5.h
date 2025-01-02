#ifndef PROST_T5_H
#define PROST_T5_H

#include <string>
#include <vector>

struct llama_model;
struct llama_context;

class ProstT5 {
public:
    ProstT5(const std::string& model_file, std::string & device);
    ~ProstT5();

    static std::vector<std::string> getDevices();
    
    std::string predict(const std::string& aa);
    void perf();

    llama_model * model;
    llama_context * ctx;
};


#endif // PROST_T5_H