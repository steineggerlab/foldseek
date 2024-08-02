#ifndef PROSTT5_H
#define PROSTT5_H

#include <string>
#include <vector>
#include <memory>

class ProstT5 {
public:
    ProstT5(
        const std::string& model_dir,
        const std::string& device_str,
        const std::vector<int>& device_indices,
        size_t inter_threads,
        size_t intra_threads,
        long max_queued_batches,
        int cpu_core_offset,
        const std::string& cpu_compute_type = "",
        const std::string& cuda_compute_type = "",
        unsigned int seed = 0
    );

    ~ProstT5();

    std::vector<std::string> predict(const std::vector<std::string>& input_seqs);

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

#endif
