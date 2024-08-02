#include <fstream>
#include <iostream>
#include <chrono>
#include <vector>
#include <string>

#include <prostt5.h>
#include <cxxopts.hpp>

int main(int argc, char* argv[]) {
  cxxopts::Options cmd_options("ct2-translator", "CTranslate2 translator client");
  cmd_options.custom_help("--model <directory> [OPTIONS]");

  cmd_options.add_options("General")
    ("h,help", "Display available options.")
    ("seed", "Seed value of the random generators.",
     cxxopts::value<unsigned int>()->default_value("0"))
    ;

  cmd_options.add_options("Device")
    ("inter_threads", "Maximum number of CPU translations to run in parallel.",
     cxxopts::value<size_t>()->default_value("1"))
    ("intra_threads", "Number of computation threads (set to 0 to use the default value).",
     cxxopts::value<size_t>()->default_value("0"))
    ("device", "Device to use (can be cpu, cuda, auto).",
     cxxopts::value<std::string>()->default_value("cpu"))
    ("device_index", "Comma-separated list of device IDs to use.",
     cxxopts::value<std::vector<int>>()->default_value("0"))
    ("cpu_core_offset", "Pin worker threads to CPU cores starting from this offset.",
     cxxopts::value<int>()->default_value("-1"))
    ;

  cmd_options.add_options("Model")
    ("model", "Path to the CTranslate2 model directory.", cxxopts::value<std::string>())
    ("compute_type", "The type used for computation: default, auto, float32, float16, bfloat16, int16, int8, int8_float32, int8_float16, or int8_bfloat16",
     cxxopts::value<std::string>()->default_value("default"))
    ("cuda_compute_type", "Computation type on CUDA devices (overrides compute_type)",
     cxxopts::value<std::string>()->default_value(""))
    ("cpu_compute_type", "Computation type on CPU devices (overrides compute_type)",
     cxxopts::value<std::string>()->default_value(""))
    ;

  cmd_options.add_options("Data")
    ("src", "Path to the source file (read from the standard input if not set).",
     cxxopts::value<std::string>())
    ("out", "Path to the output file (write to the standard output if not set).",
     cxxopts::value<std::string>())
    ("batch_size", "Size of the batch to forward into the model at once.",
     cxxopts::value<size_t>()->default_value("32"))
    ("read_batch_size", "Size of the batch to read at once (defaults to batch_size).",
     cxxopts::value<size_t>()->default_value("0"))
    ("max_queued_batches", "Maximum number of batches to load in advance (set -1 for unlimited, 0 for an automatic value).",
     cxxopts::value<long>()->default_value("0"))
    ("batch_type", "Batch type (can be examples, tokens).",
     cxxopts::value<std::string>()->default_value("examples"))
    ("max_input_length", "Truncate inputs after this many tokens (set 0 to disable).",
     cxxopts::value<size_t>()->default_value("1024"))
    ;

  auto args = cmd_options.parse(argc, argv);
  if (args.count("help")) {
    std::cerr << cmd_options.help() << std::endl;
    return 0;
  }

  ProstT5 prost_t5(
      args["model"].as<std::string>(),
      args["device"].as<std::string>(),
      args["compute_type"].as<std::string>(),
      args["device_index"].as<std::vector<int>>(),
      args["inter_threads"].as<size_t>(),
      args["intra_threads"].as<size_t>(),
      args["max_queued_batches"].as<long>(),
      args["cpu_core_offset"].as<int>(),
      args["cpu_compute_type"].as<std::string>(),
      args["cuda_compute_type"].as<std::string>(),
      args.count("seed") ? args["seed"].as<unsigned int>() : 0
  );

  auto start = std::chrono::high_resolution_clock::now();

  std::vector<std::string> input_seqs = {args["src"].as<std::string>()};
  std::vector<std::string> predictions = prost_t5.predict(input_seqs);
  for (const auto& prediction : predictions) {
      std::cout << prediction << std::endl;
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;

  std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;

  return 0;
}
