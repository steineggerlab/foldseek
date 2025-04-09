/**
 * File: discretizer.cpp
 * Project: foldcomp
 * Created: 2021-02-05 13:41:54
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Functions for discretizing float values and restoring them
 * ---
 * Last Modified: 2022-12-09 15:42:34
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "discretizer.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream> // IWYU pragma: keep
#include <iostream>

Discretizer::Discretizer(const std::vector<float>& values, unsigned int nb):
    n_bin(nb) {
    if (values.size() == 0) {
        return;
    }
    // Get min & max
    this->min = *std::min_element(values.begin(), values.end());
    this->max = *std::max_element(values.begin(), values.end());
    // Calculate factors
    this->disc_f = this->n_bin / (this->max - this->min);
    this->cont_f = (this->max - this->min) / this->n_bin;
}

void Discretizer::set_continuous_values(const std::vector<float>& values) {
    this->min = *std::min_element(values.begin(), values.end());
    this->max = *std::max_element(values.begin(), values.end());
    // Calculate factors
    this->disc_f = this->n_bin / (this->max - this->min);
    this->cont_f = (this->max - this->min) / this->n_bin;
}

std::vector<unsigned int> Discretizer::discretize(const std::vector<float>& continuous_values) {
    std::vector<float>::const_iterator it;
    unsigned int tmp_disc_value;
    std::vector<unsigned int> discretizedValues;
    discretizedValues.reserve(continuous_values.size());
    for (it = continuous_values.cbegin(); it != continuous_values.cend(); it++) {
        tmp_disc_value = (unsigned int)((*it - min) * (this->disc_f) + 0.5);
        discretizedValues.push_back(tmp_disc_value);
    }
    return discretizedValues;
}

unsigned int Discretizer::discretize(float continuous_value) {
    return (continuous_value - this->min) * (this->disc_f);
}

std::vector<float> Discretizer::continuize(const std::vector<unsigned int>& discrete_values) {
    std::vector<unsigned int>::const_iterator it;
    std::vector<float> output;
    output.reserve(discrete_values.size());
    for (it = discrete_values.cbegin(); it != discrete_values.cend(); it++) {
        float tmp_cont_value = (*it * this->cont_f) + this->min;
        output.push_back(tmp_cont_value);
    }
    return output;
}

float Discretizer::continuize(unsigned int discrete_value) {
    return (discrete_value * this->cont_f) + this->min;
}

DiscParams Discretizer::get_param() {
    DiscParams params;
    params.min = this->min;
    params.max = this->max;
    params.n_bin = this->n_bin;
    params.disc_f = this->disc_f;
    params.cont_f = this->cont_f;
    return params;
}

// Methods for tests

void Discretizer::print() {
    std::cout << "MIN: " << this->min << std::endl;
    std::cout << "MAX: " << this->max << std::endl;
    std::cout << "N_BIN: " << this->n_bin << std::endl;
    std::cout << "DISC_F: " << this->disc_f << std::endl;
    std::cout << "CONT_F: " << this->cont_f << std::endl;
}

void Discretizer::write_to_file(std::string filename) {
    std::ofstream fout(filename);
    fout << "#MIN:" << this->min << "\n";
    fout << "#MAX:" << this->max << "\n";
    fout << "#N_BIN:" << this->n_bin << "\n";
    fout << "#DISC_F:" << this->disc_f << "\n";
    fout << "#CONT_F:" << this->cont_f << "\n";
    fout << "ORIGINAL_VALUES,DISCRETIZED_VALUES\n";
    fout.close();
}

float Discretizer::average_error(const std::vector<float>& continuous_values) {
    std::vector<unsigned int> discretized_values = this->discretize(continuous_values);
    std::vector<float> restored = this->continuize(discretized_values);
    float sum = 0;
    for (size_t i = 0; i < continuous_values.size(); i++) {
        sum += std::abs(continuous_values[i] - restored[i]);
    }
    return sum / continuous_values.size();
}

float Discretizer::max_error(const std::vector<float>& continuous_values) {
    std::vector<unsigned int> discretized_values = this->discretize(continuous_values);
    std::vector<float> restored = this->continuize(discretized_values);
    float max = 0;
    for (size_t i = 0; i < continuous_values.size(); i++) {
        if (std::abs(continuous_values[i] - restored[i]) > max) {
            max = std::abs(continuous_values[i] - restored[i]);
        }
    }
    return max;
}
