/**
 * File: discretizer.h
 * Project: foldcomp
 * Created: 2021-02-05 13:41:21
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Functions for discretizing float values and restoring them
 * ---
 * Last Modified: 2022-07-20 01:56:21
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#pragma once
#include <string> // IWYU pragma: keep
#include <vector>

#define MIN_ANGLE -180.0
#define MAX_ANLGE 180.0

struct DiscParams {
    /* data */
    // Discretizer - min, max, n_bin, df, cf
    float min;
    float max;
    int n_bin;
    float disc_f;
    float cont_f;
};

/**
 * @brief Discretize float vector
 *
 */
class Discretizer {
private:
    /* data */
public:
    Discretizer() = default;
    Discretizer(const Discretizer&) = default;
    Discretizer(Discretizer&&) = default;
    Discretizer& operator=(const Discretizer&) = default;
    Discretizer& operator=(Discretizer&&) = default;

    /**
     * @brief Construct a new Discretizer object (start with values)
     *
     * @param values a float vector to discretize
     * @param nb an int representing the number of bins
     */
    Discretizer(const std::vector<float>& values, unsigned int nb);
    /**
     * @brief Construct a new Discretizer object (without values)
     *
     * @param min_ the minimum value
     * @param max_ the maximum value
     * @param nb   an int representing the number of bins
     * @param df   a float factor for discretization
     * @param cf   a float factor for continuization
     */
    Discretizer(float min_, float max_, unsigned int nb, float df, float cf):
        min(min_), max(max_), n_bin(nb), disc_f(df), cont_f(cf){};

    float min;
    float max;
    unsigned int n_bin;
    float disc_f; // discrete factor:
    float cont_f; // continous factor:

    // methods
    void set_continuous_values(const std::vector<float>& values);


    std::vector<unsigned int> discretize(const std::vector<float>& continuous_values);
    unsigned int discretize(float continuous_value);

    std::vector<float> continuize(const std::vector<unsigned int>& discrete_values);
    float continuize(unsigned int discrete_value);

    DiscParams get_param();

    // Methods for tests
    void print();
    void write_to_file(std::string filename);
    float average_error(const std::vector<float>& continuous_values);
    float max_error(const std::vector<float>& continuous_values);
};

class FixedAngleDiscretizer: public Discretizer {
public:
    FixedAngleDiscretizer(int nb) {
        this->n_bin = nb;
        this->min = MIN_ANGLE;
        this->max = MAX_ANLGE;
        this->disc_f = this->n_bin / (this->max - this->min);
        this->cont_f = (this->max - this->min) / this->n_bin;
    };
    FixedAngleDiscretizer(const std::vector<float>& /* values */, unsigned int nb) {
        this->n_bin = nb;
        this->min = MIN_ANGLE;
        this->max = MAX_ANLGE;
        this->disc_f = this->n_bin / (this->max - this->min);
        this->cont_f = (this->max - this->min) / this->n_bin;
    };
    ~FixedAngleDiscretizer(){};
};

// 2022-06-08 13:40:05
// DONE: REWROTE THE DISCRETIZER TO REMOVE VECTORS IN IT
// CHECK AND FIND ALL PARTS WHERE DISCRETIZER IS USED