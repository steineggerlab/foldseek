/**
 * File: torsion_angle.cpp
 * Project: foldcomp
 * Created: 2021-01-13 10:34:34
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code calculates torsion angles from the atomic coordinates.
 * Reference:
 *     1) "torsion.xyz.R" in Bio3D R package
 *     https://rdrr.io/cran/bio3d/man/torsion.xyz.html
 *     2) "XYZ.h" in pdbtools
 *     https://github.com/realbigws/PDB_Tool
 * ---
 * Last Modified: 2022-09-13 15:15:46
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "torsion_angle.h"

#include "atom_coordinate.h"

#include <cstddef>
#include <iostream>
#include <string>

 // Temp function
void print3DFloatVec(std::string name, const std::vector<float>& input) {
    std::cout << name << ": " << input[0] << ", " << input[1] << ", " << input[2] << std::endl;
}

std::vector<float> getTorsionFromXYZ(
    const std::vector<AtomCoordinate>& coordinates, int atm_inc
) {
    std::vector<float3d> coord = extractCoordinates(coordinates);
    return getTorsionFromXYZ(coord, atm_inc);
}

/**
 * @brief Get the torsion from xyz object
 *
 * @param coordinates a vector of float vectors
 * @param atm_inc an integer
 * @return std::vector<std::vector<float>>
 */
std::vector<float> getTorsionFromXYZ(
    const std::vector<float3d>& coordinates, int atm_inc = 1
) {
    std::vector<float> torsion_vector;
    for (size_t i = 0; i < (coordinates.size() - 3); i += atm_inc) {
        // 00. Get 4 atom coordinates
        float3d atm_1 = coordinates[i + 0];
        float3d atm_2 = coordinates[i + 1];
        float3d atm_3 = coordinates[i + 2];
        float3d atm_4 = coordinates[i + 3];
        // 01. Obtain vectors from coordinates
        float3d d1 {
            (atm_2.x - atm_1.x), (atm_2.y - atm_1.y), (atm_2.z - atm_1.z)
        };
        float3d d2{
            (atm_3.x - atm_2.x), (atm_3.y - atm_2.y), (atm_3.z - atm_2.z)
        };
        float3d d3{
            (atm_4.x - atm_3.x), (atm_4.y - atm_3.y), (atm_4.z - atm_3.z)
        };
        // 02. Calculate cross product
        float3d u1 = crossProduct(d1, d2);
        float3d u2 = crossProduct(d2, d3);

        // 03. Get cosine theta of u1 & u2
        float cos_torsion = getCosineTheta(u1, u2);
        // 04. Calculate arc cosine
        float torsion_angle;
        if (std::isnan(acos(cos_torsion))) {
            if (cos_torsion < 0) {
                torsion_angle = 180.0;
            } else {
                torsion_angle = 0.0;
            }
        } else {
            torsion_angle = acos(cos_torsion) * 180.0 / M_PI;
        }

        // 05. check torsion.xyz.R line 37
        // IMPORTANT: torsion_angle should be minus in specific case;
        // the vector on the plane beta can be represented as
        // cross product of u2 & d2
        float3d plane_beta_vec = crossProduct(u2, d2);
        // determinate of u1 & plane_beta_vec
        if ((u1.x * plane_beta_vec.x) + (u1.y * plane_beta_vec.y) + (u1.z * plane_beta_vec.z) < 0) {
            torsion_angle = -1 * torsion_angle;
        }
        torsion_vector.push_back(torsion_angle);
    }
    return torsion_vector;
}


/**
 * WARNING: TEMPORARY FUNCTION
 * @brief Fill an array of double with float vector
 *
 * @param fv
 * @param output
 */
void float3dVectorToDoubleArray(const std::vector<float>& fv, double output[3]) {
    for (int i = 0; i < 3; i++) {
        output[i] = static_cast<double>(fv[i]);
    }
}


/**
 * @brief Save torsion angles to the output
 *
 * @param output_path
 * @param torsion
 */
void writeTorsionAngles(const std::string& file_path, const std::vector<float>& torsion) {
    std::ofstream output;
    output.open(file_path, std::ios_base::out);
    for (float angle : torsion) {
        output << angle << "\n";
    }
    output.close();
}

std::vector<float> readTorsionAngles(const std::string& file_path) {
    std::ifstream input_file;
    std::string line;
    std::vector<float> output;
    float angle;
    input_file.open(file_path, std::ios_base::in);
    while(getline(input_file, line)) {
        angle = std::stof(line);
        output.push_back(angle);
    }
    input_file.close();
    return output;
}

// Encode backbone by encoding dihedrals by 2B each into intervals of 2 pi / 65536
// Current, the degree is not in radian and 360 will be used instead of 2 pi.

/**
 * @brief Encode float-based torsion angles to short
 *
 * @param torsions A float vector of torsion angles
 * @param n_bits Number of bits for encoding
 * @return std::vector<short>
 */
std::vector<short> encodeTorsionAnglesToShort(
    const std::vector<float>& torsions, unsigned int n_bits
) {
    std::vector<short> output;
    short s_angle;
    for (float f_angle: torsions) {
        s_angle = (short)(f_angle * pow(2, n_bits) / 360);
        output.push_back(s_angle);
    }
    return output;
}

/**
 * @brief Decode short-encoded torsion angles to float
 *
 * @param encoded_torsions A short vector of torsion angles
 * @param n_bits Number of bits for encoding
 * @return std::vector<float>
 */
std::vector<float> decodeEncodedTorsionAngles(
    const std::vector<short>& encoded_torsions, unsigned int n_bits
) {
    std::vector<float> output;
    float f_angle;
    for (short s_angle: encoded_torsions) {
        f_angle = (float)(s_angle);
        f_angle = f_angle * 360.0 / pow(2, n_bits);
        output.push_back(f_angle);
    }
    return output;
}

