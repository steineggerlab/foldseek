/**
 * File: torsion_angle.h
 * Project: foldcomp
 * Created: 2021-01-13 10:35:56
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Functions for torsion angle calculation
 * ---
 * Last Modified: 2022-09-13 15:14:10
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#pragma once
#include "float3d.h"

#include <iosfwd>
#include <string> // IWYU pragma: keep
#include <vector>

class AtomCoordinate;

std::vector<float> getTorsionFromXYZ(
    const std::vector<float3d>& coordinates, int atm_inc
);

std::vector<float> getTorsionFromXYZ(
    const std::vector<AtomCoordinate>& coordinates, int atm_inc
);

std::vector<float> getTorsionFromXYZ(float3d coordinates[4], int atm_inc);

void float3dVectorToDoubleArray(const std::vector<float>& fv, double output[3]);

void writeTorsionAngles(const std::string& file_path, const std::vector<float>& torsion);
std::vector<float> readTorsionAngles(const std::string& file_path);


// 2021-02-02 13:02:02 - encode and decode torsion angles to short

std::vector<short> encodeTorsionAnglesToShort(
    const std::vector<float>& torsions, unsigned int n_bits = 16
);
std::vector<float> decodeEncodedTorsionAngles(
    const std::vector<short>& encoded_torsions, unsigned int n_bits = 16
);