/**
 * File: nerf.h
 * Project: foldcomp
 * Created: 2021-01-11 16:42:32
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Contributor: Milot Mirdita (milot@mirdita.de)
 * Description:
 *     This file contains NeRF (Natural Extension of Reference Frame) algorithm
 *     implementation. NeRF is a method to calculate upcoming atoms' coordinates
 *     with three preceding atoms and bond information.
 *     Reference: https://benjamin.computer/posts/2018-03-16-mres-part2.html
 * ---
 * Last Modified: 2022-09-29 17:56:31
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

// COMMENT FROM REFERENCE
// # TODO - PROLINE has different lengths which we should take into account
// # TODO - A_TO_C angle differs by + / -5 degrees
#pragma once
#include "float3d.h"

#include <map>
#include <string>
#include <vector>

class AminoAcid;
class AtomCoordinate;

class Nerf {
private:
    // private members
public:
    /* data */
    const std::map<std::string, float> bond_lengths {
        // TODO: if residue is proline, apply different bond length
        // 2021-01-19 15:20:42 NOTE: TEST FORE SERINE
        {"N_TO_CA", 1.4581f}, {"PRO_N_TO_CA", 1.353f}, // proline has different lengths
        {"CA_TO_C", 1.5281f}, {"C_TO_N", 1.3311f},
        {"C_TO_O", 1.23f}
    };
    // bond angles are in radian
    const std::map<std::string, float> bond_angles {
        {"N_TO_CA", 121.3822f}, {"CA_TO_C", 111.2812f}, {"C_TO_N", 116.6429f},
        {"C_TO_O", 120.5117f}
    };
    //NOTE: bond_order??
    // PDB file is ordered from N terminal to C terminal
    const std::vector<std::string> bond_order { "C_TO_N", "N_TO_CA", "CA_TO_C" };

    float3d place_atom(
        const float3d* prev_atoms,
        float bond_length, float bond_angle, float torsion_angle
    );

    std::vector<AtomCoordinate> reconstructWithAAMaps(
        const std::vector<AtomCoordinate>& original_atoms,
        const std::map <std::string, std::vector<std::string>>& prev_atom_map,
        const std::map<std::string, float>& aa_torsion_angles,
        const std::map<std::string, float>& aa_bond_lengths,
        const std::map<std::string, float>& aa_bond_angles
    );

    std::vector<AtomCoordinate> reconstructWithDynamicAngles(
        const std::vector<AtomCoordinate>& original_atoms,
        const std::vector<float>& torsion_angles,
        const std::vector<float>& atom_bond_lengths,
        const std::vector<float>& atom_bond_angles
    );

    std::vector<AtomCoordinate> reconstructWithDynamicAngles(
        const std::vector<AtomCoordinate>& original_atoms,
        const std::vector<float>& torsion_angles,
        const std::vector<float>& atom_bond_angles
    );

    std::vector<AtomCoordinate> reconstructWithBreaks(
        const std::vector<AtomCoordinate>& original_atoms,
        const std::vector<float>& torsion_angles,
        const std::vector<float>& atom_bond_angles,
        const std::vector<int>& break_indices
    );

    std::vector<AtomCoordinate> reconstructWithReversed(
        std::vector<AtomCoordinate> original_atoms,
        std::vector<float> torsion_angles,
        std::vector<float> atom_bond_angles
    );

    std::vector<AtomCoordinate> reconstructAminoAcid(
        const std::vector<AtomCoordinate>& original_atoms,
        const std::vector<float>& torsion_angles,
        const AminoAcid& aa
    );

    void writeInfoForChecking(
        const std::vector<AtomCoordinate>& coord_list,
        const std::string& output_path,
        const std::string sep = ","
    );
    void writeCoordinatesBinary(
        const std::vector<AtomCoordinate>& coord_list,
        const std::string& output_path
    );
    std::vector<float> getBondAngles(const std::vector<AtomCoordinate>& original_atoms);
    std::vector<float> getBondLengths(const std::vector<AtomCoordinate>& original_atoms);

    // 2021-07-12 17:18:16
    // Method to identify breaks from AtomCoordinate vector
    std::vector<int> identifyBreaks(
        const std::vector<AtomCoordinate>& original_atoms, float cutoff = 2
    );

};
