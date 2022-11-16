/**
 * File: nerf.cpp
 * Project: foldcomp
 * Created: 2021-01-11 16:42:27
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This file contains NeRF (Natural Extension of Reference Frame) algorithm
 *     implementation. NeRF is a method to calculate upcoming atoms' coordinates
 *     with three preceding atoms and bond information.
 *     Reference: https://benjamin.computer/posts/2018-03-16-mres-part2.html
 * ---
 * Last Modified: 2022-09-29 17:47:07
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "nerf.h"

#include "amino_acid.h"
#include "atom_coordinate.h"

#include <algorithm>
#include <cstddef>
#include <fstream>

/**
 * @brief Calculate the position of next atom with torsion angle.
 *
 * @details
 * This code is a part of FoldU project. this function calculates the
 * position of the next atom with NeRF algorithm.
 * Angles are in radians, lengths in angstroms
 * @param prev_atoms XYZ coordinates of 3 previous atoms
 * @param bond_length a float, the length of current bond
 * @param bond_angle a float, the angle of current bond
 * @param torsion_angle a float, the torsion angle of current bond.
 * @return std::vector<float>
 */
float3d Nerf::place_atom(
    const float3d* prev_atoms,
    float bond_length, float bond_angle,
    float torsion_angle
){
    // 00. Get 3 atom coordinates
    float3d atm_a = prev_atoms[0];
    float3d atm_b = prev_atoms[1];
    float3d atm_c = prev_atoms[2];

    // 01. Obtain vectors from coordinates

    float3d ab {
        (atm_b.x - atm_a.x), (atm_b.y - atm_a.y), (atm_b.z - atm_a.z)
    };
    float3d bc {
        (atm_c.x - atm_b.x), (atm_c.y - atm_b.y), (atm_c.z - atm_b.z)
    };
    float bc_norm = norm(bc);

    // n2 - unit vector of direction same with d2
    float3d bcn {
        (bc.x / bc_norm), (bc.y / bc_norm), (bc.z / bc_norm)
    };
    // Current atom
    bond_angle = bond_angle * M_PI / 180.0;
    torsion_angle = torsion_angle * M_PI / 180.0;

    float3d curr_atm {
        (-1 * bond_length * cosf(bond_angle)), // x
        (bond_length * cosf(torsion_angle) * sinf(bond_angle)), // y
        (bond_length * sinf(torsion_angle) * sinf(bond_angle)) // z
    };

    // 02. Calculate cross product
    // TODO: check if we can reduce the count of norm calculation
    float3d n = crossProduct(ab, bcn);
    float n_norm = norm(n);
    n.x = n.x / n_norm;
    n.y = n.y / n_norm;
    n.z = n.z / n_norm;
    float3d nbc = crossProduct(n, bcn);
    float3d m[3] {
        {bcn.x, nbc.x, n.x},
        {bcn.y, nbc.y, n.y},
        {bcn.z, nbc.z, n.z}
    };
    //
    float3d atm_d {0.0, 0.0, 0.0};
    atm_d.x += (m[0].x * curr_atm.x);
    atm_d.x += (m[0].y * curr_atm.y);
    atm_d.x += (m[0].z * curr_atm.z);
    atm_d.y += (m[1].x * curr_atm.x);
    atm_d.y += (m[1].y * curr_atm.y);
    atm_d.y += (m[1].z * curr_atm.z);
    atm_d.z += (m[2].x * curr_atm.x);
    atm_d.z += (m[2].y * curr_atm.y);
    atm_d.z += (m[2].z * curr_atm.z);

    // 0X. Add the coordinates of atom c to atom d --> absolute coordinate
    atm_d.x += atm_c.x;
    atm_d.y += atm_c.y;
    atm_d.z += atm_c.z;

    return atm_d;
}

std::vector<AtomCoordinate> Nerf::reconstructAminoAcid(
    const std::vector<AtomCoordinate>& original_atoms,
    const std::vector<float>& torsion_angles,
    const AminoAcid& aa
) {
    // save three first atoms
    std::vector<AtomCoordinate> reconstructed_atoms = {
        original_atoms[0], original_atoms[1], original_atoms[2]
    };

    AtomCoordinate curr_atom(
        "", original_atoms[0].residue, original_atoms[0].chain,
        original_atoms[2].atom_index + 1, original_atoms[0].residue_index,
        0.0f, 0.0f, 0.0f
    );

    int total = aa.atoms.size();
    for (int i = 0; i < (total - 3); i++) {
        // Get current atom's info
        curr_atom.atom_index = reconstructed_atoms[i + 2].atom_index + 1;
        curr_atom.atom = aa.atoms[i + 3];

        const auto& side_chain = aa.sideChain.at(curr_atom.atom);

        // Fill prev_atoms
        const AtomCoordinate& prev_atom0 = findFirstAtom(reconstructed_atoms, side_chain[0]);
        const AtomCoordinate& prev_atom1 = findFirstAtom(reconstructed_atoms, side_chain[1]);
        const AtomCoordinate& prev_atom2 = findFirstAtom(reconstructed_atoms, side_chain[2]);

        float3d prev_coord[3];
        extractCoordinates(prev_coord, prev_atom0, prev_atom1, prev_atom2);
        // Bond name & Angle name
        std::string curr_bond_name = prev_atom2.atom + "_" + curr_atom.atom;
        std::string curr_angle_name = prev_atom1.atom + "_" + prev_atom2.atom + "_" + curr_atom.atom;

        float curr_bond_length = aa.bondLengths.at(curr_bond_name);
        float curr_bond_angle = aa.bondAngles.at(curr_angle_name);

        float3d curr_coord = place_atom(
            prev_coord, curr_bond_length, curr_bond_angle, torsion_angles[i]
        );

        reconstructed_atoms.emplace_back(
            curr_atom.atom, curr_atom.residue, curr_atom.chain,
            curr_atom.atom_index, curr_atom.residue_index,
            curr_coord.x, curr_coord.y, curr_coord.z
        );
    }
    return reconstructed_atoms;
}

std::vector<AtomCoordinate> Nerf::reconstructWithAAMaps(
    const std::vector<AtomCoordinate>& original_atoms,
    const std::map<std::string, std::vector<std::string>>& prev_atom_map,
    const std::map<std::string, float>& aa_torsion_angles,
    const std::map<std::string, float>& aa_bond_lengths,
    const std::map<std::string, float>& aa_bond_angles
) {
    std::map<std::string, AtomCoordinate> original_atom_map;
    for (const auto& atm : original_atoms) {
        original_atom_map[atm.atom] = atm;
    }

    std::vector<AtomCoordinate> reconstructed_atoms;
    reconstructed_atoms.reserve(original_atoms.size());

    // save three first atoms
    reconstructed_atoms.push_back(original_atoms[0]);
    reconstructed_atoms.push_back(original_atoms[1]);
    reconstructed_atoms.push_back(original_atoms[2]);

    int total = original_atoms.size();
    for (int i = 0; i < (total - 3); i++) {
        const AtomCoordinate& curr_atom = original_atoms[i + 3];
        const std::vector<std::string>& prev_atom_names = prev_atom_map.at(curr_atom.atom);
        std::vector<AtomCoordinate> prev_atoms;
        for (const auto& atm_name : prev_atom_names) {
            prev_atoms.push_back(original_atom_map[atm_name]);
        }
        std::vector<float3d> prev_coord_vec = extractCoordinates(prev_atoms);
        float3d prev_coord[3];
        for (int j = 0; j < 3; j++) {
            prev_coord[j] = prev_coord_vec[j];
        }

        std::string curr_bond_length_key = prev_atoms[2].atom + "_" + curr_atom.atom;
        std::string curr_bond_angle_key = prev_atoms[1].atom + "_" + prev_atoms[2].atom + "_" + curr_atom.atom;
        std::string curr_torsion_angle_key = prev_atoms[0].atom + "_" + prev_atoms[1].atom +
            "_" + prev_atoms[2].atom + "_" + curr_atom.atom;

        float curr_bond_length = aa_bond_lengths.at(curr_bond_length_key);
        float curr_bond_angle = aa_bond_angles.at(curr_bond_angle_key);
        // Found error here! 2021-10-15 01:03:03
        float curr_torsion_angle = aa_torsion_angles.at(curr_torsion_angle_key);
        float3d curr_coord = place_atom(
            prev_coord, curr_bond_length, curr_bond_angle, curr_torsion_angle
        );
        reconstructed_atoms.emplace_back(
            curr_atom.atom, curr_atom.residue, curr_atom.chain,
            curr_atom.atom_index, curr_atom.residue_index,
            curr_coord.x, curr_coord.y, curr_coord.z
        );
    }
    return reconstructed_atoms;
}




/**
 * @brief Reconstruct structure from torsion angles using atom-specific bond info
 *
 * @param original_atoms std::vector<AtomCoordinate>
 * @param torsion_angles
 * @param atom_bond_angles
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> Nerf::reconstructWithDynamicAngles(
    const std::vector<AtomCoordinate>& original_atoms,
    const std::vector<float>& torsion_angles,
    const std::vector<float>& atom_bond_lengths,
    const std::vector<float>& atom_bond_angles
) {
    std::vector<AtomCoordinate> reconstructed_atoms;
    reconstructed_atoms.reserve(original_atoms.size());

    // save three first atoms
    reconstructed_atoms.push_back(original_atoms[0]);
    reconstructed_atoms.push_back(original_atoms[1]);
    reconstructed_atoms.push_back(original_atoms[2]);

    int total = original_atoms.size();
    for (int i = 0; i < (total - 3); i++) {
        const AtomCoordinate& curr_atom = original_atoms[i + 3];
        float3d prev_coord[3];
        extractCoordinates(prev_coord, reconstructed_atoms[i], reconstructed_atoms[i + 1],
                reconstructed_atoms[i + 2]);
        std::string curr_bond_name = original_atoms[i + 2].atom + "_TO_" + curr_atom.atom;
        float curr_bond_length = atom_bond_lengths[i + 2];
        float curr_bond_angle = atom_bond_angles[i + 1];
        float3d curr_coord = place_atom(prev_coord, curr_bond_length, curr_bond_angle,
                            torsion_angles[i]);
        reconstructed_atoms.emplace_back(
            curr_atom.atom, curr_atom.residue, curr_atom.chain,
            curr_atom.atom_index, curr_atom.residue_index,
            curr_coord.x, curr_coord.y, curr_coord.z
        );
    }
    return reconstructed_atoms;

}

std::vector<AtomCoordinate> Nerf::reconstructWithDynamicAngles(
    const std::vector<AtomCoordinate>& original_atoms,
    const std::vector<float>& torsion_angles,
    const std::vector<float>& atom_bond_angles
) {
    std::vector<AtomCoordinate> reconstructed_atoms;
    reconstructed_atoms.reserve(original_atoms.size());

    // save three first atoms
    reconstructed_atoms.push_back(original_atoms[0]);
    reconstructed_atoms.push_back(original_atoms[1]);
    reconstructed_atoms.push_back(original_atoms[2]);

    int total = original_atoms.size();
    for (int i = 0; i < (total - 3); i++) {
        const AtomCoordinate& curr_atom = original_atoms[i + 3];
        float3d prev_coord[3];
        extractCoordinates(prev_coord, reconstructed_atoms[i], reconstructed_atoms[i + 1],
                    reconstructed_atoms[i + 2]);
        std::string curr_bond_name = original_atoms[i + 2].atom + "_TO_" + curr_atom.atom;
        float curr_bond_length = this->bond_lengths.at(curr_bond_name);
        float curr_bond_angle = atom_bond_angles[i + 1];
        float3d curr_coord = place_atom(prev_coord, curr_bond_length, curr_bond_angle,
                                torsion_angles[i]);
        reconstructed_atoms.emplace_back(
            curr_atom.atom, curr_atom.residue, curr_atom.chain,
            curr_atom.atom_index, curr_atom.residue_index,
            curr_coord.x, curr_coord.y, curr_coord.z
        );
    }
    return reconstructed_atoms;
}


std::vector<AtomCoordinate> Nerf::reconstructWithBreaks(
    const std::vector<AtomCoordinate>& original_atoms,
    const std::vector<float>& torsion_angles,
    const std::vector<float>& atom_bond_angles,
    const std::vector<int>& break_indices
) {
    int total = original_atoms.size();
    std::vector<AtomCoordinate> reconstructed_atoms;
    reconstructed_atoms.reserve(original_atoms.size());

    // save three first atoms
    for (size_t k = 0; k < break_indices.size(); k++) {
        int breakpoint = break_indices[k];
        int next_breakpoint;
        if (k != (break_indices.size() - 1)) {
            next_breakpoint = break_indices[k + 1];
        } else {
            next_breakpoint = total;
        }
        reconstructed_atoms.push_back(original_atoms[breakpoint]);
        reconstructed_atoms.push_back(original_atoms[breakpoint + 1]);
        reconstructed_atoms.push_back(original_atoms[breakpoint + 2]);

        for (int i = breakpoint; i < (next_breakpoint - 3); i++) {
            const AtomCoordinate& curr_atom = original_atoms[i + 3];
            float3d prev_coord[3];
            extractCoordinates(prev_coord, reconstructed_atoms[i], reconstructed_atoms[i + 1],
                          reconstructed_atoms[i + 2]);
            std::string curr_bond_name = original_atoms[i + 2].atom + "_TO_" + curr_atom.atom;
            float curr_bond_length = this->bond_lengths.at(curr_bond_name);
            float curr_bond_angle = atom_bond_angles[i + 1];
            float3d curr_coord = place_atom(prev_coord, curr_bond_length, curr_bond_angle,
                torsion_angles[i]);
            reconstructed_atoms.emplace_back(
                curr_atom.atom, curr_atom.residue, curr_atom.chain,
                curr_atom.atom_index, curr_atom.residue_index,
                curr_coord.x, curr_coord.y, curr_coord.z
            );
        }
    }
    return reconstructed_atoms;
}

/**
 * @brief Reconstruct the peptide in a reversed order
 * @param original_atoms
 * @param torsion_angles
 * @param atom_bond_angles
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> Nerf::reconstructWithReversed(
    std::vector<AtomCoordinate> original_atoms,
    std::vector<float> torsion_angles,
    std::vector<float> atom_bond_angles
) {

    std::reverse(original_atoms.begin(), original_atoms.end());
    std::reverse(torsion_angles.begin(), torsion_angles.end());
    std::reverse(atom_bond_angles.begin(), atom_bond_angles.end());

    // save three first atoms
    std::vector<AtomCoordinate> reconstructed_atoms = {
        original_atoms[0], original_atoms[1], original_atoms[2]
    };

    int total = original_atoms.size();
    for (int i = 0; i < (total - 3); i++) {
        const AtomCoordinate& curr_atom = original_atoms[i + 3];
        float3d prev_coord[3];
        extractCoordinates(prev_coord, reconstructed_atoms[i], reconstructed_atoms[i + 1],
                      reconstructed_atoms[i + 2] );
        std::string curr_bond_name = curr_atom.atom  + "_TO_" +original_atoms[i + 2].atom;
        float curr_bond_length = this->bond_lengths.at(curr_bond_name);
        float curr_bond_angle = atom_bond_angles[i + 1];
        float3d curr_coord = place_atom(prev_coord, curr_bond_length, curr_bond_angle,
            torsion_angles[i]);
        reconstructed_atoms.emplace_back(
            curr_atom.atom, curr_atom.residue, curr_atom.chain,
            curr_atom.atom_index, curr_atom.residue_index,
            curr_coord.x, curr_coord.y, curr_coord.z
        );
    }

    // reverse the output
    std::reverse(reconstructed_atoms.begin(), reconstructed_atoms.end());

    return reconstructed_atoms;
}


/**
 * @brief Write coordinate and bond info into csv file
 *
 * @param coord_list std::vector<AtomCoordinate>
 * @param output_path std::string, path for output
 * @param sep std::string, separation character (default: comma)
 */
void Nerf::writeInfoForChecking(
    const std::vector<AtomCoordinate>& coord_list,
    const std::string& output_path,
    const std::string sep
) {
    // setup file
    std::ofstream outfile;
    const int total = coord_list.size();
    // Things to write
    // CSV HEADER:
    //  AtomIndex,Atom,ResidueIndex,Residue,X,Y,Z,BondName,BondAngle,BondLength

    AtomCoordinate curr_atm;
    AtomCoordinate prev_atm, next_atm;
    std::string bond_name;
    std::string NA = "NA";
    float bond_angle, bond_length;
    std::string header = "AtomIndex,Atom,ResidueIndex,Residue,Chain,X,Y,Z,BondName,BondAngle,BondLength";

    outfile.open(output_path, std::ios::out);
    // Insert header
    outfile << header << "\n";

    for (int i = 0; i < total; i++) {
        //
        curr_atm = coord_list[i];

        if (i == 0) {
            outfile << curr_atm.atom_index << sep << curr_atm.atom << sep;
            outfile << curr_atm.residue_index << sep << curr_atm.residue << sep;
            outfile << curr_atm.chain << sep;
            outfile << curr_atm.coordinate.x << sep << curr_atm.coordinate.y << sep;
            outfile << curr_atm.coordinate.z << sep;
            outfile << NA << sep << NA << sep << NA << "\n";
        } else {
            // WARNING:
            // FIXME: FIX THIS --> curr_atm to next atm
            prev_atm = coord_list[i - 1];
            bond_name = prev_atm.atom + "_" + curr_atm.atom;
            bond_length = distance(prev_atm.coordinate, curr_atm.coordinate);
            // check error conditions
            // 01. same atom name with previous atom (ex: CA-CA bond)
            // --> SKIP PRINT
            if (prev_atm.atom == curr_atm.atom) {
                continue;
            }
            // print coordinate
            outfile << curr_atm.atom_index << sep << curr_atm.atom << sep;
            outfile << curr_atm.residue_index << sep << curr_atm.residue << sep;
            outfile << curr_atm.chain << sep;
            outfile << curr_atm.coordinate.x << sep << curr_atm.coordinate.y << sep;
            outfile << curr_atm.coordinate.z << sep;
            // 02. different chain
            if (prev_atm.chain != curr_atm.chain) {
                outfile << NA << sep << NA << NA << "\n";
                continue;
            }

            if (i == (total-1)) {
                outfile << bond_name << sep << NA << sep << bond_length << "\n";
            } else {
                next_atm = coord_list[i + 1];
                // handle error condition 01.
                if (next_atm.atom == curr_atm.atom) {
                    next_atm = coord_list[i + 2];
                }
                bond_angle = angle(
                    prev_atm.coordinate, curr_atm.coordinate, next_atm.coordinate
                );
                outfile << bond_name << sep << bond_angle << sep;
                outfile << bond_length << "\n";
            }
        }
    }
    outfile.close();
}


void Nerf::writeCoordinatesBinary(
    const std::vector<AtomCoordinate>& coord_list,
    const std::string& output_path
) {
    // setup file
    std::ofstream outfile;
    const int total = coord_list.size();
    // Things to write
    // CSV HEADER:
    //  AtomIndex,Atom,ResidueIndex,Residue,X,Y,Z,BondName,BondAngle,BondLength

    AtomCoordinate curr_atm;
    // std::string header = "AtomIndex,Atom,ResidueIndex,Residue,Chain,X,Y,Z";

    outfile.open(output_path, std::ios::out | std::ios::binary);

    for (int i = 0; i < total; i++) {
        curr_atm = coord_list[i];
        outfile.write((char*)&curr_atm.atom, 1);
        outfile.write((char*)&curr_atm.coordinate.x, 4);
        outfile.write((char*)&curr_atm.coordinate.y, 4);
        outfile.write((char*)&curr_atm.coordinate.z, 4);
    }
    outfile.close();
}



std::vector<float> Nerf::getBondAngles(const std::vector<AtomCoordinate>& original_atoms) {
    // define variables
    std::vector<float> output;
    output.reserve(original_atoms.size());
    const int total = original_atoms.size();
    for (int i = 1; i < (total - 1); i++) {
        const AtomCoordinate& curr_atm = original_atoms[i];
        const AtomCoordinate& prev_atm = original_atoms[i - 1];
        const AtomCoordinate& next_atm = original_atoms[i + 1];
        output.push_back(angle(prev_atm.coordinate, curr_atm.coordinate, next_atm.coordinate));
    }
    // ADD LAST ANGLE
    return output;
}

std::vector<float> Nerf::getBondLengths(const std::vector<AtomCoordinate>& original_atoms) {
    // define variables
    std::vector<float> output;
    output.reserve(original_atoms.size());
    const int total = original_atoms.size();
    for (int i = 0; i < (total - 1); i++) {
        //
        const AtomCoordinate& curr_atm = original_atoms[i];
        const AtomCoordinate& next_atm = original_atoms[i + 1];
        float bond_length = distance(curr_atm.coordinate, next_atm.coordinate);
        output.push_back(bond_length);
    }
    return output;
}

/**
 * @brief Identify breaks in chain and return the indices of break point
 *
 * @param original_atoms A vector of AtomCoordinate
 * @param cutoff Standard to define a break. Default value = 3
 * @return std::vector<int>
 */
std::vector<int> Nerf::identifyBreaks(const std::vector<AtomCoordinate>& original_atoms, float cutoff) {
    // define variables
    std::vector<int> output = {0};
    std::vector<float> bondLengths = this->getBondLengths(original_atoms);

    // No need to start with 0
    for (size_t i = 1; i < bondLengths.size(); i++) {
        // If length is bigger than the cutoff, add breakpoint
        if (bondLengths[i] > cutoff) {
            // bondLenght[i] = dist between original_atoms[i], [i+1]
            // append i+1 to the output
            if (i <= (original_atoms.size() - 3)) {
                output.push_back(i + 1);
                // std::cout << i + 1 << std::endl;
            }
        }
    }
    return output;
}

// std::vector<float> compute_positions(
//     std::vector<AtomCoordinate> coord_list
// ) {

//     std::vector<float> output {0.0, 0.0};
//     return output;
// }
