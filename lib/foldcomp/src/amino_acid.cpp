/**
 * File: amino_acid.cpp
 * Project: foldcomp
 * Created: 2021-08-18 23:20:01
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Geometric information of amino acids.
 * ---
 * Last Modified: 2022-09-29 17:11:56
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

#include "amino_acid.h"

// 2021-08-18 23:20:21

int writeAminoAcidMapToFile(std::map<std::string, AminoAcid>& aa_map, std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    std::map<std::string, AminoAcid>::iterator it;
    std::map<std::string, float>::iterator inner_it;
    outfile << "aa_name,type,variable,value\n";

    for (it = aa_map.begin(); it != aa_map.end(); it++) {
        // Print bond lengths
        for (inner_it = it->second.bondLengths.begin(); inner_it != it->second.bondLengths.end(); inner_it++) {
            outfile << it->first << ",BL," << inner_it->first << "," << inner_it->second << "\n";
        }
        // Print bond angles
        for (inner_it = it->second.bondAngles.begin(); inner_it != it->second.bondAngles.end(); inner_it++) {
            outfile << it->first << ",BA," << inner_it->first << "," << inner_it->second << "\n";
        }
        // Print torsion angles
        for (inner_it = it->second.torsionAngles.begin(); inner_it != it->second.torsionAngles.end(); inner_it++) {
            outfile << it->first << ",TA," << inner_it->first << "," << inner_it->second << "\n";
        }
    }
    outfile.close();
    return 0;
}

std::vector<std::string> getAminoAcidList(void) {
    std::vector<std::string> aa_list = {
        // Sorted alphabetically
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    };
    return aa_list;
}