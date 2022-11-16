/**
 * File: sidechain.cpp
 * Project: foldcomp
 * Created: 2021-08-19 19:14:16
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Functions for handling side-chain atoms of amino acids
 * ---
 * Last Modified: 2022-09-13 15:15:23
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "sidechain.h"

#include "amino_acid.h"
#include "atom_coordinate.h"
#include "float3d.h"
#include "nerf.h"
#include "torsion_angle.h"
#include "utility.h"

#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <utility>


std::vector<AtomCoordinate> reconstructSideChain(
    std::vector<AtomCoordinate> backbone, AtomCoordinate /* c0 */
) {
    std::vector<AtomCoordinate> output;
    AminoAcid AA;
    std::map<std::string, AminoAcid> AAS = AA.AminoAcids();
    AtomCoordinate temp = backbone[0];
    AminoAcid currResidue = AAS[temp.residue];
    Nerf nerf;
    output = nerf.reconstructWithAAMaps(
        backbone, currResidue.sideChain, currResidue.torsionAngles,
        currResidue.bondLengths, currResidue.bondAngles
    );

    std::map<std::string, float> actualBondLengths = calculateBondLengths(backbone, currResidue);
    std::map<std::string, float> actualBondAngles = calculateBondAngles(backbone, currResidue);
    std::map<std::string, float> actualTorsionAngles = calculateTorsionAngles(backbone, currResidue);
    // compareMap(currResidue.bondLengths, actualBondLengths, "PeptideBuilder", "Actual");
    // compareMap(currResidue.bondAngles, actualBondAngles, "PeptideBuilder", "Actual");
    // compareMap(currResidue.torsionAngles, actualTorsionAngles, "PeptideBuilder", "Actual");
    return output;
}


std::vector<AtomCoordinate> reconstructSideChainFromCalculatedInfo(
    std::vector<AtomCoordinate> backbone
) {
    std::vector<AtomCoordinate> output;
    AminoAcid AA;
    std::map<std::string, AminoAcid> AAS = AA.AminoAcids();
    AtomCoordinate temp = backbone[0];
    AminoAcid currResidue = AAS[temp.residue];
    Nerf nerf;

    std::map<std::string, float> actualBondLengths = calculateBondLengths(backbone, currResidue);
    std::map<std::string, float> actualBondAngles = calculateBondAngles(backbone, currResidue);
    std::map<std::string, float> actualTorsionAngles = calculateTorsionAngles(backbone, currResidue);

    output = nerf.reconstructWithAAMaps(
        backbone, currResidue.sideChain, actualTorsionAngles,
        actualBondLengths, actualBondAngles
    );
    return output;
}

std::map<std::string, float3d> convertAtomCoordinateVectorToMap(
    const std::vector<AtomCoordinate>& input
) {
    std::map<std::string, float3d> output;
    for (const auto& ac : input) {
        output.emplace(ac.atom, ac.coordinate);
    }
    return output;
}

std::map<std::string, float> calculateBondLengths(
    std::vector<AtomCoordinate> originalAtoms,  AminoAcid AA
) {
    std::map<std::string, float> output;
    std::map<std::string, float3d> origAtomMap = convertAtomCoordinateVectorToMap(originalAtoms);
    // Calculating bond angle works for side chain - 2021-10-01 14:56:20
    for (const auto& atom : AA.sideChainAtoms) {
        const std::vector<std::string>& prev_atoms = AA.sideChain.at(atom);
        float3d prev_atom = origAtomMap.at(prev_atoms[2]);
        float3d curr_atom = origAtomMap.at(atom);
        std::string k = prev_atoms[2] + "_" + atom;
        float d = distance(prev_atom, curr_atom);
        output[k] = d;
    }
    return output;
}

std::map<std::string, float> calculateBondAngles(
    std::vector<AtomCoordinate> originalAtoms, AminoAcid AA
) {
    std::map<std::string, float> output;
    std::map<std::string, float3d> origAtomMap = convertAtomCoordinateVectorToMap(originalAtoms);
    for (const auto& atom : AA.sideChainAtoms) {
        if (atom == "N") { continue; }
        const std::vector<std::string>& prev_atoms = AA.sideChain[atom];
        float3d prev_atom1 = origAtomMap[prev_atoms[1]];
        float3d prev_atom2 = origAtomMap[prev_atoms[2]];
        float3d curr_atom = origAtomMap[atom];
        std::string k = prev_atoms[1] + "_" + prev_atoms[2] + "_" + atom;
        float a = angle(prev_atom1, prev_atom2, curr_atom);
        output[k] = a;
    }
    return output;
}

std::map<std::string, float> calculateTorsionAngles(
    std::vector<AtomCoordinate> originalAtoms, AminoAcid AA
) {
    std::map<std::string, float> output;
    std::map<std::string, float3d> origAtomMap = convertAtomCoordinateVectorToMap(originalAtoms);
    for (const auto& atom : AA.sideChainAtoms) {
        if (atom == "N" || atom == "CA") { continue; }
        const std::vector<std::string>& prev_atoms = AA.sideChain.at(atom);
        float3d prev_atom1 = origAtomMap.at(prev_atoms[0]);
        float3d prev_atom2 = origAtomMap.at(prev_atoms[1]);
        float3d prev_atom3 = origAtomMap.at(prev_atoms[2]);
        float3d curr_atom  = origAtomMap.at(atom);
        std::string k = prev_atoms[0] + "_" + prev_atoms[1] + "_" + prev_atoms[2] + "_" + atom;
        std::vector<float3d> coordinates = {prev_atom1, prev_atom2, prev_atom3, curr_atom};
        std::vector<float> vt = getTorsionFromXYZ(coordinates, 1);
        output[k] = vt[0];
    }
    return output;
}

float3d findFirstAtomCoords(const std::vector<AtomCoordinate>& atoms, std::string atom_name) {
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.atom == atom_name) {
            return curr_atm.coordinate;
        }
    }
    return {0, 0, 0};
}

std::vector<float> calculateTorsionAnglesInResidue(
    const std::vector<AtomCoordinate>& originalAtoms, const AminoAcid& AA
) {
    std::vector<float> output;
    output.reserve(AA.sideChainAtoms.size());
    for (const auto& atom : AA.sideChainAtoms) {
        if (atom == "N" || atom == "CA") { continue; }
        const std::vector<std::string>& prev_atom_names = AA.sideChain.at(atom);

        float3d prev_atom1 = findFirstAtomCoords(originalAtoms, prev_atom_names[0]);
        float3d prev_atom2 = findFirstAtomCoords(originalAtoms, prev_atom_names[1]);
        float3d prev_atom3 = findFirstAtomCoords(originalAtoms, prev_atom_names[2]);
        float3d curr_atom  = findFirstAtomCoords(originalAtoms, atom);

        std::vector<float3d> coordinates = {prev_atom1, prev_atom2, prev_atom3, curr_atom};
        std::vector<float> vt = getTorsionFromXYZ(coordinates, 1);
        output.push_back(vt[0]);
    }
    return output;
}

std::vector< std::vector<float> > calculateSideChainTorsionAnglesPerResidue(
    std::vector<AtomCoordinate>& originalAtoms, const std::map<std::string, AminoAcid>& AAmap
) {
    std::vector<std::vector<AtomCoordinate>> atomByResidue = splitAtomByResidue(originalAtoms);
    std::vector<std::vector<float>> output;
    output.reserve(atomByResidue.size());
    for (const auto& residue : atomByResidue) {
        output.emplace_back(calculateTorsionAnglesInResidue(residue, AAmap.at(residue[0].residue)));
    }
    return output;
}

void compareMap(
    std::map<std::string, float> m1, std::map<std::string, float> m2,
    std::string m1Name, std::string m2Name
) {
    float value;
    for (std::map<std::string, float>::iterator it = m1.begin(); it != m1.end(); ++it) {
        std::cout << it->first << " - " << m1Name << " : " << it->second << " / ";
        value = m2[it->first];
        std::cout << m2Name << " : " << value << std::endl;
    }
}

std::vector<AtomCoordinate> reconstructSideChain(
    std::vector<std::vector<float>> /* backboneCoordinates */,
    std::string /* residueName */
) {
    std::vector<AtomCoordinate> output;
    return output;
}

int calculateSideChainInfo(
    std::vector<AtomCoordinate>& peptideCoordinates,
    std::map<std::string, AminoAcid>& AAS
) {
    // Calculate side chain information of all amino acid in the peptide
    std::vector<AtomCoordinate> currResidue;
    size_t total = peptideCoordinates.size();
    AtomCoordinate* currAtom;
    int currResInd = peptideCoordinates[0].residue_index;
    std::string currRes = peptideCoordinates[0].residue;
    std::map<std::string, float> resCountMap = {
        {"ALA", 0.0f},{"ARG", 0.0f},{"ASN", 0.0f},{"ASP", 0.0f},{"CYS", 0.0f},
        {"GLN", 0.0f},{"GLU", 0.0f},{"GLY", 0.0f},{"HIS", 0.0f},{"ILE", 0.0f},
        {"LEU", 0.0f},{"LYS", 0.0f},{"MET", 0.0f},{"PHE", 0.0f},{"PRO", 0.0f},
        {"SER", 0.0f},{"THR", 0.0f},{"TRP", 0.0f},{"TYR", 0.0f},{"VAL", 0.0f}
    };
    std::map<std::string, float> currBondLengths;
    std::map<std::string, float> currBondAngles;
    std::map<std::string, float> currTorsionAngles;

    int success;
    // TODO: do something with success, remove next line
    (void)(success);
    for (size_t i = 0; i < total; i++) {
        currAtom = &(peptideCoordinates[i]);
        if ((currAtom->residue_index == currResInd) && (i != (total - 1))) {
            currResidue.push_back( *currAtom );
        // 2021-10-25 14:34:21 - NO PROBLEM UNTIL HERE
        } else if (currAtom->residue_index != currResInd) {
            currBondLengths = calculateBondLengths(currResidue, AAS[currRes]);
            currBondAngles = calculateBondAngles(currResidue, AAS[currRes]);
            currTorsionAngles = calculateTorsionAngles(currResidue, AAS[currRes]);
            // Add calculated values total
            success = addMap(&(AAS[currRes].bondLengths), &currBondLengths);
            success = addMap(&(AAS[currRes].bondAngles), &currBondAngles);
            success = addMap(&(AAS[currRes].torsionAngles), &currTorsionAngles);
            // Go to the next residue
            resCountMap[currRes] += 1;
            currResInd = currAtom->residue_index;
            currRes = currAtom->residue;
            currResidue.clear();
            currResidue.push_back( *currAtom );
        } else if (i == (total - 1)) {
            currResidue.push_back( *currAtom );
            currBondLengths = calculateBondLengths(currResidue, AAS[currRes]);
            currBondAngles = calculateBondAngles(currResidue, AAS[currRes]);
            currTorsionAngles = calculateTorsionAngles(currResidue, AAS[currRes]);
            // Add calculated values total
            success = addMap(&(AAS[currRes].bondLengths), &currBondLengths);
            success = addMap(&(AAS[currRes].bondAngles), &currBondAngles);
            success = addMap(&(AAS[currRes].torsionAngles), &currTorsionAngles);
            resCountMap[currRes] += 1;
            currResidue.clear();
        }
    }
    // Divide the sum with the count of residues to calculate average
    for (const auto& aa : AAS) {
        success = divideMapWithConst(&(AAS[aa.first].bondLengths), resCountMap[aa.first]);
        success = divideMapWithConst(&(AAS[aa.first].bondAngles), resCountMap[aa.first]);
        success = divideMapWithConst(&(AAS[aa.first].torsionAngles), resCountMap[aa.first]);
    }
    // Calculate the average value and save it to each map
    return 0;
}

int reconstructSideChainFromPeptide(
    std::vector<AtomCoordinate>& peptideCoordinates,
    std::map<std::string, AminoAcid>& AAS,
    std::vector<AtomCoordinate>& output
) {
    //2021-11-01 23:30:42 Using codepilot to fill in the empty function
    // TODO: Change & fill in the details
    // Start from the first residue
    std::vector<AtomCoordinate> currResidue;
    std::vector<AtomCoordinate> backbone;
    std::vector<AtomCoordinate> tempResidue;

    size_t total = peptideCoordinates.size();
    int currResInd = peptideCoordinates[0].residue_index;
    std::string currRes = peptideCoordinates[0].residue;
    int success = 0;
    Nerf nerf;

    std::map<std::string, float> currResTorsionAngles;

    // Need to calculate the bond lengths, angles and torsion angles from the peptide first
    // Calculate side chain information of all amino acid in the peptide
    // Iterate atom by atom
    for (size_t i = 0; i < total; i++) {
        // Handle backbone atoms
        if (peptideCoordinates[i].atom == "N" ||
            peptideCoordinates[i].atom == "CA" ||
            peptideCoordinates[i].atom == "C") {
            // Save to a vector
            backbone.push_back(peptideCoordinates[i]);
            // Backbone will be handled later
        }
        if (peptideCoordinates[i].atom == "OXT") {
            continue;
        }
        // Handle sidechain atoms
        // Calculate side chain information of all amino acid in the peptide
        if ((peptideCoordinates[i].residue_index == currResInd) && (i != (total - 1))) {
            // If the atom is in the same residue as the previous atom, add it to the current residue
            currResidue.push_back(peptideCoordinates[i]);
        } else if (peptideCoordinates[i].residue_index != currResInd) {
            // If the atom is in a different residue, calculate the side chain information
            currResTorsionAngles = calculateTorsionAngles(currResidue, AAS[currRes]);
            // 2021-11-09 23:30:19 DONE: NEED A RECONSTRUCTION METHOD WITH INPUT OF TORSION ANGLE.
            tempResidue = residueReconstruction(currResidue, AAS[currRes], currResTorsionAngles);
            for (auto atm : tempResidue) {
                output.push_back(atm);
            }
            // Go to the next residue
            currResidue.clear();
            currResidue.push_back(peptideCoordinates[i]);
            currResInd = peptideCoordinates[i].residue_index;
            currRes = peptideCoordinates[i].residue;
        } else if (i == (total - 1)) {
            currResidue.push_back(peptideCoordinates[i]);
            // If the atom is in a different residue, calculate the side chain information
            currResTorsionAngles = calculateTorsionAngles(currResidue, AAS[currRes]);
            // 2021-11-09 23:30:19 DONE: NEED A RECONSTRUCTION METHOD WITH INPUT OF TORSION ANGLE.
            tempResidue = residueReconstruction(currResidue, AAS[currRes], currResTorsionAngles);
            for (auto atm : tempResidue) {
                output.push_back(atm);
            }
        }
    }
    // Backbone reconstruction
    std::vector<float> backboneTorsionAngles;
    std::vector<float> backboneBondAngles;
    std::vector<float> backboneBondLengths;

    backboneTorsionAngles = getTorsionFromXYZ(backbone, 1);
    backboneBondAngles = nerf.getBondAngles(backbone);
    backboneBondLengths = nerf.getBondLengths(backbone);

    std::vector<int> backboneBreaks = nerf.identifyBreaks(backbone);

    std::vector<AtomCoordinate> reconstructedBackbone = nerf.reconstructWithBreaks(
        backbone, backboneTorsionAngles, backboneBondAngles, backboneBreaks
    );

    // Save torsion angles
    // success = saveTorsionAngles(peptideCoordinates, AAS);
    return success;
}

std::vector<AtomCoordinate> residueReconstruction(
    std::vector<AtomCoordinate>& currResidue,
    AminoAcid& currAA,
    std::map<std::string, float>& currResTorsionAngles
) {
    // Get bond angles and bond lengths from currAA
    const std::map<std::string, float>& bondAngles = currAA.bondAngles;
    const std::map<std::string, float>& bondLengths = currAA.bondLengths;

    Nerf nerf;
    std::vector<AtomCoordinate> output = nerf.reconstructWithAAMaps(
        currResidue, currAA.sideChain,
        currResTorsionAngles, bondLengths, bondAngles
    );
    return output;
}

// std::vector<AtomCoordinate> residueReconstruction(
//     std::vector<AtomCoordinate>& currResidue,
//     std::string residueName,
//     std::vector<float>& currResTorsionAngles
// ){

// }

int saveTorsionAngles(
    std::vector<AtomCoordinate>& /* peptideCoordinates */,
    std::map<std::string, AminoAcid>& AAS
) {
    // TODO: This is not a valid code. Need to fix it
    int success = 0;
    // Calculate torsion angles
    // for (auto aa : AAS) {
    //     success = calculateTorsionAngles(peptideCoordinates, aa.second);
    // }

    // Open a file
    std::ofstream outFile;
    outFile.open("torsionAngles.txt");
    // Write the header
    outFile << "Residue\t";
    for (auto aa : AAS) {
        outFile << aa.first << "\t";
    }
    outFile << "\n";
    // Write the data
    for (auto aa : AAS) {
        outFile << aa.first << "\t";
        for (auto angle : AAS[aa.first].torsionAngles) {
            outFile << angle.second << "\t";
        }
        outFile << "\n";
    }
    outFile.close();
    // Save the torsion angles to the file
    return success;
}

void checkEmptyAtomsInResidue(
    std::vector<AtomCoordinate>& Residue, AminoAcid& AA
){
    // Generated with copilot -> NEED TO BE CHANGED
    // Check if the residue is empty
    if (Residue.size() == 0) {
        std::cout << "ERROR: Empty residue" << std::endl;
        return;
    }
    // Check if the residue is a valid amino acid
    if (AA.sideChain.size() == 0) {
        std::cout << "ERROR: Invalid amino acid" << std::endl;
        return;
    }
    // Check if the residue is a valid amino acid
    if (AA.sideChain.size() != Residue.size()) {
        std::cout << "ERROR: Invalid amino acid" << std::endl;
        return;
    }
}

std::map<std::string, std::vector< std::vector<float> > > groupSideChainTorsionByResidue(
    std::vector< std::vector<float> >& sideChainTorsionAngles,
    std::vector<std::string>& residueNames,
    const std::map<std::string, AminoAcid>& AAS
) {
    std::map<std::string, std::vector< std::vector<float> > > output;
    // Initialize the map
    for (const auto& aa : AAS) {
        std::vector< std::vector<float> > temp;
        output[aa.first] = temp;
    }
    //
    std::string currResidue;
    for (size_t i = 0; i < sideChainTorsionAngles.size(); i++) {
        currResidue = residueNames[i];
        output[currResidue].push_back(sideChainTorsionAngles[i]);
    }
    return output;
}

std::vector<float> getSpecificTorsionAngle(
    std::map<std::string, std::vector< std::vector<float> > > sideChainTorsionMap,
    std::string& residueName, int index
) {
    const std::vector<std::vector<float>>& temp = sideChainTorsionMap[residueName];
    std::vector<float> output;
    output.resize(temp.size());
    for (size_t i = 0; i < temp.size(); i++) {
        output.push_back(temp[i][index]);
    }
    if (output.size() == 0) {
        output = std::vector<float>(1, 0);
    }
    return output;
}