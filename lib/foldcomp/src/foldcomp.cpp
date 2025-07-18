/**
 * File: foldcomp.cpp
 * Project: foldcomp
 * Created: 2021-02-04 13:31:52
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This file contains main data structures for torsion angle compression and
 *     functions for handling them.
 * ---
 * Last Modified: 2024-08-08 19:43:31
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "foldcomp.h"

#include "amino_acid.h"
#include "sidechain.h"
#include "torsion_angle.h"
#include "utility.h"

#include <algorithm>
#include <bitset>
#include <utility>

// Changed at 2022-04-14 16:02:30
/**
 * @brief Convert a BackboneChain to a byte array (8 bytes)
 *
 * @param res a BackboneChain
 * @return char*
 */
int convertBackboneChainToBytes(BackboneChain& res, char* output) {
    int flag = 0;
    // 00-04: residue, 05-07: OMEGA first 3 bits
    output[0] = ((res.residue << 3) | ((res.omega & 0x07FF) >> 8));
    // 08-15: OMEGA last 8 bits
    output[1] = res.omega & 0x00FF;
    // 16-23: PSI first 8 bits
    output[2] = ((res.psi & 0x0FFF) >> 4);
    // 24-27: PSI last 4 bits, 28-31: PHI first 4 bits
    output[3] = ((res.psi & 0x000F) << 4) | ((res.phi & 0x0FFF) >> 8);
    // 32-39: PHI last 8 bits
    output[4] = res.phi & 0x00FF;
    // 40-47: CA_C_N
    output[5] = res.ca_c_n_angle;
    // 48-55: C_N_CA
    output[6] = res.c_n_ca_angle;
    // 56-63: N_CA_C
    output[7] = res.n_ca_c_angle;
    return flag;
}

/**
 * @brief Read a byte array and convert it to a BackboneChain.
 * This function is used for reading compressed residue from a file.
 * @param bytes a byte array (8 bytes) which encodes a BackboneChain
 * @return BackboneChain
 */
BackboneChain convertBytesToBackboneChain(char* bytes) {
    BackboneChain res;
    // 00-04: residue
    res.residue = ((bytes[0] & 0xF8) >> 3);
    // 05-07: OMEGA first 3 bits, 08-15: OMEGA last 8 bits
    res.omega = (unsigned int)(((bytes[0] & 0x0007) << 8) | (bytes[1] & 0x00FF));
    // 16-23: PSI first 8 bits, 24-27: PSI last 4 bits
    res.psi = (unsigned int)(((bytes[2] & 0x00FF) << 4) | (bytes[3] & 0x00FF) >> 4);
    // 28-31: PHI first 4 bits, 32-39: PHI last 8 bits
    res.phi = (unsigned int)(((bytes[3] & 0x000F) << 8) | (bytes[4] & 0x00FF));
    // 40-47: CA_C_N
    res.ca_c_n_angle = bytes[5];
    // 48-55: C_N_CA
    res.c_n_ca_angle = bytes[6];
    // 56-63: N_CA_C
    res.n_ca_c_angle = bytes[7];
    return res;
}

// NOTE:
BackboneChain newBackboneChain(
    char residue, unsigned int phi, unsigned int psi, unsigned int omega,
    unsigned int n_ca_c_angle, unsigned int ca_c_n_angle, unsigned int c_n_ca_angle
) {
    // unsigned int r = convertOneLetterCodeToInt(residue);
    BackboneChain res;
    res.residue = convertOneLetterCodeToInt(residue);
    res.ca_c_n_angle = ca_c_n_angle;
    res.c_n_ca_angle = c_n_ca_angle;
    res.n_ca_c_angle = n_ca_c_angle;
    res.psi = psi;
    res.omega = omega;
    res.phi = phi;
    return res;
}

BackboneChain newBackboneChain(
    unsigned int bResidue, unsigned int phi, unsigned int psi, unsigned int omega,
    unsigned int n_ca_c_angle, unsigned int ca_c_n_angle, unsigned int c_n_ca_angle
) {
    BackboneChain res;
    res.residue = bResidue;
    res.ca_c_n_angle = ca_c_n_angle;
    res.c_n_ca_angle = c_n_ca_angle;
    res.n_ca_c_angle = n_ca_c_angle;
    res.psi = psi;
    res.omega = omega;
    res.phi = phi;
    return res;
}


// WARNING:
/**
 * @brief Convert a BackboneChain to DecompressedBackboneChain
 *
 * @description As the angles are short-encoded in the compressed format,
 * this function converts the short-encoded angles to the float.
 * @param bb
 * @param header
 * @return DecompressedBackboneChain
 */
DecompressedBackboneChain decompressBackboneChain(
    const BackboneChain& bb, const CompressedFileHeader& header
) {
    DecompressedBackboneChain output;
    output.residue = convertIntToOneLetterCode(bb.residue);
    output.phi = _continuize(bb.phi, header.mins[0], header.cont_fs[0]);
    output.psi = _continuize(bb.psi, header.mins[1], header.cont_fs[1]);
    output.omega = _continuize(bb.omega, header.mins[2], header.cont_fs[2]);
    output.n_ca_c_angle = _continuize(bb.n_ca_c_angle, header.mins[3], header.cont_fs[3]);
    output.ca_c_n_angle = _continuize(bb.ca_c_n_angle, header.mins[4], header.cont_fs[4]);
    output.c_n_ca_angle = _continuize(bb.c_n_ca_angle, header.mins[5], header.cont_fs[5]);
    return output;
}


/**
 * @brief Convert a vectorof BackboneChain to DecompressedBackboneChain vector
 *
 * @param bbv
 * @param header
 * @return std::vector<DecompressedBackboneChain>
 */
std::vector<DecompressedBackboneChain> decompressBackboneChain(
    const std::vector<BackboneChain>& bbv, const CompressedFileHeader& header
) {
    std::vector<DecompressedBackboneChain> output;
    output.reserve(bbv.size());
    for (const auto& bb : bbv) {
        output.push_back(decompressBackboneChain(bb, header));
    }
    return output;
}

float _continuize(unsigned int input, float min, float cont_f) {
    float output = min + ((float)input * cont_f);
    return output;
}

/**
 * @brief Reconstruct backbone atoms from compressed backbone info
 *
 * @param prev_atoms std::vector<AtomCoordinates> of 3 previous atoms
 * @param backbone std::vector<BackboneChain>
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> reconstructBackboneAtoms(
    const std::vector<AtomCoordinate>& prevAtoms,
    const std::vector<BackboneChain>& backbone,
    CompressedFileHeader& header
) {
    Nerf nerf;
    // Save first three atoms
    std::vector<AtomCoordinate> reconstructedAtoms = {
        prevAtoms[0], prevAtoms[1], prevAtoms[2]
    };
    int total = backbone.size();

    std::vector<DecompressedBackboneChain> deBackbone = decompressBackboneChain(backbone, header);
    int currAtomIndex = prevAtoms[2].atom_index + 1;
    int currResidueIndex = prevAtoms[2].residue_index + 1;

    // Iterate through backbone
    // Should put N, CA, C in this loop
    for (int i = 0; i < (total - 1); i++) {
        const AtomCoordinate& prevAtom1 = reconstructedAtoms[i*3];
        const AtomCoordinate& prevAtom2 = reconstructedAtoms[i*3 + 1];
        const AtomCoordinate& prevAtom3 = reconstructedAtoms[i*3 + 2];
        float3d prevCoords[3];
        extractCoordinates(prevCoords, prevAtom1, prevAtom2, prevAtom3);
        // Convert char (deBackbone[i].residue) to string (currResidue)
        std::string currResidue = getThreeLetterCode(deBackbone[i + 1].residue);
        std::string currChain = prevAtom1.chain;

        // Place N
        float3d currNCoord = nerf.place_atom(
            prevCoords, C_TO_N_DIST, deBackbone[i].ca_c_n_angle, deBackbone[i].psi
        );
        // Place CA
        prevCoords[0] = prevCoords[1];
        prevCoords[1] = prevCoords[2];
        prevCoords[2] = currNCoord;
        float3d currCACoord;
        if (deBackbone[i].residue != 'P') {
            currCACoord = nerf.place_atom(
                prevCoords, N_TO_CA_DIST, deBackbone[i].c_n_ca_angle, deBackbone[i].omega
            );
        } else {
            currCACoord = nerf.place_atom(
                prevCoords, PRO_N_TO_CA_DIST, deBackbone[i].c_n_ca_angle, deBackbone[i].omega
            );
        }

        // Place C
        currAtomIndex++;
        prevCoords[0] = prevCoords[1];
        prevCoords[1] = prevCoords[2];
        prevCoords[2] = currCACoord;
        float3d currCCoord = nerf.place_atom(
            prevCoords, CA_TO_C_DIST, deBackbone[i].n_ca_c_angle, deBackbone[i].phi
        );

        reconstructedAtoms.emplace_back(
            "N", currResidue, currChain,
            currAtomIndex, currResidueIndex,
            currNCoord.x, currNCoord.y, currNCoord.z
        );
        currAtomIndex++;
        reconstructedAtoms.emplace_back(
            "CA", currResidue, currChain,
            currAtomIndex, currResidueIndex,
            currCACoord.x, currCACoord.y, currCACoord.z
        );
        currAtomIndex++;
        reconstructedAtoms.emplace_back(
            "C", currResidue, currChain,
            currAtomIndex, currResidueIndex,
            currCCoord.x, currCCoord.y, currCCoord.z
        );
        // Increment Residue Index
        currResidueIndex++;
        currAtomIndex++;

    }
    return reconstructedAtoms;
}

int reconstructBackboneReverse(
    std::vector<AtomCoordinate>& atom, std::vector< std::vector<float> >& lastCoords,
    std::vector<float>& torsion_angles, Nerf& nerf
) {
    std::vector<AtomCoordinate> atomBack = atom;
    // Last atoms
    atomBack[atomBack.size() - 3].coordinate.x = lastCoords[0][0];
    atomBack[atomBack.size() - 3].coordinate.y = lastCoords[0][1];
    atomBack[atomBack.size() - 3].coordinate.z = lastCoords[0][2];
    atomBack[atomBack.size() - 2].coordinate.x = lastCoords[1][0];
    atomBack[atomBack.size() - 2].coordinate.y = lastCoords[1][1];
    atomBack[atomBack.size() - 2].coordinate.z = lastCoords[1][2];
    atomBack[atomBack.size() - 1].coordinate.x = lastCoords[2][0];
    atomBack[atomBack.size() - 1].coordinate.y = lastCoords[2][1];
    atomBack[atomBack.size() - 1].coordinate.z = lastCoords[2][2];

    std::vector<float> bond_angles = nerf.getBondAngles(atom);

    std::vector<AtomCoordinate> atomBackward = nerf.reconstructWithReversed(
        atomBack, torsion_angles, bond_angles
    );

    atomBack = weightedAverage(atom, atomBackward);
    atom = atomBack;
    return 0;
}


int discretizeSideChainTorsionAngles(
    std::vector< std::vector<float> >& torsionPerResidue,
    std::vector<std::string>& residueNames,
    const std::map<std::string, AminoAcid>& AAS,
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap,
    std::vector<unsigned int>& output
) {
    // Declare
    int success = 0;
    std::string currResidue;
    int currResidueTorsionNum;
    float min, cont_f;
    std::vector<float> currTorsion;
    std::vector<unsigned int> currTorsionDiscretized;
    float* min_arr;
    float* cont_f_arr;

    std::map<std::string, std::vector< std::vector<float> > > sideChainTorsionMap;
    sideChainTorsionMap = groupSideChainTorsionByResidue(torsionPerResidue, residueNames, AAS);

    // Fill in Discretizer map
    for (const auto& sc : sideChainTorsionMap) {
        currResidue = sc.first;
        currResidueTorsionNum = getSideChainTorsionNum(currResidue);
        min_arr = getMinPointerFromSideChainDiscretizers(currResidue, scDiscretizers);
        cont_f_arr = getContFFromSideChainDiscretizers(currResidue, scDiscretizers);
        for (int i = 0; i < currResidueTorsionNum; i++) {
            // Get current torsion angle vector
            currTorsion = getSpecificTorsionAngle(sideChainTorsionMap, currResidue, i);
            // Discretize the torsion angles
            // Discretizer currDiscretizer = Discretizer(currTorsion, pow(2, NUM_BITS_SIDECHAIN) - 1);
            // TODO: TESTING FIXEDANGLEDISCRETIZER
            FixedAngleDiscretizer currDiscretizer = FixedAngleDiscretizer(pow(2, NUM_BITS_TEMP) - 1);

            // Save min and cont_f to SideChainDiscretizers
            min = currDiscretizer.min;
            cont_f = currDiscretizer.cont_f;
            min_arr[i] = min;
            cont_f_arr[i] = cont_f;
            // Save Discretizer to map
            scDiscretizersMap[currResidue][i] = currDiscretizer;
        }
    }

    Discretizer torsionDisc;
    unsigned int torsionDiscretized;
    // Discretize the torsion angles and make the result into a flattend vector
    for (size_t i = 0; i < torsionPerResidue.size(); i++) {
        currResidue = residueNames[i];
        currResidueTorsionNum = getSideChainTorsionNum(currResidue);
        for (int j = 0; j < currResidueTorsionNum; j++) {
            torsionDisc = scDiscretizersMap[currResidue][j];
            torsionDiscretized = torsionDisc.discretize(torsionPerResidue[i][j]);
            // Append to flattened vector
            output.push_back(torsionDiscretized);
        }
    }

    return success;
}

int continuizeSideChainTorsionAngles(
    std::vector<unsigned int>& torsionDiscretized,
    std::vector<std::string>& residueNames,
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap,
    std::vector< std::vector<float> >& output
) {

    scDiscretizersMap = initializeSideChainDiscMap();
    int success = fillSideChainDiscretizerMap(scDiscretizers, scDiscretizersMap);
    std::vector< std::vector<float> > torsionPerResidue;
    std::vector<float> currTorsionVector;
    FixedAngleDiscretizer currDiscretizer(pow(2, NUM_BITS_TEMP) - 1);
    // Iterate
    int currIndex = 0;
    for (size_t i = 0; i < residueNames.size(); i++) {
        std::string currResidue = residueNames[i];
        int currResidueTorsionNum = getSideChainTorsionNum(currResidue);
        currTorsionVector.clear();
        currTorsionVector.resize(currResidueTorsionNum);
        for (int j = 0; j < currResidueTorsionNum; j++) {
            // Get current torsion angle vector
            //currTorsion = scDiscretizersMap[currResidue][j].continuize(torsionDiscretized[currIndex]);
            float currTorsion = currDiscretizer.continuize(torsionDiscretized[currIndex]);
            currTorsionVector[j] = currTorsion;
            currIndex++;
        }
        torsionPerResidue.push_back(currTorsionVector);
    }
    output = torsionPerResidue;
    return success;
}

int fillSideChainDiscretizerMap(
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap
) {
    // Declare
    int success = 0;
    std::string currResidue;
    int currResidueTorsionNum;
    float min, cont_f;
    // Iterate through map
    for (const auto& sc : scDiscretizersMap) {
        currResidue = sc.first;
        currResidueTorsionNum = getSideChainTorsionNum(currResidue);
        for (int i = 0; i < currResidueTorsionNum; i++) {
            // Get current torsion angle vector
            min = getMinPointerFromSideChainDiscretizers(currResidue, scDiscretizers)[i];
            cont_f = getContFFromSideChainDiscretizers(currResidue, scDiscretizers)[i];
            // Save min and cont_f to SideChainDiscretizers
            scDiscretizersMap[currResidue][i].min = min;
            scDiscretizersMap[currResidue][i].cont_f = cont_f;
        }
    }
    return success;
}

int Foldcomp::_discretizeSideChainTorsionAngles(
    std::vector< std::vector<float> >& input,
    std::vector<unsigned int>& output
) {
    // Declare
    int success = 0;
    success = discretizeSideChainTorsionAngles(
        input, this->residueThreeLetter, this->AAS,
        this->sideChainDisc, this->sideChainDiscMap, output
    );
    this->header.nSideChainTorsion = output.size();
    this->nSideChainTorsion = output.size();
    return success;
}

int Foldcomp::_continuizeSideChainTorsionAngles(
    std::vector<unsigned int>& input, std::vector< std::vector<float> >& output
) {
    // Declare
    int success = 0;
    success = continuizeSideChainTorsionAngles(
        input, this->residueThreeLetter, this->sideChainDisc, this->sideChainDiscMap, output
    );
    return success;
}

/**
 * @brief [TEMP] Print the bit representation of the variables in the
 * compressed residue
 *
 * @param res
 */
void printCompressedResidue(BackboneChain& res) {
    // Print the variables in the header
    std::cout << "SIZE: " << sizeof(res) << std::endl;
    std::cout << "residue: " << res.residue << std::endl;
    std::cout << "phi: " << res.phi << std::endl;
    std::cout << "psi: " << res.psi << std::endl;
    std::cout << "omega: " << res.omega << std::endl;
    std::cout << "n_ca_c_angle: " << res.n_ca_c_angle << std::endl;
    std::cout << "ca_c_n_angle: " << res.ca_c_n_angle << std::endl;
    std::cout << "c_n_ca_angle: " << res.c_n_ca_angle << std::endl;
    std::cout << "CONVERTED BYTE ARRAY: ";
    std::bitset<8> bits;
    char* byteArray = new char[8];
    // int flag = convertBackboneChainToBytes(res, byteArray);
    for (int i = 0; i < 8; i++) {
        bits = byteArray[i];
        std::cout << bits << " ";
    }
    delete[] byteArray;
    std::cout << std::endl;
}

int Foldcomp::preprocess(const tcb::span<AtomCoordinate>& atoms) {
    int success = 0;
    this->compressedBackBone.resize(atoms.size());
    this->compressedSideChain.resize(atoms.size());

    // Discretize
    // Extract backbone
    this->backbone = filterBackbone(atoms);

    // this->title = this->strTitle.c_str();
    this->lenTitle = this->strTitle.size();

    this->nResidue = this->backbone.size() / 3;
    this->nBackbone = this->backbone.size();
    this->nAtom = atoms.size();
    this->idxResidue = atoms[0].residue_index;
    this->idxAtom = atoms[0].atom_index;
    this->chain = atoms[0].chain.c_str()[0];
    this->firstResidue = getOneLetterCode(atoms[0].residue);
    this->lastResidue = getOneLetterCode(atoms[atoms.size() - 1].residue);

    // Anchor atoms
    this->_setAnchor(atoms);

    if (atoms[atoms.size() - 1].atom == "OXT") {
        this->hasOXT = 1;
        this->OXT = atoms[atoms.size() - 1];
        this->OXT_coords = atoms[atoms.size() - 1].coordinate;
    } else {
        this->hasOXT = 0;
        this->OXT = AtomCoordinate();
        this->OXT_coords = {0.0, 0.0, 0.0};
    }

    std::vector<float> backboneTorsion = getTorsionFromXYZ(this->backbone, 1);
    this->backboneTorsionAngles = backboneTorsion;
    // Split backbone into phi, psi, omega
    // Calculate phi, psi, omega
    for (size_t i = 0; i < backboneTorsion.size(); i += 3) {
        this->psi.push_back(backboneTorsion[i]);
        this->omega.push_back(backboneTorsion[i + 1]);
        this->phi.push_back(backboneTorsion[i + 2]);
    }

    std::vector<float> backboneBondAngles = this->nerf.getBondAngles(backbone);
    this->backboneBondAngles = backboneBondAngles;
    // Split bond angles into three parts
    for (size_t i = 1; i < backboneBondAngles.size(); i++) {
        if (i % 3 == 0) {
            this->n_ca_c_angle.push_back(backboneBondAngles[i]);
        } else if (i % 3 == 1) {
            this->ca_c_n_angle.push_back(backboneBondAngles[i]);
        } else {
            this->c_n_ca_angle.push_back(backboneBondAngles[i]);
        }
    }

    // Discretize
    this->phiDisc = Discretizer(this->phi, pow(2, NUM_BITS_PHI_PSI) - 1);
    this->phiDiscretized = this->phiDisc.discretize(this->phi);
    this->omegaDisc = Discretizer(this->omega, pow(2, NUM_BITS_OMEGA) - 1);
    this->omegaDiscretized = this->omegaDisc.discretize(this->omega);
    this->psiDisc = Discretizer(this->psi, pow(2, NUM_BITS_PHI_PSI) - 1);
    this->psiDiscretized = this->psiDisc.discretize(this->psi);
    this->n_ca_c_angleDisc = Discretizer(this->n_ca_c_angle, pow(2, NUM_BITS_BOND) - 1);
    this->n_ca_c_angleDiscretized = this->n_ca_c_angleDisc.discretize(this->n_ca_c_angle);
    this->ca_c_n_angleDisc = Discretizer(this->ca_c_n_angle, pow(2, NUM_BITS_BOND) - 1);
    this->ca_c_n_angleDiscretized = this->ca_c_n_angleDisc.discretize(this->ca_c_n_angle);
    this->c_n_ca_angleDisc = Discretizer(this->c_n_ca_angle, pow(2, NUM_BITS_BOND) - 1);
    this->c_n_ca_angleDiscretized = this->c_n_ca_angleDisc.discretize(this->c_n_ca_angle);

    // Set residue names
    this->residueThreeLetter = getResidueNameVector(atoms);

    // Set Discretizer for side chain
    this->sideChainDiscMap = initializeSideChainDiscMap();
    // Calculate side chain info

    this->sideChainAnglesPerResidue = calculateSideChainTorsionAnglesPerResidue(atoms, this->AAS);

    // Discretize side chain
    //this->_discretizeSideChainTorsionAngles(this->sideChainAnglesPerResidue, this->sideChainAnglesDiscretized);
    FixedAngleDiscretizer sideChainDiscretizer = FixedAngleDiscretizer(pow(2, NUM_BITS_TEMP) - 1);
    for (size_t i = 0; i < this->sideChainAnglesPerResidue.size(); i++) {
        for (size_t j = 0; j < this->sideChainAnglesPerResidue[i].size(); j++) {
            unsigned int temp = sideChainDiscretizer.discretize(this->sideChainAnglesPerResidue[i][j]);
            this->sideChainAnglesDiscretized.push_back(temp);
        }
    }
    this->nSideChainTorsion = this->sideChainAnglesDiscretized.size();

    // Get tempFactors
    // 2022-08-31 16:28:30 - Changed to save one tempFactor per residue
    for (size_t i = 0; i < atoms.size(); i++) {
        if (atoms[i].atom == "CA") {
            this->tempFactors.push_back(atoms[i].tempFactor);
        }
    }
    // Discretize
    this->tempFactorsDisc = Discretizer(this->tempFactors, pow(2, NUM_BITS_TEMP) - 1);
    this->tempFactorsDiscretized = this->tempFactorsDisc.discretize(this->tempFactors);

    // Get header
    this->header = this->get_header();

    // Mark as processed
    this->isPreprocessed = true;

    return success;
}


std::vector<BackboneChain> Foldcomp::compress(const tcb::span<AtomCoordinate>& atoms) {
    std::vector<BackboneChain> output;
    // TODO: convert the atom coordinate vector into a vector of compressed residue
    // CURRENT VERSION - 2022-01-10 15:34:21
    // IGNORE BREAKS
    if (!this->isPreprocessed) {
        this->preprocess(atoms);
    }
    this->prevAtoms = {atoms[0], atoms[1], atoms[2]};

    AtomCoordinate currN;
    char currResCode;

    this->lastAtoms = {this->backbone[this->backbone.size() - 1],
                       this->backbone[this->backbone.size() - 2],
                       this->backbone[this->backbone.size() - 3]};

    // Need to extract backbone atoms in a separate vector

    BackboneChain res;
    for (int i = 0; i < (this->nResidue - 1); i++) {
        currN = this->backbone[i * 3];
        currResCode = getOneLetterCode(currN.residue);
        this->residues.push_back(currResCode);
        res.residue = convertOneLetterCodeToInt(currResCode);
        res.psi = this->psiDiscretized[i];
        res.omega = this->omegaDiscretized[i];
        res.phi = this->phiDiscretized[i];
        res.n_ca_c_angle = this->n_ca_c_angleDiscretized[i];
        res.ca_c_n_angle = this->ca_c_n_angleDiscretized[i];
        res.c_n_ca_angle = this->c_n_ca_angleDiscretized[i];
        output.push_back(res);
    }
    currN = this->backbone[(this->nResidue - 1) * 3];
    currResCode = getOneLetterCode(currN.residue);
    this->residues.push_back(currResCode);
    res.residue = convertOneLetterCodeToInt(currResCode);
    res.psi = 0; res.omega = 0; res.phi = 0;
    res.n_ca_c_angle = 0; res.ca_c_n_angle = 0; res.c_n_ca_angle = 0;
    output.push_back(res);
    this->compressedBackBone = output;

    this->isCompressed = true;
    return output;
}

int _restoreResidueNames(
    std::vector<BackboneChain>& compressedBackbone,
    CompressedFileHeader& /* header */,
    std::vector<char>& residueOneLetter,
    std::vector<std::string>& residueThreeLetter
) {
    int success = 0;
    std::string threeLetterCode;
    char oneLetterCode;
    residueThreeLetter.clear();
    for (size_t i = 0; i < compressedBackbone.size(); i++) {
        oneLetterCode = convertIntToOneLetterCode(compressedBackbone[i].residue);
        threeLetterCode = convertIntToThreeLetterCode(compressedBackbone[i].residue);
        residueOneLetter.push_back(oneLetterCode);
        residueThreeLetter.push_back(threeLetterCode);
    }
    return success;
}

int Foldcomp::_restoreDiscretizer(int angleType) {
    int success = 0;
    std::vector<unsigned int> temp;
    temp.reserve(this->compressedBackBone.size());
    for (const auto& cb : this->compressedBackBone) {
        switch (angleType) {
        case 0: // Phi
            temp.push_back(cb.phi);
            break;
        case 1:
            temp.push_back(cb.psi);
            break;
        case 2:
            temp.push_back(cb.omega);
            break;
        case 3:
            temp.push_back(cb.n_ca_c_angle);
            break;
        case 4:
            temp.push_back(cb.ca_c_n_angle);
            break;
        case 5:
            temp.push_back(cb.c_n_ca_angle);
            break;
        default:
            break;
        }
    }
    // Set
    switch (angleType) {
    case 0: // Phi
        this->phiDiscretized = temp;
        this->phiDisc.min = this->header.mins[0];
        this->phiDisc.cont_f = this->header.cont_fs[0];
        break;
    case 1: // Psi
        this->psiDiscretized = temp;
        this->psiDisc.min = this->header.mins[1];
        this->psiDisc.cont_f = this->header.cont_fs[1];
        break;
    case 2: // Omega
        this->omegaDiscretized = temp;
        this->omegaDisc.min = this->header.mins[2];
        this->omegaDisc.cont_f = this->header.cont_fs[2];
        break;
    case 3: // N-CA-C
        this->n_ca_c_angleDiscretized = temp;
        this->n_ca_c_angleDisc.min = this->header.mins[3];
        this->n_ca_c_angleDisc.cont_f = this->header.cont_fs[3];
        break;
    case 4: // CA-C-N
        this->ca_c_n_angleDiscretized = temp;
        this->ca_c_n_angleDisc.min = this->header.mins[4];
        this->ca_c_n_angleDisc.cont_f = this->header.cont_fs[4];
        break;
    case 5: // C-N-CA
        this->c_n_ca_angleDiscretized = temp;
        this->c_n_ca_angleDisc.min = this->header.mins[5];
        this->c_n_ca_angleDisc.cont_f = this->header.cont_fs[5];
        break;
    default:
        break;
    }

    return success;
}

// TODO: Write a function to restore prev atoms (AtomCoordinate) from coordinates (float)
// 2022-02-17 23:29:22
/**
 * @brief Restore previous atoms from the header & coordinates
 *
 * @param coords
 * @return int
 */
int Foldcomp::_restoreAtomCoordinate(float* coords) {
    int success = 0;
    std::string firstResidue = getThreeLetterCode(this->header.firstResidue);
    // TODO: Fix the chain alphabet according to the real data
    // 2022-03-04 14:30:58
    // IMPORTANT: WARNING: TODO: Atom indexing is not correct right now
    // SIDE CHAIN ATOMS SHOULD GE CONSIDERED WHEN INDEXING
    // convert char to string
    std::string chain = std::string(1, this->header.chain);
    AtomCoordinate prevN = AtomCoordinate(
        "N", firstResidue, chain, this->header.idxAtom, this->header.idxResidue,
        coords[0], coords[1], coords[2]
    );
    AtomCoordinate prevCA = AtomCoordinate(
        "CA", firstResidue, chain, this->header.idxAtom + 1, this->header.idxResidue,
        coords[3], coords[4], coords[5]
    );
    AtomCoordinate prevC = AtomCoordinate(
        "C", firstResidue, chain, this->header.idxAtom + 2, this->header.idxResidue,
        coords[6], coords[7], coords[8]
    );

    // Check if this->prevAtoms is empty
    if (this->prevAtoms.size() == 0) {
        this->prevAtoms.push_back(prevN);
        this->prevAtoms.push_back(prevCA);
        this->prevAtoms.push_back(prevC);
    } else {
        // If not empty, it will be updated
        this->prevAtoms[0] = prevN;
        this->prevAtoms[1] = prevCA;
        this->prevAtoms[2] = prevC;
    }

    return success;
}

int Foldcomp::_getAnchorNum(int threshold) {
    int nAnchor = 0;
    nAnchor = this->nResidue / threshold;
    return nAnchor;
}

void Foldcomp::_setAnchor(const tcb::span<AtomCoordinate>& atomCoordinates) {
    this->nInnerAnchor = this->_getAnchorNum(this->anchorThreshold);
    this->nAllAnchor = this->nInnerAnchor + 2; // Start and end
    // Set the anchor points - residue index
    this->anchorIndices.clear();
    int interval = this->nResidue / (this->nAllAnchor - 1);
    for (int i = 0; i < this->nAllAnchor - 1; i++) {
        this->anchorIndices.push_back(i * interval);
    }
    this->anchorIndices.push_back(this->nResidue - 1);
    //
    std::vector<int> anchorResidueIndices;
    for (size_t i = 0; i < this->anchorIndices.size(); i++) {
        anchorResidueIndices.push_back(this->anchorIndices[i] + this->idxResidue);
    }
    this->anchorAtoms = getAtomsWithResidueIndex(atomCoordinates, anchorResidueIndices);
}

std::vector<float> Foldcomp::checkTorsionReconstruction() {
    // Continuize torsion angles
    this->phi = this->phiDisc.continuize(this->phiDiscretized);
    this->psi = this->psiDisc.continuize(this->psiDiscretized);
    this->omega = this->omegaDisc.continuize(this->omegaDiscretized);
    // Append psi, omega, phi to torsion angles
    std::vector<float> output;
    output.reserve(this->phi.size());
    for (size_t i = 0; i < this->phi.size(); i++) {
        output.push_back(this->psi[i]);
        output.push_back(this->omega[i]);
        output.push_back(this->phi[i]);
    }
    return output;
}

int Foldcomp::decompress(std::vector<AtomCoordinate>& atom) {
    // 2022-11-15 11:47:49 - Removed defining new vectors for TAs & BAs
    int success;

    // Continuize torsion angles
    this->phi = this->phiDisc.continuize(this->phiDiscretized);
    this->psi = this->psiDisc.continuize(this->psiDiscretized);
    this->omega = this->omegaDisc.continuize(this->omegaDiscretized);

    // Append psi, omega, phi to torsion angles
    for (size_t i = 0; i < (this->phi.size() - 1); i++) {
        this->backboneTorsionAngles.push_back(this->psi[i]);
        this->backboneTorsionAngles.push_back(this->omega[i]);
        this->backboneTorsionAngles.push_back(this->phi[i]);
    }

    // Continuize bond angles
    this->n_ca_c_angle = this->n_ca_c_angleDisc.continuize(this->n_ca_c_angleDiscretized);
    this->ca_c_n_angle = this->ca_c_n_angleDisc.continuize(this->ca_c_n_angleDiscretized);
    this->c_n_ca_angle = this->c_n_ca_angleDisc.continuize(this->c_n_ca_angleDiscretized);
    // Append n_ca_c_angle, ca_c_n_angle, c_n_ca_angle to bond angles
    for (size_t i = 0; i < this->n_ca_c_angle.size(); i++) {
        this->backboneBondAngles.push_back(this->ca_c_n_angle[i]);
        this->backboneBondAngles.push_back(this->c_n_ca_angle[i]);
        this->backboneBondAngles.push_back(this->n_ca_c_angle[i]);
    }

    // Get the residue vector
    // TODO: EXTRACT RESIDUE NAMES FROM COMPRESSED BACKBONE
    success = _restoreResidueNames(this->compressedBackBone, this->header, this->residues, this->residueThreeLetter);

    //TODO: FILL IN THIS FUNCTION
    // nerf.reconstruct();
    std::vector<AtomCoordinate> prevForAnchor;
    std::vector<AtomCoordinate> atomByAnchor;
    for (int i = 0; i < this->nAllAnchor - 1; i++) {
        if (i == 0) {
            prevForAnchor = this->prevAtoms;
        }
        // std::vector<int>   sub(&data[100000],&data[101000]);
        // MANUALLY CHECK THAT THE INDICES ARE WITHIN BOUNDS
        int maxIndex = (int)this->compressedBackBone.size() - 1;
        size_t firstIndex = std::min(this->anchorIndices[i], maxIndex);
        size_t lastIndex = std::min(this->anchorIndices[i + 1] + 1, maxIndex);
        std::vector<BackboneChain> subBackbone(
            &this->compressedBackBone[firstIndex],
            &this->compressedBackBone[lastIndex]
        );
        // LAST RESIDUE SHOULD BE INCLUDED
        if (i == (this->nAllAnchor - 2)) {
            subBackbone.push_back(this->compressedBackBone.back());
        }
        atomByAnchor = reconstructBackboneAtoms(prevForAnchor, subBackbone, this->header);

        // Subset torsion_angles
        maxIndex = (int)this->backboneTorsionAngles.size() - 1;
        firstIndex = std::min(this->anchorIndices[i] * 3, maxIndex);
        lastIndex = std::min(this->anchorIndices[i + 1] * 3, maxIndex);
        std::vector<float> subTorsionAngles(
            &this->backboneTorsionAngles[firstIndex], &this->backboneTorsionAngles[lastIndex]
        );
        // LAST ANGLE SHOULD BE APPENDED
        if (i == (this->nAllAnchor - 2)) {
            subTorsionAngles.push_back(this->backboneTorsionAngles.back());
        }

        success = reconstructBackboneReverse(
            atomByAnchor, this->anchorCoordinates[i], subTorsionAngles, this->nerf
        );
        // Append atomByAnchor to atom
        if (i != this->nAllAnchor - 2) {
            atom.insert(atom.end(), atomByAnchor.begin(), atomByAnchor.end() - 3);
        } else {
            atom.insert(atom.end(), atomByAnchor.begin(), atomByAnchor.end());
        }
        // Update prevForAnchor - last 3 atoms of atomByAnchor
        prevForAnchor = std::vector<AtomCoordinate>(
            atomByAnchor.end() - 3, atomByAnchor.end()
        );
    }

    // Reconstruct sidechain
    std::vector<std::vector<AtomCoordinate>> backBonePerResidue = splitAtomByResidue(atom);
    std::string currResidue = getThreeLetterCode(this->header.firstResidue);
    std::vector<AtomCoordinate> fullResidue;

    success = this->_continuizeSideChainTorsionAngles(
        this->sideChainAnglesDiscretized, this->sideChainAnglesPerResidue
    );

    for (size_t i = 0; i < backBonePerResidue.size(); i++) {
        if (i != 0){
            currResidue = backBonePerResidue[i][0].residue;
        }
        fullResidue = nerf.reconstructAminoAcid(
            backBonePerResidue[i], this->sideChainAnglesPerResidue[i], AAS.at(currResidue)
        );
        if (this->useAltAtomOrder) {
            _reorderAtoms(fullResidue, AAS.at(currResidue));
        }
        backBonePerResidue[i] = fullResidue;
    }
    // Flatten backBonePerResidue & set tempFactor
    atom.clear();
    // Prepare tempFactor
    std::vector<float> tempFactors = this->tempFactorsDisc.continuize(this->tempFactorsDiscretized);
    this->tempFactors = tempFactors;
    // Iterate over each residue
    for (size_t i = 0; i < backBonePerResidue.size(); i++) {
        for (size_t j = 0; j < backBonePerResidue[i].size(); j++) {
            backBonePerResidue[i][j].tempFactor = tempFactors[i];
            atom.push_back(backBonePerResidue[i][j]);
        }
    }
    // Reindex atom index of atom
    if (this->hasOXT) {
        // Set OXT tempFactor
        this->OXT.tempFactor = tempFactors.back();
        atom.push_back(this->OXT);
    }
    setAtomIndexSequentially(atom, this->header.idxAtom);

    return success;
}

int Foldcomp::read(std::istream & file) {
    int success;
    // Open file in reading binary mode

    // Check the file starts with magic number
    char mNum[MAGICNUMBER_LENGTH];
    file.read(mNum, MAGICNUMBER_LENGTH);
    for (int i = 0; i < MAGICNUMBER_LENGTH; i++) {
        if (mNum[i] != MAGICNUMBER[i]) {
            return -1;
        }
    }
    // Read the header
    file.read(reinterpret_cast<char*>(&this->header), sizeof(this->header));
    this->read_header(this->header);
    // Read anchorIndices
    this->anchorIndices.resize(this->nAllAnchor);
    // Read int vector
    file.read(reinterpret_cast<char*>(&this->anchorIndices[0]), sizeof(int) * this->nAllAnchor);
    // Read the title
    this->strTitle = std::string(this->header.lenTitle, '\0');
    file.read(&this->strTitle[0], sizeof(char) * this->header.lenTitle);

    // Read the prev atoms
    // In the file, only xyz coordinates are stored
    // So, we need to reconstruct the atomcoordinate from the xyz coordinates & the information from the header
    float prevAtomCoords[9];
    file.read(reinterpret_cast<char*>(prevAtomCoords), sizeof(prevAtomCoords));

    // ANCHOR ATOMS
    if (this->header.nAnchor > 2) {
        float innerAnchorCoords[3];
        for (int i = 0; i < (this->header.nAnchor - 2); i++) {
            std::vector<std::vector<float>> innerAnchorCoord;
            for (int j = 0; j < 3; j++) {
                file.read(reinterpret_cast<char*>(innerAnchorCoords), sizeof(innerAnchorCoords));
                innerAnchorCoord.push_back({
                    innerAnchorCoords[0],
                    innerAnchorCoords[1],
                    innerAnchorCoords[2]
                });
            }
            this->anchorCoordinates.push_back(innerAnchorCoord);
        }
    }

    // LAST ATOMS
    float lastAtomCoords[9];
    file.read(reinterpret_cast<char*>(lastAtomCoords), sizeof(lastAtomCoords));
    for (int i = 0; i < 3; i++) {
        this->lastAtomCoordinates.push_back({ lastAtomCoords[i*3], lastAtomCoords[i*3 + 1], lastAtomCoords[i*3+ 2] });
    }
    this->anchorCoordinates.push_back(this->lastAtomCoordinates);

    file.read(&this->hasOXT, sizeof(char));
    float oxtCoords[3];
    file.read(reinterpret_cast<char*>(oxtCoords), sizeof(oxtCoords));
    this->OXT_coords = { oxtCoords[0], oxtCoords[1], oxtCoords[2] };
    this->OXT = AtomCoordinate(
        "OXT", getThreeLetterCode(this->header.lastResidue), std::string(1, this->header.chain),
        (int)this->header.nAtom, (int)this->header.nResidue, this->OXT_coords
    );

    // Read sidechain discretizers
    // file.read(reinterpret_cast<char*>(&this->sideChainDisc), sizeof(SideChainDiscretizers));time

    // Read the array of backbone bytes
    // Read 8 bytes at a time
    // Converted data will be saved at this->compressedBackBone
    // if this->compressedBackBone is not empty, then it will be overwritten
    // check if compressedBackBone is empty
    if (this->compressedBackBone.size() == 0) {
        this->compressedBackBone.resize(this->header.nResidue);
    } else {
        this->compressedBackBone.clear();
        this->compressedBackBone.resize(this->header.nResidue);
    }
    char* buffer = new char[8];
    for (int i = 0; i < this->header.nResidue; i++) {
        file.read(buffer, 8);
        // Convert buffer to compressedBackBone using convertBytesToBackboneChain
        this->compressedBackBone[i] = convertBytesToBackboneChain(buffer);
    }
    delete[] buffer;

    // Read sidechain
    int encodedSideChainSize = this->header.nSideChainTorsion;
    if (encodedSideChainSize % 2 == 1) {
        encodedSideChainSize++;
    }
    encodedSideChainSize /= 2;
    unsigned char* encodedSideChain = new unsigned char[this->header.nSideChainTorsion];
    // file.read(reinterpret_cast<char*>(encodedSideChain), sizeof(encodedSideChain));
    // Decode sidechain
    // success = decodeSideChainTorsionVector(
    //     encodedSideChain, this->header.nSideChainTorsion,
    //     this->sideChainAnglesDiscretized
    // );
    file.read(reinterpret_cast<char*>(encodedSideChain), this->header.nSideChainTorsion);
    // Read char array to unsigned int vector sideChainAnglesDiscretized
    unsigned int temp;
    for (size_t i = 0; i < this->header.nSideChainTorsion; i++) {
        // Convert char to unsigned int
        temp = encodedSideChain[i];
        this->sideChainAnglesDiscretized.push_back(temp);
    }
    delete[] encodedSideChain;

    // Read temperature factor
    // Read discretizer for temperature factor
    file.read(reinterpret_cast<char*>(&this->tempFactorsDisc.min), sizeof(float));
    file.read(reinterpret_cast<char*>(&this->tempFactorsDisc.cont_f), sizeof(float));

    unsigned char* encodedTempFactors = new unsigned char[this->header.nResidue];
    file.read(reinterpret_cast<char*>(encodedTempFactors), this->header.nResidue);

    unsigned int tempFactor;
    for (int i = 0; i < this->header.nResidue; i++) {
        tempFactor = (unsigned int)encodedTempFactors[i];
        this->tempFactorsDiscretized.push_back(tempFactor);
    }
    delete[] encodedTempFactors;

    success = _restoreAtomCoordinate(prevAtomCoords);
    if (success != 0) {
        return -2;
    }
    for (int i = 0; i < 6; i++) {
        success = _restoreDiscretizer(i);
    }
    // Close file
    return success;
}

int Foldcomp::writeStream(std::ostream& os) {
    int flag = 0;
    // Write magic number
    os.write(MAGICNUMBER, MAGICNUMBER_LENGTH);
    // Write header
    os.write((char*)&this->header, sizeof(CompressedFileHeader));
    // Write anchorIndices
    for (size_t i = 0; i < this->anchorIndices.size(); i++) {
        os.write((char*)&this->anchorIndices[i], sizeof(int));
    }
    // Write title
    os.write(this->strTitle.c_str(), this->strTitle.length());

    // 2022-08-08 19:15:30 - Changed to write all anchor atoms
    // TODO: NEED TO BE CHECKED
    for (const auto& anchors : this->anchorAtoms) {
        for (int i = 0; i < 3; i++) {
            os.write((char*)&anchors[i].coordinate.x, sizeof(float));
            os.write((char*)&anchors[i].coordinate.y, sizeof(float));
            os.write((char*)&anchors[i].coordinate.z, sizeof(float));
        }
    }

    os.write(&this->hasOXT, sizeof(char));
    os.write((char*)&this->OXT_coords.x, sizeof(float));
    os.write((char*)&this->OXT_coords.y, sizeof(float));
    os.write((char*)&this->OXT_coords.z, sizeof(float));

    // END OF BACKBONE METADATA
    // // Write sidechain discretizers
    // os.write((char*)&this->sideChainDisc, sizeof(SideChainDiscretizers));

    // Write the compressed backbone
    char* buffer = new char[8];
    for (size_t i = 0; i < this->compressedBackBone.size(); i++) {
        flag = convertBackboneChainToBytes(this->compressedBackBone[i], buffer);
        os.write(buffer, 8);
    }
    delete[] buffer;
    // 2022-06-07 18:44:09 - Check memory leak
    // Write side chain torsion angles
    int encodedSideChainSize = this->nSideChainTorsion;
    if (encodedSideChainSize % 2 == 1) {
        encodedSideChainSize++;
    }
    encodedSideChainSize /= 2;
    // char* charSideChainTorsion = encodeSideChainTorsionVector(this->sideChainAnglesDiscretized);
    // os.write(charSideChainTorsion, encodedSideChainSize);
    // TODO: TESTING
    // Convert unsigned int to char array
    unsigned char* charSideChainTorsion = new unsigned char[this->nSideChainTorsion];
    // Get array of unsigned int from sideChainAnglesDiscretized and convert to char array
    for (int i = 0; i < this->nSideChainTorsion; i++) {
        // convert unsigned int to char
        charSideChainTorsion[i] = (unsigned char)this->sideChainAnglesDiscretized[i];
    }
    os.write((char*)charSideChainTorsion, this->sideChainAnglesDiscretized.size());
    delete[] charSideChainTorsion;

    // Write temperature factors
    // Disc
    os.write((char*)&this->tempFactorsDisc.min, sizeof(float));
    os.write((char*)&this->tempFactorsDisc.cont_f, sizeof(float));

    unsigned char* charTempFactors = new unsigned char[this->header.nResidue];
    for (int i = 0; i < this->header.nResidue; i++) {
        charTempFactors[i] = (unsigned char)this->tempFactorsDiscretized[i];
    }
    os.write((char*)charTempFactors, this->tempFactorsDiscretized.size());
    delete[] charTempFactors;
    return flag;
}

int Foldcomp::write(std::string filename) {
    // Open in binary & writing mode
    std::ofstream outfile(filename, std::ios::out | std::ios::binary);
    if (!outfile) {
        return -1;
    }
    return writeStream(outfile);
}

#ifdef FOLDCOMP_EXECUTABLE
// 2022-08-29 15:42:29 TAR format support
int Foldcomp::writeTar(mtar_t& tar, std::string filename, size_t size) {
    int flag = 0;
    mtar_write_file_header(&tar, filename.c_str(), size);
    // Magic number
    mtar_write_data(&tar, MAGICNUMBER, MAGICNUMBER_LENGTH);
    // Write header
    mtar_write_data(&tar, &this->header, sizeof(CompressedFileHeader));
    // Write anchorIndices
    for (size_t i = 0; i < this->anchorIndices.size(); i++) {
        mtar_write_data(&tar, &this->anchorIndices[i], sizeof(int));
    }
    // Write title
    mtar_write_data(&tar, this->strTitle.c_str(), this->strTitle.length());
    // Write anchor atoms
    for (const auto& anchors : this->anchorAtoms) {
        for (int i = 0; i < 3; i++) {
            mtar_write_data(&tar, &anchors[i].coordinate.x, sizeof(float));
            mtar_write_data(&tar, &anchors[i].coordinate.y, sizeof(float));
            mtar_write_data(&tar, &anchors[i].coordinate.z, sizeof(float));
        }
    }
    // Write hasOXT
    mtar_write_data(&tar, &this->hasOXT, sizeof(char));
    // Write OXT_coords
    mtar_write_data(&tar, &this->OXT_coords.x, sizeof(float));
    mtar_write_data(&tar, &this->OXT_coords.y, sizeof(float));
    mtar_write_data(&tar, &this->OXT_coords.z, sizeof(float));
    // Write sideChainDisc
    // mtar_write_data(&tar, &this->sideChainDisc, sizeof(SideChainDiscretizers));
    // Write the compressed backbone
    char* buffer = new char[8];
    for (size_t i = 0; i < this->compressedBackBone.size(); i++) {
        flag = convertBackboneChainToBytes(this->compressedBackBone[i], buffer);
        mtar_write_data(&tar, buffer, 8);
    }
    delete[] buffer;
    // Write side chain torsion angles
    int encodedSideChainSize = this->nSideChainTorsion;
    if (encodedSideChainSize % 2 == 1) {
        encodedSideChainSize++;
    }
    encodedSideChainSize /= 2;

    unsigned char* charSideChainTorsion = new unsigned char[this->nSideChainTorsion];
    // Get array of unsigned int from sideChainAnglesDiscretized and convert to char array
    for (int i = 0; i < this->nSideChainTorsion; i++) {
        // convert unsigned int to char
        charSideChainTorsion[i] = (unsigned char)this->sideChainAnglesDiscretized[i];
    }
    mtar_write_data(&tar, charSideChainTorsion, this->sideChainAnglesDiscretized.size());
    delete[] charSideChainTorsion;
    // Write temperature factors
    // Disc
    mtar_write_data(&tar, &this->tempFactorsDisc.min, sizeof(float));
    mtar_write_data(&tar, &this->tempFactorsDisc.cont_f, sizeof(float));
    // Convert unsigned int to char array
    unsigned char* charTempFactors = new unsigned char[this->header.nResidue];
    // Get array of unsigned int from tempFactorsDiscretized and convert to char array
    for (int i = 0; i < this->header.nResidue; i++) {
        // convert unsigned int to char
        charTempFactors[i] = (unsigned char)this->tempFactorsDiscretized[i];
    }
    mtar_write_data(&tar, charTempFactors, this->tempFactorsDiscretized.size());
    delete[] charTempFactors;
    return flag;
}
#endif

size_t Foldcomp::getSize() {
    // Calculate the size of the compressed format
    size_t size = 0;
    // Magic number
    size += MAGICNUMBER_LENGTH;
    // Header
    size += sizeof(CompressedFileHeader);
    // Anchor indices
    size += sizeof(int) * this->anchorIndices.size();
    // Title
    size += this->strTitle.length();
    // Anchor atoms
    size += sizeof(float) * 3 * 3 * this->anchorAtoms.size();
    // OXT
    size += sizeof(char);
    size += sizeof(float) * 3;
    // Backbone
    size += sizeof(char) * 8 * this->compressedBackBone.size();
    // Side chain torsion angles
    size += sizeof(unsigned char) * this->sideChainAnglesDiscretized.size();
    // Temperature factors
    size += sizeof(float) * 2;
    size += sizeof(unsigned char) * this->tempFactorsDiscretized.size();
    return size;
}


// Functions to extract temperature factors only
int Foldcomp::continuizeTempFactors() {
    this->tempFactors = this->tempFactorsDisc.continuize(this->tempFactorsDiscretized);
    return 0;
}

int Foldcomp::writeFASTALike(std::ostream& os, const std::string& data) {
    // Output format
    // >title
    // 95461729... 889 // plddt of all residues converted to 1 decimal place or
    // MKLLSKPR... YVK // amino acid sequence
    // Write title
    os << ">" << this->strTitle << "\n" << data << "\n";
    return 0;
}

int Foldcomp::writeTSV(std::ostream& os, const std::string& data) {
    // Write title
    os << this->strTitle << "\t" << this->nResidue << "\t" << data << "\n";
    return 0;
}   

int Foldcomp::writeTorsionAngles(std::string filename) {
    int flag = 0;
    std::ofstream outfile(filename, std::ios::out);
    // Write header
    outfile << "index,phi,psi,omega" << std::endl;
    // Write backbone torsion angles
    for (size_t i = 0; i < this->phi.size(); i++) {
        outfile << i << "," << this->phi[i] << "," << this->psi[i] << "," << this->omega[i] << std::endl;
    }
    // Done
    outfile.close();
    return flag;
}

/**
 * @brief Extract information from the compressed file and write to a FASTA-like file
 *
 * @param filename
 * @param type 0: plddt, 1: sequence
 * @return int
 */
int Foldcomp::extract(std::string& data, int type, int digits) {
    int flag = 0;
    if (type == 0) {
        // Extract temperature factors
        this->continuizeTempFactors();
        if (digits < 1) {
            digits = 1;
        } else if (digits > 4) {
            digits = 4;
        }
        // Reserve string size
        int reserving;
        if (digits == 1) {
            reserving = this->tempFactors.size();
        } else if (digits == 2) {
            // 2 digits, comma separated
            reserving = this->tempFactors.size() * 3;
        } else if (digits == 3) {
            // 3 digits, one decimal place, comma separated
            reserving = this->tempFactors.size() * 5;
        } else {
            // 4 digits, two decimal places, comma separated
            reserving = this->tempFactors.size() * 6;
        }
        data.reserve(reserving);

        for (size_t i = 0; i < this->tempFactors.size(); i++) {
            float clamped = (std::clamp(this->tempFactors[i], 0.0f, 100.0f));
            char digit1 = (char)(clamped / 10.0f) + '0';
            char digit2 = (char)((int)clamped % 10) + '0';

            data.append(1, digit1);
            if (digits > 1) {
                data.append(1, digit2);
            }

            if (digits >= 3) {
                char digit3 = (char)((int)(clamped * 10.0f) % 10) + '0';
                data.append(1, '.');
                data.append(1, digit3);
            }

            if (digits == 4) {
                char digit4 = (char)((int)(clamped * 100.0f) % 10) + '0';
                data.append(1, digit4);
            }

            if (digits > 1 && i != this->tempFactors.size() - 1) {
                // Add comma separator
                data.append(1, ',');
            }
        }
    } else if (type == 1) {
        // Extract sequence
        data.reserve(this->header.nResidue);
        for (int i = 0; i < this->header.nResidue; i++) {
            char res = convertIntToOneLetterCode(this->compressedBackBone[i].residue);
            data.append(1, res);
        }
    }
    return flag;
}


//
CompressedFileHeader Foldcomp::get_header() {
    CompressedFileHeader header;
    // counts
    header.nResidue = this->nResidue;
    header.nAtom = this->nAtom;
    header.idxResidue = this->idxResidue;
    header.idxAtom = this->idxAtom;
    header.nAnchor = this->nAllAnchor;
    header.nSideChainTorsion = this->nSideChainTorsion;
    header.firstResidue = this->firstResidue;
    header.lastResidue = this->lastResidue;
    header.lenTitle = this->lenTitle;
    header.chain = this->chain;
    // discretizer parameters
    header.mins[0] = this->phiDisc.min;
    header.mins[1] = this->psiDisc.min;
    header.mins[2] = this->omegaDisc.min;
    header.mins[3] = this->n_ca_c_angleDisc.min;
    header.mins[4] = this->ca_c_n_angleDisc.min;
    header.mins[5] = this->c_n_ca_angleDisc.min;
    header.cont_fs[0] = this->phiDisc.cont_f;
    header.cont_fs[1] = this->psiDisc.cont_f;
    header.cont_fs[2] = this->omegaDisc.cont_f;
    header.cont_fs[3] = this->n_ca_c_angleDisc.cont_f;
    header.cont_fs[4] = this->ca_c_n_angleDisc.cont_f;
    header.cont_fs[5] = this->c_n_ca_angleDisc.cont_f;
    return header;
}

int Foldcomp::read_header(CompressedFileHeader& header) {
    this->nResidue = header.nResidue;
    this->nAtom = header.nAtom;
    this->idxResidue = header.idxResidue;
    this->idxAtom = header.idxAtom;
    this->nAllAnchor = header.nAnchor;
    this->nSideChainTorsion = header.nSideChainTorsion;
    this->firstResidue = header.firstResidue;
    this->lastResidue = header.lastResidue;
    this->lenTitle = header.lenTitle;
    //
    this->chain = header.chain;
    // discretizer parameters
    this->phiDisc.min = header.mins[0];
    this->psiDisc.min = header.mins[1];
    this->omegaDisc.min = header.mins[2];
    this->n_ca_c_angleDisc.min = header.mins[3];
    this->ca_c_n_angleDisc.min = header.mins[4];
    this->c_n_ca_angleDisc.min = header.mins[5];
    this->phiDisc.cont_f = header.cont_fs[0];
    this->psiDisc.cont_f = header.cont_fs[1];
    this->omegaDisc.cont_f = header.cont_fs[2];
    this->n_ca_c_angleDisc.cont_f = header.cont_fs[3];
    this->ca_c_n_angleDisc.cont_f = header.cont_fs[4];
    this->c_n_ca_angleDisc.cont_f = header.cont_fs[5];
    return 0;
}

void Foldcomp::print(int length) {
    // Print the header
    std::cout << "[Header]" << std::endl;
    std::cout << "nResidue: " << this->header.nResidue << std::endl;
    std::cout << "nAtom: " << this->header.nAtom << std::endl;
    std::cout << "idxResidue: " << this->header.idxResidue << std::endl;
    std::cout << "idxAtom: " << this->header.idxAtom << std::endl;
    std::cout << "nSideChainTorsion: " << this->header.nSideChainTorsion << std::endl;
    std::cout << "mins: " << std::endl;
    for (int i = 0; i < 6; i++) {
        std::cout << this->header.mins[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "cont_fs: " << std::endl;
    for (int i = 0; i < 6; i++) {
        std::cout << this->header.cont_fs[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "--------------------" << std::endl;

    // Print the prevAtoms
    std::cout << "[PrevAtoms]" << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << "Atom " << this->prevAtoms[i].atom << ": " << std::endl;
        std::cout << "x: " << this->prevAtoms[i].coordinate.x << std::endl;
        std::cout << "y: " << this->prevAtoms[i].coordinate.y << std::endl;
        std::cout << "z: " << this->prevAtoms[i].coordinate.z << std::endl;
    }
    std::cout << "--------------------" << std::endl;

    // Print the first element of compressedBackBone
    std::cout << "[CompressedBackbone]" << std::endl;
    for (int i = 0; i < length; i++) {
        std::cout << "Residue: " << this->compressedBackBone[i].residue << std::endl;
        std::cout << "phi-disc: " << this->compressedBackBone[i].phi << " / ";
        std::cout << this->phiDisc.continuize(this->compressedBackBone[i].phi) << std::endl;
        std::cout << "psi-disc: " << this->compressedBackBone[i].psi << " / ";
        std::cout << this->psiDisc.continuize(this->compressedBackBone[i].psi) << std::endl;
        std::cout << "omega-disc: " << this->compressedBackBone[i].omega << " / ";
        std::cout << this->omegaDisc.continuize(this->compressedBackBone[i].omega) << std::endl;
        std::cout << "n_ca_c_angle-disc: " << this->compressedBackBone[i].n_ca_c_angle << " / ";
        std::cout << this->n_ca_c_angleDisc.continuize(this->compressedBackBone[i].n_ca_c_angle) << std::endl;
        std::cout << "ca_c_n_angle-disc: " << this->compressedBackBone[i].ca_c_n_angle << " / ";
        std::cout << this->ca_c_n_angleDisc.continuize(this->compressedBackBone[i].ca_c_n_angle) << std::endl;
        std::cout << "c_n_ca_angle-disc: " << this->compressedBackBone[i].c_n_ca_angle << " / ";
        std::cout << this->c_n_ca_angleDisc.continuize(this->compressedBackBone[i].c_n_ca_angle) << std::endl;
    }
}

/*
void Foldcomp::printSideChainTorsion(std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    //
    outfile << "ResidueInd,Residue,Type,Key,RawVal,DiscVal,DiscMin,DiscContF,ReconVal,Diff\n";
    int movingIndex = 0;
    std::vector< std::vector<AtomCoordinate> > atomByResidue = splitAtomByResidue(this->rawAtoms);
    std::map<std::string, float> currBondAngle;
    std::map<std::string, float> currBondLength;
    std::map<std::string, float> currTorsionAngle;

    for (int i = 0; i < this->nResidue; i++) {
        currBondAngle = calculateBondAngles(atomByResidue[i], this->AAS.at(this->residueThreeLetter[i]));
        currBondLength = calculateBondLengths(atomByResidue[i], this->AAS.at(this->residueThreeLetter[i]));
        for (const auto& bl: currBondLength) {
            outfile << i << "," << this->residueThreeLetter[i] << ",BondLength,";
            outfile << bl.first << "," << bl.second << ",NA,NA,NA,";
            outfile << this->AAS.at(this->residueThreeLetter[i]).bondLengths.at(bl.first) << ",";
            outfile << bl.second - this->AAS.at(this->residueThreeLetter[i]).bondLengths.at(bl.first) << "\n";
        }
        for (const auto& ba : currBondAngle) {
            outfile << i << "," << this->residueThreeLetter[i] << ",BondAngle,";
            outfile << ba.first << "," << ba.second << ",NA,NA,NA," << this->AAS.at(this->residueThreeLetter[i]).bondAngles.at(ba.first) << ",";
            outfile << ba.second - this->AAS.at(this->residueThreeLetter[i]).bondAngles.at(ba.first) << "\n";
        }
        for (size_t j = 0; j < this->sideChainAnglesPerResidue[i].size(); j++) {
            outfile << i << "," << this->residueThreeLetter[i] << ",TorsionAngle,";
            outfile << j << "," << this->sideChainAnglesPerResidue[i][j] << ",";
            outfile << this->sideChainAnglesDiscretized[movingIndex] << ",";
            outfile << this->sideChainDiscMap[this->residueThreeLetter[i]][j].min << ",";
            outfile << this->sideChainDiscMap[this->residueThreeLetter[i]][j].cont_f << ",";
            outfile << this->sideChainDiscMap[this->residueThreeLetter[i]][j].continuize(this->sideChainAnglesDiscretized[movingIndex]) << ",";
            outfile << this->sideChainAnglesPerResidue[i][j] - this->sideChainDiscMap[this->residueThreeLetter[i]][j].continuize(this->sideChainAnglesDiscretized[movingIndex]) << "\n";
            movingIndex++;
        }
    }
    outfile.close();
}
*/

/**
 * @brief Checks the input file read is valid or not.
 *        This method is expected to be called after the input file is read.
 * @return int Error code. 0 if no error.
 */
ValidityError Foldcomp::checkValidity() {
    // Check size
    bool hasCorrectNumResidue = (this->header.nResidue == this->compressedBackBone.size());
    bool hasCorrectNumSideChain = (this->header.nSideChainTorsion == this->sideChainAnglesDiscretized.size());
    bool hasCorrectNumTempfactor = (this->header.nResidue == this->tempFactorsDiscretized.size());
    // Check vectors are not empty
    // For backbone, just check the torsion angles
    bool emptyBackbone = std::all_of(
        this->compressedBackBone.begin(),
        this->compressedBackBone.end(),
        [](BackboneChain bb) {
            return bb.phi == 0 && bb.psi == 0 && bb.omega == 0;
        }
    );
    bool emptySideChain = std::all_of(
        this->sideChainAnglesDiscretized.begin(),
        this->sideChainAnglesDiscretized.end(),
        [](unsigned int i){return i == 0;}
    );
    bool emptyTempFactor = std::all_of(
        this->tempFactorsDiscretized.begin(),
        this->tempFactorsDiscretized.end(),
        [](unsigned int i){return i == 0;}
    );
    // Return
    if (!hasCorrectNumResidue) {
        return E_BACKBONE_COUNT_MISMATCH;
    } else if (!hasCorrectNumSideChain) {
        return E_SIDECHAIN_COUNT_MISMATCH;
    } else if (!hasCorrectNumTempfactor) {
        return E_TEMP_FACTOR_COUNT_MISMATCH;
    } else if (emptyBackbone) {
        return E_EMPTY_BACKBONE_ANGLE;
    } else if (emptySideChain) {
        return E_EMPTY_SIDECHAIN_ANGLE;
    } else if (emptyTempFactor) {
        return E_EMPTY_TEMP_FACTOR;
    } else {
        return SUCCESS;
    }
}

void printValidityError(ValidityError err, std::string& filename) {
    // Print error message to stderr with filename
    switch (err) {
        case E_BACKBONE_COUNT_MISMATCH:
            std::clog << "[Error] Number of backbone angles does not match header: " << filename << std::endl;
            break;
        case E_SIDECHAIN_COUNT_MISMATCH:
            std::clog << "[Error] Number of sidechain angles does not match header: " << filename << std::endl;
            break;
        case E_TEMP_FACTOR_COUNT_MISMATCH:
            std::clog << "[Error] Number of temperature factors does not match header: " << filename << std::endl;
            break;
        case E_EMPTY_BACKBONE_ANGLE:
            std::clog << "[Error] All backbone angles are empty: " << filename << std::endl;
            break;
        case E_EMPTY_SIDECHAIN_ANGLE:
            std::clog << "[Error] All sidechain angles are empty: " << filename << std::endl;
            break;
        case E_EMPTY_TEMP_FACTOR:
            std::clog << "[Error] All temperature factors are empty: " << filename << std::endl;
            break;
        case SUCCESS:
            break;
        default:
            std::clog << "[Error] Unknown error: " << filename << std::endl;
            break;
    }
}

void _reorderAtoms(std::vector<AtomCoordinate>& atoms, const AminoAcid& aa) {
    std::vector<AtomCoordinate> newAtoms = atoms;
    for (size_t i = 0; i < atoms.size(); i++) {
        if (atoms[i].atom == aa.altAtoms[i]) {
            continue;
        } else {
            for (size_t j = 0; j < aa.altAtoms.size(); j++) {
                if (atoms[i].atom == aa.altAtoms[j]) {
                    newAtoms[j] = atoms[i];
                }
            }
        }
    }
    atoms = newAtoms;
}


// Sidechain

char* encodeSideChainTorsionVector(std::vector<unsigned int> vector) {
    //
    size_t size = vector.size();
    size_t newSize = size;
    if (size % 2 == 1) {
        newSize++;
    }
    newSize /= 2;
    char* output = new char[size];
    char temp;
    for (size_t i = 0; i < newSize; i++) {
        // First 4 bits
        temp = vector[i * 2] & 0x0F;
        temp = temp << 4;
        if (i * 2 + 1 < vector.size()) {
            temp |= (vector[i * 2 + 1] & 0x0F);
        } else {
            temp |= 0x0F;
        }
        output[i] = temp;
    }
    return output;
}

int decodeSideChainTorsionVector(char* input, int nTorsion, std::vector<unsigned int>& vector) {
    int size = nTorsion;
    if (size % 2 == 1) {
        size++;
    }
    size /= 2;
    unsigned int temp, first, second;
    if (!vector.empty()) {
        vector.clear();
    }
    for (int i = 0; i < size; i++) {
        temp = input[i];
        first = temp >> 4;
        second = temp & 0x0F;
        vector.push_back(first);
        if (i * 2 + 1 < nTorsion) {
            vector.push_back(second);
        }
    }
    return size;
}

unsigned char* encodeDiscretizedTempFactors(std::vector<unsigned int> vector) {
    unsigned char* output = new unsigned char[vector.size()];
    for (size_t i = 0; i < vector.size(); i++) {
        output[i] = (unsigned char)vector[i];
    }
    return output;
}

int decodeDiscretizedTempFactors(unsigned char* input, int size, std::vector<unsigned int>& vector) {
    int out = 0;
    unsigned int temp;
    // If vector is not empty, clear it
    if (!vector.empty()) {
        vector.clear();
    }
    for (int i = 0; i < size; i++) {
        temp = (unsigned int)input[i];
        vector.push_back(temp);
    }
    return out;
}

std::map<std::string, std::vector<Discretizer> > initializeSideChainDiscMap() {
    // Initialize the map
    std::map<std::string, std::vector<Discretizer> > discMap;
    // Get Amino acids list
    std::vector<std::string> aaNames = getAminoAcidList();
    int numTorsion = 0;
    // We can access to the specific angle discretizer with AA name and torsion index
    // ex) discmap["ALA"][0]
    for (size_t i = 0; i < aaNames.size(); i++) {
        numTorsion = getSideChainTorsionNum(aaNames[i]);
        discMap[aaNames[i]] = std::vector<Discretizer>(numTorsion);
    }
    return discMap;
}

float* getMinPointerFromSideChainDiscretizers(
    std::string& residue, SideChainDiscretizers& scDiscretizers
) {
    if (residue == "ALA") {
        return scDiscretizers.ala_min;
    } else if (residue == "ARG") {
        return scDiscretizers.arg_min;
    } else if (residue == "ASN") {
        return scDiscretizers.asn_min;
    } else if (residue == "ASP") {
        return scDiscretizers.asp_min;
    } else if (residue == "CYS") {
        return scDiscretizers.cys_min;
    } else if (residue == "GLN") {
        return scDiscretizers.gln_min;
    } else if (residue == "GLU") {
        return scDiscretizers.glu_min;
    } else if (residue == "GLY") {
        return scDiscretizers.gly_min;
    } else if (residue == "HIS") {
        return scDiscretizers.his_min;
    } else if (residue == "ILE") {
        return scDiscretizers.ile_min;
    } else if (residue == "LEU") {
        return scDiscretizers.leu_min;
    } else if (residue == "LYS") {
        return scDiscretizers.lys_min;
    } else if (residue == "MET") {
        return scDiscretizers.met_min;
    } else if (residue == "PHE") {
        return scDiscretizers.phe_min;
    } else if (residue == "PRO") {
        return scDiscretizers.pro_min;
    } else if (residue == "SER") {
        return scDiscretizers.ser_min;
    } else if (residue == "THR") {
        return scDiscretizers.thr_min;
    } else if (residue == "TRP") {
        return scDiscretizers.trp_min;
    } else if (residue == "TYR") {
        return scDiscretizers.tyr_min;
    } else if (residue == "VAL") {
        return scDiscretizers.val_min;
    } else {
        return NULL;
    }
}

float* getContFFromSideChainDiscretizers(
    std::string& residue, SideChainDiscretizers& scDiscretizers
) {
    if (residue == "ALA") {
        return scDiscretizers.ala_cont_fs;
    } else if (residue == "ARG") {
        return scDiscretizers.arg_cont_fs;
    } else if (residue == "ASN") {
        return scDiscretizers.asn_cont_fs;
    } else if (residue == "ASP") {
        return scDiscretizers.asp_cont_fs;
    } else if (residue == "CYS") {
        return scDiscretizers.cys_cont_fs;
    } else if (residue == "GLN") {
        return scDiscretizers.gln_cont_fs;
    } else if (residue == "GLU") {
        return scDiscretizers.glu_cont_fs;
    } else if (residue == "GLY") {
        return scDiscretizers.gly_cont_fs;
    } else if (residue == "HIS") {
        return scDiscretizers.his_cont_fs;
    } else if (residue == "ILE") {
        return scDiscretizers.ile_cont_fs;
    } else if (residue == "LEU") {
        return scDiscretizers.leu_cont_fs;
    } else if (residue == "LYS") {
        return scDiscretizers.lys_cont_fs;
    } else if (residue == "MET") {
        return scDiscretizers.met_cont_fs;
    } else if (residue == "PHE") {
        return scDiscretizers.phe_cont_fs;
    } else if (residue == "PRO") {
        return scDiscretizers.pro_cont_fs;
    } else if (residue == "SER") {
        return scDiscretizers.ser_cont_fs;
    } else if (residue == "THR") {
        return scDiscretizers.thr_cont_fs;
    } else if (residue == "TRP") {
        return scDiscretizers.trp_cont_fs;
    } else if (residue == "TYR") {
        return scDiscretizers.tyr_cont_fs;
    } else if (residue == "VAL") {
        return scDiscretizers.val_cont_fs;
    } else {
        return NULL;
    }
}

int getSideChainTorsionNum(std::string residue) {
    int out = 0;
    if (residue == "ALA") {
        out = 2;
    } else if (residue == "ARG") {
        out = 8;
    } else if (residue == "ASN") {
        out = 5;
    } else if (residue == "ASP") {
        out = 5;
    } else if (residue == "CYS") {
        out = 3;
    } else if (residue == "GLN") {
        out = 6;
    } else if (residue == "GLU") {
        out = 6;
    } else if (residue == "GLY") {
        out = 1;
    } else if (residue == "HIS") {
        out = 7;
    } else if (residue == "ILE") {
        out = 5;
    } else if (residue == "LEU") {
        out = 5;
    } else if (residue == "LYS") {
        out = 6;
    } else if (residue == "MET") {
        out = 5;
    } else if (residue == "PHE") {
        out = 8;
    } else if (residue == "PRO") {
        out = 4;
    } else if (residue == "SER") {
        out = 3;
    } else if (residue == "THR") {
        out = 4;
    } else if (residue == "TRP") {
        out = 11;
    } else if (residue == "TYR") {
        out = 9;
    } else if (residue == "VAL") {
        out = 4;
    } else {
        out = 0;
    }
    return out;
}

const std::map<std::string, AminoAcid> Foldcomp::AAS = AminoAcid::AminoAcids();
