/**
 * File: foldcomp.h
 * Project: foldcomp
 * Created: 2021-02-02 14:04:40
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Contributor: Milot Mirdita (milot@mirdita.de)
 * Description:
 *     This file contains main data structures for torsion angle compression and
 *     functions for handling them.
 * ---
 * Last Modified: 2024-08-08 19:37:12
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#pragma once
#include "float3d.h"

#include "atom_coordinate.h"
#include "discretizer.h"
#include "nerf.h"
#include "tcbspan.h"

#ifdef FOLDCOMP_EXECUTABLE
// TAR format handling - only for executable
#include "microtar/microtar.h"
#endif

#include <cstdint>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// forward declaration
class AminoAcid;

// CONSTANTS
#define NUM_TYPE_OF_ANGLES 6
#define MAGICNUMBER_LENGTH 4
#define MAGICNUMBER "FCMP"
// NUMBER OF BITS FOR ENCODING
#define NUM_BITS_PHI_PSI 12
#define NUM_BITS_OMEGA 11
#define NUM_BITS_BOND 8
#define NUM_BITS_RESIDUE 5
#define NUM_BITS_TEMP 8
#define NUM_BITS_SIDECHAIN 4

#define N_TO_CA_DIST 1.4581
#define CA_TO_C_DIST 1.5281
#define C_TO_N_DIST 1.3311
#define PRO_N_TO_CA_DIST 1.353

#define DEFAULT_ANCHOR_THRESHOLD 25

// ERROR CODES FOR CHECKING VALIDITY
enum ValidityError {
    SUCCESS = 0,
    E_BACKBONE_COUNT_MISMATCH,
    E_SIDECHAIN_COUNT_MISMATCH,
    E_TEMP_FACTOR_COUNT_MISMATCH,
    E_EMPTY_BACKBONE_ANGLE,
    E_EMPTY_SIDECHAIN_ANGLE,
    E_EMPTY_TEMP_FACTOR
};

// NOTE: ORDER OF BOND ANGLE: CA_C_N, C_N_CA, N_CA_C
// NOTE: THE ORDER OF TORSION ANGLE IS PSI->OMEGA->PHI
struct BackboneChain {
    // TOTAL BITS: 64 bits = 8 bytes
    uint64_t residue: NUM_BITS_RESIDUE;     // 5 bits
    uint64_t omega : NUM_BITS_OMEGA;        // 11 bits
    uint64_t psi : NUM_BITS_PHI_PSI;        // 12 bits
    uint64_t phi : NUM_BITS_PHI_PSI;        // 12 bits
    uint64_t ca_c_n_angle : NUM_BITS_BOND;  // 8 bits
    uint64_t c_n_ca_angle : NUM_BITS_BOND;  // 8 bits
    uint64_t n_ca_c_angle: NUM_BITS_BOND;   // 8 bits
};
static_assert(sizeof(BackboneChain) == 8, "BackboneChain must remain compatible");

// IDEA: Split backbone chain header & general header??
// First residue string should be saved in the header
// TODO: SAVE FIRST RESIDUE
struct BackboneChainHeader {
    uint16_t nResidue;   // 16 bits
    uint16_t nAtom;      // 16 bits
    uint16_t idxResidue; // 16 bits
    uint16_t idxAtom;    // 16 bits
    // char firstResidue; // 1 byte = 8 bits --> WILL BE APPLIED AFTER CHANGING BACKBONE CHAIN
    float prevAtoms[9];  // 9 * 4 byte
};
static_assert(sizeof(BackboneChainHeader) == 44, "BackboneChainHeader must remain compatible");

struct DecompressedBackboneChain {
    /* data */
    char residue;
    float n_ca_c_angle;
    float ca_c_n_angle;
    float c_n_ca_angle;
    float phi;
    float psi;
    float omega;
};

// IDEA: Implement offset array
// An array of offsets?
// struct FileOffset {
//     unsigned int firstBackboneStart: 32;
//     unsigned int sideChainStart: 32;
// };
// struct CompressedFileOffset {
//     unsigned int firstBackboneStart : 32;
//     unsigned int sideChainStart : 32;
// };

struct CompressedFileHeader {
    // Backbone
    uint16_t nResidue;   // 16 bits
    uint16_t nAtom;      // 16 bits
    uint16_t idxResidue; // 16 bits
    uint16_t idxAtom;    // 16 bits
    // Anchor points
    uint8_t nAnchor;     // 8 bits
    char chain;          // 8 bits + 8 bits padding TODO
    // Sidechain
    uint32_t nSideChainTorsion; // 32 bits
    char firstResidue;   // 8 bits
    char lastResidue;    // 8 bits + 8 bits padding TODO
    uint32_t lenTitle;   // 32 bits
    // Discretizer for backbone chain
    float mins[6];       // 6 * 32 bits
    float cont_fs[6];    // 6 * 32 bits
};
static_assert(sizeof(CompressedFileHeader) == 72, "CompressedFileHeader must remain compatible");

struct SideChainDiscretizers {
    float ala_min[2];
    float ala_cont_fs[2];
    float arg_min[8];
    float arg_cont_fs[8];
    float asn_min[5];
    float asn_cont_fs[5];
    float asp_min[5];
    float asp_cont_fs[5];
    float cys_min[3];
    float cys_cont_fs[3];
    float gln_min[6];
    float gln_cont_fs[6];
    float glu_min[6];
    float glu_cont_fs[6];
    float gly_min[1];
    float gly_cont_fs[1];
    float his_min[7];
    float his_cont_fs[7];
    float ile_min[5];
    float ile_cont_fs[5];
    float leu_min[5];
    float leu_cont_fs[5];
    float lys_min[6];
    float lys_cont_fs[6];
    float met_min[5];
    float met_cont_fs[5];
    float phe_min[6];
    float phe_cont_fs[6];
    float pro_min[4];
    float pro_cont_fs[4];
    float ser_min[3];
    float ser_cont_fs[3];
    float thr_min[4];
    float thr_cont_fs[4];
    float trp_min[11];
    float trp_cont_fs[11];
    float tyr_min[9];
    float tyr_cont_fs[9];
    float val_min[4];
    float val_cont_fs[4];
};

// Conversion
uint32_t convertCompressedResidueToFirst4Bytes(BackboneChain& res);
uint32_t convertCompressedResidueToSecond4Bytes(BackboneChain& res);
int convertBackboneChainToBytes(BackboneChain& res, char* output);
BackboneChain convertBytesToBackboneChain(char* bytes);
BackboneChain newBackboneChain(
    char residue, unsigned int phi, unsigned int psi, unsigned int omega,
    unsigned int n_ca_c_angle, unsigned int ca_c_n_angle, unsigned int c_n_ca_angle
);
BackboneChain newBackboneChain(
    unsigned int bResidue, unsigned int phi, unsigned int psi, unsigned int omega,
    unsigned int n_ca_c_angle, unsigned int ca_c_n_angle, unsigned int c_n_ca_angle
);

// TODO: Change header
DecompressedBackboneChain decompressBackboneChain(
    const BackboneChain& bb, const CompressedFileHeader& header
);
std::vector<DecompressedBackboneChain> decompressBackboneChain(
    const std::vector<BackboneChain>& bbv, const CompressedFileHeader& header
);

float _continuize(unsigned int input, float min, float cont_f);

// Reconstruct
std::vector<AtomCoordinate> reconstructBackboneAtoms(
    const std::vector<AtomCoordinate>& prevAtoms,
    const std::vector<BackboneChain>& backbone,
    const CompressedFileHeader& header
);

std::vector<AtomCoordinate> reconstructSidechainAtoms(
    std::vector<AtomCoordinate>& backBoneAtoms,
    std::vector<AtomCoordinate>& wholeAtoms,
    SideChainDiscretizers& sideChainDisc
);

int discretizeSideChainTorsionAngles(
    std::vector< std::vector<float> >& torsionPerResidue,
    std::vector<std::string>& residueNames,
    std::map<std::string, AminoAcid>& AAS,
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap,
    std::vector<unsigned int>& output
);


int continuizeSideChainTorsionAngles(
    std::vector<unsigned int>& torsionDiscretized,
    std::vector<std::string>& residueNames,
    std::map<std::string, AminoAcid>& AAS,
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap,
    std::vector< std::vector<float> >& output
);

float* getContFFromSideChainDiscretizers(
    std::string& residue, SideChainDiscretizers& scDiscretizers
);
float* getMinPointerFromSideChainDiscretizers(
    std::string& residue, SideChainDiscretizers& scDiscretizers
);
std::map<std::string, std::vector<Discretizer> > initializeSideChainDiscMap();
unsigned char* encodeDiscretizedTempFactors(std::vector<unsigned int> vector);
int decodeDiscretizedTempFactors(unsigned char* input, int size, std::vector<unsigned int>& vector);
char* encodeSideChainTorsionVector(std::vector<unsigned int> vector);
int decodeSideChainTorsionVector(char* input, int nTorsion, std::vector<unsigned int>& vector);
int getSideChainTorsionNum(std::string residue);
int fillSideChainDiscretizerMap(
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap
);

void _reorderAtoms(std::vector<AtomCoordinate>& atoms, const AminoAcid& aa);

// Print
void printCompressedResidue(BackboneChain& res);
void printValidityError(ValidityError err, std::string& filename);

struct FloatArrayWithDisc {
    unsigned short size;
    float min;
    float cont_f;
    float* array;
};

class Foldcomp {
private:
    /* private methods */
    int _restoreDiscretizer(int angleType);
    int _restoreAtomCoordinate(float* coords);
    int _preprocessBackbone();
    int _preprocessSideChain();

    // Anchor
    int _getAnchorNum(int threshold);
    void _setAnchor(const tcb::span<AtomCoordinate>& atomCoordinates);
    std::vector<AtomCoordinate> _getAnchorAtoms(bool includeStartAndEnd = true);

    int _discretizeSideChainTorsionAngles(
        std::vector< std::vector<float> >& input, std::vector<unsigned int>& output
    );

    int _continuizeSideChainTorsionAngles(
        std::vector<unsigned int>& input, std::vector< std::vector<float> >& output
    );

public:
    bool isPreprocessed = false;
    bool isCompressed = false;
    bool backwardReconstruction = true;
    bool useAltAtomOrder = false;
    // Number of atoms & residues
    int nResidue = 0;
    int nAtom = 0;
    int nBackbone = 0;
    int nSideChainTorsion = 0;
    int nInnerAnchor = 0;
    int nAllAnchor = 0;
    int anchorThreshold = DEFAULT_ANCHOR_THRESHOLD;
    // Indices for residue & atom
    int idxResidue = 0;
    int idxAtom = 0;
    char chain;
    char firstResidue;
    char lastResidue;
    char hasOXT = 1;

    // Metadata
    std::string strTitle;
    // const char* title;
    int lenTitle;
    std::map<std::string, std::string> strMetadata;
    std::map<std::string, std::vector<float> > floatMetadata;

    // Header
    CompressedFileHeader header;

    // Vectors for atoms
    std::vector<AtomCoordinate> prevAtoms; // 3 atoms
    std::vector<AtomCoordinate> lastAtoms; // 3 atoms
    std::vector< std::vector<float> > lastAtomCoordinates; // 3 atoms
    std::vector<AtomCoordinate> backbone;
    // 2022-08-05 20:47:57 - Anchors for reducing RMSD
    std::vector< std::vector<AtomCoordinate> > anchorAtoms;
    std::vector< std::vector< std::vector<float> > > anchorCoordinates;
    std::vector<int> anchorIndices;

    std::vector<BackboneChain> compressedBackBone;
    std::vector<unsigned int> compressedSideChain;
    std::vector<char> residues;
    std::vector<std::string> residueThreeLetter;
    AtomCoordinate OXT;
    float3d OXT_coords;

    // Angles
    std::vector<float> backboneTorsionAngles;
    std::vector<float> backboneBondAngles;
    std::vector<float> psi;
    std::vector<float> omega;
    std::vector<float> phi;
    std::vector<float> n_ca_c_angle;
    std::vector<float> ca_c_n_angle;
    std::vector<float> c_n_ca_angle;
    Discretizer psiDisc;
    std::vector<unsigned int> psiDiscretized;
    Discretizer omegaDisc;
    std::vector<unsigned int> omegaDiscretized;
    Discretizer phiDisc;
    std::vector<unsigned int> phiDiscretized;
    Discretizer n_ca_c_angleDisc;
    std::vector<unsigned int> n_ca_c_angleDiscretized;
    Discretizer ca_c_n_angleDisc;
    std::vector<unsigned int> ca_c_n_angleDiscretized;
    Discretizer c_n_ca_angleDisc;
    std::vector<unsigned int> c_n_ca_angleDiscretized;
    const static std::map<std::string, AminoAcid> AAS;
    Nerf nerf;
    // Sidechain angles
    std::vector<float> sideChainAngles;
    std::vector< std::vector<float> > sideChainAnglesPerResidue;
    std::vector<unsigned int> sideChainAnglesDiscretized;
    SideChainDiscretizers sideChainDisc;
    std::map<std::string, std::vector<Discretizer> > sideChainDiscMap;
    int encodedSideChainSize;
    // Temperature factors / pLDDT
    std::vector<float> tempFactors;
    std::vector<unsigned int> tempFactorsDiscretized;
    Discretizer tempFactorsDisc;

     // methods
    int preprocess(const tcb::span<AtomCoordinate>& atoms);
    std::vector<BackboneChain> compress(const tcb::span<AtomCoordinate>& atoms);
    int decompress(std::vector<AtomCoordinate>& atoms);
    int read(std::istream & filename);
    int writeStream(std::ostream& os);
    int write(std::string filename);
    // Read & write for tar files
#ifdef FOLDCOMP_EXECUTABLE
// FOLDCOMP PYTHON API DOESN'T HANDLE TAR FORMAT
    // int readTar(mtar_t& tar, std::string filename, size_t size);
    int writeTar(mtar_t& tar, std::string filename, size_t size);
#endif

    int reconstruct(std::vector<AtomCoordinate>& atoms, int mode);
    CompressedFileHeader get_header();
    int read_header(CompressedFileHeader& header);
    size_t getSize();
    // methods for getting plddt (tempFactors) or amino acid sequence
    int continuizeTempFactors();
    int writeFASTALike(std::ostream& os, const std::string& data);
    int writeTSV(std::ostream& os, const std::string& data);
    int writeTorsionAngles(std::string filename);
    int extract(std::string& data, int type, int digits);

    // temporary method for testing
    std::vector<float> checkTorsionReconstruction();
    void print(int length = 5);
    // void printSideChainTorsion(std::string filename);
    // Check validity
    ValidityError checkValidity();
};
