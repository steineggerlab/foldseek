#ifndef FOLDSEEK_PDBWRITER_H
#define FOLDSEEK_PDBWRITER_H

#include <cstdio>

// three letter code of a one letter amino acid, anything outside A-Z becomes UNK
inline const char* threeLetterResidue(char aa) {
    static constexpr const char* lookup[26] = { "ALA", "ASX", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "XLE", "LYS", "LEU", "MET", "ASN", "PYL", "PRO", "GLN", "ARG", "SER", "THR", "SEC", "VAL", "TRP", "UNK", "TYR", "GLX" };
    // make AA upper case
    aa = aa & ~0x20;
    if (aa < 'A' || aa > 'Z') {
        return "UNK";
    }
    return lookup[aa - 'A'];
}

// hybrid-36 (http://cci.lbl.gov/hybrid_36/) format for more than 99999 atoms or 9999 residues
inline void writeHybrid36(char* out, int width, int value) {
    int decimalLimit = 1;
    int base36Offset = 10;
    for (int i = 0; i < width; ++i) {
        decimalLimit *= 10;
    }
    for (int i = 1; i < width; ++i) {
        base36Offset *= 36;
    }
    int base = 10;
    if (value >= decimalLimit) {
        base = 36;
        value = value - decimalLimit + base36Offset;
    }
    // decimal numbers are right aligned with spaces, hybrid-36 numbers always fill the width
    for (int i = width - 1; i >= 0; --i) {
        out[i] = (i == width - 1 || value > 0) ? "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[value % base] : ' ';
        value /= base;
    }
    out[width] = '\0';
}

// %8.3f overflows below -999.999 and above 9999.999, drop decimals to keep the eight columns
// returns true if decimals had to be dropped
inline bool writeCoordinate(char* out, float value) {
    for (int precision = 3; precision > 0; --precision) {
        if (snprintf(out, 9, "%8.*f", precision, value) <= 8) {
            return precision < 3;
        }
    }
    snprintf(out, 9, "%8.0f", value);
    return true;
}

// writes a C-alpha ATOM record into out (at least 82 bytes), chain takes one or two characters.
// overflow is set if a number did not fit its PDB columns.
inline int writeCalphaAtom(char* out, int serial, const char* residue, const char* chain, int resSeq,
                           float x, float y, float z, bool& overflow) {
    char serialStr[6], resSeqStr[5], xStr[9], yStr[9], zStr[9];
    overflow |= serial > 99999 || resSeq > 9999;
    overflow |= writeCoordinate(xStr, x);
    overflow |= writeCoordinate(yStr, y);
    overflow |= writeCoordinate(zStr, z);
    writeHybrid36(serialStr, 5, serial);
    writeHybrid36(resSeqStr, 4, resSeq);
    return snprintf(out, 82, "ATOM  %s  CA  %s%2.2s%s    %s%s%s  1.00  0.00           C\n",
                    serialStr, residue, chain, resSeqStr, xStr, yStr, zStr);
}

// writes the TER record closing a chain into out (at least 82 bytes)
inline int writeTer(char* out, int serial, const char* residue, const char* chain, int resSeq) {
    char serialStr[6], resSeqStr[5];
    writeHybrid36(serialStr, 5, serial);
    writeHybrid36(resSeqStr, 4, resSeq);
    return snprintf(out, 82, "TER   %s      %s%2.2s%s\n", serialStr, residue, chain, resSeqStr);
}

#endif
