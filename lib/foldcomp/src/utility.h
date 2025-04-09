/**
 * File: utility.h
 * Project: foldcomp
 * Created: 2021-01-05 14:27:35
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Utility functions
 * ---
 * Last Modified: Fri Mar 03 2023
 * Modified By: Hyunbin Kim
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#pragma once
#include <cstring>
#include <fstream> // IWYU pragma: keep
#include <iostream>
#include <map>
#include <vector>

#ifdef _WIN32
#define NOMINMAX
#ifdef _MSC_VER
#include <dirent.h>
#else
#include "windows/dirent.h"
#endif // _MSC_VER
#include <io.h>
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#else
#include <dirent.h>
#include <sys/mman.h>
#endif // _WIN32

template<typename T>
std::vector<T> vectorSlice(std::vector<T> const& v, int m, int n) {
    auto first = v.cbegin() + m;
    auto last = v.cbegin() + n + 1;

    std::vector<T> vec(first, last);
    return vec;
}

template<typename T>
/**
 * @brief Add the values of map2 to map1
 *
 * @param map1
 * @param map2
 * @return int
 */
int addMap(std::map<std::string, T>* m1, std::map<std::string, T>* m2) {
    size_t m1Size = m1->size();
    if (m1Size == 0) {
        // if m1 is empty, fill in with m2
        for (auto const& x : *m2) {
            (*m1)[x.first] = x.second;
        }
    } else {
        for (auto const& x : *m2) {
            // if the key from map2 is in map1
            if (m1->find(x.first) != m1->end()) {
                // add the values to map1
                (*m1)[x.first] += x.second;
            }
        }
    }
    return 0;
}

//
template <typename T>
int divideMapWithConst(std::map<std::string, T>* m, T c) {
    for (auto const& x: *m) {
        (*m)[x.first] = (x.second / c);
    }
    return 0;
}

template <typename T>
int printMapToFile(std::map<std::string, T>* m, std::string fileName) {
    std::ofstream outFile;
    outFile.open(fileName, std::ios::out);

    for (auto const& x : *m) {
        outFile << x.first << "," << (*m)[x.first] << "\n";
    }
    outFile.close();
    return 0;
}

template <typename T>
void printVector(std::vector<T> v) {
    // Open brackets
    std::cout << "[";
    // Print each element
    for (auto const& x : v) {
        std::cout << x << " ";
    }
    // Close brackets
    std::cout << "]" << std::endl;
}

// #ifdef FOLDCOMP_EXECUTABLE
// Get all files in a directory using dirent.h
std::vector<std::string> getFilesInDirectory(const std::string& dir, bool recursive);
char *file_map(FILE *file, ssize_t *size, int extra_flags = 0);
int file_unmap(char* mem, ssize_t size);
// #endif
std::string baseName(const std::string& path);

std::string getFileWithoutExt(const std::string& file);
std::pair<std::string, std::string> getFileParts(const std::string& file);
bool isCompressible(std::pair<std::string, std::string>& fileParts);
bool stringStartsWith(const std::string& prefix, const std::string& str, const size_t offset = 0);
bool stringEndsWith(const std::string& suffix, const std::string& str);
std::vector<std::string> stringSplit(const std::string& str, const std::string& sep);

// One letter code to Three letter code
char getOneLetterCode(std::string threeLetterCode);
std::string getThreeLetterCode(char oneLetterCode);

// Int to One letter code
char convertIntToOneLetterCode(unsigned int aab);
unsigned int convertOneLetterCodeToInt(char oneLetterCode);

// Int to Three letter code
std::string convertIntToThreeLetterCode(unsigned int aab);
unsigned int convertThreeLetterCodeToInt(std::string threeLetterCode);

// CONSTANTS FOR AMINO ACID NAMING
#define AA_ALA_INT 0
#define AA_ALA_STR "ALA"
#define AA_ALA_CHAR 'A'
#define AA_ARG_INT 1
#define AA_ARG_STR "ARG"
#define AA_ARG_CHAR 'R'
#define AA_ASN_INT 2
#define AA_ASN_STR "ASN"
#define AA_ASN_CHAR 'N'
#define AA_ASP_INT 3
#define AA_ASP_STR "ASP"
#define AA_ASP_CHAR 'D'
#define AA_CYS_INT 4
#define AA_CYS_STR "CYS"
#define AA_CYS_CHAR 'C'
#define AA_GLN_INT 5
#define AA_GLN_STR "GLN"
#define AA_GLN_CHAR 'Q'
#define AA_GLU_INT 6
#define AA_GLU_STR "GLU"
#define AA_GLU_CHAR 'E'
#define AA_GLY_INT 7
#define AA_GLY_STR "GLY"
#define AA_GLY_CHAR 'G'
#define AA_HIS_INT 8
#define AA_HIS_STR "HIS"
#define AA_HIS_CHAR 'H'
#define AA_ILE_INT 9
#define AA_ILE_STR "ILE"
#define AA_ILE_CHAR 'I'
#define AA_LEU_INT 10
#define AA_LEU_STR "LEU"
#define AA_LEU_CHAR 'L'
#define AA_LYS_INT 11
#define AA_LYS_STR "LYS"
#define AA_LYS_CHAR 'K'
#define AA_MET_INT 12
#define AA_MET_STR "MET"
#define AA_MET_CHAR 'M'
#define AA_PHE_INT 13
#define AA_PHE_STR "PHE"
#define AA_PHE_CHAR 'F'
#define AA_PRO_INT 14
#define AA_PRO_STR "PRO"
#define AA_PRO_CHAR 'P'
#define AA_SER_INT 15
#define AA_SER_STR "SER"
#define AA_SER_CHAR 'S'
#define AA_THR_INT 16
#define AA_THR_STR "THR"
#define AA_THR_CHAR 'T'
#define AA_TRP_INT 17
#define AA_TRP_STR "TRP"
#define AA_TRP_CHAR 'W'
#define AA_TYR_INT 18
#define AA_TYR_STR "TYR"
#define AA_TYR_CHAR 'Y'
#define AA_VAL_INT 19
#define AA_VAL_STR "VAL"
#define AA_VAL_CHAR 'V'
#define AA_ASX_INT 20
#define AA_ASX_STR "ASX"
#define AA_ASX_CHAR 'B'
#define AA_GLX_INT 21
#define AA_GLX_STR "GLX"
#define AA_GLX_CHAR 'Z'
#define AA_STOP_INT 22
#define AA_STOP_STR "STP"
#define AA_STOP_CHAR '*'
#define AA_UNK_INT 23
#define AA_UNK_STR "UNK"
#define AA_UNK_CHAR 'X'
#define NUM_ALL_AA_CODES 24
#define NUM_VALID_AA_CODES 20