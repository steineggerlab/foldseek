/**
 * File: utility.cpp
 * Project: foldcomp
 * Created: 2021-01-05 14:29:04
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Utility functions
 * ---
 * Last Modified: Fri Mar 03 2023
 * Modified By: Hyunbin Kim
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
// TODO: TRIM WHITE SPACES FROM PDB LINE
#include "utility.h"

#include <cstdlib>
#include <fstream>
#include <string>

#include <errno.h>
#include <sys/stat.h>

// Get all files in a directory using dirent.h
void getdir(const std::string& dir, bool recursive, std::vector<std::string>& files) {
    struct stat st;
    if (stat(dir.c_str(), &st) == 0) {
        if (S_ISREG(st.st_mode) || S_ISLNK(st.st_mode)) {
            files.emplace_back(dir);
            return;
        }
    }
    DIR* dp;
    struct dirent* dirp;
    if ((dp = opendir(dir.c_str())) == NULL) {
        perror("error");
        exit(errno);
    }
    std::vector<std::string> dirs;
    while ((dirp = readdir(dp)) != NULL) {
        unsigned int type = dirp->d_type;
        if (type == DT_UNKNOWN) {
            // stat the file to determine its real type
            struct stat st;
            std::string path = dir + "/" + dirp->d_name;
            if (stat(path.c_str(), &st) == 0) {
                if (S_ISDIR(st.st_mode)) {
                    type = DT_DIR;
                } else if (S_ISREG(st.st_mode)) {
                    type = DT_REG;
                } else if (S_ISLNK(st.st_mode)) {
                    type = DT_LNK;
                } else {
                    continue;
                }
            }
        }
        if (type == DT_DIR) {
            if (recursive && std::string(dirp->d_name) != "." && std::string(dirp->d_name) != "..") {
                dirs.emplace_back(dir + "/" + dirp->d_name);
            } else {
                continue;
            }
        }
        if (type == DT_REG || type == DT_LNK) {
            files.emplace_back(dir + "/" + dirp->d_name);
        }
    }
    closedir(dp);

    for (std::vector<std::string>::const_iterator it = dirs.begin(); it != dirs.end(); ++it) {
        getdir(*it, recursive, files);
    }
}
std::vector<std::string> getFilesInDirectory(const std::string& dir, bool recursive) {
    std::vector<std::string> files;
    getdir(dir, recursive, files);
    return files;
}

char *file_map(FILE *file, ssize_t *size, int extra_flags) {
    struct stat sb;
    fstat(fileno(file), &sb);
    *size = sb.st_size;

    int fd = fileno(file);
#ifdef _MSC_VER
    HANDLE handle = (HANDLE)_get_osfhandle(fd);
    void* mapping = CreateFileMapping(handle, NULL, PAGE_READONLY, 0, 0, NULL);
    DWORD offsetLow  = DWORD(0 & 0xFFFFFFFF);
    DWORD offsetHigh = DWORD(0 >> 32);
    return (char *)MapViewOfFile(mapping, FILE_MAP_READ, offsetHigh, offsetLow, sb.st_size);
#else
    return (char *)mmap(NULL, (size_t)(*size), PROT_READ, MAP_PRIVATE | extra_flags, fd, 0);
#endif
}

int file_unmap(char *mem, ssize_t size) {
#ifdef _MSC_VER
    return UnmapViewOfFile(mem);
#else
    return munmap(mem, size);
#endif
}


std::string baseName(const std::string& path) {
    return path.substr(path.find_last_of("/\\") + 1);
}

std::string getFileWithoutExt(const std::string& file) {
    size_t basePos = file.find_last_of("/\\");
    basePos = basePos == std::string::npos ? 0 : basePos + 1;
    size_t extStart = file.substr(basePos).find_last_of(".");
    return extStart == std::string::npos ? file : file.substr(0, basePos + extStart);
}

std::pair<std::string, std::string> getFileParts(const std::string& file) {
    size_t basePos = file.find_last_of("/\\");
    basePos = basePos == std::string::npos ? 0 : basePos + 1;
    size_t extStart = file.substr(basePos).find_last_of(".");
    if (extStart == std::string::npos) {
        return std::make_pair(file, "");
    }
    return std::make_pair(file.substr(0, basePos + extStart), file.substr(basePos + extStart + 1));
}

// Allowable extensions: pdb, cif, pdb.gz, cif.gz
bool isCompressible(std::pair<std::string, std::string>& fileParts) {
    std::string ext = fileParts.second;
    if (ext == "pdb" || ext == "cif") {
        return true;
    } else if (ext == "gz") {
        std::pair<std::string, std::string> fileParts2 = getFileParts(fileParts.first);
        if (fileParts2.second == "pdb" || fileParts2.second == "cif") {
            return true;
        }
    }
    return false;
}

bool stringEndsWith(const std::string& suffix, const std::string& str) {
    if (str.length() < suffix.length()) {
        return false;
    }

    return (!str.compare(str.length() - suffix.length(), suffix.length(), suffix));
}

bool stringStartsWith(const std::string& prefix, const std::string& str, const size_t offset) {
    if (str.length() < prefix.length()) {
        return false;
    }
    return (!str.compare(offset, prefix.length(), prefix));
}

#if defined(_WIN32) || defined(_WIN64)
#define strtok_r strtok_s
#endif
std::vector<std::string> stringSplit(const std::string& str, const std::string& sep) {
    std::vector<std::string> arr;

    char* cstr = strdup(str.c_str());
    const char* csep = sep.c_str();
    char* rest;
    char* current = strtok_r(cstr, csep, &rest);
    while (current != NULL) {
        arr.emplace_back(current);
        current = strtok_r(NULL, csep, &rest);
    }
    free(cstr);

    return arr;
}

/*  FUNCTIONS FOR HANDLING AMINO ACID NAME  */

char getOneLetterCode(std::string threeLetterCode) {
    char aa;
    if (threeLetterCode == AA_ALA_STR) {
        aa = AA_ALA_CHAR;
    } else if (threeLetterCode == AA_ARG_STR) {
        aa = AA_ARG_CHAR;
    } else if (threeLetterCode == AA_ASN_STR) {
        aa = AA_ASN_CHAR;
    } else if (threeLetterCode == AA_ASP_STR) {
        aa = AA_ASP_CHAR;
    } else if (threeLetterCode == AA_CYS_STR) {
        aa = AA_CYS_CHAR;
    } else if (threeLetterCode == AA_GLN_STR) {
        aa = AA_GLN_CHAR;
    } else if (threeLetterCode == AA_GLU_STR) {
        aa = AA_GLU_CHAR;
    } else if (threeLetterCode == AA_GLY_STR) {
        aa = AA_GLY_CHAR;
    } else if (threeLetterCode == AA_HIS_STR) {
        aa = AA_HIS_CHAR;
    } else if (threeLetterCode == AA_ILE_STR) {
        aa = AA_ILE_CHAR;
    } else if (threeLetterCode == AA_LEU_STR) {
        aa = AA_LEU_CHAR;
    } else if (threeLetterCode == AA_LYS_STR) {
        aa = AA_LYS_CHAR;
    } else if (threeLetterCode == AA_MET_STR) {
        aa = AA_MET_CHAR;
    } else if (threeLetterCode == AA_PHE_STR) {
        aa = AA_PHE_CHAR;
    } else if (threeLetterCode == AA_PRO_STR) {
        aa = AA_PRO_CHAR;
    } else if (threeLetterCode == AA_SER_STR) {
        aa = AA_SER_CHAR;
    } else if (threeLetterCode == AA_THR_STR) {
        aa = AA_THR_CHAR;
    } else if (threeLetterCode == AA_TRP_STR) {
        aa = AA_TRP_CHAR;
    } else if (threeLetterCode == AA_TYR_STR) {
        aa = AA_TYR_CHAR;
    } else if (threeLetterCode == AA_VAL_STR) {
        aa = AA_VAL_CHAR;
    } else if (threeLetterCode == AA_ASX_STR) {
        aa = AA_ASX_CHAR;
    } else if (threeLetterCode == AA_GLX_STR) {
        aa = AA_GLX_CHAR;
    } else if (threeLetterCode == AA_STOP_STR) {
        aa = AA_STOP_CHAR;
    } else if (threeLetterCode == AA_UNK_STR) {
        aa = AA_UNK_CHAR;
    } else {
        aa = AA_UNK_CHAR;
    }
    return aa;
}

std::string getThreeLetterCode(char oneLetterCode) {
    std::string threeLetterCode;
    if (oneLetterCode == AA_ALA_CHAR) {
        threeLetterCode = AA_ALA_STR;
    } else if (oneLetterCode == AA_ARG_CHAR) {
        threeLetterCode = AA_ARG_STR;
    } else if (oneLetterCode == AA_ASN_CHAR) {
        threeLetterCode = AA_ASN_STR;
    } else if (oneLetterCode == AA_ASP_CHAR) {
        threeLetterCode = AA_ASP_STR;
    } else if (oneLetterCode == AA_CYS_CHAR) {
        threeLetterCode = AA_CYS_STR;
    } else if (oneLetterCode == AA_GLN_CHAR) {
        threeLetterCode = AA_GLN_STR;
    } else if (oneLetterCode == AA_GLU_CHAR) {
        threeLetterCode = AA_GLU_STR;
    } else if (oneLetterCode == AA_GLY_CHAR) {
        threeLetterCode = AA_GLY_STR;
    } else if (oneLetterCode == AA_HIS_CHAR) {
        threeLetterCode = AA_HIS_STR;
    } else if (oneLetterCode == AA_ILE_CHAR) {
        threeLetterCode = AA_ILE_STR;
    } else if (oneLetterCode == AA_LEU_CHAR) {
        threeLetterCode = AA_LEU_STR;
    } else if (oneLetterCode == AA_LYS_CHAR) {
        threeLetterCode = AA_LYS_STR;
    } else if (oneLetterCode == AA_MET_CHAR) {
        threeLetterCode = AA_MET_STR;
    } else if (oneLetterCode == AA_PHE_CHAR) {
        threeLetterCode = AA_PHE_STR;
    } else if (oneLetterCode == AA_PRO_CHAR) {
        threeLetterCode = AA_PRO_STR;
    } else if (oneLetterCode == AA_SER_CHAR) {
        threeLetterCode = AA_SER_STR;
    } else if   (oneLetterCode == AA_THR_CHAR) {
        threeLetterCode = AA_THR_STR;
    } else if (oneLetterCode == AA_TRP_CHAR) {
        threeLetterCode = AA_TRP_STR;
    } else if (oneLetterCode == AA_TYR_CHAR) {
        threeLetterCode = AA_TYR_STR;
    } else if (oneLetterCode == AA_VAL_CHAR) {
        threeLetterCode = AA_VAL_STR;
    } else if (oneLetterCode == AA_ASX_CHAR) {
        threeLetterCode = AA_ASX_STR;
    } else if (oneLetterCode == AA_GLX_CHAR) {
        threeLetterCode = AA_GLX_STR;
    } else if (oneLetterCode == AA_STOP_CHAR) {
        threeLetterCode = AA_STOP_STR;
    } else if (oneLetterCode == AA_UNK_CHAR) {
        threeLetterCode = AA_UNK_STR;
    } else {
        threeLetterCode = AA_UNK_STR;
    }
    return threeLetterCode;
}

/**
 * @brief Read 5-bit encoded residue (unsigned integer) and return
 *       corresponding amino acid (char)
 *
 * @param aab
 * @return char
 */
char convertIntToOneLetterCode(unsigned int aab) {
    char aa;
    switch (aab) {
        case AA_ALA_INT:
            aa = AA_ALA_CHAR;
            break;
        case AA_ARG_INT:
            aa = AA_ARG_CHAR;
            break;
        case AA_ASN_INT:
            aa = AA_ASN_CHAR;
            break;
        case AA_ASP_INT:
            aa = AA_ASP_CHAR;
            break;
        case AA_CYS_INT:
            aa = AA_CYS_CHAR;
            break;
        case AA_GLN_INT:
            aa = AA_GLN_CHAR;
            break;
        case AA_GLU_INT:
            aa = AA_GLU_CHAR;
            break;
        case AA_GLY_INT:
            aa = AA_GLY_CHAR;
            break;
        case AA_HIS_INT:
            aa = AA_HIS_CHAR;
            break;
        case AA_ILE_INT:
            aa = AA_ILE_CHAR;
            break;
        case AA_LEU_INT:
            aa = AA_LEU_CHAR;
            break;
        case AA_LYS_INT:
            aa = AA_LYS_CHAR;
            break;
        case AA_MET_INT:
            aa = AA_MET_CHAR;
            break;
        case AA_PHE_INT:
            aa = AA_PHE_CHAR;
            break;
        case AA_PRO_INT:
            aa = AA_PRO_CHAR;
            break;
        case AA_SER_INT:
            aa = AA_SER_CHAR;
            break;
        case AA_THR_INT:
            aa = AA_THR_CHAR;
            break;
        case AA_TRP_INT:
            aa = AA_TRP_CHAR;
            break;
        case AA_TYR_INT:
            aa = AA_TYR_CHAR;
            break;
        case AA_VAL_INT:
            aa = AA_VAL_CHAR;
            break;
        case AA_ASX_INT:
            aa = AA_ASX_CHAR;
            break;
        case AA_GLX_INT:
            aa = AA_GLX_CHAR;
            break;
        case AA_STOP_INT:
            aa = AA_STOP_CHAR;
            break;
        case AA_UNK_INT:
            aa = AA_UNK_CHAR;
            break;
        default:
            aa = AA_UNK_CHAR;
            break;
    }
    return aa;
}

unsigned int convertOneLetterCodeToInt(char oneLetterCode) {
    unsigned int output;
    switch (oneLetterCode) {
        case AA_ALA_CHAR:
            output = AA_ALA_INT;
            break;
        case AA_ARG_CHAR:
            output = AA_ARG_INT;
            break;
        case AA_ASN_CHAR:
            output = AA_ASN_INT;
            break;
        case AA_ASP_CHAR:
            output = AA_ASP_INT;
            break;
        case AA_CYS_CHAR:
            output = AA_CYS_INT;
            break;
        case AA_GLN_CHAR:
            output = AA_GLN_INT;
            break;
        case AA_GLU_CHAR:
            output = AA_GLU_INT;
            break;
        case AA_GLY_CHAR:
            output = AA_GLY_INT;
            break;
        case AA_HIS_CHAR:
            output = AA_HIS_INT;
            break;
        case AA_ILE_CHAR:
            output = AA_ILE_INT;
            break;
        case AA_LEU_CHAR:
            output = AA_LEU_INT;
            break;
        case AA_LYS_CHAR:
            output = AA_LYS_INT;
            break;
        case AA_MET_CHAR:
            output = AA_MET_INT;
            break;
        case AA_PHE_CHAR:
            output = AA_PHE_INT;
            break;
        case AA_PRO_CHAR:
            output = AA_PRO_INT;
            break;
        case AA_SER_CHAR:
            output = AA_SER_INT;
            break;
        case AA_THR_CHAR:
            output = AA_THR_INT;
            break;
        case AA_TRP_CHAR:
            output = AA_TRP_INT;
            break;
        case AA_TYR_CHAR:
            output = AA_TYR_INT;
            break;
        case AA_VAL_CHAR:
            output = AA_VAL_INT;
            break;
        case AA_ASX_CHAR:
            output = AA_ASX_INT;
            break;
        case AA_GLX_CHAR:
            output = AA_GLX_INT;
            break;
        case AA_STOP_CHAR:
            output = AA_STOP_INT;
            break;
        case AA_UNK_CHAR:
            output = AA_UNK_INT;
            break;
        default:
            output = AA_UNK_INT;
            break;
    }
    return output;
}

std::string convertIntToThreeLetterCode(unsigned int aab) {
    std::string threeLetterCode;
    switch (aab) {
        case AA_ALA_INT:
            threeLetterCode = AA_ALA_STR;
            break;
        case AA_ARG_INT:
            threeLetterCode = AA_ARG_STR;
            break;
        case AA_ASN_INT:
            threeLetterCode = AA_ASN_STR;
            break;
        case AA_ASP_INT:
            threeLetterCode = AA_ASP_STR;
            break;
        case AA_CYS_INT:
            threeLetterCode = AA_CYS_STR;
            break;
        case AA_GLN_INT:
            threeLetterCode = AA_GLN_STR;
            break;
        case AA_GLU_INT:
            threeLetterCode = AA_GLU_STR;
            break;
        case AA_GLY_INT:
            threeLetterCode = AA_GLY_STR;
            break;
        case AA_HIS_INT:
            threeLetterCode = AA_HIS_STR;
            break;
        case AA_ILE_INT:
            threeLetterCode = AA_ILE_STR;
            break;
        case AA_LEU_INT:
            threeLetterCode = AA_LEU_STR;
            break;
        case AA_LYS_INT:
            threeLetterCode = AA_LYS_STR;
            break;
        case AA_MET_INT:
            threeLetterCode = AA_MET_STR;
            break;
        case AA_PHE_INT:
            threeLetterCode = AA_PHE_STR;
            break;
        case AA_PRO_INT:
            threeLetterCode = AA_PRO_STR;
            break;
        case AA_SER_INT:
            threeLetterCode = AA_SER_STR;
            break;
        case AA_THR_INT:
            threeLetterCode = AA_THR_STR;
            break;
        case AA_TRP_INT:
            threeLetterCode = AA_TRP_STR;
            break;
        case AA_TYR_INT:
            threeLetterCode = AA_TYR_STR;
            break;
        case AA_VAL_INT:
            threeLetterCode = AA_VAL_STR;
            break;
        case AA_ASX_INT:
            threeLetterCode = AA_ASX_STR;
            break;
        case AA_GLX_INT:
            threeLetterCode = AA_GLX_STR;
            break;
        case AA_STOP_INT:
            threeLetterCode = AA_STOP_STR;
            break;
        case AA_UNK_INT:
            threeLetterCode = AA_UNK_STR;
            break;
        default:
            threeLetterCode = AA_UNK_STR;
            break;
    }
    return threeLetterCode;
}

unsigned int convertThreeLetterCodeToInt(std::string threeLetterCode) {
    unsigned int output;
    if (threeLetterCode == AA_ALA_STR) {
        output = AA_ALA_INT;
    } else if (threeLetterCode == AA_ARG_STR) {
        output = AA_ARG_INT;
    } else if (threeLetterCode == AA_ASN_STR) {
        output = AA_ASN_INT;
    } else if (threeLetterCode == AA_ASP_STR) {
        output = AA_ASP_INT;
    } else if (threeLetterCode == AA_CYS_STR) {
        output = AA_CYS_INT;
    } else if (threeLetterCode == AA_GLN_STR) {
        output = AA_GLN_INT;
    } else if (threeLetterCode == AA_GLU_STR) {
        output = AA_GLU_INT;
    } else if (threeLetterCode == AA_GLY_STR) {
        output = AA_GLY_INT;
    } else if (threeLetterCode == AA_HIS_STR) {
        output = AA_HIS_INT;
    } else if (threeLetterCode == AA_ILE_STR) {
        output = AA_ILE_INT;
    } else if (threeLetterCode == AA_LEU_STR) {
        output = AA_LEU_INT;
    } else if (threeLetterCode == AA_LYS_STR) {
        output = AA_LYS_INT;
    } else if (threeLetterCode == AA_MET_STR) {
        output = AA_MET_INT;
    } else if (threeLetterCode == AA_PHE_STR) {
        output = AA_PHE_INT;
    } else if (threeLetterCode == AA_PRO_STR) {
        output = AA_PRO_INT;
    } else if (threeLetterCode == AA_SER_STR) {
        output = AA_SER_INT;
    } else if (threeLetterCode == AA_THR_STR) {
        output = AA_THR_INT;
    } else if (threeLetterCode == AA_TRP_STR) {
        output = AA_TRP_INT;
    } else if (threeLetterCode == AA_TYR_STR) {
        output = AA_TYR_INT;
    } else if (threeLetterCode == AA_VAL_STR) {
        output = AA_VAL_INT;
    } else if (threeLetterCode == AA_ASX_STR) {
        output = AA_ASX_INT;
    } else if (threeLetterCode == AA_GLX_STR) {
        output = AA_GLX_INT;
    } else if (threeLetterCode == AA_STOP_STR) {
        output = AA_STOP_INT;
    } else if (threeLetterCode == AA_UNK_STR) {
        output = AA_UNK_INT;
    } else {
        output = AA_UNK_INT;
    }
    return output;
}
