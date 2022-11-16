/**
 * File: main.cpp
 * Project: foldcomp
 * Created: 2021-12-23 17:44:53
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Contributor: Milot Mirdita (milot@mirdita.de), Martin Steinegger (themartinsteinegger@gmail.com)
 * Description:
 *     This code contains main function for "foldcomp".
 *     Foldcomp is a fast lossy compression algorithm for protein structure.
 *     It encodes torsion angles with optimal number of bits and reconstruct
 *     3D coordinates from the encoded angles.
 * Usage:
 *    foldcomp compress input.pdb output.fcz
 *    foldcomp decompress input.fcz output.pdb
 * ---
 * Last Modified: 2022-10-07 17:25:34
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
// Headers in the project
#include "atom_coordinate.h"
#include "foldcomp.h"
#include "structure_reader.h"
#include "utility.h"

// Standard libraries
#include <cstring>
#include <fstream> // IWYU pragma: keep
#include <getopt.h>
#include <iostream>
#include <sstream> // IWYU pragma: keep
#include <string>
#include <vector>

#include <sys/stat.h>

// OpenMP for parallelization
#ifdef OPENMP
#include <omp.h>
#endif

#ifdef HAVE_GCS
#include "google/cloud/storage/client.h"
#endif

static int use_alt_order = 0;
static int anchor_residue_threshold = DEFAULT_ANCHOR_THRESHOLD;
static int save_as_tar = 0;
static int ext_mode = 0;
static int ext_merge = 1;

int print_usage(void) {
    std::cout << "Usage: foldcomp compress <pdb_file> [<fcz_file>]" << std::endl;
    std::cout << "       foldcomp compress [-t number] <pdb_dir> [<fcz_dir>]" << std::endl;
    std::cout << "       foldcomp decompress <fcz_file|tar> [<pdb_file>]" << std::endl;
    std::cout << "       foldcomp decompress [-t number] <fcz_dir|tar> [<pdb_dir>]" << std::endl;
    std::cout << "       foldcomp extract [--plddt|--amino-acid] <fcz_file> [<fasta_file>]" << std::endl;
    std::cout << "       foldcomp extract [--plddt|--amino-acid] [-t number] <fcz_dir|tar> [<fasta_dir>]" << std::endl;
    std::cout << "       foldcomp check <fcz_file>" << std::endl;
    std::cout << "       foldcomp check [-t number] <fcz_dir|tar>" << std::endl;
    std::cout << "       foldcomp rmsd <pdb1|cif1> <pdb2|cif2>" << std::endl;
    std::cout << " -h, --help           print this help message" << std::endl;
    std::cout << " -t, --threads        threads for (de)compression of folders/tar files [default=1]" << std::endl;
    std::cout << " -a, --alt            use alternative atom order [default=false]" << std::endl;
    std::cout << " -b, --break          interval size to save absolute atom coordinates [default=" << anchor_residue_threshold << "]" << std::endl;
    std::cout << " -z, --tar            save as tar file [default=false]" << std::endl;
    std::cout << " --plddt              extract pLDDT score (only for extraction mode)" << std::endl;
    std::cout << " --fasta              extract amino acid sequence (only for extraction mode)" << std::endl;
    std::cout << " --no-merge           do not merge output files (only for extraction mode)" << std::endl;
    return 0;
}

int compress(std::string input, std::string output) {
    StructureReader reader;
    reader.load(input);
    std::vector<AtomCoordinate> atomCoordinates;
    reader.readAllAtoms(atomCoordinates);
    if (atomCoordinates.size() == 0) {
        std::cout << "[Error] No atoms found in the input file: " << input << std::endl;
        return 1;
    }
    std::string title = reader.title;

    std::vector<BackboneChain> compData;
    Foldcomp compRes = Foldcomp();
    // Convert title to char
    compRes.strTitle = title;
    compRes.anchorThreshold = anchor_residue_threshold;
    compData = compRes.compress(atomCoordinates);
    // Write compressed data to file
    if (compRes.write(output) != 0) {
        std::cout << "[Error] Writing file: " << output << std::endl;
        return -1;
    }
    // DEBUGGING
    // Nerf nerf;
    // nerf.writeInfoForChecking(atomCoordinates, "BEFORE_COMPRESSION.csv");
    // clear memory
    // atomCoordinates.clear();
    // compData.clear();
    return 0;
}


int compressFromBuffer(const std::string& content, const std::string& output, std::string& name) {
    StructureReader reader;
    reader.loadFromBuffer(content.c_str(), content.size(), name);
    std::vector<AtomCoordinate> atomCoordinates;
    reader.readAllAtoms(atomCoordinates);
    if (atomCoordinates.size() == 0) {
        std::cout << "[Error] No atoms found in the input" << std::endl;
        return 1;
    }
    std::string title = reader.title;

    std::vector<BackboneChain> compData;
    Foldcomp compRes = Foldcomp();
    // Convert title to char
    compRes.strTitle = name;
    compRes.anchorThreshold = anchor_residue_threshold;
    compData = compRes.compress(atomCoordinates);
    // Write compressed data to file
    compRes.write(output + "/" + name + ".fcz");
    // DEBUGGING
    // Nerf nerf;
    // nerf.writeInfoForChecking(atomCoordinates, "BEFORE_COMPRESSION.csv");
    // clear memory
    // atomCoordinates.clear();
    // compData.clear();
    return 0;
}


int compressWithoutWriting(Foldcomp& compRes, std::string input) {
    StructureReader reader;
    reader.load(input);
    std::vector<AtomCoordinate> atomCoordinates;
    reader.readAllAtoms(atomCoordinates);
    if (atomCoordinates.size() == 0) {
        std::cout << "[Error] No atoms found in the input file: " << input << std::endl;
        return 1;
    }
    std::string title = reader.title;

    std::vector<BackboneChain> compData;
    // Convert title to char
    compRes.strTitle = title;
    compRes.anchorThreshold = anchor_residue_threshold;
    compData = compRes.compress(atomCoordinates);
    return 0;
}

int compressFromBufferWithoutWriting(Foldcomp& compRes, const std::string& content, std::string& name) {
    StructureReader reader;
    reader.loadFromBuffer(content.c_str(), content.size(), name);
    std::vector<AtomCoordinate> atomCoordinates;
    reader.readAllAtoms(atomCoordinates);
    if (atomCoordinates.size() == 0) {
        std::cout << "[Error] No atoms found in the input" << std::endl;
        return 1;
    }
    std::string title = name;
    std::vector<BackboneChain> compData;
    // Convert title to char
    compRes.strTitle = name;
    compRes.anchorThreshold = anchor_residue_threshold;
    compData = compRes.compress(atomCoordinates);
    return 0;
}

int decompress(std::istream &file, std::string output) {
    int flag = 0;
    Foldcomp compRes = Foldcomp();
    flag = compRes.read(file);
    if (flag != 0) {
        std::cerr << "[Error] Reading" << std::endl;
        return 1;
    }
    std::vector<AtomCoordinate> atomCoordinates;
    compRes.useAltAtomOrder = use_alt_order;
    flag = compRes.decompress(atomCoordinates);
    if (flag != 0) {
        std::cerr << "[Error] decompressing compressed data." << std::endl;
        return 1;
    }
    // Write decompressed data to file
    flag = writeAtomCoordinatesToPDBFile(atomCoordinates, compRes.strTitle, output);

    return flag;
}

int extract(std::istream& file, std::string output) {
    int flag = 0;
    Foldcomp compRes = Foldcomp();
    flag = compRes.read(file);
    if (flag != 0) {
        std::cerr << "[Error] reading" << std::endl;
        return 1;
    }
    std::vector<std::string> data;
    compRes.extract(data, ext_mode);
    compRes.writeFASTALike(output, data);
    return 0;
}

int check(std::istream& file, std::string& filename) {
    int flag = 0;
    Foldcomp compRes = Foldcomp();
    flag = compRes.read(file);
    if (flag != 0) {
        std::cerr << "[Error] reading file: " << filename << std::endl;
        return 1;
    }
    ValidityError err;
    err = compRes.checkValidity();
    printValidityError(err, filename);
    return flag;
}

int rmsd(std::string pdb1, std::string pdb2) {
    // Read
    StructureReader reader1;
    reader1.load(pdb1);
    std::vector<AtomCoordinate> atomCoordinates1;
    reader1.readAllAtoms(atomCoordinates1);
    StructureReader reader2;
    reader2.load(pdb2);
    std::vector<AtomCoordinate> atomCoordinates2;
    reader2.readAllAtoms(atomCoordinates2);
    // Check
    if (atomCoordinates1.size() == 0) {
        std::cerr << "[Error] No atoms found in the input file: " << pdb1 << std::endl;
        return 1;
    }
    if (atomCoordinates2.size() == 0) {
        std::cerr << "[Error] No atoms found in the input file: " << pdb2 << std::endl;
        return 1;
    }
    if (atomCoordinates1.size() != atomCoordinates2.size()) {
        std::cerr << "[Error] The number of atoms in the two files are different." << std::endl;
        return 1;
    }
    std::vector<AtomCoordinate> backbone1 = filterBackbone(atomCoordinates1);
    std::vector<AtomCoordinate> backbone2 = filterBackbone(atomCoordinates2);
    // Print
    std::cout << pdb1 << '\t' << pdb2 << '\t';
    std::cout << backbone1.size() / 3 << '\t' << atomCoordinates1.size() << '\t';
    std::cout << RMSD(backbone1, backbone2) << '\t';
    std::cout << RMSD(atomCoordinates1, atomCoordinates2) << std::endl;
    return 0;
}

int main(int argc, char* const *argv) {
    if (argc < 3) {
        return print_usage();
    }

    int flag = 0;
    int option_index = 0;
    int num_threads = 1;
    int has_output = 0;

    // Mode - non-optional argument
    enum {
        COMPRESS,
        DECOMPRESS,
        COMPRESS_MULTIPLE,
        COMPRESS_MULTIPLE_GCS,
        DECOMPRESS_MULTIPLE,
        DECOMPRESS_MULTIPLE_TAR,
        EXTRACT,
        EXTRACT_MULTIPLE,
        EXTRACT_MULTIPLE_TAR,
        CHECK,
        CHECK_MULTIPLE,
        CHECK_MULTIPLE_TAR,
        RMSD
    } mode = COMPRESS;

    // Define command line options
    static struct option long_options[] = {
            {"help",          no_argument,          0, 'h'},
            {"alt",           no_argument,          0, 'a'},
            {"tar",           no_argument,          0, 'z'},
            {"plddt",         no_argument,  &ext_mode,  0 },
            {"fasta",         no_argument,  &ext_mode,  1 },
            {"no-merge",      no_argument, &ext_merge,  0 },
            {"threads", required_argument,          0, 't'},
            {"break",   required_argument,          0, 'b'},
            {0,                         0,          0,  0 }
    };

    // Parse command line options with getopt_long
    flag = getopt_long(argc, argv, "hazt:b:c:", long_options, &option_index);

    while (flag != -1) {
        switch (flag) {
            case 'h':
                return print_usage();
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'a':
                use_alt_order = 1;
                break;
            case 'z':
                save_as_tar = 1;
                break;
            case 'b':
                anchor_residue_threshold = atoi(optarg);
                break;
            case '?':
                return print_usage();
            default:
                break;
        }
        flag = getopt_long(argc, argv, "hazt:b:", long_options, &option_index);
    }

    // Parse non-option arguments
    // argv[optind]: MODE
    // argv[optind + 1]: INPUT
    // argv[optind + 2]: OUTPUT (optional)

    if ((optind + 1) >= argc) {
        std::cerr << "[Error] Not enough arguments." << std::endl;
        return print_usage();
    }

    struct stat st;
    int fileExists = stat(argv[optind + 1], &st);
    // get mode from command line
    if (strcmp(argv[optind], "compress") == 0) {
        // Check argv[2] is file, directory, or gcs URI
        // TODO: stdin support
        // TODO: COMPRESS_MULTIPLE_TAR
        // If gcs URI, mode = COMPRESS_MULTIPLE_GCS
        // If directory, mode = COMPRESS_MULTIPLE
        // If file, mode = COMPRESS
#ifdef HAVE_GCS
        if ((optind + 1) < argc && stringStartsWith("gcs://", argv[optind + 1])) {
            mode = COMPRESS_MULTIPLE_GCS;
            fileExists = 0;
        } else
#endif
        if (st.st_mode & S_ISDIR(st.st_mode)) {
            mode = COMPRESS_MULTIPLE;
        } else {
            mode = COMPRESS;
        }
    } else if (strcmp(argv[optind], "decompress") == 0) {
        // Check argv[2] is file or directory
        // If directory, mode = DECOMPRESS_MULTIPLE
        // If file, mode = DECOMPRESS
        char *end = strrchr(argv[optind + 1], '.');
        if (st.st_mode & S_ISDIR(st.st_mode)) {
            mode = DECOMPRESS_MULTIPLE;
        } else if (strcmp(end, ".tar") == 0) {
            mode = DECOMPRESS_MULTIPLE_TAR;
        } else {
            mode = DECOMPRESS;
        }
    } else if (strcmp(argv[optind], "extract") == 0) {
        // Check argv[2] is file or directory
        // If directory, mode = EXTRACT_MULTIPLE
        // If file, mode = EXTRACT
        char* end = strrchr(argv[optind + 1], '.');
        if (st.st_mode & S_ISDIR(st.st_mode)) {
            mode = EXTRACT_MULTIPLE;
        } else if (strcmp(end, ".tar") == 0) {
            mode = EXTRACT_MULTIPLE_TAR;
        } else {
            mode = EXTRACT;
        }
    } else if (strcmp(argv[optind], "check") == 0){
        char* end = strrchr(argv[optind + 1], '.');
        if (st.st_mode & S_ISDIR(st.st_mode)) {
            mode = CHECK_MULTIPLE;
        }
        else if (strcmp(end, ".tar") == 0) {
            mode = CHECK_MULTIPLE_TAR;
        }
        else {
            mode = CHECK;
        }
    } else if (strcmp(argv[optind], "rmsd") == 0) {
        // NO MULTIPLE MODE FOR RMSD
        mode = RMSD;
    } else {
        return print_usage();
    }
    // Error if no input file given
    if (mode != COMPRESS_MULTIPLE_GCS && fileExists == -1) {
        std::cerr << "[Error] " << argv[optind + 1] << " does not exist." << std::endl;
        return 1;
    }

    std::string input = argv[optind + 1];
    std::string output;
    if (argc == optind + 3) {
        has_output = 1;
        output = argv[optind + 2];
    }

    // check if mode is compress or decompress
    if (mode == COMPRESS) {
        // compress a single file
        if (!has_output) {
            output = getFileWithoutExt(input) + ".fcz";
        }
        std::cout << "Compressing " << input << " to " << output << std::endl;
        compress(input, output);
        flag = 0;
    } else if (mode == DECOMPRESS) {
        // decompress a single file
        if (!has_output) {
            output = getFileWithoutExt(input) + "_fcz.pdb";
        }
        std::ifstream inputFile(input, std::ios::binary);
        // Check if file is open
        if (!inputFile.is_open()) {
            std::cout << "[Error] Could not open file " << input << std::endl;
            return -1;
        }

        std::cout << "Decompressing " << input << " to " << output << std::endl;
        decompress(inputFile, output);
        inputFile.close();
        flag = 0;
    } else if (mode == EXTRACT) {
        // In extract mode, specific information is directly extracted from the fcz file
        if (!has_output) {
            if (ext_mode == 0) {
                output = getFileWithoutExt(input) + ".plddt.txt";
                std::cout << "Extracting PLDDT from " << input << " to " << output << std::endl;
            } else if (ext_mode == 1) {
                output = getFileWithoutExt(input) + ".fasta";
                std::cout << "Extracting amino acid sequence from " << input << " to " << output << std::endl;
            }
        }
        std::ifstream inputFile(input, std::ios::binary);
        // Check if file is open
        if (!inputFile.is_open()) {
            std::cout << "[Error] Could not open file " << input << std::endl;
            return -1;
        }
        extract(inputFile, output);
        inputFile.close();
        flag = 0;
    } else if (mode == CHECK){
        // Check if the file is a valid fcz file
        std::ifstream inputFile(input, std::ios::binary);
        std::clog << "Checking " << input << std::endl;
        if (!inputFile.is_open()) {
            std::cerr << "[Error] Could not open file " << input << std::endl;
            return -1;
        }
        flag = check(inputFile, input);
        inputFile.close();
    } else if (mode == RMSD) {
        // Calculate RMSD between two PDB files
        rmsd(input, output);
        flag = 0;
    } else if (mode == COMPRESS_MULTIPLE) {
        // compress multiple files
        if (input[input.length() - 1] != '/') {
            input += "/";
        }
        if (!has_output) {
            output = input.substr(0, input.length() - 1) + "_fcz/";
        } else {
            if (stringEndsWith(output, ".tar")) {
                save_as_tar = 1;
            }
        }

        if (output[output.length() - 1] != '/' && !save_as_tar) {
            output += "/";
        }
        // Check output directory exists or not
        if (!save_as_tar) {
            if (stat(output.c_str(), &st) == -1) {
#if defined(_WIN32) || defined(_WIN64)
                _mkdir(output.c_str());
#else
                mkdir(output.c_str(), 0755);
#endif
            }
        }
        // Get all files in input directory
        std::string file;
        std::string inputFile;
        std::string outputFile;
        std::cout << "Compressing files in " << input;
        std::cout << " using " << num_threads << " threads" << std::endl;
        if (save_as_tar) {
            std::cout << "Output tar file: " << output << std::endl;
        } else {
            std::cout << "Output directory: " << output << std::endl;
        }
        std::vector<std::string> files = getFilesInDirectory(input);
        // Parallelize
        if (!save_as_tar) {
#pragma omp parallel num_threads(num_threads)
            {
#pragma omp for
                for (size_t i = 0; i < files.size(); i++) {
                    std::string file = files[i];
                    std::string inputFile = input + file;
                    std::string outputFile = output + getFileWithoutExt(file) + ".fcz";
                    compress(inputFile, outputFile);
                }
            }
        } else {
            mtar_t tar;
            std::string tarFile = output.substr(0, output.length() - 1) + ".tar";
            mtar_open(&tar, tarFile.c_str(), "w");
#pragma omp parallel num_threads(num_threads)
            {
#pragma omp for
                for (size_t i = 0; i < files.size(); i++) {
                    std::string file = files[i];
                    std::string inputFile = input + file;
                    std::string outputFile = output + getFileWithoutExt(file) + ".fcz";
                    Foldcomp compRes = Foldcomp();
                    compressWithoutWriting(compRes, inputFile);
#pragma omp critical
                {
                    compRes.writeTar(tar, outputFile, compRes.getSize());
                }
                }
            }
            mtar_finalize(&tar);
            mtar_close(&tar);
        }
    } else if (mode == COMPRESS_MULTIPLE_GCS) {
        // compress multiple files from gcs
#ifdef HAVE_GCS
        if (!has_output) {
            std::cerr << "Please specify output directory" << std::endl;
            return 1;
        }
        if (output[output.length() - 1] != '/') {
            output += "/";
        }
        // Check output directory exists or not
        if (stat(output.c_str(), &st) == -1) {
#if defined(_WIN32) || defined(_WIN64)
            _mkdir(output.c_str());
#else
            mkdir(output.c_str(), 0755);
#endif
        }
        // Get all files in input directory
        std::cout << "Compressing files in " << input;
        std::cout << " using " << num_threads << " threads" << std::endl;
        std::cout << "Output directory: " << output << std::endl;
        // Parallelize
        namespace gcs = ::google::cloud::storage;
        auto options = google::cloud::Options{}
            .set<gcs::ConnectionPoolSizeOption>(num_threads)
            .set<google::cloud::storage_experimental::HttpVersionOption>("2.0");
        auto client = gcs::Client(options);
        std::vector<std::string> parts = stringSplit(input, "/");
        if (parts.size() == 1) {
            std::cerr << "Invalid gcs URI" << std::endl;
            return 1;
        }
        std::string bucket_name = parts[1];

        // Filter for splitting input into 10 different processes
        //char filter = parts[2][0];
        int num_tar = num_threads;
        mtar_t tarArray[num_tar];
        std::vector<std::string> tarFiles;
        for (int i = 0; i < num_tar; i++) {
            std::string tarFile = output + "AF2_Uniprot_foldcomp." + std::to_string(i) + ".tar";
            tarFiles.push_back(tarFile);
            mtar_open(&tarArray[i], tarFile.c_str(), "w");
        }

#pragma omp parallel num_threads(num_threads)
        {
#pragma omp single
            // Get object list from gcs bucket
            for (auto&& object_metadata : client.ListObjects(bucket_name, gcs::Projection::NoAcl(), gcs::MaxResults(100000))) {
                std::string obj_name = object_metadata->name();
                // Set zero padding for ID with 4 digits
#pragma omp task firstprivate(obj_name)
                {
                    // Filter for splitting input into 10 different processes
                    // bool skipFilter = filter != '\0' && obj_name.length() >= 9 && obj_name[8] == filter;
                    bool skipFilter = true;
                    bool allowedSuffix = stringEndsWith(".cif", obj_name) || stringEndsWith(".pdb", obj_name);
                    if (skipFilter && allowedSuffix) {
                        auto reader = client.ReadObject(bucket_name, obj_name);
                        if (!reader.status().ok()) {
                            std::cerr << "Could not read object " << obj_name << std::endl;
                        } else {

                            std::string contents{ std::istreambuf_iterator<char>{reader}, {} };
                            Foldcomp compRes = Foldcomp();
                            std::string outputFile = output + getFileWithoutExt(obj_name) + ".fcz";
                            compressFromBufferWithoutWriting(compRes, contents, obj_name);
                            //compRes.writeTar(tar, outputFile, compRes.getSize());
                            int tar_id = 0;
#ifdef OPENMP
                            tar_id = omp_get_thread_num();
#endif

                            compRes.writeTar(tarArray[tar_id], outputFile, compRes.getSize());
                        }
                    }

                }
            }
            // Close tar
#pragma omp taskwait
            for (int i = 0; i < num_tar; i++) {
                mtar_finalize(&tarArray[i]);
                mtar_close(&tarArray[i]);
            }

        }
#endif
        flag = 0;
    } else if (mode == DECOMPRESS_MULTIPLE || mode == DECOMPRESS_MULTIPLE_TAR) {
        if (!has_output) {
            output = input.substr(0, input.length() - 1) + "_pdb/";
        }
        if (output[output.length() - 1] != '/') {
            output += "/";
        }
        // Check output directory exists or not
        if (stat(output.c_str(), &st) == -1) {
#if defined(_WIN32) || defined(_WIN64)
            _mkdir(output.c_str());
#else
            mkdir(output.c_str(), 0755);
#endif
        }
        if (mode == DECOMPRESS_MULTIPLE) {
            // decompress multiple files
            if (input[input.length() - 1] != '/') {
                input += "/";
            }

            // Get all files in input directory
            std::cout << "Decompressing files in " << input;
            std::cout << " using " << num_threads << " threads" << std::endl;
            std::cout << "Output directory: " << output << std::endl;
            std::vector<std::string> files = getFilesInDirectory(input);
#pragma omp parallel num_threads(num_threads)
            {
#pragma omp for
                for (size_t i = 0; i < files.size(); i++) {
                    std::string inputFile = input + files[i];
                    std::ifstream input(inputFile, std::ios::binary);
                    // Check if file is open
                    if (!input.is_open()) {
                        std::cout << "[Error] Could not open file " << inputFile << std::endl;
                        continue;
                    }
                    std::string outputFile = output + getFileWithoutExt(files[i]) + ".pdb";
                    decompress(input, outputFile);
                    input.close();
                }
            }
            flag = 0;
        } else if (mode == DECOMPRESS_MULTIPLE_TAR) {
            mtar_t tar;
            if (mtar_open(&tar, argv[optind + 1], "r") != MTAR_ESUCCESS) {
                std::cerr << "[Error] open tar " << argv[optind + 1] << " failed." << std::endl;
                return 1;
            }
#pragma omp parallel shared(tar) num_threads(num_threads)
            {
                bool proceed = true;
                mtar_header_t header;
                size_t bufferSize = 1024 * 1024;
                char *dataBuffer = (char *) malloc(bufferSize);
                std::string name;
                while (proceed) {
                    bool writeEntry = true;
#pragma omp critical
                    {
                        if (mtar_read_header(&tar, &header) != MTAR_ENULLRECORD) {
                            //TODO GNU tar has special blocks for long filenames
                            name = header.name;
                            if (header.size > bufferSize) {
                                bufferSize = header.size * 1.5;
                                dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                            }
                            if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                                std::cerr << "[Error] reading tar entry " << name << " failed." << std::endl;
                                writeEntry = false;
                                proceed = false;
                            } else {
                                writeEntry = true;
                                proceed = true;
                            }
                            mtar_next(&tar);
                            writeEntry = (header.type == MTAR_TREG) ? writeEntry : false;
                        } else {
                            proceed = false;
                            writeEntry = false;
                        }
                    } // end read in
                    if (proceed && writeEntry) {
                        std::istringstream input(std::string(dataBuffer, header.size));
                        std::string name_clean = name.substr(name.find_last_of("/\\") + 1);
                        std::string outputFile = output + name_clean + ".pdb";
                        decompress(input, outputFile);
                    }
                } // end while loop
            } // end openmp
            flag = 0;
        }
    } else if (mode == EXTRACT_MULTIPLE || mode == EXTRACT_MULTIPLE_TAR) {
            // extract multiple files
            if (mode == EXTRACT_MULTIPLE) {
                if (input[input.length() - 1] != '/') {
                    input += "/";
                }
                if (!has_output) {
                    if (ext_mode == 0) {
                        output = input.substr(0, input.length() - 1) + "_plddt/";
                    }
                    else if (ext_mode == 1) {
                        output = input.substr(0, input.length() - 1) + "_fasta/";
                    }
                }
                if (output[output.length() - 1] != '/') {
                    output += "/";
                }
                // Check output directory exists or not
                if (stat(output.c_str(), &st) == -1) {
#if defined(_WIN32) || defined(_WIN64)
                    _mkdir(output.c_str());
#else
                    mkdir(output.c_str(), 0755);
#endif
                }
                // Get all files in input directory
                std::string defaultOutputFile = "";
                if (ext_mode == 0){
                    defaultOutputFile = output + "plddt.txt";
                } else if (ext_mode == 1) {
                    defaultOutputFile = output + "aa.fasta";
                }
                std::ofstream defaultOutput(defaultOutputFile, std::ios::out);
                if (ext_merge == 0) {
                    // close and delete defaultOutput
                    defaultOutput.close();
                }
                std::cout << "Extracting files in " << input;
                std::cout << " using " << num_threads << " threads" << std::endl;
                std::cout << "Output directory: " << output << std::endl;
                std::vector<std::string> files = getFilesInDirectory(input);
#pragma omp parallel num_threads(num_threads)
                {
#pragma omp for
                    for (size_t i = 0; i < files.size(); i++) {
                        std::string inputFile = input + files[i];
                        std::ifstream input(inputFile, std::ios::binary);
                        // Check if file is open
                        if (!input.is_open()) {
                            std::cout << "[Error] Could not open file " << inputFile << std::endl;
                            continue;
                        }
                        std::string outputFile;
                        if (ext_merge == 1) {
                            std::vector<std::string> data;
                            Foldcomp compRes = Foldcomp();
                            compRes.read(input);
                            compRes.extract(data, ext_mode);
                            #pragma omp critical
                            {
                                defaultOutput << ">" << compRes.strTitle << "\n";
                                for (size_t j = 0; j < data.size(); j++) {
                                    defaultOutput << data[j];
                                }
                                defaultOutput << "\n";
                            }
                        } else {
                            if (ext_mode == 0) {
                                // output file extension is ".plddt.txt"
                                outputFile = output + getFileWithoutExt(files[i]) + ".plddt.txt";
                            }
                            else if (ext_mode == 1) {
                                outputFile = output + getFileWithoutExt(files[i]) + ".fasta";
                            }
                            extract(input, outputFile);
                        }
                        input.close();
                    }
                }
                flag = 0;
                if (ext_merge == 1) {
                    // close and delete defaultOutput
                    defaultOutput.close();
                }
            } else if (mode == EXTRACT_MULTIPLE_TAR) {
                mtar_t tar;
                if (mtar_open(&tar, input.c_str(), "r") != MTAR_ESUCCESS) {
                    std::cerr << "[Error] Open tar " << input << " failed." << std::endl;
                    return 1;
                }
                std::string defaultOutputFile = "";
                if (!has_output) {
                    if (ext_mode == 0) {
                        output = getFileWithoutExt(input) + "_plddt/";
                        defaultOutputFile = output + "plddt.txt";
                    }
                    else if (ext_mode == 1) {
                        output = getFileWithoutExt(input) + "_fasta/";
                        defaultOutputFile = output + "aa.fasta";
                    }
                }

                // Check output directory exists or not
                if (stat(output.c_str(), &st) == -1) {
#if defined(_WIN32) || defined(_WIN64)
                    _mkdir(output.c_str());
#else
                    mkdir(output.c_str(), 0755);
#endif
                }
                std::ofstream defaultOutput(defaultOutputFile, std::ios::out);
                if (ext_merge == 0) {
                    // close and delete defaultOutput
                    defaultOutput.close();
                }
                std::cout << "Extracting files in " << input << " using " << num_threads << " threads" << std::endl;
                // TAR READING PART BY MARTIN STEINEGGER
#pragma omp parallel shared(tar) num_threads(num_threads)
                {
                    bool proceed = true;
                    mtar_header_t header;
                    size_t bufferSize = 1024 * 1024;
                    char* dataBuffer = (char*)malloc(bufferSize);
                    std::string name;
                    while (proceed) {
                        bool writeEntry = true;
#pragma omp critical
                        {
                            if (mtar_read_header(&tar, &header) != MTAR_ENULLRECORD) {
                                //TODO GNU tar has special blocks for long filenames
                                name = header.name;
                                if (header.size > bufferSize) {
                                    bufferSize = header.size * 1.5;
                                    dataBuffer = (char*)realloc(dataBuffer, bufferSize);
                                }
                                if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                                    std::cerr << "[Error] Reading tar entry " << name << " failed." << std::endl;
                                    writeEntry = false;
                                    proceed = false;
                                }
                                else {
                                    writeEntry = true;
                                    proceed = true;
                                }
                                mtar_next(&tar);
                                writeEntry = (header.type == MTAR_TREG) ? writeEntry : false;
                            }
                            else {
                                proceed = false;
                                writeEntry = false;
                            }
                        } // end read in
                        if (proceed && writeEntry) {
                            std::istringstream input(std::string(dataBuffer, header.size));
                            std::string name_clean = name.substr(name.find_last_of("/\\") + 1);
                            std::string outputFile = "";

                            if (ext_merge == 1) {
                                std::vector<std::string> data;
                                Foldcomp compRes = Foldcomp();
                                compRes.read(input);
                                compRes.extract(data, ext_mode);
                            #pragma omp critical
                                {
                                    defaultOutput << ">" << compRes.strTitle << "\n";
                                    for (size_t j = 0; j < data.size(); j++) {
                                        defaultOutput << data[j];
                                    }
                                    defaultOutput << "\n";
                                }
                            } else {
                                if (ext_mode == 0) {
                                    // output file extension is ".plddt.txt"
                                    outputFile = output + name_clean + ".plddt.txt";
                                }
                                else if (ext_mode == 1) {
                                    outputFile = output + name_clean + ".fasta";
                                }
                                extract(input, outputFile);
                            }
                        }
                    } // end while loop
                } // end openmp
                flag = 0;
                if (ext_merge == 1) {
                    // close and delete defaultOutput
                    defaultOutput.close();
                }
            }
    } else if (mode == CHECK_MULTIPLE) {
        if (input[input.length() - 1] != '/') {
            input += "/";
        }
        std::vector<std::string> files = getFilesInDirectory(input);
        std::clog << "Checking files in " << input << " using " << num_threads << " threads" << std::endl;
#pragma omp parallel num_threads(num_threads)
        {
#pragma omp for
            for (size_t i = 0; i < files.size(); i++) {
                std::string inputFile = input + files[i];
                std::ifstream input(inputFile, std::ios::binary);
                // Check if file is open
                if (!input.is_open()) {
                    std::cerr << "[Error] Could not open file " << inputFile << std::endl;
                    continue;
                }
                check(input, inputFile);
                input.close();
            }
        }
        flag = 0;
    } else if (mode == CHECK_MULTIPLE_TAR) {
        mtar_t tar;
        if (mtar_open(&tar, input.c_str(), "r") != MTAR_ESUCCESS) {
            std::cerr << "[Error] open tar " << input << " failed." << std::endl;
            return 1;
        }
        std::clog << "Checking files in " << input << " using " << num_threads << " threads" << std::endl;
        // TAR READING PART BY MARTIN STEINEGGER
#pragma omp parallel shared(tar) num_threads(num_threads)
        {
            bool proceed = true;
            mtar_header_t header;
            size_t bufferSize = 1024 * 1024;
            char* dataBuffer = (char*)malloc(bufferSize);
            std::string name;
            while (proceed) {
                bool writeEntry = true;
#pragma omp critical
                {
                    if (mtar_read_header(&tar, &header) != MTAR_ENULLRECORD) {
                        //TODO GNU tar has special blocks for long filenames
                        name = header.name;
                        if (header.size > bufferSize) {
                            bufferSize = header.size * 1.5;
                            dataBuffer = (char*)realloc(dataBuffer, bufferSize);
                        }
                        if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                            std::cerr << "[Error] Reading tar entry " << name << " failed." << std::endl;
                            writeEntry = false;
                            proceed = false;
                        }
                        else {
                            writeEntry = true;
                            proceed = true;
                        }
                        mtar_next(&tar);
                        writeEntry = (header.type == MTAR_TREG) ? writeEntry : false;
                    }
                    else {
                        proceed = false;
                        writeEntry = false;
                    }
                } // end read in
                if (proceed && writeEntry) {
                    std::istringstream input(std::string(dataBuffer, header.size));
                    std::string name_clean = name.substr(name.find_last_of("/\\") + 1);
                    check(input, name_clean);
                }
            } // end while loop
        } // end openmp
        flag = 0;
    } else {
        std::cerr << "Invalid mode." << std::endl;
        return 1;
    }    // Print log
    if (mode != RMSD) {
        std::clog << "Done." << std::endl;
    }
    return flag;
}
