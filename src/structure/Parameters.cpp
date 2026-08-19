#include "Parameters.hpp"

#include "InputProbe.hpp"

void print_help(){
    std::cout << "Usage: StrucTTY <query...> [OPTIONS]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -m, --mode <MODE>       Color mode:\n";
    std::cout << "                            protein (default), chain, rainbow,\n";
    std::cout << "                            plddt, interface, conservation, align,\n";
    std::cout << "                            align-fs (Foldseek alignment only),\n";
    std::cout << "                            align-near (nearest-neighbour distance)\n";
    std::cout << "  -c, --chains <FILE>     Show only selected chains (see example/chainfile)\n";
    std::cout << "  -s, --structure         Show secondary structure (alpha helix, beta sheet)\n";
    std::cout << "  --msa <FILE>            MSA file for conservation score (FASTA/A3M)\n";
    std::cout << "  -fst, --foldseek-target <PATH>\n";
    std::cout << "                          Target source for Foldseek hits: Foldseek DB,\n";
    std::cout << "                          structure directory, structure file, or 'auto'\n";
    std::cout << "                          ('auto' downloads hits from public DBs)\n";
    std::cout << "  -fsr, --foldseek-result <FILE>\n";
    std::cout << "                          Foldseek result: m8 (12/17/21/29 columns) or\n";
    std::cout << "                          multimer _report (14 columns)\n";
    std::cout << "                          -fst and -fsr must be given together\n";
    std::cout << "  -fm, --foldmason <FILE> FoldMason result (JSON or FASTA MSA)\n";
    std::cout << "  -n, --nopanel           Hide info panel\n";
    std::cout << "  -b, --benchmark         Benchmark mode (measure FPS/latency)\n";
    std::cout << "  --help                  Show this help message\n";
    std::cout << "\nSupported inputs (4 kinds; detected automatically):\n";
    std::cout << "  1. structure file       .pdb / .cif / .ent (+ .gz)\n";
    std::cout << "  2. structure directory  a directory of those files. As the query, every\n";
    std::cout << "                          accession in -fsr is looked up inside it (]/[ walks them)\n";
    std::cout << "  3. Foldseek DB          base path of a DB built from structures\n";
    std::cout << "                          (needs <db>_ca; sequence-derived DBs have none)\n";
    std::cout << "  4. Foldseek result      m8 (12/17/21/29 columns) or\n";
    std::cout << "                          multimer _report (14 columns)\n";
    std::cout << "  Sequence FASTA is NOT supported -- it carries no 3D coordinates.\n";
    std::cout << "\nRecipes:\n";
    std::cout << "  StrucTTY query.cif\n";
    std::cout << "  StrucTTY query.cif target.cif -m aligned\n";
    std::cout << "  StrucTTY query.cif -fst targetDB   -fsr result.m8   -m aligned\n";
    std::cout << "  StrucTTY query.cif -fst pdb_dir/   -fsr result.m8\n";
    std::cout << "  StrucTTY query.cif -fst auto       -fsr result.m8\n";
    std::cout << "  StrucTTY query_dir/ -fst targetDB  -fsr result.m8   (multi-query: ]/[)\n";
    std::cout << "  StrucTTY queryDB   -fst targetDB   -fsr result.m8   (same, from a DB)\n";
    std::cout << "  StrucTTY queryDB   -fst targetDB   -fsr out_report  (multimer)\n";
}

Parameters::Parameters(int argc, char* argv[]) {
    arg_okay = true;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_help();
            std::exit(0);
        }
    }
    
    if (argc <= 1) {
        std::cerr << "Need input file dir" << std::endl;
        arg_okay = false;
        return;
    }

    for (int i = 1; i < argc; i++) {
        try {
            if (!strcmp(argv[i], "-m") || !strcmp(argv[i], "--mode")) {
                if (i + 1 < argc) {
                    std::string val(argv[i + 1]);
                    std::transform(val.begin(), val.end(), val.begin(), ::tolower); // to lowercase
                    if (val == "chain" || val == "rainbow" || val == "protein" ||
                        val == "plddt" || val == "interface" || val == "conservation" ||
                        val == "align" || val == "align-fs" || val == "align-near") {
                        mode = val;
                        i++;
                    } else {
                        throw std::runtime_error("Error: Invalid value for --mode. Use 'protein', 'chain', 'rainbow', 'plddt', 'interface', 'conservation', 'align', 'align-fs', or 'align-near'.");
                    }
                } else {
                    throw std::runtime_error("Error: Missing value for -m / --mode.");
                }
            } else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--chains")) {
                if (i + 1 < argc) {
                    chainfile = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for -c / --chains.");
                }
            }
            else if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "--structure")) {
                show_structure = true;
            } else if (!strcmp(argv[i], "-n") || !strcmp(argv[i], "--nopanel")) {
                no_panel = true;
            } else if (fs::exists(argv[i]) &&
                       (fs::is_regular_file(argv[i]) || fs::is_directory(argv[i])) &&
                       in_file.size() < 9){
                // 디렉터리 query: -fsr 의 query accession 들을 이 디렉터리에서 찾아
                // 체인마다 하나씩 순회한다(foldseek 뷰어가 query DB 로 하는 것과 같은 구성).
                in_file.push_back(argv[i]);
            } else if (!strcmp(argv[i], "--msa")) {
                if (i + 1 < argc) {
                    msa_file = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for --msa.");
                }
            } else if (!strcmp(argv[i], "-fst") || !strcmp(argv[i], "--foldseek-target")) {
                if (i + 1 < argc) {
                    foldseek_target = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for -fst / --foldseek-target.");
                }
            } else if (!strcmp(argv[i], "-fsr") || !strcmp(argv[i], "--foldseek-result")) {
                if (i + 1 < argc) {
                    foldseek_result = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for -fsr / --foldseek-result.");
                }
            } else if (!strcmp(argv[i], "--foldmason") || !strcmp(argv[i], "-fm")) {
                if (i + 1 < argc) {
                    foldmason_file = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for --foldmason / -fm.");
                }
            } else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--benchmark")) {
                benchmark_mode = true;
                show_structure = true;
            } else {
                throw std::runtime_error("Error: Unknown parameter: " + std::string(argv[i]));
            }
        }       
        catch (const std::exception& e) {
            std::cerr << "Wrong input parameters: " << e.what() << std::endl;
            std::cerr << "Error at argument: " << argv[i] << std::endl;
            arg_okay = false;
            return;
        }
    }
    if (!input_probe::validate_inputs(in_file, foldseek_target, foldseek_result)) {
        arg_okay = false;
        return;
    }
    return;
}

void Parameters::print_args() {
    cout << "Input parameters >> " << endl;
    cout << "  in_file: " << endl;
    for (size_t i = 0; i < in_file.size(); i++) {
        std::cout << "\t" << in_file[i] << '\n';
    }
    cout << "  mode: " << mode << endl;
    cout << "  chainfile: " << chainfile << endl;
    cout << "  show_structure: " << show_structure << endl;
    cout << "  benchmark_mode: " << benchmark_mode << endl;
    if (!foldseek_target.empty() || !foldseek_result.empty()) {
        cout << "  foldseek_target: " << foldseek_target << endl;
        cout << "  foldseek_result: " << foldseek_result << endl;
    }

    cout << "\n";
    return;
}
