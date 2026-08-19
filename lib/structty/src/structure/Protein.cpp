#include "Protein.hpp"

static gemmi::Structure read_structure(const std::string& path) {
    gemmi::Structure st = gemmi::read_structure_file(path);
    st.remove_empty_chains();
    st.merge_chain_parts();
    return st;
}

// target_chains 는 "-"(전체 허용) 또는 콤마로 구분된 체인 ID 목록이다.
// 토큰 단위 정확 비교 — substring 비교였을 때는 필터 "C" 가 체인 "C-2" 에도 걸렸다.
static inline bool chain_ok(const std::string& target, const std::string& cid) {
    if (target == "-") return true;
    size_t start = 0;
    while (true) {
        const size_t comma = target.find(',', start);
        const size_t end = (comma == std::string::npos) ? target.size() : comma;
        if (target.compare(start, end - start, cid) == 0) return true;
        if (comma == std::string::npos) return false;
        start = comma + 1;
    }
}

Protein::Protein(const std::string& in_file_, const std::string& target_chains_, const bool& show_structure_) {
    structureMaker = StructureMaker();

    in_file = in_file_;
    target_chains = target_chains_;
    show_structure = show_structure_;
    cx = cy = cz = scale = 0.0;
}

Protein::Protein(const std::string& name, bool show_structure_val) {
    structureMaker = StructureMaker();

    in_file = name;
    target_chains = "-";
    show_structure = show_structure_val;
    cx = cy = cz = scale = 0.0;
}

Protein::~Protein() {
}

std::map<std::string, std::vector<Atom>>& Protein::get_atoms() {
    return screen_atoms;  
}

std::map<std::string, int> Protein::get_residue_count() {
    return chain_res_count;
}

std::map<std::string, int> Protein::get_chain_length() {
    std::map<std::string, int> result;
    for (const auto& [chainID, atoms] : init_atoms) {
        result[chainID] = atoms.size();
    }
    return result;
}

int Protein::get_chain_length(std::string chainID) {
    if (screen_atoms.count(chainID)) {
        return screen_atoms[chainID].size();
    }
    return 0;
}

int Protein::get_length() {
    int total_atoms = 0;
    for (const auto& [chainID, atoms] : init_atoms) {
        total_atoms += atoms.size();
    }
    return total_atoms;
}



BoundingBox& Protein::get_bounding_box() {
     return bounding_box; 
}

void Protein::set_scale(float scale_) { 
    scale = scale_;
    cx = 0.5f * (bounding_box.min_x + bounding_box.max_x);
    cy = 0.5f * (bounding_box.min_y + bounding_box.max_y);
    cz = 0.5f * (bounding_box.min_z + bounding_box.max_z);
    ssPredictor.set_scale(1.0f/scale);
}    

bool Protein::is_ss_in_file(const std::string& in_file) {
    try {
        gemmi::Structure st = read_structure(in_file);
        return (!st.helices.empty() || !st.sheets.empty());
    } catch (...) {
        return false;
    }
}

void Protein::set_bounding_box() {
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            bounding_box.min_x = std::min(bounding_box.min_x, atom.x);
            bounding_box.min_y = std::min(bounding_box.min_y, atom.y);
            bounding_box.min_z = std::min(bounding_box.min_z, atom.z);
            bounding_box.max_x = std::max(bounding_box.max_x, atom.x);
            bounding_box.max_y = std::max(bounding_box.max_y, atom.y);
            bounding_box.max_z = std::max(bounding_box.max_z, atom.z);
        }
    }
}         

void Protein::count_seqres(const std::string& file) {
    // std::cout << "  count SEQRES\n";
    chain_res_count.clear();

    try {
        gemmi::Structure st = read_structure(file);

        for (const gemmi::Entity& ent : st.entities) {
            int len = (int)ent.full_sequence.size();
            if (len <= 0) continue;

            for (const std::string& cname : ent.subchains) {
                if (cname.empty()) continue;
                std::string cid = cname;
                chain_res_count[cid] = len;
            }
        }
    } catch (...) {
        std::cerr << "Error: gemmi SEQRES read failed\n";
    }
}

void Protein::load_init_atoms(const std::string& in_file, 
                              const std::string& target_chains,
                              const std::vector<std::tuple<std::string, int, std::string, int, char>>& ss_info, 
                              float * vectorpointers , bool yesUT) {
    // std::cout << "  load atoms\n";
    init_atoms.clear();

    gemmi::Structure st = read_structure(in_file);
    gemmi::Model& model = st.first_model();

    for (gemmi::Chain& chain : model.chains) {
        std::string cid = chain.name.empty() ? "?" : chain.name;
        if (!chain_ok(target_chains, cid))
            continue;

        for (gemmi::Residue& res : chain.residues) {
            const gemmi::Atom* ca = res.get_ca();
            if (!ca) continue;
            if (!res.seqid.num.has_value()) continue;

            int resn = (int)res.seqid.num;
            float x = (float)ca->pos.x;
            float y = (float)ca->pos.y;
            float z = (float)ca->pos.z;

            Atom a(x, y, z);
            a.bfactor        = (float)ca->b_iso;
            a.residue_number = resn;
            a.residue_name   = res.name;

            for (auto& t : ss_info) {
                std::string sc; int s; std::string ec; int e; char type;
                std::tie(sc, s, ec, e, type) = t;

                if (cid == sc && resn >= s && resn <= e) {
                    a.set_structure(type);
                    break;
                }
            }
            init_atoms[cid].push_back(a);
        }
    }
}

void Protein::load_init_atoms(const std::string& in_file, 
                              const std::string& target_chains, float * vectorpointers, bool yesUT) {
    // std::cout << "  load atoms\n";
    init_atoms.clear();

    gemmi::Structure st = read_structure(in_file);
    gemmi::Model& model = st.first_model();

    for (gemmi::Chain& chain : model.chains) {
        std::string cid = chain.name.empty() ? "?" : chain.name;
        if (!chain_ok(target_chains, cid))
            continue;

        for (gemmi::Residue& res : chain.residues) {
            const gemmi::Atom* ca = res.get_ca();
            if (!ca) continue;

            float x = (float)ca->pos.x;
            float y = (float)ca->pos.y;
            float z = (float)ca->pos.z;

            Atom a(x, y, z, 'x');
            a.bfactor      = (float)ca->b_iso;
            a.residue_name = res.name;
            if (res.seqid.num.has_value())
                a.residue_number = (int)res.seqid.num;
            init_atoms[cid].push_back(a);
        }
    }
}

void Protein::load_ss_info(const std::string& in_file,
                               const std::string& target_chains,
                               std::vector<std::tuple<std::string,int,std::string,int,char>>& ss_info)
{
    // std::cout << "  load SS info\n";
    ss_info.clear();

    gemmi::Structure st = read_structure(in_file);

    // Helix → H
    for (const gemmi::Helix& h : st.helices) {
        auto beg = h.start;
        auto end = h.end;

        if (!beg.res_id.seqid.num.has_value() ||
            !end.res_id.seqid.num.has_value())
            continue;

        std::string bc = beg.chain_name.empty() ? std::string("?") : beg.chain_name;
        std::string ec = end.chain_name.empty() ? std::string("?") : end.chain_name;
        if (!chain_ok(target_chains, bc)) continue;

        int bs = (int)beg.res_id.seqid.num;
        int es = (int)end.res_id.seqid.num;

        ss_info.emplace_back(bc, bs, ec, es, 'H');
    }

    // Sheet → S
    for (const gemmi::Sheet& sheet : st.sheets) {
        for (const gemmi::Sheet::Strand& s : sheet.strands) {
            auto beg = s.start;
            auto end = s.end;

            if (!beg.res_id.seqid.num.has_value() ||
                !end.res_id.seqid.num.has_value())
                continue;

            std::string bc = beg.chain_name.empty() ? std::string("?"): beg.chain_name;
            std::string ec = end.chain_name.empty() ? std::string("?") : end.chain_name;
            if (!chain_ok(target_chains, bc)) continue;

            int bs = (int)beg.res_id.seqid.num;
            int es = (int)end.res_id.seqid.num;

            ss_info.emplace_back(bc, bs, ec, es, 'S');
        }
    }
}

std::ostream& operator<<(std::ostream& os, const std::tuple<std::string, int, std::string, int, char>& t) {
    os << "("
       << std::get<0>(t) << ", "
       << std::get<1>(t) << ", "
       << std::get<2>(t) << ", "
       << std::get<3>(t) << ", "
       << std::get<4>(t) << ")";
    return os;
}

void Protein::load_data(float * vectorpointers, bool yesUT) {    
    // pdb
    if (in_file.find(".pdb") != std::string::npos || in_file.find(".cif") != std::string::npos) {
        if (show_structure){
            if (is_ss_in_file(in_file)){
                std::vector<std::tuple<std::string, int, std::string, int, char>> ss_info;
                load_ss_info(in_file, target_chains, ss_info);
                load_init_atoms(in_file, target_chains, ss_info, vectorpointers, yesUT);
            }
            else{
                load_init_atoms(in_file, target_chains, vectorpointers, yesUT);
                ssPredictor.run(init_atoms);
            }
        }
        else{
            load_init_atoms(in_file, target_chains, vectorpointers, yesUT);
        }
        
        if (init_atoms.empty()) {
            std::cerr << "Error: input PDB file is empty." << std::endl;
            return;
        }
        
        if (show_structure){
            structureMaker.calculate_ss_points(init_atoms, screen_atoms);
        }
        else{ screen_atoms = init_atoms; }
        count_seqres(in_file);
    }

    // others
    else{
        std::cerr << "Error: input file format is not supported." << std::endl;
        return;
    }
    std::cout << std::endl;
}

static const char* aa_to_three(char one_letter) {
    static const char* table[26] = {
        "ALA","ASX","CYS","ASP","GLU","PHE","GLY","HIS","ILE",
        "XLE","LYS","LEU","MET","ASN","PYL","PRO","GLN","ARG",
        "SER","THR","SEC","VAL","TRP","XAA","TYR","GLX"
    };
    int idx = std::toupper(static_cast<unsigned char>(one_letter)) - 'A';
    if (idx < 0 || idx >= 26) return "UNK";
    return table[idx];
}

void Protein::load_from_ca(const std::vector<float>& coords, size_t n_residues,
                            const std::string& aa_seq) {
    init_atoms.clear();
    ca_only_ = true;

    for (size_t j = 0; j < n_residues; ++j) {
        float x = coords[j];
        float y = coords[j + n_residues];
        float z = coords[j + 2 * n_residues];

        Atom a(x, y, z, 'x');
        a.residue_number = static_cast<int>(j + 1);
        if (j < aa_seq.size()) {
            a.residue_name = aa_to_three(aa_seq[j]);
        } else {
            a.residue_name = "UNK";
        }
        init_atoms["A"].push_back(a);
    }

    if (init_atoms.empty()) {
        std::cerr << "Error: Foldseek DB entry is empty." << std::endl;
        return;
    }

    if (show_structure) {
        ssPredictor.run(init_atoms);
        structureMaker.calculate_ss_points(init_atoms, screen_atoms);
    } else {
        screen_atoms = init_atoms;
    }

    chain_res_count.clear();
    chain_res_count["A"] = static_cast<int>(n_residues);
}

void Protein::set_rotate(int x_rotate, int y_rotate, int z_rotate){
    const float PI = 3.14159265359;
    // const float UNIT = 12;
    const float UNIT = 48;

    if (x_rotate != 0) {
        float values[9] = {1, 0, 0,
                   0, cos(x_rotate * PI / UNIT), -sin(x_rotate * PI / UNIT), 
                   0, sin(x_rotate * PI / UNIT), cos(x_rotate * PI / UNIT)};

        float* rotate_mat = new float[9];
        for (int i = 0; i < 9; i++) {
            rotate_mat[i] = values[i];
        }
        do_rotation(rotate_mat);
    }
    else if (y_rotate != 0) {
        float values[9] = {cos(y_rotate * PI / UNIT), 0, sin(y_rotate * PI / UNIT),
                   0, 1, 0, 
                   -sin(y_rotate * PI / UNIT), 0, cos(y_rotate * PI / UNIT)};

        float* rotate_mat = new float[9];
        for (int i = 0; i < 9; i++) {
            rotate_mat[i] = values[i];
        }
        do_rotation(rotate_mat);
    }
    else if (z_rotate != 0) {
        float values[9] = {cos(z_rotate * PI / UNIT), -sin(z_rotate * PI / UNIT), 0,
                   sin(z_rotate * PI / UNIT), cos(z_rotate * PI / UNIT), 0, 
                   0, 0, 1};

        float* rotate_mat = new float[9];
        for (int i = 0; i < 9; i++) {
            rotate_mat[i] = values[i];
        }
        do_rotation(rotate_mat);
    }

}


void Protein::set_shift(float shift_x, float shift_y, float shift_z) { 
    float* shift_mat = new float[3];
    shift_mat[0] = shift_x;
    shift_mat[1] = shift_y;
    shift_mat[2] = shift_z;
    do_shift(shift_mat);
}

void Protein::do_naive_rotation(float * rotate_mat) {
    float avgx = 0;
    float avgy = 0;
    float avgz = 0;
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            avgx = atom.x;
            avgy = atom.y;
            avgz = atom.z;
            atom.x  = avgx * rotate_mat[0] + avgy * rotate_mat[1]+ avgz * rotate_mat[2];
            atom.y  = avgx * rotate_mat[3] + avgy * rotate_mat[4]+ avgz * rotate_mat[5];
            atom.z  = avgx * rotate_mat[6] + avgy * rotate_mat[7]+ avgz * rotate_mat[8];
        }
    }
}
void Protein::do_rotation(float * rotate_mat) {
    float avgx = 0;
    float avgy = 0;
    float avgz = 0;
    int num = 0;
    // for (auto& [chainID, chain_atoms] : screen_atoms) {
    //     for (Atom& atom : chain_atoms) {
    //         avgx += atom.x;
    //         avgy += atom.y;
    //         avgz += atom.z;
    //         num += 1;   
    //     }
    // }
    // avgx /= num;
    // avgy /= num;
    // avgz /= num;

    float minx,maxx,miny,maxy,minz,maxz;
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            if (num == 0) {
                minx = atom.x;
                maxx = atom.x;
                miny = atom.y;
                maxy = atom.y;
                minz = atom.z;
                maxz = atom.z;
            } else {
                if (atom.x <minx) {
                    minx = atom.x;
                }
                if (atom.x > maxx) {
                    maxx = atom.x;
                }
                if (atom.y <miny) {
                    miny = atom.y;
                }
                if (atom.y > maxy) {
                    maxy = atom.y;
                }
                if (atom.z <minz) {
                    minz = atom.z;
                }
                if (atom.z > maxz) {
                    maxz = atom.z;
                }
            }
            num++;
        }
    }
    avgx = (minx + maxx) /2;
    avgy = (miny + maxy) /2;
    avgz = (minz + maxz) /2;
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            float x = atom.x;
            float y = atom.y;
            float z = atom.z;
            atom.x = avgx + rotate_mat[0] * (x - avgx) + rotate_mat[1] * (y - avgy) + rotate_mat[2] * (z - avgz);
            atom.y = rotate_mat[3] * (x - avgx) + avgy + rotate_mat[4] * (y - avgy) + rotate_mat[5] * (z - avgz);
            atom.z = rotate_mat[6] * (x - avgx) + rotate_mat[7] * (y - avgy) +  avgz + rotate_mat[8] * (z - avgz);
        }
    }
}

void Protein::do_shift(float* shift_mat) {
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            atom.x += shift_mat[0];
            atom.y += shift_mat[1];
            atom.z += shift_mat[2];
        }
    }
}


void Protein::do_scale(float scale) {
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            atom.x *= scale;
            atom.y *= scale;
            atom.z *= scale;
        }
    }
}

// 기능 1: Interface Region ─────────────────────────────────────────────────
// 두 체인 간 CA-CA 거리 < threshold 인 잔기들을 init_atoms에서 is_interface = true 로 표시
void Protein::compute_interface_pair(const std::string& chain_A,
                                     const std::string& chain_B,
                                     float threshold) {
    float thr2 = threshold * threshold;

    auto& atomsA = init_atoms[chain_A];
    auto& atomsB = init_atoms[chain_B];

    for (Atom& a : atomsA) {
        for (const Atom& b : atomsB) {
            float dx = a.x - b.x;
            float dy = a.y - b.y;
            float dz = a.z - b.z;
            if (dx*dx + dy*dy + dz*dz < thr2) {
                a.is_interface = true;
                break;
            }
        }
    }

    for (Atom& b : atomsB) {
        for (const Atom& a : atomsA) {
            float dx = b.x - a.x;
            float dy = b.y - a.y;
            float dz = b.z - a.z;
            if (dx*dx + dy*dy + dz*dz < thr2) {
                b.is_interface = true;
                break;
            }
        }
    }
}

// init_atoms.is_interface를 screen_atoms으로 전파 (공통 헬퍼: bool 필드용)
// - residue_number >= 0 (coil/junction): residue_number로 직접 매핑
// - residue_number == -1 (H/S geometry): screen_atoms 내 인덱스 기준 가장 가까운 coil 원자 값 사용
//   (좌표계 불일치 방지: scaled screen_atoms와 unscaled init_atoms의 3D 비교는 부정확)
static void sync_bool_field_to_screen(
    std::map<std::string, std::vector<Atom>>& screen_atoms,
    const std::map<std::string, std::vector<Atom>>& init_atoms,
    bool Atom::* field)
{
    for (auto& [chainID, screen_chain] : screen_atoms) {
        auto init_it = init_atoms.find(chainID);
        if (init_it == init_atoms.end()) continue;
        const std::vector<Atom>& init_chain = init_it->second;

        int n = (int)screen_chain.size();

        // 1단계: residue_number >= 0 (coil/junction) 원자를 init_atoms에서 매핑
        //         동시에 coil 원자의 (index, value) 리스트 수집 (H/S geometry용)
        struct CoilEntry { int idx; bool val; };
        std::vector<CoilEntry> coil_entries;
        coil_entries.reserve(n / 2);

        for (int k = 0; k < n; ++k) {
            Atom& sa = screen_chain[k];
            if (sa.residue_number >= 0) {
                bool found_val = sa.*field;  // 기본값 유지
                for (const Atom& ia : init_chain) {
                    if (ia.residue_number == sa.residue_number) {
                        found_val = ia.*field;
                        break;
                    }
                }
                sa.*field = found_val;
                coil_entries.push_back({k, found_val});
            }
        }

        // 2단계: residue_number == -1 (H/S geometry) 원자 → 인덱스 기준 최근접 coil 원자
        if (!coil_entries.empty()) {
            for (int k = 0; k < n; ++k) {
                Atom& sa = screen_chain[k];
                if (sa.residue_number < 0) {
                    int best_dist = std::numeric_limits<int>::max();
                    bool best_val = false;
                    for (const CoilEntry& ce : coil_entries) {
                        int d = (k > ce.idx) ? (k - ce.idx) : (ce.idx - k);
                        if (d < best_dist) {
                            best_dist = d;
                            best_val  = ce.val;
                        }
                    }
                    sa.*field = best_val;
                }
            }
        }
    }
}

void Protein::sync_interface_to_screen() {
    sync_bool_field_to_screen(screen_atoms, init_atoms, &Atom::is_interface);
}

void Protein::sync_aligned_to_screen() {
    sync_bool_field_to_screen(screen_atoms, init_atoms, &Atom::is_aligned);
}

// 기능 4: UT transform을 init_atoms에도 적용 (screen_atoms와 동일한 rotation+shift)
void Protein::apply_ut_to_init_atoms(const float* U, const float* T) {
    for (auto& [chainID, chain_atoms] : init_atoms) {
        for (Atom& atom : chain_atoms) {
            float x = atom.x;
            float y = atom.y;
            float z = atom.z;
            atom.x = U[0]*x + U[1]*y + U[2]*z + T[0];
            atom.y = U[3]*x + U[4]*y + U[5]*z + T[1];
            atom.z = U[6]*x + U[7]*y + U[8]*z + T[2];
        }
    }
}

// 기능 4: nearest-neighbor 기반 정렬 잔기 표시 (UT transform 이미 적용된 init_atoms 사용)
// this와 other 두 단백질 모두 is_aligned 설정 후 screen_atoms에 반영
void Protein::compute_aligned_regions_nn(Protein& other, float threshold) {
    float thr2 = threshold * threshold;

    // query (this) 측: 최근접 target CA < threshold이면 is_aligned = true
    for (auto& [cid_a, chain_a] : init_atoms) {
        for (Atom& a : chain_a) {
            for (const auto& [cid_b, chain_b] : other.init_atoms) {
                for (const Atom& b : chain_b) {
                    float dx = a.x - b.x;
                    float dy = a.y - b.y;
                    float dz = a.z - b.z;
                    if (dx*dx + dy*dy + dz*dz < thr2) {
                        a.is_aligned = true;
                        break;
                    }
                }
                if (a.is_aligned) break;
            }
        }
    }

    // target (other) 측: 최근접 query CA < threshold이면 is_aligned = true
    for (auto& [cid_b, chain_b] : other.init_atoms) {
        for (Atom& b : chain_b) {
            for (const auto& [cid_a, chain_a] : init_atoms) {
                for (const Atom& a : chain_a) {
                    float dx = b.x - a.x;
                    float dy = b.y - a.y;
                    float dz = b.z - a.z;
                    if (dx*dx + dy*dy + dz*dz < thr2) {
                        b.is_aligned = true;
                        break;
                    }
                }
                if (b.is_aligned) break;
            }
        }
    }

    sync_aligned_to_screen();
    other.sync_aligned_to_screen();
}

// 기능 4: alignment string 기반 정렬 잔기 표시 (Foldseek qaln/taln 사용)
// qaln[i] != '-' and taln[i] != '-' 인 위치의 CA 쌍 거리 < threshold이면 양쪽 is_aligned = true
void Protein::compute_aligned_regions_from_aln(Protein& other,
                                               const std::string& qaln,
                                               const std::string& taln,
                                               int q_start,
                                               int t_start,
                                               float threshold,
                                               bool skip_distance_check) {
    if (qaln.size() != taln.size()) return;
    float thr2 = threshold * threshold;

    // is_aligned 플래그 초기화 (multihit 전환 시 이전 hit 잔류 방지)
    for (auto& [cid, chain] : init_atoms) {
        for (Atom& a : chain) {
            a.is_aligned = false;
        }
    }
    for (auto& [cid, chain] : other.init_atoms) {
        for (Atom& a : chain) {
            a.is_aligned = false;
        }
    }

    // query CA flat 목록 (chain 이름 순 = map 순)
    std::vector<Atom*> query_cas;
    for (auto& [cid, chain] : init_atoms) {
        for (Atom& a : chain) {
            query_cas.push_back(&a);
        }
    }

    // target CA flat 목록
    std::vector<Atom*> target_cas;
    for (auto& [cid, chain] : other.init_atoms) {
        for (Atom& a : chain) {
            target_cas.push_back(&a);
        }
    }

    // qaln/taln 은 q_start/t_start 잔기(1-based)에서 시작한다. 0 이하 값은 1로 취급.
    int q_idx = (q_start > 1) ? (q_start - 1) : 0;
    int t_idx = (t_start > 1) ? (t_start - 1) : 0;
    const int q_size = (int)query_cas.size();
    const int t_size = (int)target_cas.size();

    for (size_t i = 0; i < qaln.size(); ++i) {
        bool q_gap = (qaln[i] == '-');
        bool t_gap = (taln[i] == '-');

        if (!q_gap && !t_gap) {
            if (q_idx < q_size && t_idx < t_size) {
                Atom* qa = query_cas[q_idx];
                Atom* ta = target_cas[t_idx];
                if (skip_distance_check) {
                    // alignment string만으로 is_aligned 설정 (좌표 비교 생략)
                    qa->is_aligned = true;
                    ta->is_aligned = true;
                } else {
                    float dx = qa->x - ta->x;
                    float dy = qa->y - ta->y;
                    float dz = qa->z - ta->z;
                    if (dx*dx + dy*dy + dz*dz < thr2) {
                        qa->is_aligned = true;
                        ta->is_aligned = true;
                    }
                }
            }
        }

        if (!q_gap) q_idx++;
        if (!t_gap) t_idx++;
    }

    sync_aligned_to_screen();
    other.sync_aligned_to_screen();
}

// 기능 5: float 필드를 init_atoms → screen_atoms로 전파하는 공통 헬퍼
// - residue_number >= 0 (coil/junction): residue_number로 직접 매핑
// - residue_number == -1 (H/S geometry): 인덱스 기준 최근접 coil 원자 값 사용
static void sync_float_field_to_screen(
    std::map<std::string, std::vector<Atom>>& screen_atoms,
    const std::map<std::string, std::vector<Atom>>& init_atoms,
    float Atom::* field,
    float default_val = -1.0f)
{
    for (auto& [chainID, screen_chain] : screen_atoms) {
        auto init_it = init_atoms.find(chainID);
        if (init_it == init_atoms.end()) continue;
        const std::vector<Atom>& init_chain = init_it->second;

        int n = (int)screen_chain.size();

        struct CoilEntry { int idx; float val; };
        std::vector<CoilEntry> coil_entries;
        coil_entries.reserve(n / 2 + 1);

        // 1단계: coil/junction (residue_number >= 0) → init_atoms에서 직접 매핑
        for (int k = 0; k < n; ++k) {
            Atom& sa = screen_chain[k];
            if (sa.residue_number >= 0) {
                float found_val = default_val;
                for (const Atom& ia : init_chain) {
                    if (ia.residue_number == sa.residue_number) {
                        found_val = ia.*field;
                        break;
                    }
                }
                sa.*field = found_val;
                coil_entries.push_back({k, found_val});
            }
        }

        // 2단계: H/S geometry (residue_number == -1) → 인덱스 기준 최근접 coil 원자 값
        if (!coil_entries.empty()) {
            for (int k = 0; k < n; ++k) {
                Atom& sa = screen_chain[k];
                if (sa.residue_number < 0) {
                    int best_dist = std::numeric_limits<int>::max();
                    float best_val = default_val;
                    for (const CoilEntry& ce : coil_entries) {
                        int d = (k > ce.idx) ? (k - ce.idx) : (ce.idx - k);
                        if (d < best_dist) {
                            best_dist = d;
                            best_val  = ce.val;
                        }
                    }
                    sa.*field = best_val;
                }
            }
        }
    }
}

void Protein::sync_conservation_to_screen() {
    sync_float_field_to_screen(screen_atoms, init_atoms, &Atom::conservation_score, -1.0f);
}

// 기능 5: 0-based 인덱스 순서로 conservation_score를 init_atoms에 적용 후 screen_atoms에 전파
// scores 벡터 길이는 init_atoms 전체 잔기 수와 일치해야 정확하다.
// map 순서(체인 ID 사전순) → 각 체인 내 순서(읽은 순서)로 순회하여 적용한다.
void Protein::apply_conservation_scores(const std::vector<float>& scores) {
    int score_len = (int)scores.size();
    if (score_len == 0) return;

    // 전체 잔기 수 계산
    int total_residues = 0;
    for (auto& [cid, chain] : init_atoms) {
        total_residues += (int)chain.size();
    }

    if (total_residues <= score_len) {
        // 기존 방식: scores가 전체 잔기를 커버하면 순차 적용
        int idx = 0;
        for (auto& [cid, chain] : init_atoms) {
            for (Atom& a : chain) {
                if (idx < score_len) {
                    a.conservation_score = scores[idx];
                }
                idx++;
            }
        }
    } else {
        // Homodimer 대응: scores가 전체 잔기보다 짧으면 체인별로 반복 적용
        for (auto& [cid, chain] : init_atoms) {
            int idx = 0;
            for (Atom& a : chain) {
                if (idx < score_len) {
                    a.conservation_score = scores[idx];
                }
                idx++;
            }
        }
    }

    sync_conservation_to_screen();
}

// 공개 진입점: 모든 체인 쌍에 대해 interface를 계산하고 screen_atoms에 반영
void Protein::compute_interface(float threshold) {
    // 체인 ID 목록 수집
    std::vector<std::string> chain_ids;
    chain_ids.reserve(init_atoms.size());
    for (const auto& [cid, _] : init_atoms) {
        chain_ids.push_back(cid);
    }

    // 체인이 1개 이하면 inter-chain interface 없음
    if (chain_ids.size() < 2) return;

    // 모든 체인 쌍 계산
    for (size_t i = 0; i < chain_ids.size(); ++i) {
        for (size_t j = i + 1; j < chain_ids.size(); ++j) {
            compute_interface_pair(chain_ids[i], chain_ids[j], threshold);
        }
    }

    // screen_atoms으로 전파
    sync_interface_to_screen();
}