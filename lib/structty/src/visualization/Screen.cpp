#include "Screen.hpp"
#include "Terminal.hpp"
#include "Common.hpp"
#include <cstring>       // strncpy, memset
#include <unordered_set>
#include <string>
#include <limits>
#include <filesystem>

const float FOV = 90.0f;
const float PI  = 3.14159265359f;
const int MAX_STRUCT_NUM = 9;

Screen::Screen(const int& width, const int& height, const bool& show_structure,
               const std::string& mode)
    : renderer_(width, height, mode, show_structure)
{
    screen_width = width;
    screen_height = height;
    screen_show_structure = show_structure;
    screen_mode = mode;
    aspect_ratio = (float)screen_width / screen_height;
    zoom_level = 2.0f;

    camera = new Camera(width, height, mode);
    panel  = new Panel(width, mode, show_structure);
}

Screen::~Screen() {
    for (Protein* p : data) delete p;
    data.clear();
    free_tmatrix();
    delete camera;
    delete panel;
}

void Screen::free_tmatrix() {
    if (!vectorpointer) return;
    for (size_t i = 0; i < vectorpointer_len_; i++) delete[] vectorpointer[i];
    delete[] vectorpointer;
    vectorpointer = nullptr;
    vectorpointer_len_ = 0;
}

static float compute_scene_radius_from_render_positions(const std::vector<Protein*>& data) {
    float max_r2 = 0.0f;
    bool any = false;

    for (auto* p : data) {
        for (const auto& [chainID, chain_atoms] : p->get_atoms()) {
            for (const auto& atom : chain_atoms) {
                float* pos = atom.get_position();
                float x = pos[0], y = pos[1], z = pos[2];
                float r2 = x*x + y*y + z*z;
                if (r2 > max_r2) max_r2 = r2;
                any = true;
            }
        }
    }
    if (!any) return 1.0f;
    return std::sqrt(max_r2);
}


std::string Screen::chain_from_accession(const std::string& accession,
                                         const std::string& file_path) {
    if (accession.empty() || file_path.empty()) return "-";
    // foldseek createdb 는 체인마다 `<파일 stem>_<auth chain>` 엔트리를 만든다
    // (예: 1nzy-assembly1_A, 1nzy-assembly1_C-2).
    const std::string stem = std::filesystem::path(file_path).stem().string();
    if (stem.empty()) return "-";
    const std::string prefix = stem + "_";
    if (accession.size() <= prefix.size()) return "-";
    if (accession.compare(0, prefix.size(), prefix) != 0) return "-";
    const std::string chain = accession.substr(prefix.size());
    return chain.empty() ? "-" : chain;
}

std::string Screen::apply_accession_chain(int idx, const std::string& accession,
                                          const std::string& file_path) {
    if (idx < 0 || idx >= (int)chainVec.size()) return "";
    if (chainVec[idx] != "-") return "";  // 사용자 --chains 지정이 우선
    const std::string chain = chain_from_accession(accession, file_path);
    if (chain == "-") return "";
    chainVec[idx] = chain;
    return chain;
}

void Screen::set_protein(const std::string& in_file, int ii, const bool& show_structure) {
    Protein* protein = new Protein(in_file, chainVec.at(ii), show_structure);
    data.push_back(protein);
    pan_x.push_back(0.0f);
    pan_y.push_back(0.0f);

}

bool Screen::load_query_into_data0(const std::string& accession,
                                   const bool& show_structure) {
    std::vector<float> coords;
    std::string aa_seq;
    size_t n_res = query_db_reader_.read_entry(accession, coords, aa_seq);
    if (n_res == 0) {
        std::cerr << "Warning: query accession not found in query DB: " << accession << "\n";
        return false;
    }

    Protein* query_protein = new Protein(accession, show_structure);
    query_protein->load_from_ca(coords, n_res, aa_seq);
    data.push_back(query_protein);
    pan_x.push_back(0.0f);
    pan_y.push_back(0.0f);
    // chainVec was sized by set_chainfile(); ensure an entry exists for data[0].
    if ((int)chainVec.size() < (int)data.size()) chainVec.push_back("-");
    return true;
}

int Screen::load_complex_chains(const std::string& complex,
                                const std::vector<std::string>& chains,
                                bool from_db, FoldseekDBReader& reader,
                                const std::string& source_path, bool source_is_dir,
                                const bool& show_structure) {
    int loaded = 0;
    if (from_db) {
        for (const std::string& ch : chains) {
            if (load_chain_into_data(reader, complex + "_" + ch, show_structure)) {
                mm_chain_labels_.push_back(complex + "_" + ch);
                loaded++;
            }
        }
        return loaded;
    }

    std::string file_path = source_path;
    if (source_is_dir || source_path.empty()) {
        std::string status_msg;
        file_path = PDBDownloader::resolve_target_file(complex, source_path, status_msg);
    }
    if (file_path.empty()) {
        std::cerr << "Warning: no structure file for complex " << complex << "\n";
        return 0;
    }

    for (const std::string& ch : chains) {
        Protein* p = new Protein(file_path, ch, show_structure);
        float tvec[3] = {0.f, 0.f, 0.f};
        p->load_data(tvec, false);
        if (p->get_chain_length(ch) == 0) {
            std::cerr << "Warning: chain " << ch << " not found in " << file_path << "\n";
            delete p;
            continue;
        }
        data.push_back(p);
        pan_x.push_back(0.0f);
        pan_y.push_back(0.0f);
        chainVec.push_back(ch);
        mm_chain_labels_.push_back(complex + "_" + ch);
        loaded++;
    }
    return loaded;
}

bool Screen::set_query_from_db(const std::string& query_db_path,
                               const std::string& accession,
                               const bool& show_structure) {
    if (!query_db_reader_.open(query_db_path)) {
        std::cerr << "Warning: Failed to open query Foldseek DB: " << query_db_path << "\n";
        return false;
    }
    std::unordered_set<std::string> q_ids = { accession };
    if (!query_db_reader_.prepare(q_ids)) {
        std::cerr << "Warning: Failed to prepare query DB for accession: " << accession << "\n";
        return false;
    }
    return load_query_into_data0(accession, show_structure);
}

void Screen::set_query_nav(const std::vector<std::string>& query_ids,
                           const std::map<std::string, std::vector<FoldseekHit>>& hits_by_query,
                           const std::string& query_db_path,
                           const std::string& target_db_path,
                           const bool& show_structure) {
    query_ids_ = query_ids;
    hits_by_query_ = hits_by_query;
    query_db_path_ = query_db_path;
    target_db_path_ = target_db_path;
    multi_query_show_structure_ = show_structure;

    // Open + prepare the query DB once for all query accessions.
    if (!query_db_reader_.is_open()) {
        if (!query_db_reader_.open(query_db_path)) {
            std::cerr << "Warning: Failed to open query Foldseek DB: " << query_db_path << "\n";
            return;
        }
    }
    std::unordered_set<std::string> q_set(query_ids.begin(), query_ids.end());
    query_db_reader_.prepare(q_set);

    // Open the target DB once (targets are read per-hit in load_next_hit).
    if (!target_db_path.empty() && !fs_db_reader_.is_open()) {
        open_foldseek_db(target_db_path);
    }

    current_query_idx_ = -1;
    activate_query(0);
}

void Screen::set_query_nav_from_file(const std::vector<std::string>& query_ids,
                                    const std::map<std::string, std::vector<FoldseekHit>>& hits_by_query,
                                    const std::string& query_source,
                                    const std::string& target_db_path,
                                    const bool& show_structure,
                                    const bool& source_is_directory) {
    query_ids_ = query_ids;
    hits_by_query_ = hits_by_query;
    query_file_ = query_source;
    query_file_is_dir_ = source_is_directory;
    target_db_path_ = target_db_path;
    multi_query_show_structure_ = show_structure;

    // Targets still come from a Foldseek DB when one was given; the query side is
    // the only difference from set_query_nav().
    if (!target_db_path.empty() && !fs_db_reader_.is_open()) {
        open_foldseek_db(target_db_path);
    }

    current_query_idx_ = -1;
    activate_query(0);
}

void Screen::activate_query(int idx) {
    if (idx < 0 || idx >= (int)query_ids_.size()) return;
    current_query_idx_ = idx;
    const std::string& acc = query_ids_[idx];

    // Tear down the current scene (query + any loaded targets).
    for (Protein* p : data) delete p;
    data.clear();
    pan_x.clear();
    pan_y.clear();
    chainVec.clear();
    chainVec.push_back("-");
    if (panel) panel->reset_entries();

    // Load the new query structure as data[0]. With a structure file as the query the
    // same file is re-read with the chain the accession names; otherwise it comes from
    // the (already prepared) query DB.
    if (!query_file_.empty()) {
        std::string query_path = query_file_;
        if (query_file_is_dir_) {
            // accession(`<stem>_<chain>`) 으로 디렉터리에서 원본 파일을 찾는다 —
            // target 쪽과 같은 규칙(뒤쪽 `_<chain>` 을 떼며 탐색).
            std::string ignored;
            query_path = PDBDownloader::resolve_target_file(acc, query_file_, ignored);
            if (query_path.empty()) {
                std::cerr << "Warning: no structure file for query '" << acc
                          << "' in " << query_file_ << std::endl;
                return;
            }
        }
        const std::string chain = chain_from_accession(acc, query_path);
        chainVec[0] = chain;   // "-" 면 전체 체인
        set_protein(query_path, 0, multi_query_show_structure_);
    } else if (!load_query_into_data0(acc, multi_query_show_structure_)) {
        return;
    }

    // Re-normalize for this query: set_tmatrix re-sizes vectorpointer, and
    // normalize_proteins() recomputes norm_scale / centroid from data[0] and
    // re-adds the query panel entry (entries were cleared above).
    set_tmatrix();
    normalize_proteins();
    update_total_len_ca();

    // Set this query's hit list and load its first target (which applies the
    // per-hit superposition using the freshly computed norm_scale / centroid).
    const auto it = hits_by_query_.find(acc);
    if (it != hits_by_query_.end() && !it->second.empty()) {
        set_foldseek_hits(it->second);
        prepare_foldseek_db(it->second);
        current_hit_idx = -1;
        load_next_hit(+1);
    } else {
        foldseek_hits.clear();
        current_hit_idx = -1;
        if (panel) {
            FoldseekHitInfo fi;
            fi.valid = true;
            fi.query_idx = current_query_idx_ + 1;
            fi.total_queries = (int)query_ids_.size();
            panel->set_foldseek_hit_info(fi);
        }
    }
}

void Screen::switch_query(int delta) {
    if ((int)query_ids_.size() < 2) return;
    int new_idx = current_query_idx_ + delta;
    new_idx = std::max(0, std::min((int)query_ids_.size() - 1, new_idx));
    if (new_idx == current_query_idx_) return;
    activate_query(new_idx);
}

// ── Step 7 (D6): multimer `_report` 경로 ──────────────────────────────────────

static std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= s.size()) {
        const size_t comma = s.find(',', start);
        const size_t end = (comma == std::string::npos) ? s.size() : comma;
        out.push_back(s.substr(start, end - start));
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    if (s.empty()) out.clear();
    return out;
}

bool Screen::load_chain_into_data(FoldseekDBReader& reader,
                                  const std::string& accession,
                                  const bool& show_structure) {
    std::vector<float> coords;
    std::string aa_seq;
    size_t n_res = reader.read_entry(accession, coords, aa_seq);
    if (n_res == 0) {
        std::cerr << "Warning: complex chain not found in DB: " << accession << "\n";
        return false;
    }
    Protein* p = new Protein(accession, show_structure);
    p->load_from_ca(coords, n_res, aa_seq);
    data.push_back(p);
    pan_x.push_back(0.0f);
    pan_y.push_back(0.0f);
    if ((int)chainVec.size() < (int)data.size()) chainVec.push_back("-");
    return true;
}

void Screen::normalize_complex() {
    // 모든 query 체인을 하나의 공유 centroid/scale 로 정규화 (상대 위치 보존, grid 없음).
    global_bb = BoundingBox();
    for (auto* p : data) {
        p->set_bounding_box();
        global_bb = global_bb + p->get_bounding_box();
    }
    float max_ext = std::max(global_bb.max_x - global_bb.min_x,
                             global_bb.max_y - global_bb.min_y);
    max_ext = std::max(max_ext, global_bb.max_z - global_bb.min_z);
    float scale = (max_ext > 0.f) ? (2.0f / max_ext) : 1.0f;

    float gx = 0.5f * (global_bb.min_x + global_bb.max_x);
    float gy = 0.5f * (global_bb.min_y + global_bb.max_y);
    float gz = 0.5f * (global_bb.min_z + global_bb.max_z);
    float shift[3] = { -gx, -gy, -gz };
    for (auto* p : data) {
        p->set_scale(scale);
        p->do_shift(shift);
        p->do_scale(scale);
    }

    norm_scale = scale;
    norm_cx = gx; norm_cy = gy; norm_cz = gz;
    rot_pivot_[0] = rot_pivot_[1] = rot_pivot_[2] = 0.f;
    for (size_t i = 0; i < pan_x.size(); i++) { pan_x[i] = 0.0f; pan_y[i] = 0.0f; }
    yesUT = true;             // overlay 프레임 (grid layout 비활성)
    depth_calibrated = false;

    float radius = compute_scene_radius_from_render_positions(data);
    focal_offset = std::clamp(2.5f * radius + 1.0f, 2.0f, 8.0f);
}

void Screen::transform_target_chain(int idx, const float U[9], const float T[3]) {
    if (idx < 0 || idx >= (int)data.size() || !data[idx]) return;
    Protein* tp = data[idx];
    tp->set_bounding_box();
    float t_cx = tp->cx, t_cy = tp->cy, t_cz = tp->cz;
    // query 와 동일한 norm_scale 로 우선 중심+스케일 (per-chain centroid 는 아래 공식에서 상쇄)
    tp->set_scale(norm_scale);
    float t_shift[3] = { -t_cx, -t_cy, -t_cz };
    tp->do_shift(t_shift);
    tp->do_scale(norm_scale);

    // complex U/T(Å, target→query) → norm 공간 T. 공식은 load_next_hit 29컬럼 경로와 동일:
    // result = (U*orig + T - norm_c) * norm_scale (per-chain t_cen 상쇄).
    float Utc[3] = {
        U[0]*t_cx + U[1]*t_cy + U[2]*t_cz,
        U[3]*t_cx + U[4]*t_cy + U[5]*t_cz,
        U[6]*t_cx + U[7]*t_cy + U[8]*t_cz
    };
    float Tn[3] = {
        (Utc[0] + T[0] - norm_cx) * norm_scale,
        (Utc[1] + T[1] - norm_cy) * norm_scale,
        (Utc[2] + T[2] - norm_cz) * norm_scale
    };
    float Ta[3] = { T[0], T[1], T[2] };
    apply_foldseek_transform(idx, U, Tn, Ta);
}

std::pair<int,int> Screen::mm_complex_range(int complex_idx) const {
    const int dsz = (int)data.size();
    if (complex_idx == 0) {
        int hi = std::min(mm_query_chain_count_ - 1, dsz - 1);
        return {0, hi};
    } else {
        if (mm_query_chain_count_ >= dsz) return {dsz, dsz - 1}; // empty range: no target loaded
        return {mm_query_chain_count_, dsz - 1};
    }
}

void Screen::set_multimer_report(const std::vector<MultimerHit>& hits,
                                 const std::string& query_source,
                                 bool query_is_db,
                                 bool query_is_dir,
                                 const std::string& target_db_path,
                                 const bool& show_structure) {
    multimer_mode_ = true;
    multi_query_show_structure_ = show_structure;
    mm_query_source_ = query_source;
    mm_query_is_db_ = query_is_db;
    mm_query_is_dir_ = query_is_dir;
    query_db_path_ = query_is_db ? query_source : std::string();
    target_db_path_ = target_db_path;

    // query complex 별로 그룹화 (등장 순서 보존)
    for (const auto& h : hits) {
        auto it = mm_hits_by_query_.find(h.qComplex);
        if (it == mm_hits_by_query_.end()) {
            mm_query_complexes_.push_back(h.qComplex);
            mm_hits_by_query_.emplace(h.qComplex, std::vector<MultimerHit>{ h });
        } else {
            it->second.push_back(h);
        }
    }
    if (mm_query_complexes_.empty()) { multimer_mode_ = false; return; }

    // query 가 DB 면 열고 모든 query 체인 accession 을 인덱싱한다.
    if (mm_query_is_db_) {
        if (!query_db_reader_.open(query_source)) {
            std::cerr << "Warning: Failed to open query complex DB: " << query_source << "\n";
            multimer_mode_ = false;
            return;
        }
        std::unordered_set<std::string> q_acc;
        for (const auto& qc : mm_query_complexes_)
            for (const auto& ch : mm_hits_by_query_[qc].front().qChains)
                q_acc.insert(qc + "_" + ch);
        query_db_reader_.prepare(q_acc);
    } else if (!mm_query_is_dir_) {
        // 구조 파일 하나로는 그 파일이 담은 complex 만 볼 수 있다.
        const std::string stem = std::filesystem::path(mm_query_source_).stem().string();
        std::vector<std::string> kept;
        for (const std::string& qc : mm_query_complexes_) {
            if (qc == stem) kept.push_back(qc);
        }
        if (kept.empty()) {
            std::cerr << "Error: none of the query complexes in the report match "
                      << mm_query_source_ << "\n"
                      << "       Pass the directory holding them, or a Foldseek DB." << std::endl;
            multimer_mode_ = false;
            return;
        }
        if (kept.size() != mm_query_complexes_.size()) {
            std::cerr << "Query complexes in " << mm_query_source_ << ": " << kept.size()
                      << " of " << mm_query_complexes_.size()
                      << " (pass a directory or a Foldseek DB to walk them all)" << std::endl;
        }
        mm_query_complexes_ = kept;
    }

    // target complex DB 열고 모든 target 체인 accession 인덱싱
    if (!target_db_path.empty() && fs_db_reader_.open(target_db_path)) {
        std::unordered_set<std::string> t_acc;
        for (const auto& h : hits)
            for (const auto& ch : h.tChains)
                t_acc.insert(h.tComplex + "_" + ch);
        fs_db_reader_.prepare(t_acc);
    }

    mm_current_query_idx_ = -1;
    activate_multimer_query(0);
}

void Screen::activate_multimer_query(int idx) {
    if (idx < 0 || idx >= (int)mm_query_complexes_.size()) return;
    mm_current_query_idx_ = idx;
    const std::string& qc = mm_query_complexes_[idx];

    // scene teardown
    for (Protein* p : data) delete p;
    data.clear();
    pan_x.clear();
    pan_y.clear();
    chainVec.clear();
    mm_chain_labels_.clear();
    if (panel) panel->reset_entries();

    // query complex 의 모든 체인 로드 (체인 목록은 같은 query 의 모든 hit 가 동일 → front 사용)
    const std::vector<std::string>& qChains = mm_hits_by_query_[qc].front().qChains;
    load_complex_chains(qc, qChains, mm_query_is_db_, query_db_reader_,
                        mm_query_source_, mm_query_is_dir_, multi_query_show_structure_);
    mm_query_chain_count_ = (int)data.size();
    if (mm_query_chain_count_ == 0) return;

    set_tmatrix();
    normalize_complex();
    update_total_len_ca();

    mm_current_hit_idx_ = -1;
    load_multimer_hit(+1);
}

void Screen::load_multimer_hit(int delta) {
    if (mm_current_query_idx_ < 0) return;
    const std::string& qc = mm_query_complexes_[mm_current_query_idx_];
    auto it = mm_hits_by_query_.find(qc);
    if (it == mm_hits_by_query_.end() || it->second.empty()) return;
    const std::vector<MultimerHit>& hits = it->second;

    int new_idx;
    if (mm_current_hit_idx_ < 0) {
        new_idx = (delta >= 0) ? 0 : (int)hits.size() - 1;
    } else {
        new_idx = mm_current_hit_idx_ + delta;
        new_idx = std::max(0, std::min((int)hits.size() - 1, new_idx));
        if (new_idx == mm_current_hit_idx_) return;
    }
    mm_current_hit_idx_ = new_idx;
    const MultimerHit& h = hits[new_idx];

    // 이전 target 체인 제거 (query 체인 보존)
    while ((int)data.size() > mm_query_chain_count_) {
        delete data.back();
        data.pop_back();
        pan_x.pop_back();
        pan_y.pop_back();
    }
    while ((int)chainVec.size() > mm_query_chain_count_) chainVec.pop_back();
    while ((int)mm_chain_labels_.size() > mm_query_chain_count_) mm_chain_labels_.pop_back();

    // target complex 체인 로드 + complex U/T 적용
    const int before = (int)data.size();
    load_complex_chains(h.tComplex, h.tChains, fs_db_reader_.is_open(), fs_db_reader_,
                        fs_db_path, true, screen_show_structure);
    if (h.has_transform) {
        for (int i = before; i < (int)data.size(); i++) {
            transform_target_chain(i, h.U, h.T);
        }
    }

    if (is_aligned_mode(screen_mode)) {
        // 멀티머 _report 에는 정렬 문자열이 없어 거리 판정만 가능하다. align-fs 는
        // 진입 전에 거부되므로(structty.cpp) 여기서는 폴백 사실만 패널에 남긴다.
        for (int qi = 0; qi < mm_query_chain_count_; qi++) {
            for (int ti = mm_query_chain_count_; ti < (int)data.size(); ti++) {
                if (data[qi] && data[ti]) data[qi]->compute_aligned_regions_nn(*data[ti], 4.0f);
            }
        }
        if (panel) panel->set_align_method("nearest-nbr");
    }

    // 패널 재구성: query+target 전체 entry + hit 정보
    if (panel) {
        const std::vector<std::string> q_tms = split_csv(h.qChainTms);
        const std::vector<std::string> t_tms = split_csv(h.tChainTms);
        panel->reset_entries();
        for (size_t i = 0; i < data.size(); i++) {
            std::string label = (i < mm_chain_labels_.size())
                                ? mm_chain_labels_[i] : data[i]->get_file_name();
            const bool is_query = (int)i < mm_query_chain_count_;
            const std::vector<std::string>& tms = is_query ? q_tms : t_tms;
            const size_t ci = is_query ? i : i - (size_t)mm_query_chain_count_;
            const size_t expect = is_query ? (size_t)mm_query_chain_count_
                                           : data.size() - (size_t)mm_query_chain_count_;
            if (tms.size() == expect && ci < tms.size()) {
                label += "  TM " + tms[ci];
            }
            panel->add_panel_info(label,
                                  data[i]->get_chain_length(),
                                  data[i]->get_residue_count(),
                                  is_query ? 0 : 1, (int)i);
        }
        FoldseekHitInfo fi;
        fi.valid = true;
        fi.multimer = true;
        fi.current_idx = mm_current_hit_idx_ + 1;
        fi.total_hits = (int)hits.size();
        fi.target = h.tComplex;
        fi.qtmscore = h.qTMScore;
        fi.ttmscore = h.tTMScore;
        fi.qcov = h.qComplexCov;
        fi.tcov = h.tComplexCov;
        fi.interface_lddt = h.interfaceLddt;
        fi.ass_id = h.assId;
        fi.query_idx = mm_current_query_idx_ + 1;
        fi.total_queries = (int)mm_query_complexes_.size();
        panel->set_foldseek_hit_info(fi);
    }
}

void Screen::switch_multimer_query(int delta) {
    if ((int)mm_query_complexes_.size() < 2) return;
    int new_idx = mm_current_query_idx_ + delta;
    new_idx = std::max(0, std::min((int)mm_query_complexes_.size() - 1, new_idx));
    if (new_idx == mm_current_query_idx_) return;
    activate_multimer_query(new_idx);
}

void Screen::set_tmatrix() {
    free_tmatrix();  // release previous allocation (multi-query re-setup)
    size_t filenum = data.size();
    vectorpointer = new float*[filenum];
    for (size_t i = 0; i < filenum; i++) {
        vectorpointer[i] = new float[3];
        vectorpointer[i][0] = 0;
        vectorpointer[i][1] = 0;
        vectorpointer[i][2] = 0;
    }
    vectorpointer_len_ = filenum;
}

void Screen::set_chainfile(const std::string& chainfile, int filesize) {
    for (int i = 0; i < filesize; i++) {
        chainVec.push_back("-");
    }
    if (chainfile.empty()) return;

    std::ifstream file(chainfile);
    if (!file.is_open()) {
        std::cerr << "Failed to open chainfile\n";
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        int index;
        std::string chainlist;
        iss >> index >> chainlist;
        if (index >= filesize) {
            std::cout << "Index in chain file should be 0..filenum-1\n";
            continue;
        }
        chainVec[index] = chainlist;
    }
    file.close();
}

void Screen::normalize_proteins() {
    for (size_t i = 0; i < data.size(); i++) {
        auto* p = data[i];
        // Cα-only proteins (loaded from a Foldseek DB via load_from_ca, e.g. the
        // query-from-DB path) already have init_atoms/screen_atoms populated and
        // have no plaintext file to parse — skip the gemmi load_data step.
        if (!p->is_ca_only()) {
            p->load_data(vectorpointer[i], yesUT);
        }
        // 적용된 체인 필터를 파일명 뒤에 표기 (좁으면 Panel 이 잘라낸다)
        std::string label = p->get_file_name();
        if (i < chainVec.size() && chainVec[i] != "-") {
            label += " [" + chainVec[i] + "]";
        }
        panel->add_panel_info(label,
                              p->get_chain_length(),
                              p->get_residue_count());
    }

    global_bb = BoundingBox();
    for (auto* p : data) {
        p->set_bounding_box();
        global_bb = global_bb + p->get_bounding_box();
    }

    float max_ext = std::max(global_bb.max_x - global_bb.min_x,
                             global_bb.max_y - global_bb.min_y);
    max_ext = std::max(max_ext, global_bb.max_z - global_bb.min_z);
    float scale = (max_ext > 0.f) ? (2.0f / max_ext) : 1.0f;

    // hit 변환 공식(load_next_hit / apply_hit_transform / transform_target_chain)은
    // "query 는 norm_c* 만큼 shift 된 뒤 norm_scale 로 스케일됐다" 를 전제한다.
    // 따라서 실제로 shift 에 사용한 값을 그대로 norm_c* 로 남긴다.
    float applied_cx = 0.0f;
    float applied_cy = 0.0f;
    float applied_cz = 0.0f;

    {
        for (auto* p : data) {
            // set_scale() 이 bounding box 중심을 cx/cy/cz 에 채운다. 따라서 center_shift 는
            // set_scale() 이후에 계산해야 한다. (이전에는 먼저 캡처해서 생성자 초기값 0 이
            // 들어가 중심 이동이 무효였고, 체인 단위 입력이 화면 밖으로 밀렸다)
            p->set_scale(scale);
            float center_shift[3] = { -p->cx, -p->cy, -p->cz };
            p->do_shift(center_shift);
            p->do_scale(scale);
        }

        // 단백질마다 자기 centroid 로 shift 한다. hit 변환의 기준은 query(data[0])
        // 이므로 그 값을 남긴다. set_scale() 이 이미 실행됐으므로 data[0]->cx 는
        // shift 에 사용한 centroid 와 같은 값이다.
        if (!data.empty() && data[0]) {
            applied_cx = data[0]->cx;
            applied_cy = data[0]->cy;
            applied_cz = data[0]->cz;
        }
    }

    depth_calibrated = false;

    // 기능 3: 정규화 파라미터 저장 (hit 탐색 시 동일 스케일 적용)
    norm_scale = scale;
    norm_cx = applied_cx;
    norm_cy = applied_cy;
    norm_cz = applied_cz;

    // 정규화 후 씬 중심이 원점이므로 회전 pivot 은 원점이다
    // (normalize_complex() 의 멀티머 경로와 동일).
    rot_pivot_[0] = 0.f;
    rot_pivot_[1] = 0.f;
    rot_pivot_[2] = 0.f;

    float radius = compute_scene_radius_from_render_positions(data);

    focal_offset = std::clamp(2.5f * radius + 1.0f, 2.0f, 8.0f);

    int n = (int)pan_x.size();
    int cols, rows;
    switch (n) {
        case 1:  cols=1; rows=1; break;
        case 2:  cols=2; rows=1; break;
        case 3:  cols=3; rows=1; break;
        case 4:  cols=2; rows=2; break;
        case 5:  cols=3; rows=2; break;
        case 6:  cols=3; rows=2; break;
        case 7:  cols=3; rows=3; break;
        case 8:  cols=3; rows=3; break;
        default: cols=3; rows=(n+2)/3; break;
    }

    if (n > 1) {
        int max_dim = std::max(cols, rows);
        float step      = (max_dim == 2) ? 0.75f : 0.5f;
        float foc_scale = (max_dim == 2) ? 0.8f  : 0.6f;
        focal_offset *= max_dim * foc_scale;
        for (int i = 0; i < n; i++) {
            int col = i % cols;
            int row = i / cols;
            pan_x[i] =  (col - (cols - 1) / 2.0f) * step;
            pan_y[i] = -((row - (rows - 1) / 2.0f) * step);
        }
    } else {
        for (int i = 0; i < n; i++) {
            pan_x[i] = 0.0f;
            pan_y[i] = 0.0f;
        }
    }
}

void Screen::calibrate_depth_baseline_first_view() {
    const float nearPlane = 0.05f;
    float fovRads = 1.0f / std::tan((FOV / zoom_level) * 0.5f / 180.0f * PI);

    float minz = std::numeric_limits<float>::infinity();
    float maxz = -std::numeric_limits<float>::infinity();
    bool any = false;

    for (size_t ii = 0; ii < data.size(); ii++) {
        Protein* target = data[ii];

        for (const auto& [chainID, chain_atoms] : target->get_atoms()) {
            if (chain_atoms.empty()) continue;
            int num_atoms = target->get_chain_length(chainID);

            for (int i = 0; i < num_atoms; ++i) {
                float* position = chain_atoms[i].get_position();
                float x = position[0];
                float y = position[1];
                // 카메라 Z 부호 통일: Renderer::project_and_fill() 과 동일한 관례 사용
                // (회전행렬은 +Z가 시청자 쪽인 표준 오른손좌표계를 가정하므로 부호 반전 필요).
                float z = -position[2] + focal_offset;

                if (z < nearPlane) continue;

                float projectedX = (x / z) * fovRads + pan_x[ii];
                float projectedY = (y / z) * fovRads + pan_y[ii];

                int screenX = (int)((projectedX + 1.0f) * 0.5f * screen_width);
                int screenY = (int)((1.0f - projectedY) * 0.5f * screen_height);

                if (screenX < 0 || screenX >= screen_width || screenY < 0 || screenY >= screen_height) continue;

                minz = std::min(minz, z);
                maxz = std::max(maxz, z);
                any = true;
            }
        }
    }

    if (!any) {
        depth_base_min_z = focal_offset;
        depth_base_max_z = focal_offset + 1.0f;
    } else {
        const float minSpan = 0.3f;
        if ((maxz - minz) < minSpan) maxz = minz + minSpan;
        depth_base_min_z = minz;
        depth_base_max_z = maxz;
    }

    depth_calibrated = true;
}

int Screen::intern_chain(const std::string& id) {
    auto it = chain_intern_.find(id);
    if (it != chain_intern_.end()) return it->second;
    int idx = (int)chain_names_.size();
    chain_names_.push_back(id);
    chain_intern_.emplace(id, idx);
    return idx;
}

const std::string& Screen::chain_name(int idx) const {
    static const std::string empty;
    if (idx < 0 || idx >= (int)chain_names_.size()) return empty;
    return chain_names_[idx];
}

std::vector<RenderAtom> Screen::to_render_atoms() {
    std::vector<RenderAtom> result;
    result.reserve(50000);
    for (size_t ii = 0; ii < data.size(); ++ii) {
        Protein* target = data[ii];
        for (const auto& [chainID, chain_atoms] : target->get_atoms()) {
            if (chain_atoms.empty()) continue;
            for (const auto& a : chain_atoms) {
                RenderAtom ra;
                ra.x = a.x; ra.y = a.y; ra.z = a.z;
                ra.structure          = a.get_structure();
                ra.bfactor            = a.bfactor;
                ra.is_interface       = a.is_interface;
                ra.is_aligned         = a.is_aligned;
                ra.conservation_score = a.conservation_score;
                ra.residue_number     = a.residue_number;
                strncpy(ra.residue_name, a.residue_name.c_str(), 3);
                ra.residue_name[3]    = '\0';
                ra.chain_id           = intern_chain(chainID);
                if (multimer_mode_) {
                    ra.protein_index  = ((int)ii < mm_query_chain_count_) ? 0 : 1;
                    ra.chain_color_id = (int)ii;
                } else {
                    ra.protein_index  = (int)ii;
                }
                ra.pan_x              = pan_x[ii];
                ra.pan_y              = pan_y[ii];
                result.push_back(std::move(ra));
            }
        }
    }
    return result;
}


void Screen::draw_screen(bool no_panel) {
    const bool prof = (bm && bm->enabled);
    auto t0 = Benchmark::clock::now();

    calibrate_depth_baseline_first_view();
    auto t_cal = Benchmark::clock::now();

    renderer_.set_depth_params(focal_offset, zoom_level, depth_base_min_z, depth_base_max_z);

    std::vector<RenderAtom> render_atoms = to_render_atoms();
    auto t_tra = Benchmark::clock::now();

    renderer_.render(render_atoms);
    auto t_ren = Benchmark::clock::now();

    Terminal::Size term_sz = Terminal::get_size();
    int rows = term_sz.rows;
    int cols = term_sz.cols;
    int panel_cols = std::min(cols, screen_width);

    // compact_level 자동 결정: 패널이 터미널 높이의 40% 이하가 되도록
    int max_panel_h = rows * 2 / 5;  // 40%
    int compact = 0;
    int panel_h = panel->get_height_for_width(panel_cols, 0);
    if (panel_h > max_panel_h) {
        compact = 1;
        panel_h = panel->get_height_for_width(panel_cols, 1);
    }
    if (panel_h > max_panel_h) {
        compact = 2;
        panel_h = panel->get_height_for_width(panel_cols, 2);
    }
    if (panel_h > max_panel_h) {
        compact = 3;
        panel_h = panel->get_height_for_width(panel_cols, 3);
    }
    if (panel_h > rows) panel_h = rows;

    int offset = 0;
    if (!no_panel) {
        offset += panel_h;
    }
    if (offset > rows) offset = rows;

    auto t_layout = Benchmark::clock::now();

    Terminal::cursor_home();
    print_screen(offset);
    auto t_print = Benchmark::clock::now();

    int start_row = rows;
    if (!no_panel) {
        start_row -= panel_h;
    }
    if (start_row < 0) start_row = 0;

    if (!no_panel){
        panel->draw_panel(start_row, 0, panel_h, panel_cols, compact);
    }
    // 기능 6: 패널 위치 저장 (hover 갱신 시 부분 재렌더링에 사용)
    last_panel_h         = no_panel ? 0 : panel_h;
    last_panel_start_row = no_panel ? rows : start_row;
    last_panel_cols      = panel_cols;
    last_no_panel        = no_panel;
    fflush(stdout);

    auto t1 = Benchmark::clock::now();
    int64_t render_dt_ms = Benchmark::ms_since(t0, t1);

    if (prof) {
        // 서브페이즈 계측 (µs). dt_ms 컬럼에 µs 값, num_ca 컬럼에 보조 수치.
        bm->log("ph_calibrate", -1, Benchmark::us_since(t0,       t_cal),    total_len_ca);
        bm->log("ph_torender",  -1, Benchmark::us_since(t_cal,    t_tra),    total_len_ca);
        bm->log("ph_render",    -1, Benchmark::us_since(t_tra,    t_ren),    total_len_ca);
        // render 내부 3분할 (µs): clear / project_and_fill / zbuffer_resolve
        bm->log("ph_r_clear",   -1, renderer_.get_last_us_clear(), total_len_ca);
        bm->log("ph_r_fill",    -1, renderer_.get_last_us_fill(),  total_len_ca);
        bm->log("ph_r_zbuf",    -1, renderer_.get_last_us_zbuf(),  total_len_ca);
        bm->log("ph_layout",    -1, Benchmark::us_since(t_ren,    t_layout), total_len_ca);
        bm->log("ph_print",     -1, Benchmark::us_since(t_layout, t_print),  total_len_ca);
        bm->log("ph_panel",     -1, Benchmark::us_since(t_print,  t1),       total_len_ca);
        // project_and_fill 이 생성한 RenderPoint 총 개수 (num_ca 컬럼 재활용).
        bm->log("points",       -1, (int64_t)renderer_.get_last_point_count(), total_len_ca);

        if (!ttff_logged) {
            bm->log("ttff", -1, Benchmark::ms_since(bm->t0, t1), total_len_ca, (int64_t)data.size());
            ttff_logged = true;
        }
        bm->mark_frame_end(render_dt_ms, total_len_ca, (int64_t)data.size());
    }
}

void Screen::print_screen_braille(int y_offset) {
    // Render logicalPixels (2*W x 4*H) to terminal using Unicode Braille characters.
    // Each terminal cell covers a 2-wide x 4-tall sub-pixel block:
    //
    //   subcol=0  subcol=1          Braille dot numbering (bit positions):
    //   subrow=0: dot1(0)  dot4(3)
    //   subrow=1: dot2(1)  dot5(4)
    //   subrow=2: dot3(2)  dot6(5)
    //   subrow=3: dot7(6)  dot8(7)
    //
    // Unicode Braille: U+2800 + bitmask
    static const int dot_bits[2][4] = {
        {0, 1, 2, 6},  // left column  (subcol=0)
        {3, 4, 5, 7}   // right column (subcol=1)
    };

    Terminal::Size sz = Terminal::get_size();
    int rows = sz.rows;
    int cols = sz.cols;

    const int logical_w = screen_width * 2;
    const int logical_h = screen_height * 4;

    const auto& pixels = renderer_.get_pixels();

    // Accumulate ANSI output into a single buffer, then fwrite once per frame
    std::string out;
    out.reserve(static_cast<size_t>(screen_width) * screen_height * 32);

    for (int ty = 0; ty < screen_height; ++ty) {
        int row = ty - (y_offset / 2) - 3;
        if (row < 0) continue;
        if (row >= rows) break;

        // Erase this terminal row before drawing — prevents ghosting from prior frames
        char row_clear[32];
        int rn = snprintf(row_clear, sizeof(row_clear), "\033[%d;1H\033[K", row + 1);
        if (rn > 0) out.append(row_clear, static_cast<size_t>(rn));

        int max_cols = std::min(screen_width, cols);

        for (int tx = 0; tx < max_cols; ++tx) {
            int   bitmask       = 0;
            int   best_color_id = 0;
            float best_depth    = std::numeric_limits<float>::infinity();

            for (int sc = 0; sc < 2; ++sc) {
                for (int sr = 0; sr < 4; ++sr) {
                    int lx = tx * 2 + sc;
                    int ly = ty * 4 + sr;
                    if (lx >= logical_w || ly >= logical_h) continue;

                    const RenderPoint& lp = pixels[ly * logical_w + lx];
                    if (lp.color_id > 0) {
                        bitmask |= (1 << dot_bits[sc][sr]);
                        if (lp.depth < best_depth) {
                            best_depth    = lp.depth;
                            best_color_id = lp.color_id;
                        }
                    }
                }
            }

            if (bitmask > 0 && best_color_id > 0) {
                // Cursor position (1-based)
                char pos[24];
                int pn = snprintf(pos, sizeof(pos), "\033[%d;%dH", row + 1, tx + 1);
                if (pn > 0) out.append(pos, static_cast<size_t>(pn));

                // ANSI fg colour
                out += Palettes::palette_to_ansi_fg_str(best_color_id);

                // Braille UTF-8: U+2800+bitmask encoded as 3-byte UTF-8
                out += (char)0xE2;
                out += (char)(0xA0 | (bitmask >> 6));
                out += (char)(0x80 | (bitmask & 0x3F));

                // Reset colour
                out += "\033[0m";
            }
            // Empty cells: already cleared by \033[K at row start, no extra output needed
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
}

void Screen::print_screen(int y_offset) {
    print_screen_braille(y_offset);
}

void Screen::set_zoom_level(float zoom){
    if ((zoom_level + zoom > 0.5)&&(zoom_level + zoom < 50)){
        float f_old = 1.0f / std::tan((FOV / zoom_level) * 0.5f / 180.0f * PI);
        zoom_level += zoom;
        float f_new = 1.0f / std::tan((FOV / zoom_level) * 0.5f / 180.0f * PI);
        float r = f_new / f_old;
        for (size_t i = 0; i < pan_x.size(); i++) {
            pan_x[i] *= r;
            pan_y[i] *= r;
        }
    }
}

// 기능 6: 마우스 커서 위치의 잔기 정보 검색 → 패널 hover 섹션 부분 갱신
void Screen::update_hover_info(int mx, int my) {
    // 터미널 좌표 → 스크린 좌표 변환
    // print_screen() 공식: row = screen_i - (y_offset/2) - 3
    // 역변환: screen_i = terminal_row + (panel_h/2) + 3
    int panel_offset = last_panel_h / 2 + 3;

    const RenderPoint* best = nullptr;
    float best_depth = std::numeric_limits<float>::infinity();

    {
        int ty_screen = my + panel_offset;
        if (ty_screen >= 0 && ty_screen < screen_height && mx >= 0 && mx < screen_width) {
            const int logical_w = screen_width * 2;
            const int logical_h = screen_height * 4;
            const auto& pixels = renderer_.get_pixels();
            for (int sc = 0; sc < 2; ++sc) {
                for (int sr = 0; sr < 4; ++sr) {
                    int lx = mx * 2 + sc;
                    int ly = ty_screen * 4 + sr;
                    if (lx >= logical_w || ly >= logical_h) continue;
                    const RenderPoint& lp = pixels[ly * logical_w + lx];
                    if (lp.residue_number >= 0 && lp.depth < best_depth) {
                        best_depth = lp.depth;
                        best = &lp;
                    }
                }
            }
        }
    }

    if (best) {
        panel->set_hover_residue(chain_name(best->chainID), best->residue_name,
                                 best->residue_number, best->structure,
                                 best->bfactor, best->conservation_score);
    } else {
        panel->clear_hover_residue();
    }

    // 패널이 보이는 경우에만 hover 섹션을 부분 갱신
    if (!last_no_panel && last_panel_h > 0) {
        int hover_row = panel->get_last_hover_row();
        if (hover_row >= 0) {
            panel->draw_hover_section(hover_row, last_panel_cols);
            fflush(stdout);
        }
    }
}

bool Screen::handle_input(bool& needs_redraw) {
    int key = Terminal::read_key();
    return handle_input_impl(key, needs_redraw);
}

bool Screen::handle_input(int key) {
    bool dummy = true;
    return handle_input_impl(key, dummy);
}

bool Screen::handle_input_impl(int key, bool& needs_redraw) {
    needs_redraw = true;  // 기본: 전체 재렌더링 필요

    if (key == Terminal::KEY_MOUSE) {
        Terminal::MouseEvent ev;
        if (Terminal::read_mouse(ev)) {
            update_hover_info(ev.x, ev.y);
        }
        needs_redraw = false;  // 마우스 이동은 전체 재렌더링 불필요
        return true;
    }

    bool keep_show = true;

    auto pan_step_x = 2.0f * 4.0f / screen_width;
    auto pan_step_y = 2.0f * 2.0f / screen_height;

    auto apply_pan = [&](int idx, float dx, float dy){
        if (idx < 0 || idx >= (int)pan_x.size()) return;
        pan_x[idx] += dx;
        pan_y[idx] += dy;
    };

    if (bm && bm->enabled) bm->mark_event(key);
    Terminal::flush_input();
    switch(key){
        // select protein / complex
        case 48:
            structNum = -1;
            break;
        case 49:
        case 50:
        case 51:
        case 52:
        case 53:
        case 54:
        case 55:
        case 56:
        case 57:
            if (multimer_mode_) {
                // 1 = query complex (idx 0), 2 = target complex (idx 1); 3+ ignored
                int requested = key - 49;
                if (requested < 2) structNum = requested;
            } else {
                if (key - 48 <= (int)data.size()) structNum = key - 49;
            }
            break;
        // A, a (minus x-axis)
        case 65:
        case 97:
            if (multimer_mode_ && structNum >= 0) {
                auto [lo, hi] = mm_complex_range(structNum);
                for (int i = lo; i <= hi; i++) apply_pan(i, -pan_step_x, 0.0f);
            } else if (structNum != -1) {
                apply_pan(structNum, -pan_step_x, 0.0f);
            } else {
                for (int i = 0; i < (int)data.size(); i++) apply_pan(i, -pan_step_x, 0.0f);
            }
            break;
        // D, d (plus x-axis)
        case 68:
        case 100:
            if (multimer_mode_ && structNum >= 0) {
                auto [lo, hi] = mm_complex_range(structNum);
                for (int i = lo; i <= hi; i++) apply_pan(i, +pan_step_x, 0.0f);
            } else if (structNum != -1) {
                apply_pan(structNum, +pan_step_x, 0.0f);
            } else {
                for (int i = 0; i < (int)data.size(); i++) apply_pan(i, +pan_step_x, 0.0f);
            }
            break;
        // S, s (minus y-axis)
        case 83:
        case 115:
            if (multimer_mode_ && structNum >= 0) {
                auto [lo, hi] = mm_complex_range(structNum);
                for (int i = lo; i <= hi; i++) apply_pan(i, 0.0f, -pan_step_y);
            } else if (structNum != -1) {
                apply_pan(structNum, 0.0f, -pan_step_y);
            } else {
                for (int i = 0; i < (int)data.size(); i++) apply_pan(i, 0.0f, -pan_step_y);
            }
            break;
        // W, w (plus y-axis)
        case 87:
        case 119:
            if (multimer_mode_ && structNum >= 0) {
                auto [lo, hi] = mm_complex_range(structNum);
                for (int i = lo; i <= hi; i++) apply_pan(i, 0.0f, +pan_step_y);
            } else if (structNum != -1) {
                apply_pan(structNum, 0.0f, +pan_step_y);
            } else {
                for (int i = 0; i < (int)data.size(); i++) apply_pan(i, 0.0f, +pan_step_y);
            }
            break;

        // X, x (rotate x-centered)
        case 88:
        case 120:
            if (multimer_mode_ && structNum >= 0) {
                auto [lo, hi] = mm_complex_range(structNum);
                for (int i = lo; i <= hi; i++) data[i]->set_rotate(1, 0, 0);
            } else if (structNum != -1) {
                data[structNum]->set_rotate(1, 0, 0);
            } else if (yesUT) {
                float c = cos(PI / 48.0f), s = sin(PI / 48.0f);
                float m[9] = {1,0,0, 0,c,-s, 0,s,c};
                float neg_piv[3] = { -rot_pivot_[0], -rot_pivot_[1], -rot_pivot_[2] };
                for (auto* p : data) {
                    p->do_shift(neg_piv);
                    p->do_naive_rotation(m);
                    p->do_shift(rot_pivot_);
                }
            } else {
                for (int i = 0; i < (int)data.size(); i++) data[i]->set_rotate(1, 0, 0);
            }
            break;
        // Y, y (rotate y-centered)
        case 89:
        case 121:
            if (multimer_mode_ && structNum >= 0) {
                auto [lo, hi] = mm_complex_range(structNum);
                for (int i = lo; i <= hi; i++) data[i]->set_rotate(0, 1, 0);
            } else if (structNum != -1) {
                data[structNum]->set_rotate(0, 1, 0);
            } else if (yesUT) {
                float c = cos(PI / 48.0f), s = sin(PI / 48.0f);
                float m[9] = {c,0,s, 0,1,0, -s,0,c};
                float neg_piv[3] = { -rot_pivot_[0], -rot_pivot_[1], -rot_pivot_[2] };
                for (auto* p : data) {
                    p->do_shift(neg_piv);
                    p->do_naive_rotation(m);
                    p->do_shift(rot_pivot_);
                }
            } else {
                for (int i = 0; i < (int)data.size(); i++) data[i]->set_rotate(0, 1, 0);
            }
            break;
        // Z, z (rotate z-centered)
        case 90:
        case 122:
            if (multimer_mode_ && structNum >= 0) {
                auto [lo, hi] = mm_complex_range(structNum);
                for (int i = lo; i <= hi; i++) data[i]->set_rotate(0, 0, 1);
            } else if (structNum != -1) {
                data[structNum]->set_rotate(0, 0, 1);
            } else if (yesUT) {
                float c = cos(PI / 48.0f), s = sin(PI / 48.0f);
                float m[9] = {c,-s,0, s,c,0, 0,0,1};
                float neg_piv[3] = { -rot_pivot_[0], -rot_pivot_[1], -rot_pivot_[2] };
                for (auto* p : data) {
                    p->do_shift(neg_piv);
                    p->do_naive_rotation(m);
                    p->do_shift(rot_pivot_);
                }
            } else {
                for (int i = 0; i < (int)data.size(); i++) data[i]->set_rotate(0, 0, 1);
            }
            break;

        // F, f (zoom out)
        case 70:
        case 102:
            set_zoom_level(-0.3);
            break;   
        // R, R (zoom in)
        case 82:
        case 114:
            set_zoom_level(0.3);
            break;   

        // C, c (camera)
        case 67:
        case 99:     
        {
            camera->screenshot(renderer_.get_pixels(), screen_width * 2, screen_height * 4);
            break;
        }
        // N, n (next Foldseek hit / target complex)
        case 78:
        case 110:
            if (multimer_mode_) load_multimer_hit(+1);
            else if (!foldseek_hits.empty()) load_next_hit(+1);
            break;

        // P, p (prev Foldseek hit / target complex)
        case 80:
        case 112:
            if (multimer_mode_) load_multimer_hit(-1);
            else if (!foldseek_hits.empty()) load_next_hit(-1);
            break;

        // ] (next query — Step 5 multi-query / Step 7 multimer nav)
        case 93:
            if (multimer_mode_) switch_multimer_query(+1);
            else switch_query(+1);
            break;

        // [ (prev query)
        case 91:
            if (multimer_mode_) switch_multimer_query(-1);
            else switch_query(-1);
            break;

        // Q, q
        case 81:
        case 113:
            keep_show = false;
            break;

        default:
            break;
    }

    return keep_show;
}

// 기능 1: 로드된 모든 Protein에 대해 inter-chain interface를 계산
void Screen::compute_interface_all(float threshold) {
    for (Protein* p : data) {
        if (p) p->compute_interface(threshold);
    }
}

// 기능 5: 지정 protein에 conservation scores 적용
void Screen::apply_msa_conservation(int protein_idx, const std::vector<float>& scores) {
    if (protein_idx >= 0 && protein_idx < (int)data.size() && data[protein_idx]) {
        data[protein_idx]->apply_conservation_scores(scores);
    }
}

// 기능 4: 로드된 모든 Protein 쌍에 대해 nearest-neighbor 기반 정렬 잔기를 계산
void Screen::compute_aligned_all(float threshold) {
    if (data.size() < 2) {
        // 단일 단백질인 경우: 전체 잔기를 aligned로 표시
        for (Protein* p : data) {
            if (p) p->compute_aligned_regions_nn(*p, threshold);
        }
    } else {
        for (size_t i = 0; i < data.size(); ++i) {
            for (size_t j = i + 1; j < data.size(); ++j) {
                if (data[i] && data[j]) {
                    data[i]->compute_aligned_regions_nn(*data[j], threshold);
                }
            }
        }
    }
    // 정렬 문자열이 아니라 거리 기반 판정임을 패널에 정직하게 표시한다.
    // (기존에는 단백질이 1개일 때 early return 하면서 라벨이 누락됐다)
    if (panel) panel->set_align_method("nearest-nbr");
}

// 기능 4: -fs 기반 — Foldseek hit의 U/T transform을 지정 protein에 적용
// T_norm: 정규화 공간 T (screen_atoms 이동용)
// T_angstrom: Å 공간 T (init_atoms 갱신용). nullptr이면 init_atoms 미갱신.
//   → compute_aligned_regions_from_aln은 init_atoms 기준 거리 비교를 하므로
//     반드시 Å 공간 T를 전달해야 aligned 색상이 올바르게 동작함.
void Screen::apply_foldseek_transform(int protein_idx, const float* U_flat,
                                      const float* T_norm, const float* T_angstrom) {
    if (protein_idx < 0 || protein_idx >= (int)data.size() || !data[protein_idx]) return;
    // screen_atoms에 적용 (시각적 정렬)
    data[protein_idx]->do_naive_rotation(const_cast<float*>(U_flat));
    data[protein_idx]->do_shift(const_cast<float*>(T_norm));
    // init_atoms에 Å 공간 T로 적용 (거리 비교 기준)
    if (T_angstrom) {
        data[protein_idx]->apply_ut_to_init_atoms(U_flat, T_angstrom);
    }
    yesUT = true;

    // pan 초기화 — overlay 모드로 전환 (grid layout 해제)
    for (size_t i = 0; i < pan_x.size(); i++) {
        pan_x[i] = 0.0f;
        pan_y[i] = 0.0f;
    }

    // focal_offset 재계산 (bounding box 변경 반영)
    float radius = compute_scene_radius_from_render_positions(data);
    focal_offset = std::clamp(2.5f * radius + 1.0f, 2.0f, 8.0f);
    depth_calibrated = false;
}

// 기능 4: -fs 기반 — alignment string으로 aligned 잔기 계산 (protein0 vs protein1)
void Screen::compute_aligned_from_aln(const std::string& qaln, const std::string& taln,
                                      int q_start, int t_start,
                                      float threshold, bool skip_distance_check) {
    if ((int)data.size() < 2 || !data[0] || !data[1]) return;
    data[0]->compute_aligned_regions_from_aln(*data[1], qaln, taln,
                                             q_start, t_start,
                                             threshold, skip_distance_check);
}

// 기능 4: 패널에 정렬 방식 표시 설정
void Screen::set_align_method(const std::string& method) {
    if (panel) panel->set_align_method(method);
}

// ── 기능 3: Foldseek hit 탐색 ──────────────────────────────────────────────────

// 3×3 대칭 행렬 Jacobi 고유값 분해
// a: 입력 (수정됨), d: 고유값, v: 고유벡터 (열 = 고유벡터)
static void jacobi3_sym(float a[3][3], float d[3], float v[3][3]) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) v[i][j] = (i == j) ? 1.f : 0.f;
    for (int iter = 0; iter < 100; ++iter) {
        int p = 0, q = 1;
        float max_val = std::abs(a[0][1]);
        if (std::abs(a[0][2]) > max_val) { max_val = std::abs(a[0][2]); p = 0; q = 2; }
        if (std::abs(a[1][2]) > max_val) { max_val = std::abs(a[1][2]); p = 1; q = 2; }
        if (max_val < 1e-9f) break;
        float diff = a[q][q] - a[p][p];
        float t;
        if (std::abs(diff) < 1e-9f * std::abs(a[p][q])) {
            t = 1.f;
        } else {
            float phi = diff / (2.f * a[p][q]);
            t = (phi > 0) ? (1.f / (phi + std::sqrt(1.f + phi*phi)))
                           : (1.f / (phi - std::sqrt(1.f + phi*phi)));
        }
        float c = 1.f / std::sqrt(1.f + t*t);
        float s = t * c;
        a[p][p] -= t * a[p][q];
        a[q][q] += t * a[p][q];
        a[p][q] = a[q][p] = 0.f;
        int r3 = 3 - p - q;
        float apr = a[p][r3], aqr = a[q][r3];
        a[p][r3] = a[r3][p] =  c * apr - s * aqr;
        a[q][r3] = a[r3][q] =  s * apr + c * aqr;
        for (int i = 0; i < 3; i++) {
            float vp = v[i][p], vq = v[i][q];
            v[i][p] = c * vp - s * vq;
            v[i][q] = s * vp + c * vq;
        }
    }
    for (int i = 0; i < 3; i++) d[i] = a[i][i];
}

static float mat3_det_flat(const float m[9]) {
    return m[0]*(m[4]*m[8]-m[5]*m[7])
          -m[1]*(m[3]*m[8]-m[5]*m[6])
          +m[2]*(m[3]*m[7]-m[4]*m[6]);
}

// Kabsch 알고리즘: U*Q[i]+T ≈ P[i] 최소화
// P: 참조 좌표 (Nx3 flat), Q: 회전 대상 좌표 (Nx3 flat), N: 쌍 수
// 결과: U[9] (row-major 3×3), T[3]
static void kabsch(const std::vector<float>& P, const std::vector<float>& Q,
                   int N, float U[9], float T[3]) {
    // 단위 행렬로 초기화
    for (int i = 0; i < 9; i++) U[i] = (i % 4 == 0) ? 1.f : 0.f;
    T[0] = T[1] = T[2] = 0.f;
    if (N < 3) return;

    // 무게중심
    float Pc[3] = {}, Qc[3] = {};
    for (int i = 0; i < N; i++) {
        Pc[0] += P[i*3]; Pc[1] += P[i*3+1]; Pc[2] += P[i*3+2];
        Qc[0] += Q[i*3]; Qc[1] += Q[i*3+1]; Qc[2] += Q[i*3+2];
    }
    for (int k = 0; k < 3; k++) { Pc[k] /= N; Qc[k] /= N; }

    // H = Σ (Q[i]-Qc)(P[i]-Pc)^T  [3×3]
    float H[3][3] = {};
    for (int i = 0; i < N; i++) {
        float q[3] = { Q[i*3]-Qc[0], Q[i*3+1]-Qc[1], Q[i*3+2]-Qc[2] };
        float p[3] = { P[i*3]-Pc[0], P[i*3+1]-Pc[1], P[i*3+2]-Pc[2] };
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                H[r][c] += q[r] * p[c];
    }

    // H^T H (대칭) → 고유값 분해로 V (오른쪽 특이벡터) 계산
    float HtH[3][3] = {};
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            for (int k = 0; k < 3; k++)
                HtH[r][c] += H[k][r] * H[k][c];

    float d[3];
    float Vmat[3][3];
    jacobi3_sym(HtH, d, Vmat);

    // 고유값 내림차순 정렬
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2-i; j++) {
            if (d[j] < d[j+1]) {
                std::swap(d[j], d[j+1]);
                for (int k = 0; k < 3; k++) std::swap(Vmat[k][j], Vmat[k][j+1]);
            }
        }
    }

    // 특이값
    float sigma[3];
    for (int i = 0; i < 3; i++) sigma[i] = (d[i] > 0.f) ? std::sqrt(d[i]) : 0.f;

    // 왼쪽 특이벡터: Umat[:,i] = H * Vmat[:,i] / sigma[i]
    float Umat[3][3] = {};
    for (int i = 0; i < 3; i++) {
        if (sigma[i] < 1e-7f) continue;
        for (int r = 0; r < 3; r++) {
            for (int k = 0; k < 3; k++) Umat[r][i] += H[r][k] * Vmat[k][i];
            Umat[r][i] /= sigma[i];
        }
    }
    // 退화 처리: 세 번째 열 = 교차곱
    if (sigma[2] < 1e-7f) {
        Umat[0][2] = Umat[1][0]*Umat[2][1] - Umat[2][0]*Umat[1][1];
        Umat[1][2] = Umat[2][0]*Umat[0][1] - Umat[0][0]*Umat[2][1];
        Umat[2][2] = Umat[0][0]*Umat[1][1] - Umat[1][0]*Umat[0][1];
    }

    // det(V*U^T) 반사 보정
    float VUt[9];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            VUt[r*3+c] = 0;
            for (int k = 0; k < 3; k++) VUt[r*3+c] += Vmat[r][k] * Umat[c][k];
        }
    float det_sign = (mat3_det_flat(VUt) < 0) ? -1.f : 1.f;

    // R = V * diag(1,1,det_sign) * U^T
    float diag[3] = {1.f, 1.f, det_sign};
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            U[r*3+c] = 0;
            for (int k = 0; k < 3; k++) U[r*3+c] += Vmat[r][k] * diag[k] * Umat[c][k];
        }

    // T = Pc - R*Qc
    for (int r = 0; r < 3; r++) {
        T[r] = Pc[r];
        for (int c = 0; c < 3; c++) T[r] -= U[r*3+c] * Qc[c];
    }
}

void Screen::set_foldseek_hits(const std::vector<FoldseekHit>& hits) {
    foldseek_hits = hits;
    current_hit_idx = -1;
    // 패널에 총 hit 수 표시
    if (panel && !hits.empty()) {
        FoldseekHitInfo fi;
        fi.valid = true;
        fi.current_idx = 0;
        fi.total_hits = (int)hits.size();
        fi.query_idx = (current_query_idx_ >= 0) ? current_query_idx_ + 1 : 0;
        fi.total_queries = (int)query_ids_.size();
        panel->set_foldseek_hit_info(fi);
    }
}

void Screen::set_fs_db_path(const std::string& path) {
    fs_db_path = path;
}

bool Screen::open_foldseek_db(const std::string& db_base_path) {
    if (fs_db_reader_.open(db_base_path)) {
        std::cout << "  Foldseek DB opened: " << db_base_path << "\n";
        return true;
    }
    std::cerr << "Warning: Failed to open Foldseek DB: " << db_base_path << "\n";
    return false;
}

bool Screen::prepare_foldseek_db(const std::vector<FoldseekHit>& hits) {
    if (!fs_db_reader_.is_open()) return false;
    std::unordered_set<std::string> target_ids;
    for (const auto& hit : hits) {
        target_ids.insert(hit.target);
    }
    return fs_db_reader_.prepare(target_ids);
}

void Screen::load_next_hit(int delta) {
    if (foldseek_hits.empty()) return;

    int new_idx;
    if (current_hit_idx < 0) {
        // 초기 로드: delta > 0 이면 0번, delta < 0 이면 마지막
        new_idx = (delta >= 0) ? 0 : (int)foldseek_hits.size() - 1;
    } else {
        new_idx = current_hit_idx + delta;
        new_idx = std::max(0, std::min((int)foldseek_hits.size() - 1, new_idx));
        if (new_idx == current_hit_idx) return;
    }
    current_hit_idx = new_idx;

    const FoldseekHit& hit = foldseek_hits[current_hit_idx];

    std::string status_msg;
    FoldseekHitInfo fs_info;
    fs_info.valid      = true;
    fs_info.current_idx = current_hit_idx + 1;
    fs_info.total_hits  = (int)foldseek_hits.size();
    fs_info.target      = hit.target;
    fs_info.evalue      = hit.evalue;
    fs_info.prob        = hit.prob;
    fs_info.lddt        = hit.lddt;
    fs_info.qtmscore    = hit.qtmscore;
    fs_info.ttmscore    = hit.ttmscore;
    fs_info.query_idx     = (current_query_idx_ >= 0) ? current_query_idx_ + 1 : 0;
    fs_info.total_queries = (int)query_ids_.size();

    // 기존 target protein (index 1+) 제거
    while ((int)data.size() > 1) {
        delete data.back();
        data.pop_back();
        pan_x.pop_back();
        pan_y.pop_back();
    }
    while ((int)chainVec.size() > 1) chainVec.pop_back();

    Protein* target_protein = nullptr;
    bool loaded_from_db = false;

    // ★ DB 직접 읽기 경로
    if (fs_db_reader_.is_open()) {
        std::vector<float> coords;
        std::string aa_seq;
        size_t n_res = fs_db_reader_.read_entry(hit.target, coords, aa_seq);
        if (n_res > 0) {
            chainVec.push_back("-");
            target_protein = new Protein(hit.target, screen_show_structure);
            target_protein->load_from_ca(coords, n_res, aa_seq);
            data.push_back(target_protein);
            pan_x.push_back(0.0f);
            pan_y.push_back(0.0f);
            loaded_from_db = true;
            status_msg = "loaded from DB";
            if (screen_mode == "plddt") {
                status_msg += " (pLDDT not available)";
            }
        }
    }

    // Fallback: 기존 PDBDownloader 경로
    if (!loaded_from_db) {
        std::string file_path = PDBDownloader::resolve_target_file(
            hit.target, fs_db_path, status_msg);

        if (file_path.empty()) {
            fs_info.status_msg = status_msg;
            if (panel) panel->set_foldseek_hit_info(fs_info);
            return;
        }

        DBType db_type = PDBDownloader::detect_db_type(hit.target);
        std::string chain_filter = PDBDownloader::extract_chain(hit.target, db_type);
        if (chain_filter == "-") {
            // foldseek 체인 accession(`<stem>_<chain>`)은 detect_db_type 이 Unknown 으로
            // 보기 때문에 extract_chain 이 "-" 를 준다. 해석된 파일명으로 한 번 더 시도.
            // 어느 체인이 걸렸는지는 패널의 `[B-2]` 라벨에 나오므로 로그는 남기지 않는다
            // (hit 을 넘길 때마다 찍히면 화면이 밀린다).
            const std::string from_acc = chain_from_accession(hit.target, file_path);
            if (from_acc != "-") {
                chain_filter = from_acc;
            }
        }
        chainVec.push_back(chain_filter);

        target_protein = new Protein(file_path, chain_filter, screen_show_structure);
        data.push_back(target_protein);
        pan_x.push_back(0.0f);
        pan_y.push_back(0.0f);

        float tvec[3] = {0.f, 0.f, 0.f};
        target_protein->load_data(tvec, false);
    }

    fs_info.status_msg = status_msg;
    if (panel) panel->set_foldseek_hit_info(fs_info);

    // 패널 entry 갱신 (적용된 체인 필터를 파일명 뒤에 표기)
    if (panel) {
        std::string label = target_protein->get_file_name();
        if (chainVec.size() > 1 && chainVec[1] != "-") {
            label += " [" + chainVec[1] + "]";
        }
        panel->update_entry(1, label,
                            target_protein->get_chain_length(),
                            target_protein->get_residue_count());
    }

    // query와 동일한 norm_scale 로 target 정규화
    target_protein->set_bounding_box();
    target_protein->set_scale(norm_scale);
    float t_cx = target_protein->cx;
    float t_cy = target_protein->cy;
    float t_cz = target_protein->cz;
    float t_shift[3] = { -t_cx, -t_cy, -t_cz };
    target_protein->do_shift(t_shift);
    target_protein->do_scale(norm_scale);

    // U/T transform 계산
    float U[9] = {1,0,0, 0,1,0, 0,0,1};
    float T[3]  = {0,0,0};
    float T_ang[3] = {0,0,0};  // Å 공간 T (init_atoms 갱신용)
    bool computed_transform = false;
    std::string align_method_str;

    if (hit.has_transform) {
        // 29컬럼 포맷: U/T를 정규화 공간으로 변환
        // target_norm = (target_orig - t_centroid) * norm_scale
        // alns_norm   = (alns_orig   - q_centroid) * norm_scale
        // U_norm      = U_fs (회전은 스케일 불변)
        // T_norm      = (U_fs*t_centroid + T_fs - q_centroid) * norm_scale
        const float* Uf = hit.U;
        const float* Tf = hit.T;
        float Utc[3] = {
            Uf[0]*t_cx + Uf[1]*t_cy + Uf[2]*t_cz,
            Uf[3]*t_cx + Uf[4]*t_cy + Uf[5]*t_cz,
            Uf[6]*t_cx + Uf[7]*t_cy + Uf[8]*t_cz
        };
        for (int i = 0; i < 9; i++) U[i] = Uf[i];
        T[0] = (Utc[0] + Tf[0] - norm_cx) * norm_scale;
        T[1] = (Utc[1] + Tf[1] - norm_cy) * norm_scale;
        T[2] = (Utc[2] + Tf[2] - norm_cz) * norm_scale;
        // init_atoms용: Foldseek 원본 Å 공간 T (hit.T)를 그대로 사용
        T_ang[0] = Tf[0]; T_ang[1] = Tf[1]; T_ang[2] = Tf[2];
        computed_transform = true;
        align_method_str = "aln-string";

    } else if (hit.is_alis_format && !hit.alns.empty() && hit.has_aln) {
        // alis 21컬럼 포맷: Kabsch SVD로 U/T 역산
        // alns는 원래 query frame 좌표 → 정규화: (alns - q_centroid) * norm_scale
        // target CA 좌표는 이미 정규화됨 (do_shift + do_scale 후)

        // target CA atom 플랫 리스트 (정규화된 좌표)
        std::vector<std::array<float,3>> target_cas;
        for (const auto& [cid, chain] : target_protein->get_ca_atoms()) {
            for (const auto& atom : chain) {
                target_cas.push_back({(atom.x - t_cx) * norm_scale,
                                      (atom.y - t_cy) * norm_scale,
                                      (atom.z - t_cz) * norm_scale});
            }
        }

        // taln 순회: 비-갭 위치마다 쌍 수집
        std::vector<float> P_norm, Q_norm;
        int aln_idx = 0;  // alns 인덱스 (3 float/잔기)
        int t_seq_idx = hit.tstart - 1;  // target CA 0-based 인덱스

        for (size_t ai = 0; ai < hit.taln.size(); ai++) {
            if (hit.taln[ai] == '-') continue;
            if (t_seq_idx < (int)target_cas.size() &&
                aln_idx * 3 + 2 < (int)hit.alns.size()) {
                // 정규화된 alns
                P_norm.push_back((hit.alns[aln_idx*3]   - norm_cx) * norm_scale);
                P_norm.push_back((hit.alns[aln_idx*3+1] - norm_cy) * norm_scale);
                P_norm.push_back((hit.alns[aln_idx*3+2] - norm_cz) * norm_scale);
                // 정규화된 target CA
                Q_norm.push_back(target_cas[t_seq_idx][0]);
                Q_norm.push_back(target_cas[t_seq_idx][1]);
                Q_norm.push_back(target_cas[t_seq_idx][2]);
            }
            aln_idx++;
            t_seq_idx++;
        }

        int N = (int)std::min(P_norm.size(), Q_norm.size()) / 3;
        if (N >= 3) {
            kabsch(P_norm, Q_norm, N, U, T);
            // BUG-A 2단계: kabsch는 정규화 공간에서 계산됨.
            // init_atoms는 Å 공간이므로 T_Å = T_norm/norm_scale + q_centroid - U*t_centroid
            {
                const float q_cen[3] = {norm_cx, norm_cy, norm_cz};
                const float t_cen[3] = {t_cx, t_cy, t_cz};
                for (int r = 0; r < 3; r++) {
                    T_ang[r] = T[r] / norm_scale + q_cen[r];
                    for (int c = 0; c < 3; c++) T_ang[r] -= U[r*3+c] * t_cen[c];
                }
            }
            computed_transform = true;
            align_method_str = "kabsch-alns";
        }
    } else if (hit.has_aln && !hit.is_alis_format) {
        // 작업 1-A: 17컬럼 포맷 — qaln/taln 기반 Kabsch SVD
        // qaln/taln에서 양쪽 비-갭 위치의 CA 쌍을 수집 → 정규화 공간에서 Kabsch

        // query CA flat 리스트 (정규화된 screen_atoms)
        std::vector<std::array<float,3>> query_cas;
        for (const auto& [cid, chain] : data[0]->get_ca_atoms()) {
            for (const auto& atom : chain) {
                query_cas.push_back({(atom.x - norm_cx) * norm_scale,
                                      (atom.y - norm_cy) * norm_scale,
                                      (atom.z - norm_cz) * norm_scale});
            }
        }

        // target CA flat 리스트 (정규화된 screen_atoms)
        std::vector<std::array<float,3>> target_cas;
        for (const auto& [cid, chain] : target_protein->get_ca_atoms()) {
            for (const auto& atom : chain) {
                target_cas.push_back({(atom.x - t_cx) * norm_scale,
                                       (atom.y - t_cy) * norm_scale,
                                       (atom.z - t_cz) * norm_scale});
            }
        }

        // qaln/taln 순회하여 aligned CA 쌍 수집
        std::vector<float> P_norm, Q_norm;
        int q_idx = hit.qstart - 1;  // 0-based
        int t_idx = hit.tstart - 1;  // 0-based
        const int q_size = (int)query_cas.size();
        const int t_size = (int)target_cas.size();

        for (size_t ai = 0; ai < hit.qaln.size() && ai < hit.taln.size(); ai++) {
            bool q_gap = (hit.qaln[ai] == '-');
            bool t_gap = (hit.taln[ai] == '-');

            if (!q_gap && !t_gap) {
                if (q_idx < q_size && t_idx < t_size) {
                    P_norm.push_back(query_cas[q_idx][0]);
                    P_norm.push_back(query_cas[q_idx][1]);
                    P_norm.push_back(query_cas[q_idx][2]);
                    Q_norm.push_back(target_cas[t_idx][0]);
                    Q_norm.push_back(target_cas[t_idx][1]);
                    Q_norm.push_back(target_cas[t_idx][2]);
                }
            }

            if (!q_gap) q_idx++;
            if (!t_gap) t_idx++;
        }

        int N = (int)std::min(P_norm.size(), Q_norm.size()) / 3;
        if (N >= 3) {
            kabsch(P_norm, Q_norm, N, U, T);
            {
                const float q_cen[3] = {norm_cx, norm_cy, norm_cz};
                const float t_cen2[3] = {t_cx, t_cy, t_cz};
                for (int r = 0; r < 3; r++) {
                    T_ang[r] = T[r] / norm_scale + q_cen[r];
                    for (int c = 0; c < 3; c++) T_ang[r] -= U[r*3+c] * t_cen2[c];
                }
            }
            computed_transform = true;
            align_method_str = "kabsch-qaln";
        }
    }

    // 작업 1-B: fallback — 12컬럼 등 transform 정보 없는 경우 전체 CA 순서 매칭 Kabsch
    if (!computed_transform) {
        std::vector<std::array<float,3>> query_cas;
        for (const auto& [cid, chain] : data[0]->get_ca_atoms()) {
            for (const auto& atom : chain) {
                query_cas.push_back({(atom.x - norm_cx) * norm_scale,
                                      (atom.y - norm_cy) * norm_scale,
                                      (atom.z - norm_cz) * norm_scale});
            }
        }

        std::vector<std::array<float,3>> target_cas;
        for (const auto& [cid, chain] : target_protein->get_ca_atoms()) {
            for (const auto& atom : chain) {
                target_cas.push_back({(atom.x - t_cx) * norm_scale,
                                       (atom.y - t_cy) * norm_scale,
                                       (atom.z - t_cz) * norm_scale});
            }
        }

        int N = (int)std::min(query_cas.size(), target_cas.size());
        if (N >= 3) {
            std::vector<float> P_norm, Q_norm;
            P_norm.reserve(N * 3);
            Q_norm.reserve(N * 3);
            for (int i = 0; i < N; i++) {
                P_norm.push_back(query_cas[i][0]);
                P_norm.push_back(query_cas[i][1]);
                P_norm.push_back(query_cas[i][2]);
                Q_norm.push_back(target_cas[i][0]);
                Q_norm.push_back(target_cas[i][1]);
                Q_norm.push_back(target_cas[i][2]);
            }
            kabsch(P_norm, Q_norm, N, U, T);
            {
                const float q_cen[3] = {norm_cx, norm_cy, norm_cz};
                const float t_cen2[3] = {t_cx, t_cy, t_cz};
                for (int r = 0; r < 3; r++) {
                    T_ang[r] = T[r] / norm_scale + q_cen[r];
                    for (int c = 0; c < 3; c++) T_ang[r] -= U[r*3+c] * t_cen2[c];
                }
            }
            computed_transform = true;
            align_method_str = "kabsch-seq";
        }
    }

    if (computed_transform) {
        apply_foldseek_transform(1, U, T, T_ang);
    }

    // 정렬 방식 패널 갱신
    fs_info.align_method = align_method_str;
    if (panel) panel->set_foldseek_hit_info(fs_info);

    // align 계열 모드일 때 is_aligned 계산
    if (is_aligned_mode(screen_mode)) {
        if (screen_mode == "align-near") {
            // 거리 판정을 명시적으로 요구한 경우: 정렬 문자열이 있어도 쓰지 않는다.
            compute_aligned_all();
        } else if (hit.has_aln) {
            // qaln/taln 은 hit.qstart/hit.tstart 잔기에서 시작한다
            compute_aligned_from_aln(hit.qaln, hit.taln, hit.qstart, hit.tstart, 5.0f, true);
            set_align_method("aln-string");
        } else if (screen_mode == "align-fs") {
            // 폴백 금지. 진입 전에 거부되지만, 혼합 포맷에 대비해 아무것도 칠하지 않는다.
            set_align_method("none");
        } else {
            // 정렬 문자열이 없는 포맷(12컬럼 등): 최근접 이웃 폴백.
            // compute_aligned_all() 이 패널 라벨을 "nearest-nbr" 로 설정한다.
            compute_aligned_all();
        }
    }

    depth_calibrated = false;
}

// ── 기능 3: 이미 로드된 target에 hit transform 적용 ─────────────────────────
void Screen::apply_hit_transform(int target_protein_idx, const FoldseekHit& hit) {
    if (target_protein_idx < 0 || target_protein_idx >= (int)data.size() || !data[target_protein_idx]) return;
    if (!data[0]) return;

    Protein* target_protein = data[target_protein_idx];
    float t_cx = target_protein->cx;
    float t_cy = target_protein->cy;
    float t_cz = target_protein->cz;

    float U[9] = {1,0,0, 0,1,0, 0,0,1};
    float T[3]  = {0,0,0};
    float T_ang[3] = {0,0,0};
    bool computed_transform = false;

    if (hit.has_transform) {
        // 29컬럼 포맷: U/T 직접 사용
        const float* Uf = hit.U;
        const float* Tf = hit.T;
        float Utc[3] = {
            Uf[0]*t_cx + Uf[1]*t_cy + Uf[2]*t_cz,
            Uf[3]*t_cx + Uf[4]*t_cy + Uf[5]*t_cz,
            Uf[6]*t_cx + Uf[7]*t_cy + Uf[8]*t_cz
        };
        for (int i = 0; i < 9; i++) U[i] = Uf[i];
        T[0] = (Utc[0] + Tf[0] - norm_cx) * norm_scale;
        T[1] = (Utc[1] + Tf[1] - norm_cy) * norm_scale;
        T[2] = (Utc[2] + Tf[2] - norm_cz) * norm_scale;
        T_ang[0] = Tf[0]; T_ang[1] = Tf[1]; T_ang[2] = Tf[2];
        computed_transform = true;
    } else if (hit.is_alis_format && !hit.alns.empty() && hit.has_aln) {
        // 21컬럼 alis 포맷: Kabsch SVD
        std::vector<std::array<float,3>> target_cas;
        for (const auto& [cid, chain] : target_protein->get_ca_atoms()) {
            for (const auto& atom : chain) {
                target_cas.push_back({(atom.x - t_cx) * norm_scale,
                                       (atom.y - t_cy) * norm_scale,
                                       (atom.z - t_cz) * norm_scale});
            }
        }

        std::vector<float> P_norm, Q_norm;
        int aln_idx = 0;
        int t_seq_idx = hit.tstart - 1;

        for (size_t ai = 0; ai < hit.taln.size(); ai++) {
            if (hit.taln[ai] == '-') continue;
            if (t_seq_idx < (int)target_cas.size() &&
                aln_idx * 3 + 2 < (int)hit.alns.size()) {
                P_norm.push_back((hit.alns[aln_idx*3]   - norm_cx) * norm_scale);
                P_norm.push_back((hit.alns[aln_idx*3+1] - norm_cy) * norm_scale);
                P_norm.push_back((hit.alns[aln_idx*3+2] - norm_cz) * norm_scale);
                Q_norm.push_back(target_cas[t_seq_idx][0]);
                Q_norm.push_back(target_cas[t_seq_idx][1]);
                Q_norm.push_back(target_cas[t_seq_idx][2]);
            }
            aln_idx++;
            t_seq_idx++;
        }

        int N = (int)std::min(P_norm.size(), Q_norm.size()) / 3;
        if (N >= 3) {
            kabsch(P_norm, Q_norm, N, U, T);
            const float q_cen[3] = {norm_cx, norm_cy, norm_cz};
            const float t_cen[3] = {t_cx, t_cy, t_cz};
            for (int r = 0; r < 3; r++) {
                T_ang[r] = T[r] / norm_scale + q_cen[r];
                for (int c = 0; c < 3; c++) T_ang[r] -= U[r*3+c] * t_cen[c];
            }
            computed_transform = true;
        }
    } else if (hit.has_aln && !hit.is_alis_format) {
        // 17컬럼: qaln/taln 기반 Kabsch
        std::vector<std::array<float,3>> query_cas;
        for (const auto& [cid, chain] : data[0]->get_ca_atoms()) {
            for (const auto& atom : chain) {
                query_cas.push_back({(atom.x - norm_cx) * norm_scale,
                                      (atom.y - norm_cy) * norm_scale,
                                      (atom.z - norm_cz) * norm_scale});
            }
        }
        std::vector<std::array<float,3>> target_cas;
        for (const auto& [cid, chain] : target_protein->get_ca_atoms()) {
            for (const auto& atom : chain) {
                target_cas.push_back({(atom.x - t_cx) * norm_scale,
                                       (atom.y - t_cy) * norm_scale,
                                       (atom.z - t_cz) * norm_scale});
            }
        }

        std::vector<float> P_norm, Q_norm;
        int q_idx = hit.qstart - 1;
        int t_idx = hit.tstart - 1;
        const int q_size = (int)query_cas.size();
        const int t_size = (int)target_cas.size();

        for (size_t ai = 0; ai < hit.qaln.size() && ai < hit.taln.size(); ai++) {
            bool q_gap = (hit.qaln[ai] == '-');
            bool t_gap = (hit.taln[ai] == '-');
            if (!q_gap && !t_gap) {
                if (q_idx < q_size && t_idx < t_size) {
                    P_norm.push_back(query_cas[q_idx][0]);
                    P_norm.push_back(query_cas[q_idx][1]);
                    P_norm.push_back(query_cas[q_idx][2]);
                    Q_norm.push_back(target_cas[t_idx][0]);
                    Q_norm.push_back(target_cas[t_idx][1]);
                    Q_norm.push_back(target_cas[t_idx][2]);
                }
            }
            if (!q_gap) q_idx++;
            if (!t_gap) t_idx++;
        }

        int N = (int)std::min(P_norm.size(), Q_norm.size()) / 3;
        if (N >= 3) {
            kabsch(P_norm, Q_norm, N, U, T);
            const float q_cen[3] = {norm_cx, norm_cy, norm_cz};
            const float t_cen2[3] = {t_cx, t_cy, t_cz};
            for (int r = 0; r < 3; r++) {
                T_ang[r] = T[r] / norm_scale + q_cen[r];
                for (int c = 0; c < 3; c++) T_ang[r] -= U[r*3+c] * t_cen2[c];
            }
            computed_transform = true;
        }
    }

    // fallback: 전체 CA 순서 매칭
    if (!computed_transform) {
        std::vector<std::array<float,3>> query_cas;
        for (const auto& [cid, chain] : data[0]->get_ca_atoms()) {
            for (const auto& atom : chain) {
                query_cas.push_back({(atom.x - norm_cx) * norm_scale,
                                      (atom.y - norm_cy) * norm_scale,
                                      (atom.z - norm_cz) * norm_scale});
            }
        }
        std::vector<std::array<float,3>> target_cas;
        for (const auto& [cid, chain] : target_protein->get_ca_atoms()) {
            for (const auto& atom : chain) {
                target_cas.push_back({(atom.x - t_cx) * norm_scale,
                                       (atom.y - t_cy) * norm_scale,
                                       (atom.z - t_cz) * norm_scale});
            }
        }

        int N = (int)std::min(query_cas.size(), target_cas.size());
        if (N >= 3) {
            std::vector<float> P_norm, Q_norm;
            P_norm.reserve(N * 3);
            Q_norm.reserve(N * 3);
            for (int i = 0; i < N; i++) {
                P_norm.push_back(query_cas[i][0]);
                P_norm.push_back(query_cas[i][1]);
                P_norm.push_back(query_cas[i][2]);
                Q_norm.push_back(target_cas[i][0]);
                Q_norm.push_back(target_cas[i][1]);
                Q_norm.push_back(target_cas[i][2]);
            }
            kabsch(P_norm, Q_norm, N, U, T);
            const float q_cen[3] = {norm_cx, norm_cy, norm_cz};
            const float t_cen2[3] = {t_cx, t_cy, t_cz};
            for (int r = 0; r < 3; r++) {
                T_ang[r] = T[r] / norm_scale + q_cen[r];
                for (int c = 0; c < 3; c++) T_ang[r] -= U[r*3+c] * t_cen2[c];
            }
            computed_transform = true;
        }
    }

    if (computed_transform) {
        apply_foldseek_transform(target_protein_idx, U, T, T_ang);
    }

    // align 계열 모드일 때 is_aligned 계산
    if (is_aligned_mode(screen_mode)) {
        if (screen_mode == "align-near") {
            compute_aligned_all();
        } else if (hit.has_aln) {
            data[0]->compute_aligned_regions_from_aln(
                *data[target_protein_idx], hit.qaln, hit.taln,
                hit.qstart, hit.tstart, 5.0f, true);
            set_align_method("aln-string");
        } else if (screen_mode == "align-fs") {
            set_align_method("none");
        } else {
            compute_aligned_all();
        }
    }

    depth_calibrated = false;
}

// ── 기능 8: FoldMason MSA 기반 superposition ─────────────────────────────────

void Screen::set_foldmason(std::unique_ptr<FoldMasonParser> parser) {
    foldmason_parser = std::move(parser);
}

void Screen::set_foldmason_panel_info(const FoldMasonInfo& info) {
    if (panel) panel->set_foldmason_info(info);
}

void Screen::apply_foldmason_superposition(int query_protein_idx, int target_protein_idx,
                                           int fm_query_entry_idx, int fm_target_entry_idx) {
    if (!foldmason_parser) return;
    if (query_protein_idx < 0 || query_protein_idx >= (int)data.size() || !data[query_protein_idx]) return;
    if (target_protein_idx < 0 || target_protein_idx >= (int)data.size() || !data[target_protein_idx]) return;

    const auto& entries = foldmason_parser->get_entries();
    if (fm_query_entry_idx >= (int)entries.size() || fm_target_entry_idx >= (int)entries.size()) return;

    auto pairs = foldmason_parser->build_aligned_pairs(fm_query_entry_idx, fm_target_entry_idx);

    // query/target protein CA atoms 플랫화 (screen_atoms 순서)
    std::vector<std::array<float,3>> query_cas;
    for (auto& [cid, chain] : data[query_protein_idx]->get_atoms()) {
        for (const auto& a : chain) query_cas.push_back({a.x, a.y, a.z});
    }
    std::vector<std::array<float,3>> target_cas;
    for (auto& [cid, chain] : data[target_protein_idx]->get_atoms()) {
        for (const auto& a : chain) target_cas.push_back({a.x, a.y, a.z});
    }

    const int q_size = (int)query_cas.size();
    const int t_size = (int)target_cas.size();

    // P = query (참조), Q = target (회전 대상)
    std::vector<float> P_flat, Q_flat;

    if (!pairs.empty()) {
        // MSA aligned pairs 사용
        for (const auto& [ref_res, oth_res] : pairs) {
            if (ref_res >= q_size || oth_res >= t_size) continue;
            P_flat.push_back(query_cas[ref_res][0]);
            P_flat.push_back(query_cas[ref_res][1]);
            P_flat.push_back(query_cas[ref_res][2]);
            Q_flat.push_back(target_cas[oth_res][0]);
            Q_flat.push_back(target_cas[oth_res][1]);
            Q_flat.push_back(target_cas[oth_res][2]);
        }
    } else {
        // 작업 1-C: fallback — pairs가 비어있으면 전체 CA 순서 매칭
        int fallback_n = std::min(q_size, t_size);
        for (int i = 0; i < fallback_n; i++) {
            P_flat.push_back(query_cas[i][0]);
            P_flat.push_back(query_cas[i][1]);
            P_flat.push_back(query_cas[i][2]);
            Q_flat.push_back(target_cas[i][0]);
            Q_flat.push_back(target_cas[i][1]);
            Q_flat.push_back(target_cas[i][2]);
        }
    }

    int N = (int)std::min(P_flat.size(), Q_flat.size()) / 3;
    if (N < 3) return;

    float U[9], T[3];
    kabsch(P_flat, Q_flat, N, U, T);

    // BUG-A 2단계: kabsch는 정규화 screen_atoms 공간에서 계산됨.
    // init_atoms는 Å 공간이므로 T_Å = T_norm/norm_scale + q_centroid - U*t_centroid
    float T_ang[3];
    {
        const float q_cen[3] = {norm_cx, norm_cy, norm_cz};
        const float t_cen[3] = {data[target_protein_idx]->cx,
                                 data[target_protein_idx]->cy,
                                 data[target_protein_idx]->cz};
        for (int r = 0; r < 3; r++) {
            T_ang[r] = T[r] / norm_scale + q_cen[r];
            for (int c = 0; c < 3; c++) T_ang[r] -= U[r*3+c] * t_cen[c];
        }
    }
    apply_foldseek_transform(target_protein_idx, U, T, T_ang);

    // align 계열 모드일 때 is_aligned 잔기 설정
    // MSA aa strings을 qaln/taln으로 사용 (gap 형식 동일)
    if (screen_mode == "align-near") {
        compute_aligned_all();
    } else if (is_aligned_mode(screen_mode)) {
        // MSA aa 문자열은 서열 전체를 덮으므로 시작 오프셋은 1, 1
        data[query_protein_idx]->compute_aligned_regions_from_aln(
            *data[target_protein_idx],
            entries[fm_query_entry_idx].aa,
            entries[fm_target_entry_idx].aa,
            1, 1, 5.0f, true);
        set_align_method("msa-col");
    }

    depth_calibrated = false;
}