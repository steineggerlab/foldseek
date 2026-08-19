#pragma once
#include "Protein.hpp"
#include "Atom.hpp"
#include "RenderPoint.hpp"
#include "Palette.hpp"
#include "Camera.hpp"
#include "Panel.hpp"
#include "Benchmark.hpp"
#include "Renderer.hpp"
#include "../structure/FoldseekParser.hpp"
#include "../structure/PDBDownloader.hpp"
#include "../structure/FoldMasonParser.hpp"
#include "../structure/FoldseekDBReader.hpp"
#include "../structure/MultimerReportParser.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>   // clamp, max
#include <limits>      // numeric_limits
#include <map>
#include <unordered_map>
#include <memory>
#include <string>

class Screen {
public:
    Screen(const int& width, const int& height, const bool& show_structure,
           const std::string& mode);
    ~Screen();

    // 기능 6: 마우스 hover — main loop용 (needs_redraw 반환)
    bool handle_input(bool& needs_redraw);
    // 벤치마크용: needs_redraw 없이 키 직접 처리
    bool handle_input(int key);


    void set_protein(const std::string& in_file, int ii, const bool& show_structure);

    // Foldseek m8 accession(`<파일 stem>_<체인>`)에서 체인을 도출해 idx 번 입력 파일의
    // 체인 필터로 쓴다. 사용자가 --chains 로 이미 지정한 값("-" 이 아니면)은 건드리지 않는다.
    // 반환: 실제로 적용한 체인 ID (적용하지 않았으면 빈 문자열)
    std::string apply_accession_chain(int idx, const std::string& accession,
                                     const std::string& file_path);

    // accession 이 "<파일 stem>_" 로 시작하면 그 뒤를 체인 ID 로 반환, 아니면 "-"
    static std::string chain_from_accession(const std::string& accession,
                                            const std::string& file_path);

    // Step 4 (D3/D4): load the query structure (data[0]) from a Foldseek query
    // tmp DB by accession, instead of parsing a plaintext CLI path. Returns true
    // on success. Used for folder/tar/gz query inputs.
    bool set_query_from_db(const std::string& query_db_path,
                           const std::string& accession,
                           const bool& show_structure);

    // Step 5 (D7/D8/D9): register multi-query navigation state and activate the
    // first query. query_ids is in .m8 order; hits_by_query groups hits per query.
    // Owns the full per-query setup (load query, normalize, target hits, first hit)
    // so query switching reuses one code path. target_db_path may be empty.
    void set_query_nav(const std::vector<std::string>& query_ids,
                       const std::map<std::string, std::vector<FoldseekHit>>& hits_by_query,
                       const std::string& query_db_path,
                       const std::string& target_db_path,
                       const bool& show_structure);
    // delta=+1 (next, ']'), delta=-1 (prev, '['). No-op if <2 queries.
    // 구조 파일 query 용 multi-chain 내비게이션: query_ids 는 이 파일의 체인 accession
    // (`<stem>_<chain>`) 이고, 전환 때마다 같은 파일을 다른 체인 필터로 다시 읽는다.
    // query_source 가 파일이면 그 파일을, 디렉터리면 accession 이름으로 그 안에서 찾은
    // 파일을 매 전환마다 다시 읽는다(체인 필터는 accession 에서 뽑는다).
    void set_query_nav_from_file(const std::vector<std::string>& query_ids,
                                 const std::map<std::string, std::vector<FoldseekHit>>& hits_by_query,
                                 const std::string& query_source,
                                 const std::string& target_db_path,
                                 const bool& show_structure,
                                 const bool& source_is_directory = false);
    void switch_query(int delta);
    int  query_count() const { return (int)query_ids_.size(); }

    // Step 7 (D6): multimer `_report` 경로. complex 단위로 query/target 전체 체인을
    // 로드하고 complex U/T 로 target 을 query frame 에 겹침. query complex 간 ]/[ 이동,
    // target complex hit 간 n/p 순회. query/target DB 는 complex DB(체인별 엔트리).
    void set_multimer_report(const std::vector<MultimerHit>& hits,
                             const std::string& query_source,
                             bool query_is_db,
                             bool query_is_dir,
                             const std::string& target_db_path,
                             const bool& show_structure);

    void normalize_proteins();

    void set_tmatrix();
    void set_chainfile(const std::string& chainfile, int filesize);
    void set_zoom_level(float zoom);

    // 기능 1: 로드된 모든 Protein에 대해 inter-chain interface를 계산
    void compute_interface_all(float threshold = 8.0f);

    // 기능 4: 로드된 모든 Protein 쌍에 대해 nearest-neighbor 기반 정렬 잔기를 계산
    void compute_aligned_all(float threshold = 10.0f);

    // 기능 5: 지정 protein에 conservation scores 적용
    void apply_msa_conservation(int protein_idx, const std::vector<float>& scores);

    // 기능 4: -fs 기반 — Foldseek hit의 U/T transform을 지정 protein에 적용
    // U_flat: row-major 3x3 (9 elements)
    // T_norm: 정규화 공간 T (screen_atoms do_shift용)
    // T_angstrom: Å 공간 T (init_atoms용). nullptr이면 init_atoms 미갱신
    void apply_foldseek_transform(int protein_idx, const float* U_flat, const float* T_norm,
                                  const float* T_angstrom = nullptr);

    // 기능 4: -fs 기반 — alignment string으로 aligned 잔기 계산 (protein0 vs protein1)
    // q_start/t_start: Foldseek qstart/tstart (1-based). 전체 서열 정렬은 1, 1.
    void compute_aligned_from_aln(const std::string& qaln, const std::string& taln,
                                  int q_start = 1, int t_start = 1,
                                  float threshold = 5.0f,
                                  bool skip_distance_check = false);

    // 기능 4: 패널에 정렬 방식 표시 설정 ("aln-string" or "nearest-nbr")
    void set_align_method(const std::string& method);

    // 기능 3: Foldseek hit 탐색
    void set_foldseek_hits(const std::vector<FoldseekHit>& hits);
    void set_fs_db_path(const std::string& path);
    void load_next_hit(int delta);  // delta=+1: next, delta=-1: prev, delta=0: first

    // Foldseek DB 직접 읽기 (v2: 선택적 스캔)
    bool open_foldseek_db(const std::string& db_base_path);
    bool prepare_foldseek_db(const std::vector<FoldseekHit>& hits);

    // 기능 3: 이미 로드된 target protein에 hit의 U/T transform 적용
    void apply_hit_transform(int target_protein_idx, const FoldseekHit& hit);

    // 기능 8: FoldMason MSA 기반 superposition + aligned region 설정
    void set_foldmason(std::unique_ptr<FoldMasonParser> parser);
    void apply_foldmason_superposition(int query_protein_idx, int target_protein_idx,
                                       int fm_query_entry_idx, int fm_target_entry_idx);
    void set_foldmason_panel_info(const FoldMasonInfo& info);

    void draw_screen(bool no_panel);

    // 기능 6: 마우스 hover — 현재 커서 위치의 잔기 정보를 패널에 반영
    void update_hover_info(int mx, int my);

    void set_benchmark(Benchmark* b) { bm = b; }
    
    void update_total_len_ca() {
        int64_t sum = 0;
        for (auto* p : data) {
            if (!p) continue;
            sum += (int64_t)p->get_length();  // Cα-only length라고 가정
        }
        total_len_ca = sum;
    }

private:
    int screen_width, screen_height;
    float aspect_ratio;
    bool screen_show_structure;
    bool yesUT = false;
    std::string screen_mode;
    int structNum = -1;

    // Braille sub-pixel rendering (logical pixels now owned by renderer_)
    Renderer renderer_;
    void print_screen_braille(int y_offset);

    // Chain id 인터닝: RenderPoint/RenderAtom 이 std::string 대신 int index 를 들고
    // 다니도록 하기 위한 string↔int 매핑 (단계 3 에서 사용). chain id 는 multi-char
    // 이므로 char 로 축약 불가 — int index + 역참조 테이블로 무손실 유지.
    std::vector<std::string> chain_names_;              // index → chain 문자열
    std::unordered_map<std::string, int> chain_intern_; // chain 문자열 → index
    int intern_chain(const std::string& id);            // 없으면 등록 후 index 반환
    const std::string& chain_name(int idx) const;       // index → 문자열 (범위 밖이면 빈 문자열)

    float focal_offset = 3.0f;
    float zoom_level;

    std::vector<float> pan_x;
    std::vector<float> pan_y;
    std::vector<std::string> chainVec;
    float** vectorpointer = nullptr;

    std::vector<Protein*> data;

    BoundingBox global_bb;
    Camera* camera = nullptr;
    Panel* panel = nullptr;

    bool depth_calibrated = false;
    float depth_base_min_z = 0.0f;
    float depth_base_max_z = 1.0f;

    Benchmark* bm = nullptr;
    bool ttff_logged = false;

    // 기능 6: 마우스 hover — 패널 위치 정보 (hover 갱신 시 필요)
    int last_panel_h = 0;
    int last_panel_start_row = 0;
    int last_panel_cols = 0;
    bool last_no_panel = false;

    // 기능 3: Foldseek hit 탐색
    std::vector<FoldseekHit> foldseek_hits;
    int current_hit_idx = -1;
    std::string fs_db_path;

    // Foldseek DB 직접 읽기
    FoldseekDBReader fs_db_reader_;

    // Step 4: query 구조를 query tmp DB 에서 읽기 (target 용 fs_db_reader_ 와 별도)
    FoldseekDBReader query_db_reader_;

    // Step 5: multi-query 내비게이션 상태
    std::vector<std::string> query_ids_;                            // .m8 순서
    std::map<std::string, std::vector<FoldseekHit>> hits_by_query_; // query별 hit 그룹
    int  current_query_idx_ = -1;
    std::string query_db_path_;
    // 비어있지 않으면 query 를 DB 가 아니라 이 구조 파일에서 읽는다(체인 필터는 accession 에서)
    std::string query_file_;
    // query_file_ 이 디렉터리인가 — accession 마다 그 안에서 파일을 찾는다
    bool query_file_is_dir_ = false;
    std::string target_db_path_;
    bool multi_query_show_structure_ = false;

    // 현재 query_ids_[idx] 를 활성화: scene teardown → query 로드 → 정규화 → target hit 셋업
    void activate_query(int idx);
    // query_db_reader_ 에서 accession 읽어 data[0] 로 push (set_query_from_db/activate_query 공유)
    bool load_query_into_data0(const std::string& accession, const bool& show_structure);

    // Step 7: multimer 상태 + 헬퍼
    bool multimer_mode_ = false;
    std::vector<std::string> mm_query_complexes_;                       // 등장 순서(유니크)
    std::map<std::string, std::vector<MultimerHit>> mm_hits_by_query_;  // query complex별 hit
    int mm_current_query_idx_ = -1;
    int mm_current_hit_idx_ = -1;
    int mm_query_chain_count_ = 0;   // data[0..count) = 현재 query complex 체인
    // query 소스: DB 면 query_db_reader_, 아니면 구조 파일 또는 그 디렉터리
    std::string mm_query_source_;
    bool mm_query_is_db_ = true;
    bool mm_query_is_dir_ = false;
    // 패널에 쓸 체인 이름(파일에서 읽으면 파일 경로가 아니라 accession 을 보여준다)
    std::vector<std::string> mm_chain_labels_;
    // query complex idx 활성화: teardown → query 체인 전체 로드 → 정규화 → 첫 target complex
    void activate_multimer_query(int idx);
    // 현재 query complex 의 target complex hit 간 순회(n/p) / query complex 간 이동(]/[)
    void load_multimer_hit(int delta);
    void switch_multimer_query(int delta);
    // 공유 centroid/scale 로 현재 data(=query 체인들) 정규화 (grid 없이 overlay 프레임)
    void normalize_complex();
    // complex DB 에서 accession 체인을 읽어 data 끝에 push. 실패 시 false.
    bool load_chain_into_data(FoldseekDBReader& reader, const std::string& accession,
                              const bool& show_structure);
    // complex 의 체인들을 소스 종류에 맞게 로드한다. DB 면 리더에서, 그 외에는 구조 파일에서
    // 읽는다(디렉터리·auto 는 complex 이름으로 한 번만 해석). 로드한 체인 수를 돌려준다.
    int load_complex_chains(const std::string& complex,
                            const std::vector<std::string>& chains,
                            bool from_db, FoldseekDBReader& reader,
                            const std::string& source_path, bool source_is_dir,
                            const bool& show_structure);
    // 이미 로드된 target 체인(data[idx])에 complex U/T(Å) 적용 + query frame 정규화
    void transform_target_chain(int idx, const float U[9], const float T[3]);
    // set_tmatrix 재호출 시 이전 vectorpointer 해제 (leak 방지)
    void free_tmatrix();
    // multimer mode에서 complex_idx(0=query, 1=target) → data[] 범위 {lo, hi} (inclusive)
    std::pair<int,int> mm_complex_range(int complex_idx) const;
    size_t vectorpointer_len_ = 0;

    // 기능 8: FoldMason MSA
    std::unique_ptr<FoldMasonParser> foldmason_parser;
    float norm_scale = 1.0f;   // normalize_proteins() 에서 저장
    float norm_cx = 0.0f;      // 정규화 시 실제로 shift 에 사용한 centroid
    float norm_cy = 0.0f;      //   = query(data[0]) 자기 centroid
    float norm_cz = 0.0f;
    float rot_pivot_[3] = {0.f, 0.f, 0.f}; // yesUT 회전 pivot: 정규화 후 씬 중심이 원점이므로 항상 {0,0,0}

    void calibrate_depth_baseline_first_view();

    std::vector<RenderAtom> to_render_atoms();
    void print_screen(int panel_lines);

    // 기능 6: handle_input 내부 구현 (두 public 오버로드 공유)
    bool handle_input_impl(int key, bool& needs_redraw);

    int64_t total_len_ca = 0;
};