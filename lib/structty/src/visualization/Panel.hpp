#pragma once
#include "Palette.hpp"
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <cstring>

struct Entry {
    std::string file_name;
    std::map<std::string, int> chain_atom_info;
    std::map<std::string, int> chain_residue_info;
    // -1 이면 entry 순번을 쓴다. 멀티머는 complex 단위 색 + 전역 체인 순번을 넘긴다.
    int color_group = -1;
    int chain_color_base = -1;
};

// 기능 8: FoldMason MSA 정보 (패널 표시용)
struct FoldMasonInfo {
    bool        valid        = false;
    int         entry_count  = 0;
    std::string align_method;  // "msa-col" or "-"
};

// 기능 3: Foldseek hit 정보 (패널 표시용)
struct FoldseekHitInfo {
    bool     valid        = false;
    int      current_idx  = 0;   // 1-based
    int      total_hits   = 0;
    int      query_idx    = 0;   // 1-based (Step 5: multi-query nav); 0 = single query
    int      total_queries = 0;  // >1 enables query i/N display + ]/[ hint
    std::string target;
    float    evalue       = 0.0f;
    float    prob         = -1.0f;
    float    lddt         = -1.0f;
    float    qtmscore     = -1.0f;
    std::string align_method;   // "aln-string" or "nearest-nbr"
    std::string status_msg;     // 다운로드 상태 / 오류 메시지
    // 멀티머 `_report` 전용: e-value·정렬 문자열이 없는 대신 complex 지표를 보여준다
    bool     multimer     = false;
    float    ttmscore     = -1.0f;
    float    qcov         = -1.0f;
    float    tcov         = -1.0f;
    std::string interface_lddt;
    int      ass_id       = -1;
};

class Panel {
public:
    Panel(int width, const std::string& mode, bool show_structure = false);

    void add_panel_info(const std::string& file_name,
                        const std::map<std::string, int>& chain_info,
                        const std::map<std::string, int>& chain_residue_info,
                        int color_group = -1, int chain_color_base = -1);

    // Step 5: clear all protein entries (used when switching query in multi-query mode)
    void reset_entries() { entries.clear(); }

    // 기능 3: 두 번째 항목(target protein) entry 갱신
    void update_entry(int idx, const std::string& file_name,
                      const std::map<std::string, int>& chain_info,
                      const std::map<std::string, int>& chain_residue_info);

    int get_height() const;
    int get_height_for_width(int max_cols, int compact_level = 0) const;

    void draw_panel(int start_row, int start_col,
                    int max_rows, int max_cols,
                    int compact_level = 0) const;

    // 기능 4: 정렬 방식 표시 ("nearest-nbr" or "aln-string")
    void set_align_method(const std::string& method);

    // 기능 3: Foldseek hit info 섹션
    void set_foldseek_hit_info(const FoldseekHitInfo& info);
    void clear_foldseek_hit_info();
    int  get_foldseek_section_height() const;

    // 기능 6: 마우스 hover — Residue Info 섹션
    void set_hover_residue(const std::string& chainID,
                           const char* residue_name,
                           int residue_number,
                           char structure,
                           float bfactor,
                           float conservation_score);
    void clear_hover_residue();

    // Residue Info 섹션의 줄 수 (항상 고정)
    int get_residue_section_height() const;

    // draw_panel()에서 기록된 Residue Info 시작 행
    int get_last_hover_row() const;

    // hover 섹션만 부분 갱신
    void draw_hover_section(int hover_start_row, int max_cols) const;

    // 기능 8: FoldMason MSA 섹션
    void set_foldmason_info(const FoldMasonInfo& info);
    void clear_foldmason_info();
    int  get_foldmason_section_height() const;

private:
    std::vector<Entry> entries;

    int panel_width;
    std::string panel_mode;
    bool panel_show_structure = false;
    std::string align_method;

    // 기능 3: foldseek hit info 상태
    FoldseekHitInfo fs_hit_info;

    // 기능 8: foldmason info 상태
    FoldMasonInfo fm_info;

    // 기능 6: hover 상태
    mutable int last_hover_row = -1;  // draw_panel()에서 기록된 Residue Info 시작 행
    bool hover_valid = false;
    std::string hover_chain;
    char hover_residue_name[4] = {};
    int hover_residue_number = -1;
    char hover_structure = 0;
    float hover_bfactor = 0.0f;
    float hover_conservation = -1.0f;
};
