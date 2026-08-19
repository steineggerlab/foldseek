#include "Panel.hpp"
#include "Terminal.hpp"
#include "Common.hpp"
#include <cstring>  // strncpy
#include <cstdio>   // snprintf, fwrite, fputc, fputs

Panel::Panel(int width, const std::string& mode, bool show_structure)
    : panel_width(width), panel_mode(mode), panel_show_structure(show_structure) {}

void Panel::add_panel_info(const std::string& file_name, 
                           const std::map<std::string, int>& chain_info, 
                           const std::map<std::string, int>& chain_residue_info,
                           int color_group, int chain_color_base) {
    entries.push_back(Entry{
        file_name,
        chain_info,
        chain_residue_info,
        color_group,
        chain_color_base
    });
}

void Panel::set_align_method(const std::string& method) {
    align_method = method;
}

// 기능 3: entry 갱신 (target protein 교체 시)
void Panel::update_entry(int idx, const std::string& file_name,
                          const std::map<std::string, int>& chain_info,
                          const std::map<std::string, int>& chain_residue_info) {
    if (idx < 0) return;
    if (idx < (int)entries.size()) {
        entries[idx] = Entry{file_name, chain_info, chain_residue_info};
    } else {
        while ((int)entries.size() < idx) {
            entries.push_back(Entry{});
        }
        entries.push_back(Entry{file_name, chain_info, chain_residue_info});
    }
}

// 기능 3: Foldseek hit info 섹션
void Panel::set_foldseek_hit_info(const FoldseekHitInfo& info) {
    fs_hit_info = info;
}

void Panel::clear_foldseek_hit_info() {
    fs_hit_info = FoldseekHitInfo{};
}

int Panel::get_foldseek_section_height() const {
    if (!fs_hit_info.valid && fs_hit_info.total_hits == 0) return 0;
    return fs_hit_info.multimer ? 6 : 8;
}

// 기능 8: FoldMason MSA 섹션 ------------------------------------------------

void Panel::set_foldmason_info(const FoldMasonInfo& info) {
    fm_info = info;
}

void Panel::clear_foldmason_info() {
    fm_info = FoldMasonInfo{};
}

int Panel::get_foldmason_section_height() const {
    if (!fm_info.valid) return 0;
    return 3;  // "FoldMason MSA" + "Entries: N" + "Align: X"
}

// 기능 6: Residue Info hover -----------------------------------------------

void Panel::set_hover_residue(const std::string& chainID,
                               const char* residue_name,
                               int residue_number,
                               char structure,
                               float bfactor,
                               float conservation_score) {
    hover_valid          = true;
    hover_chain          = chainID;
    strncpy(hover_residue_name, residue_name, 3);
    hover_residue_name[3] = '\0';
    hover_residue_number  = residue_number;
    hover_structure       = structure;
    hover_bfactor         = bfactor;
    hover_conservation    = conservation_score;
}

void Panel::clear_hover_residue() {
    hover_valid = false;
}

int Panel::get_last_hover_row() const {
    return last_hover_row;
}

int Panel::get_residue_section_height() const {
    // 항상 고정: header + chain + res + ss (= 4 기본)
    // + pLDDT 줄 (plddt 모드일 때) + Cons 줄 (conservation 모드일 때)
    int h = 4;
    if (panel_mode == "plddt")        h += 1;
    if (panel_mode == "conservation") h += 1;
    return h;
}

void Panel::draw_hover_section(int hover_start_row, int max_cols) const {
    // hover_start_row: "Residue Info" 헤더 행 (separator는 이미 draw_panel이 그렸으므로 제외)
    // 이 함수는 Residue Info 섹션만 다시 그린다 (bottom border 제외)
    int r = hover_start_row;
    int right_limit = max_cols - 1;

    auto clear_ln = [&](int rr) {
        Terminal::move_cursor(rr, 0);
        Terminal::clear_to_eol();
        Terminal::move_cursor(rr, 0);
    };

    auto put_text = [&](int rr, const char* s) {
        (void)rr;
        int len = (int)std::strlen(s);
        int k   = std::min(len, right_limit);
        if (k > 0) fwrite(s, 1, (size_t)k, stdout);
    };

    // "Residue Info" 헤더
    clear_ln(r);
    put_text(r, "Residue Info");
    ++r;

    // Chain: X  or  Chain: -
    clear_ln(r);
    if (hover_valid) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "Chain: %s", hover_chain.empty() ? "-" : hover_chain.c_str());
        put_text(r, buf);
    } else {
        put_text(r, "Chain: -");
    }
    ++r;

    // Res: GLU 42  or  Res: -
    clear_ln(r);
    if (hover_valid && hover_residue_number >= 0) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "Res:   %s %d",
                      (hover_residue_name[0] ? hover_residue_name : "?"),
                      hover_residue_number);
        put_text(r, buf);
    } else {
        put_text(r, "Res:   -");
    }
    ++r;

    // SS: Helix / Sheet / Coil  or  SS: -
    clear_ln(r);
    if (hover_valid) {
        const char* ss_str = "Coil";
        if      (hover_structure == 'H') ss_str = "Helix";
        else if (hover_structure == 'S') ss_str = "Sheet";
        char buf[32];
        std::snprintf(buf, sizeof(buf), "SS:    %s", ss_str);
        put_text(r, buf);
    } else {
        put_text(r, "SS:    -");
    }
    ++r;

    // pLDDT: 87.3  (plddt 모드일 때만)
    if (panel_mode == "plddt") {
        clear_ln(r);
        if (hover_valid) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "pLDDT: %.1f", hover_bfactor);
            put_text(r, buf);
        } else {
            put_text(r, "pLDDT: -");
        }
        ++r;
    }

    // Cons: 0.82  (conservation 모드일 때만)
    if (panel_mode == "conservation") {
        clear_ln(r);
        if (hover_valid && hover_conservation >= 0.0f) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "Cons:  %.2f", hover_conservation);
            put_text(r, buf);
        } else {
            put_text(r, "Cons:  -");
        }
        ++r;
    }
}

// 기능 6 끝 -----------------------------------------------------------------

int Panel::get_height() const {
    int lines = 0;
    lines += 3;
    if (is_aligned_mode(panel_mode) && !align_method.empty()) {
        lines += 1;  // "Align: nearest-nbr" or "Align: aln-string"
    }
    for (const auto& entry : entries) {
        lines += 1;
        int n = (int)entry.chain_atom_info.size();
        int chain_lines = (n == 0) ? 1 : ((n + 2) / 3); // 3 per line
        lines += chain_lines;
        lines += 1;
    }
    // 기능 3: Foldseek hit info 섹션
    int fs_h = get_foldseek_section_height();
    if (fs_h > 0) {
        lines += 1;   // separator
        lines += fs_h;
    }
    // 기능 8: FoldMason MSA 섹션
    int fm_h = get_foldmason_section_height();
    if (fm_h > 0) {
        lines += 1;   // separator
        lines += fm_h;
    }
    // 기능 6: separator + Residue Info 섹션
    lines += 1;                             // separator (---)
    lines += get_residue_section_height();  // 고정 줄 수
    lines += 1;                             // bottom border
    return lines;
}
void Panel::draw_panel(int start_row, int start_col,
                       int max_rows, int max_cols,
                       int compact_level) const {
    const int num_protein_colors = 9;
    const int num_chain_colors   = 15;
    if (max_rows <= 0 || max_cols <= 0) return;

    const int top    = start_row;
    const int left   = start_col;
    const int bottom = start_row + max_rows; // exclusive
    const int right  = start_col + max_cols; // exclusive

    const int right_limit = right - 1;

    auto in_rows = [&](int rr){ return rr >= top && rr < bottom; };
    auto remain_cols = [&](int x){ return right_limit - x; };

    auto clear_line = [&](int rr){
        if (!in_rows(rr)) return;
        Terminal::move_cursor(rr, left);
        Terminal::clear_to_eol();
        Terminal::move_cursor(rr, left);
    };

    auto put_n = [&](int& rr, int& x, const char* s, int n){
        (void)rr;
        if (!in_rows(rr)) return;
        int rem = remain_cols(x);
        if (rem <= 0 || n <= 0) return;
        int k = std::min(rem, n);
        fwrite(s, 1, (size_t)k, stdout);
        x += k;
    };

    auto put_str = [&](int& rr, int& x, const std::string& s){
        put_n(rr, x, s.c_str(), (int)s.size());
    };

    auto put_cstr = [&](int& rr, int& x, const char* s){
        put_n(rr, x, s, (int)std::strlen(s));
    };

    auto put_indent = [&](int& rr, int& x){
        put_str(rr, x, "  ");
    };

    int r = start_row;

    // Top border
    clear_line(r);
    {
        int x = left;
        put_str(r, x, "*");
        int mid = std::max(0, std::min(panel_width - 2, max_cols - 2));
        put_str(r, x, std::string(mid, '='));
        if (max_cols >= 2) put_str(r, x, "*");
    }
    ++r;
    if (!in_rows(r)) return;

    // Help line
    clear_line(r);
    {
        int x = left;
        put_str(r, x, "W A S D : ^ < v >  ");
        put_str(r, x, "R F : Zoom In/Out  ");
        put_str(r, x, "X Y Z : Rotate X, Y, Z axis  ");
        put_str(r, x, "C : Screenshot  ");
        put_str(r, x, "Q : Quit");
    }
    ++r;
    if (!in_rows(r)) return;

    // Separator
    clear_line(r);
    {
        Terminal::move_cursor(r, left);
        int w = std::min(panel_width, max_cols);
        w = std::min(w, max_cols - 1);
        for (int i = 0; i < w; ++i) fputc('-', stdout);
    }
    ++r;
    if (!in_rows(r)) return;

    // 기능 4: align 계열 모드일 때 정렬 방식 표시
    if (is_aligned_mode(panel_mode) && !align_method.empty()) {
        if (!in_rows(r)) return;
        clear_line(r);
        int x = left;
        put_cstr(r, x, "Align: ");
        put_str(r, x, align_method);
        ++r;
        if (!in_rows(r)) return;
    }

    // Body
    int file_idx = 0;
    for (const auto& entry : entries) {
        if (!in_rows(r)) break;

        const std::string& file_name = entry.file_name;
        const auto& chain_info       = entry.chain_atom_info;

        const int color_idx = (entry.color_group >= 0) ? entry.color_group : file_idx;

        int protein_pair = 0;
        if (panel_mode == "protein") {
            protein_pair = (color_idx % num_protein_colors) + 1;  // pairs 1-9
        } else if (is_aligned_mode(panel_mode)) {
            protein_pair = (color_idx % num_protein_colors) + 101;  // pairs 101-109
        } else if (panel_mode == "chain" && entry.chain_color_base >= 0) {
            // entry 하나가 체인 하나인 경우(멀티머)에는 이름 줄에 그 체인 색을 쓴다
            protein_pair = 21 + (entry.chain_color_base % num_chain_colors);
        }

        // file name line
        clear_line(r);
        {
            int x = left;
            if (protein_pair > 0) fputs(Palettes::palette_to_ansi_fg_str(protein_pair).c_str(), stdout);
            put_str(r, x, file_name);
            if (protein_pair > 0) fputs("\033[0m", stdout);
        }
        ++r;
        if (!in_rows(r)) break;

        // chain lines (compact_level < 3)
        if (compact_level <= 1) {
            // Level 0, 1: 기존 chain 상세 표시
            clear_line(r);
            Terminal::move_cursor(r, left);
            int x = left;
            put_indent(r, x);

            int count = 0;
            for (const auto& [chainID, length] : chain_info) {
                if (!in_rows(r)) break;

                int residue_cnt = 0;
                auto itC = entry.chain_residue_info.find(chainID);
                if (itC != entry.chain_residue_info.end()) residue_cnt = itC->second;

                char buf[64];
                int token_len = std::snprintf(buf, sizeof(buf), "%s: %d (%d)  ",
                                              chainID.c_str(), residue_cnt, length);
                if (token_len < 0) token_len = 0;
                if (token_len >= (int)sizeof(buf)) token_len = (int)sizeof(buf) - 1;

                if (x + token_len > right_limit) {
                    ++r;
                    if (!in_rows(r)) break;
                    clear_line(r);
                    Terminal::move_cursor(r, left);
                    x = left;
                    put_indent(r, x);
                }

                int chain_pair = 0;
                if (panel_mode == "chain") {
                    const int base = (entry.chain_color_base >= 0)
                                     ? entry.chain_color_base : file_idx * 10;
                    chain_pair = 21 + ((base + count) % num_chain_colors);  // pairs 21-35
                }

                int pair_to_use = (panel_mode == "protein" || is_aligned_mode(panel_mode)) ? protein_pair : chain_pair;

                if (pair_to_use > 0) fputs(Palettes::palette_to_ansi_fg_str(pair_to_use).c_str(), stdout);
                put_n(r, x, buf, token_len);
                if (pair_to_use > 0) fputs("\033[0m", stdout);

                ++count;
            }
            ++r;
        } else if (compact_level == 2) {
            // Level 2: "N chains" 한 줄 요약
            clear_line(r);
            {
                int x = left;
                put_indent(r, x);
                char buf[32];
                std::snprintf(buf, sizeof(buf), "%d chains", (int)chain_info.size());
                put_cstr(r, x, buf);
            }
            ++r;
        }
        // Level 3: chain 정보 생략 (파일명만 표시)

        // blank lines
        if (compact_level == 0) {
            // Level 0: blank line 2줄
            if (!in_rows(r)) break;
            clear_line(r);
            ++r;
        } else if (compact_level == 1) {
            // Level 1: blank line 없음 (++r은 chain 블록에서 이미 수행)
        }
        // Level 2, 3: blank line 없음

        ++file_idx;
    }

    if (!in_rows(r)) return;

    // 기능 3: Foldseek hit info 섹션
    {
        int fs_h = get_foldseek_section_height();
        if (fs_h > 0) {
            // separator
            clear_line(r);
            {
                Terminal::move_cursor(r, left);
                int w = std::min(panel_width, max_cols);
                w = std::min(w, max_cols - 1);
                for (int i = 0; i < w; ++i) fputc('-', stdout);
            }
            ++r;
            if (!in_rows(r)) return;

            const FoldseekHitInfo& fi = fs_hit_info;

            // Line 1: "Foldseek Hits" (+ "Q[i/N]" when multiple queries)
            clear_line(r);
            {
                int x = left;
                put_cstr(r, x, "Foldseek Hits");
                if (fi.total_queries > 1) {
                    char qbuf[32];
                    std::snprintf(qbuf, sizeof(qbuf), " Q[%d/%d]",
                                  fi.query_idx, fi.total_queries);
                    put_cstr(r, x, qbuf);
                }
            }
            ++r; if (!in_rows(r)) return;

            // Line 2: "[X / Y]" or hit count
            clear_line(r);
            {
                int x = left;
                char buf[32];
                std::snprintf(buf, sizeof(buf), "[%d / %d]", fi.current_idx, fi.total_hits);
                put_cstr(r, x, buf);
            }
            ++r; if (!in_rows(r)) return;

            // Line 3: "Target: ..."
            clear_line(r);
            {
                int x = left;
                put_cstr(r, x, "Target: ");
                if (!fi.target.empty()) put_str(r, x, fi.target);
                else put_cstr(r, x, "-");
            }
            ++r; if (!in_rows(r)) return;

            if (fi.multimer) {
                clear_line(r);
                {
                    int x = left;
                    char buf[80];
                    std::snprintf(buf, sizeof(buf), "TM  q %.3f / t %.3f   Cov q %.3f / t %.3f",
                                  fi.qtmscore, fi.ttmscore, fi.qcov, fi.tcov);
                    put_cstr(r, x, buf);
                }
                ++r; if (!in_rows(r)) return;

                clear_line(r);
                {
                    int x = left;
                    put_cstr(r, x, "iLDDT ");
                    put_str(r, x, fi.interface_lddt.empty() ? "-" : fi.interface_lddt);
                    char buf[32];
                    std::snprintf(buf, sizeof(buf), "   assId %d", fi.ass_id);
                    put_cstr(r, x, buf);
                }
                ++r; if (!in_rows(r)) return;
            } else {
            // Line 4: "E-val:  ..."
            clear_line(r);
            {
                int x = left;
                char buf[32];
                std::snprintf(buf, sizeof(buf), "E-val:  %.2e", fi.evalue);
                put_cstr(r, x, buf);
            }
            ++r; if (!in_rows(r)) return;

            // Line 5: prob or status_msg
            clear_line(r);
            {
                int x = left;
                if (!fi.status_msg.empty()) {
                    put_str(r, x, fi.status_msg);
                } else if (fi.prob >= 0.0f) {
                    char buf[32];
                    std::snprintf(buf, sizeof(buf), "Prob:   %.3f", fi.prob);
                    put_cstr(r, x, buf);
                } else {
                    put_cstr(r, x, "-");
                }
            }
            ++r; if (!in_rows(r)) return;

            // Line 6: lDDT and TM scores
            clear_line(r);
            {
                int x = left;
                char buf[80];
                int n = 0;
                if (fi.lddt >= 0.0f) {
                    n = std::snprintf(buf, sizeof(buf), "lDDT:   %.3f", fi.lddt);
                }
                if (fi.qtmscore >= 0.0f && n < (int)sizeof(buf)) {
                    if (fi.ttmscore >= 0.0f) {
                        std::snprintf(buf + n, sizeof(buf) - n, "%sTM  q %.3f / t %.3f",
                                      n > 0 ? "   " : "", fi.qtmscore, fi.ttmscore);
                    } else {
                        std::snprintf(buf + n, sizeof(buf) - n, "%sTM  q %.3f",
                                      n > 0 ? "   " : "", fi.qtmscore);
                    }
                    n = (int)std::strlen(buf);
                }
                if (n > 0) {
                    put_cstr(r, x, buf);
                } else {
                    put_cstr(r, x, "-");
                }
            }
            ++r; if (!in_rows(r)) return;

            // Line 7: superposition method (how U/T was obtained). The colour
            // source is the "Align:" line above and can differ from this one.
            clear_line(r);
            {
                int x = left;
                if (!fi.align_method.empty()) {
                    put_cstr(r, x, "Superpose: ");
                    put_str(r, x, fi.align_method);
                } else {
                    put_cstr(r, x, "Superpose: -");
                }
            }
            ++r; if (!in_rows(r)) return;
            }

            // Line 8: nav hint ( ]/[ query hint when multiple queries )
            clear_line(r);
            {
                int x = left;
                if (fi.total_queries > 1) {
                    put_cstr(r, x, "[N]/[P] hit  ]/[ qry");
                } else {
                    put_cstr(r, x, "[N]ext  [P]rev");
                }
            }
            ++r; if (!in_rows(r)) return;
        }
    }

    // 기능 8: FoldMason MSA 섹션
    {
        int fm_h = get_foldmason_section_height();
        if (fm_h > 0) {
            // separator
            clear_line(r);
            {
                Terminal::move_cursor(r, left);
                int w = std::min(panel_width, max_cols);
                w = std::min(w, max_cols - 1);
                for (int i = 0; i < w; ++i) fputc('-', stdout);
            }
            ++r; if (!in_rows(r)) return;

            const FoldMasonInfo& fi = fm_info;

            // Line 1: "FoldMason MSA"
            clear_line(r);
            { int x = left; put_cstr(r, x, "FoldMason MSA"); }
            ++r; if (!in_rows(r)) return;

            // Line 2: "Entries: N"
            clear_line(r);
            {
                int x = left;
                char buf[32];
                std::snprintf(buf, sizeof(buf), "Entries: %d", fi.entry_count);
                put_cstr(r, x, buf);
            }
            ++r; if (!in_rows(r)) return;

            // Line 3: "Align: msa-col" or "Align: -"
            clear_line(r);
            {
                int x = left;
                put_cstr(r, x, "Align: ");
                if (!fi.align_method.empty()) {
                    put_str(r, x, fi.align_method);
                } else {
                    put_cstr(r, x, "-");
                }
            }
            ++r; if (!in_rows(r)) return;
        }
    }

    // 기능 6: Residue Info 섹션 separator
    clear_line(r);
    {
        Terminal::move_cursor(r, left);
        int w = std::min(panel_width, max_cols);
        w = std::min(w, max_cols - 1);
        for (int i = 0; i < w; ++i) fputc('-', stdout);
    }
    ++r;
    if (!in_rows(r)) return;

    // 기능 6: Residue Info 섹션 (draw_hover_section과 동일한 내용)
    // draw_panel()은 전체 패널을 그리므로 여기서도 Residue Info를 그린다.
    // hover가 없을 때는 모두 "-" 표시.
    last_hover_row = r;  // hover 부분 갱신 시 정확한 행 사용
    draw_hover_section(r, max_cols);
    r += get_residue_section_height();

    if (!in_rows(r)) return;

    // Bottom border
    clear_line(r);
    {
        int x = left;
        put_str(r, x, "*");
        int mid = std::max(0, std::min(panel_width - 2, max_cols - 2));
        put_str(r, x, std::string(mid, '='));
        if (max_cols >= 2) put_str(r, x, "*");
    }
}

static int count_wrapped_lines_for_chaininfo(
    const std::map<std::string,int>& chain_info,
    const std::map<std::string,int>& chain_residue_info,
    int avail_cols,
    int indent_cols = 2
){
    if (avail_cols <= indent_cols) return 1;

    int lines = 1;
    int x = indent_cols;

    for (const auto& [chainID, length] : chain_info) {
        int residue_cnt = 0;
        auto itC = chain_residue_info.find(chainID);
        if (itC != chain_residue_info.end()) residue_cnt = itC->second;

        char buf[128];
        std::snprintf(buf, sizeof(buf), "%s: %d (%d)  ", chainID.c_str(), residue_cnt, length);

        int token_len = (int)std::strlen(buf);

        if (token_len > avail_cols - indent_cols) {
            if (x > indent_cols) { lines++; x = indent_cols; }
            x = indent_cols + std::min(token_len, avail_cols - indent_cols);
            continue;
        }

        if (x + token_len > avail_cols) {
            lines++;
            x = indent_cols;
        }
        x += token_len;
    }
    return lines;
}
int Panel::get_height_for_width(int max_cols, int compact_level) const {
    int lines = 0;

    lines += 3; // Top border + Help line + Separator

    if (is_aligned_mode(panel_mode) && !align_method.empty()) {
        lines += 1;  // "Align: ..." line
    }

    int avail_cols = max_cols;
    if (avail_cols < 1) avail_cols = 1;

    for (const auto& entry : entries) {
        lines += 1; // file name line

        if (compact_level <= 1) {
            // Level 0, 1: chain 정보 표시
            lines += count_wrapped_lines_for_chaininfo(
                entry.chain_atom_info, entry.chain_residue_info,
                /*avail_cols=*/avail_cols,
                /*indent_cols=*/2
            );
        } else if (compact_level == 2) {
            // Level 2: "N chains" 한 줄 요약
            lines += 1;
        }
        // Level 3: 파일명만 (chain 정보 없음)

        // blank lines
        if (compact_level == 0) {
            lines += 2;
        } else if (compact_level == 1) {
            lines += 1;
        }
        // Level 2, 3: blank line 없음
    }

    // 기능 3: Foldseek hit info 섹션
    int fs_h = get_foldseek_section_height();
    if (fs_h > 0) {
        lines += 1;   // separator
        lines += fs_h;
    }
    // 기능 8: FoldMason MSA 섹션
    int fm_h = get_foldmason_section_height();
    if (fm_h > 0) {
        lines += 1;   // separator
        lines += fm_h;
    }
    // 기능 6: separator + Residue Info 섹션 + bottom border
    lines += 1;                             // separator (---)
    lines += get_residue_section_height();  // 고정 줄 수
    lines += 1;                             // Bottom border
    return lines;
}

