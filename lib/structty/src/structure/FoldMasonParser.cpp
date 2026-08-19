#include "FoldMasonParser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cctype>

// ─── 간단한 JSON 파서 헬퍼 ───────────────────────────────────────────────────

static void skip_whitespace(const std::string& s, size_t& i) {
    while (i < s.size() && std::isspace((unsigned char)s[i])) i++;
}

// 현재 위치가 '"' 일 때 JSON 문자열을 읽어 반환 (탈출 시퀀스 미지원)
static std::string parse_json_string(const std::string& s, size_t& i) {
    if (i >= s.size() || s[i] != '"') return "";
    i++; // skip opening quote
    std::string result;
    while (i < s.size() && s[i] != '"') {
        if (s[i] == '\\') { i++; if (i < s.size()) result += s[i]; }
        else result += s[i];
        i++;
    }
    if (i < s.size()) i++; // skip closing quote
    return result;
}

// 현재 위치에서 JSON 숫자를 float으로 읽어 반환
static float parse_json_float(const std::string& s, size_t& i) {
    size_t start = i;
    if (i < s.size() && (s[i] == '-' || s[i] == '+')) i++;
    while (i < s.size() && (std::isdigit((unsigned char)s[i]) || s[i] == '.' ||
                             s[i] == 'e' || s[i] == 'E' || s[i] == '-' || s[i] == '+')) i++;
    if (i == start) return 0.0f;
    return std::stof(s.substr(start, i - start));
}

// JSON에서 "key" 를 찾아 이동 (현재 위치부터 탐색)
static bool find_json_key(const std::string& s, size_t& i, const std::string& key) {
    std::string needle = "\"" + key + "\"";
    size_t pos = s.find(needle, i);
    if (pos == std::string::npos) return false;
    i = pos + needle.size();
    skip_whitespace(s, i);
    if (i < s.size() && s[i] == ':') i++;
    skip_whitespace(s, i);
    return true;
}

// ─── FoldMasonParser::load_json ─────────────────────────────────────────────

bool FoldMasonParser::load_json(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "FoldMasonParser: cannot open " << path << std::endl;
        return false;
    }
    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());
    f.close();

    entries.clear();
    scores.clear();
    msa_lddt = -1.0f;
    format = "json";

    size_t i = 0;

    // entries 배열 파싱
    if (!find_json_key(content, i, "entries")) {
        std::cerr << "FoldMasonParser: 'entries' key not found in " << path << std::endl;
        return false;
    }
    skip_whitespace(content, i);
    if (i >= content.size() || content[i] != '[') return false;
    i++; // skip '['

    while (i < content.size()) {
        skip_whitespace(content, i);
        if (i >= content.size()) break;
        if (content[i] == ']') { i++; break; }
        if (content[i] != '{') { i++; continue; }

        FoldMasonEntry entry;
        size_t entry_end = i;
        // 중첩 중괄호 매칭으로 entry 블록 끝 찾기
        int depth = 0;
        for (size_t j = i; j < content.size(); j++) {
            if (content[j] == '{') depth++;
            else if (content[j] == '}') { depth--; if (depth == 0) { entry_end = j; break; } }
        }
        std::string entry_str = content.substr(i, entry_end - i + 1);

        // name
        size_t ei = 1;
        if (find_json_key(entry_str, ei, "name")) {
            entry.name = parse_json_string(entry_str, ei);
        }
        // aa
        ei = 0;
        if (find_json_key(entry_str, ei, "aa")) {
            entry.aa = parse_json_string(entry_str, ei);
        }
        // ss
        ei = 0;
        if (find_json_key(entry_str, ei, "ss")) {
            entry.ss = parse_json_string(entry_str, ei);
        }
        // ca — 배열 of 배열 [[x,y,z], ...]
        ei = 0;
        if (find_json_key(entry_str, ei, "ca")) {
            skip_whitespace(entry_str, ei);
            if (ei < entry_str.size() && entry_str[ei] == '[') {
                ei++; // outer [
                while (ei < entry_str.size()) {
                    skip_whitespace(entry_str, ei);
                    if (ei >= entry_str.size() || entry_str[ei] == ']') { ei++; break; }
                    if (entry_str[ei] == '[') {
                        ei++; // inner [
                        std::array<float,3> coord;
                        skip_whitespace(entry_str, ei);
                        coord[0] = parse_json_float(entry_str, ei);
                        skip_whitespace(entry_str, ei);
                        if (ei < entry_str.size() && entry_str[ei] == ',') ei++;
                        skip_whitespace(entry_str, ei);
                        coord[1] = parse_json_float(entry_str, ei);
                        skip_whitespace(entry_str, ei);
                        if (ei < entry_str.size() && entry_str[ei] == ',') ei++;
                        skip_whitespace(entry_str, ei);
                        coord[2] = parse_json_float(entry_str, ei);
                        skip_whitespace(entry_str, ei);
                        if (ei < entry_str.size() && entry_str[ei] == ']') ei++;
                        entry.ca.push_back(coord);
                    } else if (entry_str[ei] == ',') {
                        ei++;
                    } else {
                        ei++;
                    }
                }
            }
        }

        entries.push_back(std::move(entry));
        i = entry_end + 1;

        // skip comma
        skip_whitespace(content, i);
        if (i < content.size() && content[i] == ',') i++;
    }

    // scores 배열 파싱
    size_t si = 0;
    if (find_json_key(content, si, "scores")) {
        skip_whitespace(content, si);
        if (si < content.size() && content[si] == '[') {
            si++;
            while (si < content.size()) {
                skip_whitespace(content, si);
                if (si >= content.size() || content[si] == ']') { si++; break; }
                if (content[si] == ',') { si++; continue; }
                if (std::isdigit((unsigned char)content[si]) || content[si] == '-' || content[si] == '.') {
                    scores.push_back(parse_json_float(content, si));
                } else {
                    si++;
                }
            }
        }
    }

    // statistics.msaLDDT 파싱
    size_t stat_i = 0;
    if (find_json_key(content, stat_i, "statistics")) {
        if (find_json_key(content, stat_i, "msaLDDT")) {
            msa_lddt = parse_json_float(content, stat_i);
        }
    }

    return !entries.empty();
}

// ─── FoldMasonParser::load_fasta ─────────────────────────────────────────────

bool FoldMasonParser::load_fasta(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "FoldMasonParser: cannot open " << path << std::endl;
        return false;
    }

    entries.clear();
    scores.clear();
    msa_lddt = -1.0f;
    format = "fasta";

    std::string line;
    FoldMasonEntry cur;
    bool has_entry = false;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (has_entry) entries.push_back(std::move(cur));
            cur = FoldMasonEntry{};
            // name: trim '>' prefix and whitespace
            cur.name = line.substr(1);
            size_t p = cur.name.find_first_of(" \t");
            if (p != std::string::npos) cur.name = cur.name.substr(0, p);
            has_entry = true;
        } else {
            if (has_entry) cur.aa += line;
        }
    }
    if (has_entry) entries.push_back(std::move(cur));

    return !entries.empty();
}

// ─── msa_length ─────────────────────────────────────────────────────────────

int FoldMasonParser::msa_length() const {
    if (entries.empty()) return 0;
    return (int)entries[0].aa.size();
}

// ─── build_query_col_map ────────────────────────────────────────────────────

std::vector<int> FoldMasonParser::build_query_col_map(int ref_idx) const {
    if (ref_idx < 0 || ref_idx >= (int)entries.size()) return {};
    const std::string& aa = entries[ref_idx].aa;
    std::vector<int> col_map;
    for (int col = 0; col < (int)aa.size(); col++) {
        if (aa[col] != '-') col_map.push_back(col);
    }
    return col_map;
}

// ─── build_aligned_pairs ────────────────────────────────────────────────────

std::vector<std::pair<int,int>> FoldMasonParser::build_aligned_pairs(
        int ref_idx, int other_idx) const {
    if (ref_idx < 0 || ref_idx >= (int)entries.size()) return {};
    if (other_idx < 0 || other_idx >= (int)entries.size()) return {};
    const std::string& ref_aa = entries[ref_idx].aa;
    const std::string& oth_aa = entries[other_idx].aa;
    if (ref_aa.size() != oth_aa.size()) return {};

    std::vector<std::pair<int,int>> pairs;
    int ref_res = 0, oth_res = 0;
    for (int col = 0; col < (int)ref_aa.size(); col++) {
        bool ref_gap = (ref_aa[col] == '-');
        bool oth_gap = (oth_aa[col] == '-');
        if (!ref_gap && !oth_gap) {
            pairs.push_back({ref_res, oth_res});
        }
        if (!ref_gap) ref_res++;
        if (!oth_gap) oth_res++;
    }
    return pairs;
}

// ─── compute_column_entropy ─────────────────────────────────────────────────

std::vector<float> FoldMasonParser::compute_column_entropy(bool use_ss) const {
    if (entries.empty()) return {};
    int ncols = msa_length();
    std::vector<float> result(ncols, 0.0f);

    for (int col = 0; col < ncols; col++) {
        int counts[128] = {};
        int total = 0;
        for (const auto& e : entries) {
            const std::string& seq = use_ss ? e.ss : e.aa;
            if (col >= (int)seq.size()) continue;
            char c = seq[col];
            if (c == '-' || c == '\0') continue;
            counts[(unsigned char)c]++;
            total++;
        }
        if (total == 0) { result[col] = 0.0f; continue; }

        float H = 0.0f;
        for (int ci = 0; ci < 128; ci++) {
            if (counts[ci] == 0) continue;
            float p = (float)counts[ci] / (float)total;
            H -= p * std::log2(p);
        }
        // conservation = 1 - H / log2(20)
        static const float log2_20 = std::log2(20.0f);
        float conservation = 1.0f - H / log2_20;
        if (conservation < 0.0f) conservation = 0.0f;
        if (conservation > 1.0f) conservation = 1.0f;
        result[col] = conservation;
    }
    return result;
}
