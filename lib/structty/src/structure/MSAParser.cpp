#include "MSAParser.hpp"

#include <fstream>
#include <cmath>
#include <array>
#include <algorithm>

// A3M 정규화: 소문자(삽입, insertion) 제거 → 쿼리 기준 정렬 길이만 남긴다.
// 대문자 알파벳과 갭('-')만 유지한다.
static std::string normalize_a3m(const std::string& seq) {
    std::string result;
    result.reserve(seq.size());
    for (char c : seq) {
        if (c == '-' || (c >= 'A' && c <= 'Z')) {
            result += c;
        }
        // 소문자(삽입) → 제거
    }
    return result;
}

bool MSAParser::load(const std::string& filepath) {
    std::ifstream ifs(filepath);
    if (!ifs.is_open()) {
        return false;
    }

    sequence_names.clear();
    sequences.clear();
    query_sequence.clear();
    msa_pos_to_query_idx.clear();

    std::string line;

    // 첫 번째 비빈 줄로 형식 판별
    while (std::getline(ifs, line)) {
        if (!line.empty()) break;
    }

    if (line.empty()) return false;

    // Stockholm 형식: 미지원
    if (line.find("# STOCKHOLM") == 0) {
        return false;
    }

    // FASTA/A3M: 첫 문자가 '>'
    if (line[0] != '>') {
        return false;
    }

    // FASTA/A3M 파싱
    std::string cur_name = line.substr(1);
    std::string cur_seq;

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!cur_seq.empty()) {
                sequence_names.push_back(cur_name);
                sequences.push_back(normalize_a3m(cur_seq));
            }
            cur_name = line.substr(1);
            cur_seq.clear();
        } else {
            cur_seq += line;
        }
    }

    // 마지막 시퀀스 저장
    if (!cur_seq.empty()) {
        sequence_names.push_back(cur_name);
        sequences.push_back(normalize_a3m(cur_seq));
    }

    if (sequences.empty()) return false;

    query_sequence = sequences[0];
    return true;
}

std::vector<float> MSAParser::compute_conservation() {
    if (sequences.empty()) return {};

    int query_len = (int)query_sequence.size();
    if (query_len == 0) return {};

    std::vector<float> conservation(query_len, -1.0f);

    const float log2_20 = std::log2(20.0f);

    for (int pos = 0; pos < query_len; ++pos) {
        // 각 위치에서 아미노산 빈도 카운트 (갭 제외)
        std::array<int, 26> counts{};
        int total = 0;

        for (const std::string& seq : sequences) {
            if (pos >= (int)seq.size()) continue;
            char c = seq[pos];
            if (c >= 'A' && c <= 'Z') {
                counts[static_cast<int>(c - 'A')]++;
                total++;
            }
            // '-' (갭) 은 빈도 계산에서 제외
        }

        if (total == 0) {
            conservation[pos] = -1.0f;
            continue;
        }

        // Shannon entropy: H(i) = -Σ_aa [f(aa,i) * log2(f(aa,i))]
        float H = 0.0f;
        for (int aa = 0; aa < 26; ++aa) {
            if (counts[aa] == 0) continue;
            float f = static_cast<float>(counts[aa]) / static_cast<float>(total);
            H -= f * std::log2(f);
        }

        // conservation(i) = 1 - H(i) / log2(20)
        float score = 1.0f - H / log2_20;
        // 수치 오차 클램핑
        conservation[pos] = std::max(0.0f, std::min(1.0f, score));
    }

    return conservation;
}
