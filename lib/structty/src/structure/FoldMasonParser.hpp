#pragma once
#include <string>
#include <vector>
#include <array>
#include <utility>

struct FoldMasonEntry {
    std::string name;
    std::string aa;   // gap-포함 아미노산 정렬 (MSA 열 길이 동일)
    std::string ss;   // gap-포함 3Di 구조 알파벳 (옵션, JSON만)
    std::vector<std::array<float,3>> ca;  // Cα 좌표 (gap 제외 잔기 수와 동일, JSON만)
};

class FoldMasonParser {
public:
    bool load_json(const std::string& path);
    bool load_fasta(const std::string& path);

    const std::vector<FoldMasonEntry>& get_entries() const { return entries; }
    const std::vector<float>& get_scores() const { return scores; }
    float get_msa_lddt() const { return msa_lddt; }
    int msa_length() const;
    int entry_count() const { return (int)entries.size(); }

    // ref_idx 잔기(gap 제외) → MSA 열 인덱스 매핑
    std::vector<int> build_query_col_map(int ref_idx = 0) const;

    // 두 entry 간 gap-free 공통 열의 잔기 쌍 추출 (superposition용)
    // 반환: {ref_residue_idx, other_residue_idx} 쌍 목록 (0-based, gap 제외 순서)
    std::vector<std::pair<int,int>> build_aligned_pairs(int ref_idx, int other_idx) const;

    // 열별 Shannon entropy 기반 conservation 계산 (1=완전 보존, 0=최소 보존)
    // use_ss=true이면 ss 필드 사용, false이면 aa 사용
    std::vector<float> compute_column_entropy(bool use_ss = false) const;

private:
    std::vector<FoldMasonEntry> entries;
    std::vector<float> scores;
    float msa_lddt = -1.0f;
    std::string format;  // "json" 또는 "fasta"
};
