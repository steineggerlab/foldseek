#pragma once

#include <string>
#include <vector>

// 기능 5: MSA Conservation Score
// MMseqs2 easy-msa 결과(.m8) 또는 FASTA/A3M 형식의 MSA 파일을 파싱하고
// 각 쿼리 잔기 위치의 conservation score(Shannon entropy 기반)를 계산한다.
class MSAParser {
public:
    // filepath: FASTA/A3M 형식의 MSA 파일 경로 (.m8 포함)
    // 반환: 성공 시 true, 실패(파일 없음·형식 불명) 시 false
    bool load(const std::string& filepath);

    // 각 쿼리 위치의 conservation score를 계산하여 반환
    // 반환 벡터 길이 = 쿼리 서열 길이
    // 값 범위: [0.0, 1.0] (0=가변, 1=보존), -1.0 = 정보 없음(갭만 존재)
    std::vector<float> compute_conservation();

private:
    std::vector<std::string> sequence_names;
    std::vector<std::string> sequences;   // A3M 정규화 후 (소문자 삽입 제거)
    std::string query_sequence;           // sequences[0] (쿼리)
    std::vector<int> msa_pos_to_query_idx;
};
