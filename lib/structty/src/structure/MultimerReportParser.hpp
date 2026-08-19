#pragma once

#include <string>
#include <vector>

// D6: Foldseek multimer 결과 `_report` 파일(14컬럼 tsv) 파서.
// `.m8` 와 달리 complex(assembly) 단위 1행 = (query complex, target complex) 매칭이며
// u/t 는 complex 전체 superposition (target → query frame). 컬럼 순서(spec §4-4):
//   1 qComplexName 2 tComplexName 3 qChains 4 tChains 5 qTMScore 6 tTMScore
//   7 u(회전 9개,쉼표) 8 t(이동 3개,쉼표) 9 qComplexCov 10 tComplexCov
//   11 qChainTms 12 tChainTms 13 interfaceLddtScore 14 assId

struct MultimerHit {
    std::string qComplex;                 // query complex 이름 (DB lookup base)
    std::string tComplex;                 // target complex 이름
    std::vector<std::string> qChains;     // query 체인 ID 목록 (쉼표 분리)
    std::vector<std::string> tChains;     // target 체인 ID 목록
    float qTMScore = -1.0f;
    float tTMScore = -1.0f;

    // complex 단위 superposition: target_query_frame = U * target + T  (row-major)
    float U[9] = {1,0,0, 0,1,0, 0,0,1};
    float T[3] = {0,0,0};
    bool has_transform = false;

    float qComplexCov = -1.0f;
    float tComplexCov = -1.0f;
    std::string qChainTms;                // 체인별 TM (쉼표 구분, qChains 순서)
    std::string tChainTms;                // 체인별 TM (쉼표 구분, tChains 순서)
    std::string interfaceLddt;            // 문자열(다중 값 가능)
    int assId = -1;
};

class MultimerReportParser {
public:
    // _report 파일 로드. 성공 시 true.
    bool load(const std::string& filepath);

    const std::vector<MultimerHit>& get_hits() const { return hits_; }
    int hit_count() const { return (int)hits_.size(); }

private:
    std::vector<MultimerHit> hits_;
};
