#pragma once

#include <string>
#include <vector>

// 기능 3/4: Foldseek 결과 파일(.m8) 파싱
// 지원 포맷:
//   12 cols: 기본 m8 (qaln/taln/U/T 없음) — foldseek easy-search 기본 출력
//   17 cols: 축소 형식 (qaln/taln 포함, U/T 없음)
//            = foldseek --format-output query,target,fident,alnlen,mismatch,gapopen,
//              qstart,qend,tstart,tend,evalue,bits,lddt,qtmscore,ttmscore,qaln,taln
//            foldseek `--view-structty` 가 뷰어에 넘기는 viewer m8 포맷.
//            컬럼 순서가 이와 다르면 lddt/qtmscore 가 뒤바뀌고 qaln/taln 을 숫자로 읽는다.
//            생성에는 검색 단계의 backtrace(`-a`)가 필요하다.
//   21 cols: alis 포맷 (alns, tseq, taxid, taxname 포함, U/T 없음) ← Foldseek easy-search 기본 출력
//   29 cols: 전체 형식 (U(9), T(3), qaln, taln 포함)

struct FoldseekHit {
    std::string query;
    std::string target;   // 첫 번째 공백 이전 토큰만 저장
    float fident    = 0.0f;
    int   alnlen    = 0;

    // cols 6-9 (qstart, qend, tstart, tend) — 잔기 번호 매핑에 필요
    int   qstart    = 1;
    int   qend      = 0;
    int   tstart    = 1;  // taln의 첫 번째 non-gap 잔기 번호
    int   tend      = 0;

    // 21컬럼: col[10]=prob, col[11]=evalue  /  29컬럼: col[10]=evalue
    float prob      = -1.0f;  // alis 21컬럼 포맷: 구조 유사도 확률 (0~1)
    float evalue    = 0.0f;
    float bits      = 0.0f;

    int   qlength   = 0;      // alis 21컬럼: col[13] = 쿼리 전체 길이
    int   tlength   = 0;      // alis 21컬럼: col[14] = 타겟 전체 길이

    // 29컬럼 전용 필드
    float lddt      = -1.0f;
    float qtmscore  = -1.0f;
    float ttmscore  = -1.0f;

    // U: row-major 3x3 rotation matrix (flat, 9 elements)
    // T: translation vector (3 elements)
    // 29컬럼: Foldseek에서 직접 제공 / alis 21컬럼: load_next_hit()에서 Kabsch SVD로 채워짐
    float U[9] = {};
    float T[3] = {};
    bool has_transform = false;

    std::string qaln;   // 쿼리 alignment 문자열 (gap='-')
    std::string taln;   // 타겟 alignment 문자열 (gap='-')
    bool has_aln = false;

    // alis 21컬럼 전용 필드
    std::vector<float> alns;  // 정렬된 타겟 CA 좌표 (query frame 기준, 3N floats)
    std::string tseq;         // 타겟 전체 서열
    std::string taxid;        // taxonomy ID
    std::string taxname;      // taxonomy 이름 (내부 공백 포함 가능)
    bool is_alis_format = false;  // 21컬럼 alis 포맷 여부
};

class FoldseekParser {
public:
    // filepath: Foldseek easy-search 출력 파일(.m8 형식, 탭 구분, 헤더 없음)
    // 반환: 성공 시 true, 실패 시 false
    bool load(const std::string& filepath);

    const std::vector<FoldseekHit>& get_hits() const { return hits; }
    int hit_count() const { return (int)hits.size(); }

private:
    std::vector<FoldseekHit> hits;

    // 컬럼 수 기반 포맷 자동 감지 후 파싱
    bool parse_line(const std::vector<std::string>& cols, FoldseekHit& hit, int fmt);
};
