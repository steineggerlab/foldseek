#pragma once

#include <string>
#include <vector>

// 입력 경로의 종류 판별 — CLI 검증(`-fst`/`-fsr` 쌍, 좌표 없는 입력 차단)에서 쓴다.
// StrucTTY 가 다루는 입력은 4종뿐이다:
//   구조 파일(PDB/mmCIF) · 구조 파일 디렉터리 · Foldseek DB · Foldseek 결과(m8 / 멀티머 _report)
// 서열 FASTA 는 좌표가 없어 렌더할 수 없다(foldseek createdb 의 ProstT5 경로는 `_ss` 만
// 만들고 `_ca` 를 만들지 않는다). 따라서 별도 종류로 판별해 조기에 거부한다.
//
// 이름 충돌 방지: 이 라이브러리는 foldseek/mmseqs 에 정적 링크되므로
// `probe`/`db_has_ca` 같은 일반적인 이름을 전역에 두지 않고 namespace 안에 둔다.
namespace input_probe {

enum class InputKind {
    Structure,      // .pdb/.cif/.ent (+.gz)
    StructureDir,   // 구조 파일 디렉터리
    FoldseekDB,     // <path>.dbtype 존재
    ResultM8,       // 탭 12/17/21(+)/29 컬럼
    ResultReport,   // 탭 14 컬럼 (멀티머 _report)
    SequenceFasta,  // 첫 비어있지 않은 줄이 '>' — 미지원
    Unknown
};

// 판별 순서: Foldseek DB → 디렉터리 → 구조 확장자 → FASTA → 탭 컬럼 수 → Unknown
// 열기 실패·빈 파일·해석 불가는 예외 없이 Unknown 을 돌려준다(프로젝트 관례).
InputKind probe(const std::string& path);

// Foldseek DB 에 Cα 좌표 부분(`<db>_ca`)이 있는가.
// 서열 FASTA 유래 DB 나 `--index-exclude 2` 로 만든 DB 는 이 파일이 없다.
bool db_has_ca(const std::string& db_path);

// 첫 비주석(`#`)·비공백 줄의 탭 필드 수. 열기 실패나 해당 줄이 없으면 0.
int tsv_column_count(const std::string& path);

// 로그·에러 메시지용 짧은 이름
const char* kind_name(InputKind kind);

// 렌더할 수 없는 입력 조합을 걸러낸다. CLI(Parameters)와 라이브러리 진입점
// (structty::run)이 같은 검사를 쓰도록 여기에 둔다. 실패 시 이유를 stderr 로
// 출력하고 false 를 돌려준다.
bool validate_inputs(const std::vector<std::string>& query,
                     const std::string& target, const std::string& result);

}  // namespace input_probe
