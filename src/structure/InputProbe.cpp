#include "InputProbe.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace input_probe {

namespace {

// 소문자 사본 — 확장자 비교용
std::string to_lower(const std::string& s) {
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });
    return out;
}

bool ends_with(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    return s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// 구조 파일 확장자. Protein.cpp:271 이 ".pdb"/".cif" 부분 문자열로 판정하는 것과
// 어긋나지 않게 같은 두 확장자를 포함하고, gemmi 가 읽는 .ent 와 .gz 도 받는다.
bool has_structure_extension(const std::string& path) {
    std::string p = to_lower(path);
    if (ends_with(p, ".gz")) {
        p.erase(p.size() - 3);
    }
    return ends_with(p, ".pdb") || ends_with(p, ".cif") || ends_with(p, ".ent");
}

bool exists_file(const std::string& path) {
    std::error_code ec;
    return fs::is_regular_file(path, ec);
}

bool is_directory(const std::string& path) {
    std::error_code ec;
    return fs::is_directory(path, ec);
}

// 첫 비주석·비공백 줄. 없으면 빈 문자열.
std::string first_data_line(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) return std::string();

    std::string line;
    while (std::getline(ifs, line)) {
        const size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;   // 공백만
        if (line[first] == '#') continue;           // 주석
        return line;
    }
    return std::string();
}

// 컬럼 수 → 결과 파일 종류.
// FoldseekParser::load 의 포맷 감지와 같은 규칙이어야 한다:
//   12 / 17 / 21 / 29 를 각각 인식하고, taxname 에 탭이 섞일 수 있으므로
//   21 초과도 alis(21)로 받아준다(FoldseekParser.cpp `if (ncols > 21)`).
// 멀티머 `_report` 는 14컬럼으로 이 집합과 겹치지 않는다
// (MultimerReportParser.hpp 의 14컬럼 레이아웃).
InputKind kind_from_columns(int ncols) {
    if (ncols == 14) return InputKind::ResultReport;
    if (ncols == 12 || ncols == 17 || ncols >= 21) return InputKind::ResultM8;
    return InputKind::Unknown;
}

// 디렉터리는 대표 파일 1개만 확인한다(전수 검사는 큰 DB 디렉터리에서 비싸다).
// 결정적 결과를 위해 이름 정렬 후 첫 정규 파일을 고른다.
std::string representative_file(const std::string& dir) {
    std::error_code ec;
    std::vector<std::string> names;
    for (const auto& entry : fs::directory_iterator(dir, ec)) {
        if (ec) return std::string();
        std::error_code file_ec;
        if (!entry.is_regular_file(file_ec)) continue;
        names.push_back(entry.path().string());
    }
    if (names.empty()) return std::string();
    return *std::min_element(names.begin(), names.end());
}

}  // namespace

int tsv_column_count(const std::string& path) {
    const std::string line = first_data_line(path);
    if (line.empty()) return 0;

    int fields = 1;
    for (const char c : line) {
        if (c == '\t') ++fields;
    }
    return fields;
}

bool db_has_ca(const std::string& db_path) {
    return exists_file(db_path + "_ca.dbtype");
}

InputKind probe(const std::string& path) {
    if (path.empty()) return InputKind::Unknown;

    // (a) Foldseek DB — 데이터 파일 옆에 `<path>.dbtype` 이 함께 있다
    if (exists_file(path + ".dbtype")) return InputKind::FoldseekDB;

    // (b) 디렉터리 — 대표 파일로 구조 디렉터리인지 판단
    if (is_directory(path)) {
        const std::string sample = representative_file(path);
        if (sample.empty()) return InputKind::Unknown;
        if (has_structure_extension(sample)) return InputKind::StructureDir;
        const std::string line = first_data_line(sample);
        if (!line.empty() && line[0] == '>') return InputKind::SequenceFasta;
        return InputKind::Unknown;
    }

    if (!exists_file(path)) return InputKind::Unknown;

    // (c) 구조 확장자 — 내용을 읽기 전에 확장자로 끊는다.
    //     .gz 는 평문으로 읽을 수 없고, mmCIF 본문에 탭이 섞이면
    //     컬럼 수 판정이 오작동할 수 있다.
    if (has_structure_extension(path)) return InputKind::Structure;

    const std::string line = first_data_line(path);
    if (line.empty()) return InputKind::Unknown;

    // (d) 서열 FASTA — 첫 데이터 줄이 '>' 로 시작
    if (line[0] == '>') return InputKind::SequenceFasta;

    // (e) Foldseek 결과 — 탭 컬럼 수로 m8 / _report 구분
    return kind_from_columns(tsv_column_count(path));
}

const char* kind_name(const InputKind kind) {
    switch (kind) {
        case InputKind::Structure:     return "structure file";
        case InputKind::StructureDir:  return "structure directory";
        case InputKind::FoldseekDB:    return "Foldseek DB";
        case InputKind::ResultM8:      return "Foldseek m8";
        case InputKind::ResultReport:  return "Foldseek multimer report";
        case InputKind::SequenceFasta: return "sequence FASTA";
        case InputKind::Unknown:       break;
    }
    return "unknown";
}

namespace {

// 좌표가 없는 입력을 렌더 시작 전에 거른다.
// 서열 FASTA: foldseek createdb 의 ProstT5 경로는 3Di(`_ss`)만 예측하고 `_ca` 를 만들지 않는다.
// `_ca` 없는 DB: 서열 유래이거나 `--index-exclude 2` 로 만든 DB.
bool check_has_coordinates(const std::string& path, const std::string& role) {
    const InputKind kind = probe(path);

    if (kind == InputKind::SequenceFasta) {
        std::cerr << "Error: " << role << " '" << path << "' is a sequence FASTA.\n"
                  << "       StrucTTY needs 3D coordinates: pass PDB/mmCIF files, a directory\n"
                  << "       of them, or a Foldseek DB built from structures.\n"
                  << "       (foldseek createdb --prostt5-model predicts 3Di but writes no _ca)"
                  << std::endl;
        return false;
    }
    if (kind == InputKind::FoldseekDB && !db_has_ca(path)) {
        std::cerr << "Error: Foldseek DB '" << path << "' has no C-alpha coordinates ("
                  << path << "_ca is missing).\n"
                  << "       It was built from sequences, or with --index-exclude 2."
                  << std::endl;
        return false;
    }
    return true;
}

}  // namespace

bool validate_inputs(const std::vector<std::string>& query,
                     const std::string& target, const std::string& result) {
    if (query.empty()) {
        std::cerr << "Error: Need input file dir" << std::endl;
        return false;
    }

    // 쌍 규칙: target 소스 없이 결과만 주면 hit 구조를 읽을 수 없고,
    // 결과 없이 target 만 주면 어떤 hit 을 띄울지 알 수 없다.
    if (!target.empty() && result.empty()) {
        std::cerr << "Error: -fst / --foldseek-target given without -fsr / --foldseek-result.\n"
                  << "       Both must be given together." << std::endl;
        return false;
    }
    if (target.empty() && !result.empty()) {
        std::cerr << "Error: -fsr / --foldseek-result given without -fst / --foldseek-target.\n"
                  << "       Both must be given together. Use -fst auto to download hit structures."
                  << std::endl;
        return false;
    }

    // "auto" 는 예약값이라 경로 판별을 건너뛴다(다운로드 모드).
    for (const std::string& q : query) {
        if (!check_has_coordinates(q, "query")) {
            return false;
        }
    }
    if (!target.empty() && target != "auto" && !check_has_coordinates(target, "-fst target")) {
        return false;
    }

    if (result.empty()) {
        return true;
    }

    const InputKind result_kind = probe(result);
    if (result_kind != InputKind::ResultM8 && result_kind != InputKind::ResultReport) {
        const int ncols = tsv_column_count(result);
        std::cerr << "Error: -fsr '" << result << "' is not a Foldseek result file.\n"
                  << "       Expected tab-separated 12/17/21/29 columns (m8) or 14 columns"
                  << " (multimer _report),\n       ";
        if (ncols == 0) {
            std::cerr << "but the file could not be read or has no data rows.";
        } else {
            std::cerr << "but found " << ncols << " columns.";
        }
        std::cerr << std::endl;
        return false;
    }

    return true;
}

}  // namespace input_probe
