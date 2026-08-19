#include "PDBDownloader.hpp"

#include <cstdlib>   // getenv
#include <cstdio>    // popen, pclose
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>  // stat, mkdir
#include <cerrno>

// ─────────────────────────────────────────────────────────────────────────────
// 내부 헬퍼
// ─────────────────────────────────────────────────────────────────────────────

static bool file_exists_nonempty(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    return st.st_size > 0;
}

static bool mkdir_p(const std::string& path) {
    // 단순 단일-단계 mkdir (캐시 디렉토리는 고정 깊이)
    struct stat st;
    if (stat(path.c_str(), &st) == 0) return true;
#if defined(_WIN32)
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0755) == 0 || errno == EEXIST;
#endif
}

// 홈 디렉토리 가져오기
static std::string get_home_dir() {
    const char* home = std::getenv("HOME");
    if (home) return std::string(home);
    return "/tmp";
}

// 캐시 루트 경로 (~/.cache/structty/pdb/)
static std::string cache_root() {
    return get_home_dir() + "/.cache/structty/pdb";
}

// 경로에 디렉토리가 있는지 확인 후 없으면 생성 (최대 3 레벨)
static void ensure_cache_dir() {
    std::string home = get_home_dir();
    mkdir_p(home + "/.cache");
    mkdir_p(home + "/.cache/structty");
    mkdir_p(home + "/.cache/structty/pdb");
}

// 문자열 소문자 변환
static std::string to_lower(const std::string& s) {
    std::string r = s;
    for (char& c : r) c = (char)std::tolower((unsigned char)c);
    return r;
}

// 정규표현식 없이 단순 패턴 매칭 헬퍼
static bool starts_with(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.substr(0, prefix.size()) == prefix;
}

static bool contains(const std::string& s, const std::string& sub) {
    return s.find(sub) != std::string::npos;
}

static bool is_alnum_upper(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9');
}

static bool is_alnum_lower_or_digit(char c) {
    return (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9');
}

// ─────────────────────────────────────────────────────────────────────────────
// PDBDownloader 구현
// ─────────────────────────────────────────────────────────────────────────────

std::string PDBDownloader::extract_target_id(const std::string& raw_target) {
    size_t pos = raw_target.find(' ');
    if (pos == std::string::npos) return raw_target;
    return raw_target.substr(0, pos);
}

DBType PDBDownloader::detect_db_type(const std::string& target_id) {
    if (target_id.empty()) return DBType::Unknown;

    // 패턴 8: TED — AF-..._TED\d+ (AlphaFold보다 먼저 검사)
    if (starts_with(target_id, "AF-") && contains(target_id, "_TED")) {
        // _TED 뒤가 숫자인지 확인
        size_t ted_pos = target_id.rfind("_TED");
        if (ted_pos != std::string::npos) {
            std::string after = target_id.substr(ted_pos + 4);
            bool all_digit = !after.empty();
            for (char c : after) if (c < '0' || c > '9') { all_digit = false; break; }
            if (all_digit) return DBType::TED;
        }
    }

    // 패턴 5: BFVD_Local — _unrelaxed_rank_001_alphafold2 포함
    if (contains(target_id, "_unrelaxed_rank_001_alphafold2")) {
        return DBType::BFVD_Local;
    }

    // 패턴 7: GMGCL — GMGC10.\d+_...\...
    if (starts_with(target_id, "GMGC10.")) {
        return DBType::GMGCL;
    }

    // 패턴 9: BFMD — BFD\d+_.+ (BFD/BFMD Foldseek DB, 예: BFD100_12345678)
    if (starts_with(target_id, "BFD") && target_id.size() > 3) {
        // BFD 뒤에 숫자가 오고 그 뒤 '_' 이후 내용이 있는 형태
        size_t i = 3;
        while (i < target_id.size() && target_id[i] >= '0' && target_id[i] <= '9') i++;
        if (i > 3 && i < target_id.size() && target_id[i] == '_') {
            return DBType::BFMD;
        }
    }

    // 패턴 2: AlphaFold — AF-[A-Z0-9]+-F\d+-model_v\d+
    if (starts_with(target_id, "AF-")) {
        // AF-XXXXX-F\d+-model_v\d+
        size_t dash1 = target_id.find('-', 3);
        if (dash1 != std::string::npos) {
            size_t dash2 = target_id.find('-', dash1 + 1);
            if (dash2 != std::string::npos) {
                std::string frag = target_id.substr(dash1 + 1, dash2 - dash1 - 1);
                if (!frag.empty() && frag[0] == 'F' &&
                    target_id.find("model_v", dash2) != std::string::npos) {
                    return DBType::AlphaFold;
                }
            }
        }
    }

    // 패턴 3: ESMAtlas30 — MGYP + 12자리 숫자
    if (starts_with(target_id, "MGYP") && target_id.size() == 16) {
        bool all_digit = true;
        for (size_t i = 4; i < 16; i++) {
            if (target_id[i] < '0' || target_id[i] > '9') { all_digit = false; break; }
        }
        if (all_digit) return DBType::ESMAtlas30;
    }

    // 패턴 1a: PDB100 — [0-9][a-z0-9]{3}-assembly.*\.cif\.gz_[A-Za-z0-9]+
    // 예: 3a0d-assembly1.cif.gz_A
    if (target_id.size() >= 6 &&
        target_id[0] >= '0' && target_id[0] <= '9' &&
        is_alnum_lower_or_digit(target_id[1]) &&
        is_alnum_lower_or_digit(target_id[2]) &&
        is_alnum_lower_or_digit(target_id[3]) &&
        target_id[4] == '-' &&
        contains(target_id, ".cif.gz_")) {
        return DBType::PDB;
    }

    // 패턴 1: PDB — [0-9][a-z0-9]{3}_[A-Za-z0-9]+
    if (target_id.size() >= 6 &&
        target_id[0] >= '0' && target_id[0] <= '9' &&
        is_alnum_lower_or_digit(target_id[1]) &&
        is_alnum_lower_or_digit(target_id[2]) &&
        is_alnum_lower_or_digit(target_id[3]) &&
        target_id[4] == '_' &&
        target_id.size() >= 6) {
        return DBType::PDB;
    }

    // 패턴 4: CATH50 — [0-9][a-z0-9]{3}[A-Za-z][0-9]{2} (6글자)
    if (target_id.size() == 6 &&
        target_id[0] >= '0' && target_id[0] <= '9' &&
        is_alnum_lower_or_digit(target_id[1]) &&
        is_alnum_lower_or_digit(target_id[2]) &&
        is_alnum_lower_or_digit(target_id[3]) &&
        std::isalpha((unsigned char)target_id[4]) &&
        target_id[5] >= '0' && target_id[5] <= '9') {
        return DBType::CATH50;
    }

    // 패턴 6: BFVD_Official — [A-Z0-9]{6,10}(_[0-9]+)?
    // UniProt accession 형식 (대문자+숫자, 6~10자, 선택적 _숫자 접미사)
    {
        std::string base = target_id;
        size_t underscore = target_id.rfind('_');
        if (underscore != std::string::npos) {
            std::string suffix = target_id.substr(underscore + 1);
            bool all_digit = !suffix.empty();
            for (char c : suffix) if (c < '0' || c > '9') { all_digit = false; break; }
            if (all_digit) base = target_id.substr(0, underscore);
        }
        if (base.size() >= 6 && base.size() <= 10) {
            bool all_upper = true;
            for (char c : base) if (!is_alnum_upper(c)) { all_upper = false; break; }
            if (all_upper) return DBType::BFVD_Official;
        }
    }

    return DBType::Unknown;
}

std::string PDBDownloader::extract_chain(const std::string& target_id, DBType db_type) {
    if (db_type != DBType::PDB) return "-";

    // PDB100 형식: 3a0d-assembly1.cif.gz_A
    size_t cifgz_pos = target_id.find(".cif.gz_");
    if (cifgz_pos != std::string::npos) {
        return target_id.substr(cifgz_pos + 8);  // ".cif.gz_" = 8글자
    }

    // "2xyz_B" → "B"
    // "1a0n_MODEL_1_B" → "B" (마지막 '_' 이후)
    size_t pos = target_id.find('_');
    if (pos == std::string::npos) return "-";

    std::string after = target_id.substr(pos + 1);
    // MODEL_N_X 형태 처리: 마지막 '_' 이후
    size_t last_us = after.rfind('_');
    if (last_us != std::string::npos) {
        after = after.substr(last_us + 1);
    }
    return after.empty() ? "-" : after;
}

std::string PDBDownloader::extract_pdb_id(const std::string& target_id) {
    // PDB100 형식: 3a0d-assembly1.cif.gz_A
    size_t dash_pos = target_id.find('-');
    if (dash_pos == 4) return target_id.substr(0, 4);

    // "2xyz_B" → "2xyz"
    size_t pos = target_id.find('_');
    if (pos == std::string::npos) return target_id.substr(0, 4);
    return target_id.substr(0, std::min(pos, (size_t)4));
}

std::string PDBDownloader::extract_uniprot_id(const std::string& target_id, DBType db_type) {
    if (db_type == DBType::BFVD_Local) {
        // "A0A0N7HVG9_unrelaxed_rank_001_alphafold2_..." → "A0A0N7HVG9"
        size_t pos = target_id.find('_');
        if (pos != std::string::npos) return target_id.substr(0, pos);
        return target_id;
    }
    if (db_type == DBType::BFVD_Official) {
        // "A0A345AIN9_1" → "A0A345AIN9", "A0A345AIN9" → "A0A345AIN9"
        size_t underscore = target_id.rfind('_');
        if (underscore != std::string::npos) {
            std::string suffix = target_id.substr(underscore + 1);
            bool all_digit = !suffix.empty();
            for (char c : suffix) if (c < '0' || c > '9') { all_digit = false; break; }
            if (all_digit) return target_id.substr(0, underscore);
        }
        return target_id;
    }
    return target_id;
}

std::string PDBDownloader::get_download_url(const std::string& target_id, DBType db_type) {
    switch (db_type) {
        case DBType::PDB: {
            std::string pdb_id = extract_pdb_id(target_id);
            // RCSB는 단일 체인 제공 안 함 — 전체 CIF 다운로드
            return "https://files.rcsb.org/download/" + pdb_id + ".cif";
        }
        case DBType::AlphaFold: {
            // AF-P12345-F1-model_v4 형태 그대로 URL에 사용
            return "https://alphafold.ebi.ac.uk/files/" + target_id + ".cif";
        }
        case DBType::ESMAtlas30: {
            return "https://api.esmatlas.com/fetchPredictedStructure/" + target_id;
        }
        case DBType::CATH50: {
            return "https://www.cathdb.info/version/v4_3_0/api/rest/id/" + target_id + ".pdb";
        }
        case DBType::BFVD_Local:
        case DBType::BFVD_Official: {
            // UniProt accession 추출 후 BFVD API 사용
            std::string uniprot = extract_uniprot_id(target_id, db_type);
            return "https://bfvd.steineggerlab.workers.dev/pdb/" + uniprot + ".pdb";
        }
        case DBType::TED: {
            // TED API: target ID 그대로 사용 (accept 헤더는 resolve_target_file에서 처리)
            return "https://ted.cathdb.info/api/v1/files/" + target_id + ".pdb";
        }
        // 다운로드 URL이 없는 DB들
        case DBType::GMGCL:
        case DBType::BFMD:
        case DBType::Unknown:
        default:
            return "";
    }
}

std::string PDBDownloader::get_cache_path(const std::string& target_id, DBType db_type) {
    ensure_cache_dir();
    std::string root = cache_root();

    switch (db_type) {
        case DBType::PDB: {
            std::string pdb_id = extract_pdb_id(target_id);
            return root + "/" + pdb_id + ".cif";
        }
        case DBType::AlphaFold:
            return root + "/" + target_id + ".cif";
        case DBType::ESMAtlas30:
            return root + "/" + target_id + ".pdb";
        case DBType::CATH50:
            return root + "/" + target_id + ".pdb";
        case DBType::BFVD_Local:
        case DBType::BFVD_Official: {
            // UniProt accession 기준으로 캐시 (전체 target ID는 너무 김)
            std::string uniprot = extract_uniprot_id(target_id, db_type);
            return root + "/" + uniprot + ".pdb";
        }
        case DBType::TED:
            return root + "/" + target_id + ".pdb";
        default:
            return root + "/" + target_id + ".pdb";
    }
}

// <db_path>/<name> 에 붙일 확장자 후보를 순서대로 시도한다.
// 시도 순서: .pdb, .cif, .ent, 각각의 .gz, 확장자 없음 — 모두 소문자 이름으로도 한 번 더.
static std::string find_by_name(const std::string& name, const std::string& db_path) {
    const std::string base       = db_path + "/" + name;
    const std::string base_lower = db_path + "/" + to_lower(name);

    const std::vector<std::string> candidates = {
        base + ".pdb",       base + ".cif",       base + ".ent",
        base + ".pdb.gz",    base + ".cif.gz",    base + ".ent.gz",
        base,                // 확장자 없는 경우
        base_lower + ".pdb", base_lower + ".cif", base_lower + ".ent",
        base_lower + ".pdb.gz", base_lower + ".cif.gz", base_lower + ".ent.gz",
        base_lower,
    };

    for (const std::string& cand : candidates) {
        if (file_exists_nonempty(cand)) return cand;
    }
    return "";
}

std::string PDBDownloader::find_in_db_path(const std::string& target_id,
                                            const std::string& db_path) {
    if (db_path.empty()) return "";

    const std::string direct = find_by_name(target_id, db_path);
    if (!direct.empty()) return direct;

    // foldseek createdb 는 멀티머를 체인마다 쪼개 `<파일 stem>_<auth chain>` accession 을
    // 만든다(예: 1dci-assembly1_B-2). 디렉터리에는 원본 파일 하나(`1dci-assembly1.cif`)만
    // 있으므로 뒤쪽 `_<chain>` 을 하나씩 떼어내며 다시 찾는다. 어느 체인을 그릴지는
    // Screen::chain_from_accession() 이 찾아낸 파일 stem 과 accession 을 비교해 정한다.
    // stem 자체에 `_` 가 있을 수 있어(`my_protein_A`) 뒤에서부터 순서대로 자른다.
    size_t cut = target_id.rfind('_');
    while (cut != std::string::npos && cut > 0) {
        const std::string found = find_by_name(target_id.substr(0, cut), db_path);
        if (!found.empty()) return found;
        if (cut == 0) break;
        cut = target_id.rfind('_', cut - 1);
    }
    return "";
}

bool PDBDownloader::download_file(const std::string& url, const std::string& dest_path,
                                  const std::string& extra_header) {
    if (url.empty() || dest_path.empty()) return false;

    // curl 우선, 없으면 wget
    // -f: fail silently on HTTP errors, -s: silent, -L: follow redirects
    std::string cmd;
    if (extra_header.empty()) {
        cmd = "curl -f -s -L -o \"" + dest_path + "\" \"" + url + "\" 2>/dev/null"
              " || wget -q -O \"" + dest_path + "\" \"" + url + "\" 2>/dev/null";
    } else {
        cmd = "curl -f -s -L -H \"" + extra_header + "\" -o \"" + dest_path + "\" \"" + url + "\" 2>/dev/null"
              " || wget -q --header=\"" + extra_header + "\" -O \"" + dest_path + "\" \"" + url + "\" 2>/dev/null";
    }

    FILE* p = popen(cmd.c_str(), "r");
    if (!p) return false;
    pclose(p);

    return file_exists_nonempty(dest_path);
}

std::string PDBDownloader::get_no_url_message(DBType db_type,
                                               const std::string& target_id) {
    switch (db_type) {
        case DBType::GMGCL:
            return "GMGCL: no download URL available. Web-server only DB.";
        case DBType::BFMD:
            return "BFMD: no download URL. Use --db-path.";
        default:
            return "File not found: " + target_id;
    }
}

std::string PDBDownloader::resolve_target_file(const std::string& target_id,
                                               const std::string& db_path,
                                               std::string& status_msg) {
    status_msg.clear();
    DBType db_type = detect_db_type(target_id);

    // 1. db_path가 있으면 거기서 먼저 탐색
    if (!db_path.empty()) {
        std::string found = find_in_db_path(target_id, db_path);
        if (!found.empty()) return found;
    }

    // 2. 캐시에 이미 있는지 확인
    std::string cache_path = get_cache_path(target_id, db_type);
    if (file_exists_nonempty(cache_path)) return cache_path;

    // 3. 다운로드 URL이 있으면 다운로드 시도
    std::string url = get_download_url(target_id, db_type);
    if (!url.empty()) {
        // TED API는 accept 헤더 필요
        std::string header;
        if (db_type == DBType::TED) {
            header = "accept: application/octet-stream";
        }
        bool ok = download_file(url, cache_path, header);
        if (ok) return cache_path;
        // 실패 이유 메시지
        if (db_type == DBType::ESMAtlas30) {
            status_msg = "ESMAtlas API rate limit, retry later";
        } else {
            status_msg = "Download failed: " + target_id;
        }
        return "";
    }

    // 4. URL 없는 DB: 안내 메시지
    status_msg = get_no_url_message(db_type, target_id);
    return "";
}
