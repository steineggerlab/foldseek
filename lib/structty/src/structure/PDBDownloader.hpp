#pragma once

#include <string>

// plan.md 기능 3: Foldseek 결과 target ID 패턴 감지 및 구조 파일 다운로드
//
// 지원하는 DB 종류 (target ID 패턴 기반):
//   PDB          — ^[0-9][a-z0-9]{3}_[A-Za-z0-9]+$
//   AlphaFold    — ^AF-[A-Z0-9]+-F[0-9]+-model_v[0-9]+$
//   ESMAtlas30   — ^MGYP[0-9]{12}$
//   CATH50       — ^[0-9][a-z0-9]{3}[A-Za-z][0-9]{2}$
//   BFVD_Local   — _unrelaxed_rank_001_alphafold2 포함 → UniProt 추출 → bfvd.steineggerlab.workers.dev
//   BFVD_Official— ^[A-Z0-9]{6,10}(_[0-9]+)?$ → UniProt 추출 → bfvd.steineggerlab.workers.dev
//   GMGCL        — ^GMGC10\.[0-9_]+\..+$ (다운로드 URL 없음)
//   TED          — ^AF-.*_TED[0-9]+$ → ted.cathdb.info (accept 헤더 필요)
//   BFMD         — ^BFD[0-9]+_.+$ (BFD/BFMD DB, 다운로드 URL 없음, --db-path 필요)
//   Unknown      — 그 외

enum class DBType {
    PDB,
    AlphaFold,
    ESMAtlas30,
    CATH50,
    BFVD_Local,
    BFVD_Official,
    GMGCL,
    TED,
    BFMD,
    Unknown
};

class PDBDownloader {
public:
    // target 컬럼 raw 값에서 첫 공백 이전 토큰 추출
    static std::string extract_target_id(const std::string& raw_target);

    // target ID 패턴으로 DB 종류 감지
    static DBType detect_db_type(const std::string& target_id);

    // PDB 타겟에서 체인 ID 추출 (예: "2xyz_B" → "B", "1a0n_MODEL_1_B" → "B")
    // AlphaFold/ESMAtlas 등 단일 체인은 "-" 반환
    static std::string extract_chain(const std::string& target_id, DBType db_type);

    // PDB 타겟에서 4자리 PDB ID 추출 (예: "2xyz_B" → "2xyz")
    static std::string extract_pdb_id(const std::string& target_id);

    // BFVD target ID에서 UniProt accession 추출
    //   BFVD_Local : "A0A0N7HVG9_unrelaxed_rank_001_..." → "A0A0N7HVG9"
    //   BFVD_Official: "A0A345AIN9_1" → "A0A345AIN9", "A0A345AIN9" → "A0A345AIN9"
    static std::string extract_uniprot_id(const std::string& target_id, DBType db_type);

    // 다운로드 URL 생성 (URL이 없는 DB는 빈 문자열 반환)
    static std::string get_download_url(const std::string& target_id, DBType db_type);

    // 캐시 파일 경로 반환 (~/.cache/structty/pdb/{id}.{ext})
    // BFVD: UniProt accession 기준 파일명 사용
    // 디렉토리가 없으면 생성한다
    static std::string get_cache_path(const std::string& target_id, DBType db_type);

    // db_path 내에서 target 파일 탐색 (여러 확장자 시도)
    // 발견 시 경로 반환, 없으면 빈 문자열
    static std::string find_in_db_path(const std::string& target_id,
                                       const std::string& db_path);

    // 파일 다운로드 (curl → wget 순서로 시도, popen 기반)
    // dest_path: 저장 경로, url: 다운로드 URL
    // extra_header: HTTP 헤더 추가 (예: "accept: application/octet-stream"), 빈 문자열이면 생략
    // 성공 시 true, 실패 시 false
    static bool download_file(const std::string& url, const std::string& dest_path,
                              const std::string& extra_header = "");

    // target에 대한 로컬 파일 경로 결정:
    //   1. db_path가 있으면 거기서 탐색
    //   2. 캐시에 이미 있으면 캐시 경로 반환
    //   3. 다운로드 URL이 있으면 캐시에 다운로드
    //   4. 모두 실패하면 빈 문자열 + status_msg 설정
    // status_msg: 실패 이유나 특수 안내 메시지 (패널 표시용)
    static std::string resolve_target_file(const std::string& target_id,
                                           const std::string& db_path,
                                           std::string& status_msg);

    // DB 종류에 따른 "다운로드 불가" 안내 메시지
    static std::string get_no_url_message(DBType db_type, const std::string& target_id);
};
