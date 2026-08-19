#include <iostream>
#include <clocale>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <memory>
#include <map>
#include <vector>
#include <string>
#include "structty.h"
#include "Terminal.hpp"
#include "Common.hpp"
#include "Screen.hpp"
#include "MSAParser.hpp"
#include "FoldseekParser.hpp"
#include "MultimerReportParser.hpp"
#include "FoldMasonParser.hpp"
#include "InputProbe.hpp"
#include "PDBDownloader.hpp"
#include "Benchmark.hpp"

namespace fs = std::filesystem;

namespace structty {

namespace {

// 결과 파일 첫 데이터 행의 target 컬럼(2번째). `-fst auto` 가 그 accession 을 공개 DB
// 패턴으로 해석할 수 있는지 미리 알리는 데만 쓴다. 실패하면 빈 문자열.
std::string first_target_accession(const std::string& result_path) {
    std::ifstream ifs(result_path);
    if (!ifs.is_open()) return std::string();

    std::string line;
    while (std::getline(ifs, line)) {
        const size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;

        const size_t tab1 = line.find('\t');
        if (tab1 == std::string::npos) return std::string();
        const size_t tab2 = line.find('\t', tab1 + 1);
        const size_t end = (tab2 == std::string::npos) ? line.size() : tab2;
        return line.substr(tab1 + 1, end - tab1 - 1);
    }
    return std::string();
}

// 기능 4: align 계열 모드는 정렬 문자열(qaln/taln)이 필요하다. 12컬럼 m8 처럼 정렬
// 문자열이 없으면 align 은 최근접 이웃 판정으로 대체하고, align-fs 는 진행하지 않는다.
// 진행해도 되면 true.
bool check_alignment_available(const std::vector<FoldseekHit>& hits,
                               const std::string& mode,
                               const std::string& result_path) {
    if (mode != "align" && mode != "align-fs") return true;

    for (const FoldseekHit& hit : hits) {
        if (hit.has_aln) return true;
    }

    const bool strict = (mode == "align-fs");
    std::cerr << (strict ? "Error: -m align-fs needs alignment strings (qaln/taln), "
                         : "Warning: -m align needs alignment strings (qaln/taln), ")
              << "but none are present in " << result_path << ".\n"
              << (strict ? "       Use -m align-near for distance-based colouring, "
                           "or regenerate the result with:\n"
                         : "         Falling back to nearest-neighbour (10 A) colouring "
                           "-- the panel shows 'nearest-nbr'.\n"
                           "         Regenerate the result with:\n")
              << "           foldseek convertalis <queryDB> <targetDB> <resultDB> out.m8 \\\n"
              << "             --format-output query,target,fident,alnlen,mismatch,gapopen,"
              << "qstart,qend,tstart,tend,evalue,bits,lddt,qtmscore,ttmscore,qaln,taln\n"
              << "         (the search itself must run with -a so backtraces exist)"
              << std::endl;
    return !strict;
}

}  // namespace

bool run(const RunOptions& opts) {
    if (!input_probe::validate_inputs(opts.input_files, opts.foldseek_target,
                                      opts.foldseek_result)) {
        return false;
    }

    setlocale(LC_ALL, "");
    Terminal::enter_raw_mode();

    Terminal::Size term_sz = Terminal::get_size();
    int term_rows = term_sz.rows;
    int term_cols = term_sz.cols;

    Screen screen(term_cols, term_rows, opts.show_structure, opts.mode);

    Benchmark bm;
    const bool bench = opts.benchmark;
    using BenchClock = Benchmark::clock;
    BenchClock::time_point t_load0;

    if (bench) {
        fs::path in_file(opts.input_files[0]);
        bm.start(std::filesystem::current_path().string() + "/structty_bench_"
                 + current_timestamp() + "_" + in_file.stem().string() + ".csv");
        std::cout << "CWD = " << std::filesystem::current_path() << std::endl;
        screen.set_benchmark(&bm);
        t_load0 = Benchmark::clock::now();
    }

    screen.set_chainfile(opts.chains_file, (int)opts.input_files.size());

    // -fst / -fsr 로 들어온 입력의 종류를 판별해 씬 경로를 정한다.
    // 과거의 --db / --db-path / --query-db / --report-format 조합을 대체한다:
    //   query 가 Foldseek DB   → query-from-DB (구 --query-db)
    //   -fst 가 Foldseek DB    → CA 직접 읽기 (구 --db)
    //   -fst 가 디렉터리·구조  → accession 이름으로 로컬 탐색 (구 --db-path)
    //   -fst auto              → PDBDownloader 로 내려받기
    //   -fsr 이 14컬럼 _report → 멀티머 경로 (구 --report-format)
    using input_probe::InputKind;
    const std::string& fs_result   = opts.foldseek_result;
    const bool         target_auto = (opts.foldseek_target == "auto");
    const InputKind result_kind = fs_result.empty()
                                    ? InputKind::Unknown : input_probe::probe(fs_result);
    const InputKind query_kind  = opts.input_files.empty()
                                    ? InputKind::Unknown : input_probe::probe(opts.input_files[0]);
    const InputKind target_kind = (opts.foldseek_target.empty() || target_auto)
                                    ? InputKind::Unknown : input_probe::probe(opts.foldseek_target);

    const std::string query_db  = (query_kind == InputKind::FoldseekDB)
                                    ? opts.input_files[0] : std::string();
    const std::string target_db = (target_kind == InputKind::FoldseekDB)
                                    ? opts.foldseek_target : std::string();
    // PDBDownloader::find_in_db_path() 는 디렉터리를 받는다. 구조 파일 하나를 -fst 로 준
    // 경우에는 그 파일이 있는 디렉터리를 탐색 대상으로 삼는다.
    std::string target_dir;
    if (target_kind == InputKind::StructureDir) {
        target_dir = opts.foldseek_target;
    } else if (target_kind == InputKind::Structure) {
        target_dir = fs::path(opts.foldseek_target).parent_path().string();
        if (target_dir.empty()) target_dir = ".";
    }

    // 해석된 입력 1줄 요약 — 자동 판별 결과를 눈으로 확인할 수 있어야 한다
    // (판별이 틀리면 씬 경로 자체가 달라지므로).
    if (!fs_result.empty() || !opts.foldseek_target.empty()) {
        // stderr 로 낸다: stdout 은 렌더 출력(raw write)과 std::cout 이 섞여 있어서
        // 여기서 cout 을 flush 하면 프레임 앞부분 개행 위치가 달라진다.
        std::cerr << "Foldseek input: query=" << input_probe::kind_name(query_kind)
                  << ", target=";
        if (target_auto) {
            std::cerr << "auto (download)";
        } else {
            std::cerr << input_probe::kind_name(target_kind);
        }
        std::cerr << ", result=" << input_probe::tsv_column_count(fs_result) << " cols ("
                  << input_probe::kind_name(result_kind) << "), mode=" << opts.mode
                  << std::endl;

        // auto 는 공개 DB accession 패턴만 해석한다. foldseek createdb 가 만든
        // 체인 accession(`<stem>_<chain>`)은 어느 패턴에도 맞지 않아 URL 이 없다.
        if (target_auto) {
            const std::string acc = first_target_accession(fs_result);
            if (!acc.empty() && PDBDownloader::detect_db_type(acc) == DBType::Unknown) {
                std::cerr << "Warning: -fst auto cannot resolve hit accession '" << acc
                          << "' to a public database URL.\n"
                          << "         Pass a Foldseek DB or a structure directory to -fst"
                          << " instead; the panel shows the per-hit reason." << std::endl;
            }
        }
    }

    // Query-from-DB (D3/D4) + multi-query nav (D7/D8/D9): when the query input is
    // a Foldseek DB, read the query structure(s) from it instead of parsing a
    // structure file. Handles folder/tar/gz query inputs that gemmi's single-file
    // reader cannot open. Hits are grouped by the .m8 query column so ]/[ can
    // navigate between queries.
    // Step 7 (D6): multimer `_report` 경로. 14컬럼 tsv 를 complex 단위로 파싱해
    // query/target complex 전체 체인을 로드하고 complex U/T 로 겹침. query/target DB 는
    // complex DB(체인별 엔트리). -fsr 이 _report 일 때만 진입(.m8 경로와 배타적).
    bool multimer_report = false;
    if (result_kind == InputKind::ResultReport && opts.mode == "align-fs") {
        std::cerr << "Error: -m align-fs needs alignment strings, but a multimer report "
                  << "carries none.\n"
                  << "       Use -m align-near, or run a monomer search with qaln/taln."
                  << std::endl;
        return false;
    }
    if (result_kind == InputKind::ResultReport) {
        MultimerReportParser mr_parser;
        if (mr_parser.load(fs_result) && mr_parser.hit_count() > 0) {
            if (!target_db.empty()) {
                screen.open_foldseek_db(target_db);
            } else {
                screen.set_fs_db_path(target_dir);
            }
            screen.set_multimer_report(mr_parser.get_hits(), opts.input_files[0],
                                       query_kind == InputKind::FoldseekDB,
                                       query_kind == InputKind::StructureDir,
                                       target_db, opts.show_structure);
            multimer_report = true;
        } else {
            std::cerr << "Warning: failed to parse multimer report: "
                      << fs_result << std::endl;
        }
    }

    bool query_from_db = false;
    std::vector<std::string> query_ids;                              // .m8 순서
    std::map<std::string, std::vector<FoldseekHit>> hits_by_query;   // query별 hit 그룹
    if (!multimer_report && !query_db.empty() && result_kind == InputKind::ResultM8) {
        FoldseekParser q_parser;
        if (q_parser.load(fs_result) && q_parser.hit_count() > 0) {
            for (const FoldseekHit& h : q_parser.get_hits()) {
                auto it = hits_by_query.find(h.query);
                if (it == hits_by_query.end()) {
                    query_ids.push_back(h.query);
                    hits_by_query.emplace(h.query, std::vector<FoldseekHit>{ h });
                } else {
                    it->second.push_back(h);
                }
            }
            query_from_db = !query_ids.empty();
            if (query_from_db &&
                !check_alignment_available(q_parser.get_hits(), opts.mode, fs_result)) {
                return false;
            }
        }
    }

    // query 가 Foldseek DB 인데 위 두 경로 어느 쪽도 서지 않으면 더 진행할 수 없다.
    // plaintext 경로는 gemmi 로 구조 파일을 여는데, DB 데이터 파일은 열 수 없다.
    if (!query_db.empty() && !multimer_report && !query_from_db) {
        std::cerr << "Error: query '" << query_db << "' is a Foldseek DB, but no usable hits "
                  << "were read from " << (fs_result.empty() ? "(no -fsr given)" : fs_result)
                  << ".\n       A Foldseek DB query needs a matching -fsr result file."
                  << std::endl;
        Terminal::exit_raw_mode();
        return false;
    }

    if (multimer_report) {
        // Scene already set up by set_multimer_report() above (complex chains
        // loaded + superposed). Nothing more to do here.
    } else if (query_from_db) {
        // Screen owns the full per-query setup (load query from DB, normalize,
        // open target DB, load first hit). Query switching via ]/[ reuses the
        // same path. The workflow handoff passes only the query as input_files,
        // so plaintext targets are not loaded here.
        screen.set_query_nav(query_ids, hits_by_query, query_db,
                             target_db, opts.show_structure);
    } else {
    // plaintext 경로: 체인 필터를 구조 로드 전에 정해야 하므로 m8 을 먼저 읽는다.
    // foldseek createdb 는 멀티머를 체인마다 쪼개므로(`<stem>_<chain>`), m8 accession 이
    // 가리키는 체인만 남겨야 qaln/taln 인덱스가 파일의 CA 순서와 맞는다.
    FoldseekParser fs_nav_parser;
    bool fs_hits_loaded = false;
    if (!fs_result.empty()) {
        fs_hits_loaded = fs_nav_parser.load(fs_result) &&
                         fs_nav_parser.hit_count() > 0;
        if (!fs_hits_loaded) {
            std::cerr << "Warning: no usable hits in Foldseek result: "
                      << fs_result << std::endl;
        }

        // target 소스가 없다는 경고는 필요 없어졌다 — -fsr 은 -fst 와 쌍으로만 받으므로
        // (Parameters 검증 2) 결과 파일이 있으면 target 소스도 항상 있다.

        if (fs_hits_loaded && !check_alignment_available(fs_nav_parser.get_hits(),
                                                        opts.mode, fs_result)) {
            return false;
        }
    }

    // easy-search 는 query 디렉터리 전체의 hit 을 m8 하나에 쓴다. 지금 열어둔 구조와
    // 무관한 query 의 hit 까지 [N]ext 로 순회하면 엉뚱한 구조가 겹쳐 뜨므로, query
    // accession 이 이 파일을 가리키는 hit 만 남긴다
    // (accession = `<파일 stem>` 또는 foldseek createdb 의 `<파일 stem>_<chain>`).
    std::vector<FoldseekHit> nav_hits;
    if (fs_hits_loaded) {
        std::string qstem = opts.input_files.empty()
                                ? std::string()
                                : fs::path(opts.input_files[0]).stem().string();
        // `x.cif.gz` 는 stem() 이 `x.cif` 를 주므로 한 번 더 벗긴다
        if (qstem.find('.') != std::string::npos) {
            qstem = fs::path(qstem).stem().string();
        }
        for (const FoldseekHit& hit : fs_nav_parser.get_hits()) {
            const bool same_query =
                !qstem.empty() &&
                (hit.query == qstem ||
                 (hit.query.size() > qstem.size() &&
                  hit.query.compare(0, qstem.size(), qstem) == 0 &&
                  hit.query[qstem.size()] == '_'));
            if (same_query) nav_hits.push_back(hit);
        }
        const size_t total = fs_nav_parser.get_hits().size();
        if (nav_hits.empty()) {
            // accession 이름 규칙이 다른 m8(직접 만든 결과 등)에서는 전부 유지한다.
            nav_hits = fs_nav_parser.get_hits();
        } else if (nav_hits.size() != total) {
            std::cerr << "Foldseek hits for query '" << qstem << "': " << nav_hits.size()
                      << " of " << total << " (hits of other queries in " << fs_result
                      << " are skipped)" << std::endl;
        }
    }

    // 구조 파일 query 라도 체인이 여러 개면 `]`/`[` 로 체인을 넘기고 각 체인의 hit 을
    // `N`/`P` 로 볼 수 있게 한다. 체인 전환 때마다 같은 파일을 다른 체인 필터로 다시 읽는다.
    // 사용자가 `-c` 로 체인을 직접 고르거나 target 파일을 추가로 준 경우, MSA/FoldMason 을
    // 쓰는 경우는 아래 단일 체인 경로를 그대로 둔다(각각 별도 씬 구성이 필요하다).
    bool file_query_nav = false;
    const bool query_is_dir = (query_kind == InputKind::StructureDir);
    if (fs_hits_loaded && opts.input_files.size() == 1 && opts.chains_file.empty() &&
        opts.foldmason_file.empty() && opts.msa_file.empty()) {
        std::vector<std::string> chain_ids;                              // m8 등장 순서
        std::map<std::string, std::vector<FoldseekHit>> hits_by_chain;
        for (const FoldseekHit& hit : nav_hits) {
            auto it = hits_by_chain.find(hit.query);
            if (it == hits_by_chain.end()) {
                chain_ids.push_back(hit.query);
                hits_by_chain.emplace(hit.query, std::vector<FoldseekHit>{ hit });
            } else {
                it->second.push_back(hit);
            }
        }
        // 체인이 하나뿐이면 기존 경로가 더 단순하다(씬 구성이 동일하다).
        // 디렉터리 query 는 하나여도 이 경로로 간다 — 파일을 accession 으로 찾아야 하므로.
        if (chain_ids.size() >= 2 || (query_is_dir && !chain_ids.empty())) {
            std::cerr << "Foldseek queries in " << opts.input_files[0] << ": "
                      << chain_ids.size()
                      << (query_is_dir ? " (]/[ switches query chain, N/P walks its hits)"
                                       : " chains (]/[ switches chain, N/P walks its hits)")
                      << std::endl;
            screen.set_fs_db_path(target_dir);
            screen.set_query_nav_from_file(chain_ids, hits_by_chain, opts.input_files[0],
                                           target_db, opts.show_structure, query_is_dir);
            file_query_nav = true;
        }
    }

    if (!file_query_nav) {
    if (fs_hits_loaded) {
        // query(0번 입력): m8 query 컬럼이 이 파일의 체인을 가리키면 그 체인만 남긴다.
        // 그 체인 하나만 로드하므로, **같은 파일의 다른 체인이 낸 hit 은 버려야 한다** —
        // 남겨두면 체인 B 의 정렬을 체인 A 위에 겹치게 된다.
        std::string query_acc;
        if (!opts.input_files.empty()) {
            for (const FoldseekHit& hit : nav_hits) {
                const std::string applied =
                    screen.apply_accession_chain(0, hit.query, opts.input_files[0]);
                if (!applied.empty()) {
                    query_acc = hit.query;
                    std::cerr << "Chain filter from Foldseek query '" << hit.query
                              << "': " << opts.input_files[0] << " -> chain "
                              << applied << std::endl;
                    break;
                }
            }
        }
        if (!query_acc.empty()) {
            std::vector<std::string> other_chains;
            std::vector<FoldseekHit> chain_hits;
            for (const FoldseekHit& hit : nav_hits) {
                if (hit.query == query_acc) {
                    chain_hits.push_back(hit);
                } else if (std::find(other_chains.begin(), other_chains.end(), hit.query)
                           == other_chains.end()) {
                    other_chains.push_back(hit.query);
                }
            }
            if (!other_chains.empty()) {
                std::cerr << "Foldseek hits for '" << query_acc << "': " << chain_hits.size()
                          << " (skipping " << (nav_hits.size() - chain_hits.size())
                          << " hits of";
                for (size_t oi = 0; oi < other_chains.size(); oi++) {
                    std::cerr << (oi == 0 ? " " : ", ") << other_chains[oi];
                }
                std::cerr << ")\n"
                          << "         Pass the Foldseek query DB as the query to walk every"
                          << " chain with ]/[." << std::endl;
            }
            nav_hits = std::move(chain_hits);
        }
        const std::vector<FoldseekHit>& all_hits = nav_hits;
        // target(1번 이후 입력): m8 target 컬럼 기준으로 동일 처리
        for (int ti = 1; ti < (int)opts.input_files.size(); ti++) {
            for (const FoldseekHit& hit : all_hits) {
                const std::string applied =
                    screen.apply_accession_chain(ti, hit.target, opts.input_files[ti]);
                if (!applied.empty()) {
                    std::cerr << "Chain filter from Foldseek target '" << hit.target
                              << "': " << opts.input_files[ti] << " -> chain "
                              << applied << std::endl;
                    break;
                }
            }
        }
    }

    for (int i = 0; i < (int)opts.input_files.size(); i++) {
        screen.set_protein(opts.input_files[i], i, opts.show_structure);
    }
    screen.set_tmatrix();
    screen.normalize_proteins();
    screen.update_total_len_ca();

    // 기능 1: interface 모드일 때 inter-chain interface 계산 (threshold=8.0Å)
    if (opts.mode == "interface") {
        screen.compute_interface_all();
    }

    // 기능 4: align 계열 모드 — foldseek/FoldMason 결과가 없는 경우 (nearest-neighbor fallback)
    if (is_aligned_mode(opts.mode) && fs_result.empty() && opts.foldmason_file.empty()) {
        if (opts.mode == "align-fs") {
            std::cerr << "Error: -m align-fs needs a Foldseek result (-fsr) or a FoldMason "
                      << "alignment (-fm).\n"
                      << "       Use -m align-near to colour by distance instead."
                      << std::endl;
            return false;
        }
        screen.compute_aligned_all();
    }

    // Foldseek DB 직접 읽기 모드 (-fst 가 Foldseek DB 일 때)
    if (!target_db.empty()) {
        screen.open_foldseek_db(target_db);
    }

    // 기능 3: Foldseek hit 탐색 설정 (-fsr 결과가 있을 때)
    if (!fs_result.empty()) {
        if (fs_hits_loaded) {
            if ((int)opts.input_files.size() > 1) {
                // 작업 3-A: 다중 타겟 — m8 hit을 CLI target 파일명 기준으로 필터링
                std::vector<std::string> target_stems;
                for (int i = 1; i < (int)opts.input_files.size(); i++) {
                    fs::path p(opts.input_files[i]);
                    target_stems.push_back(p.stem().string());
                }

                const std::vector<FoldseekHit>& all_hits = nav_hits;
                std::vector<FoldseekHit> filtered_hits;
                for (const FoldseekHit& hit : all_hits) {
                    for (const std::string& stem : target_stems) {
                        if (hit.target.find(stem) != std::string::npos ||
                            stem.find(hit.target) != std::string::npos) {
                            filtered_hits.push_back(hit);
                            break;
                        }
                    }
                }
                screen.set_foldseek_hits(filtered_hits);
                screen.set_fs_db_path(target_dir);
                screen.prepare_foldseek_db(filtered_hits);

                // 작업 3-B: 각 target protein에 매칭 hit의 transform 적용
                for (int ti = 1; ti < (int)opts.input_files.size(); ti++) {
                    fs::path tp(opts.input_files[ti]);
                    std::string tstem = tp.stem().string();

                    for (const FoldseekHit& hit : filtered_hits) {
                        if (hit.target.find(tstem) != std::string::npos ||
                            tstem.find(hit.target) != std::string::npos) {
                            screen.apply_hit_transform(ti, hit);
                            break;
                        }
                    }
                }
            } else {
                // 단일 입력: 기존 동작 유지
                screen.set_foldseek_hits(nav_hits);
                screen.set_fs_db_path(target_dir);
                screen.prepare_foldseek_db(nav_hits);
                screen.load_next_hit(+1);
            }
        }
    }

    // 기능 5: conservation 모드일 때 MSA 파일 로드 및 conservation score 계산
    if (opts.mode == "conservation" && !opts.msa_file.empty()) {
        MSAParser msa_parser;
        if (msa_parser.load(opts.msa_file)) {
            std::vector<float> scores = msa_parser.compute_conservation();
            screen.apply_msa_conservation(0, scores);
        } else {
            std::cerr << "Warning: Failed to load MSA file: " << opts.msa_file << std::endl;
        }
    }

    // 기능 8: FoldMason MSA 기반 superposition + conservation
    if (!opts.foldmason_file.empty()) {
        auto fm_parser = std::make_unique<FoldMasonParser>();
        const std::string& fm_path = opts.foldmason_file;
        bool fm_loaded = false;

        std::string ext;
        {
            size_t dot = fm_path.rfind('.');
            if (dot != std::string::npos) ext = fm_path.substr(dot + 1);
            for (auto& c : ext) c = (char)std::tolower((unsigned char)c);
        }

        if (ext == "json") {
            fm_loaded = fm_parser->load_json(fm_path);
        } else {
            fm_loaded = fm_parser->load_fasta(fm_path);
        }

        if (!fm_loaded || fm_parser->entry_count() == 0) {
            std::cerr << "Warning: Failed to load FoldMason file: " << fm_path << std::endl;
        } else {
            // entry 매칭: 파일명(확장자 제외) 기준으로 탐색, 실패 시 0/1 순서 가정
            int fm_query_idx = 0;
            int fm_target_idx = 1;
            const auto& entries = fm_parser->get_entries();
            if ((int)opts.input_files.size() >= 1) {
                fs::path qp(opts.input_files[0]);
                std::string qstem = qp.stem().string();
                for (int ei = 0; ei < fm_parser->entry_count(); ei++) {
                    if (entries[ei].name.find(qstem) != std::string::npos) {
                        fm_query_idx = ei; break;
                    }
                }
            }
            if ((int)opts.input_files.size() >= 2) {
                fs::path tp(opts.input_files[1]);
                std::string tstem = tp.stem().string();
                for (int ei = 0; ei < fm_parser->entry_count(); ei++) {
                    if (ei == fm_query_idx) continue;
                    if (entries[ei].name.find(tstem) != std::string::npos) {
                        fm_target_idx = ei; break;
                    }
                }
            }

            int total_entries = fm_parser->entry_count();

            // conservation 색상 (단일/다중 구조 모두)
            if (opts.mode == "conservation") {
                std::vector<float> entropy = fm_parser->compute_column_entropy(false);
                auto col_map = fm_parser->build_query_col_map(fm_query_idx);
                std::vector<float> mapped_scores(col_map.size(), 0.0f);
                for (int ri = 0; ri < (int)col_map.size(); ri++) {
                    int col = col_map[ri];
                    if (col < (int)entropy.size()) mapped_scores[ri] = entropy[col];
                }
                screen.apply_msa_conservation(0, mapped_scores);
            }

            // 패널 정보 구성
            FoldMasonInfo fm_info;
            fm_info.valid = true;
            fm_info.entry_count = total_entries;

            // 두 구조 superposition (두 번째 PDB 있을 때만)
            bool do_superposition = ((int)opts.input_files.size() >= 2 && total_entries >= 2);
            if (do_superposition) {
                fm_info.align_method = (opts.mode == "align-near") ? "nearest-nbr"
                                     : (is_aligned_mode(opts.mode) ? "msa-col" : "-");
            } else {
                fm_info.align_method = "-";
            }

            screen.set_foldmason(std::move(fm_parser));
            screen.set_foldmason_panel_info(fm_info);

            if (do_superposition) {
                screen.apply_foldmason_superposition(0, 1, fm_query_idx, fm_target_idx);
            }
        }
    }
    } // end if (!file_query_nav)
    } // end else (non-query-from-DB scene setup)

    if (bench) {
        auto t_load1 = Benchmark::clock::now();
        bm.log("load", -1, Benchmark::ms_since(t_load0, t_load1));

        const std::vector<int> script = {
            'X','Y','Z','A','D','W','S','R','F'
        };

        const int warmup = 200;
        const int events = 2000;

        // Warmup run (not measured)
        bool old_enabled = bm.enabled;
        bm.enabled = false;
        for (int i = 0; i < warmup; i++) {
            screen.draw_screen(opts.no_panel);
            screen.handle_input(script[i % script.size()]);
        }
        bm.enabled = old_enabled;

        // Measured run
        bm.log("bench_begin", -1, 0.0);
        for (int i = 0; i < events; i++) {
            screen.draw_screen(opts.no_panel);
            screen.handle_input(script[i % script.size()]);
        }
        bm.log("bench_end", -1, 0.0);
    } else {
        bool running = true;
        bool needs_redraw = true;
        while (running) {
            if (needs_redraw) {
                screen.draw_screen(opts.no_panel);
            }
            // KEY_MOUSE 이벤트 시 needs_redraw=false: 패널 부분 갱신만 수행
            // 키보드 이벤트 시 needs_redraw=true: 전체 재렌더링
            running = screen.handle_input(needs_redraw);
        }
    }

    Terminal::exit_raw_mode();
    return true;
}

} // namespace structty
