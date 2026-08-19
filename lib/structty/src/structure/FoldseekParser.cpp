#include "FoldseekParser.hpp"

#include <fstream>
#include <sstream>

// 탭으로 구분된 행을 컬럼 벡터로 분리
static std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> cols;
    std::istringstream ss(line);
    std::string col;
    while (std::getline(ss, col, '\t')) {
        cols.push_back(col);
    }
    return cols;
}

// 쉼표로 구분된 float 문자열 파싱
static std::vector<float> parse_csv_floats(const std::string& s) {
    std::vector<float> result;
    std::istringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        try { result.push_back(std::stof(tok)); } catch (...) {}
    }
    return result;
}

static float safe_stof(const std::string& s, float default_val = 0.0f) {
    try { return std::stof(s); } catch (...) { return default_val; }
}

static int safe_stoi(const std::string& s, int default_val = 0) {
    try { return std::stoi(s); } catch (...) { return default_val; }
}

// target 컬럼에서 첫 번째 공백 이전 토큰만 추출
// PDB m8 결과에 "3a0d_A Crystal Structure of..." 형태로 description이 붙는 경우 대응
static std::string extract_target_id(const std::string& raw) {
    size_t pos = raw.find(' ');
    if (pos == std::string::npos) return raw;
    return raw.substr(0, pos);
}

bool FoldseekParser::parse_line(const std::vector<std::string>& cols,
                                FoldseekHit& hit, int fmt) {
    if ((int)cols.size() < fmt) return false;

    hit.query  = cols[0];
    hit.target = extract_target_id(cols[1]);
    hit.fident = safe_stof(cols[2]);
    hit.alnlen = safe_stoi(cols[3]);
    // cols[4]=mismatch, cols[5]=gapopen
    hit.qstart = safe_stoi(cols[6], 1);
    hit.qend   = safe_stoi(cols[7], 0);
    hit.tstart = safe_stoi(cols[8], 1);
    hit.tend   = safe_stoi(cols[9], 0);

    if (fmt == 12) {
        hit.evalue = safe_stof(cols[10], 0.0f);
        hit.bits   = safe_stof(cols[11], 0.0f);
        hit.has_transform = false;
        hit.has_aln       = false;
        return true;
    }

    if (fmt == 21) {
        // alis 포맷: col[10]=prob, col[11]=evalue, col[12]=bits
        //            col[13]=qlength, col[14]=tlength
        //            col[15]=qaln, col[16]=taln
        //            col[17]=alns (comma-separated CA coords, query frame 기준)
        //            col[18]=tseq, col[19]=taxid, col[20]=taxname
        hit.prob    = safe_stof(cols[10], -1.0f);
        hit.evalue  = safe_stof(cols[11], 0.0f);
        hit.bits    = safe_stof(cols[12], 0.0f);
        hit.qlength = safe_stoi(cols[13], 0);
        hit.tlength = safe_stoi(cols[14], 0);
        hit.qaln    = cols[15];
        hit.taln    = cols[16];
        hit.alns    = parse_csv_floats(cols[17]);
        hit.tseq    = cols[18];
        hit.taxid   = cols[19];
        // taxname은 마지막 컬럼 — 탭이 없는 경우 이미 전체가 들어있음
        hit.taxname = cols[20];
        hit.has_aln       = (!hit.qaln.empty() && !hit.taln.empty());
        hit.has_transform = false;  // Kabsch로 load_next_hit()에서 채워짐
        hit.is_alis_format = true;
        return true;
    }

    // fmt 17 or 29: cols[12]=lddt, [13]=qtmscore, [14]=ttmscore
    hit.lddt      = safe_stof(cols[12], -1.0f);
    hit.qtmscore  = safe_stof(cols[13], -1.0f);
    hit.ttmscore  = safe_stof(cols[14], -1.0f);
    hit.evalue    = safe_stof(cols[10], 0.0f);
    hit.bits      = safe_stof(cols[11], 0.0f);

    if (fmt == 17) {
        hit.qaln = cols[15];
        hit.taln = cols[16];
        hit.has_aln       = true;
        hit.has_transform = false;
        return true;
    }

    // fmt == 29: cols[15..23]=U (9 values), cols[24..26]=T (3 values)
    //            cols[27]=qaln, cols[28]=taln
    for (int i = 0; i < 9; ++i) {
        hit.U[i] = safe_stof(cols[15 + i]);
    }
    for (int i = 0; i < 3; ++i) {
        hit.T[i] = safe_stof(cols[24 + i]);
    }
    hit.qaln = cols[27];
    hit.taln = cols[28];
    hit.has_transform = true;
    hit.has_aln       = true;
    return true;
}

bool FoldseekParser::load(const std::string& filepath) {
    hits.clear();

    std::ifstream ifs(filepath);
    if (!ifs.is_open()) {
        return false;
    }

    int fmt = 0;  // detected format (12, 17, 21, or 29)
    std::string line;

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::vector<std::string> cols = split_tabs(line);
        int ncols = (int)cols.size();

        // 첫 번째 데이터 행에서 포맷 감지
        if (fmt == 0) {
            if (ncols == 12) {
                fmt = 12;
            } else if (ncols == 17) {
                fmt = 17;
            } else if (ncols == 21) {
                fmt = 21;
            } else if (ncols == 29) {
                fmt = 29;
            } else {
                // taxname에 탭이 포함될 수 있으므로 21+N 컬럼도 alis로 처리
                if (ncols > 21) {
                    fmt = 21;
                    // 마지막 컬럼들을 taxname으로 합치기
                    // cols는 이미 split됐으므로 idx 20 이후를 합친다
                    while ((int)cols.size() > 21) {
                        cols[20] += "\t" + cols[21];
                        cols.erase(cols.begin() + 21);
                    }
                    ncols = 21;
                } else {
                    return false;
                }
            }
        }

        // alis 포맷에서 taxname에 탭이 있을 경우 처리
        if (fmt == 21 && (int)cols.size() > 21) {
            while ((int)cols.size() > 21) {
                cols[20] += "\t" + cols[21];
                cols.erase(cols.begin() + 21);
            }
        }

        FoldseekHit hit;
        if (parse_line(cols, hit, fmt)) {
            hits.push_back(std::move(hit));
        }
    }

    return !hits.empty();
}
