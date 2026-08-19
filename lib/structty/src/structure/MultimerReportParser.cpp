#include "MultimerReportParser.hpp"

#include <fstream>
#include <sstream>
#include <cstdlib>

namespace {

// 구분자로 문자열 분리 (빈 토큰 포함, 단 마지막 trailing 빈 토큰은 보존).
std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == delim) { out.push_back(cur); cur.clear(); }
        else cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

// 쉼표 분리 float 목록 파싱. expected 개수만큼 out 에 채움(부족하면 false).
bool parse_floats(const std::string& s, float* out, int expected) {
    std::vector<std::string> toks = split(s, ',');
    if ((int)toks.size() < expected) return false;
    for (int i = 0; i < expected; i++) {
        out[i] = std::strtof(toks[i].c_str(), nullptr);
    }
    return true;
}

} // namespace

bool MultimerReportParser::load(const std::string& filepath) {
    std::ifstream in(filepath);
    if (!in.is_open()) return false;

    hits_.clear();
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (!line.empty() && line.back() == '\r') line.pop_back();

        std::vector<std::string> cols = split(line, '\t');
        if (cols.size() < 14) continue;  // 14컬럼 미만은 스킵

        MultimerHit h;
        h.qComplex = cols[0];
        h.tComplex = cols[1];
        h.qChains  = split(cols[2], ',');
        h.tChains  = split(cols[3], ',');
        h.qTMScore = std::strtof(cols[4].c_str(), nullptr);
        h.tTMScore = std::strtof(cols[5].c_str(), nullptr);
        h.has_transform = parse_floats(cols[6], h.U, 9) &&
                          parse_floats(cols[7], h.T, 3);
        h.qComplexCov = std::strtof(cols[8].c_str(), nullptr);
        h.tComplexCov = std::strtof(cols[9].c_str(), nullptr);
        h.qChainTms = cols[10];
        h.tChainTms = cols[11];
        h.interfaceLddt = cols[12];
        h.assId = std::atoi(cols[13].c_str());

        hits_.push_back(std::move(h));
    }
    return !hits_.empty();
}
