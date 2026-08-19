#include "FoldseekDBReader.hpp"
#include <cstring>
#include <cstdlib>
#include <iostream>

FoldseekDBReader::FoldseekDBReader()
    : ca_data_file_(nullptr), seq_data_file_(nullptr), opened_(false) {}

FoldseekDBReader::~FoldseekDBReader() {
    close();
}

bool FoldseekDBReader::open(const std::string& basePath) {
    close();
    base_path_ = basePath;

    // data file 핸들만 열기 (파싱 없음)
    ca_data_file_ = std::fopen((basePath + "_ca").c_str(), "rb");
    if (!ca_data_file_) {
        std::cerr << "FoldseekDBReader: failed to open " << basePath << "_ca\n";
        return false;
    }

    seq_data_file_ = std::fopen(basePath.c_str(), "rb");
    if (!seq_data_file_) {
        std::cerr << "FoldseekDBReader: failed to open " << basePath << "\n";
        std::fclose(ca_data_file_);
        ca_data_file_ = nullptr;
        return false;
    }

    opened_ = true;
    return true;
}

void FoldseekDBReader::close() {
    if (ca_data_file_) {
        std::fclose(ca_data_file_);
        ca_data_file_ = nullptr;
    }
    if (seq_data_file_) {
        std::fclose(seq_data_file_);
        seq_data_file_ = nullptr;
    }
    lookup_.clear();
    ca_index_.clear();
    seq_index_.clear();
    opened_ = false;
}

bool FoldseekDBReader::is_open() const {
    return opened_;
}

// .lookup 순차 스캔: target_ids에 매칭되는 accession만 저장
bool FoldseekDBReader::scan_lookup(const std::string& path,
                                    const std::unordered_set<std::string>& target_ids) {
    FILE* fp = std::fopen(path.c_str(), "r");
    if (!fp) {
        std::cerr << "FoldseekDBReader: failed to open " << path << "\n";
        return false;
    }

    lookup_.clear();
    lookup_.reserve(target_ids.size());

    char line[4096];
    while (std::fgets(line, sizeof(line), fp)) {
        // 형식: <internal_id>\t<accession>\t<file_number>
        // 첫 번째 탭 찾기
        char* tab1 = std::strchr(line, '\t');
        if (!tab1) continue;

        // internal_id 파싱
        *tab1 = '\0';
        unsigned long internal_id = std::strtoul(line, nullptr, 10);

        // 두 번째 탭 찾기
        char* acc_start = tab1 + 1;
        char* tab2 = std::strchr(acc_start, '\t');
        if (tab2) *tab2 = '\0';

        // newline 제거
        size_t acc_len = std::strlen(acc_start);
        if (acc_len > 0 && acc_start[acc_len - 1] == '\n') {
            acc_start[acc_len - 1] = '\0';
        }

        std::string accession(acc_start);

        if (target_ids.count(accession)) {
            lookup_[accession] = static_cast<unsigned int>(internal_id);

            // 모든 target을 찾으면 조기 종료
            if (lookup_.size() == target_ids.size()) break;
        }
    }

    std::fclose(fp);
    return true;
}

// .index 순차 스캔: target_internal_ids에 매칭되는 id만 저장
bool FoldseekDBReader::scan_index(const std::string& path,
                                   const std::unordered_set<unsigned int>& target_internal_ids,
                                   std::unordered_map<unsigned int, IndexEntry>& index_out) {
    FILE* fp = std::fopen(path.c_str(), "r");
    if (!fp) {
        std::cerr << "FoldseekDBReader: failed to open " << path << "\n";
        return false;
    }

    index_out.clear();
    index_out.reserve(target_internal_ids.size());

    char line[4096];
    while (std::fgets(line, sizeof(line), fp)) {
        // 형식: <internal_id>\t<offset>\t<length>
        char* tab1 = std::strchr(line, '\t');
        if (!tab1) continue;

        *tab1 = '\0';
        unsigned long id = std::strtoul(line, nullptr, 10);
        unsigned int internal_id = static_cast<unsigned int>(id);

        if (!target_internal_ids.count(internal_id)) continue;

        char* off_start = tab1 + 1;
        char* tab2 = std::strchr(off_start, '\t');
        if (!tab2) continue;

        *tab2 = '\0';
        size_t offset = static_cast<size_t>(std::strtoull(off_start, nullptr, 10));
        size_t length = static_cast<size_t>(std::strtoull(tab2 + 1, nullptr, 10));

        index_out[internal_id] = {offset, length};

        // 모든 target을 찾으면 조기 종료
        if (index_out.size() == target_internal_ids.size()) break;
    }

    std::fclose(fp);
    return true;
}

bool FoldseekDBReader::prepare(const std::unordered_set<std::string>& target_ids) {
    if (!opened_) return false;
    if (target_ids.empty()) return true;

    std::cout << "  Scanning Foldseek DB for " << target_ids.size() << " targets...\n";

    // 단계 1: .lookup 스캔
    if (!scan_lookup(base_path_ + ".lookup", target_ids)) return false;
    std::cout << "    lookup: " << lookup_.size() << "/" << target_ids.size() << " found\n";

    if (lookup_.empty()) return true;

    // lookup_에서 얻은 internal_id 집합 구축
    std::unordered_set<unsigned int> internal_ids;
    internal_ids.reserve(lookup_.size());
    for (const auto& pair : lookup_) {
        internal_ids.insert(pair.second);
    }

    // 단계 2: _ca.index 스캔
    if (!scan_index(base_path_ + "_ca.index", internal_ids, ca_index_)) return false;
    std::cout << "    ca_index: " << ca_index_.size() << " entries\n";

    // 단계 3: .index 스캔 (AA 서열)
    if (!scan_index(base_path_ + ".index", internal_ids, seq_index_)) return false;
    std::cout << "    seq_index: " << seq_index_.size() << " entries\n";

    std::cout << "  Foldseek DB ready.\n";
    return true;
}

bool FoldseekDBReader::read_raw(FILE* file, size_t offset, size_t length, std::vector<char>& buf) {
    buf.resize(length);
    if (std::fseek(file, static_cast<long>(offset), SEEK_SET) != 0) return false;
    size_t read_bytes = std::fread(buf.data(), 1, length, file);
    return (read_bytes == length);
}

void FoldseekDBReader::decode_ca(const char* data, size_t entryLen, size_t chainLen,
                                  std::vector<float>& out) {
    out.resize(chainLen * 3);

    // Uncompressed: entryLen >= chainLen * 3 * sizeof(float)
    if (entryLen >= chainLen * 3 * sizeof(float)) {
        std::memcpy(out.data(), data, chainLen * 3 * sizeof(float));
        return;
    }

    // Compressed: Planar delta decoding (Coordinate16.h 로직 이식)
    const char* ptr = data;

    // X 차원
    int32_t start = 0;
    std::memcpy(&start, ptr, sizeof(int32_t));
    ptr += sizeof(int32_t);
    out[0] = start / 1000.0f;
    int32_t diffSum = 0;
    for (size_t i = 1; i < chainLen; ++i) {
        int16_t intDiff = 0;
        std::memcpy(&intDiff, ptr, sizeof(int16_t));
        ptr += sizeof(int16_t);
        diffSum += intDiff;
        out[i] = (start + diffSum) / 1000.0f;
    }

    // Y 차원
    diffSum = 0;
    std::memcpy(&start, ptr, sizeof(int32_t));
    ptr += sizeof(int32_t);
    out[chainLen] = start / 1000.0f;
    for (size_t i = chainLen + 1; i < 2 * chainLen; ++i) {
        int16_t intDiff = 0;
        std::memcpy(&intDiff, ptr, sizeof(int16_t));
        ptr += sizeof(int16_t);
        diffSum += intDiff;
        out[i] = (start + diffSum) / 1000.0f;
    }

    // Z 차원
    diffSum = 0;
    std::memcpy(&start, ptr, sizeof(int32_t));
    ptr += sizeof(int32_t);
    out[2 * chainLen] = start / 1000.0f;
    for (size_t i = 2 * chainLen + 1; i < 3 * chainLen; ++i) {
        int16_t intDiff = 0;
        std::memcpy(&intDiff, ptr, sizeof(int16_t));
        ptr += sizeof(int16_t);
        diffSum += intDiff;
        out[i] = (start + diffSum) / 1000.0f;
    }
}

size_t FoldseekDBReader::read_entry(const std::string& accession,
                                     std::vector<float>& coords,
                                     std::string& aa_seq) {
    if (!opened_) return 0;

    // accession → internal_id
    auto it = lookup_.find(accession);
    if (it == lookup_.end()) return 0;
    unsigned int internal_id = it->second;

    // AA 서열 읽기 (chainLen 결정용)
    auto seq_it = seq_index_.find(internal_id);
    if (seq_it == seq_index_.end()) return 0;
    const IndexEntry& seq_entry = seq_it->second;
    if (seq_entry.length == 0) return 0;

    std::vector<char> seq_buf;
    if (!read_raw(seq_data_file_, seq_entry.offset, seq_entry.length, seq_buf)) return 0;

    // null/newline 제거
    size_t seq_len = seq_entry.length;
    while (seq_len > 0 && (seq_buf[seq_len - 1] == '\0' || seq_buf[seq_len - 1] == '\n')) {
        --seq_len;
    }
    aa_seq.assign(seq_buf.data(), seq_len);
    size_t chainLen = aa_seq.size();
    if (chainLen == 0) return 0;

    // Cα 좌표 읽기
    auto ca_it = ca_index_.find(internal_id);
    if (ca_it == ca_index_.end()) return 0;
    const IndexEntry& ca_entry = ca_it->second;
    if (ca_entry.length == 0) return 0;

    std::vector<char> ca_buf;
    if (!read_raw(ca_data_file_, ca_entry.offset, ca_entry.length, ca_buf)) return 0;

    // 디코딩
    decode_ca(ca_buf.data(), ca_entry.length, chainLen, coords);

    return chainLen;
}
