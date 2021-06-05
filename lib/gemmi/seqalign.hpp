// Simple pairwise sequence alignment.
//
// Code in this file is based on and derived from files ksw2_gg.c and ksw2.h
// from https://github.com/lh3/ksw2, which is under the MIT license.
// The original code, written by Heng Li, has more features and has more
// efficient variants that use SSE instructions.

#ifndef GEMMI_SEQALIGN_HPP_
#define GEMMI_SEQALIGN_HPP_

#include <cstdint>
#include <algorithm> // for reverse
#include <map>
#include <string>
#include <vector>

namespace gemmi {

struct AlignmentScoring {
  int match = 1;
  int mismatch = -1;
  int gapo = -1;
  int gape = -1;
  std::vector<std::int8_t> score_matrix;
  std::vector<std::string> matrix_encoding;
};

struct AlignmentResult {
  struct Item {
    std::uint32_t value;
    char op() const { return "MID"[value & 0xf]; }
    std::uint32_t len() const { return value >> 4; }
  };
  int score = 0;
  int match_count = 0;
  std::string match_string;
  std::vector<Item> cigar;

  std::string cigar_str() const {
    std::string s;
    for (Item item : cigar) {
      s += std::to_string(item.len());
      s += item.op();
    }
    return s;
  }

  // 1=query, 2=target, other=shorter
  std::size_t input_length(int which) const {
    std::size_t counters[3] = {0, 0, 0};
    for (Item item : cigar)
      counters[item.value & 0xf] += item.len();
    if (which == 1 || which == 2)
      return counters[0] + counters[which];
    return counters[0] + std::min(counters[1], counters[2]);
  }
  double calculate_identity(int which=0) const {
    return 100. * match_count / input_length(which);
  }

  // In the backtrack matrix, value p[] has the following structure:
  //   bit 0-2: which type gets the max - 0 for H, 1 for E, 2 for F
  //   bit 3/0x08: 1 if a continuation on the E state
  //   bit 4/0x10: 1 if a continuation on the F state
  void backtrack_to_cigar(const std::uint8_t *p, int i, int j) {
    i--;
    int j0 = j--;
    int state = 0;
    while (i >= 0 && j >= 0) {
      // at the beginning of the loop, _state_ tells us which state to check
      // if requesting the H state, find state one maximizes it.
      uint32_t tmp = p[(std::size_t)i * j0 + j];
      if (state == 0 || (tmp & (1 << (state + 2))) == 0)
        state = tmp & 7;
      if (state == 0) { // match
        push_cigar(0, 1);
        --i;
        --j;
      } else if (state == 1) { // deletion
        push_cigar(2, 1);
        --i;
      } else { // insertion
        push_cigar(1, 1);
        --j;
      }
    }
    if (i >= 0)
      push_cigar(2, i + 1); // first deletion
    else if (j >= 0)
      push_cigar(1, j + 1); // first insertion
    std::reverse(cigar.begin(), cigar.end());
  }


  void count_matches(const std::vector<std::uint8_t>& query,
                     const std::vector<std::uint8_t>& target) {
    match_count = 0;
    size_t pos1 = 0, pos2 = 0;
    for (Item item : cigar)
      if (item.op() == 'M') {
        for (uint32_t i = 0; i < item.len(); ++i)
          if (query[pos1++] == target[pos2++]) {
            ++match_count;
            match_string += '|';
          } else {
            match_string += '.';
          }
      } else if (item.op() == 'I') {
        pos1 += item.len();
        match_string.append(item.len(), ' ');
      } else /*item.op() == 'D'*/ {
        pos2 += item.len();
        match_string.append(item.len(), ' ');
      }
  }

  std::string add_gaps(const std::string& s, unsigned which) const {
    std::string out;
    size_t pos = 0;
    for (Item item : cigar) {
      bool show = (item.value & 0xf) == 0 || (item.value & 0xf) == which;
      for (uint32_t i = 0; i < item.len(); ++i)
        out += show ? s.at(pos++) : '-';
    }
    return out;
  }

  std::string formatted(const std::string& a, const std::string& b) const {
    std::string r;
    r.reserve((match_string.size() + 1) * 3);
    r += add_gaps(a, 1);
    r += '\n';
    r += match_string;
    r += '\n';
    r += add_gaps(b, 2);
    r += '\n';
    return r;
  }

  // op: 0=match/mismatch, 1=insertion, 2=deletion
  void push_cigar(std::uint32_t op, int len) {
    if (cigar.empty() || op != (cigar.back().value & 0xf))
      cigar.push_back({len<<4 | op});
    else
      cigar.back().value += len<<4;
  }
};

// All values in query and target must be less then m.
// free_gapo marks positions in target where gap opening is free.
inline
AlignmentResult align_sequences(const std::vector<std::uint8_t>& query,
                                const std::vector<std::uint8_t>& target,
                                const std::vector<bool>& free_gapo,
                                std::uint8_t m,
                                const AlignmentScoring& scoring) {
  // generate the query profile
  std::int8_t *query_profile = new std::int8_t[query.size() * m];
  {
    std::uint32_t mat_size = (std::uint32_t) scoring.matrix_encoding.size();
    std::int32_t i = 0;
    for (std::uint8_t k = 0; k < m; ++k)
      for (std::uint8_t q : query)
        if (k < mat_size && q < mat_size)
          query_profile[i++] = scoring.score_matrix[k * m + q];
        else
          query_profile[i++] = (k == q ? scoring.match : scoring.mismatch);
  }

  struct eh_t { std::int32_t h, e; };
  eh_t *eh = new eh_t[query.size() + 1];
  std::int32_t gape = scoring.gape;
  std::int32_t gapoe = scoring.gapo + gape;

  // fill the first row
  {
    std::int32_t gap0 = !free_gapo.empty() && free_gapo[0] ? gape : gapoe;
    eh[0].h = 0;
    eh[0].e = gap0 + gapoe;
    for (std::int32_t j = 1; j <= (std::int32_t)query.size(); ++j) {
      eh[j].h = gap0 + gape * (j - 1);
      eh[j].e = gap0 + gapoe + gape * j;
    }
  }

  // backtrack matrix; in each cell: f<<4|e<<2|h
  std::uint8_t *z = new std::uint8_t[query.size() * target.size()];
  // DP loop
  for (std::int32_t i = 0; i < (std::int32_t)target.size(); ++i) {
    std::uint8_t target_item = target[i];
    std::int8_t *scores = &query_profile[target_item * query.size()];
    std::uint8_t *zi = &z[i * query.size()];
    std::int32_t h1 = gapoe + gape * i;
    std::int32_t f = gapoe + gapoe + gape * i;
    std::int32_t gapx = i+1 < (std::int32_t)free_gapo.size() && free_gapo[i+1]
                        ? gape : gapoe;
    for (std::size_t j = 0; j < query.size(); ++j) {
      // At the beginning of the loop:
      //  eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
      // Cells are computed in the following order:
      //   H(i,j)   = max{H(i-1,j-1) + S(i,j), E(i,j), F(i,j)}
      //   E(i+1,j) = max{H(i,j)+gapo, E(i,j)} + gape
      //   F(i,j+1) = max{H(i,j)+gapo, F(i,j)} + gape
      eh_t *p = &eh[j];
      std::int32_t h = p->h;
      std::int32_t e = p->e;
      p->h = h1;
      h += scores[j];
      std::uint8_t direction = 0;
      if (h < e) {
        direction = 1;  // deletion
        h = e;
      }
      if (h <= f) {
        direction = 2;  // insertion
        h = f;
      }
      h1 = h;

      h += gapoe;
      e += gape;
      if (e > h)
        direction |= 0x08;
      else
        e = h;

      h = h1 + gapx;
      p->e = e;
      f += gape;
      if (f > h)
        direction |= 0x10;
      else
        f = h;

      // z[i,j] keeps h for the current cell and e/f for the next cell
      zi[j] = direction;
    }
    eh[query.size()].h = h1;
    eh[query.size()].e = -0x40000000; // -infinity
  }

  AlignmentResult result;
  result.score = eh[query.size()].h;
  delete [] query_profile;
  delete [] eh;
  result.backtrack_to_cigar(z, (int)target.size(), (int)query.size());
  delete [] z;
  result.count_matches(query, target);
  return result;
}

inline
AlignmentResult align_string_sequences(const std::vector<std::string>& query,
                                       const std::vector<std::string>& target,
                                       const std::vector<bool>& free_gapo,
                                       const AlignmentScoring& scoring) {
  std::map<std::string, std::uint8_t> encoding;
  for (const std::string& res_name : scoring.matrix_encoding)
    encoding.emplace(res_name, (std::uint8_t)encoding.size());
  for (const std::string& s : query)
    encoding.emplace(s, (std::uint8_t)encoding.size());
  for (const std::string& s : target)
    encoding.emplace(s, (std::uint8_t)encoding.size());
  if (encoding.size() > 255)
    return AlignmentResult();
  std::vector<std::uint8_t> encoded_query(query.size());
  for (size_t i = 0; i != query.size(); ++i)
    encoded_query[i] = encoding.at(query[i]);
  std::vector<std::uint8_t> encoded_target(target.size());
  for (size_t i = 0; i != target.size(); ++i)
    encoded_target[i] = encoding.at(target[i]);
  return align_sequences(encoded_query, encoded_target,
                         free_gapo, (std::uint8_t)encoding.size(), scoring);
}

inline AlignmentScoring prepare_blosum62_scoring() {
  AlignmentScoring s;
  s.match = 1;
  s.mismatch = -4;
  s.gapo = -10;  // BLAST uses BLOSUM-62 with gap cost (10,1)
  s.gape = -1;
  s.score_matrix = {
    4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,
   -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,
   -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,
   -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,
    0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,
   -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,
   -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,
    0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,
   -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,
   -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,
   -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,
   -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,
   -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,
   -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,
   -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,
    1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,
    0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,
   -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,
   -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,
    0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,
  };
  s.matrix_encoding = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
  };
  return s;
}

} // namespace gemmi
#endif
