#include "SSPredictor.hpp"
#include <algorithm>

std::vector<char> SSPredictor::compute_breaks(const std::vector<Atom>& A) {
    const size_t n = A.size();
    std::vector<char> br(n, 0);
    if (n < 2) return br;
    for (size_t i = 0; i + 1 < n; ++i) {
        if (dist(A[i], A[i+1]) * scale > break_gap) {
            br[i] = 1;
        }
    }
    return br;
}

void SSPredictor::vote(const std::vector<Atom>& A,
                       const std::vector<char>& is_break,
                       std::vector<int>& h_score,
                       std::vector<int>& e_score)
{
    const size_t n = A.size();
    if (n < 4) return;

    auto same_segment = [&](size_t i, size_t j) -> bool {
        if (i > j) std::swap(i, j);
        for (size_t k = i; k < j; ++k) if (is_break[k]) return false;
        return true;
    };

    for (size_t i = 0; i + 3 < n; ++i) {
        if (!same_segment(i, i+3)) continue;

        float d13 = dist(A[i], A[i+2]) * scale;
        float d14 = dist(A[i], A[i+3]) * scale;
        float tau = torsion_deg(A[i], A[i+1], A[i+2], A[i+3]);
        float at = std::fabs(tau);

        int h = 0, e = 0;

        if (d13 >= d13_helix_min && d13 <= d13_helix_max) h++;
        if (d14 >= d14_helix_min && d14 <= d14_helix_max) h++;
        if (at >= tors_helix_abs_min && at <= tors_helix_abs_max) h++;

        if (d13 >= d13_beta_min) e++;
        if (d14 >= d14_beta_min) e++;
        if (at >= tors_beta_abs_min) e++;

        if (h >= 2) {
            h_score[i+1] += 1;
            h_score[i+2] += 1;
        }
        if (e >= 2) {
            e_score[i+1] += 1;
            e_score[i+2] += 1;
        }
    }
}

std::vector<char> SSPredictor::label_from_scores(const std::vector<int>& h_score,
                                                 const std::vector<int>& e_score)
{
    const size_t n = h_score.size();
    std::vector<char> lab(n, 'x'); // default
    for (size_t i = 0; i < n; ++i) {
        int h = h_score[i], e = e_score[i];
        if (h >= vote_threshold && h > e) lab[i] = 'H';
        else if (e >= vote_threshold && e > h) lab[i] = 'S';
        else lab[i] = 'x';
    }
    return lab;
}

void SSPredictor::squash_short_segments(std::vector<char>& lab, char target, int min_len) {
    const size_t n = lab.size();
    size_t i = 0;
    while (i < n) {
        if (lab[i] != target) { ++i; continue; }
        size_t j = i;
        while (j < n && lab[j] == target) ++j;
        int len = static_cast<int>(j - i);
        if (len < min_len) {
            for (size_t k = i; k < j; ++k) lab[k] = 'x';
        }
        i = j;
    }
}

void SSPredictor::smooth_labels(const std::vector<char>& is_break,
                                std::vector<char>& lab)
{
    const size_t n = lab.size();
    if (n == 0) return;

    auto not_across_break = [&](size_t a, size_t b) {
        if (a > b) std::swap(a, b);
        for (size_t k = a; k < b; ++k) if (is_break[k]) return false;
        return true;
    };

    // remove island (size <=smooth_island)
    int K = smooth_island;
    if (K >= 1) {
        for (char t : {'H', 'S'}) {
            size_t i = 0;
            while (i < n) {
                if (lab[i] != t) { ++i; continue; }
                size_t j = i;
                while (j < n && lab[j] == t) ++j;
                int len = static_cast<int>(j - i);

                if (len <= K) {
                    bool across_break = false;
                    if (i > 0 && !not_across_break(i-1, i)) across_break = true;
                    if (j < n && !not_across_break(j-1, j)) across_break = true;
                    if (!across_break) {
                        for (size_t k = i; k < j; ++k) lab[k] = 'x';
                    }
                }
                i = j;
            }
        }
    }

    // min len filter
    squash_short_segments(lab, 'H', helix_min_len);
    squash_short_segments(lab, 'S',  beta_min_len);
}

void SSPredictor::run_chain(std::vector<Atom>& chain_atoms) {
    const size_t n = chain_atoms.size();
    if (n == 0) return;

    // break position
    std::vector<char> is_break = compute_breaks(chain_atoms);

    // vote
    std::vector<int> h_score(n, 0), e_score(n, 0);
    vote(chain_atoms, is_break, h_score, e_score);

    // label -> smoothing
    std::vector<char> lab = label_from_scores(h_score, e_score);
    smooth_labels(is_break, lab);

    // result
    for (size_t i = 0; i < n; ++i) {
        chain_atoms[i].set_structure(lab[i]);
    }
}

void SSPredictor::run(std::map<std::string, std::vector<Atom>>& atoms) {
    std:: cout << "  predict secondary structure\n";    
    for (auto& chain : atoms) {
        auto& chain_atoms = chain.second;
        run_chain(chain_atoms);
    }
}
