#pragma once
#include <vector>
#include <map>
#include <cmath>
#include <cstddef>
#include "Atom.hpp"  

class SSPredictor {
public:
    float scale = 1.0f;
    float d13_helix_min = 4.8f, d13_helix_max = 6.0f; 
    float d14_helix_min = 4.8f, d14_helix_max = 6.0f; 
    float d13_beta_min  = 6.4f;                       
    float d14_beta_min  = 8.5f;                     

    float tors_helix_abs_min = 35.0f, tors_helix_abs_max = 75.0f;  
    float tors_beta_abs_min  = 110.0f;                              

    int   helix_min_len = 4;
    int   beta_min_len  = 3;

    float break_gap = 4.8f;


    int   vote_threshold = 2;

    int   smooth_island = 1;

    void run(std::map<std::string, std::vector<Atom>>& atoms);

    void run_chain(std::vector<Atom>& chain_atoms);

    void set_scale(float scale_) { scale = scale_; }

private:
    static inline float dist(const Atom& a, const Atom& b) {
        float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    static inline void cross(const float a[3], const float b[3], float out[3]) {
        out[0] = a[1]*b[2] - a[2]*b[1];
        out[1] = a[2]*b[0] - a[0]*b[2];
        out[2] = a[0]*b[1] - a[1]*b[0];
    }

    static inline float dot(const float a[3], const float b[3]) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }

    static inline void sub(const Atom& p, const Atom& q, float out[3]) {
        out[0] = p.x - q.x; out[1] = p.y - q.y; out[2] = p.z - q.z;
    }

    static inline void normalize(float v[3]) {
        float n = std::sqrt(dot(v, v));
        if (n > 1e-8f) { v[0]/=n; v[1]/=n; v[2]/=n; }
    }

    static float torsion_deg(const Atom& a, const Atom& b, const Atom& c, const Atom& d) {
        float b1[3], b2[3], b3[3];
        sub(b, a, b1); sub(c, b, b2); sub(d, c, b3);

        float n1[3], n2[3];
        cross(b1, b2, n1);
        cross(b2, b3, n2);

        float b2n[3] = { b2[0], b2[1], b2[2] };
        normalize(b2n);

        float m1[3];
        cross(n1, b2n, m1);

        float x = dot(n1, n2);
        float y = dot(m1, n2);

        float angle = std::atan2(y, x) * 180.0f / 3.14159265358979323846f;
        return angle;
    }

    // break mask: if nearest residue is far, break. true
    std::vector<char> compute_breaks(const std::vector<Atom>& A);

    // vote based on the distance/torsion
    void vote(const std::vector<Atom>& A,
              const std::vector<char>& is_break,
              std::vector<int>& h_score,
              std::vector<int>& e_score);

    // score -> label
    std::vector<char> label_from_scores(const std::vector<int>& h_score,
                                        const std::vector<int>& e_score);

    // smoothing (remove island, min len filter)
    void smooth_labels(const std::vector<char>& is_break,
                       std::vector<char>& lab);

    // sequenctial region length check
    static void squash_short_segments(std::vector<char>& lab, char target, int min_len);
};
