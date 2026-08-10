// LoLAlign homology probability: 39-feature MLP that emits P(homolog) for a LoLAlign hit.
// Originally referred to as the "v20" classifier (homology_classifier_cv0_v20.pt).
//
// Features + constants are a 1:1 port of nn_classifier/phase13_afdb_features_predict.py.
// Adjacency semantics match nn_classifier/phase12_build_query_adjacency.py:
//   undirected edge (i,j) iff CA distance < 20 A; degree = #neighbors (excl. self).
// Coordinates are foldseek's flat per-axis CA arrays (x[0..N), y[0..N), z[0..N)).
//
// Set the env var LOLALIGN_PROB_DEBUG=1 to print the 39-vector + prob per call.
#ifndef LOLALIGN_PROBABILITY_H
#define LOLALIGN_PROBABILITY_H

#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "LoLAlignProbabilityWeights.h"

class LoLAlignProbability {
public:
    static constexpr float  RADIUS = 20.0f;
    static constexpr float  ANF_C  = 5.0f;
    static constexpr int    NB    = 7;   // ANF histogram bins
    static constexpr int    NSEP6 = 6;   // seq-sep bins
    static constexpr int    NSEP5 = 5;   // merged seq-sep bins (tamax)

    LoLAlignProbability() : qLen_(0) {}

    // Precompute everything that depends only on the query (degree, per-residue 6-bin
    // seq-sep histogram, edge list). Call once per query, before the target loop.
    void prepareQuery(int queryLen, const float* qx, const float* qy, const float* qz) {
        qLen_ = queryLen;
        buildAdjacency(queryLen, qx, qy, qz, qDeg_, qEdgeS_, qEdgeD_, qEdgeBin_);
        buildResHist6(queryLen, qEdgeS_, qEdgeD_, qEdgeBin_, qHist6_);
    }

    // Compute P(homolog) for one alignment.
    //  tx,ty,tz    : target CA coords (flat per-axis, len targetLen)
    //  qAln,tAln   : 0-based aligned residue indices (len nAln), query/target
    //  evalue,bits : raw LoL scores (res.eval, res.score)
    //  lddt        : LDDTCalculator avgLddtScore for this alignment
    //  alnlen      : M+I+D positions; gapopen: #gap runs; numGaps: I+D positions
    float predict(int targetLen, const float* tx, const float* ty, const float* tz,
                  const int* qAln, const int* tAln, int nAln,
                  double evalue, double bits, double lddt,
                  int alnlen, int gapopen, int numGaps) {
        std::vector<int> tDeg, tEdgeS, tEdgeD, tEdgeBin;
        buildAdjacency(targetLen, tx, ty, tz, tDeg, tEdgeS, tEdgeD, tEdgeBin);
        std::vector<float> tHist6;
        buildResHist6(targetLen, tEdgeS, tEdgeD, tEdgeBin, tHist6);

        float feat[39];
        // --- BASE scalars (0..7) -------------------------------------------------
        feat[0] = (float)evalue;
        feat[1] = (float)std::log(std::max(evalue, 1e-300));
        feat[2] = (float)bits;
        feat[3] = (float)(sign(bits) * std::log1p(std::fabs(bits)));
        feat[4] = (float)lddt;
        feat[5] = (float)std::log((double)alnlen + 1.0);
        feat[6] = (float)((double)gapopen / std::max(alnlen, 1));
        feat[7] = (float)((double)numGaps / std::max(alnlen, 1));
        // --- q_anf (8..14), t_anf (15..21) --------------------------------------
        float h[NB];
        anfHist(qLen_, qDeg_, qEdgeS_, qEdgeD_, qAln, nAln, h);
        for (int i = 0; i < NB; ++i) feat[8 + i] = h[i];
        anfHist(targetLen, tDeg, tEdgeS, tEdgeD, tAln, nAln, h);
        for (int i = 0; i < NB; ++i) feat[15 + i] = h[i];
        // --- tamax (22..26) -----------------------------------------------------
        float tam[NSEP5];
        tamaxFeat(qHist6_, qLen_, tHist6, targetLen, qAln, tAln, nAln, tam);
        for (int i = 0; i < NSEP5; ++i) feat[22 + i] = tam[i];
        // --- qsep (27..32) ------------------------------------------------------
        float qs[NSEP6];
        qsepFeat(qHist6_, qLen_, qAln, nAln, qs);
        for (int i = 0; i < NSEP6; ++i) feat[27 + i] = qs[i];
        // --- qfrac (33..38) -----------------------------------------------------
        float qf[NSEP6];
        qfracFeat(qLen_, qEdgeS_, qEdgeD_, qEdgeBin_, qAln, nAln, qf);
        for (int i = 0; i < NSEP6; ++i) feat[33 + i] = qf[i];

        return forward(feat);
    }

private:
    int qLen_;
    std::vector<int>   qDeg_, qEdgeS_, qEdgeD_, qEdgeBin_;
    std::vector<float> qHist6_;   // qLen_ * NSEP6, row-major

    static inline double sign(double x) { return (x > 0) - (x < 0); }

    static inline int sepBin(int s) {            // = #{e in {5,13,31,51,101} : e <= |s|}
        int a = s < 0 ? -s : s, b = 0;
        if (a >= 5)   ++b;
        if (a >= 13)  ++b;
        if (a >= 31)  ++b;
        if (a >= 51)  ++b;
        if (a >= 101) ++b;
        return b;
    }

    // Undirected edge list (i<j, CA dist<20) + degree, from flat per-axis coords.
    static void buildAdjacency(int N, const float* x, const float* y, const float* z,
                               std::vector<int>& deg, std::vector<int>& es,
                               std::vector<int>& ed, std::vector<int>& eb) {
        deg.assign(N, 0); es.clear(); ed.clear(); eb.clear();
        const float R2 = RADIUS * RADIUS;
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                float dx = x[i] - x[j], dy = y[i] - y[j], dz = z[i] - z[j];
                if (dx * dx + dy * dy + dz * dz < R2) {   // dist < 20
                    ++deg[i]; ++deg[j];
                    es.push_back(i); ed.push_back(j); eb.push_back(sepBin(i - j));
                }
            }
        }
    }

    // Per-residue 6-bin seq-sep histogram, row-normalized. out: N*NSEP6.
    static void buildResHist6(int N, const std::vector<int>& es, const std::vector<int>& ed,
                              const std::vector<int>& eb, std::vector<float>& out) {
        out.assign((size_t)N * NSEP6, 0.0f);
        for (size_t e = 0; e < es.size(); ++e) {
            int b = eb[e];
            out[(size_t)es[e] * NSEP6 + b] += 1.0f;
            out[(size_t)ed[e] * NSEP6 + b] += 1.0f;
        }
        for (int i = 0; i < N; ++i) {
            float* row = &out[(size_t)i * NSEP6];
            float tot = 0.0f;
            for (int b = 0; b < NSEP6; ++b) tot += row[b];
            if (tot > 0.0f) for (int b = 0; b < NSEP6; ++b) row[b] /= tot;
        }
    }

    // ANF: histogram of (aligned-neighbors / (degree + C)) over aligned residues w/ deg>0.
    static void anfHist(int N, const std::vector<int>& deg,
                        const std::vector<int>& es, const std::vector<int>& ed,
                        const int* aln, int nAln, float* out) {
        for (int i = 0; i < NB; ++i) out[i] = 0.0f;
        if (deg.empty() || es.empty()) return;
        std::vector<char> amask(N, 0);
        for (int k = 0; k < nAln; ++k) { int a = aln[k]; if (a >= 0 && a < N) amask[a] = 1; }
        std::vector<int> nAlnNb(N, 0);
        for (size_t e = 0; e < es.size(); ++e) {
            int s = es[e], d = ed[e];
            nAlnNb[s] += amask[d];
            nAlnNb[d] += amask[s];
        }
        static const float EDGES[NB + 1] =
            {0.0f, 0.25f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0001f};
        float tot = 0.0f;
        for (int i = 0; i < N; ++i) {
            if (!amask[i] || deg[i] <= 0) continue;
            float frac = (float)nAlnNb[i] / ((float)deg[i] + ANF_C);
            int b = NB - 1;                          // last bin closed on the right
            for (int k = 0; k < NB; ++k) {
                if (frac >= EDGES[k] && frac < EDGES[k + 1]) { b = k; break; }
            }
            out[b] += 1.0f; tot += 1.0f;
        }
        if (tot > 0.0f) for (int i = 0; i < NB; ++i) out[i] /= tot;
    }

    static inline void to5(const float* h6, float* h5) {
        h5[0] = h6[0]; h5[1] = h6[1]; h5[2] = h6[2];
        h5[3] = h6[3] + h6[4]; h5[4] = h6[5];
    }

    // tamax: per-bin abs-max of (to5(qHist6[qi]) - to5(tHist6[ti])) over aligned pairs.
    static void tamaxFeat(const std::vector<float>& qH6, int qN,
                          const std::vector<float>& tH6, int tN,
                          const int* qAln, const int* tAln, int nAln, float* out) {
        for (int i = 0; i < NSEP5; ++i) out[i] = 0.0f;
        bool any = false;
        float q5[NSEP5], t5[NSEP5];
        for (int k = 0; k < nAln; ++k) {
            int qi = qAln[k], ti = tAln[k];
            if (qi < 0 || qi >= qN || ti < 0 || ti >= tN) continue;
            to5(&qH6[(size_t)qi * NSEP6], q5);
            to5(&tH6[(size_t)ti * NSEP6], t5);
            for (int b = 0; b < NSEP5; ++b) {
                float a = std::fabs(q5[b] - t5[b]);
                if (a > out[b]) out[b] = a;
            }
            any = true;
        }
        if (!any) for (int i = 0; i < NSEP5; ++i) out[i] = 0.0f;
    }

    // qsep: mean of query per-residue 6-bin histogram over aligned query residues.
    static void qsepFeat(const std::vector<float>& qH6, int qN,
                         const int* qAln, int nAln, float* out) {
        for (int i = 0; i < NSEP6; ++i) out[i] = 0.0f;
        int cnt = 0;
        for (int k = 0; k < nAln; ++k) {
            int qi = qAln[k];
            if (qi < 0 || qi >= qN) continue;
            const float* row = &qH6[(size_t)qi * NSEP6];
            for (int b = 0; b < NSEP6; ++b) out[b] += row[b];
            ++cnt;
        }
        if (cnt > 0) for (int b = 0; b < NSEP6; ++b) out[b] /= (float)cnt;
    }

    // qfrac: per seq-sep bin, fraction of query contacts that are aligned-aligned.
    static void qfracFeat(int qN, const std::vector<int>& es, const std::vector<int>& ed,
                          const std::vector<int>& eb, const int* qAln, int nAln, float* out) {
        for (int i = 0; i < NSEP6; ++i) out[i] = 0.0f;
        if (es.empty()) return;
        std::vector<char> mask(qN, 0);
        for (int k = 0; k < nAln; ++k) { int a = qAln[k]; if (a >= 0 && a < qN) mask[a] = 1; }
        double all_cnt[NSEP6] = {0}, aln_cnt[NSEP6] = {0};
        for (size_t e = 0; e < es.size(); ++e) {
            int b = eb[e]; char ms = mask[es[e]], md = mask[ed[e]];
            all_cnt[b] += (double)ms + (double)md;
            if (ms && md) aln_cnt[b] += 2.0;
        }
        for (int b = 0; b < NSEP6; ++b)
            out[b] = all_cnt[b] > 0.0 ? (float)(aln_cnt[b] / all_cnt[b]) : 0.0f;
    }

    // (x-mean)/std -> 39->32 relu -> 32->16 relu -> 16->1 -> sigmoid.
    static float forward(const float* feat) {
        float x[lolprobw::N_FEAT];
        for (int i = 0; i < lolprobw::N_FEAT; ++i)
            x[i] = (feat[i] - lolprobw::FEAT_MEAN[i]) / lolprobw::FEAT_STD[i];

        float h0[lolprobw::H0];
        for (int o = 0; o < lolprobw::H0; ++o) {
            float a = lolprobw::B0[o];
            const float* w = &lolprobw::W0[(size_t)o * lolprobw::N_FEAT];
            for (int i = 0; i < lolprobw::N_FEAT; ++i) a += w[i] * x[i];
            h0[o] = a > 0.0f ? a : 0.0f;
        }
        float h1[lolprobw::H1];
        for (int o = 0; o < lolprobw::H1; ++o) {
            float a = lolprobw::B1[o];
            const float* w = &lolprobw::W1[(size_t)o * lolprobw::H0];
            for (int i = 0; i < lolprobw::H0; ++i) a += w[i] * h0[i];
            h1[o] = a > 0.0f ? a : 0.0f;
        }
        float logit = lolprobw::B2;
        for (int i = 0; i < lolprobw::H1; ++i) logit += lolprobw::W2[i] * h1[i];
        float prob = 1.0f / (1.0f + std::exp(-logit));

        if (getenv("LOLALIGN_PROB_DEBUG")) {
            fprintf(stderr, "LOLALIGN_PROB_FEAT");
            for (int i = 0; i < lolprobw::N_FEAT; ++i) fprintf(stderr, "\t%.6g", feat[i]);
            fprintf(stderr, "\tprob=%.6f\n", prob);
        }
        return prob;
    }
};

#endif
