#include "Renderer.hpp"
#include "Palette.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

Renderer::Renderer(int width, int height, const std::string& mode, bool show_structure)
    : width_(width), height_(height), mode_(mode), show_structure_(show_structure)
{
    logical_pixels_.resize(width_ * 2 * height_ * 4);
    // per-point 문자열 비교를 피하기 위해 모드를 enum 으로 1회 캐시.
    if      (mode_ == "chain")        mode_id_ = Mode::Chain;
    else if (mode_ == "rainbow")      mode_id_ = Mode::Rainbow;
    else if (mode_ == "plddt")        mode_id_ = Mode::Plddt;
    else if (mode_ == "interface")    mode_id_ = Mode::Interface;
    else if (mode_ == "align" || mode_ == "align-fs" || mode_ == "align-near")
                                      mode_id_ = Mode::Aligned;
    else if (mode_ == "conservation") mode_id_ = Mode::Conservation;
    else                              mode_id_ = Mode::Protein;
}

void Renderer::set_depth_params(float focal_offset, float zoom_level,
                                float depth_min_z, float depth_max_z) {
    focal_offset_     = focal_offset;
    zoom_level_       = zoom_level;
    depth_base_min_z_ = depth_min_z;
    depth_base_max_z_ = depth_max_z;
}

void Renderer::render(const std::vector<RenderAtom>& atoms) {
    t_enter_ = std::chrono::steady_clock::now();  // 계측: clear 이전 시각
    clear();
    auto t_clear = std::chrono::steady_clock::now();

    last_point_count_ = 0;
    project_and_fill(atoms);  // 생성 즉시 logical_pixels_ 에 z-test 기록 (중간 벡터 없음)
    auto t_fill = std::chrono::steady_clock::now();

    // 계측: clear / fill. zbuffer 는 fill 에 흡수되어 0.
    last_us_clear_ = std::chrono::duration_cast<std::chrono::microseconds>(t_clear - t_enter_).count();
    last_us_fill_  = std::chrono::duration_cast<std::chrono::microseconds>(t_fill  - t_clear).count();
    last_us_zbuf_  = 0;
}

const std::vector<RenderPoint>& Renderer::get_pixels() const {
    return logical_pixels_;
}

int Renderer::get_logical_width()  const { return width_  * 2; }
int Renderer::get_logical_height() const { return height_ * 4; }

void Renderer::clear() {
    logical_pixels_.assign(width_ * 2 * height_ * 4, RenderPoint());
}

int Renderer::compute_depth_band(float z) const {
    float range = depth_base_max_z_ - depth_base_min_z_;
    if (range <= 0.0f) return 1;
    float t = (z - depth_base_min_z_) / range;
    if      (t < 0.25f) return 0;
    else if (t < 0.60f) return 1;
    else                return 2;
}


int Renderer::color_from_style(const PlotStyle& s, int band) const {
    const RenderAtom& a = *s.a;

    // show_structure 는 protein 모드에서만 SS 색으로 override (원본 assign_colors_impl 과 동일).
    // helix/sheet 는 전용 near/far 색(47-50), coil 은 protein near/far 배열을 재사용해 3밴드.
    if (show_structure_ && mode_id_ == Mode::Protein) {
        if (a.structure == 'H') {
            if (band == 0) return 47;
            if (band == 1) return 41;
            return 48;
        }
        if (a.structure == 'S') {
            if (band == 0) return 49;
            if (band == 1) return 42;
            return 50;
        }
        const int coil_idx = s.protein_idx % 9;
        if (band == 0) return coil_idx + 120;  // PROTEIN_NEAR
        if (band == 1) return coil_idx + 11;   // PROTEIN_DIM
        return coil_idx + 85;                  // PROTEIN_COIL_FAR (DIM 보다 한 단계 어둡다)
    }

    switch (mode_id_) {
        case Mode::Protein: {
            int idx = s.protein_idx % 9;
            if (band == 0) return idx + 120;
            if (band == 1) return idx + 1;
            return idx + 200;
        }
        case Mode::Chain: {
            int ci = (a.chain_color_id >= 0 ? a.chain_color_id
                                            : s.protein_idx * 10 + s.chain_color_idx) % 15;
            if (band == 0) return 130 + ci;
            if (band == 1) return 21 + ci;
            return 145 + ci;
        }
        case Mode::Rainbow: {
            int ci = (int)(s.rainbow_frac * 20.0f);
            if (ci < 0) ci = 0;
            if (ci > 19) ci = 19;
            if (band == 0) return ci + 160;
            if (band == 1) return ci + 51;
            return ci + 180;
        }
        case Mode::Plddt: {
            float plddt = a.bfactor;
            int base = (plddt >= 90) ? 0 : (plddt >= 70) ? 1 : (plddt >= 50) ? 2 : 3;
            if (band == 0) return 209 + base;
            if (band == 1) return 71 + base;
            return 213 + base;
        }
        case Mode::Interface: {
            if (band == 0) return a.is_interface ? 237 : 238;
            if (band == 1) return a.is_interface ? 43 : 44;
            return a.is_interface ? 239 : 240;
        }
        case Mode::Aligned: {
            const int idx = s.protein_idx % 9;
            if (a.is_aligned) {
                if (band == 0) return idx + 241;  // aligned near
                if (band == 1) return idx + 101;  // aligned mid (bright)
                return idx + 251;                 // aligned far
            }
            if (band == 0) return 111;  // non-aligned near
            if (band == 1) return 110;  // non-aligned mid
            return 250;                 // non-aligned far
        }
        case Mode::Conservation: {
            float score = a.conservation_score;
            if (score < 0) {
                // 미채점 잔기: 중립 회색 3단 (10 / 36 / 20)
                if (band == 0) return 10;
                if (band == 1) return 36;
                return 20;
            }
            int idx = (int)(score * 9.0f);
            if (idx < 0) idx = 0;
            if (idx > 9) idx = 9;
            if (band == 0) return 217 + idx;
            if (band == 1) return 75 + idx;
            return 227 + idx;
        }
    }
    return 0;
}

void Renderer::plot(int x, int y, float depth, const PlotStyle& s) {
    ++last_point_count_;  // 계측: 생성 시도한 점 수
    const int logical_w = width_  * 2;
    const int logical_h = height_ * 4;
    if (x < 0 || x >= logical_w || y < 0 || y >= logical_h) return;

    RenderPoint& p = logical_pixels_[y * logical_w + x];
    if (depth >= p.depth) return;  // 뒤쪽/동일 → 버림 (frontmost 유지 = 기존 zbuffer 시맨틱)

    int band = compute_depth_band(depth);  // 통과한 점만 band/color 계산 (오버드로 절약)
    const RenderAtom& a = *s.a;
    p.x = x; p.y = y; p.depth = depth; p.pixel = ' ';
    p.color_id           = color_from_style(s, band);
    p.chainID            = a.chain_id;
    p.structure          = a.structure;
    p.depth_band         = band;
    p.bfactor            = a.bfactor;
    p.is_interface       = a.is_interface;
    p.is_aligned         = a.is_aligned;
    p.conservation_score = a.conservation_score;
    p.residue_number     = a.residue_number;
    std::strncpy(p.residue_name, a.residue_name, 3);
    p.residue_name[3]    = '\0';
}

void Renderer::draw_line_impl(int x1, int x2, int y1, int y2,
                              float z1, float z2,
                              const PlotStyle& s, int half) {
    int dx = x2 - x1;
    int dy = y2 - y1;
    float dz = z2 - z1;

    int steps = std::max(std::abs(dx), std::abs(dy));
    if (steps == 0) steps = 1;

    float xInc = (float)dx / steps;
    float yInc = (float)dy / steps;
    float zInc = dz / steps;

    float fx = (float)x1, fy = (float)y1, fz = z1;
    for (int i = 0; i <= steps; ++i) {
        int ix = (int)fx, iy = (int)fy;
        for (int oy = -half; oy <= half; oy++) {
            for (int ox = -half; ox <= half; ox++) {
                if (ox != 0 && oy != 0) continue;
                plot(ix + ox, iy + oy, fz, s);  // plot 이 경계검사 + z-test 수행
            }
        }
        fx += xInc; fy += yInc; fz += zInc;
    }
}

void Renderer::project_and_fill(const std::vector<RenderAtom>& atoms) {
    const float nearPlane = 0.05f;
    float fovRads = 1.0f / std::tan((FOV_ / zoom_level_) * 0.5f / 180.0f * PI_);

    const int logical_w = width_  * 2;
    const int logical_h = height_ * 4;

    size_t n = atoms.size();
    size_t protein_start = 0;

    while (protein_start < n) {
        int cur_protein = atoms[protein_start].protein_index;
        size_t protein_end = protein_start;
        while (protein_end < n && atoms[protein_end].protein_index == cur_protein)
            ++protein_end;
        int protein_atoms = (int)(protein_end - protein_start);

        int chain_color_idx = 0;  // protein 내 chain 순번 (chain 모드)
        size_t chain_start = protein_start;
        while (chain_start < protein_end) {
            int cur_chain = atoms[chain_start].chain_id;
            size_t chain_end = chain_start;
            while (chain_end < protein_end && atoms[chain_end].chain_id == cur_chain)
                ++chain_end;

            int num_atoms = (int)(chain_end - chain_start);
            int prevScreenX = -1, prevScreenY = -1;
            float prevZ = -1.0f;

            for (int i = 0; i < num_atoms; ++i) {
                const RenderAtom& a = atoms[chain_start + i];
                float x = a.x;
                float y = a.y;
                // 카메라 Z 부호 통일: 회전행렬(Protein::set_rotate)은 표준 오른손좌표계(+Z가 시청자 쪽)를
                // 가정하므로, 카메라 공간에서는 로컬 Z가 클수록(=시청자에 가까울수록) near가 되어야 함.
                float z = -a.z + focal_offset_;

                if (z < nearPlane) {
                    prevScreenX = prevScreenY = -1;
                    prevZ = -1.0f;
                    continue;
                }

                char structure = a.structure;
                float projectedX = (x / z) * fovRads + a.pan_x;
                float projectedY = (y / z) * fovRads + a.pan_y;
                int screenX = (int)((projectedX + 1.0f) * 0.5f * logical_w);
                int screenY = (int)((1.0f - projectedY) * 0.5f * logical_h);

                // rainbow: protein 내 원자 진행도 (원본 point-index 기반 → atom-index 기반, deviation)
                int local = (int)((chain_start + i) - protein_start);
                float frac = (protein_atoms > 1) ? (float)local / (protein_atoms - 1) : 0.0f;
                PlotStyle style{ &a, cur_protein, chain_color_idx, frac };

                if (prevScreenX != -1 && prevScreenY != -1) {
                    if (structure != 'H' && structure != 'S') {
                        const RenderAtom& P0 = atoms[chain_start + std::max(0, i - 2)];
                        const RenderAtom& P1 = atoms[chain_start + (i - 1)];
                        const RenderAtom& P2 = atoms[chain_start + i];
                        const RenderAtom& P3 = atoms[chain_start + std::min(num_atoms - 1, i + 1)];
                        int cr_prevX = prevScreenX, cr_prevY = prevScreenY;
                        float cr_prevZ = prevZ;
                        for (int cr = 1; cr <= 4; ++cr) {
                            float t = cr * 0.25f;
                            float t2 = t * t, t3 = t2 * t;
                            float cx = 0.5f * ((-P0.x + 3*P1.x - 3*P2.x + P3.x)*t3
                                             + ( 2*P0.x - 5*P1.x + 4*P2.x - P3.x)*t2
                                             + (-P0.x + P2.x)*t + 2*P1.x);
                            float cy = 0.5f * ((-P0.y + 3*P1.y - 3*P2.y + P3.y)*t3
                                             + ( 2*P0.y - 5*P1.y + 4*P2.y - P3.y)*t2
                                             + (-P0.y + P2.y)*t + 2*P1.y);
                            // 카메라 Z 부호 통일: 위 단일-원자 투영과 동일하게 블렌드된 로컬 Z 전체를
                            // 음수화한 뒤 focal_offset_ 를 더함 (focal_offset_ 자체는 부호 반전 대상 아님).
                            float cz = -(0.5f * ((-P0.z + 3*P1.z - 3*P2.z + P3.z)*t3
                                             + ( 2*P0.z - 5*P1.z + 4*P2.z - P3.z)*t2
                                             + (-P0.z + P2.z)*t + 2*P1.z)) + focal_offset_;
                            if (cz < nearPlane) break;
                            float cpX = (cx / cz) * fovRads + a.pan_x;
                            float cpY = (cy / cz) * fovRads + a.pan_y;
                            int cX = (int)((cpX + 1.0f) * 0.5f * logical_w);
                            int cY = (int)((1.0f - cpY) * 0.5f * logical_h);
                            draw_line_impl(cr_prevX, cX, cr_prevY, cY, cr_prevZ, cz, style, 1);
                            cr_prevX = cX; cr_prevY = cY; cr_prevZ = cz;
                        }
                    } else {
                        draw_line_impl(prevScreenX, screenX, prevScreenY, screenY, prevZ, z, style, 1);
                    }
                }

                if (screenX >= 0 && screenX < logical_w && screenY >= 0 && screenY < logical_h) {
                    if (structure == 'x') {
                        plot(screenX, screenY,     z, style);
                        plot(screenX, screenY - 1, z, style);
                        plot(screenX, screenY + 1, z, style);
                    } else {
                        for (int oy = -1; oy <= 1; oy++)
                            for (int ox = -1; ox <= 1; ox++) {
                                if (ox != 0 && oy != 0) continue;
                                plot(screenX + ox, screenY + oy, z, style);
                            }
                    }
                }

                prevScreenX = screenX;
                prevScreenY = screenY;
                prevZ = z;
            }

            chain_start = chain_end;
            ++chain_color_idx;
        }

        protein_start = protein_end;
    }
}
