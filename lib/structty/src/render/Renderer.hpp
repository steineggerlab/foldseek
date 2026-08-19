#pragma once
#include <vector>
#include <string>
#include <chrono>
#include <cstdint>
#include "RenderPoint.hpp"

struct RenderAtom {
    float x, y, z;
    char  structure;
    float bfactor;
    bool  is_interface;
    bool  is_aligned;
    float conservation_score;
    int   residue_number;
    char  residue_name[4];
    int   chain_id = -1;   // 인터닝된 chain index. was std::string
    int   chain_color_id = -1;  // >=0 이면 chain 모드 색을 이 값으로 정한다(멀티머: 전체 체인 순번)
    int   protein_index;
    float pan_x, pan_y;
};

class Renderer {
public:
    Renderer(int width, int height, const std::string& mode, bool show_structure);

    void set_depth_params(float focal_offset, float zoom_level,
                          float depth_min_z, float depth_max_z);
    void render(const std::vector<RenderAtom>& atoms);

    const std::vector<RenderPoint>& get_pixels() const;
    int get_logical_width()  const;
    int get_logical_height() const;

    // 계측용: 직전 render() 에서 project_and_fill 이 생성한 RenderPoint 총 개수
    // (coil→catmull 점 폭발 가설 검증용).
    size_t get_last_point_count() const { return last_point_count_; }

    // 계측용: 직전 render() 내부 3분할 소요시간 (µs).
    int64_t get_last_us_clear() const { return last_us_clear_; }
    int64_t get_last_us_fill()  const { return last_us_fill_;  }
    int64_t get_last_us_zbuf()  const { return last_us_zbuf_;  }

private:
    enum class Mode { Protein, Chain, Rainbow, Plddt, Interface, Aligned, Conservation };

    // 점 하나를 logical_pixels_ 에 기록할 때 필요한 원자 단위 색/속성 컨텍스트.
    struct PlotStyle {
        const RenderAtom* a;      // 원자 속성(bfactor/interface/aligned/conservation/residue/structure/chain)
        int   protein_idx;
        int   chain_color_idx;    // protein 내 chain 순번 (chain 모드 색)
        float rainbow_frac;       // 0..1, protein 내 원자 진행도 (rainbow 모드 색)
    };

    int         width_, height_;
    std::string mode_;
    Mode        mode_id_ = Mode::Protein;
    bool        show_structure_;
    float       focal_offset_     = 3.0f;
    float       zoom_level_       = 2.0f;
    float       depth_base_min_z_ = 0.0f;
    float       depth_base_max_z_ = 1.0f;

    std::vector<RenderPoint> logical_pixels_;
    size_t                   last_point_count_ = 0;  // 계측용: 생성 시도한 점 수
    std::chrono::steady_clock::time_point t_enter_;  // 계측용: render 진입 시각
    int64_t                  last_us_clear_ = 0;     // 계측용
    int64_t                  last_us_fill_  = 0;     // 계측용
    int64_t                  last_us_zbuf_  = 0;     // 계측용(직접 rasterize 로 0)

    static constexpr float FOV_ = 90.0f;
    static constexpr float PI_  = 3.14159265359f;

    void clear();
    int  compute_depth_band(float z) const;
    int  color_from_style(const PlotStyle& s, int band) const;
    // 점 하나를 logical_pixels_ 에 z-test 후 직접 기록 (통과 시에만 color 계산).
    void plot(int x, int y, float depth, const PlotStyle& s);
    void draw_line_impl(int x1, int x2, int y1, int y2,
                        float z1, float z2,
                        const PlotStyle& s, int half);
    void project_and_fill(const std::vector<RenderAtom>& atoms);
};
