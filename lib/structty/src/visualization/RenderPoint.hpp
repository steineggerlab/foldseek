#pragma once
#include <limits>
#include <string>

struct RenderPoint {
    int         x = 0;
    int         y = 0;
    float       depth = std::numeric_limits<float>::infinity();
    char        pixel = ' ';
    int         color_id = 0;
    int         chainID = -1;   // 인터닝된 chain index (Screen::chain_name 로 역참조). was std::string
    char        structure = 0;

    // 기능 1: interface region
    bool        is_interface = false;

    // 기능 4: 정렬 구조
    bool        is_aligned = false;

    // 기능 2: pLDDT
    float       bfactor = 0.0f;

    // 기능 5: conservation
    float       conservation_score = -1.0f;

    // 기능 6: 잔기 정보
    // std::string 금지: RenderPoint는 매 프레임 대량 생성·소멸하므로 heap allocation 방지
    int         residue_number = -1;
    char        residue_name[4] = {};  // 최대 3글자(GLU, ALA 등) + null terminator

    // 기능 7: depth fog (braille 모드 원근감)
    int         depth_band = 0;  // 0=near(vivid), 1=mid(normal), 2=far(dim)
};
