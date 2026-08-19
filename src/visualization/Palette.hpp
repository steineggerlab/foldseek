#pragma once
#include <array>
#include <cstdint>
#include <string>

struct RGBA { uint8_t r,g,b,a; };
static_assert(sizeof(RGBA) == 4, "RGBA must be 4 bytes");

namespace Palettes {
    inline const std::array<int, 20> RAINBOW = {
        196, 202, 208, 214, 220, 226, 190, 154, 118, 82,
        49,  51,  45,  39,  33,  27,  21,  93,  129, 201
    };

    // 9 named protein colors (vivid), pairs 1-9
    inline const std::array<int, 9> PROTEIN_COLORS = {
        100,  //  1  olive        #878700
         80,  //  2  turquoise    #5fd7d7
         27,  //  3  navy         #005fff
        129,  //  4  purple       #af00ff
        213,  //  5  pink         #ff87ff
        209,  //  6  coral        #ff875f
        130,  //  7  brown        #af5f00
        214,  //  8  orange       #ffaf00
        160,  //  9  red          #d70000
    };

    // Dimmed counterparts for coil atoms in protein+-s mode, pairs 11-19
    inline const std::array<int, 9> PROTEIN_DIM_COLORS = {
         58,  // 11  olive dim     #5f5f00
         30,  // 12  turquoise dim #008787
         18,  // 13  navy dim      #000087
         54,  // 14  purple dim    #5f0087
        132,  // 15  pink dim      #af5f87
        131,  // 16  coral dim     #af5f5f
         94,  // 17  brown dim     #875f00
        172,  // 18  orange dim    #d78700
         88,  // 19  red dim       #870000
    };

    // Aligned 강조용 밝은 protein colors, pairs 101-109
    inline const std::array<int, 9> PROTEIN_BRIGHT_COLORS = {
        148,  // 101  olive bright    #afd700
        123,  // 102  turquoise bright #87d7ff
         69,  // 103  navy bright     #5f87ff
        171,  // 104  purple bright   #d75fff
        219,  // 105  pink bright     #ffafff
        216,  // 106  coral bright    #ffaf87
        178,  // 107  brown bright    #d7af00
        221,  // 108  orange bright   #ffd75f
        203,  // 109  red bright      #ff5f5f
    };

    // Non-aligned 통일 회색, pair 110
    inline constexpr int ALIGNED_NONALIGNED_DIM = 240;  // xterm-240 #585858

    // Depth fog: protein near (bright) colors, pairs 120-128
    // 기존 대비 +80~120 RGB 차이 — 매우 밝은 변형
    inline const std::array<int, 9> PROTEIN_NEAR_COLORS = {
        184,  // 120  olive near      #d7d700 (215,215,0)
        159,  // 121  turquoise near  #afffff (175,255,255)
        105,  // 122  navy near       #8787ff (135,135,255)
        177,  // 123  purple near     #d787ff (215,135,255)
        225,  // 124  pink near       #ffd7ff (255,215,255)
        223,  // 125  coral near      #ffd7af (255,215,175)
        220,  // 126  brown near      #ffd700 (255,215,0)
        229,  // 127  orange near     #ffffaf (255,255,175)
        210,  // 128  red near        #ff8787 (255,135,135)
    };

    // Depth fog: protein far colors (dark but hue-retaining), pairs 200-208
    inline const std::array<int, 9> PROTEIN_FAR_COLORS = {
         58,  // 200  olive far       #5f5f00 (95,95,0)
         30,  // 201  turquoise far   #008787 (0,135,135)
         18,  // 202  navy far        #000087 (0,0,135)
         54,  // 203  purple far      #5f0087 (95,0,135)
         96,  // 204  pink far        #875f87 (135,95,135)
        130,  // 205  coral far       #af5f00 (175,95,0)
         94,  // 206  brown far       #875f00 (135,95,0)
        136,  // 207  orange far      #af8700 (175,135,0)
         88,  // 208  red far         #870000 (135,0,0)
    };

    // Depth fog: chain near colors (brighter), pairs 130-144
    inline const std::array<int, 15> CHAIN_NEAR_COLORS = {
        184,  // 130  olive near      #d7d700
        159,  // 131  turquoise near  #afffff
        105,  // 132  navy near       #8787ff
        177,  // 133  purple near     #d787ff
        225,  // 134  pink near       #ffd7ff
        223,  // 135  coral near      #ffd7af
        220,  // 136  brown near      #ffd700
        229,  // 137  orange near     #ffffaf
        210,  // 138  red near        #ff8787
         87,  // 139  teal near       #5fffff
        192,  // 140  lime near       #d7ff87
        213,  // 141  magenta near    #ff87ff
        228,  // 142  gold near       #ffff87
        189,  // 143  lavender near   #d7d7ff
        217,  // 144  salmon near     #ffafaf
    };

    // Depth fog: chain far colors (dark but hue-retaining), pairs 145-159
    inline const std::array<int, 15> CHAIN_FAR_COLORS = {
         58,  // 145  olive far       #5f5f00
         30,  // 146  turquoise far   #008787
         18,  // 147  navy far        #000087
         54,  // 148  purple far      #5f0087
         96,  // 149  pink far        #875f87
        130,  // 150  coral far       #af5f00
         94,  // 151  brown far       #875f00
        136,  // 152  orange far      #af8700
         88,  // 153  red far         #870000
         23,  // 154  teal far        #005f5f
         64,  // 155  lime far        #5f8700
         90,  // 156  magenta far     #870087
        136,  // 157  gold far        #af8700
         60,  // 158  lavender far    #5f5f87
        131,  // 159  salmon far      #af5f5f
    };

    // Depth fog: rainbow near colors (brighter), pairs 160-179
    inline const std::array<int, 20> RAINBOW_NEAR = {
        210, 216, 222, 228, 229, 230, 192, 157, 157, 120,
         85,  87,  81,  75, 105,  63,  57,  99, 177, 213
    };

    // Depth fog: rainbow far colors (dark hue-retaining), pairs 180-199
    inline const std::array<int, 20> RAINBOW_FAR = {
        124, 130, 136, 172, 178, 142, 106,  70,  34,  28,
         23,  24,  24,  18,  18,  18,  54,  54,  90, 127
    };

    // Depth fog: pLDDT near (brighter, hue-retaining), pairs 209-212
    // 26->69 #5f87ff, 81->117 #87d7ff, 220->229 #ffffaf, 209->216 #ffaf87
    inline const std::array<int, 4> PLDDT_NEAR = {69, 117, 229, 216};
    // Depth fog: pLDDT far (dark, hue-retaining), pairs 213-216
    // 26->24 #005f87, 81->31 #0087af, 220->58 #5f5f00, 209->95 #875f5f
    inline const std::array<int, 4> PLDDT_FAR = {24, 31, 58, 95};

    // Depth fog: conservation near (brighter), pairs 217-226
    inline const std::array<int, 10> CONSERVATION_NEAR = {57, 63, 69, 75, 85, 157, 192, 229, 221, 210};
    // Depth fog: conservation far (dark hue-retaining), pairs 227-236
    inline const std::array<int, 10> CONSERVATION_FAR = {18, 18, 19, 24, 28, 30, 64, 58, 94, 88};

    // Depth fog: interface near/far, pairs 237-238 (near), 239-240 (far)
    inline constexpr int INTERFACE_NEAR_COLOR     = 213;  // xterm-213 #ff87ff (brighter magenta)
    inline constexpr int INTERFACE_DIM_NEAR_COLOR = 243;  // xterm-243 #767676 (lighter gray)
    inline constexpr int INTERFACE_FAR_COLOR      =  90;  // xterm-90  #870087 (dark magenta)
    inline constexpr int INTERFACE_DIM_FAR_COLOR  = 236;  // xterm-236 #303030 (dark gray, intentional)

    // Depth fog: aligned near/far, pairs 241-249 (near), 250 (far dim)
    // aligned bright near = even brighter protein colors
    inline const std::array<int, 9> ALIGNED_BRIGHT_NEAR = {
        190,  // 241  olive aligned near   #d7ff00
        159,  // 242  turquoise aligned near #afffff
        111,  // 243  navy aligned near    #8787ff
        183,  // 244  purple aligned near  #d7afff
        231,  // 245  pink aligned near    #ffffff
        224,  // 246  coral aligned near   #ffd7d7
        227,  // 247  brown aligned near   #ffff5f
        230,  // 248  orange aligned near  #ffffd7
        211,  // 249  red aligned near     #ff875f
    };
    // aligned dim far = darker gray
    inline constexpr int ALIGNED_NONALIGNED_FAR = 235;  // xterm-235 #262626

    // Depth fog: aligned far, pairs 251-259.
    // PROTEIN_BRIGHT(101-109) 보다 어둡고 PROTEIN_FAR(200-208) 보다는 밝게 —
    // 멀리 있어도 "정렬됨" 이 구분되도록 색상은 유지하고 명도만 낮춘다.
    inline const std::array<int, 9> ALIGNED_BRIGHT_FAR = {
        106,  // 251  olive aligned far     #87af00
         73,  // 252  turquoise aligned far #5fafaf
         26,  // 253  navy aligned far      #005fd7
         98,  // 254  purple aligned far    #875fd7
        176,  // 255  pink aligned far      #d787d7
        173,  // 256  coral aligned far     #d7875f
        137,  // 257  brown aligned far     #af875f
        172,  // 258  orange aligned far    #d78700
        167,  // 259  red aligned far       #d75f5f
    };

    // Non-aligned near (240 mid / 235 far 보다 밝은 회색), pair 111
    inline constexpr int ALIGNED_NONALIGNED_NEAR = 245;  // xterm-245 #8a8a8a

    // Depth fog: secondary structure colors for `-s` + protein mode.
    // mid 는 기존 pair 41(helix, xterm-226) / 42(sheet, xterm-51) 를 그대로 쓴다.
    inline constexpr int SS_HELIX_NEAR = 228;  //  47  #ffff87
    inline constexpr int SS_HELIX_FAR  = 142;  //  48  #afaf00
    inline constexpr int SS_SHEET_NEAR =  87;  //  49  #5fffff
    inline constexpr int SS_SHEET_FAR  =  37;  //  50  #00afaf

    // Depth fog: `-s` protein 모드의 coil far, pairs 85-93.
    // PROTEIN_DIM(11-19, coil mid) 보다 한 단계 더 어둡다. DIM 과 PROTEIN_FAR(200-208)
    // 은 6개 슬롯에서 값이 같아 coil 이 2단계로 뭉쳤기 때문에 전용 배열을 둔다.
    inline const std::array<int, 9> PROTEIN_COIL_FAR = {
        236,  //  85  olive coil far      #303030 (큐브에 더 어두운 올리브가 없어 중립 dark)
         23,  //  86  turquoise coil far  #005f5f
         17,  //  87  navy coil far       #00005f
         53,  //  88  purple coil far     #5f005f
         96,  //  89  pink coil far       #875f87
         95,  //  90  coral coil far      #875f5f
         58,  //  91  brown coil far      #5f5f00
        130,  //  92  orange coil far     #af5f00
         52,  //  93  red coil far        #5f0000
    };

    // Conservation 미채점(score < 0) 잔기: 중립 회색 3단, pairs 10 / 36 / 20
    inline constexpr int CONSERVATION_UNSCORED_NEAR = 249;  //  10  #b2b2b2
    inline constexpr int CONSERVATION_UNSCORED_MID  = 244;  //  36  #808080
    inline constexpr int CONSERVATION_UNSCORED_FAR  = 238;  //  20  #444444

    // Interface region colors, pairs 43-44
    // pair 43: interface residue (bright magenta, xterm-201)
    // pair 44: non-interface dim  (olive dim, xterm-58)
    inline constexpr int INTERFACE_COLOR     = 201;
    inline constexpr int INTERFACE_DIM_COLOR = 58;

    // Aligned region colors, pairs 45-46
    // pair 45: aligned residue (bright green, xterm-46)
    // pair 46: non-aligned dim (olive dim, xterm-58)
    inline constexpr int ALIGNED_COLOR     = 46;
    inline constexpr int ALIGNED_DIM_COLOR = 58;

    // pLDDT colors: pairs 71-74 (matches AlphaFold DB confidence legend)
    // 71: Very High (>=90) xterm-26  #005fd7 (dark/royal blue,  official #0053D6)
    // 72: Confident  (70-90) xterm-81  #5fd7ff (light sky blue,  official #65CBF3)
    // 73: Low        (50-70) xterm-220 #ffd700 (golden yellow,   official #FFDB13)
    // 74: Very Low   (<50)   xterm-209 #ff875f (salmon orange,   official #FF7D45)
    inline const std::array<int, 4> PLDDT_COLORS = {26, 81, 220, 209};

    // Conservation gradient: pairs 75-84 (blue→red, 10 steps, variable→conserved)
    inline const std::array<int, 10> CONSERVATION_COLORS = {21, 27, 33, 39, 49, 118, 190, 226, 214, 196};

    // 15 chain colors (9 given + 6 extended), pairs 21-35
    inline const std::array<int, 15> CHAIN_COLORS = {
        100,  // 21  olive
         80,  // 22  turquoise
         27,  // 23  navy
        129,  // 24  purple
        213,  // 25  pink
        209,  // 26  coral
        130,  // 27  brown
        214,  // 28  orange
        160,  // 29  red
         37,  // 30  teal         #00afaf
        118,  // 31  lime         #87ff00
        201,  // 32  magenta      #ff00ff
        220,  // 33  gold         #ffd700
        183,  // 34  lavender     #d7afff
        210,  // 35  salmon       #ff8787
    };

    inline constexpr RGBA ID2RGBA[256] = {
        RGBA{0, 0, 0, 255}, // 0
        RGBA{128, 0, 0, 255}, // 1
        RGBA{0, 128, 0, 255}, // 2
        RGBA{128, 128, 0, 255}, // 3
        RGBA{0, 0, 128, 255}, // 4
        RGBA{128, 0, 128, 255}, // 5
        RGBA{0, 128, 128, 255}, // 6
        RGBA{192, 192, 192, 255}, // 7
        RGBA{128, 128, 128, 255}, // 8
        RGBA{255, 0, 0, 255}, // 9
        RGBA{0, 255, 0, 255}, // 10
        RGBA{255, 255, 0, 255}, // 11
        RGBA{0, 0, 255, 255}, // 12
        RGBA{255, 0, 255, 255}, // 13
        RGBA{0, 255, 255, 255}, // 14
        RGBA{255, 255, 255, 255}, // 15
        RGBA{0, 0, 0, 255}, // 16
        RGBA{0, 0, 95, 255}, // 17
        RGBA{0, 0, 135, 255}, // 18
        RGBA{0, 0, 175, 255}, // 19
        RGBA{0, 0, 215, 255}, // 20
        RGBA{0, 0, 255, 255}, // 21
        RGBA{0, 95, 0, 255}, // 22
        RGBA{0, 95, 95, 255}, // 23
        RGBA{0, 95, 135, 255}, // 24
        RGBA{0, 95, 175, 255}, // 25
        RGBA{0, 95, 215, 255}, // 26
        RGBA{0, 95, 255, 255}, // 27
        RGBA{0, 135, 0, 255}, // 28
        RGBA{0, 135, 95, 255}, // 29
        RGBA{0, 135, 135, 255}, // 30
        RGBA{0, 135, 175, 255}, // 31
        RGBA{0, 135, 215, 255}, // 32
        RGBA{0, 135, 255, 255}, // 33
        RGBA{0, 175, 0, 255}, // 34
        RGBA{0, 175, 95, 255}, // 35
        RGBA{0, 175, 135, 255}, // 36
        RGBA{0, 175, 175, 255}, // 37
        RGBA{0, 175, 215, 255}, // 38
        RGBA{0, 175, 255, 255}, // 39
        RGBA{0, 215, 0, 255}, // 40
        RGBA{0, 215, 95, 255}, // 41
        RGBA{0, 215, 135, 255}, // 42
        RGBA{0, 215, 175, 255}, // 43
        RGBA{0, 215, 215, 255}, // 44
        RGBA{0, 215, 255, 255}, // 45
        RGBA{0, 255, 0, 255}, // 46
        RGBA{0, 255, 95, 255}, // 47
        RGBA{0, 255, 135, 255}, // 48
        RGBA{0, 255, 175, 255}, // 49
        RGBA{0, 255, 215, 255}, // 50
        RGBA{0, 255, 255, 255}, // 51
        RGBA{95, 0, 0, 255}, // 52
        RGBA{95, 0, 95, 255}, // 53
        RGBA{95, 0, 135, 255}, // 54
        RGBA{95, 0, 175, 255}, // 55
        RGBA{95, 0, 215, 255}, // 56
        RGBA{95, 0, 255, 255}, // 57
        RGBA{95, 95, 0, 255}, // 58
        RGBA{95, 95, 95, 255}, // 59
        RGBA{95, 95, 135, 255}, // 60
        RGBA{95, 95, 175, 255}, // 61
        RGBA{95, 95, 215, 255}, // 62
        RGBA{95, 95, 255, 255}, // 63
        RGBA{95, 135, 0, 255}, // 64
        RGBA{95, 135, 95, 255}, // 65
        RGBA{95, 135, 135, 255}, // 66
        RGBA{95, 135, 175, 255}, // 67
        RGBA{95, 135, 215, 255}, // 68
        RGBA{95, 135, 255, 255}, // 69
        RGBA{95, 175, 0, 255}, // 70
        RGBA{95, 175, 95, 255}, // 71
        RGBA{95, 175, 135, 255}, // 72
        RGBA{95, 175, 175, 255}, // 73
        RGBA{95, 175, 215, 255}, // 74
        RGBA{95, 175, 255, 255}, // 75
        RGBA{95, 215, 0, 255}, // 76
        RGBA{95, 215, 95, 255}, // 77
        RGBA{95, 215, 135, 255}, // 78
        RGBA{95, 215, 175, 255}, // 79
        RGBA{95, 215, 215, 255}, // 80
        RGBA{95, 215, 255, 255}, // 81
        RGBA{95, 255, 0, 255}, // 82
        RGBA{95, 255, 95, 255}, // 83
        RGBA{95, 255, 135, 255}, // 84
        RGBA{95, 255, 175, 255}, // 85
        RGBA{95, 255, 215, 255}, // 86
        RGBA{95, 255, 255, 255}, // 87
        RGBA{135, 0, 0, 255}, // 88
        RGBA{135, 0, 95, 255}, // 89
        RGBA{135, 0, 135, 255}, // 90
        RGBA{135, 0, 175, 255}, // 91
        RGBA{135, 0, 215, 255}, // 92
        RGBA{135, 0, 255, 255}, // 93
        RGBA{135, 95, 0, 255}, // 94
        RGBA{135, 95, 95, 255}, // 95
        RGBA{135, 95, 135, 255}, // 96
        RGBA{135, 95, 175, 255}, // 97
        RGBA{135, 95, 215, 255}, // 98
        RGBA{135, 95, 255, 255}, // 99
        RGBA{135, 135, 0, 255}, // 100
        RGBA{135, 135, 95, 255}, // 101
        RGBA{135, 135, 135, 255}, // 102
        RGBA{135, 135, 175, 255}, // 103
        RGBA{135, 135, 215, 255}, // 104
        RGBA{135, 135, 255, 255}, // 105
        RGBA{135, 175, 0, 255}, // 106
        RGBA{135, 175, 95, 255}, // 107
        RGBA{135, 175, 135, 255}, // 108
        RGBA{135, 175, 175, 255}, // 109
        RGBA{135, 175, 215, 255}, // 110
        RGBA{135, 175, 255, 255}, // 111
        RGBA{135, 215, 0, 255}, // 112
        RGBA{135, 215, 95, 255}, // 113
        RGBA{135, 215, 135, 255}, // 114
        RGBA{135, 215, 175, 255}, // 115
        RGBA{135, 215, 215, 255}, // 116
        RGBA{135, 215, 255, 255}, // 117
        RGBA{135, 255, 0, 255}, // 118
        RGBA{135, 255, 95, 255}, // 119
        RGBA{135, 255, 135, 255}, // 120
        RGBA{135, 255, 175, 255}, // 121
        RGBA{135, 255, 215, 255}, // 122
        RGBA{135, 255, 255, 255}, // 123
        RGBA{175, 0, 0, 255}, // 124
        RGBA{175, 0, 95, 255}, // 125
        RGBA{175, 0, 135, 255}, // 126
        RGBA{175, 0, 175, 255}, // 127
        RGBA{175, 0, 215, 255}, // 128
        RGBA{175, 0, 255, 255}, // 129
        RGBA{175, 95, 0, 255}, // 130
        RGBA{175, 95, 95, 255}, // 131
        RGBA{175, 95, 135, 255}, // 132
        RGBA{175, 95, 175, 255}, // 133
        RGBA{175, 95, 215, 255}, // 134
        RGBA{175, 95, 255, 255}, // 135
        RGBA{175, 135, 0, 255}, // 136
        RGBA{175, 135, 95, 255}, // 137
        RGBA{175, 135, 135, 255}, // 138
        RGBA{175, 135, 175, 255}, // 139
        RGBA{175, 135, 215, 255}, // 140
        RGBA{175, 135, 255, 255}, // 141
        RGBA{175, 175, 0, 255}, // 142
        RGBA{175, 175, 95, 255}, // 143
        RGBA{175, 175, 135, 255}, // 144
        RGBA{175, 175, 175, 255}, // 145
        RGBA{175, 175, 215, 255}, // 146
        RGBA{175, 175, 255, 255}, // 147
        RGBA{175, 215, 0, 255}, // 148
        RGBA{175, 215, 95, 255}, // 149
        RGBA{175, 215, 135, 255}, // 150
        RGBA{175, 215, 175, 255}, // 151
        RGBA{175, 215, 215, 255}, // 152
        RGBA{175, 215, 255, 255}, // 153
        RGBA{175, 255, 0, 255}, // 154
        RGBA{175, 255, 95, 255}, // 155
        RGBA{175, 255, 135, 255}, // 156
        RGBA{175, 255, 175, 255}, // 157
        RGBA{175, 255, 215, 255}, // 158
        RGBA{175, 255, 255, 255}, // 159
        RGBA{215, 0, 0, 255}, // 160
        RGBA{215, 0, 95, 255}, // 161
        RGBA{215, 0, 135, 255}, // 162
        RGBA{215, 0, 175, 255}, // 163
        RGBA{215, 0, 215, 255}, // 164
        RGBA{215, 0, 255, 255}, // 165
        RGBA{215, 95, 0, 255}, // 166
        RGBA{215, 95, 95, 255}, // 167
        RGBA{215, 95, 135, 255}, // 168
        RGBA{215, 95, 175, 255}, // 169
        RGBA{215, 95, 215, 255}, // 170
        RGBA{215, 95, 255, 255}, // 171
        RGBA{215, 135, 0, 255}, // 172
        RGBA{215, 135, 95, 255}, // 173
        RGBA{215, 135, 135, 255}, // 174
        RGBA{215, 135, 175, 255}, // 175
        RGBA{215, 135, 215, 255}, // 176
        RGBA{215, 135, 255, 255}, // 177
        RGBA{215, 175, 0, 255}, // 178
        RGBA{215, 175, 95, 255}, // 179
        RGBA{215, 175, 135, 255}, // 180
        RGBA{215, 175, 175, 255}, // 181
        RGBA{215, 175, 215, 255}, // 182
        RGBA{215, 175, 255, 255}, // 183
        RGBA{215, 215, 0, 255}, // 184
        RGBA{215, 215, 95, 255}, // 185
        RGBA{215, 215, 135, 255}, // 186
        RGBA{215, 215, 175, 255}, // 187
        RGBA{215, 215, 215, 255}, // 188
        RGBA{215, 215, 255, 255}, // 189
        RGBA{215, 255, 0, 255}, // 190
        RGBA{215, 255, 95, 255}, // 191
        RGBA{215, 255, 135, 255}, // 192
        RGBA{215, 255, 175, 255}, // 193
        RGBA{215, 255, 215, 255}, // 194
        RGBA{215, 255, 255, 255}, // 195
        RGBA{255, 0, 0, 255}, // 196
        RGBA{255, 0, 95, 255}, // 197
        RGBA{255, 0, 135, 255}, // 198
        RGBA{255, 0, 175, 255}, // 199
        RGBA{255, 0, 215, 255}, // 200
        RGBA{255, 0, 255, 255}, // 201
        RGBA{255, 95, 0, 255}, // 202
        RGBA{255, 95, 95, 255}, // 203
        RGBA{255, 95, 135, 255}, // 204
        RGBA{255, 95, 175, 255}, // 205
        RGBA{255, 95, 215, 255}, // 206
        RGBA{255, 95, 255, 255}, // 207
        RGBA{255, 135, 0, 255}, // 208
        RGBA{255, 135, 95, 255}, // 209
        RGBA{255, 135, 135, 255}, // 210
        RGBA{255, 135, 175, 255}, // 211
        RGBA{255, 135, 215, 255}, // 212
        RGBA{255, 135, 255, 255}, // 213
        RGBA{255, 175, 0, 255}, // 214
        RGBA{255, 175, 95, 255}, // 215
        RGBA{255, 175, 135, 255}, // 216
        RGBA{255, 175, 175, 255}, // 217
        RGBA{255, 175, 215, 255}, // 218
        RGBA{255, 175, 255, 255}, // 219
        RGBA{255, 215, 0, 255}, // 220
        RGBA{255, 215, 95, 255}, // 221
        RGBA{255, 215, 135, 255}, // 222
        RGBA{255, 215, 175, 255}, // 223
        RGBA{255, 215, 215, 255}, // 224
        RGBA{255, 215, 255, 255}, // 225
        RGBA{255, 255, 0, 255}, // 226
        RGBA{255, 255, 95, 255}, // 227
        RGBA{255, 255, 135, 255}, // 228
        RGBA{255, 255, 175, 255}, // 229
        RGBA{255, 255, 215, 255}, // 230
        RGBA{255, 255, 255, 255}, // 231
        RGBA{8, 8, 8, 255}, // 232
        RGBA{18, 18, 18, 255}, // 233
        RGBA{28, 28, 28, 255}, // 234
        RGBA{38, 38, 38, 255}, // 235
        RGBA{48, 48, 48, 255}, // 236
        RGBA{58, 58, 58, 255}, // 237
        RGBA{68, 68, 68, 255}, // 238
        RGBA{78, 78, 78, 255}, // 239
        RGBA{88, 88, 88, 255}, // 240
        RGBA{98, 98, 98, 255}, // 241
        RGBA{108, 108, 108, 255}, // 242
        RGBA{118, 118, 118, 255}, // 243
        RGBA{128, 128, 128, 255}, // 244
        RGBA{138, 138, 138, 255}, // 245
        RGBA{148, 148, 148, 255}, // 246
        RGBA{158, 158, 158, 255}, // 247
        RGBA{168, 168, 168, 255}, // 248
        RGBA{178, 178, 178, 255}, // 249
        RGBA{188, 188, 188, 255}, // 250
        RGBA{198, 198, 198, 255}, // 251
        RGBA{208, 208, 208, 255}, // 252
        RGBA{218, 218, 218, 255}, // 253
        RGBA{228, 228, 228, 255}, // 254
        RGBA{238, 238, 238, 255} // 255
    };

    // color_id (logical palette index / ncurses pair number) → xterm-256 foreground color number.
    // Mirrors Camera::renderPoint2image() so ANSI and PNG output use identical colors.
    inline int palette_to_xterm256_fg(int cid) {
        if      (cid == 10)                return CONSERVATION_UNSCORED_NEAR;
        else if (cid == 20)                return CONSERVATION_UNSCORED_FAR;
        else if (cid == 36)                return CONSERVATION_UNSCORED_MID;
        else if (cid == 47)                return SS_HELIX_NEAR;
        else if (cid == 48)                return SS_HELIX_FAR;
        else if (cid == 49)                return SS_SHEET_NEAR;
        else if (cid == 50)                return SS_SHEET_FAR;
        else if (cid >= 85  && cid <= 93)  return PROTEIN_COIL_FAR[cid - 85];
        else if (cid == 111)               return ALIGNED_NONALIGNED_NEAR;
        else if (cid >= 251 && cid <= 259) return ALIGNED_BRIGHT_FAR[cid - 251];
        else if (cid >= 1   && cid <= 9)   return PROTEIN_COLORS[cid - 1];
        else if (cid >= 11  && cid <= 19)  return PROTEIN_DIM_COLORS[cid - 11];
        else if (cid >= 21  && cid <= 35)  return CHAIN_COLORS[cid - 21];
        else if (cid == 41)                return 226;
        else if (cid == 42)                return 51;
        else if (cid == 43)                return INTERFACE_COLOR;
        else if (cid == 44)                return INTERFACE_DIM_COLOR;
        else if (cid == 45)                return ALIGNED_COLOR;
        else if (cid == 46)                return ALIGNED_DIM_COLOR;
        else if (cid >= 51  && cid <= 70)  return RAINBOW[cid - 51];
        else if (cid >= 71  && cid <= 74)  return PLDDT_COLORS[cid - 71];
        else if (cid >= 75  && cid <= 84)  return CONSERVATION_COLORS[cid - 75];
        else if (cid >= 101 && cid <= 109) return PROTEIN_BRIGHT_COLORS[cid - 101];
        else if (cid == 110)               return ALIGNED_NONALIGNED_DIM;
        else if (cid >= 120 && cid <= 128) return PROTEIN_NEAR_COLORS[cid - 120];
        else if (cid >= 130 && cid <= 144) return CHAIN_NEAR_COLORS[cid - 130];
        else if (cid >= 145 && cid <= 159) return CHAIN_FAR_COLORS[cid - 145];
        else if (cid >= 160 && cid <= 179) return RAINBOW_NEAR[cid - 160];
        else if (cid >= 180 && cid <= 199) return RAINBOW_FAR[cid - 180];
        else if (cid >= 200 && cid <= 208) return PROTEIN_FAR_COLORS[cid - 200];
        else if (cid >= 209 && cid <= 212) return PLDDT_NEAR[cid - 209];
        else if (cid >= 213 && cid <= 216) return PLDDT_FAR[cid - 213];
        else if (cid >= 217 && cid <= 226) return CONSERVATION_NEAR[cid - 217];
        else if (cid >= 227 && cid <= 236) return CONSERVATION_FAR[cid - 227];
        else if (cid == 237)               return INTERFACE_NEAR_COLOR;
        else if (cid == 238)               return INTERFACE_DIM_NEAR_COLOR;
        else if (cid == 239)               return INTERFACE_FAR_COLOR;
        else if (cid == 240)               return INTERFACE_DIM_FAR_COLOR;
        else if (cid >= 241 && cid <= 249) return ALIGNED_BRIGHT_NEAR[cid - 241];
        else if (cid == 250)               return ALIGNED_NONALIGNED_FAR;
        else                               return 231;  // fallback white
    }

    // color_id → ANSI xterm-256 foreground escape sequence "\033[38;5;Nm"
    inline std::string palette_to_ansi_fg_str(int color_id) {
        return "\033[38;5;" + std::to_string(palette_to_xterm256_fg(color_id)) + "m";
    }

    // ANSI attribute reset
    inline const char* palette_to_ansi_reset() {
        return "\033[0m";
    }
}
