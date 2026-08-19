#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "Camera.hpp"

Camera::Camera(const int width, const int height, const std::string mode){
    camera_width = width;
    camera_height = height;
    camera_mode = mode;
    
    camera_dir = get_home_dir() + "/Pictures/StrucTTY_screenshot/";

    try {
        std::filesystem::create_directories(camera_dir);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Failed to create screenshot directory: " << e.what() << std::endl;
    }
}

void Camera::screenshot(const std::vector<RenderPoint>& pixels, int width, int height){
    std::vector<RGBA> screenImage;
    screenImage.assign(width * height * height_duplicate, RGBA{0,0,0,0});
    renderPoint2image(pixels, width, height, screenImage);
    if(save_image(screenImage, width, height)){
        std::cout << "\nScreenshot saved in " << camera_dir << std::endl;
    } else{
        std::cout << "\nFailed to save screenshot in " << camera_dir << std::endl;
    }
    std::this_thread::sleep_for(std::chrono::seconds(1));
}

void Camera::renderPoint2image(const std::vector<RenderPoint>& pixels, int width, int height,
                               std::vector<RGBA>& screenImage){
    float min_depth = std::numeric_limits<float>::max();
    float max_depth = std::numeric_limits<float>::lowest();
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int d = 0; d < height_duplicate; d++){
                int cid = pixels[y * width + x].color_id;

                if (cid <= 0) {
                } else {
                    int xterm;
                    if      (cid == 10)                xterm = Palettes::CONSERVATION_UNSCORED_NEAR;
                    else if (cid == 20)                xterm = Palettes::CONSERVATION_UNSCORED_FAR;
                    else if (cid == 36)                xterm = Palettes::CONSERVATION_UNSCORED_MID;
                    else if (cid == 47)                xterm = Palettes::SS_HELIX_NEAR;
                    else if (cid == 48)                xterm = Palettes::SS_HELIX_FAR;
                    else if (cid == 49)                xterm = Palettes::SS_SHEET_NEAR;
                    else if (cid == 50)                xterm = Palettes::SS_SHEET_FAR;
                    else if (cid >= 85 && cid <= 93)   xterm = Palettes::PROTEIN_COIL_FAR[cid - 85];
                    else if (cid == 111)               xterm = Palettes::ALIGNED_NONALIGNED_NEAR;
                    else if (cid >= 251 && cid <= 259) xterm = Palettes::ALIGNED_BRIGHT_FAR[cid - 251];
                    else if (cid >= 1  && cid <= 9)    xterm = Palettes::PROTEIN_COLORS[cid - 1];
                    else if (cid >= 11 && cid <= 19)   xterm = Palettes::PROTEIN_DIM_COLORS[cid - 11];
                    else if (cid >= 21 && cid <= 35)   xterm = Palettes::CHAIN_COLORS[cid - 21];
                    else if (cid == 41)                xterm = 226;  // yellow helix
                    else if (cid == 42)                xterm = 51;   // cyan sheet
                    else if (cid == 43)                xterm = Palettes::INTERFACE_COLOR;
                    else if (cid == 44)                xterm = Palettes::INTERFACE_DIM_COLOR;
                    else if (cid == 45)                xterm = Palettes::ALIGNED_COLOR;
                    else if (cid == 46)                xterm = Palettes::ALIGNED_DIM_COLOR;
                    else if (cid >= 51 && cid <= 70)   xterm = Palettes::RAINBOW[cid - 51];
                    else if (cid >= 71 && cid <= 74)   xterm = Palettes::PLDDT_COLORS[cid - 71];
                    else if (cid >= 75 && cid <= 84)   xterm = Palettes::CONSERVATION_COLORS[cid - 75];
                    else if (cid >= 101 && cid <= 109) xterm = Palettes::PROTEIN_BRIGHT_COLORS[cid - 101];
                    else if (cid == 110)               xterm = Palettes::ALIGNED_NONALIGNED_DIM;
                    else if (cid >= 120 && cid <= 128) xterm = Palettes::PROTEIN_NEAR_COLORS[cid - 120];
                    else if (cid >= 130 && cid <= 144) xterm = Palettes::CHAIN_NEAR_COLORS[cid - 130];
                    else if (cid >= 145 && cid <= 159) xterm = Palettes::CHAIN_FAR_COLORS[cid - 145];
                    else if (cid >= 160 && cid <= 179) xterm = Palettes::RAINBOW_NEAR[cid - 160];
                    else if (cid >= 180 && cid <= 199) xterm = Palettes::RAINBOW_FAR[cid - 180];
                    else if (cid >= 200 && cid <= 208) xterm = Palettes::PROTEIN_FAR_COLORS[cid - 200];
                    else if (cid >= 209 && cid <= 212) xterm = Palettes::PLDDT_NEAR[cid - 209];
                    else if (cid >= 213 && cid <= 216) xterm = Palettes::PLDDT_FAR[cid - 213];
                    else if (cid >= 217 && cid <= 226) xterm = Palettes::CONSERVATION_NEAR[cid - 217];
                    else if (cid >= 227 && cid <= 236) xterm = Palettes::CONSERVATION_FAR[cid - 227];
                    else if (cid == 237)               xterm = Palettes::INTERFACE_NEAR_COLOR;
                    else if (cid == 238)               xterm = Palettes::INTERFACE_DIM_NEAR_COLOR;
                    else if (cid == 239)               xterm = Palettes::INTERFACE_FAR_COLOR;
                    else if (cid == 240)               xterm = Palettes::INTERFACE_DIM_FAR_COLOR;
                    else if (cid >= 241 && cid <= 249) xterm = Palettes::ALIGNED_BRIGHT_NEAR[cid - 241];
                    else if (cid == 250)               xterm = Palettes::ALIGNED_NONALIGNED_FAR;
                    else                               xterm = 231;  // fallback white
                    screenImage[((y * height_duplicate) + d) * width + x] =
                        Palettes::ID2RGBA[xterm];
                }
            }

            min_depth = std::min(min_depth, pixels[y * width + x].depth);
            if (pixels[y * width + x].depth != std::numeric_limits<float>::infinity())
                max_depth = std::max(max_depth, pixels[y * width + x].depth);
        }
    }
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int d = 0; d < height_duplicate; d++){
                screenImage[((y * height_duplicate) + d) * width + x].a = get_alpha_from_depth(pixels[y * width + x].depth, min_depth, max_depth);
            }
        }
    }

    return;
}

bool Camera::save_image(std::vector<RGBA>& screenImage, int width, int height){
    std::string screenshot_dir = camera_dir + current_timestamp() + ".png";
    const unsigned char* bytes = reinterpret_cast<const unsigned char*>(screenImage.data());
    int ok = stbi_write_png(screenshot_dir.c_str(), width, height * height_duplicate, 4, bytes, width * 4);
    if (!ok) { return false; }
    else { return true; }
}

int Camera::get_alpha_from_depth(float z, float min_z, float max_z) {
    float zn = (z - min_z) / (max_z - min_z);

    if (zn < 0.1) return 255;
    else if (zn < 0.26) return 218;
    else if (zn < 0.42) return 182;
    else if (zn < 0.58) return 145;
    else if (zn < 0.72) return 109;
    else if (zn < 0.9) return 73;
    else if (zn < 1.0) return 36;
    else return 0;
}