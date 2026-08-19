#pragma once
#include "Protein.hpp"
#include "Atom.hpp"
#include "Palette.hpp"
#include "Camera.hpp"
#include "Common.hpp"
#include "RenderPoint.hpp"
#include <ctime>
#include <thread>
#include <chrono>  
#include <vector>
#include <iostream>
#include <filesystem> 
#include <unordered_map>
#include <iomanip>
#include <sstream>

inline std::string get_home_dir() {
    if (const char* home = std::getenv("HOME")) return std::string(home);
    throw std::runtime_error("HOME is not set");
}

class Camera {
public:
    Camera(const int width, const int height, const std::string mode);
    ~Camera() = default;

    void screenshot(const std::vector<RenderPoint>& pixels, int width, int height);
    void renderPoint2image(const std::vector<RenderPoint>& pixels, int width, int height,
                           std::vector<RGBA>& screenImage);
    bool save_image(std::vector<RGBA>& screenImage, int width, int height);
    int get_alpha_from_depth(float z, float min_z, float max_z);

private:
    float focal_offset = 10.0f;

    std::string camera_dir;
    int camera_width, camera_height;
    int screenshot_idx;

    int height_duplicate = 2;

    std::string camera_mode = "default";
};