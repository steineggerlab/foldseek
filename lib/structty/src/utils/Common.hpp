#pragma once
#include <ctime>
#include <chrono>  
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

static std::string current_timestamp() {
    using namespace std::chrono;

    auto now = system_clock::now();
    auto t = system_clock::to_time_t(now);
    auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

    std::tm tm{};
    #if defined(_WIN32)
        localtime_s(&tm, &t);
    #else
        localtime_r(&t, &tm);
    #endif


    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d_%H%M%S")
        << "_" << std::setw(3) << std::setfill('0') << ms.count();

    return oss.str();
}

// align 계열 색 모드: 자동(align) · foldseek 정렬 구간만(align-fs) · 거리 기반(align-near)
static bool is_aligned_mode(const std::string& mode) {
    return mode == "align" || mode == "align-fs" || mode == "align-near";
}
