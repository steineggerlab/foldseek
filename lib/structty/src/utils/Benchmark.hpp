// Benchmark.hpp
#pragma once
#include <chrono>
#include <cstdint>
#include <fstream>
#include <string>

struct Benchmark {
    using clock = std::chrono::steady_clock;

    bool enabled = false;
    std::ofstream out;

    clock::time_point t0;
    clock::time_point last_fps_t;
    uint64_t frames_since_last = 0;

    // last input event (for input-to-frame latency)
    bool has_event = false;
    int last_key = -1;
    clock::time_point last_event_t;

    void start(const std::string& path) {
        enabled = true;
        out.open(path, std::ios::out);
        t0 = clock::now();
        last_fps_t = t0;
        frames_since_last = 0;
        has_event = false;
        last_key = -1;

        // CSV header
        out << "t_ms,tag,key,dt_ms,num_ca,num_file\n";
        out.flush();
    }

    static inline int64_t ms_since(clock::time_point a, clock::time_point b) {
        return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count();
    }

    // 서브페이즈 계측용: ms 는 소형 구조에서 0 으로 뭉개지므로 µs 해상도 사용.
    static inline int64_t us_since(clock::time_point a, clock::time_point b) {
        return std::chrono::duration_cast<std::chrono::microseconds>(b - a).count();
    }

    void log(const char* tag, int key, int64_t dt_ms, int64_t num_ca = 0, int64_t num_file = 0) {
        if (!enabled) return;
        auto now = clock::now();
        out << ms_since(t0, now) << "," << tag << "," << key << "," << dt_ms
            << "," << num_ca << "," << num_file << "\n";
    }

    void mark_event(int key) {
        if (!enabled) return;
        has_event = true;
        last_key = key;
        last_event_t = clock::now();
    }

    // call once per completed frame flush
    void mark_frame_end(int64_t render_dt_ms, int64_t total_ca = 0, int64_t structs = 1) {
        if (!enabled) return;

        // FPS bookkeeping (optional log cadence)
        frames_since_last++;
        auto now = clock::now();
        auto elapsed_ms = ms_since(last_fps_t, now);
        if (elapsed_ms >= 1000) {
            double fps = (double)frames_since_last * 1000.0 / (double)elapsed_ms;
            log("fps", -1, (int64_t)(fps + 0.5), total_ca, structs);
            frames_since_last = 0;
            last_fps_t = now;
        }

        // input-to-frame latency: event -> this frame end
        if (has_event) {
            int64_t lat_ms = ms_since(last_event_t, now);
            log("latency", last_key, lat_ms, total_ca, structs);
            has_event = false; // consume one event
        }

        // per-frame render time
        log("frame", -1, render_dt_ms, total_ca, structs);
        out.flush();  // 부분 실행(중단) 시에도 계측 데이터 보존
    }
};
