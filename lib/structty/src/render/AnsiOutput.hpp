#pragma once
#include <string>
#include <vector>
#include "RenderPoint.hpp"

class AnsiOutput {
public:
    // Convert a logical pixel buffer to an ANSI escape string.
    // logical_width  = Renderer::get_logical_width()  (2 * terminal_cols)
    // logical_height = Renderer::get_logical_height() (4 * terminal_rows)
    // Each terminal cell covers a 2-wide × 4-tall sub-pixel block (Braille rendering).
    static std::string to_ansi_string(
        const std::vector<RenderPoint>& pixels,
        int logical_width,
        int logical_height
    );

    // Write to_ansi_string() result directly to stdout.
    // Calls setlocale(LC_ALL, "") to ensure UTF-8 Braille output is decoded correctly.
    static void print_to_stdout(
        const std::vector<RenderPoint>& pixels,
        int logical_width,
        int logical_height
    );
};
