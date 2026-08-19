# --------------------------------------------------
# Third-party (submodules)
# --------------------------------------------------

# gemmi 0.5.1 is header-only — it ships no consumable library target (its
# add_subdirectory only builds the gemmi *program* executables). So instead of
# add_subdirectory(gemmi), define an INTERFACE target exposing its headers.
# Guard: when embedded in a parent (e.g. Foldseek) that already provides a
# gemmi_cpp target, skip and reuse theirs.
if(NOT TARGET gemmi_cpp)
    add_library(gemmi_cpp INTERFACE)
    # SYSTEM: gemmi 0.5.1 headers emit warnings (e.g. fileutil.hpp
    # -Wignored-attributes) we don't own; treat them as system headers.
    target_include_directories(gemmi_cpp SYSTEM INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/gemmi/include)
    target_compile_features(gemmi_cpp INTERFACE cxx_std_17)
    # gemmi 0.5.1 headers (atox.hpp, cifdoc.hpp) use std::uint8_t without
    # including <cstdint>; GCC 13+ no longer pulls it in transitively, so
    # force-include it ahead of every gemmi header (GCC/Clang only — all of
    # StrucTTY's targets: Linux/macOS/mingw).
    target_compile_options(gemmi_cpp INTERFACE -include cstdint)
    if(ZLIB_FOUND)
        target_link_libraries(gemmi_cpp INTERFACE ZLIB::ZLIB)
    else()
        # No system zlib: expose gemmi's bundled zlib headers so gz.hpp compiles.
        target_include_directories(gemmi_cpp SYSTEM INTERFACE
            ${CMAKE_CURRENT_SOURCE_DIR}/gemmi/third_party/zlib)
    endif()
    add_library(gemmi::gemmi_cpp ALIAS gemmi_cpp)
endif()

# stb_image_write: single public-domain header, used for PNG screenshots
add_library(stb INTERFACE)

target_include_directories(stb
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/stb
)
