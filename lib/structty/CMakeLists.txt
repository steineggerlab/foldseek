cmake_minimum_required(VERSION 3.15)
project(StrucTTY LANGUAGES C CXX)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set OFF when embedding via add_subdirectory (e.g. Foldseek) to skip the
# StrucTTY executable target. The structty library is always built.
option(STRUCTTY_BUILD_APP  "Build StrucTTY terminal application" ON)
# Set ON to also build the render_test integration test binary.
option(BUILD_RENDER_TEST   "Build render_test integration test binary" OFF)

find_package(ZLIB)  # optional — gemmi bundles its own zlib when not found

# Third-party deps are built in lib/
add_subdirectory(lib)

# Render library (ncurses-free static lib; Renderer.cpp + AnsiOutput.cpp)
add_subdirectory(src/render)

# -----------------------------------------------------------------------
# structty static library — interactive viewer, no ncurses dependency.
# Foldseek embeds this via add_subdirectory(structty) + link structty.
# -----------------------------------------------------------------------
file(GLOB_RECURSE STRUCTTY_LIB_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/structty.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/structure/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/visualization/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/*.cpp"
)

file(GLOB_RECURSE STRUCTTY_LIB_HEADERS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/include/structty.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/structure/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/visualization/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/render/*.hpp"
)

add_library(structty STATIC
    ${STRUCTTY_LIB_SOURCES}
    ${STRUCTTY_LIB_HEADERS}
)

target_include_directories(structty
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/utils
        ${CMAKE_CURRENT_SOURCE_DIR}/src/structure
        ${CMAKE_CURRENT_SOURCE_DIR}/src/visualization
        ${CMAKE_CURRENT_SOURCE_DIR}/src/render
)

target_link_libraries(structty
    PUBLIC
        structty_render
        gemmi::gemmi_cpp
        stb
)
if(ZLIB_FOUND)
    target_link_libraries(structty PUBLIC ZLIB::ZLIB)
endif()

# -----------------------------------------------------------------------
# StrucTTY executable — thin wrapper around structty::run()
# -----------------------------------------------------------------------
if (STRUCTTY_BUILD_APP)
    add_executable(StrucTTY ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
    target_link_libraries(StrucTTY PRIVATE structty)
endif()

if (BUILD_RENDER_TEST)
    add_executable(render_test example/render_test.cpp)
    target_link_libraries(render_test PRIVATE structty_render gemmi::gemmi_cpp)
endif()
