add_library(structty_render STATIC
    Renderer.cpp
    AnsiOutput.cpp
)

# Use generator expressions so include paths are correct in both build-tree
# (add_subdirectory) and install-tree (find_package) contexts.
target_include_directories(structty_render
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/visualization>
        $<INSTALL_INTERFACE:include>
)

# No ncurses linkage — intentional; consumers verify with:
#   nm -u libstructty_render.a | grep -i ncurses
# ── Install rules ──────────────────────────────────────────────────────────────
include(GNUInstallDirs)

install(TARGETS structty_render
    EXPORT  structty_renderTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Public headers: structty_render.h + render headers + RenderPoint (used by both)
install(FILES
    ${PROJECT_SOURCE_DIR}/include/structty_render.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Renderer.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AnsiOutput.hpp
    ${PROJECT_SOURCE_DIR}/src/visualization/RenderPoint.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# CMake package files for find_package(structty_render) support
install(EXPORT structty_renderTargets
    FILE      structty_renderTargets.cmake
    NAMESPACE structty_render::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/structty_render
)

install(FILES
    ${PROJECT_SOURCE_DIR}/cmake/structty_renderConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/structty_render
)
