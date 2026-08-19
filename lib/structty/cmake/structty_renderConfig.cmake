# structty_renderConfig.cmake
#
# Supports: find_package(structty_render REQUIRED)
# After finding: target_link_libraries(my_target PRIVATE structty_render::structty_render)
#
# Consumers must also link gemmi::gemmi_cpp (or an equivalent that provides
# Gemmi symbols) because structty_render's static archive references them.
# Foldseek already links Gemmi, so no extra work is needed there.

cmake_minimum_required(VERSION 3.15)
include(CMakeFindDependencyMacro)
find_dependency(ZLIB)

include("${CMAKE_CURRENT_LIST_DIR}/structty_renderTargets.cmake")
