file(GLOB_RECURSE STRUCT_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/structure/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/visualization/*.cpp"
)

file(GLOB_RECURSE STRUCT_HEADERS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/structure/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/visualization/*.hpp"
)

add_library(structty_core STATIC
    ${STRUCT_SOURCES}
    ${STRUCT_HEADERS}
)

target_include_directories(structty_core
    PUBLIC     
        ${CMAKE_CURRENT_SOURCE_DIR}/utils
        ${CMAKE_CURRENT_SOURCE_DIR}/structure
        ${CMAKE_CURRENT_SOURCE_DIR}/visualization
)

target_link_libraries(structty_core
    PUBLIC
        structty_utils   # gemmi + stb
)
