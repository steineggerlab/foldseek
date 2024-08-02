add_library(
  prostt5 STATIC
    prostt5.h
    prostt5.cpp
)
target_link_libraries(prostt5 PRIVATE ${PROJECT_NAME})

if (BUILD_CLI)
  add_executable(prostt5cli
    prostt5_cli.cc)
  target_include_directories(prostt5cli
    PRIVATE 
      ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/cxxopts/include
      ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(prostt5cli PRIVATE ${PROJECT_NAME} prostt5)
  install(
    TARGETS prostt5cli
    DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()