cmake_minimum_required(VERSION 2.8)

project(ig_component_splitter CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${IGREC_MAIN_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

link_libraries(input ${COMMON_LIBRARIES})
link_libraries(boost_program_options)

add_executable(ig_component_splitter
  ig_component_splitter.cpp
  utils.cpp
  config.cpp
  component.cpp)
target_link_libraries(ig_component_splitter build_info)
