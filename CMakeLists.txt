cmake_minimum_required(VERSION 3.13)
project(linfin)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

add_subdirectory(src)
add_subdirectory(libs)
