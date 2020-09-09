cmake_minimum_required(VERSION 3.9)

project(Fredholm2D_test CXX)

set(LIB_FREDHOLM_2D_INCLUDE_DIR ../src)

find_library(LIB_FREDHOLM_2D  Fredholm2D ../)

find_package(OpenMP REQUIRED)

add_executable(${PROJECT_NAME} test.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC ${LIB_FREDHOLM_2D_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${LIB_FREDHOLM_2D} OpenMP::OpenMP_CXX)