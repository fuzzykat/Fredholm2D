cmake_minimum_required(VERSION 3.9)

project(Fredholm2D CXX)

add_library(${PROJECT_NAME} SHARED 
    src/FE2DSolver.cpp  
    src/GaussNodes.cpp  
    src/LAOperations.cpp  
    src/LASSolver.cpp  
    src/Spline2D.cpp)
    
find_package(OpenMP REQUIRED)

target_include_directories(${PROJECT_NAME} PRIVATE src)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)