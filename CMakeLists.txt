cmake_minimum_required(VERSION 3.27)
project(ConservationLawsProject)

set(CMAKE_CXX_STANDARD 20)

find_package (Eigen3 3.4 REQUIRED)
find_package (Ceres REQUIRED)

add_executable(ConservationLawsProject main.cpp
        Solver.cpp
        Geometry.cpp
        initial.cpp
        AnalyticalSolution.cpp
        Godunov.cpp
        Godunov.hpp)

target_link_libraries(ConservationLawsProject Eigen3::Eigen Ceres::ceres)
