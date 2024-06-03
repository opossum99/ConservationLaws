cmake_minimum_required(VERSION 3.27)
project(HPC_Project)

set(CMAKE_CXX_STANDARD 20)

find_package (Eigen3 3.4 REQUIRED)

add_executable(HPC_Project main.cpp
        Solver.cpp
        initial.cpp)

target_link_libraries(HPC_Project Eigen3::Eigen)
