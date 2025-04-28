# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(hw1 VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(hw1
    src/main.c
    src/algebra.c)

# Include directories (if needed)
target_include_directories(hw1 PRIVATE inc)

