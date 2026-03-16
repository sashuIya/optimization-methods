cmake_minimum_required(VERSION 3.14)
project(OptimizationMethods VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Source directory
include_directories(src)

# Main executable
add_executable(task3 
    src/main.cc 
    src/Vector.cc 
    src/LennardJones.cc 
    src/Optimizers.cc
)

# Testing
enable_testing()
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(unit_tests 
    tests/test_vector.cpp 
    tests/test_physics.cpp 
    tests/test_optimizers.cpp
    src/Vector.cc 
    src/LennardJones.cc 
    src/Optimizers.cc
)
target_link_libraries(unit_tests gtest_main)

include(GoogleTest)
gtest_discover_tests(unit_tests)
