cmake_minimum_required(VERSION 3.14)
project(OptimizationMethods VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Source and Include directories
include_directories(include)

# Main executable
add_executable(task3 src/main.cpp 
    src/Vector.cpp 
    src/LennardJones.cpp 
    src/Optimizers.cpp
)

# Testing
enable_testing()
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
)
# For Windows: prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(unit_tests 
    tests/test_vector.cpp 
    tests/test_physics.cpp 
    tests/test_optimizers.cpp
    src/Vector.cpp 
    src/LennardJones.cpp 
    src/Optimizers.cpp
)
target_link_libraries(unit_tests gtest_main)

include(GoogleTest)
gtest_discover_tests(unit_tests)
