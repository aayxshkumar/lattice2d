cmake_minimum_required(VERSION 3.24)
project(lattice2d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(LATTICE_BUILD_TESTS "Build unit tests" ON)
# e.g. -DLATTICE_SANITIZERS=address,undefined
set(LATTICE_SANITIZERS "" CACHE STRING "Comma-separated -fsanitize= list")

if(LATTICE_SANITIZERS)
  add_compile_options(-fsanitize=${LATTICE_SANITIZERS} -fno-omit-frame-pointer -fno-sanitize-recover=all)
  add_link_options(-fsanitize=${LATTICE_SANITIZERS})
endif()

# No fused multiply-add contraction: results then depend only on the source
# and libm, not on whether the compiler chose to fuse a*b+c.
add_compile_options(-ffp-contract=off)
set(LATTICE_WARNINGS -Wall -Wextra -Wpedantic -Werror)

add_library(lattice_core src/shape.cpp src/collision.cpp src/world.cpp src/scenes.cpp)
target_include_directories(lattice_core PUBLIC include)
target_compile_options(lattice_core PRIVATE ${LATTICE_WARNINGS})

add_executable(lattice_sim tools/lattice_sim.cpp)
target_link_libraries(lattice_sim PRIVATE lattice_core)
target_compile_options(lattice_sim PRIVATE ${LATTICE_WARNINGS})

if(LATTICE_BUILD_TESTS)
  include(FetchContent)
  # Vendored tarball when present (offline machines), else GitHub; same hash.
  set(_gtest_tarball ${CMAKE_CURRENT_SOURCE_DIR}/third_party/googletest-1.15.2.tar.gz)
  if(EXISTS ${_gtest_tarball})
    set(_gtest_url ${_gtest_tarball})
  else()
    set(_gtest_url https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz)
  endif()
  FetchContent_Declare(googletest
    URL ${_gtest_url}
    URL_HASH SHA256=7b42b4d6ed48810c5362c265a17faebe90dc2373c885e5216439d37927f02926
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  enable_testing()
  add_executable(lattice_tests tests/shape_test.cpp tests/collision_test.cpp tests/world_test.cpp)
  target_link_libraries(lattice_tests PRIVATE lattice_core GTest::gtest_main)
  target_compile_options(lattice_tests PRIVATE ${LATTICE_WARNINGS})
  include(GoogleTest)
  gtest_discover_tests(lattice_tests DISCOVERY_MODE PRE_TEST DISCOVERY_TIMEOUT 60)
  add_test(NAME physical_checks COMMAND lattice_sim --check)
endif()
