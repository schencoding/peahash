if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif ()

project(dash)
cmake_minimum_required(VERSION 3.10)
#include(FetchContent)
#include(ExternalProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long -fPIC -march=native")
set(CMAKE_ENABLE_COMPILE_COMMANDS "ON")

if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
  message("-- Release mode, all optimizations enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -Ofast")
endif ()

##################### PMDK ####################
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(PMDK REQUIRED)

include_directories(build/_deps/epoch_reclaimer-src)

option(USE_PMEM "enable persistent memory support" ON)

set(libs_to_link
        pthread gflags)

if (USE_PMEM MATCHES "ON")
  message(STATUS "persistent memory support enabled, going to build with PMDK")
  add_definitions(-DPMEM)
  list(APPEND libs_to_link ${PMDK_LIBRARIES})
endif ()

if (USE_PMEM MATCHES "ON")
  add_executable(test_pmem test/test_pmem.cpp)
  target_include_directories(test_pmem PRIVATE util)
  target_link_libraries(test_pmem PRIVATE ${libs_to_link})
else ()
  add_executable(gen_data test/gen_data_graph.cpp)
  # add_executable(gen_sample test/resample.cpp)
endif ()

add_definitions(-DAVX512F)
