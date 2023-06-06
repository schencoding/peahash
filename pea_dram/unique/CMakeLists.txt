if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif ()

project(pehash)
cmake_minimum_required(VERSION 3.11)
include(FetchContent)
include(ExternalProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long -fPIC -march=native")
set(CMAKE_ENABLE_COMPILE_COMMANDS "ON")

if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
  message("-- Release mode, all optimizations enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -Ofast")
endif ()

option(USE_PMEM "enable persistent memory support" ON)

set(libs_to_link
        pthread gflags    
        )

if (USE_PMEM MATCHES "ON")
  message(STATUS "add definition PMEM")
  add_definitions(-DPMEM)
endif ()

if (USE_PMEM MATCHES "ON")
  add_executable(test_pmem test/test_pmem.cpp)
  # add_executable(example test/example.cpp)
  # add_executable(gen_data test/gen_data.cpp)
  # add_executable(resample test/resample.cpp)
  target_include_directories(test_pmem PRIVATE util)
  target_link_libraries(test_pmem PRIVATE ${libs_to_link})
  # target_link_libraries(example PRIVATE ${libs_to_link})
  add_dependencies(test_pmem PMDK)
  # add_dependencies(example PMDK)
else ()
endif ()

add_definitions(-DAVX512F)
