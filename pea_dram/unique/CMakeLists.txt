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

set(libs_to_link
        pthread gflags      
        )

add_executable(test_pmem test/test_pmem.cpp)
target_include_directories(test_pmem PRIVATE util)
target_link_libraries(test_pmem PRIVATE ${libs_to_link})


add_definitions(-DAVX512F)
