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

##################### PMDK ####################
set(PMDK_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/pmdk")
ExternalProject_Add(PMDK
        GIT_REPOSITORY file:///home/liuzhuoxuan/workspace/pmdk
        GIT_TAG 1.7.1
	BUILD_IN_SOURCE 1
        BUILD_COMMAND $(MAKE)
        PREFIX ${PMDK_PREFIX}
        CONFIGURE_COMMAND ""
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
        )

include_directories(${PMDK_PREFIX}/src/PMDK/src/include)
if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  link_directories(${PMDK_PREFIX}/src/PMDK/src/debug)
else ()
  link_directories(${PMDK_PREFIX}/src/PMDK/src/nondebug)
endif ()


##################### Epoch Reclaim ####################
FetchContent_Declare(
        epoch_reclaimer
        GIT_REPOSITORY file:///home/liuzhuoxuan/workspace/epoch-reclaimer
#        GIT_TAG master
)
if (NOT epoch_reclaimer_POPULATED)
  FetchContent_Populate(epoch_reclaimer)
  include_directories(${epoch_reclaimer_SOURCE_DIR})
endif ()

option(USE_PMEM "enable persistent memory support" ON)

set(libs_to_link
        pthread gflags)

if (USE_PMEM MATCHES "ON")
  message(STATUS "persistent memory support enabled, going to build with PMDK")
  add_definitions(-DPMEM)
  list(APPEND libs_to_link pmemobj pmem)
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
# add_definitions(-DDEBUG)
# add_definitions(-DAU)
# add_definitions(-DANA)
