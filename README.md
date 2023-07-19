# PeaHash: Performant Extendible Adaptive Hashing Index

Code for SIGMOD'23 Pea Hash Paper


## Directory
- pea_dram/unique: Support unique keys in DRAM
- pea_pmem/unique: Support unique keys in PMEM
- pea_pmem/duplicate: Support duplicate keys in PMEM

## What's included
- src/ is the main source files of pea hash.
- util/ is some util functions.
- test/ is test files.
- Compile.sh shows the usage of the test_pmem.

## Prerequisite

1. Hardware: Intel Optane Memory

2. install cmake,  g++

3. Installing PMDK is a little difficult. Please refer to the Dependencies of https://github.com/pmem/pmdk if incurring error.

## Building
CMakeLists, test/test_pmem*.cpp, compile.sh
These are important test files.

PS. If your server has no access to github.com, please refer to pea_pmem/unique/CMakeLists.txt.offline for help.


## Running benchmark

As stated in our paper, we run the tests in a single NUMA node with 24 physical CPU cores. We pin threads to physical cores compactly assuming thread ID == core ID (e.g., for a dual-socket system, we assume cores 0-23 are located in socket 0, and cores 24-47 in socket 1).  To run benchmarks, use the `test_pmem` executable in the `build` directory. It supports the following arguments:

```bash
./build/test_pmem --helpshort
Usage: 
    ./build/test_pmem [OPTION...]

-index      the index to evaluate:dash-ex/dash-lh/cceh/level (default: "dash-ex")
-op         the type of operation to execute:insert/pos/neg/delete/mixed (default: "full")
-n          the number of warm-up workload (default: 0)
-p          the number of operations(insert/search/delete) to execute (default: 20000000)
-t          the number of concurrent threads (default: 1)
-r          search ratio for mixed workload: 0.0~1.0 (default: 1.0)
-s          insert ratio for mixed workload: 0.0~1.0 (default: 0.0)
-d          delete ratio for mixed workload: 0.0~1.0 (default: 0.0)
-e          whether to register epoch in application level: 0/1 (default: 0)
-k          the type of stored keys: fixed/variable (default: "fixed")
-vl         the length of the variable length key (default: 16)
```
Check out also the `run.sh` script for example benchmarks and easy testing of the hash tables. 


# Availability and Reproducibility for SIGMOD'23
- pea_dram/unique: Fig.6
- pea_pmem/unique: Fig.4, Fig.5
- pea_pmem/duplicate: Fig.8, Fig.9, Fig.10
- pea_pmem/different_strategy: Fig.7