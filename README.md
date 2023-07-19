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

3. Since [Dash](https://github.com/baotonglu/dash) use its Customized PMDK, we install PMDK in our CMakeLists. Installing PMDK is a little difficult. Please refer to the Dependencies of https://github.com/pmem/pmdk if incurring error.



## Building
CMakeLists, test/test_pmem*.cpp, compile.sh
These are important test files.

PS. If your server has no access to github.com, please refer to pea_pmem/unique/CMakeLists.txt.offline for help.


# Availability and Reproducibility for SIGMOD'23
- pea_dram/unique: Fig.6
- pea_pmem/unique: Fig.4, Fig.5
- pea_pmem/duplicate: Fig.8, Fig.9, Fig.10
- pea_pmem/different_strategy: Fig.7