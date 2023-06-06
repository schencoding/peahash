# PeaHash: Performant Extendible Adaptive Hashing Index

Code for SIGMOD'23 Pea Hash Paper


## Directory
- pea_dram/unique: Support unique keys in DRAM
- pea_pmem/unique: Support unique keys in PMEM
- pea_pmem/duplicate: Support duplicate keys in PMEM

## What's included
- src/ is the main source files of pea hash.
- util/ is some util functions.
- test/test_pmem is an examplary test file.
- CMakeLists.txt is a template which must be changed.
- Compile.sh shows the usage of the test_pmem.

## Building
CMakeLists, test/test_pmem.cpp, compile.sh
These are important test files.
Please refer the three important files to compile your own example.


## Dependencies
No dependency is needed in Pea Hash.

If you see some strange dependencies in the example test file and CMakeLists, these are obsolete dependencies, please delete them.

If you want to use the test_pmem, 

set(libs_to_link gflags)

in CMakeLists.