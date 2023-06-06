#!/bin/sh

# rm build/CMakeCache.txt
# rm build/test_pmem
rm build/gen_data
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_PMEM=OFF -S . -B build
cd build && make -j
./gen_data