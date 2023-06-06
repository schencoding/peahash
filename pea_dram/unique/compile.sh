#!/bin/bash

index_type=(pea clht level dash-16)

# number of threads to run
thread_num=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)

op_type=(full mixed)

# prepare
rm build/test_pmem
rm build/CMakeCache.txt
rm /mnt/mypmem1/liuzhuoxuan/pmem_pea.data
rm /mnt/mypmem1/liuzhuoxuan/pmem_dash-16.data
rm /mnt/mypmem1/liuzhuoxuan/pmem_dash-512.data
rm /mnt/mypmem1/liuzhuoxuan/pmem_level.data
rm /mnt/mypmem1/liuzhuoxuan/pmem_cceh.data

# compile
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_PMEM=ON -S . -B build
cd build && make -j