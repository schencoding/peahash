#!/bin/bash

index_type=(pea dash-ex)

skew=(0.0 0.2 0.4 0.6 0.8)

# prepare
rm build/test_pmem
rm build/CMakeCache.txt
rm /mnt/mypmem1/liuzhuoxuan/pmem_pea.data
rm /mnt/mypmem1/liuzhuoxuan/pmem_dash-ex.data

# prepare
cd ../../
make clean

# compile
cmake -DCMAKE_BUILD_TYPE=Release -S ./pea_pmem/duplicate -B build
cd build && make -j