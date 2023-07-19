#!/bin/bash

index_type=(pea clht level dash-16)

# number of threads to run
thread_num=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)

op_type=(full mixed)

# prepare
cd ../../
make clean

# compile
cmake -DCMAKE_BUILD_TYPE=Release -S ./pea_dram/unique -B build
cd build && make -j