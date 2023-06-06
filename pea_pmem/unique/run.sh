#!/bin/bash

index_type=(pea dash-512 dash-16 cceh level)

# number of threads to run
thread_num=(1 2 4 8 16)

op_type=(insert recovery)

# prepare
rm build/test_pmem
rm build/CMakeCache.txt

# compile
cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_PMEM=ON -S . -B build
# cmake -DCMAKE_BUILD_TYPE=Release -DUSE_PMEM=ON -S . -B build
cd build && make -j

# execute
a=0
while(($a<1))
do
    a=`expr $a + 1`
    for k in 3   #index
    do
        for i in 0  #thread
        do
            # echo "No.$a ${index_type[$k]} 1 threads recover test begin"
            # run
            # timeout 1h numactl --cpunodebind=1 --membind=1 ./test_pmem \
            # -index ${index_type[$k]} \
            # -op ${op_type[0]}
            timeout 1h numactl --cpunodebind=1 --membind=1 ./test_pmem \
            -index ${index_type[$k]} \
            -op ${op_type[1]}
            # >> ../log/recover${a}.csv
            # clean data
            # rm /mnt/mypmem1/liuzhuoxuan/pmem_${index_type[$k]}.data
        done
    done
done    