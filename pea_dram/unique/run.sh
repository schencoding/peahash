#!/bin/bash

index_type=(pea clht level dash-16 dash-512 cceh)

# number of threads to run
thread_num=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)

op_type=(full mixed)

# prepare
rm build/test_pmem
rm build/CMakeCache.txt

# compile
cmake -DCMAKE_BUILD_TYPE=Release -S . -B build
cd build && make -j

# execute
a=0
while(($a<1))
do
    a=`expr $a + 1`
    for j in 0          #op
    do
        # for k in 0 1 2 3 4 5 #index
        for k in 0
        do
            # for i in 0 1 3 7 15  #thread
            for i in 0
            do
                echo "No.$a ${op_type[$j]} ${index_type[$k]} ${thread_num[$i]}threads test begin"
                # run
                timeout 1h ./test_pmem \
                -index ${index_type[$k]} \
                -op ${op_type[$j]} \
                -t ${thread_num[$i]} \
                # >> dram_test${a}.csv
                # clean data
            done
        done
    done
done    