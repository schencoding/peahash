#!/bin/bash

index_type=(pea dash-512 dash-16 cceh level)

# number of threads to run
thread_num=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)

op_type=(full mixed)

# prepare
make clean

# compile
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_PMEM=ON -S . -B build
cd build && make -j

# execute
a=0
while(($a<5))
do
    a=`expr $a + 1`
    for j in 0          #op
    do
        for k in 0 #index
        do
            # for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15  #thread
            for i in 0
            do
                echo "No.$a ${op_type[$j]} ${index_type[$k]} ${thread_num[$i]}threads test begin"
                # run
                timeout 1h numactl --cpunodebind=1 --membind=1 ./test_pmem \
                -index ${index_type[$k]} \
                -op ${op_type[$j]} \
                -t ${thread_num[$i]} \
                >> ../log/full${a}.csv
                # clean data
                rm /mnt/mypmem1/liuzhuoxuan/pmem_${index_type[$k]}.data
            done
        done
    done
done    