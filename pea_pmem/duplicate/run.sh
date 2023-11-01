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

# execute
a=0
while(($a<5))
do
    a=`expr $a + 1`
    for h in 0 1 2 3 4      # skew rate
    do
        for k in 0 1  #index
        do
            echo "No.$a ${index_type[$k]} test begin, skew factor ${skew[$h]}"
            # run
            timeout 1h numactl --cpunodebind=1 --membind=1 ./test_pmem \
            -index ${index_type[$k]} \
            -skew ${skew[$h]} \
            >> ../log/skew_uniq_delete${a}.csv
            # clean data
            rm /mnt/mypmem1/liuzhuoxuan/pmem_${index_type[$k]}.data
        done
    done
done    