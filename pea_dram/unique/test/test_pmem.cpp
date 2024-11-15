// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong. All rights reserved.
// Licensed under the MIT license.

/**
 * 
 * 测mix需要在search, delete中加negotiate 
 *  测full不要
 */

#include <gflags/gflags.h>
#include <immintrin.h>
#include <sys/time.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstring>
#include <mutex>
#include <thread>

#ifndef PEA_ONLY
#include "../src/CCEH/CCEH_baseline.h"
#include "../src/Dash-16/ex_finger.h"
#include "../src/Dash-512/ex_finger.h"
#include "../src/Level/level_base.h"
#include "../src/CLHT/clht_lb_res.h"
#endif
#include "../src/PeaHash/pea_hash.h"
#include "../src/PeaHash/space_manager.h"
#include "../util/key_generator.hpp"
#include "../util/uniform.hpp"

std::string pool_name = "/mnt/mypmem1/liuzhuoxuan/";
// std::string pool_name = "/home/liuzhuoxuan/workspace/peahash-test/workload/";
DEFINE_string(index, "pea",
              "the index to evaluate:pea/dash-16/dash-512/cceh/level/clht");
DEFINE_string(k, "fixed", "the type of stored keys: fixed/variable");
DEFINE_string(distribution, "uniform",
              "The distribution of the workload: uniform/skew");
DEFINE_uint64(i, 256, "the initial number of segments in pea/extendible hashing");
DEFINE_uint64(t, 1, "the number of concurrent threads");
DEFINE_uint64(n, 1000000, "the number of pre-insertion load");
DEFINE_uint64(loadType, 0, "type of pre-load integers: random (0) - range (1)");
DEFINE_uint64(p, 69000000,
              "the number of operations(insert/search/deletion) to execute");
DEFINE_string(
        op, "full",
        "which type of operation to execute:insert/pos/neg/delete/mixed/skew-all/full/skew-search");
DEFINE_double(r, 0.8, "read ratio for mixed workload:0~1.0");
DEFINE_double(s, 0.2, "insert ratio for mixed workload: 0~1.0");
DEFINE_double(d, 0, "delete ratio for mixed workload:0~1.0");
DEFINE_double(skew, 0.8, "skew factor of the workload");
DEFINE_uint32(e, 1, "whether register epoch in application level:0/1");
DEFINE_uint32(ms, 100, "#miliseconds to sample the operations");
DEFINE_uint32(vl, 16, "the length of the variable length key");
DEFINE_uint64(ps, 10ul, "The size of the memory pool (GB)");
DEFINE_uint64(ed, 5, "The frequency to enroll into the epoch");

uint64_t initCap, thread_num, load_num, operation_num;
std::string operation;
std::string distribution;
std::string key_type;
std::string index_type;
int bar_a, bar_b, bar_c;
double read_ratio, insert_ratio, delete_ratio, skew_factor;
std::mutex mtx;
std::condition_variable cv;
bool finished = false;
bool open_epoch;
uint32_t msec, var_length;
struct timeval tv1, tv2, tv3;
size_t pool_size = 1024ul * 1024ul * 1024ul * 5ul;
key_generator_t *uniform_generator;
uint64_t EPOCH_DURATION;
uint64_t load_type = 0;

struct operation_record_t {
    uint64_t number;
    uint64_t dummy[7]; /*patch to a cacheline size, avoid false sharing*/
};

operation_record_t operation_record[1024];

struct range {
    int index;
    uint64_t begin;
    uint64_t end;
    int length; /*if this is the variable length key, use this parameter to
                 indicate the length of the key*/
    void *workload;
    uint64_t random_num;
    struct timeval tv;
};

void set_affinity(uint32_t idx) {
    cpu_set_t my_set;
    CPU_ZERO(&my_set);
    CPU_SET(2 * idx + 1, &my_set);
    sched_setaffinity(0, sizeof(cpu_set_t), &my_set);
}

template<class T>
Hash<T> *InitializeIndex(int seg_num) {
    Hash<T> *eh;
    bool file_exist = false;
    gettimeofday(&tv1, NULL);
#ifndef PEA_ONLY
    if (index_type == "dash-16") {
        std::cout << "Dash-16" << std::endl;
#ifdef PREALLOC
        extendible_16::TlsTablePool<Key_t>::Initialize();
#endif
        eh = new extendible_16::Finger_EH<T>(seg_num);
    } else if (index_type == "dash-512") {
        std::cout << "Dash-512" << std::endl;
#ifdef PREALLOC
        extendible_512::TlsTablePool<Key_t>::Initialize();
#endif
        eh = new extendible_512::Finger_EH<T>(8);
    } else if (index_type == "cceh") {
        std::cout << "CCEH" << std::endl;
        eh = new cceh::CCEH<T>(seg_num);
    } else if (index_type == "level") {
        std::cout << "Level Hashing" << std::endl;
        eh = new level::LevelHashing<T>();
        int level_size = 14;
        level::initialize_level(reinterpret_cast<level::LevelHashing<T> *>(eh),
                                &level_size);
    } else if (index_type == "clht") {
        std::cout << "CLHT" << std::endl;
        eh = new clht_lb::CLHT<T>(81920);
    } else 
#endif
    if (index_type == "pea") {
        std::cout << "Pea" << std::endl;
        std::string index_pool_name = pool_name + "pmem_pea.data";
        if (FileExists(index_pool_name.c_str())) file_exist = true;
        SpaceManager::Initialize(index_pool_name.c_str(), pool_size, thread_num);
        size_t pea_size = 8;
        eh = new pea::PeaHashing<T>(pea_size, thread_num);
    } 
    if (operation == "recovery") {
        gettimeofday(&tv3, NULL);  // test end
        eh->Recovery();
        gettimeofday(&tv2, NULL);  // test end
        double duration = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
                          (double) (tv2.tv_sec - tv1.tv_sec);
        double scanning_time = (double) (tv2.tv_usec - tv3.tv_usec) / 1000000 +
                               (double) (tv2.tv_sec - tv3.tv_sec);
        std::cout << "The recovery algorithm time = " << scanning_time << std::endl;
        std::cout << "Total recovery time (open pool + recovery algorithm) = "
                  << duration << std::endl;
    }

    return eh;
}

/*generate 8-byte number and store it in the memory_region*/
void generate_8B(void *memory_region, uint64_t generate_num, bool persist,
                 key_generator_t *key_generator) {
    uint64_t *array = reinterpret_cast<uint64_t *>(memory_region);

    for (uint64_t i = 0; i < generate_num; ++i) {
        array[i] = key_generator->next_uint64();
    }
}

void skew_generate_8B(void *memory_region, uint64_t generate_num, bool persist,
                      zipfian_key_generator_t *key_generator) {
    uint64_t *array = reinterpret_cast<uint64_t *>(memory_region);

    for (uint64_t i = 0; i < generate_num; ++i) {
        array[i] = key_generator->next_uint64_skew(4999);
    }
}

/*generate 16-byte string and store it in the memory_region*/
void generate_16B(void *memory_region, uint64_t generate_num, int length,
                  bool persist, key_generator_t *key_generator) {
    string_key *var_key;
    uint64_t *_key;
    uint64_t random_num;
    char *workload = reinterpret_cast<char *>(memory_region);

    int word_num = (length / 8) + (((length % 8) != 0) ? 1 : 0);
    _key = reinterpret_cast<uint64_t *>(malloc(word_num * sizeof(uint64_t)));

    for (uint64_t i = 0; i < generate_num; ++i) {
        var_key = reinterpret_cast<string_key *>(workload +
                                                 i * (length + sizeof(string_key)));
        var_key->length = length;
        random_num = key_generator->next_uint64();
        for (int j = 0; j < word_num; ++j) {
            _key[j] = random_num;
        }
        memcpy(var_key->key, _key, length);
    }
}

template<class T>
void Load(int kv_num, Hash<T> *index, int length, void *workload) {
    std::cout << "Start load warm-up workload" << std::endl;
    if (kv_num == 0) return;
    std::string fixed("fixed");
    T *_worklod = reinterpret_cast<T *>(workload);
    T key;
    if constexpr (!std::is_pointer_v<T>) {
        for (uint64_t i = 0; i < kv_num; ++i) {
            index->Insert(_worklod[i], DEFAULT, true);
        }
    }
    std::cout << "Finish loading " << kv_num << " keys" << std::endl;
}

inline void spin_wait() {
    SUB(&bar_b, 1);
    while (LOAD(&bar_a) == 1); /*spinning*/
}

inline void end_notify(struct range *rg) {
    gettimeofday(&rg->tv, NULL);
    if (SUB(&bar_c, 1) == 0) {
        std::unique_lock<std::mutex> lck(mtx);
        finished = true;
        cv.notify_one();
    }
}

inline void end_sub() { SUB(&bar_c, 1); }

template<class T>
void concurr_insert_without_epoch(struct range *_range, Hash<T> *index) {
    set_affinity(_range->index);
    int begin = _range->begin;
    int end = _range->end;
    char *workload = reinterpret_cast<char *>(_range->workload);
    T key;

    spin_wait();
    T *key_array = reinterpret_cast<T *>(workload);
    for (uint64_t i = begin; i < end; ++i) {
        index->Insert(key_array[i], DEFAULT);
    }

    end_notify(_range);
}


template<class T>
void concurr_search_without_epoch(struct range *_range, Hash<T> *index) {
    set_affinity(_range->index);
    int begin = _range->begin;
    int end = _range->end;
    char *workload = reinterpret_cast<char *>(_range->workload);
    T key;
    uint64_t not_found = 0;

    spin_wait();

    if constexpr (!std::is_pointer_v<T>) {
        T *key_array = reinterpret_cast<T *>(workload);
        for (uint64_t i = begin; i < end; ++i) {
            if (index->Get(key_array[i]) == NONE) {
                not_found++;
            }
        }
    } else {
        T var_key;
        uint64_t string_key_size = sizeof(string_key) + _range->length;
        for (uint64_t i = begin; i < end; ++i) {
            var_key = reinterpret_cast<T>(workload + string_key_size * i);
            if (index->Get(var_key) == NONE) {
                not_found++;
            }
        }
    }
    std::cout << "not_found = " << not_found << std::endl;
    end_notify(_range);
}

template<class T>
void concurr_delete_without_epoch(struct range *_range, Hash<T> *index) {
    set_affinity(_range->index);
    int begin = _range->begin;
    int end = _range->end;
    char *workload = reinterpret_cast<char *>(_range->workload);
    T key;
    uint64_t not_found = 0;

    spin_wait();
    if constexpr (!std::is_pointer_v<T>) {
        T *key_array = reinterpret_cast<T *>(workload);
        for (uint64_t i = begin; i < end; ++i) {
            if (index->Delete(key_array[i]) == false) {
                not_found++;
            }
        }
    } else {
        T var_key;
        int string_key_size = sizeof(string_key) + _range->length;
        for (uint64_t i = begin; i < end; ++i) {
            var_key = reinterpret_cast<T>(workload + string_key_size * i);
            if (index->Delete(var_key) == false) {
                not_found++;
            }
        }
    }
    // std::cout << "not_found = " << not_found << std::endl;
    end_notify(_range);
}

template<class T>
void mixed_without_epoch(struct range *_range, Hash<T> *index) {
    set_affinity(_range->index);
    uint64_t begin = _range->begin;
    uint64_t end = _range->end;
    char *workload = reinterpret_cast<char *>(_range->workload);
    T *key_array = reinterpret_cast<T *>(_range->workload);
    T key;
    int string_key_size = sizeof(string_key) + _range->length;

    UniformRandom rng(_range->random_num);
    uint32_t random;
    uint32_t not_found = 0;

    uint32_t insert_sign = (uint32_t) (insert_ratio * 100);
    uint32_t read_sign = (uint32_t) (read_ratio * 100) + insert_sign;
    uint32_t delete_sign = (uint32_t) (delete_ratio * 100) + read_sign;

    spin_wait();

    for (uint64_t i = begin; i < end; ++i) {
        if constexpr (std::is_pointer_v<T>) { /* variable length*/
            key = reinterpret_cast<T>(workload + string_key_size * i);
        } else {
            key = key_array[i];
        }

        random = rng.next_uint32() % 100;
        if (random < insert_sign) { /*insert*/
            index->Insert(key, DEFAULT);
        } else if (random < read_sign) { /*get*/
            if (index->Get(key) == NONE) {
                not_found++;
            }
        } else { /*delete*/
            index->Delete(key);
        }
    }
    // std::cout << "not_found = " << not_found << std::endl;
    /*the last thread notify the main thread to wake up*/
    end_notify(_range);
}

template<class T>
void GeneralBench(range *rarray, Hash<T> *index, int thread_num,
                  uint64_t operation_num, std::string profile_name,
                  void (*test_func)(struct range *, Hash<T> *)) {
    std::thread *thread_array[1024];
    profile_name = profile_name + std::to_string(thread_num);
    double duration;
    finished = false;
    bar_a = 1;
    bar_b = thread_num;
    bar_c = thread_num;

    std::cout << profile_name << std::endl;
    //  System::profile(profile_name, [&]() {
    for (uint64_t i = 0; i < thread_num; ++i) {
        thread_array[i] = new std::thread(*test_func, &rarray[i], index);
    }

    while (LOAD(&bar_b) != 0);                                     // Spin
    std::unique_lock<std::mutex> lck(mtx);  // get the lock of condition variable

    gettimeofday(&tv1, NULL);
    STORE(&bar_a, 0);  // start test
    while (!finished) {
        cv.wait(lck);  // go to sleep and wait for the wake-up from child threads
    }
    gettimeofday(&tv2, NULL);  // test end

    for (int i = 0; i < thread_num; ++i) {
        thread_array[i]->join();
        delete thread_array[i];
    }

    double longest = (double) (rarray[0].tv.tv_usec - tv1.tv_usec) / 1000000 +
                     (double) (rarray[0].tv.tv_sec - tv1.tv_sec);
    double shortest = longest;
    duration = longest;

    for (int i = 1; i < thread_num; ++i) {
        double interval = (double) (rarray[i].tv.tv_usec - tv1.tv_usec) / 1000000 +
                          (double) (rarray[i].tv.tv_sec - tv1.tv_sec);
        duration += interval;
        if (shortest > interval) shortest = interval;
        if (longest < interval) longest = interval;
    }
    //std::cout << "The time difference is " << longest - shortest << std::endl;
    duration = duration / thread_num;
    printf(
            "%f,%f,%f\n",
            operation_num / duration, operation_num / shortest, operation_num / longest);
    //  });
    // std::cout << profile_name << " End" << std::endl;
}

void *GenerateWorkload(uint64_t generate_num, int length) {
    /*Since there are both positive search and negative search, it should generate
     * 2 * generate_num workload*/
    void *workload;
    if (key_type == "fixed") {
        workload = malloc(generate_num * sizeof(uint64_t));
        generate_8B(workload, generate_num, false, uniform_generator);
    } else { /*Generate the variable lengh workload*/
        std::cout << "Genereate workload for variable length key " << std::endl;
        workload = malloc(generate_num * (length + sizeof(string_key)));
        generate_16B(workload, generate_num, length, false, uniform_generator);
        std::cout << "Finish Generation" << std::endl;
    }
    return workload;
}

void *GenerateSkewWorkload(uint64_t load_num, uint64_t exist_num,
                           uint64_t non_exist_num, int length) {
    void *workload;
    if (key_type == "fixed") {
        workload =
                malloc((load_num + exist_num + non_exist_num) * sizeof(uint64_t));
        uint64_t *fixed_workload = reinterpret_cast<uint64_t *>(workload);
        /* For the warm-up workload, it is generated using uniform generator*/
        if (load_type == 1) {
            key_generator_t *range_generator = new range_key_generator_t(1);
            generate_8B(fixed_workload, load_num, false, range_generator);
            delete range_generator;
        } else {
            generate_8B(fixed_workload, load_num, false, uniform_generator);
        }

        if (exist_num) {
            zipfian_key_generator_t *skew_generator =
                    new zipfian_key_generator_t(1, exist_num, skew_factor);
            if (operation == "skew-search") {
                skew_generate_8B(fixed_workload + load_num, exist_num, false, skew_generator);
            } else {
                generate_8B(fixed_workload + load_num, exist_num, false, skew_generator);
            }
            delete skew_generator;
        }

        if (non_exist_num) {
            key_generator_t *skew_generator = new zipfian_key_generator_t(
                    exist_num + load_num, exist_num + non_exist_num + load_num,
                    skew_factor);
            generate_8B(fixed_workload + load_num + exist_num, non_exist_num, false,
                        skew_generator);
            delete skew_generator;
        }
    } else { /*Generate the variable lengh workload*/
        std::cout << "Genereate workload for variable length key " << std::endl;
        workload = malloc((load_num + exist_num + non_exist_num) *
                          (length + sizeof(string_key)));
        if (load_type == 1) {
            key_generator_t *range_generator = new range_key_generator_t(1);
            generate_16B(workload, load_num, length, false, range_generator);
            delete range_generator;
        } else {
            generate_16B(workload, load_num, length, false, uniform_generator);
        }

        char *char_workload = reinterpret_cast<char *>(workload);
        if (exist_num) {
            key_generator_t *skew_generator =
                    new zipfian_key_generator_t(1, exist_num, skew_factor);
            generate_16B(char_workload + load_num * (length + sizeof(string_key)),
                         exist_num, length, false, skew_generator);
            delete skew_generator;
        }

        if (non_exist_num) {
            key_generator_t *skew_generator = new zipfian_key_generator_t(
                    exist_num + load_num, exist_num + non_exist_num + load_num,
                    skew_factor);
            generate_16B(char_workload +
                         (load_num + exist_num) * (length + sizeof(string_key)),
                         non_exist_num, length, false, skew_generator);
            delete skew_generator;
        }
        std::cout << "Finish Generation" << std::endl;
    }
    return workload;
}

void shuffle(uint64_t *workload, uint64_t number) {
    srand(time(nullptr));
    for (size_t i = 0; i < 2 * number; i++) {
        int replace_id = rand() % number;
        uint64_t swap_one = workload[replace_id];
        workload[replace_id] = workload[0];
        workload[0] = swap_one;
    }
}

template<class T>
void Run() {
    /* Index for Finger_EH*/
    uniform_generator = new uniform_key_generator_t();
    Hash<T> *index = InitializeIndex<T>(initCap);
    uint64_t generate_num = operation_num * 2 + load_num;
    /* Generate the workload and corresponding range array*/
    // std::cout << "Generate workload" << std::endl;
    void *workload;
    if (distribution == "uniform") {
        workload = GenerateWorkload(generate_num, var_length);
    } else {
        workload = GenerateSkewWorkload(load_num, operation_num, operation_num,
                                        var_length);
    }

    void *insert_workload;
    insert_workload = workload;
    // std::cout << "Finish Generate workload" << std::endl;

    // std::cout << "load num = " << load_num << std::endl;
    Load<T>(load_num, index, var_length, insert_workload);
    void *not_used_workload;
    void *not_used_insert_workload;

    if (key_type == "fixed") {
        uint64_t *key_array = reinterpret_cast<uint64_t *>(workload);
        not_used_workload = reinterpret_cast<void *>(key_array + load_num);
        not_used_insert_workload = not_used_workload;
    }

    /* Description of the workload*/
    srand((unsigned) time(NULL));
    struct range *rarray;
    uint64_t chunk_size = operation_num / thread_num;
    rarray = reinterpret_cast<range *>(malloc(thread_num * (sizeof(range))));
    for (int i = 0; i < thread_num; ++i) {
        rarray[i].index = i;
        rarray[i].random_num = rand();
        rarray[i].begin = i * chunk_size;
        rarray[i].end = (i + 1) * chunk_size;
        rarray[i].length = var_length;
        rarray[i].workload = not_used_workload;
    }
    rarray[thread_num - 1].end = operation_num;

    /* Benchmark Phase */
    std::cout << "Comprehensive Benchmark" << std::endl;
    std::cout << "insertion start" << std::endl;
    for (int i = 0; i < thread_num; ++i) {
        rarray[i].workload = not_used_insert_workload;
    }
    if (operation != "skew-all") {
        GeneralBench<T>(rarray, index, thread_num, operation_num, "Insert",
                        &concurr_insert_without_epoch);
    }
    // index->getNumber();
    shuffle((uint64_t *) not_used_workload, operation_num);
    for (int i = 0; i < thread_num; ++i) {
        rarray[i].workload = not_used_workload;
    }
    GeneralBench<T>(rarray, index, thread_num, operation_num, "Pos_search",
                    &concurr_search_without_epoch);
    for (int i = 0; i < thread_num; ++i) {
        rarray[i].begin = operation_num + i * chunk_size;
        rarray[i].end = operation_num + (i + 1) * chunk_size;
    }
    rarray[thread_num - 1].end = 2 * operation_num;
    GeneralBench<T>(rarray, index, thread_num, operation_num, "Neg_search",
                    &concurr_search_without_epoch);

    for (int i = 0; i < thread_num; ++i) {
        rarray[i].begin = i * chunk_size;
        rarray[i].end = (i + 1) * chunk_size;
    }
    rarray[thread_num - 1].end = operation_num;
    shuffle((uint64_t *) not_used_workload, operation_num);
    GeneralBench<T>(rarray, index, thread_num, operation_num, "Delete",
                    &concurr_delete_without_epoch);

    // index->getNumber();


    /*TODO Free the workload memory*/
}

bool check_ratio() {
    int read_portion = (int) (read_ratio * 100);
    int insert_portion = (int) (insert_ratio * 100);
    int delete_portion = (int) (delete_ratio * 100);
    if ((read_portion + insert_portion + delete_portion) != 100) return false;
    return true;
}

int main(int argc, char *argv[]) {
    set_affinity(0);
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    initCap = FLAGS_i;
    thread_num = FLAGS_t;
    load_num = FLAGS_n;
    operation_num = FLAGS_p;
    key_type = FLAGS_k;
    index_type = FLAGS_index;
    distribution = FLAGS_distribution;
    load_type = FLAGS_loadType;
    std::cout << "Distribution = " << distribution << std::endl;
    std::string fixed("fixed");
    operation = FLAGS_op;
    open_epoch = FLAGS_e;
    EPOCH_DURATION = FLAGS_ed;
    msec = FLAGS_ms;
    var_length = FLAGS_vl;
    pool_size = FLAGS_ps * 1024ul * 1024ul * 1024ul; /*pool_size*/
    if (open_epoch == true)
        std::cout << "EPOCH registration in application level" << std::endl;

    read_ratio = FLAGS_r;
    insert_ratio = FLAGS_s;
    delete_ratio = FLAGS_d;
    skew_factor = FLAGS_skew;
    if (distribution == "skew")
        std::cout << "Skew theta = " << skew_factor << std::endl;

    // if (operation == "mixed") {
    //   std::cout << "Search ratio = " << read_ratio << std::endl;
    //   std::cout << "Insert ratio = " << insert_ratio << std::endl;
    //   std::cout << "Delete ratio = " << delete_ratio << std::endl;
    // }

    if (!check_ratio()) {
        std::cout << "The ratio is wrong!" << std::endl;
        return 0;
    }

    if (key_type.compare(fixed) == 0) {
        Run<uint64_t>();
    } else {
        std::cout << "Variable-length key = " << var_length << std::endl;
        Run<string_key *>();
    }
}
