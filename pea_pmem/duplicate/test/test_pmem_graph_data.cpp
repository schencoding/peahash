// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong. All rights reserved.
// Licensed under the MIT license.
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

// #include "../src/Dash/ex_finger.h"
#include "../src/Dash/ex_fingers.h"
#include "../src/PeaHash/peas_hash.h"
#include "../src/PeaHash/space_manager.h"
#include "../util/System.hpp"
#include "../util/key_generator.hpp"
#include "../util/uniform.hpp"
#include "libpmemobj.h"

#define SAMPLE_NUM  5000000L
#define MAX_LEN     250000000L

std::string insert_pth;
std::string skew_sp_pth;
std::string uniq_sp_pth;
std::string pool_name = "/mnt/mypmem1/liuzhuoxuan/";
DEFINE_string(index, "dash-ex",
              "the index to evaluate:pea/dash-ex/dash-lh/cceh/level");
DEFINE_string(k, "fixed", "the type of stored keys: fixed/variable");
DEFINE_string(distribution, "skew",
              "The distribution of the workload: uniform/skew");
DEFINE_uint64(i, 256, "the initial number of segments in extendible hashing");
DEFINE_uint64(t, 1, "the number of concurrent threads");
DEFINE_uint64(n, 10000000, "the number of pre-insertion load"); 
DEFINE_uint64(loadType, 0, "type of pre-load integers: random (0) - range (1)");
DEFINE_uint64(p, 190000000,      
              "the number of operations(insert/search/deletion) to execute");
DEFINE_string(
    op, "skew-full",
    "which type of operation to execute:insert/pos/neg/delete/mixed/full/recovery/skew-search/skew-full");
DEFINE_double(r, 1, "read ratio for mixed workload:0~1.0");
DEFINE_double(s, 0, "insert ratio for mixed workload: 0~1.0");
DEFINE_double(d, 0, "delete ratio for mixed workload:0~1.0");
DEFINE_double(skew, 0.8, "skew factor of the workload");
DEFINE_uint32(e, 0, "whether register epoch in application level:0/1");
DEFINE_uint32(ms, 100, "#miliseconds to sample the operations");
DEFINE_uint32(vl, 16, "the length of the variable length key");
DEFINE_uint64(ps, 30ul, "The size of the memory pool (GB)");
DEFINE_uint64(ed, 1000, "The frequency to enroll into the epoch");

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
uint32_t msec, value_pool_size;
struct timeval tv1, tv2, tv3;
size_t pool_size = 1024ul * 1024ul * 1024ul * 30ul;
key_generator_t *uniform_generator;
uint64_t EPOCH_DURATION;
uint64_t load_type = 0;
Value_t *value_pool;

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

template <class T>
Hash<T> *InitializeIndex(int seg_num) {
  Hash<T> *eh;
  bool file_exist = false;
  gettimeofday(&tv1, NULL);
  if (index_type == "dash-ex") {
    std::cout << "Initialize Dash-EH" << std::endl;
    std::string index_pool_name = pool_name + "pmem_dash-ex.data";
    if (FileExists(index_pool_name.c_str())) file_exist = true;
    Allocator::Initialize(index_pool_name.c_str(), pool_size);
#ifdef PREALLOC
    extendible::TlsTablePool<Key_t>::Initialize();
#endif
    eh = reinterpret_cast<Hash<T> *>(
        Allocator::GetRoot(sizeof(extendibles::Finger_EH<T>)));
    new (eh) extendibles::Finger_EH<T>(seg_num, Allocator::Get()->pm_pool_);
  } else if (index_type == "pea") {
    std::cout << "Initialize Pea Hash" << std::endl;
    std::string index_pool_name = pool_name + "pmem_pea.data";
    if (FileExists(index_pool_name.c_str())) file_exist = true;
    SpaceManager::Initialize(index_pool_name.c_str(), pool_size, thread_num);
    size_t pea_size = 8;      // pea_size = 8
    if(!file_exist) {
      eh = new peas::PeaHashing<T>(pea_size);
    }
  }

  return eh;
}

template <class T>
void Load(int kv_num, Hash<T> *index, void *workload) {
  // std::cout << "Start load warm-up workload" << std::endl;
  if (kv_num == 0) return;
  std::string fixed("fixed");
  T *_worklod = reinterpret_cast<T *>(workload);
  T key;
  if constexpr (!std::is_pointer_v<T>) {
    for (uint64_t i = 0; i < kv_num; ++i) {
      index->Insert(_worklod[i], DEFAULT);
    }
  }
#ifdef AU
  // pea::space_num = 0;
  // extendible::space_num = 0;
  // level::space_num = 0;
#endif
  // std::cout << "Finish loading " << kv_num << " keys" << std::endl;
}

inline void spin_wait() {
  SUB(&bar_b, 1);
  while (LOAD(&bar_a) == 1)
    ; /*spinning*/
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

template <class T>
void concurr_insert_without_epoch(struct range *_range, Hash<T> *index) {
  set_affinity(_range->index);
  int begin = _range->begin;
  int end = _range->end;
  char *workload = reinterpret_cast<char *>(_range->workload);
  T key;

  spin_wait();
  if constexpr (!std::is_pointer_v<T>) {
    T *key_array = reinterpret_cast<T *>(workload);
    for (uint64_t i = begin; i < end; ++i) {
      int ret = index->Insert(key_array[i], DEFAULT);
    }
  }

  end_notify(_range);
}

template <class T>
void concurr_insert(struct range *_range, Hash<T> *index) {
  set_affinity(_range->index);
  int begin = _range->begin;
  int end = _range->end;
  char *workload = reinterpret_cast<char *>(_range->workload);
  T key;
  uint64_t repeat_key = 0;

  if constexpr (!std::is_pointer_v<T>) {
    T *key_array = reinterpret_cast<T *>(workload);
    uint64_t round = (end - begin) / EPOCH_DURATION;
    uint64_t i = 0;
    spin_wait();

    while (i < round) {
      auto epoch_guard = Allocator::AquireEpochGuard();
      uint64_t _end = begin + (i + 1) * EPOCH_DURATION;
      for (uint64_t j = begin + i * EPOCH_DURATION; j < _end; ++j) {
        index->Insert(key_array[j], DEFAULT, true);
      }
      ++i;
    }

    {
      auto epoch_guard = Allocator::AquireEpochGuard();
      for (i = begin + EPOCH_DURATION * round; i < end; ++i) {
        index->Insert(key_array[i], DEFAULT, true);
      }
    }
  } else {
    T var_key;
    uint64_t round = (end - begin) / EPOCH_DURATION;
    uint64_t i = 0;
    uint64_t string_key_size = sizeof(string_key) + _range->length;

    spin_wait();
    while (i < round) {
      auto epoch_guard = Allocator::AquireEpochGuard();
      uint64_t _end = begin + (i + 1) * EPOCH_DURATION;
      for (uint64_t j = begin + i * EPOCH_DURATION; j < _end; ++j) {
        var_key = reinterpret_cast<T>(workload + string_key_size * j);
        index->Insert(var_key, DEFAULT, true);
      }
      ++i;
    }

    {
      auto epoch_guard = Allocator::AquireEpochGuard();
      for (i = begin + EPOCH_DURATION * round; i < end; ++i) {
        var_key = reinterpret_cast<T>(workload + string_key_size * i);
        index->Insert(var_key, DEFAULT, true);
      }
    }
  }

  end_notify(_range);
}

template <class T>
void concurr_search_without_epoch(struct range *_range, Hash<T> *index) {
  set_affinity(_range->index);
  int begin = _range->begin;
  int end = _range->end;
  char *workload = reinterpret_cast<char *>(_range->workload);
  T key;
  uint64_t not_found = 0;
  uint64_t ent_num = 0;
  spin_wait();

  if constexpr (!std::is_pointer_v<T>) {
    T *key_array = reinterpret_cast<T *>(workload);
    for (uint64_t i = begin; i < end; ++i) {
      uint64_t num = value_pool_size / 2;
      int retval = -1;
      retval = index->Get(key_array[i], value_pool, &num);
      if(retval == -1){
        LOG(i);
        LOG_FATAL("Error, space exceed");
      } else if(num == 0) {
        not_found++;
      } else {
        ent_num += num;
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
  std::cout << "entry number = " << ent_num << std::endl;
  end_notify(_range);
}

template <class T>
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
      index->Delete(key_array[i]);
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
  std::cout << "not_found = " << not_found << std::endl;
  end_notify(_range);
}

template <class T>
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

  uint32_t insert_sign = (uint32_t)(insert_ratio * 100);
  uint32_t read_sign = (uint32_t)(read_ratio * 100) + insert_sign;
  uint32_t delete_sign = (uint32_t)(delete_ratio * 100) + read_sign;

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
  std::cout << "not_found = " << not_found << std::endl;
  /*the last thread notify the main thread to wake up*/
  end_notify(_range);
}

template <class T>
void GeneralBench(range *rarray, Hash<T> *index, int thread_num,
                  uint64_t op_num, std::string profile_name,
                  void (*test_func)(struct range *, Hash<T> *)) {
  std::thread *thread_array[1024];
  profile_name = profile_name + std::to_string(thread_num);
  double duration;
  finished = false;
  bar_a = 1;
  bar_b = thread_num;
  bar_c = thread_num;

  std::cout << profile_name << " Begin" << std::endl;
  //  System::profile(profile_name, [&]() {
  for (uint64_t i = 0; i < thread_num; ++i) {
    thread_array[i] = new std::thread(*test_func, &rarray[i], index);
  }

  while (LOAD(&bar_b) != 0)
    ;                                     // Spin
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

  double longest = (double)(rarray[0].tv.tv_usec - tv1.tv_usec) / 1000000 +
                   (double)(rarray[0].tv.tv_sec - tv1.tv_sec);
  double shortest = longest;
  duration = longest;

  for (int i = 1; i < thread_num; ++i) {
    double interval = (double)(rarray[i].tv.tv_usec - tv1.tv_usec) / 1000000 +
                      (double)(rarray[i].tv.tv_sec - tv1.tv_sec);
    duration += interval;
    if (shortest > interval) shortest = interval;
    if (longest < interval) longest = interval;
  }
  //std::cout << "The time difference is " << longest - shortest << std::endl;
  duration = duration / thread_num;
  printf(
      "%d threads, Time = %f s, throughput = %f "
      "ops/s, fastest = %f, slowest = %f\n",
      thread_num, duration, op_num / duration, op_num / shortest,
      op_num / longest);
  
  std::cout << profile_name << " End" << std::endl;
}

void resample(const uint64_t *workload, uint64_t *sample_arr, uint64_t total_num, uint64_t sample_num) {
  srand(time(nullptr));
  for (size_t i = 0; i < sample_num; i++) {
    int resample_id = rand() % total_num;
    sample_arr[i] = workload[resample_id];
  }
}

void shuffle(uint64_t *workload, uint64_t number) {
  srand(time(nullptr));
  for (size_t i = 0; i < 2*number; i++)
  {
    int replace_id = rand() % number;
    uint64_t swap_one = workload[replace_id];
    workload[replace_id] = workload[0];
    workload[0] = swap_one;
  }
}

template <class T>
void Run() {
  /* Initialize Index for Finger_EH*/
  uniform_generator = new uniform_key_generator_t();
  Hash<T> *index = InitializeIndex<T>(initCap);
  /* Generate the workload and corresponding range array*/
  // std::cout << "Generate load workload" << std::endl;
  // void *workload = malloc(load_num * sizeof(uint64_t));
  void *insert_workload = malloc(MAX_LEN * sizeof(uint64_t));

  // std::cout << "load num = " << load_num << std::endl;
  // Load<T>(load_num, index, workload);

  FILE *ins_stream;
  size_t num_ins;
  if( (ins_stream = fopen(insert_pth.data(), "r" )) != NULL ) {
    num_ins = fread(insert_workload, sizeof(uint64_t), MAX_LEN, ins_stream);
    printf( "Number of items to insert = %ld\n", num_ins );
    fclose(ins_stream);
  } else {
    printf( "File could not be opened\n" );
  }

  // std::cout << "Finish Generate insert workload" << std::endl;

  /* Description of the workload*/
  srand((unsigned)time(NULL));
  struct range *rarray;
  rarray = reinterpret_cast<range *>(malloc(thread_num * (sizeof(range))));
  for (int i = 0; i < thread_num; ++i) {
    rarray[i].index = i;
    rarray[i].random_num = rand();
    rarray[i].begin = 0;
    rarray[i].end = num_ins;
    rarray[i].length = 8;
    rarray[i].workload = insert_workload;
  }
  rarray[thread_num - 1].end = num_ins;

  /* Benchmark Phase */
  if (operation == "skew-full") { /*do the benchmark for all single operations*/
    // std::cout << "Comprehensive Benchmark" << std::endl;
    // std::cout << "insertion start" << std::endl;    

    GeneralBench<T>(rarray, index, thread_num, num_ins, "Insert",
                      &concurr_insert_without_epoch);

    index->getNumber();

    FILE *ss_stream;
    size_t num_ss;
    if( (ss_stream = fopen(skew_sp_pth.data(), "r" )) != NULL ) {
      num_ss = fread(insert_workload, sizeof(uint64_t), SAMPLE_NUM, ss_stream);
      printf( "Number of items to search = %ld\n", num_ss );
      fclose(ss_stream);
    } else {
      printf( "File could not be opened\n" );
    } 
    rarray[0].end = num_ss;
    // skew search
    GeneralBench<T>(rarray, index, thread_num, num_ss, "skew_search",
                        &concurr_search_without_epoch);  

    FILE *us_stream;
    size_t num_us;
    if( (us_stream = fopen(uniq_sp_pth.data(), "r" )) != NULL ) {
      num_us = fread(insert_workload, sizeof(uint64_t), SAMPLE_NUM, us_stream);
      printf( "Number of items to search = %ld\n", num_us );
      fclose(us_stream);
    } else {
      printf( "File could not be opened\n" );
    }
    rarray[0].end = num_us;

    GeneralBench<T>(rarray, index, thread_num, num_us, "uniq_search",
                        &concurr_search_without_epoch); 
    
    shuffle((uint64_t *)insert_workload, num_us);

    GeneralBench<T>(rarray, index, thread_num, num_us, "Delete",
                      &concurr_delete_without_epoch);
    index->getNumber();
#ifdef ANA
    // std::cout << "1 val cnt " << val1cnt << std::endl;
    // std::cout << "2 vals cnt " << val2cnt << std::endl;
    // std::cout << ">= 3: " << valg3cnt << std::endl;
#endif

  }
  /*TODO Free the workload memory*/
  free(insert_workload);
}

bool check_ratio() {
  int read_portion = (int)(read_ratio * 100);
  int insert_portion = (int)(insert_ratio * 100);
  int delete_portion = (int)(delete_ratio * 100);
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
  value_pool_size = operation_num / 2;
  value_pool = (Value_t *)malloc(sizeof(Value_t) * value_pool_size);
  key_type = FLAGS_k;
  index_type = FLAGS_index;
  distribution = FLAGS_distribution;
  load_type = FLAGS_loadType;
  // std::cout << "Distribution = " << distribution << std::endl;
  std::string fixed("fixed");
  operation = FLAGS_op;
  open_epoch = FLAGS_e;
  EPOCH_DURATION = FLAGS_ed;
  msec = FLAGS_ms;
  pool_size = FLAGS_ps * 1024ul * 1024ul * 1024ul; /*pool_size*/
  if (open_epoch == true)
    std::cout << "EPOCH registration in application level" << std::endl;

  read_ratio = FLAGS_r;
  insert_ratio = FLAGS_s;
  delete_ratio = FLAGS_d;
  skew_factor = FLAGS_skew;
  std::cout << "Skew theta = " << skew_factor << std::endl;
  if (skew_factor == 0.0) {
    insert_pth = "../temp/data_190M_s00";    
  } else if (skew_factor == 0.2) {
    insert_pth = "../temp/data_190M_s02";        
  } else if (skew_factor == 0.4) {
    insert_pth = "../temp/data_190M_s04";    
  } else if (skew_factor == 0.6) {
    insert_pth = "../temp/data_190M_s06";        
  } else if (skew_factor == 0.8) {
    insert_pth = "../temp/data_190M_s08";        
  } else if (skew_factor == 1.1) {
    insert_pth = "../keys/soc-Slashdot0902";
  } else if (skew_factor == 1.2) {
    insert_pth = "../keys/wiki-Talk";
  } else if (skew_factor == 1.4) {
    insert_pth = "../keys/cit-Patents";
  } else if (skew_factor == 1.6) {
    insert_pth = "../keys/soc-LiveJournal1";
  } else if (skew_factor == 1.8) {
    insert_pth = "../keys/com-orkut.ungraph";
  } else {
    LOG_FATAL("Please modify file name!");
  }
  if (skew_factor < 1) {
    skew_sp_pth = insert_pth + "_sp20M";
    uniq_sp_pth = insert_pth + "_uniq_sp20M";
  } else {
    skew_sp_pth = insert_pth + "_sp5M";
    uniq_sp_pth = insert_pth + "_uniq_sp";
  }
  
  if (operation == "mixed") {
    std::cout << "Search ratio = " << read_ratio << std::endl;
    std::cout << "Insert ratio = " << insert_ratio << std::endl;
    std::cout << "Delete ratio = " << delete_ratio << std::endl;
  }

  if (!check_ratio()) {
    std::cout << "The ratio is wrong!" << std::endl;
    return 0;
  }

  if (key_type.compare(fixed) == 0) {
    Run<uint64_t>();
  }
}
