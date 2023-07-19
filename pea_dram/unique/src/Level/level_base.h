#ifndef LEVEL_HASHING_H_
#define LEVEL_HASHING_H_

/*
We do several optimization and correctness patches for level hashing, including:
(1) remove fence between storing value and storing key during insert() because
these two stores are in the same cacheline and will mot be reordered. (2) use
lock-striping rather than slot-based locking to make the lock area fitting into
the CPU cache, improving the scalability. (3) use pthread lock to ensure concurrency, original version has some problem
(4) add uniqnuess check during the insert opeartion to avoid inserting duplicate keys.
*/

#include <inttypes.h>
#include <stdint.h>
#include <sys/stat.h>

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <mutex>
#include <shared_mutex>

#include "../../util/hash.h"
#include "../../util/pair.h"
#include "../../util/utils.h"
#include "../Hash.h"
#include "spinlock.h"
#include <pthread.h>

#define ASSOC_NUM 3
#define NODE_TYPE 1000
#define LEVEL_TYPE 2000
#define LOCK_TYPE 3000
// #define COUNTING 1
//#define BATCH 1
#define UNIQUE_CHECK 1
namespace level {
#ifdef AU
uint64_t space_num = 0;
uint64_t entry_num = 0;
#endif

template <class T>
struct Entry {
  T key;
  Value_t value;
  Entry() {
    key = INVALID;
    value = NONE;
  }
};

/*
 * 1 cacheline setting: ASSOC_NUM = 3; dummy = 13
 * 2 cacheline setting: ASSOC_NUM = 7; dummy = 9
 * 4 cacheline setting: ASSOC_NUM = 15; dummy = 1
 */

template <class T>
struct Node {
  uint8_t token[ASSOC_NUM];
  char dummy[13];
  Entry<T> slot[ASSOC_NUM];
  void *operator new[](size_t size) {
    void *ret;
    posix_memalign(&ret, 64, size);
    return ret;
  }

  void *operator new(size_t size) {
    void *ret;
    posix_memalign(&ret, 64, size);
    return ret;
  }
};

typedef struct level_lock{
    spinlock s_lock[ASSOC_NUM];
} level_lock;

template <class T>
class LevelHashing : public Hash<T> {
 public:
  Node<T> *buckets[2];
  Node<T> *interim_level_buckets;
  level_lock* level_locks[2];          // Allocate a fine-grained lock for each slot
  uint64_t level_item_num[2];

  uint64_t levels;
  uint64_t addr_capacity;
  uint64_t total_capacity;
  uint64_t f_seed;
  uint64_t s_seed;
  uint32_t resize_num;
  /* lock information */
  std::atomic<int64_t> resizing_lock;
  pthread_rwlock_t *mutex;
  int nlocks;
  int locksize;
  bool resizing;

  void generate_seeds(void);
  void resize();
  int b2t_movement(uint64_t);
  uint8_t try_movement(uint64_t, uint64_t, T, Value_t);
  uint64_t F_HASH(T key);
  uint64_t S_HASH(T key);

  LevelHashing(void);
  int Get(T key, Value_t *ret_val, size_t *number){
    return 0;
  }
  LevelHashing(size_t);
  ~LevelHashing(void);
  void display_size() {
    std::cout << "The Node size is " << sizeof(Node<T>) << std::endl;
    std::cout << "The Entry size is " << sizeof(Entry<T>) << std::endl;

    if (buckets[0]) {
      printf("0x%" PRIXPTR "\n", (uintptr_t)buckets[0]);
    }

    if (buckets[1]) {
      printf("0x%" PRIXPTR "\n", (uintptr_t)buckets[1]);
    }
  }

  int Insert(T, Value_t);
#ifdef AU
  int Insert0(T, Value_t);
#endif
  int Insert(T key, Value_t value, bool is_in_epoch) {
    return Insert(key, value);
  }
  bool Delete(T);
  bool Delete(T key, bool is_in_epoch) { return Delete(key); }
  Value_t Get(T);
  Value_t Get(T key, bool flag) { return Get(key); }
  void Recovery() {
  }
  uint64_t getNumber() {
    // std::cout << "Entry Size: " << sizeof(struct Entry<T>) << std::endl;
    // std::cout << "Node Size: " << sizeof(struct Node<T>) << std::endl;
    std::cout << "First items " << level_item_num[0] << std::endl;
    std::cout << "Second items " << level_item_num[1] << std::endl;
  //  std::cout << "First load factor = "
  //            << (double)level_item_num[0] / (addr_capacity * ASSOC_NUM)
  //            << std::endl;
  //  std::cout << "Second load factor = "
  //            << (double)level_item_num[1] / (addr_capacity * ASSOC_NUM / 2)
  //            << std::endl;
#ifdef AU
    double average_utility = ((entry_num + 1.0) * entry_num ) / (2.0 * space_num);
    std::cout << "entry count = "<< entry_num << std::endl;
    std::cout << "space unit = " << space_num << std::endl;
    std::cout << "average utility = " << average_utility << std::endl;
#endif
    return 0;
  }
  size_t Capacity(void) {
    return (addr_capacity + addr_capacity / 2) * ASSOC_NUM;
  }
};

#define F_IDX(hash, capacity) (hash % (capacity / 2))
#define S_IDX(hash, capacity) ((hash % (capacity / 2)) + (capacity / 2))

template <class T>
uint64_t LevelHashing<T>::F_HASH(T key) {
    return h(&key, sizeof(Key_t), f_seed);
}

template <class T>
uint64_t LevelHashing<T>::S_HASH(T key) {
    return h(&key, sizeof(Key_t), s_seed);
}

template <class T>
void LevelHashing<T>::generate_seeds() {
  srand(time(NULL));
  do {
    f_seed = rand();
    s_seed = rand();
    f_seed = f_seed << (rand() % 63);
    s_seed = s_seed << (rand() % 63);
  } while (f_seed == s_seed);
}

template <class T>
LevelHashing<T>::LevelHashing(void) {}

template <class T>
LevelHashing<T>::~LevelHashing(void) {}

void *cache_align(void *ptr) {
  uint64_t pp = (uint64_t)ptr;
  pp += 48;
  return (void *)pp;
}

/* Initialize Function for Level Hashing*/
template <class T>
void initialize_level(LevelHashing<T> *level, void *arg) {

          /* modify*/
          level->levels = *((int *)arg);
          level->resize_num = 0;
          level->resizing_lock = 0;
          level->resizing = false;
          level->addr_capacity = pow(2, level->levels);
          level->total_capacity = pow(2, level->levels) + pow(2, level->levels - 1);
          level->locksize = 128;
          level->nlocks = (3 * level->addr_capacity / 2) / level->locksize + 1;

          level->level_locks[0] = (level_lock *)calloc(pow(2, level->levels), sizeof(level_lock));
          level->level_locks[1] = (level_lock *)calloc(pow(2, level->levels - 1), sizeof(level_lock));

          level->generate_seeds();
          level->level_item_num[0] = 0;
          level->level_item_num[1] = 0;
          level->interim_level_buckets = NULL;
          /*allocate*/

          /* Intialize pointer*/
          level->buckets[0] = (Node<T> *)calloc(level->addr_capacity, sizeof(Node<T>));
          level->buckets[1] = (Node<T> *)calloc(level->addr_capacity / 2, sizeof(Node<T>));
          level->mutex = (pthread_rwlock_t *)calloc(sizeof(pthread_rwlock_t), level->nlocks);
          for (size_t i = 0; i < level->nlocks; i++)
          {
            level->mutex[i] = PTHREAD_RWLOCK_INITIALIZER;
          }
}

template <class T>
#ifdef AU
int LevelHashing<T>::Insert0(T key, Value_t value) {
#else
  int LevelHashing<T>::Insert(T key, Value_t value) {
# endif
  RETRY:
  while (resizing_lock.load() == 1) {
    asm("nop");
  }
  uint64_t f_hash = F_HASH(key);
  uint64_t s_hash = S_HASH(key);
  uint32_t f_idx = F_IDX(f_hash, addr_capacity);
  uint32_t s_idx = S_IDX(s_hash, addr_capacity);

#ifdef UNIQUE_CHECK
  uint32_t lock_idx = f_idx / locksize;
  while (pthread_rwlock_trywrlock(&mutex[lock_idx]) != 0) {
    if (resizing == true) {
      goto RETRY;
    }
  }

  if (resizing == true) {
    pthread_rwlock_unlock(&mutex[lock_idx]);
    goto RETRY;
  }

  bool inserted = false;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < ASSOC_NUM; ++j) {
        if ((buckets[i][f_idx].token[j] == 1) &&
            (buckets[i][f_idx].slot[j].key == key)) {
          inserted = true;
          goto UNIQUE;
        }
    }

    for (int j = 0; j < ASSOC_NUM; ++j) {
      if ((buckets[i][s_idx].token[j] == 1) &&
          (buckets[i][s_idx].slot[j].key == key)) {
        inserted = true;
        goto UNIQUE;
      }
    }

    f_idx = F_IDX(f_hash, addr_capacity / 2);
    s_idx = S_IDX(s_hash, addr_capacity / 2);
  }

  UNIQUE:
  if (inserted) {
    pthread_rwlock_unlock(&mutex[lock_idx]);
    return -1; // unique check failure
  }

  pthread_rwlock_unlock(&mutex[lock_idx]);
  f_idx = F_IDX(f_hash, addr_capacity);
  s_idx = S_IDX(s_hash, addr_capacity);
#endif

  int i, j;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < ASSOC_NUM; j++) {
      spin_lock(&this->level_locks[i][f_idx].s_lock[j]);
      if (buckets[i][f_idx].token[j] == 0) {
        buckets[i][f_idx].slot[j].value = value;
        buckets[i][f_idx].slot[j].key = key;
        buckets[i][f_idx].token[j] = 1;
#ifdef COUNTING
        level_item_num[i]++;
#endif
        spin_unlock(&this->level_locks[i][f_idx].s_lock[j]);
        return 0;
      }
      spin_unlock(&this->level_locks[i][f_idx].s_lock[j]);
      spin_lock(&this->level_locks[i][s_idx].s_lock[j]);

      if (buckets[i][s_idx].token[j] == 0) {
        buckets[i][s_idx].slot[j].value = value;
        buckets[i][s_idx].slot[j].key = key;
        buckets[i][s_idx].token[j] = 1;
#ifdef COUNTING
        level_item_num[i]++;
#endif
        spin_unlock(&this->level_locks[i][s_idx].s_lock[j]);
        return 0;
      }
      spin_unlock(&this->level_locks[i][s_idx].s_lock[j]);
    }

    f_idx = F_IDX(f_hash, addr_capacity / 2);
    s_idx = S_IDX(s_hash, addr_capacity / 2);
  }

  f_idx = F_IDX(f_hash, addr_capacity);
  s_idx = S_IDX(s_hash, addr_capacity);
  int empty_loc;
  int64_t lock = 0;
  if (resizing_lock.compare_exchange_strong(lock, 1)) {
    for (i = 0; i < 2; i++) {
      if (!try_movement(f_idx, i, key, value)) {
        resizing_lock.store(0);
        return 0;
      }
      if (!try_movement(s_idx, i, key, value)) {
        resizing_lock.store(0);
        return 0;
      }
      f_idx = F_IDX(f_hash, addr_capacity / 2);
      s_idx = S_IDX(s_hash, addr_capacity / 2);
    }

    if (resize_num > 0) {
      {
        pthread_rwlock_wrlock(&mutex[f_idx / locksize]);
        empty_loc = b2t_movement(f_idx);
        if (empty_loc != -1) {
          buckets[1][f_idx].slot[empty_loc].value = value;
          buckets[1][f_idx].slot[empty_loc].key = key;
          buckets[1][f_idx].token[empty_loc] = 1;
#ifdef COUNTING
          level_item_num[1]++;
#endif
          resizing_lock.store(0);
          pthread_rwlock_unlock(&mutex[f_idx / locksize]);
          return 0;
        }
        pthread_rwlock_unlock(&mutex[f_idx / locksize]);
      }
      {
        pthread_rwlock_wrlock(&mutex[s_idx / locksize]);
        empty_loc = b2t_movement(s_idx);
        if (empty_loc != -1) {
          buckets[1][s_idx].slot[empty_loc].value = value;
          buckets[1][s_idx].slot[empty_loc].key = key;
          buckets[1][s_idx].token[empty_loc] = 1;
#ifdef COUNTING
          level_item_num[1]++;
#endif
          resizing_lock.store(0);
          pthread_rwlock_unlock(&mutex[s_idx / locksize]);
          return 0;
        }
        pthread_rwlock_unlock(&mutex[s_idx / locksize]);
      }
    }
    resize();
    resizing_lock.store(0);
  }
  goto RETRY;
}

template <class T>
void LevelHashing<T>::resize(void) {
  // std::cout << "Resizing towards levels " << levels + 1 << std::endl;
  resizing = true;
  for (int i = 0; i < nlocks; ++i) {
    pthread_rwlock_wrlock(&mutex[i]);
  }

  size_t new_addr_capacity = pow(2, levels + 1);
  free(this->level_locks[0]);
  free(this->level_locks[1]);
  this->level_locks[0] = (level_lock *)calloc(new_addr_capacity, sizeof(level_lock));
  this->level_locks[1] = (level_lock *)calloc(new_addr_capacity / 2, sizeof(level_lock));
  nlocks = (3 * 2 * addr_capacity / 2) / locksize + 1;
  pthread_rwlock_t *old_mut = mutex;
  auto new_mutex = (pthread_rwlock_t *)calloc(sizeof(pthread_rwlock_t), nlocks);
  for (size_t i = 0; i < nlocks; i++) {
    new_mutex[i] = PTHREAD_RWLOCK_INITIALIZER;
  }
  mutex = new_mutex;

  interim_level_buckets = (Node<T> *)calloc(new_addr_capacity, sizeof(Node<T>));
  if (!interim_level_buckets) {
      printf("The resizing fails: 2\n");
      exit(1);
  }

  uint64_t new_level_item_num = 0;
  uint64_t old_idx;
  for (old_idx = 0; old_idx < pow(2, levels - 1); old_idx++) {
    uint64_t i, j;
    for (i = 0; i < ASSOC_NUM; i++) {
      if (buckets[1][old_idx].token[i] == 1) {
        T key = buckets[1][old_idx].slot[i].key;
        Value_t value = buckets[1][old_idx].slot[i].value;

        uint32_t f_idx = F_IDX(F_HASH(key), new_addr_capacity);
        uint32_t s_idx = S_IDX(S_HASH(key), new_addr_capacity);

        uint8_t insertSuccess = 0;
        for (j = 0; j < ASSOC_NUM; j++) {
          if (interim_level_buckets[f_idx].token[j] == 0) {
            interim_level_buckets[f_idx].slot[j].value = value;
            interim_level_buckets[f_idx].slot[j].key = key;

            interim_level_buckets[f_idx].token[j] = 1;

            insertSuccess = 1;
#ifdef COUNTING
            new_level_item_num++;
#endif
            break;
          } else if (interim_level_buckets[s_idx].token[j] == 0) {
            interim_level_buckets[s_idx].slot[j].value = value;
            interim_level_buckets[s_idx].slot[j].key = key;

            interim_level_buckets[s_idx].token[j] = 1;

            insertSuccess = 1;
#ifdef COUNTING
            new_level_item_num++;
#endif
            break;
          }
        }

#ifndef BATCH
        buckets[1][old_idx].token[i] = 0;
#endif
      }
    }
  }


          levels++;
          resize_num++;
          free(this->buckets[1]);
          buckets[1] = buckets[0];
          buckets[0] = interim_level_buckets;
          interim_level_buckets = NULL;
#ifdef COUNTING
          level_item_num[1] = level_item_num[0];
    level_item_num[0] = new_level_item_num;
#endif
          addr_capacity = new_addr_capacity;
          total_capacity = pow(2, levels) + pow(2, levels - 1);
          free(old_mut);
          resizing = false;

  // std::cout << "Done! :Resizing towards levels " << levels << std::endl;
}

template <class T>
uint8_t LevelHashing<T>::try_movement(uint64_t idx,
                                      uint64_t level_num, T key,
                                      Value_t value) {
  uint64_t i, j, jdx;
  
    for (i = 0; i < ASSOC_NUM; i++) {
      spin_lock(&this->level_locks[level_num][idx].s_lock[i]);
      T m_key = buckets[level_num][idx].slot[i].key;
      Value_t m_value = buckets[level_num][idx].slot[i].value;
      uint64_t f_hash = F_HASH(m_key);
      uint64_t s_hash = S_HASH(m_key);
      uint64_t f_idx = F_IDX(f_hash, addr_capacity / (1 + level_num));
      uint64_t s_idx = S_IDX(s_hash, addr_capacity / (1 + level_num));

      if (f_idx == idx)
        jdx = s_idx;
      else
        jdx = f_idx;

      for (j = 0; j < ASSOC_NUM; j++) {
        spin_lock(&this->level_locks[level_num][jdx].s_lock[j]);
        if (buckets[level_num][jdx].token[j] == 0) {
          buckets[level_num][jdx].slot[j].value = m_value;
          buckets[level_num][jdx].slot[j].key = m_key;
          buckets[level_num][jdx].token[j] = 1;
          buckets[level_num][idx].token[i] = 0;
          spin_unlock(&this->level_locks[level_num][jdx].s_lock[j]);

          buckets[level_num][idx].slot[i].value = value;
          buckets[level_num][idx].slot[i].key = key;
          buckets[level_num][idx].token[i] = 1;
          spin_unlock(&this->level_locks[level_num][idx].s_lock[i]);  
#ifdef COUNTING
          level_item_num[level_num]++;
#endif
          return 0;
        }
        spin_unlock(&this->level_locks[level_num][jdx].s_lock[j]);
      }
      spin_unlock(&this->level_locks[level_num][idx].s_lock[i]);        
    }
  return 1;
}

template <class T>
int LevelHashing<T>::b2t_movement(uint64_t idx) {
  T key;
  Value_t value;
  uint64_t s_hash, f_hash;
  uint64_t s_idx, f_idx;
  uint64_t i, j;

  for (i = 0; i < ASSOC_NUM; i++) {
    key = buckets[1][idx].slot[i].key;
    value = buckets[1][idx].slot[i].value;
    f_hash = F_HASH(key);
    s_hash = S_HASH(key);
    f_idx = F_IDX(f_hash, addr_capacity);
    s_idx = S_IDX(s_hash, addr_capacity);

    for (j = 0; j < ASSOC_NUM; j++) {
      if ((idx / locksize) != (f_idx / locksize))
        pthread_rwlock_wrlock(&mutex[f_idx / locksize]);

      if (buckets[0][f_idx].token[j] == 0) {
        buckets[0][f_idx].slot[j].value = value;
        buckets[0][f_idx].slot[j].key = key;
        buckets[0][f_idx].token[j] = 1;
        buckets[1][idx].token[i] = 0;
#ifdef COUNTING
        level_item_num[0]++;
        level_item_num[1]--;
#endif
        if ((idx / locksize) != (f_idx / locksize))
          pthread_rwlock_unlock(&mutex[f_idx / locksize]);
        return i;
      }
      if ((idx / locksize) != (f_idx / locksize))
        pthread_rwlock_unlock(&mutex[f_idx / locksize]);
      if ((idx / locksize) != (s_idx / locksize))
        pthread_rwlock_wrlock(&mutex[s_idx / locksize]);

      if (buckets[0][s_idx].token[j] == 0) {
        buckets[0][s_idx].slot[j].value = value;
        buckets[0][s_idx].slot[j].key = key;
        buckets[0][s_idx].token[j] = 1;
        buckets[1][idx].token[i] = 0;
#ifdef COUNTING
        level_item_num[0]++;
        level_item_num[1]--;
#endif
        if ((idx / locksize) != (s_idx / locksize))
          pthread_rwlock_unlock(&mutex[s_idx / locksize]);
        return i;
      }
      if ((idx / locksize) != (s_idx / locksize))
        pthread_rwlock_unlock(&mutex[s_idx / locksize]);
    }
  }
  return -1;
}

template <class T>
Value_t LevelHashing<T>::Get(T key) {
  RETRY:
  while (resizing == true) {
    asm("nop");
  }
  uint64_t f_hash = F_HASH(key);
  uint64_t s_hash = S_HASH(key);
  uint32_t f_idx = F_IDX(f_hash, addr_capacity);
  uint32_t s_idx = S_IDX(s_hash, addr_capacity);
  int i = 0, j;

  for (i = 0; i < 2; i++) {
      for (j = 0; j < ASSOC_NUM; j++) {
        spin_lock(&this->level_locks[i][f_idx].s_lock[j]);
        if (buckets[i][f_idx].token[j] == 1 &&
            buckets[i][f_idx].slot[j].key == key) {
          Value_t ret = buckets[i][f_idx].slot[j].value;
          spin_unlock(&this->level_locks[i][f_idx].s_lock[j]);
          return ret;
        }
        spin_unlock(&this->level_locks[i][f_idx].s_lock[j]);
      }    
    
      for (j = 0; j < ASSOC_NUM; j++) {
          spin_lock(&this->level_locks[i][s_idx].s_lock[j]);
          if (buckets[i][s_idx].token[j] == 1 &&
              buckets[i][s_idx].slot[j].key == key) {
            Value_t ret = buckets[i][s_idx].slot[j].value;
            spin_unlock(&this->level_locks[i][s_idx].s_lock[j]);
            return ret;
          }
          spin_unlock(&this->level_locks[i][s_idx].s_lock[j]);
      }
    
    f_idx = F_IDX(f_hash, addr_capacity / 2);
    s_idx = S_IDX(s_hash, addr_capacity / 2);
  }
  return NONE;
}

template <class T>
bool LevelHashing<T>::Delete(T key) {
#ifdef COUNTING
        level_item_num[0]--;
#endif
  RETRY:
  while (resizing == true) {
    asm("nop");
  }
  uint64_t f_hash = F_HASH(key);
  uint64_t s_hash = S_HASH(key);
  uint32_t f_idx = F_IDX(f_hash, addr_capacity);
  uint32_t s_idx = S_IDX(s_hash, addr_capacity);
  int i = 0, j;

  for (i = 0; i < 2; i++) {
      for (j = 0; j < ASSOC_NUM; j++) {
          spin_lock(&this->level_locks[i][f_idx].s_lock[j]);
          if (buckets[i][f_idx].token[j] == 1 &&
              buckets[i][f_idx].slot[j].key == key) {
            buckets[i][f_idx].token[j] = 0;
            spin_unlock(&this->level_locks[i][f_idx].s_lock[j]);
            return true;
          }
          spin_unlock(&this->level_locks[i][f_idx].s_lock[j]);
      }
      for (j = 0; j < ASSOC_NUM; j++) {
          spin_lock(&this->level_locks[i][s_idx].s_lock[j]);
          if (buckets[i][s_idx].token[j] == 1 &&
              buckets[i][s_idx].slot[j].key == key) {
            buckets[i][s_idx].token[j] = 0;
            spin_unlock(&this->level_locks[i][s_idx].s_lock[j]);
            return true;
          }
          spin_unlock(&this->level_locks[i][s_idx].s_lock[j]);
      }
    f_idx = F_IDX(f_hash, addr_capacity / 2);
    s_idx = S_IDX(s_hash, addr_capacity / 2);
  }
  return false;
}
}  // namespace level
#endif  // LEVEL_HASHING_H_