// Copyright (c) 2021-2023 Institute of Computing Technology, Chinese Academy of Sciences
// Pea Hash is licensed under Mulan PSL v2.

#pragma once

#include <hash.h>
#include <pair.h>
#include <pref.h>
#include <utils.h>

#include "space_manager.h"

namespace pea {

#ifdef AU
uint64_t space_num = 0;
uint64_t entry_num = 0;
uint64_t segment_num = 0;
#endif

// 2 ^ 14 * (2 ^ 3 + 1)
#define LOCK_TABLE_SIZE 147456
#define LOCK_SIZE_PER_SEG 9


// 31:16 version; 15:8 ref-count; 7: lock bit; 6:0 thread_id;
struct lock_type {
  uint32_t thread_id : 8;
  uint32_t ref_count : 8;
  uint32_t version : 16;
};
typedef struct lock_type lock_type;
lock_type lock_map_table[LOCK_TABLE_SIZE];

static uint64_t thread_num;
static uint64_t all_approve;
uint8_t *thread_flag;
const uint8_t RUNNING = 0;
const uint8_t PROPOSE = 1;
const uint8_t APPROVE = 2;

// 31: pointer to Pea; 30: lock bit; 29:0 version
uint32_t dir_lock;

const uint32_t dir_lock_set = (1U << 30U);
constexpr size_t kNumBucket = 2048;
constexpr size_t segment_size = kNumBucket * 256;

template <class T>
struct _Pair {
  T key;
  Value_t value;
};

struct header_t {           // 16B
  uint16_t r_bitmap;        // [15:1]bitmap, 0 reserved
  uint8_t fingerprint[14];  // Compare former 14 fingerprints,
  // last entry compared by key
};

template <class T>
struct _line {  // 64B
  _Pair<T> entry[4];
};

template <class T>
union Bucket {
  _line<T> line[4];
  _Pair<T> entry[16];
  header_t header;
};

// 8B
struct seg_meta {
  uint8_t local_depth;
  uint8_t policy;
  uint16_t reserve;
  std::atomic<uint32_t> chunk_id;

  void copyFrom(seg_meta &other) {
    local_depth = other.local_depth;
    policy = other.policy;
    chunk_id.store(other.chunk_id.load(std::memory_order_acquire), std::memory_order_release);
  }
};

struct Directory {
  seg_meta _[0];
};

constexpr uint8_t SINGLE_HASH = 0;
constexpr uint8_t DOUBLE_HASH = 1;
constexpr uint8_t STASH = 2;
constexpr uint8_t POLICY_COUNT = 3;

template <class T>
struct function {
  typedef Bucket<T> bucket_t;
  // return -3: Duplicate key; -1: full; 0: success
  int (*insert)(T key, Value_t value, uint64_t hash, bucket_t *first_bucket,
                uint32_t chunk_id, seg_meta* segment, bool &is_dir_locked);
  int (*opt_insert)(T key, Value_t value, uint64_t hash, bucket_t *first_bucket,
                uint32_t chunk_id, seg_meta* segment, uint32_t old_dir_lock_value, uint32_t old_seg_lock_value);
  /*if success, then return 0, else return -1*/
  int (*del)(T key, uint64_t hash, bucket_t *first_bucket, uint32_t chunk_id);
  Value_t (*lookup)(T key, uint64_t hash, bucket_t *first_bucket,
                    uint32_t chunk_id);
  bool (*redistribute)(bucket_t *first_bucket, uint8_t run);  // unused
};

static inline uint32_t lock_hash(uint32_t chunk_idx, uint32_t bucket_idx) {
  uint32_t retval = (chunk_idx * LOCK_SIZE_PER_SEG) + (bucket_idx & (LOCK_SIZE_PER_SEG - 2)) + 1;
  return retval;
}

static inline uint32_t seg_lock_hash(uint32_t seg_idx) {
  uint32_t retval = (seg_idx * LOCK_SIZE_PER_SEG);
  return retval;
}

class hashPolicy {
 private:
  //  static constexpr auto hasher = std::hash<std::thread::id>();

 public:
  hashPolicy() = default;
  ~hashPolicy() = default;
  //---------------------- fingerprint --------------

  // the highest 8 bits represent fingerprint
  static uint8_t inline FINGERPRINT(uint64_t hash) { return hash >> 56U; }

  //--------------------bucket index 0---------------------------
  /**
   * generate the 1st index of a key
   */
  static uint32_t inline bucket_idx0(uint64_t hash) {
    return (hash << 8U) >> 53U;
  }

  /**
   * generate the 2nd index of a key
   */
  static uint32_t inline bucket_idx1(uint64_t hash) {
    return (hash << 19U) >> 53U;
  }

  //--------------------------------- bitmap ---------------------------
  /*
   * check if the corresponding bit in bitmap is set 1
   */
  static uint16_t inline get_bit_true(uint16_t bitmap, uint8_t slot) {
    return bitmap & (1U << slot);
  }

  /*
   * set the bit in bitmap false
   */
  static uint16_t inline bitmap_set0(uint16_t bitmap, uint8_t slot) {
    return bitmap & ~(1U << slot);
  }

  /*
   * set the bit in bitmap true
   */
  static uint16_t inline bitmap_set1(uint16_t bitmap, uint8_t slot) {
    return bitmap | (1U << slot);
  }

  static uint8_t inline first_empty_slot(uint16_t bitmap) {
    uint32_t mask = ~((uint32_t)bitmap | 0x1U);
    return _bit_scan_forward((int)mask);
  }

  // the reserved bit should be unset, otherwise maybe undefined
  static uint8_t inline last_empty_slot(uint16_t bitmap) {
    return _bit_scan_reverse((int)(uint16_t)~bitmap);
  }

  // ------------------ Bucket Lock -------------------------
  static inline bool try_acquire_lock(uint32_t lock_id) {
    uint32_t *lock = (uint32_t *)&lock_map_table[lock_id];
    uint32_t old_v = __atomic_load_n(lock, __ATOMIC_ACQUIRE);
    lock_type old_value = *((lock_type *)&old_v);
    if (old_value.ref_count > 0 && (old_value.thread_id != thread_id)) {
      return false;
    }
    lock_type new_value = old_value;
    new_value.thread_id = thread_id;
    new_value.ref_count ++;
    uint32_t new_v = *((uint32_t *)&new_value);
    return CAS(lock, (uint32_t *)&old_value, new_v);
  }

  // Spin wait until get lock
  static inline void get_spin_lock(uint32_t lock_id) {
    while(try_acquire_lock(lock_id) == false);
  }

  static inline void acquire_lock(uint32_t id) {
    get_spin_lock(id);
  }

  static inline void acquire_2locks(uint32_t lock_id0, uint32_t lock_id1) {
    if (lock_id1 > lock_id0) {
      get_spin_lock(lock_id0);
      get_spin_lock(lock_id1);
    } else {
      get_spin_lock(lock_id1);
      get_spin_lock(lock_id0);
    }
  }

  // increment version and release the lock
  static inline void release_lock(uint32_t lock_id) {
    uint32_t *lock = (uint32_t *)&lock_map_table[lock_id];
    uint32_t value = __atomic_load_n(lock, __ATOMIC_ACQUIRE);
    lock_type v = *((lock_type *)&value);
    if (v.ref_count > 0){
      v.ref_count--;
      v.version ++;
      uint32_t new_v = *((uint32_t *)&v);
      __atomic_store_n(lock, new_v, __ATOMIC_RELEASE);
    }
  }

  static inline uint32_t get_lock_word(uint32_t lock_id) {
    uint32_t *lock = (uint32_t *)&lock_map_table[lock_id];
    uint32_t value = __atomic_load_n(lock, __ATOMIC_ACQUIRE);
    return value;
  }

  static inline bool is_locked(uint32_t *v_ptr) {
    lock_type v = *((lock_type *)v_ptr);
    return (v.ref_count > 0);
  }

  // useless
  static inline bool is_locked_by_others(uint32_t lock_id) {
    uint32_t *lock = (uint32_t *)&lock_map_table[lock_id];
    uint32_t value = __atomic_load_n(lock, __ATOMIC_ACQUIRE);
    lock_type v = *((lock_type *)&value);
    return (v.ref_count > 0 && (v.thread_id != thread_id));
  }


  static inline void release2locks(uint32_t lock_id0, uint32_t lock_id1) {
    if (lock_id1 > lock_id0) {
      release_lock(lock_id1);
      release_lock(lock_id0);
    } else {
      release_lock(lock_id0);
      release_lock(lock_id1);
    }
  }

};

// support variable length key
template <class T>
uint64_t hash(T key) {
  return h(&key, sizeof(key));
}

template <class T>
class PeaHashing : public Hash<T> {
 public:
  PeaHashing();
  PeaHashing(size_t, uint64_t);
  ~PeaHashing();
  inline int Insert(T key, Value_t value);
  inline int PessimisticInsert(uint64_t key_hash, T key, Value_t value);
  int Insert(T key, Value_t value, bool);
  inline bool Delete(T);
  bool Delete(T, bool);
  inline Value_t Get(T);
  Value_t Get(T key, bool is_in_epoch);
  int Get(T key, Value_t *ret_val, size_t *number) { return 0; }
  void Recovery();
  uint64_t getNumber();

  void Directory_Update(uint32_t, seg_meta *, uint8_t);

  void Directory_Doubling(uint32_t, seg_meta *);

 private:
  void init_function_list();

  // g bits in the bigger endian of [31: 32 - global_depth]hash
  // then the expand bit is in the little endian
  static uint32_t inline CHUNK_IDX(uint64_t hash, uint8_t mask_g) {
    return (hash << 32U) >> mask_g;
  }

  // bit test
  static uint64_t inline expand_bit(uint64_t hash, uint8_t new_local_depth) {
    return hash & (1UL << (32UL - new_local_depth));
  }

  static uint64_t inline expand_bits(uint64_t hash, uint8_t new_depth) {
    return hash << (30U + new_depth) >> 62U;
  }

  void free4chunk(seg_meta *segments) {
    for (int i = 0; i < 4; ++i) {
      SM->freeChunk(segments[i].chunk_id.load(std::memory_order_acquire));
    }
  }

  uint8_t *global_depth;

  Directory *dir;

  SpaceManager *SM;

  function<T> class_list[POLICY_COUNT];

  seg_meta *split4(Bucket<T> *this_bucket, Bucket<T> *, uint8_t local_depth,
                   uint32_t chunk_id);

  /**
   * @param old_bucket rehashed bucket
   * @param new_depth depth of new segments
   * @param seg_ptrs array of four segments
   */
  static inline void bucket_rehash2segments_naive(Bucket<T> *old_bucket,
                                                  uint8_t new_depth,
                                                  Bucket<T> **seg_ptrs,
                                                  uint32_t i) {
    T old_key;
    uint64_t key_hash;
    Bucket<T> *new_bucket;
    _Pair<T> *old_entry;
    uint16_t new_bitmap;
    header_t *old_hdr = &old_bucket->header;
    uint16_t old_bitmap = old_hdr->r_bitmap;
    for (uint8_t j = 15; j > 0; --j) {
      old_entry = old_bucket->entry + j;
      old_key = old_entry->key;
      key_hash = hash(old_key);
      if (hashPolicy::get_bit_true(old_bitmap, j)) {
        new_bucket = seg_ptrs[expand_bits(key_hash, new_depth)] + i;
        new_bitmap = new_bucket->header.r_bitmap;
        new_bucket->entry[j] = *old_entry;
        if (j != 15)
          new_bucket->header.fingerprint[j - 1] = old_hdr->fingerprint[j - 1];
        new_bucket->header.r_bitmap = hashPolicy::bitmap_set1(new_bitmap, j);
      }
    }
  }

  static inline void stash_rehash2stashes_naive(Bucket<T> *old_bucket,
                                                uint8_t new_depth,
                                                Bucket<T> **buc_ptrs) {
    T old_key;
    uint64_t key_hash;
    Bucket<T> *new_bucket;
    _Pair<T> *old_entry;
    uint16_t new_bitmap;
    header_t *old_hdr = &old_bucket->header;
    uint16_t old_bitmap = old_hdr->r_bitmap;
    for (uint8_t j = 15; j > 0; --j) {
      old_entry = old_bucket->entry + j;
      old_key = old_entry->key;
      key_hash = hash(old_key);
      if (hashPolicy::get_bit_true(old_bitmap, j)) {
        new_bucket = buc_ptrs[expand_bits(key_hash, new_depth)];
        new_bitmap = new_bucket->header.r_bitmap;
        new_bucket->entry[j] = *old_entry;
        if (j != 15)
          new_bucket->header.fingerprint[j - 1] = old_hdr->fingerprint[j - 1];
        new_bucket->header.r_bitmap = hashPolicy::bitmap_set1(new_bitmap, j);
      }
    }
  }

  /**
   * @return true: success; false: fail
   */
  static inline bool bucket_rehash2segments_single(Bucket<T> *old_bucket,
                                                   uint8_t new_depth,
                                                   Bucket<T> **seg_ptrs) {
    T old_key;
    uint64_t key_hash;
    uint32_t bucket_id0;
    Bucket<T> *new_bucket;
    _Pair<T> *old_entry;
    header_t *old_hdr = &old_bucket->header;
    uint16_t old_bitmap = old_hdr->r_bitmap;
    for (uint8_t k = 15; k > 0; --k) {
      if (hashPolicy::get_bit_true(old_bitmap, k)) {
        old_entry = old_bucket->entry + k;
        old_key = old_entry->key;
        key_hash = hash(old_key);
        bucket_id0 = hashPolicy::bucket_idx0(key_hash);
        new_bucket = seg_ptrs[expand_bits(key_hash, new_depth)] + bucket_id0;
        uint16_t new_bitmap = new_bucket->header.r_bitmap;
        uint8_t entry_id = hashPolicy::last_empty_slot(new_bitmap | 0x8000U);
        if (entry_id == 0) {
          entry_id = hashPolicy::last_empty_slot(new_bitmap);
        }
        switch (entry_id) {
          case 0:
            return false;
          case 15:
            break;
          default:
            uint8_t fgprt = (k == 15) ? hashPolicy::FINGERPRINT(key_hash)
                                      : old_hdr->fingerprint[k - 1];
            new_bucket->header.fingerprint[entry_id - 1] = fgprt;
            break;
        }
        new_bucket->entry[entry_id] = *old_entry;
        // set bitmap
        new_bucket->header.r_bitmap =
            hashPolicy::bitmap_set1(new_bitmap, entry_id);
      }
    }
    return true;
  }

  static inline bool bucket_rehash2segments_double(Bucket<T> *old_bucket,
                                                   uint8_t new_depth,
                                                   Bucket<T> **seg_ptrs) {
    T old_key;
    uint64_t key_hash;
    uint32_t bucket_id0, bucket_id1;
    Bucket<T> *new_bucket, *first_bucket, *bucket0, *bucket1;
    _Pair<T> *old_entry;
    uint16_t bitmap0, bitmap1, new_bitmap;
    header_t *old_hdr = &old_bucket->header;
    uint16_t old_bitmap = old_hdr->r_bitmap;
    for (uint8_t k = 15; k > 0; --k) {
      if (hashPolicy::get_bit_true(old_bitmap, k)) {
        old_entry = old_bucket->entry + k;
        old_key = old_entry->key;
        key_hash = hash(old_key);
        // choose the one with fewer entries
        bucket_id0 = hashPolicy::bucket_idx0(key_hash);
        bucket_id1 = hashPolicy::bucket_idx1(key_hash);
        first_bucket = seg_ptrs[expand_bits(key_hash, new_depth)];
        bucket0 = first_bucket + bucket_id0;
        bucket1 = first_bucket + bucket_id1;
        bitmap0 = bucket0->header.r_bitmap;
        bitmap1 = bucket1->header.r_bitmap;
        // @warning the function do not consider the reserve bit
        int count0 = _popcnt32(bitmap0), count1 = _popcnt32(bitmap1);
        if (count0 < count1) {
          new_bucket = bucket0;
          new_bitmap = bitmap0;
        } else {
          new_bucket = bucket1;
          new_bitmap = bitmap1;
        }
        uint8_t entry_id = hashPolicy::last_empty_slot(new_bitmap | 0x8000U);
        if (entry_id == 0) {
          entry_id = hashPolicy::last_empty_slot(new_bitmap);
        }
        switch (entry_id) {
          case 0:
            return false;
          case 15:
            break;
          default:
            uint8_t fgprt = (k == 15) ? hashPolicy::FINGERPRINT(key_hash)
                                      : old_hdr->fingerprint[k - 1];
            new_bucket->header.fingerprint[entry_id - 1] = fgprt;
            break;
        }
        new_bucket->entry[entry_id] = *old_entry;
        // set bitmap
        new_bucket->header.r_bitmap =
            hashPolicy::bitmap_set1(new_bitmap, entry_id);
      }
    }
    return true;
  }

};

static inline bool try_get_dir_lock() {
  uint32_t old_value = __atomic_load_n(&dir_lock, __ATOMIC_ACQUIRE);
  if (old_value & dir_lock_set) {
    return false;
  }
  uint32_t new_value = old_value | dir_lock_set;
  return CAS(&dir_lock, &old_value, new_value);
}

static inline void acquire_dir_lock() {
  while (try_get_dir_lock() == false);
}

// get lock word
static inline uint32_t get_dir_lock() {
  return __atomic_load_n(&dir_lock, __ATOMIC_ACQUIRE);
}

// release lock
static inline void release_dir_lock() {
  uint32_t v = (dir_lock + 1) & 0xbfffffffU;
  __atomic_store_n(&dir_lock, v, __ATOMIC_RELEASE);
}

// release and flip lock
static inline void release_flip_lock() {
  uint32_t v = (dir_lock + 1) & 0xbfffffffU;
  asm volatile("btc $31, %1" : "+r"(v) :);
  __atomic_store_n(&dir_lock, v, __ATOMIC_RELEASE);
}

template <class T>
PeaHashing<T>::PeaHashing(size_t initCap, uint64_t t_num) {
  // @warning: InitCap should be exponential of 2
  if (initCap < 2) {
    LOG_FATAL("[ERROR] The initCap must be larger than 1");
  }

  SM = SpaceManager::Get();

  global_depth = reinterpret_cast<uint8_t *>(SM->meta0);
  dir = reinterpret_cast<Directory *>(SM->dir0);
  uint8_t G = static_cast<size_t>(log2(initCap));
#ifdef AU
  segment_num += initCap;
#endif
  /* Initialize the Directory*/
  for (int i = initCap - 1; i >= 0; --i) {
    uint32_t chunk_idx = SM->allocChunk();
    auto seg_begin = reinterpret_cast<Bucket<T> *>(SM->chunk_addr(chunk_idx));
    for (uint64_t j = 0; j < kNumBucket; ++j) {
      Bucket<T> *curr_bucket = seg_begin + j;
      memset(curr_bucket, 0, 2);
      
    }
    dir->_[i].local_depth = G;
    dir->_[i].policy = SINGLE_HASH;
    dir->_[i].chunk_id = chunk_idx;
  }
  
  *global_depth = G;

  dir_lock = 0;
  memset(lock_map_table, 0, LOCK_TABLE_SIZE * 4);
  thread_num = t_num;

  init_function_list();
}

template <class T>
PeaHashing<T>::~PeaHashing<T>() {
  delete[] thread_flag;
  std::cout << "Destructor." << std::endl;
}

template <class T>
int PeaHashing<T>::Insert(T key, Value_t value) {
  uint64_t key_hash = hash(key);
RETRY:
  uint32_t old_dir_lock_value = __atomic_load_n(&dir_lock, __ATOMIC_ACQUIRE);
  if (old_dir_lock_value & dir_lock_set) {
    goto RETRY;
  }
  uint8_t cur_g = *global_depth;                                               
  uint32_t idx = CHUNK_IDX(key_hash, 64 - cur_g);
  seg_meta* target_meta = dir->_ + idx;
  // attention visit w/o lock
  uint32_t chunk_idx = target_meta->chunk_id.load(std::memory_order_acquire);
  uint32_t seg_lock_id = seg_lock_hash(chunk_idx);
  uint32_t seg_lock_word = hashPolicy::get_lock_word(seg_lock_id);
  if (hashPolicy::is_locked(&seg_lock_word)) goto RETRY;

  uint8_t policy = target_meta->policy;
  /*if (L < cur_g || policy != STASH) {
    is_dir_locked = false;
    if (old_dir_lock_value != __atomic_load_n(&dir_lock, __ATOMIC_ACQUIRE)){
      hashPolicy::release_lock(seg_lock_hash(chunk_idx));
      goto RETRY;
    } TOO VAGUE
  }*/
  auto seg_begin = reinterpret_cast<Bucket<T> *>(SM->chunk_addr(chunk_idx));
  // insert to segment
  int ret = class_list[policy].opt_insert(key, value, key_hash, seg_begin,
                                      chunk_idx, target_meta, old_dir_lock_value, seg_lock_word);
  switch (ret) {
  case -4:
    goto RETRY;
  case -1:
    return PessimisticInsert(key_hash, key, value);
  default:
    return ret;
  }
}


/**
 * @return -3: Duplicate key exist!; 0: Success.
 */
template <class T>
int PeaHashing<T>::PessimisticInsert(uint64_t key_hash, T key, Value_t value) {
RETRY:
  bool is_dir_locked = false;
  acquire_dir_lock();
  is_dir_locked = true;
  uint8_t cur_g = *global_depth;                                               
  uint32_t idx = CHUNK_IDX(key_hash, 64 - cur_g);
  seg_meta* target_meta = dir->_ + idx;
  // attention visit w/o lock
  uint32_t chunk_idx = target_meta->chunk_id.load(std::memory_order_acquire);
  if (hashPolicy::try_acquire_lock(seg_lock_hash(chunk_idx)) == false) {
    release_dir_lock();
    goto RETRY;
  }
  uint8_t policy = target_meta->policy;
  uint8_t L = target_meta->local_depth;
  if (L < cur_g || policy != STASH) {
    release_dir_lock();
    is_dir_locked = false;
  }
  if (chunk_idx != target_meta->chunk_id.load(std::memory_order_acquire)) {
    LOG("[INFO] Insert seg meta changed");
    if (is_dir_locked) {
      release_dir_lock();
    }
    goto RETRY;    
  }
  auto seg_begin = reinterpret_cast<Bucket<T> *>(SM->chunk_addr(chunk_idx));
  // insert to segment
  int ret = class_list[policy].insert(key, value, key_hash, seg_begin,
                                      chunk_idx, target_meta, is_dir_locked);
  if (ret == -1) {
    goto RETRY;
  } else if (ret == -2) {
    // ret = -2 stash insert failed

    // 1. allocate 2 new chunks, redistribute in it
    uint32_t stash_c_id = SM->get_stash_cid(chunk_idx);
    auto stash_chunk =
        reinterpret_cast<Bucket<T> *>(SM->chunk_addr(stash_c_id));
    uint32_t stash_b_id = (chunk_idx << 1U) % kNumBucket;
    Bucket<T> *oldSB0 = stash_chunk + stash_b_id;
    seg_meta *new_segments = split4(seg_begin, oldSB0, L, chunk_idx);
    if (is_dir_locked) {
      // Directory doubling and change hash pointer to the new
      // directory(Atomic Write)
      // get every segment lock to ensure all segment are updated
      // the most simple way is to add lock in each segment lock in the hash table
      // it will lock 2 times in current seg
      for (size_t i = 0; i < LOCK_TABLE_SIZE; i += LOCK_SIZE_PER_SEG) {
        hashPolicy::acquire_lock(i);
      }
      Directory_Doubling(idx, new_segments);
      // free all seg locks
      for (int32_t i = LOCK_TABLE_SIZE - LOCK_SIZE_PER_SEG; i >= 0; i -= LOCK_SIZE_PER_SEG) {
        hashPolicy::release_lock(i);
      }

    } else {
      // Directory updating and new seg_lock is still acquired
      Directory_Update(idx, new_segments, L);
    }
    // free bucket locks in the segment; lock in this segment should not interfere lock in other segments
    uint32_t last_lock_id = chunk_idx * LOCK_SIZE_PER_SEG + LOCK_SIZE_PER_SEG - 1;
    for (size_t i = 0; i < LOCK_SIZE_PER_SEG - 1; i++) {
      hashPolicy::release_lock(last_lock_id - i);
    }
    hashPolicy::release_lock(seg_lock_hash(chunk_idx));
    if (is_dir_locked){
      release_flip_lock();
    }
    // 3. free old chunk, free new_segments pair
    // crash: redo from 3
    delete[] new_segments;
    _mm_stream_si64(reinterpret_cast<long long int *>(oldSB0), 0);
    _mm_stream_si64(reinterpret_cast<long long int *>(oldSB0 + 1), 0);
    SM->freeChunk(chunk_idx);
    goto RETRY;
  } else {
    return ret;
  }
}

template <class T>
seg_meta *PeaHashing<T>::split4(Bucket<T> *this_bucket, Bucket<T> *oldSB0,
                                uint8_t local_depth, uint32_t chunk_id) {
  // Allocate 4 segments
  uint32_t i;
  uint8_t j, k;
  uint8_t new_depth = local_depth + 2;

  auto *ret_val = new seg_meta[4];
  auto **seg_ptrs = new Bucket<T> *[4];

  uint32_t c_ids[4];
  SM->alloc4chunk(c_ids);
  for (j = 0; j < 4; j++) {
    ret_val[j].local_depth = new_depth;
    ret_val[j].policy = SINGLE_HASH;
    ret_val[j].chunk_id.store(c_ids[j], std::memory_order_release);
    seg_ptrs[j] = (Bucket<T> *)malloc(segment_size);
    for (i = 0; i < kNumBucket; ++i) {
      Bucket<T> *cur_bucket = seg_ptrs[j] + i;
      memset(cur_bucket, 0, 2);
    }
  }

  // Stash bucket also rehash into corresponding ones, if fail, allocate stash
  // buckets
  Bucket<T> *oldSB1 = oldSB0 + 1;
  bucket_rehash2segments_single(oldSB0, new_depth, seg_ptrs);
  if (!bucket_rehash2segments_single(oldSB1, new_depth, seg_ptrs)) {
    goto DOUBLE_REHASH;
  }
  for (i = 0; i < kNumBucket; ++i) {
    Bucket<T> *old_bucket = this_bucket + i;
    if (!bucket_rehash2segments_single(old_bucket, new_depth, seg_ptrs)) {
      goto DOUBLE_REHASH;
    }
  }
  // success
  bool suc0, suc1;
  goto FINAL;

DOUBLE_REHASH:
  // clear segments
  for (j = 0; j < 4; j++) {
    ret_val[j].policy = DOUBLE_HASH;
    for (i = 0; i < kNumBucket; ++i) {
      Bucket<T> *cur_bucket = seg_ptrs[j] + i;
      memset(cur_bucket, 0, 2);
    }
  }
  // double_rehash
  for (i = 0; i < kNumBucket; ++i) {
    Bucket<T> *old_bucket = this_bucket + i;
    bucket_rehash2segments_naive(old_bucket, new_depth, seg_ptrs, i);
  }
  // stash rehash to double segments
  suc0 = bucket_rehash2segments_double(oldSB0, new_depth, seg_ptrs);
  suc1 = bucket_rehash2segments_double(oldSB1, new_depth, seg_ptrs);
  if (!(suc0 & suc1)) {
    // STASH REHASH
    auto **new_stash_buckets0 = new Bucket<T> *[4];
    auto **new_stash_buckets1 = new Bucket<T> *[4];
    for (j = 0; j < 4; ++j) {
      ret_val[j].policy = STASH;
      uint32_t new_chunk_id = ret_val[j].chunk_id.load(std::memory_order_acquire);
      SM->allocStashChunkId(new_chunk_id);
      uint32_t new_stash_c_id = SM->stash_chunk_id_array[new_chunk_id / stashUnitsPChunk];
      auto new_stash_chunk =
          reinterpret_cast<Bucket<T> *>(SM->chunk_addr(new_stash_c_id));
      uint32_t new_stash_b_id = (new_chunk_id << 1U) % kNumBucket;
      new_stash_buckets0[j] = new_stash_chunk + new_stash_b_id;
      new_stash_buckets1[j] = new_stash_buckets0[j] + 1;
    }
    stash_rehash2stashes_naive(oldSB0, new_depth, new_stash_buckets0);
    stash_rehash2stashes_naive(oldSB1, new_depth, new_stash_buckets1);
    delete[] new_stash_buckets1;
    delete[] new_stash_buckets0;
  }

FINAL:

  for (j = 0; j < 4; j++) {
    Bucket<T> *from_seg = seg_ptrs[j];
    auto to_seg =
        reinterpret_cast<Bucket<T> *>(SM->chunk_addr(ret_val[j].chunk_id.load(std::memory_order_acquire)));
#ifdef AVX512F
    for (i = 0; i < kNumBucket; ++i) {
      _line<T> *from_line = from_seg[i].line;
      _line<T> *to_line = to_seg[i].line;
      for (k = 0; k < 4; ++k) {
        __m512i new_line0 = _mm512_loadu_si512(from_line + k);
        _mm512_stream_si512((__m512i *)(to_line + k), new_line0);
      }
    }
#else
    memcpy(to_seg, from_seg, kNumBucket * 256);
#endif
  }

  for (j = 0; j < 4; ++j) {
    free(seg_ptrs[j]);
  }
  delete[] seg_ptrs;
  return ret_val;
}

template <class T>
void PeaHashing<T>::Directory_Update(uint32_t idx, seg_meta *new_s,
                                     uint8_t old_local_d) {
  //  std::cout << "Dir update at index: " << idx << std::endl;
  seg_meta *dir_entry = dir->_;
  uint8_t depth_diff = *global_depth - old_local_d;
  uint32_t table00_start_idx = idx >> depth_diff << depth_diff;
  uint32_t old_chunk_size = 1UL << depth_diff;
  uint32_t new_chunk_size = old_chunk_size >> 2U;
  uint32_t other_ids = table00_start_idx + 1;
  for (uint32_t i = other_ids; i < table00_start_idx + new_chunk_size; ++i) {
    dir_entry[i].copyFrom(new_s[0]);
  }
  for (int i = 1; i < 4; ++i) {
    uint32_t start_idx = table00_start_idx + new_chunk_size * i;
    for (uint32_t j = start_idx; j < start_idx + new_chunk_size; ++j) {
      dir_entry[j].copyFrom(new_s[i]);
    }
  }
  dir_entry[table00_start_idx].copyFrom(new_s[0]);
  // for (int i = 3; i >= 0; i--) {
  //   hashPolicy::release_lock(seg_lock_hash(new_s[i].chunk_id));
  // }
}

template <class T>
void PeaHashing<T>::Directory_Doubling(uint32_t idx, seg_meta *new_s) {
  uint8_t old_G = *global_depth;
  // std::cout << "Dir to " << (uint32_t)old_G + 2
  //           << std::endl;
  seg_meta *d = dir->_;
  size_t old_capacity = 1UL << old_G;
  size_t new_capacity = old_capacity << 2U;
  // allocate dir space and modify it

  Directory *new_d;
  if (dir_lock & 0x80000000U) {
    new_d = reinterpret_cast<Directory *>(SM->dir0);
    global_depth = reinterpret_cast<uint8_t *>(SM->meta0);
  } else {
    new_d = reinterpret_cast<Directory *>(SM->dir1);
    global_depth = reinterpret_cast<uint8_t *>(SM->meta1);
  }

  for (uint32_t i = 0; i < old_capacity; ++i) {
    uint32_t new_i = i << 2U;
    for (int j = 0; j < 4; ++j) {
      new_d->_[new_i + j].copyFrom(d[i]);
    }
  }

  uint32_t idx0 = idx << 2U;
  for (int i = 0; i < 4; ++i) {
    new_d->_[idx0 + i].copyFrom(new_s[i]);
  } 
  *global_depth = old_G + 2;
  dir = new_d;
}

// TODO (this delete is mocking optimistic concurrency)
template <class T>
bool PeaHashing<T>::Delete(T key) {
  uint64_t key_hash = hash(key);
  int ret_val;
  uint32_t idx = CHUNK_IDX(key_hash, 64 - *global_depth);   //TODO atomic
  seg_meta* target_meta = dir->_ + idx;
  uint32_t chunk_idx = target_meta->chunk_id.load(std::memory_order_relaxed);

  auto seg_begin = reinterpret_cast<Bucket<T> *>(SM->chunk_addr(chunk_idx));
  ret_val = class_list[target_meta->policy].del(key, key_hash, seg_begin,
                                               chunk_idx);
  return (ret_val == 0);
}

template <class T>
Value_t PeaHashing<T>::Get(T key) {
  uint64_t key_hash = hash(key);
RETRY:
  uint32_t old_dir_lock_value = __atomic_load_n(&dir_lock, __ATOMIC_ACQUIRE); //TODO denote split
  if (old_dir_lock_value & dir_lock_set) {
    goto RETRY;
  }
  Value_t ret_val;
  uint32_t idx = CHUNK_IDX(key_hash, 64 - *global_depth);   //TODO atomic
  seg_meta* target_meta = dir->_ + idx;
  uint32_t chunk_idx = target_meta->chunk_id.load(std::memory_order_relaxed);

  auto seg_begin = reinterpret_cast<Bucket<T> *>(SM->chunk_addr(chunk_idx));
  ret_val = class_list[target_meta->policy].lookup(key, key_hash, seg_begin,
                                                  chunk_idx);
  if (old_dir_lock_value != __atomic_load_n(&dir_lock, __ATOMIC_ACQUIRE)) {
    goto RETRY;
  }
  return ret_val;
}

/*
 * According the line bitmap, which slot is empty in the line.
 * line bitmap is the binary representation of row id.
 * the first column is how many empty slots in the line (if it > 3, set 3)
 * the other three number is slot index in the line. use 255 as padding.
 */
static uint8_t empty_slots_in_line[16][4] = {
    {3, 0, 1, 2},    // 0000
    {3, 1, 2, 3},    // 0001
    {3, 0, 2, 3},    // 0010
    {2, 2, 3, 255},  // 0011

    {3, 0, 1, 3},      // 0100
    {2, 1, 3, 255},    // 0101
    {2, 0, 3, 255},    // 0110
    {1, 3, 255, 255},  // 0111

    {3, 0, 1, 2},      // 1000
    {2, 1, 2, 255},    // 1001
    {2, 0, 2, 255},    // 1010
    {1, 2, 255, 255},  // 1011

    {2, 0, 1, 255},     // 1100
    {1, 1, 255, 255},   // 1101
    {1, 0, 255, 255},   // 1110
    {0, 255, 255, 255}  // 1111
};


template <class T>
class singleHash : public hashPolicy {
 public:
  typedef Bucket<T> bucket_t;
  typedef _Pair<T> pea_entry;
  singleHash() = default;
  ~singleHash() = default;

  static Value_t lookup(T key, uint64_t hash, bucket_t *first_bucket,
                        uint32_t chunk_id) {
    uint8_t j;
    // Calculate the index of bucket
    uint32_t bucket_id = bucket_idx0(hash);
    uint32_t lock_id = lock_hash(chunk_id, bucket_id);
    bucket_t *bucket = &first_bucket[bucket_id];
    uint8_t fgprt = FINGERPRINT(hash);
    uint32_t lock_word;
    Value_t ret_val;

  RE_SEARCH:
    lock_word = get_lock_word(lock_id);
    if (is_locked(&lock_word)) goto RE_SEARCH;
    pea_entry *this_entry;
    header_t *hdr = &bucket->header;
    uint16_t bitmap = hdr->r_bitmap;

    // for fingerprint 1-14 SIMD comparison
    // a. set every byte to key_hash in a 16B register
    __m128i key_16B = _mm_set1_epi8((char)fgprt);

    // b. load meta into another 16B register
    __m128i fgprt_16B = _mm_load_si128((const __m128i *)hdr);

    // c. compare them
    __m128i cmp_res = _mm_cmpeq_epi8(key_16B, fgprt_16B);

    // d. generate a mask
    auto mask = (unsigned int)_mm_movemask_epi8(cmp_res);  // 1: same; 0: diff

    // AND bitmap
    mask = mask & ((unsigned int)bitmap << 1U);

    // search every matching candidate
    while (mask) {
      int jj = bitScan(mask) - 1;  // next candidate
      this_entry = bucket->entry + jj - 1;
      if (key == this_entry->key) {
        ret_val = this_entry->value;
        if (get_lock_word(lock_id) != lock_word) goto RE_SEARCH;
        return ret_val;
      }

      mask &= ~(0x1U << (uint16_t)jj);  // remove this bit
    }                                   // end while

    if (get_bit_true(bitmap, 15)) {
      this_entry = bucket->entry + 15;
      if (key == this_entry->key) {
        ret_val = this_entry->value;
        if (get_lock_word(lock_id) != lock_word) goto RE_SEARCH;
        return ret_val;
      }
    }

    return NONE;
  }

  static inline bool check_key_exist(bucket_t *bucket, uint16_t bitmap, T key,
                                     uint8_t fgprt) {
    pea_entry *this_entry;
    header_t *hdr = &bucket->header;
    __m128i key_16B = _mm_set1_epi8((char)fgprt);
    __m128i fgprt_16B = _mm_load_si128((const __m128i *)hdr);
    __m128i cmp_res = _mm_cmpeq_epi8(key_16B, fgprt_16B);
    auto mask = (unsigned int)_mm_movemask_epi8(cmp_res);  // 1: same; 0: diff
    mask = mask & ((unsigned int)bitmap << 1U);
    while (mask) {
      int jj = bitScan(mask) - 1;  // next candidate
      this_entry = bucket->entry + jj - 1;
      if (key == this_entry->key) {
        return true;
      }
      mask &= ~(0x1U << (uint16_t)jj);  // remove this bit
    }                                   // end while

    if (get_bit_true(bitmap, 15)) {
      this_entry = bucket->entry + 15;
      if (key == this_entry->key) {
        return true;
      }
    }
    return false;
  }

  static inline void insert_into_bucket(bucket_t *bucket, T key, Value_t value, uint8_t fgprt) {
    header_t * header = &bucket->header;
    uint16_t bitmap = header->r_bitmap;
    uint8_t slot = first_empty_slot(bitmap);
    // write entry
    bucket->entry[slot] = {.key = key, .value = value};

    uint8_t lineid = slot >> 2U;
    if (lineid == 0) {
      header->fingerprint[slot - 1] = fgprt;
      // write word 0
      // flush
      bitmap = bitmap_set1(bitmap, slot);
    } else {
      if (slot != 15) {
        header->fingerprint[slot - 1] = fgprt;
      }
      bitmap = bitmap_set1(bitmap, slot);
      uint8_t *line_empty_slots =
          empty_slots_in_line[(uint8_t)(bitmap >> (lineid * 4U)) & 0x0fU];
      for (uint8_t from = 1; from <= line_empty_slots[0]; ++from) {
        uint8_t to = line_empty_slots[from];
        // copy value
        bucket->line[lineid].entry[to] = bucket->line[0].entry[from];
        uint8_t to_slotid = lineid * 4 + to;
        if (to_slotid != 15) {
          header->fingerprint[to_slotid - 1] = header->fingerprint[from - 1];
        }
        // validate and invalidate corresponding bitmap
        bitmap = bitmap_set1(bitmap, to_slotid);
        bitmap = bitmap_set0(bitmap, from);
      }
      // after entry moving, persist the new address    
    }
    // atomic write
    header->r_bitmap = bitmap; 
  }

  /**
   * @return -3: Duplicate key exist!; -1: full; 0:
   * Success.
   * @details is_dir_locked must be false
   */
  static int insert(T key, Value_t value, uint64_t hash, bucket_t *first_bucket,
                    uint32_t chunk_id, seg_meta *target_seg, bool &is_dir_locked) {
    bucket_t *bucket;
    uint32_t bucket_id = bucket_idx0(hash);
    bucket = &first_bucket[bucket_id];
    LINE_PREF(bucket);
    uint8_t fgprt = FINGERPRINT(hash);
    uint32_t lock_id = lock_hash(chunk_id, bucket_id);
    acquire_lock(lock_id);

    uint16_t bitmap = bucket->header.r_bitmap;
    if (check_key_exist(bucket, bitmap, key, fgprt)){
      release_lock(lock_id);
      release_lock(seg_lock_hash(chunk_id));
      return -3;
    }
    if(_popcnt32(bitmap) == 15) { // full
      target_seg->policy = DOUBLE_HASH;
      release_lock(lock_id);
      release_lock(seg_lock_hash(chunk_id));
      return -1;
    } else {
      release_lock(seg_lock_hash(chunk_id));
      insert_into_bucket(bucket, key, value, fgprt);
      release_lock(lock_id);
      return 0;
    }
  }

  static int opt_insert(T key, Value_t value, uint64_t hash, bucket_t *first_bucket,
                    uint32_t chunk_id, seg_meta *target_seg, uint32_t old_dir_lock_value, uint32_t old_seg_lock_value) {
    bucket_t *bucket;
    uint32_t bucket_id = bucket_idx0(hash);
    bucket = &first_bucket[bucket_id];
    LINE_PREF(bucket);
    uint8_t fgprt = FINGERPRINT(hash);
    uint32_t lock_id = lock_hash(chunk_id, bucket_id);
    acquire_lock(lock_id);
    if (old_seg_lock_value != get_lock_word(seg_lock_hash(chunk_id))
        || chunk_id != target_seg->chunk_id.load(std::memory_order_acquire)
        || old_dir_lock_value != __atomic_load_n(&dir_lock, __ATOMIC_ACQUIRE)
        ) {
      release_lock(lock_id);
      return -4;
    }

    uint16_t bitmap = bucket->header.r_bitmap;
    if (check_key_exist(bucket, bitmap, key, fgprt)){
      release_lock(lock_id);
      return -3;
    }
    if(_popcnt32(bitmap) == 15) {
      release_lock(lock_id);
      return -1;
    } else {
      insert_into_bucket(bucket, key, value, fgprt);
      release_lock(lock_id);
      return 0;
    }
  }

  /*if success, then return 0, else return -1*/
  static int del(T key, uint64_t hash, bucket_t *first_bucket,
                 uint32_t chunk_id) {
    uint8_t j;
    bucket_t *bucket;
    uint32_t bucket_id = bucket_idx0(hash);
    bucket = &first_bucket[bucket_id];
    uint8_t fgprt = FINGERPRINT(hash);
    pea_entry *this_entry;
    uint32_t lock_id = lock_hash(chunk_id, bucket_id);
    uint32_t seg_lock_word;

    acquire_lock(lock_id);
    header_t *hdr = &bucket->header;
    uint16_t bitmap = hdr->r_bitmap;

    __m128i key_16B = _mm_set1_epi8((char)fgprt);
    __m128i fgprt_16B = _mm_load_si128((const __m128i *)hdr);
    __m128i cmp_res = _mm_cmpeq_epi8(key_16B, fgprt_16B);
    auto mask = (unsigned int)_mm_movemask_epi8(cmp_res);  // 1: same; 0: diff
    mask = mask & ((unsigned int)bitmap << 1U);
    while (mask) {
      int jj = bitScan(mask) - 1;  // next candidate
      this_entry = bucket->entry + jj - 1;
      if (key == this_entry->key) {
        hdr->r_bitmap = bitmap_set0(bitmap, jj - 1);
        release_lock(lock_id);
        return 0;
      }
      mask &= ~(0x1U << (uint16_t)jj);  // remove this bit
    }                                   // end while

    if (get_bit_true(bitmap, 15)) {
      this_entry = bucket->entry + 15;
      if (key == this_entry->key) {
        hdr->r_bitmap = bitmap_set0(bitmap, 15);
        release_lock(lock_id);
        return 0;
      }
    }
    release_lock(lock_id);
    return -1;
  }
};

template <class T>
class doubleHash : public hashPolicy {
 public:
  typedef Bucket<T> bucket_t;
  typedef _Pair<T> pea_entry;

  doubleHash() = default;

  ~doubleHash() = default;

  static Value_t inline bucket_lookup(T key, uint32_t lock_id, bucket_t *bucket,
                                      uint8_t fgprt) {
    uint8_t j;
    pea_entry *this_entry;
    Value_t ret_val;
  RE_SEARCH:
    uint32_t lock_word = get_lock_word(lock_id);
    if (is_locked(&lock_word)) goto RE_SEARCH;

    header_t *hdr = &bucket->header;
    uint16_t bitmap = hdr->r_bitmap;

    // for fingerprint 1-14 SIMD comparison
    // a. set every byte to key_hash in a 16B register
    __m128i key_16B = _mm_set1_epi8((char)fgprt);

    // b. load meta into another 16B register
    __m128i fgprt_16B = _mm_load_si128((const __m128i *)hdr);

    // c. compare them
    __m128i cmp_res = _mm_cmpeq_epi8(key_16B, fgprt_16B);

    // d. generate a mask
    auto mask = (unsigned int)_mm_movemask_epi8(cmp_res);  // 1: same; 0: diff

    // AND bitmap
    mask = mask & ((unsigned int)bitmap << 1U);

    // search every matching candidate
    while (mask) {
      int jj = bitScan(mask) - 1;  // next candidate
      this_entry = bucket->entry + jj - 1;
      if (key == this_entry->key) {
        ret_val = this_entry->value;
        if (get_lock_word(lock_id) != lock_word) goto RE_SEARCH;
        return ret_val;
      }

      mask &= ~(0x1U << (uint16_t)jj);  // remove this bit
    }                                   // end while

    if (get_bit_true(bitmap, 15)) {
      this_entry = bucket->entry + 15;
      if (key == this_entry->key) {
        ret_val = this_entry->value;
        if (get_lock_word(lock_id) != lock_word) goto RE_SEARCH;
        return ret_val;
      }
    }
    return NONE;
  }

  static int inline bucket_delete(T key, uint32_t lock_id, bucket_t *bucket,
                                  uint8_t fgprt, uint32_t chunk_id) {
    pea_entry *this_entry;
    uint8_t j;
    acquire_lock(lock_id);
    header_t *hdr = &bucket->header;
    uint16_t bitmap = hdr->r_bitmap;
    __m128i key_16B = _mm_set1_epi8((char)fgprt);
    __m128i fgprt_16B = _mm_load_si128((const __m128i *)hdr);
    __m128i cmp_res = _mm_cmpeq_epi8(key_16B, fgprt_16B);
    auto mask = (unsigned int)_mm_movemask_epi8(cmp_res);  // 1: same; 0: diff
    mask = mask & ((unsigned int)bitmap << 1U);
    while (mask) {
      int jj = bitScan(mask) - 1;  // next candidate
      this_entry = bucket->entry + jj - 1;
      if (key == this_entry->key) {
        hdr->r_bitmap = bitmap_set0(bitmap, jj - 1);        
        release_lock(lock_id);
        return 0;
      }
      mask &= ~(0x1U << (uint16_t)jj);  // remove this bit
    }

    if (get_bit_true(bitmap, 15)) {
      this_entry = bucket->entry + 15;
      if (key == this_entry->key) {
        hdr->r_bitmap = bitmap_set0(bitmap, 15);        
        release_lock(lock_id);
        return 0;
      }
    }
    release_lock(lock_id);
    return -1;
  }

  /* See hashTable.h */
  static Value_t lookup(T key, uint64_t hash, bucket_t *first_bucket,
                        uint32_t chunk_id) {
    // Calculate 2 indexes of buckets
    uint32_t bucket_id0 = bucket_idx0(hash);
    uint32_t bucket_id1 = bucket_idx1(hash);
    bucket_t *bucket0 = &first_bucket[bucket_id0];
    bucket_t *bucket1 = &first_bucket[bucket_id1];
    LINE_PREF(bucket0);
    LINE_PREF(bucket1);

    uint32_t lock_id1 = lock_hash(chunk_id, bucket_id1);
    uint32_t lock_id0 = lock_hash(chunk_id, bucket_id0);
    uint8_t fgprt = FINGERPRINT(hash);

    Value_t ret_val;
    ret_val = bucket_lookup(key, lock_id1, bucket1, fgprt);
    if (ret_val == NONE) ret_val = bucket_lookup(key, lock_id0, bucket0, fgprt);
    return ret_val;
  }

  static inline bool check_key_exist(bucket_t *bucket, uint16_t bitmap, T key,
                                     uint8_t fgprt) {
    pea_entry *this_entry;
    header_t *hdr = &bucket->header;
    __m128i key_16B = _mm_set1_epi8((char)fgprt);
    __m128i fgprt_16B = _mm_load_si128((const __m128i *)hdr);
    __m128i cmp_res = _mm_cmpeq_epi8(key_16B, fgprt_16B);
    auto mask = (unsigned int)_mm_movemask_epi8(cmp_res);  // 1: same; 0: diff
    mask = mask & ((unsigned int)bitmap << 1U);
    while (mask) {
      int jj = bitScan(mask) - 1;  // next candidate
      this_entry = bucket->entry + jj - 1;
      if (key == this_entry->key) {
        return true;
      }
      mask &= ~(0x1U << (uint16_t)jj);  // remove this bit
    }                                   // end while

    if (get_bit_true(bitmap, 15)) {
      this_entry = bucket->entry + 15;
      if (key == this_entry->key) {
        return true;
      }
    }
    return false;
  }

  static inline void insert_into_bucket(bucket_t *bucket, T key, Value_t value, uint8_t fgprt) {
    header_t * header = &bucket->header;
    uint16_t bitmap = header->r_bitmap;
    uint8_t slot = first_empty_slot(bitmap);
    // write entry
    bucket->entry[slot] = {.key = key, .value = value};

    uint8_t lineid = slot >> 2U;
    if (lineid == 0) {
      header->fingerprint[slot - 1] = fgprt;
      // write word 0
      // flush
      bitmap = bitmap_set1(bitmap, slot);
    } else {
      if (slot != 15) {
        header->fingerprint[slot - 1] = fgprt;
      }
      bitmap = bitmap_set1(bitmap, slot);
      uint8_t *line_empty_slots =
          empty_slots_in_line[(uint8_t)(bitmap >> (lineid * 4U)) & 0x0fU];
      for (uint8_t from = 1; from <= line_empty_slots[0]; ++from) {
        uint8_t to = line_empty_slots[from];
        // copy value
        bucket->line[lineid].entry[to] = bucket->line[0].entry[from];
        uint8_t to_slotid = lineid * 4 + to;
        if (to_slotid != 15) {
          header->fingerprint[to_slotid - 1] = header->fingerprint[from - 1];
        }
        // validate and invalidate corresponding bitmap
        bitmap = bitmap_set1(bitmap, to_slotid);
        bitmap = bitmap_set0(bitmap, from);
      }
      // after entry moving, persist the new address
      
      
    }
    // atomic write
    header->r_bitmap = bitmap;
    
    
  }

  /**
   * @return -3: Duplicate key exist!; -1: full; 0:
   * Success.
   * @details is_dir_locked must be false
   */
  static int insert(T key, Value_t value, uint64_t hash, bucket_t *first_bucket,
                    uint32_t chunk_id, seg_meta* target_seg, bool &is_dir_locked) {
    bucket_t *bucket0, *bucket1;
    // prefetch 8 cache line of 2 candidate buckets
    uint32_t bucket_id0 = bucket_idx0(hash);
    uint32_t bucket_id1 = bucket_idx1(hash);
    bucket0 = &first_bucket[bucket_id0];
    bucket1 = &first_bucket[bucket_id1];
    LINE_PREF(bucket0);
    LINE_PREF(bucket1);
    uint8_t fgprt = FINGERPRINT(hash);
    uint32_t lock_id0 = lock_hash(chunk_id, bucket_id0);
    uint32_t lock_id1 = lock_hash(chunk_id, bucket_id1);

    uint16_t bitmap;
    bucket_t *bucket;
    int count;
    // LOCK both of them, first lock the smaller one, then fix the bigger one
    // @attention lock depend and dead-lock
    acquire_2locks(lock_id0, lock_id1);

    uint16_t bitmap0 = bucket0->header.r_bitmap;
    uint16_t bitmap1 = bucket1->header.r_bitmap;

    // @warning the function do not consider the reserve bit
    int count0 = _popcnt32(bitmap0), count1 = _popcnt32(bitmap1);

    if (check_key_exist(bucket0, bitmap0, key, fgprt) ||
        check_key_exist(bucket1, bitmap1, key, fgprt)) {
      release2locks(lock_id0, lock_id1);
      release_lock(seg_lock_hash(chunk_id));
      return -3;
    }
    if (count0 < count1) {
      bitmap = bitmap0;
      bucket = bucket0;
      count = count0;
    } else {
      bitmap = bitmap1;
      bucket = bucket1;
      count = count1;
    }
    if(count == 15) {
      target_seg->policy = STASH;
      auto SM = SpaceManager::Get();
      // Warning acquire and release a new lock
      SM->allocStashChunkId(chunk_id);
      release2locks(lock_id0, lock_id1);
      release_lock(seg_lock_hash(chunk_id));
      return -1;
    } else {
      release_lock(seg_lock_hash(chunk_id));
      insert_into_bucket(bucket, key, value, fgprt);
      release2locks(lock_id0, lock_id1);
      return 0;
    }
  }

// opt_insert for 2Choice and 2Choice with Stash
// if change policy, goto pessi_insert
  static int opt_insert(T key, Value_t value, uint64_t hash, bucket_t *first_bucket,
                    uint32_t chunk_id, seg_meta *target_seg, uint32_t old_dir_lock_value, uint32_t old_seg_lock_value) {
    bucket_t *bucket0, *bucket1;
    // prefetch 8 cache line of 2 candidate buckets
    uint32_t bucket_id0 = bucket_idx0(hash);
    uint32_t bucket_id1 = bucket_idx1(hash);
    bucket0 = &first_bucket[bucket_id0];
    bucket1 = &first_bucket[bucket_id1];
    LINE_PREF(bucket0);
    LINE_PREF(bucket1);
    uint8_t fgprt = FINGERPRINT(hash);
    uint32_t lock_id0 = lock_hash(chunk_id, bucket_id0);
    uint32_t lock_id1 = lock_hash(chunk_id, bucket_id1);

    uint16_t bitmap;
    bucket_t *bucket;
    int count;
    // LOCK both of them, first lock the smaller one, then fix the bigger one
    // @attention lock depend and dead-lock
    acquire_2locks(lock_id0, lock_id1);
    if (old_seg_lock_value != get_lock_word(seg_lock_hash(chunk_id))
        || chunk_id != target_seg->chunk_id.load(std::memory_order_acquire)
        || old_dir_lock_value != __atomic_load_n(&dir_lock, __ATOMIC_ACQUIRE)
        ) {
      release2locks(lock_id0, lock_id1);
      return -4;
    }

    uint16_t bitmap0 = bucket0->header.r_bitmap;
    uint16_t bitmap1 = bucket1->header.r_bitmap;

    // @warning the function do not consider the reserve bit
    int count0 = _popcnt32(bitmap0), count1 = _popcnt32(bitmap1);

    if (check_key_exist(bucket0, bitmap0, key, fgprt) ||
        check_key_exist(bucket1, bitmap1, key, fgprt)) {
      release2locks(lock_id0, lock_id1);
      return -3;
    }
    if (count0 < count1) {
      bitmap = bitmap0;
      bucket = bucket0;
      count = count0;
    } else {
      bitmap = bitmap1;
      bucket = bucket1;
      count = count1;
    }
    if(count == 15) {
      release2locks(lock_id0, lock_id1);
      return -1;
    } else {
      insert_into_bucket(bucket, key, value, fgprt);
      release2locks(lock_id0, lock_id1);
      return 0;
    }
  }

  /*if success, then return 0, else return -1*/
  static int del(T key, uint64_t hash, bucket_t *first_bucket,
                 uint32_t chunk_id) {
    // Calculate 2 indexes of buckets
    uint32_t bucket_id0 = bucket_idx0(hash);
    uint32_t bucket_id1 = bucket_idx1(hash);
    bucket_t *bucket0 = &first_bucket[bucket_id0];
    bucket_t *bucket1 = &first_bucket[bucket_id1];
    LINE_PREF(bucket0);
    LINE_PREF(bucket1);
    uint8_t fgprt = FINGERPRINT(hash);

    uint32_t lock_id0 = lock_hash(chunk_id, bucket_id0);
    uint32_t lock_id1 = lock_hash(chunk_id, bucket_id1);

    int ret_val =
        bucket_delete(key, lock_id1, bucket1, fgprt, chunk_id);
    if (ret_val == -1)
      ret_val =
          bucket_delete(key, lock_id0, bucket0, fgprt, chunk_id);
    return ret_val;
  }
};

template <class T>
class doubleWithStash : public hashPolicy {
 public:
  typedef Bucket<T> bucket_t;
  typedef _Pair<T> pea_entry;
  doubleWithStash() = default;
  ~doubleWithStash() = default;

  static Value_t inline bucket_lookup(T key, uint32_t lock_id, bucket_t *bucket,
                                      uint8_t fgprt) {
    uint8_t j;
    pea_entry *this_entry;
    Value_t ret_val;
  RE_SEARCH:
    uint32_t lock_word = get_lock_word(lock_id);
    if (is_locked(&lock_word)) goto RE_SEARCH;
    header_t *hdr = &bucket->header;
    uint16_t bitmap = hdr->r_bitmap;

    // for fingerprint 1-14 SIMD comparison
    // a. set every byte to key_hash in a 16B register
    __m128i key_16B = _mm_set1_epi8((char)fgprt);

    // b. load meta into another 16B register
    __m128i fgprt_16B = _mm_load_si128((const __m128i *)hdr);

    // c. compare them
    __m128i cmp_res = _mm_cmpeq_epi8(key_16B, fgprt_16B);

    // d. generate a mask
    auto mask = (unsigned int)_mm_movemask_epi8(cmp_res);  // 1: same; 0: diff

    // AND bitmap
    mask = mask & ((unsigned int)bitmap << 1U);

    // search every matching candidate
    while (mask) {
      int jj = bitScan(mask) - 1;  // next candidate
      this_entry = bucket->entry + jj - 1;
      if (key == this_entry->key) {
        ret_val = this_entry->value;
        if (get_lock_word(lock_id) != lock_word) goto RE_SEARCH;
        return ret_val;
      }

      mask &= ~(0x1U << (uint16_t)jj);  // remove this bit
    }                                   // end while
    if (get_bit_true(bitmap, 15)) {
      this_entry = bucket->entry + 15;
      if (key == this_entry->key) {
        ret_val = this_entry->value;
        if (get_lock_word(lock_id) != lock_word) goto RE_SEARCH;
        return ret_val;
      }
    }
    return NONE;
  }

  static int inline bucket_delete(T key, bucket_t *bucket,
                                  uint8_t fgprt, uint32_t chunk_id) {
    pea_entry *this_entry;
    uint8_t j;
    header_t *hdr = &bucket->header;
    uint16_t bitmap = hdr->r_bitmap;
    __m128i key_16B = _mm_set1_epi8((char)fgprt);
    __m128i fgprt_16B = _mm_load_si128((const __m128i *)hdr);
    __m128i cmp_res = _mm_cmpeq_epi8(key_16B, fgprt_16B);
    auto mask = (unsigned int)_mm_movemask_epi8(cmp_res);  // 1: same; 0: diff
    mask = mask & ((unsigned int)bitmap << 1U);
    while (mask) {
      int jj = bitScan(mask) - 1;  // next candidate
      this_entry = bucket->entry + jj - 1;
      if (key == this_entry->key) {
        hdr->r_bitmap = bitmap_set0(bitmap, jj - 1);       
        return 0;
      }
      mask &= ~(0x1U << (uint16_t)jj);  // remove this bit
    }

    if (get_bit_true(bitmap, 15)) {
      this_entry = bucket->entry + 15;
      if (key == this_entry->key) {
        hdr->r_bitmap = bitmap_set0(bitmap, 15);
        
        
        return 0;
      }
    }
    return -1;
  }

  static Value_t lookup(T key, uint64_t hash, bucket_t *first_bucket,
                        uint32_t chunk_id) {
    // Calculate 2 indexes of buckets
    uint32_t bucket_id0 = bucket_idx0(hash);
    uint32_t bucket_id1 = bucket_idx1(hash);
    bucket_t *bucket0 = &first_bucket[bucket_id0];
    bucket_t *bucket1 = &first_bucket[bucket_id1];
    LINE_PREF(bucket0);
    LINE_PREF(bucket1);
    auto SM = SpaceManager::Get();
    uint32_t stash_c_id = SM->get_stash_cid(chunk_id);
    auto stash_chunk =
        reinterpret_cast<Bucket<T> *>(SM->chunk_addr(stash_c_id));
    uint32_t stash_b0id = (chunk_id << 1U) % kNumBucket;
    uint32_t stash_b1id = stash_b0id + 1;
    bucket_t *bucket_s0 = stash_chunk + stash_b0id;
    bucket_t *bucket_s1 = stash_chunk + stash_b1id;
    LINE_PREF(bucket_s0);
    LINE_PREF(bucket_s1);

    uint32_t lock_id0 = lock_hash(chunk_id, bucket_id0);
    uint32_t lock_id1 = lock_hash(chunk_id, bucket_id1);

    uint32_t chunk_lock = chunk_id;
    uint8_t fgprt = FINGERPRINT(hash);
    Value_t ret_val;

    ret_val = bucket_lookup(key, lock_id1, bucket1, fgprt);
    if (ret_val == NONE) ret_val = bucket_lookup(key, lock_id0, bucket0, fgprt);
    if (ret_val == NONE)
      ret_val = bucket_lookup(key, chunk_lock, bucket_s0, fgprt);
    if (ret_val == NONE)
      ret_val = bucket_lookup(key, chunk_lock, bucket_s1, fgprt);
    return ret_val;
  }

  static inline bool check_key_exist(bucket_t *bucket, uint16_t bitmap, T key,
                                     uint8_t fgprt) {
    pea_entry *this_entry;
    header_t *hdr = &bucket->header;
    __m128i key_16B = _mm_set1_epi8((char)fgprt);
    __m128i fgprt_16B = _mm_load_si128((const __m128i *)hdr);
    __m128i cmp_res = _mm_cmpeq_epi8(key_16B, fgprt_16B);
    auto mask = (unsigned int)_mm_movemask_epi8(cmp_res);  // 1: same; 0: diff
    mask = mask & ((unsigned int)bitmap << 1U);
    while (mask) {
      int jj = bitScan(mask) - 1;  // next candidate
      this_entry = bucket->entry + jj - 1;
      if (key == this_entry->key) {
        return true;
      }
      mask &= ~(0x1U << (uint16_t)jj);  // remove this bit
    }                                   // end while

    if (get_bit_true(bitmap, 15)) {
      this_entry = bucket->entry + 15;
      if (key == this_entry->key) {
        return true;
      }
    }
    return false;
  }

  static inline void insert_into_bucket(bucket_t *bucket, T key, Value_t value, uint8_t fgprt) {
    header_t * header = &bucket->header;
    uint16_t bitmap = header->r_bitmap;
    uint8_t slot = first_empty_slot(bitmap);
    // write entry
    bucket->entry[slot] = {.key = key, .value = value};

    uint8_t lineid = slot >> 2U;
    if (lineid == 0) {
      header->fingerprint[slot - 1] = fgprt;
      // write word 0
      // flush
      bitmap = bitmap_set1(bitmap, slot);
    } else {
      if (slot != 15) {
        header->fingerprint[slot - 1] = fgprt;
      }
      bitmap = bitmap_set1(bitmap, slot);
      uint8_t *line_empty_slots =
          empty_slots_in_line[(uint8_t)(bitmap >> (lineid * 4U)) & 0x0fU];
      for (uint8_t from = 1; from <= line_empty_slots[0]; ++from) {
        uint8_t to = line_empty_slots[from];
        // copy value
        bucket->line[lineid].entry[to] = bucket->line[0].entry[from];
        uint8_t to_slotid = lineid * 4 + to;
        if (to_slotid != 15) {
          header->fingerprint[to_slotid - 1] = header->fingerprint[from - 1];
        }
        // validate and invalidate corresponding bitmap
        bitmap = bitmap_set1(bitmap, to_slotid);
        bitmap = bitmap_set0(bitmap, from);
      }
      // after entry moving, persist the new address
      
      
    }
    // atomic write
    header->r_bitmap = bitmap;
    
    
  }

  /**
   * @param hash the hash value of key
   * @param first_bucket the addr of the first bucket of the chunk
   * @param chunk_id the chunk_id in the SM
   * @return -3: Duplicate key exist!; -1: full; 0:
   * Success.
   */
   static int insert(T key, Value_t value, uint64_t hash, bucket_t *first_bucket,
                    uint32_t chunk_id, seg_meta* target_seg, bool &is_dir_locked) {
    bucket_t *bucket0, *bucket1;
    uint32_t bucket_id0 = bucket_idx0(hash);
    uint32_t bucket_id1 = bucket_idx1(hash);
    bucket0 = &first_bucket[bucket_id0];
    bucket1 = &first_bucket[bucket_id1];
    LINE_PREF(bucket0);
    LINE_PREF(bucket1);
    // find the stash bucket from chunk_id
    // line prefetch(stash bucket)
    uint8_t fgprt = FINGERPRINT(hash);

    uint32_t lock_id0 = lock_hash(chunk_id, bucket_id0);
    uint32_t lock_id1 = lock_hash(chunk_id, bucket_id1);

    uint16_t bitmap;
    bucket_t *bucket;
    int count;
    // LOCK both of them, first lock the smaller one, then fix the bigger one
    // @attention lock depend and dead-lock
    acquire_2locks(lock_id0, lock_id1);

    uint16_t bitmap0 = bucket0->header.r_bitmap;
    uint16_t bitmap1 = bucket1->header.r_bitmap;

    // @warning the function do not consider the reserve bit
    int count0 = _popcnt32(bitmap0), count1 = _popcnt32(bitmap1);
    // release another bucket's lock
    if (check_key_exist(bucket0, bitmap0, key, fgprt) || check_key_exist(bucket1, bitmap1, key, fgprt) ) {
      release2locks(lock_id0, lock_id1);
      release_lock(seg_lock_hash(chunk_id));
      if (is_dir_locked) release_dir_lock();
      return -3;
    }
    if (count0 < count1) {
      bitmap = bitmap0;
      bucket = bucket0;
      count = count0;
    } else {
      bitmap = bitmap1;
      bucket = bucket1;
      count = count1;
    }
    if (count < 15) {
      release_lock(seg_lock_hash(chunk_id));
      if (is_dir_locked) release_dir_lock();
      insert_into_bucket(bucket, key, value, fgprt);
      release2locks(lock_id0, lock_id1);
      return 0;
    }
    release2locks(lock_id0, lock_id1);

    // STASH

    auto SM = SpaceManager::Get();
    uint32_t stash_c_id = SM->get_stash_cid(chunk_id);
    auto stash_chunk =
        reinterpret_cast<Bucket<T> *>(SM->chunk_addr(stash_c_id));
    uint32_t stash_b0id = (chunk_id << 1U) % kNumBucket;
    bucket_t *bucket_s0 = stash_chunk + stash_b0id;
    uint32_t stash_b1id = stash_b0id + 1;
    bucket_t *bucket_s1 = stash_chunk + stash_b1id;
    LINE_PREF(bucket_s0);
    LINE_PREF(bucket_s1);

    bitmap0 = bucket_s0->header.r_bitmap;
    bitmap1 = bucket_s1->header.r_bitmap;
    // WARNING unique check is default and melted
    if (check_key_exist(bucket_s1, bitmap1, key, fgprt) ||
        check_key_exist(bucket_s0, bitmap0, key, fgprt)) {
      release_lock(seg_lock_hash(chunk_id));
      if (is_dir_locked) release_dir_lock();
      return -3;
    }
    if (_popcnt32(bitmap0) < 15){
      bucket = bucket_s0;
    } else if (_popcnt32(bitmap1) < 15) {
      bucket = bucket_s1;
    } else {
      // both buckets are full, is_seg_locked = true
      // get every bucket lock to ensure all bucket insert in the segment are finished
      uint32_t first_lock_id = chunk_id * LOCK_SIZE_PER_SEG + 1;
      for (size_t i = 0; i < LOCK_SIZE_PER_SEG - 1; i++) {
        acquire_lock(first_lock_id + i);
      }
      return -2;
    }
    // stash buckets are not all full
    if (is_dir_locked) {
      release_dir_lock();
    }
    insert_into_bucket(bucket, key, value, fgprt);
    release_lock(seg_lock_hash(chunk_id));
    return 0;
  }
  
  /*if success, then return 0; -1: not find*/
  static int del(T key, uint64_t hash, bucket_t *first_bucket, uint32_t chunk_id) {
    // Calculate 2 indexes of buckets
    uint32_t bucket_id0 = bucket_idx0(hash);
    uint32_t bucket_id1 = bucket_idx1(hash);
    bucket_t *bucket0 = &first_bucket[bucket_id0];
    bucket_t *bucket1 = &first_bucket[bucket_id1];
    LINE_PREF(bucket0);
    LINE_PREF(bucket1);
    auto SM = SpaceManager::Get();
    uint32_t stash_c_id = SM->get_stash_cid(chunk_id);
    auto stash_chunk =
        reinterpret_cast<Bucket<T> *>(SM->chunk_addr(stash_c_id));
    uint32_t stash_b0id = (chunk_id << 1U) % kNumBucket;
    uint32_t stash_b1id = stash_b0id + 1;
    bucket_t *bucket_s0 = stash_chunk + stash_b0id;
    bucket_t *bucket_s1 = stash_chunk + stash_b1id;

    uint8_t fgprt = FINGERPRINT(hash);

    uint32_t lock_id0 = lock_hash(chunk_id, bucket_id0);
    uint32_t lock_id1 = lock_hash(chunk_id, bucket_id1);

    uint32_t chunk_lock = chunk_id;

    acquire_lock(lock_id1);
    int ret_val = bucket_delete(key, bucket1, fgprt, chunk_id);
    release_lock(lock_id1);
    if (ret_val == -1) {
      acquire_lock(lock_id0);  
      ret_val = bucket_delete(key, bucket0, fgprt, chunk_id);
      release_lock(lock_id0);
    }
    if (ret_val == -1) {
      LINE_PREF(bucket_s0);
      LINE_PREF(bucket_s1);
      acquire_lock(chunk_lock);
      ret_val = bucket_delete(key, bucket_s0, fgprt, chunk_id);
      if (ret_val == -1) {
        ret_val = bucket_delete(key, bucket_s1, fgprt, chunk_id);
      }
      release_lock(chunk_lock);
    }
    return ret_val;
  }
};

template <class T>
void PeaHashing<T>::init_function_list() {
  class_list[SINGLE_HASH].insert = singleHash<T>::insert;
  class_list[SINGLE_HASH].del = singleHash<T>::del;
  class_list[SINGLE_HASH].lookup = singleHash<T>::lookup;
  class_list[SINGLE_HASH].opt_insert = singleHash<T>::opt_insert;
  class_list[DOUBLE_HASH].insert = doubleHash<T>::insert;
  class_list[DOUBLE_HASH].del = doubleHash<T>::del;
  class_list[DOUBLE_HASH].lookup = doubleHash<T>::lookup;
  class_list[DOUBLE_HASH].opt_insert = doubleHash<T>::opt_insert;
  class_list[STASH].insert = doubleWithStash<T>::insert;
  class_list[STASH].del = doubleWithStash<T>::del;
  class_list[STASH].lookup = doubleWithStash<T>::lookup;
  class_list[STASH].opt_insert = doubleHash<T>::opt_insert; // Attention: if double hash fail goto double with stash
}

template <class T>
void PeaHashing<T>::Recovery() {
}

void check_all_lock_released() {
  // for (unsigned int i : lock_map_table) {
  //   if (i & lockSet) LOG("The concurrency control is problem!");
  // }
}

template <class T>
uint64_t PeaHashing<T>::getNumber() {
  check_all_lock_released();
  seg_meta *dir_entry = dir->_;
  uint8_t g = *global_depth;
  uint8_t depth_diff;
  uint64_t count = 0;
  uint64_t seg_count = 0;
  seg_meta *s;
  int capacity = pow(2, g);
  uint32_t SCN = SM->stash_chunk_num;
  for (int i = 0; i < capacity;) {
    s = dir_entry + i;
    depth_diff = g - s->local_depth;
    auto first_bucket =
        reinterpret_cast<Bucket<T> *>(SM->chunk_addr(s->chunk_id));
    int stash = 0;
    for (int j = i; j < i + pow(2, depth_diff); ++j) {
      if (dir_entry[j].policy == STASH) {
        stash = 1;
      }
    }
    uint32_t s_c_id_in_array = s->chunk_id / stashUnitsPChunk;
    if (stash == 1) {
      uint32_t stash_c_id = __atomic_load_n(SM->stash_chunk_id_array + s_c_id_in_array, __ATOMIC_ACQUIRE);
      auto stash_chunk =
          reinterpret_cast<Bucket<T> *>(SM->chunk_addr(stash_c_id));
      uint32_t stash_b_id = (s->chunk_id << 1U) % kNumBucket;
      Bucket<T> *bucket_s = stash_chunk + stash_b_id;
      for (int j = 0; j < 2; ++j) {
        bucket_s += j;
        header_t *check_header = &bucket_s->header;
        uint16_t bitmap = check_header->r_bitmap;
        for (uint8_t slot = 1; slot < 16; ++slot) {
          if (hashPolicy::get_bit_true(bitmap, slot)) {
            count++;
          }
        }
      }
    }
    for (uint64_t j = 0; j < kNumBucket; ++j) {
      Bucket<T> *check_bucket = first_bucket + j;
      header_t *check_header = &check_bucket->header;
      uint16_t bitmap = check_header->r_bitmap;
      for (uint8_t slot = 1; slot < 16; ++slot) {
        if (hashPolicy::get_bit_true(bitmap, slot)) {
          count++;
        }
      }
    }
    seg_count++;
    i += pow(2, depth_diff);
  }
  std::cout << "seg_count = " << seg_count << std::endl;
  std::cout << "entry_count = " << count << std::endl;

#ifdef AU
  double average_utility = (entry_num + 1.0) * entry_num / (2.0 * space_num);
  std::cout << "segment num = " << segment_num << std::endl;
  std::cout << "entry num = " << entry_num << std::endl;
  std::cout << "space unit = " << space_num << std::endl;
  std::cout << "AU = " << average_utility << std::endl;
#endif
  return count;
}

template <class T>
PeaHashing<T>::PeaHashing(void) {
  std::cout << "Reinitialize up" << std::endl;
}

template <class T>
int PeaHashing<T>::Insert(T key, Value_t value, bool is_load) {
  if (is_load) {
    return Insert(key, value);
  } else {
    LOG_FATAL("Please use another interface");
  }
}

template <class T>
bool PeaHashing<T>::Delete(T key, bool flag) {
  return Delete(key);
}

template <class T>
Value_t PeaHashing<T>::Get(T key, bool is_in_epoch) {
  return Get(key);
}

}  // namespace pea
