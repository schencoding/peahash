// Copyright (c) 2021-2023 Institute of Computing Technology, Chinese Academy of Sciences
// Pea Hash is licensed under Mulan PSL v2.
/**
 *
 * This file allocate space in persistent memory to init the hash table.
 *
 */
#pragma once

#include <pair.h>
#include "x86intrin.h"
#include <cstdlib>

std::atomic<uint8_t> gl = -1;
thread_local uint8_t thread_id = (gl++);      //TODO add some thread initialize to initialzie the thread_mask

constexpr uint32_t bucketNumPChunk = 2048;
constexpr size_t stashUnitsPChunk = bucketNumPChunk / 2;
constexpr uint8_t degreeMax = 11;
constexpr uint32_t doubleBPChunk = bucketNumPChunk / 2;
constexpr uint64_t bucket_size = 256UL;
constexpr uint64_t chunk_size = bucketNumPChunk * 256UL;
constexpr uint16_t dir_alignment = kCacheLineSize;
constexpr uint32_t stash_lock_set = 1U << 31U;
constexpr uint32_t stash_version_mask = (1U << 31U) - 1;
constexpr uint32_t bucket_LLI_size = 256 / sizeof(long long int);
constexpr uint32_t deg_mask = 0xfU;
constexpr uint32_t number_mask = ((1U << 16) - 1) << 4;
constexpr uint64_t number_64_mask = ~(((1UL << 16) - 1) << 4);
constexpr uint32_t block_id_mask = ((1U << 11) - 1) << 21;

struct chunkMeta {
  uint32_t chunk_id;
  uint32_t FEBid;  // first empty bucket id
  uint64_t free_list;
  uint32_t block_bitmap[0];
};

struct bucket_log {
  uint8_t mem[256];
};

struct ptr_log {
  uint64_t head_ptr;
  uint64_t tail_ptr;
};

struct log {
  uint32_t buc_log_cur;
  uint32_t ptr_log_cur;
  bucket_log buc_entry;
  ptr_log ptr_entry[0];
};
typedef struct log log_t;

struct SpaceManager {
 private:
  int *is_pmem;
  size_t *mapped_len;

  // DRAM variable
  uint32_t bitmap_len;
  static SpaceManager *instance_;
  chunkMeta *valueMeta[degreeMax + 1];

 public:
  // set denotes empty, unset denotes used
  uint32_t *chunk_bitmap;
  uint32_t FECid;            // id of first empty chunk
  uint32_t stash_chunk_num; 
  // 31: lock bit; 30:0 version
  uint32_t rw_lock;

  // NVM address
  void *pmem_addr;
  uint32_t *stash_chunk_id_array;
  log_t *log;
  char *meta0;
  char *dir0;
  char *meta1;
  char *dir1;

  static void Initialize(const char *pool_name, size_t pool_size,
                         uint64_t thread_num) {
    instance_ = new SpaceManager(pool_name, pool_size, thread_num);
    // std::cout << "Pool opened." << std::endl;
  }

  // Warning: use it
  void Close_pool() {
    free(pmem_addr);
    delete mapped_len;
    delete is_pmem;
    free(chunk_bitmap);
    delete instance_;
  }

  static inline char *cache_aligned(char *meta) {
    uint64_t remainder = (uint64_t)meta % kCacheLineSize;
    if (remainder != 0)
      return meta + kCacheLineSize - remainder;
    else
      return meta;
  }

  static inline char *bucket_aligned(char *meta) {
    uint64_t remainder = (uint64_t)meta % 256;
    if (remainder != 0)
      return meta + 256 - remainder;
    else
      return meta;
  }

  SpaceManager(const char *pool_name, size_t pool_size, uint64_t thread_num) {
    // 512KB allignment
    int ret = posix_memalign(&pmem_addr, 524288, pool_size);
    if (ret) {
      fprintf (stderr, "posix_memalign: %s\n", strerror (ret));
      exit(-1);
    }
    if (!pmem_addr) {
      fprintf(stderr, "pmem_map_file fails: 1\n");
      exit(1);
    }

    // initialize DRAM fields
    size_t chunk_cap = pool_size / chunk_size;
    bitmap_len = chunk_cap / 32UL;
    uint32_t stash_chunk_cap = chunk_cap / bucketNumPChunk + 1;
    chunk_bitmap = new uint32_t[bitmap_len];
    size_t dir_chunk_num = 1 + chunk_cap * 20 / chunk_size;
    chunk_bitmap[0] = 0xfffffffcU;
    if (bitmap_len > 1) memset(chunk_bitmap + 1, 0xff, bitmap_len * 4 - 4);
    FECid = 2;
    for (int i = 2; i < 2 * dir_chunk_num; ++i) {
      allocChunk();
    }
    stash_chunk_num = 0;
    rw_lock = 0;

    // initialize NVM fields
    stash_chunk_id_array = static_cast<uint32_t *>(pmem_addr);
    char *tail = (char *)stash_chunk_id_array + 4 * stash_chunk_cap;
    for (int i = 0; i <= degreeMax; ++i) {
      valueMeta[i] = reinterpret_cast<chunkMeta *>(cache_aligned(tail));
      valueMeta[i]->chunk_id = allocChunk();
      valueMeta[i]->FEBid = 0;
      valueMeta[i]->free_list = 0;
      uint16_t block_bitmap_len = (i > 5) ? 1 : 1 << (6 - i);
      memset(valueMeta[i]->block_bitmap, 0xff, block_bitmap_len * 4);
      tail = (char *)valueMeta[i] + 16 + sizeof(uint32_t) * block_bitmap_len;
    }
    // initialize LOG
    log = reinterpret_cast<log_t *>(cache_aligned(tail));
    log->buc_log_cur = 0;
    log->ptr_log_cur = 0;
    tail = (char *)log + 1400;
    meta0 = bucket_aligned(tail);
    // make sure start of chunk is 256B alignment
    dir0 = meta0 + dir_alignment;
    meta1 = meta0 + chunk_size * dir_chunk_num;
    dir1 = meta1 + dir_alignment;
  }

  void bitmap_set(uint32_t id) {
    uint32_t *bitmap = chunk_bitmap + id / 32;
    uint8_t index_32 = id % 32;
    uint32_t old_value = *bitmap;
    *bitmap = old_value & ~(1U << index_32);
  }

  static SpaceManager *Get() { return instance_; }

  bool try_get_rwlock() {
    uint32_t old_value = __atomic_load_n(&rw_lock, __ATOMIC_ACQUIRE);
    if (old_value & stash_lock_set) {
      return false;
    }
    uint32_t new_value = old_value | stash_lock_set;
    return CAS(&rw_lock, &old_value, new_value);
  }

  void acquire_rwlock() {
    while (try_get_rwlock() == false);    
  }

  void release_rwlock() {
    uint32_t v = (rw_lock + 1) & stash_version_mask;
    __atomic_store_n(&rw_lock, v, __ATOMIC_RELEASE);
  }

  uint32_t get_rwlock() {
    return __atomic_load_n(&rw_lock, __ATOMIC_ACQUIRE);
  }

  uint32_t get_stash_cid(uint32_t chunk_id) {
    uint32_t retval;
  RE_GET:
    uint32_t lock = get_rwlock();
    if (lock & stash_lock_set) {
      goto RE_GET;
    } else {
      retval = stash_chunk_id_array[chunk_id / stashUnitsPChunk];
    }
    if (lock != get_rwlock()){
      goto RE_GET;
    } else {
      return retval;
    }
  }

  /**
   * @param chunk_offset
   * @return allocate 1MB pmem
   */
  inline void *chunk_addr(uint32_t chunk_id) {
    return meta0 + (uint64_t)chunk_id * chunk_size;
  }

  //  inline void *bucket_addr(uint32_t buc_id) {
  //    uint32_t cid = chunk_table[thread_id].chunk_id;
  //    return (char *)chunk_addr(cid) + buc_id * bucket_size;
  //  }

  inline uint32_t chunk_id(char *chunk_addr) {
    return (chunk_addr - meta0) / chunk_size;
  }

  uint32_t allocChunk() {
    SCAN:
    uint32_t old_FECid = __atomic_load_n(&FECid, __ATOMIC_ACQUIRE);
    uint32_t bitword_id = old_FECid / 32U;
    uint32_t new_value;
    uint32_t *bitmap_start = chunk_bitmap + bitword_id;
    uint32_t old_value;
    int index_32;
    uint32_t *scan_bitmap = bitmap_start - 1;
    do {
      scan_bitmap++;
      old_value = __atomic_load_n(scan_bitmap, __ATOMIC_ACQUIRE);
    } while (old_value == 0);
    index_32 = _bit_scan_forward(old_value);
    // unset
    new_value = old_value & ~(1U << index_32);
    if (!CAS(scan_bitmap, &old_value, new_value)) {
      goto SCAN;
    }

    uint32_t new_FECid = (scan_bitmap - chunk_bitmap) * 32 + index_32;
    __atomic_store_n(&FECid, new_FECid + 1, __ATOMIC_RELEASE);
    return new_FECid;
  }

  void alloc4chunk(uint32_t *chunk_ids) {
    uint8_t cnt = 0;
    do {
      SCAN:
      uint32_t old_FECid = __atomic_load_n(&FECid, __ATOMIC_ACQUIRE);
      uint32_t bitword_id = old_FECid / 32U;
      uint32_t new_value, old_value;
      uint32_t *bitmap_start = chunk_bitmap + bitword_id;
      uint32_t index_32;
      uint32_t *scan_bitmap = bitmap_start - 1;
      do {
        scan_bitmap++;
        old_value = __atomic_load_n(scan_bitmap, __ATOMIC_ACQUIRE);
      } while (old_value == 0);

      uint8_t available_cnt = _popcnt32(old_value);
      uint8_t to_alloc_cnt =
          (available_cnt < 4 - cnt) ? available_cnt : 4 - cnt;
      new_value = old_value;
      for (int i = 0; i < to_alloc_cnt; ++i) {
        index_32 = _bit_scan_forward(new_value);
        chunk_ids[cnt + i] = index_32;
        new_value = new_value & ~(1U << index_32);
      }
      if (!CAS(scan_bitmap, &old_value, new_value)) {
        goto SCAN;
      }
      for (int i = 0; i < to_alloc_cnt; ++i) {
        chunk_ids[cnt + i] += (scan_bitmap - chunk_bitmap) * 32;
      }
      cnt += to_alloc_cnt;
      __atomic_store_n(&FECid, chunk_ids[cnt - 1] + 1, __ATOMIC_RELEASE);
    } while (cnt != 4);
  }

  // w clear
  void allocStashChunkId(uint32_t chunk_id) {
    acquire_rwlock();
    uint32_t id_in_stash_chunk_array = chunk_id / doubleBPChunk;
    // CAS lock until success, load old_stash_chunk_num
    if (stash_chunk_num <= id_in_stash_chunk_array) {
      for (int i = stash_chunk_num; i <= id_in_stash_chunk_array; ++i) {
        uint32_t *stash_chunk_id = stash_chunk_id_array + i;
        uint32_t new_chunk = allocChunk();
        auto stash_seg = (long long int *)chunk_addr(new_chunk);
        for (uint32_t j = 0; j < bucketNumPChunk; ++j) {
          _mm_stream_si64(stash_seg + j * bucket_LLI_size, 0);
        }
        // pmem_drain();
        *stash_chunk_id = new_chunk;
        // pmem_flush(stash_chunk_id, 64);
      }
      // atomic store array_len ++ and release lock
      stash_chunk_num = id_in_stash_chunk_array + 1;
    }
    release_rwlock();
  }

  void recoveryStashChunk(uint32_t chunk_id) {
    uint32_t id_in_stash_chunk_array = chunk_id / doubleBPChunk;
    for (int i = stash_chunk_num; i <= id_in_stash_chunk_array; ++i) {
      bitmap_set(stash_chunk_id_array[i]);
    }
  }

  void freeChunk(uint32_t free_id) {
    uint32_t *bitmap = chunk_bitmap + free_id / 32;
    uint8_t index_32 = free_id % 32;
    uint32_t old_value = __atomic_load_n(bitmap, __ATOMIC_ACQUIRE);
    uint32_t new_value;
    do {
      // set
      new_value = old_value | (1U << index_32);
    } while (!CAS(bitmap, &old_value, new_value));

    uint32_t old_FECid = __atomic_load_n(&FECid, __ATOMIC_ACQUIRE);
    do {
      if (free_id >= old_FECid) break;
    } while (!CAS(&FECid, &old_FECid, free_id));
  }

  void freeChunk(char *chunk) {
    uint32_t free_id = chunk_id(chunk);
    freeChunk(free_id);
  }

  char *getAddr(uint32_t chunk_id, uint32_t block_id, uint32_t degree) {
    return (char *)chunk_addr(chunk_id) + block_id * (bucket_size << degree);
  }

  static inline uint64_t makePtr(uint64_t chunk_id, u_int64_t block_id,
                                 uint64_t degree) {
    return (chunk_id << 32) | (block_id << 21) | degree;
  }

  uint64_t *getHdr(uint64_t val_ptr) {
    uint64_t chunk_id = val_ptr >> 32;
    uint64_t block_id = (val_ptr & (uint64_t)block_id_mask) >> 21;
    uint64_t degree = val_ptr & 0xfUL;
    char *block_hdr =
        (char *)chunk_addr(chunk_id) + block_id * (bucket_size << degree);
    return (uint64_t *)block_hdr;
  }

  /**
   * @brief  deg = 11, the free list uses 64bit that is chunk id + 4 bit deg,
   the pointer is in the former 64 bit deg < 11, the free list uses 64 bit
   value_ptr, the pointer in the former 64 bit of a block is also a value_ptr
   * @return void *: the block address; uint64_t: the pointer
   * @attention deg = 11: W/O clean the next hdr
   */
  std::pair<void *, uint64_t> allocBlock(uint8_t degree) {
    chunkMeta *local_chunk = valueMeta[degree];
    if (local_chunk->free_list >> 32U) {
      uint64_t dangling_blk = local_chunk->free_list;
      // log alloc header and itself(tail)
      log_ptr_write(dangling_blk, dangling_blk);
      uint64_t *new_hdr;
      new_hdr = getHdr(dangling_blk);
      // update free hdr
      local_chunk->free_list = *new_hdr;
      // pmem_persist(&local_chunk->free_list, 8);
      // give prior hdr
      return std::make_pair(new_hdr, dangling_blk);
    }
    uint32_t *FEBid = &local_chunk->FEBid;
    RE_ALLOC:
    uint32_t old_FEBid = __atomic_load_n(FEBid, __ATOMIC_ACQUIRE);
    // when bitmap is full
    if (old_FEBid == bucketNumPChunk >> degree) {
      local_chunk->chunk_id = allocChunk();
      local_chunk->FEBid = 0;
      uint16_t block_bitmap_len = (degree > 5) ? 1 : 1 << (6 - degree);
      memset(local_chunk->block_bitmap, 0xff, block_bitmap_len * 4);
      goto RE_ALLOC;
    }
    uint32_t bitword_id = old_FEBid / 32U;
    uint32_t new_value;
    uint32_t *bitmap_start = local_chunk->block_bitmap + bitword_id;
    uint32_t old_value;
    int index_32;
    uint32_t *scan_bitmap = bitmap_start - 1;
    do {
      scan_bitmap++;
      old_value = __atomic_load_n(scan_bitmap, __ATOMIC_ACQUIRE);
    } while (old_value == 0);
    index_32 = _bit_scan_forward(old_value);
    new_value = old_value & ~(1U << index_32);
    if (!CAS(scan_bitmap, &old_value, new_value)) {
      goto RE_ALLOC;
    }
    uint32_t new_FEBid =
        (scan_bitmap - local_chunk->block_bitmap) * 32 + index_32;
    __atomic_store_n(FEBid, new_FEBid + 1, __ATOMIC_RELEASE);
    // pmem_persist(FEBid, 4);
    uint64_t cid = local_chunk->chunk_id;
    void *ret1 = (char *)chunk_addr(cid) + new_FEBid * (bucket_size << degree);
    return std::make_pair(ret1, makePtr(cid, new_FEBid, degree));
  }

  void updateFreeHeader(uint64_t ptr, uint8_t degree) {
    valueMeta[degree]->free_list = ptr;
    // pmem_flush(&valueMeta[degree]->free_list, 8);
    // pmem_drain();
  }

  uint64_t *getFreeList(uint8_t degree) {
    return &valueMeta[degree]->free_list;
  }

  // todo, in real recovery, should know the bucket addr
  void log_buc_write(void *bucket) {
    memcpy(log->buc_entry.mem, bucket, 256);
    // pmem_flush(log->buc_entry.mem, 256);
    // pmem_drain();
    log->buc_log_cur = 1;
    // pmem_flush(&log->buc_log_cur, sizeof(uint32_t));
    // pmem_drain();
  }

  void log_buc_clean() {
    log->buc_log_cur = 0;
    // pmem_flush(&log->buc_log_cur, sizeof(uint32_t));
    // pmem_drain();
  }

  // todo, in real recovery, should know the bucket is deleting or allocating,
  // denote degree
  void log_ptr_write(uint64_t head, uint64_t tail) {
    log->ptr_entry[log->ptr_log_cur] = {head, tail};
    // pmem_flush(log->ptr_entry + log->ptr_log_cur, 16);
    // pmem_drain();
    log->ptr_log_cur++;
    // pmem_flush(&log->ptr_log_cur, 4);
    // pmem_drain();
  }

  void log_ptr_clean() {
    log->ptr_log_cur--;
    // pmem_flush(&log->ptr_log_cur, 4);
    // pmem_drain();
  }

  void log_ptr_clear() {
    log->ptr_log_cur = 0;
    // pmem_flush(&log->ptr_log_cur, 4);
    // pmem_drain();
  }
};

SpaceManager *SpaceManager::instance_ = nullptr;

static void ZAllocate(void** ptr, uint32_t alignment, size_t size) {
  int ret = posix_memalign(ptr, alignment, size);
  if (ret) {
    fprintf (stderr, "posix_memalign: %s\n",
           strerror (ret));
  }
  memset(*ptr, 0, size);
}
