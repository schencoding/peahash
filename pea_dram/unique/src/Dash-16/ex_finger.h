
// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong. All rights reserved.
// Licensed under the MIT license.
//
// Dash Extendible Hashing
// Authors:
// Baotong Lu <btlu@cse.cuhk.edu.hk>
// Xiangpeng Hao <xiangpeng_hao@sfu.ca>
// Tianzheng Wang <tzwang@sfu.ca>

#pragma once

#include <immintrin.h>
#include <omp.h>

#include <bitset>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <shared_mutex>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "../../util/hash.h"
#include "../../util/pair.h"
#include "../Hash.h"
#include "../util/utils.h"

#ifdef PMEM
#include <libpmemobj.h>
#endif

namespace extendible_16 {
// #define COUNTING 1
//#define PREALLOC 1

#ifdef AU
uint64_t space_num = 0;
uint64_t entry_num = 0;
uint64_t segment_num = 0;
uint64_t real_seg_num = 0;
# endif

template <class T>
struct _Pair {
  T key;
  Value_t value;
};

const uint32_t lockSet = ((uint32_t)1 << 31);
const uint32_t lockMask = ((uint32_t)1 << 31) - 1;
const int overflowSet = 1 << 4;
const int countMask = (1 << 4) - 1;
const uint64_t tailMask = (1UL << 56) - 1;
const uint64_t headerMask = ((1UL << 8) - 1) << 56;
const uint8_t overflowBitmapMask = (1 << 4) - 1;

constexpr size_t k_PairSize = 16;  // a k-v _Pair with a bit
constexpr size_t kNumPairPerBucket =
    14; /* it is determined by the usage of the fingerprint*/
constexpr size_t kFingerBits = 8;
constexpr size_t kMask = (1 << kFingerBits) - 1;
const constexpr size_t kNumBucket =
    64; /* the number of normal buckets in one segment 64*/
constexpr size_t stashBucket =
    2; /* the number of stash buckets in one segment*/
constexpr int allocMask = (1 << kNumPairPerBucket) - 1;
constexpr size_t bucketMask = ((1 << (int)log2(kNumBucket)) - 1);
constexpr size_t stashMask = (1 << (int)log2(stashBucket)) - 1;
constexpr uint8_t stashHighMask = ~((uint8_t)stashMask);

#define BUCKET_INDEX(hash) ((hash >> kFingerBits) & bucketMask)
#define GET_COUNT(var) ((var)&countMask)
#define GET_MEMBER(var) (((var) >> 4) & allocMask)
#define GET_INVERSE_MEMBER(var) ((~((var) >> 4)) & allocMask)
#define GET_BITMAP(var) ((var) >> 18)

inline bool var_compare(char *str1, char *str2, int len1, int len2) {
  if (len1 != len2) return false;
  return !memcmp(str1, str2, len1);
}

template <class T>
struct Bucket {
  inline int find_empty_slot() {
    if (GET_COUNT(bitmap) == kNumPairPerBucket) {
      return -1;
    }
    auto mask = ~(GET_BITMAP(bitmap));
    return __builtin_ctz(mask);
  }

  /*true indicates overflow, needs extra check in the stash*/
  inline bool test_overflow() { return overflowCount; }

  inline bool test_stash_check() { return (overflowBitmap & overflowSet); }

  inline void clear_stash_check() {
    overflowBitmap = overflowBitmap & (~overflowSet);
  }

  inline void set_indicator(uint8_t meta_hash, Bucket<T> *neighbor,
                            uint8_t pos) {
    int mask = overflowBitmap & overflowBitmapMask;
    mask = ~mask;
    auto index = __builtin_ctz(mask);

    if (index < 4) {
      finger_array[14 + index] = meta_hash;
      overflowBitmap = ((uint8_t)(1 << index) | overflowBitmap);
      overflowIndex =
          (overflowIndex & (~(3 << (index * 2)))) | (pos << (index * 2));
    } else {
      mask = neighbor->overflowBitmap & overflowBitmapMask;
      mask = ~mask;
      index = __builtin_ctz(mask);
      if (index < 4) {
        neighbor->finger_array[14 + index] = meta_hash;
        neighbor->overflowBitmap =
            ((uint8_t)(1 << index) | neighbor->overflowBitmap);
        neighbor->overflowMember =
            ((uint8_t)(1 << index) | neighbor->overflowMember);
        neighbor->overflowIndex =
            (neighbor->overflowIndex & (~(3 << (index * 2)))) |
            (pos << (index * 2));
      } else { /*overflow, increase count*/
        overflowCount++;
      }
    }
    overflowBitmap = overflowBitmap | overflowSet;
  }

  /*both clear this bucket and its neighbor bucket*/
  inline void unset_indicator(uint8_t meta_hash, Bucket<T> *neighbor, T key,
                              uint64_t pos) {
    /*also needs to ensure that this meta_hash must belongs to other bucket*/
    bool clear_success = false;
    int mask1 = overflowBitmap & overflowBitmapMask;
    for (int i = 0; i < 4; ++i) {
      if (CHECK_BIT(mask1, i) && (finger_array[14 + i] == meta_hash) &&
          (((1 << i) & overflowMember) == 0) &&
          (((overflowIndex >> (2 * i)) & stashMask) == pos)) {
        overflowBitmap = overflowBitmap & ((uint8_t)(~(1 << i)));
        overflowIndex = overflowIndex & (~(3 << (i * 2)));
        assert(((overflowIndex >> (i * 2)) & stashMask) == 0);
        clear_success = true;
        break;
      }
    }

    int mask2 = neighbor->overflowBitmap & overflowBitmapMask;
    if (!clear_success) {
      for (int i = 0; i < 4; ++i) {
        if (CHECK_BIT(mask2, i) &&
            (neighbor->finger_array[14 + i] == meta_hash) &&
            (((1 << i) & neighbor->overflowMember) != 0) &&
            (((neighbor->overflowIndex >> (2 * i)) & stashMask) == pos)) {
          neighbor->overflowBitmap =
              neighbor->overflowBitmap & ((uint8_t)(~(1 << i)));
          neighbor->overflowMember =
              neighbor->overflowMember & ((uint8_t)(~(1 << i)));
          neighbor->overflowIndex = neighbor->overflowIndex & (~(3 << (i * 2)));
          assert(((neighbor->overflowIndex >> (i * 2)) & stashMask) == 0);
          clear_success = true;
          break;
        }
      }
    }

    if (!clear_success) {
      overflowCount--;
    }

    mask1 = overflowBitmap & overflowBitmapMask;
    mask2 = neighbor->overflowBitmap & overflowBitmapMask;
    if (((mask1 & (~overflowMember)) == 0) && (overflowCount == 0) &&
        ((mask2 & neighbor->overflowMember) == 0)) {
      clear_stash_check();
    }
  }

  int unique_check(uint8_t meta_hash, T key, Bucket<T> *neighbor,
                   Bucket<T> *stash) {
    if ((check_and_get(meta_hash, key, false) != NONE) ||
        (neighbor->check_and_get(meta_hash, key, true) != NONE)) {
      return -1;
    }

    if (test_stash_check()) {
      auto test_stash = false;
      if (test_overflow()) {
        test_stash = true;
      } else {
        int mask = overflowBitmap & overflowBitmapMask;
        if (mask != 0) {
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) && (finger_array[14 + i] == meta_hash) &&
                (((1 << i) & overflowMember) == 0)) {
              test_stash = true;
              goto STASH_CHECK;
            }
          }
        }

        mask = neighbor->overflowBitmap & overflowBitmapMask;
        if (mask != 0) {
          for (int i = 0; i < 4; ++i) {
            if (CHECK_BIT(mask, i) &&
                (neighbor->finger_array[14 + i] == meta_hash) &&
                (((1 << i) & neighbor->overflowMember) != 0)) {
              test_stash = true;
              break;
            }
          }
        }
      }
      STASH_CHECK:
      if (test_stash == true) {
        for (int i = 0; i < stashBucket; ++i) {
          Bucket *curr_bucket = stash + i;
          if (curr_bucket->check_and_get(meta_hash, key, false) != NONE) {
            return -1;
          }
        }
      }
    }
    return 0;
  }

  inline int get_current_mask() {
    int mask = GET_BITMAP(bitmap) & GET_INVERSE_MEMBER(bitmap);
    return mask;
  }

  Value_t check_and_get(uint8_t meta_hash, T key, bool probe) {
    int mask = 0;
    SSE_CMP8(finger_array, meta_hash);
    if (!probe) {
      mask = mask & GET_BITMAP(bitmap) & (~GET_MEMBER(bitmap));
    } else {
      mask = mask & GET_BITMAP(bitmap) & GET_MEMBER(bitmap);
    }

    if (mask == 0) {
      return NONE;
    }

    if constexpr (std::is_pointer_v<T>) {
      /* variable-length key*/
      string_key *_key = reinterpret_cast<string_key *>(key);
      for (int i = 0; i < 14; i += 1) {
        if (CHECK_BIT(mask, i) &&
            (var_compare((reinterpret_cast<string_key *>(_[i].key))->key,
                         _key->key,
                         (reinterpret_cast<string_key *>(_[i].key))->length,
                         _key->length))) {
          return _[i].value;
        }
      }
    } else {
      /*fixed-length key*/
      /*loop unrolling*/
      for (int i = 0; i < 12; i += 4) {
        if (CHECK_BIT(mask, i) && (_[i].key == key)) {
          return _[i].value;
        }

        if (CHECK_BIT(mask, i + 1) && (_[i + 1].key == key)) {
          return _[i + 1].value;
        }

        if (CHECK_BIT(mask, i + 2) && (_[i + 2].key == key)) {
          return _[i + 2].value;
        }

        if (CHECK_BIT(mask, i + 3) && (_[i + 3].key == key)) {
          return _[i + 3].value;
        }
      }

      if (CHECK_BIT(mask, 12) && (_[12].key == key)) {
        return _[12].value;
      }

      if (CHECK_BIT(mask, 13) && (_[13].key == key)) {
        return _[13].value;
      }
    }
    return NONE;
  }

  inline void set_hash(int index, uint8_t meta_hash, bool probe) {
    finger_array[index] = meta_hash;
    uint32_t new_bitmap = bitmap | (1 << (index + 18));
    if (probe) {
      new_bitmap = new_bitmap | (1 << (index + 4));
    }
    new_bitmap += 1;
    bitmap = new_bitmap;
  }

  inline uint8_t get_hash(int index) { return finger_array[index]; }

  inline void unset_hash(int index, bool nt_flush = false) {
    uint32_t new_bitmap =
        bitmap & (~(1 << (index + 18))) & (~(1 << (index + 4)));
    assert(GET_COUNT(bitmap) <= kNumPairPerBucket);
    assert(GET_COUNT(bitmap) > 0);
    new_bitmap -= 1;
      bitmap = new_bitmap;
    }

  inline void get_lock() {
    uint32_t new_value = 0;
    uint32_t old_value = 0;
    do {
      while (true) {
        old_value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
        if (!(old_value & lockSet)) {
          old_value &= lockMask;
          break;
        }
      }
      new_value = old_value | lockSet;
    } while (!CAS(&version_lock, &old_value, new_value));
  }

  inline bool try_get_lock() {
    uint32_t v = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    if (v & lockSet) {
      return false;
    }
    auto old_value = v & lockMask;
    auto new_value = v | lockSet;
    return CAS(&version_lock, &old_value, new_value);
  }

  inline void release_lock() {
    uint32_t v = version_lock;
    __atomic_store_n(&version_lock, v + 1 - lockSet, __ATOMIC_RELEASE);
  }

  /*if the lock is set, return true*/
  inline bool test_lock_set(uint32_t &version) {
    version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    return (version & lockSet) != 0;
  }

  // test whether the version has change, if change, return true
  inline bool test_lock_version_change(uint32_t old_version) {
    auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
    return (old_version != value);
  }

  int Insert(T key, Value_t value, uint8_t meta_hash, bool probe) {
    auto slot = find_empty_slot();
    assert(slot < kNumPairPerBucket);
    if (slot == -1) {
      return -1;
    }
    _[slot].value = value;
    _[slot].key = key;
    set_hash(slot, meta_hash, probe);
    return 0;
  }

  /*if delete success, then return 0, else return -1*/
  int Delete(T key, uint8_t meta_hash, bool probe) {
    /*do the simd and check the key, then do the delete operation*/
    int mask = 0;
    SSE_CMP8(finger_array, meta_hash);
    if (!probe) {
      mask = mask & GET_BITMAP(bitmap) & (~GET_MEMBER(bitmap));
    } else {
      mask = mask & GET_BITMAP(bitmap) & GET_MEMBER(bitmap);
    }

    /*loop unrolling*/
    if constexpr (std::is_pointer_v<T>) {
      string_key *_key = reinterpret_cast<string_key *>(key);
      /*loop unrolling*/
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) &&
              (var_compare((reinterpret_cast<string_key *>(_[i].key))->key,
                           _key->key,
                           (reinterpret_cast<string_key *>(_[i].key))->length,
                           _key->length))) {
            unset_hash(i, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 1) &&
              (var_compare(
                  reinterpret_cast<string_key *>(_[i + 1].key)->key, _key->key,
                  (reinterpret_cast<string_key *>(_[i + 1].key))->length,
                  _key->length))) {
            unset_hash(i + 1, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 2) &&
              (var_compare(
                  reinterpret_cast<string_key *>(_[i + 2].key)->key, _key->key,
                  (reinterpret_cast<string_key *>(_[i + 2].key))->length,
                  _key->length))) {
            unset_hash(i + 2, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 3) &&
              (var_compare(
                  reinterpret_cast<string_key *>(_[i + 3].key)->key, _key->key,
                  (reinterpret_cast<string_key *>(_[i + 3].key))->length,
                  _key->length))) {
            unset_hash(i + 3, false);
            return 0;
          }
        }

        if (CHECK_BIT(mask, 12) &&
            (var_compare(reinterpret_cast<string_key *>(_[12].key)->key,
                         _key->key,
                         (reinterpret_cast<string_key *>(_[12].key))->length,
                         _key->length))) {
          unset_hash(12, false);
          return 0;
        }

        if (CHECK_BIT(mask, 13) &&
            (var_compare(reinterpret_cast<string_key *>(_[13].key)->key,
                         _key->key,
                         (reinterpret_cast<string_key *>(_[13].key))->length,
                         _key->length))) {
          unset_hash(13, false);
          return 0;
        }
      }

    } else {
      if (mask != 0) {
        for (int i = 0; i < 12; i += 4) {
          if (CHECK_BIT(mask, i) && (_[i].key == key)) {
            unset_hash(i, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 1) && (_[i + 1].key == key)) {
            unset_hash(i + 1, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 2) && (_[i + 2].key == key)) {
            unset_hash(i + 2, false);
            return 0;
          }

          if (CHECK_BIT(mask, i + 3) && (_[i + 3].key == key)) {
            unset_hash(i + 3, false);
            return 0;
          }
        }

        if (CHECK_BIT(mask, 12) && (_[12].key == key)) {
          unset_hash(12, false);
          return 0;
        }

        if (CHECK_BIT(mask, 13) && (_[13].key == key)) {
          unset_hash(13, false);
          return 0;
        }
      }
    }
    return -1;
  }

  int Insert_with_noflush(T key, Value_t value, uint8_t meta_hash, bool probe) {
    auto slot = find_empty_slot();
    /* this branch can be removed*/
    assert(slot < kNumPairPerBucket);
    if (slot == -1) {
      std::cout << "Cannot find the empty slot, for key " << key << std::endl;
      return -1;
    }
    _[slot].value = value;
    _[slot].key = key;
    set_hash(slot, meta_hash, probe);
    return 0;
  }

  void Insert_displace(T key, Value_t value, uint8_t meta_hash, int slot,
                       bool probe) {
    _[slot].value = value;
    _[slot].key = key;
    set_hash(slot, meta_hash, probe);
  }

  void Insert_displace_with_noflush(T key, Value_t value, uint8_t meta_hash,
                                    int slot, bool probe) {
    _[slot].value = value;
    _[slot].key = key;
    set_hash(slot, meta_hash, probe);
  }

  /* Find the displacment element in this bucket*/
  inline int Find_org_displacement() {
    uint32_t mask = GET_INVERSE_MEMBER(bitmap);
    if (mask == 0) {
      return -1;
    }
    return __builtin_ctz(mask);
  }

  /*find element that it is in the probe*/
  inline int Find_probe_displacement() {
    uint32_t mask = GET_MEMBER(bitmap);
    if (mask == 0) {
      return -1;
    }
    return __builtin_ctz(mask);
  }

  inline void resetLock() { version_lock = 0; }

  inline void resetOverflowFP() {
    overflowBitmap = 0;
    overflowIndex = 0;
    overflowMember = 0;
    overflowCount = 0;
    clear_stash_check();
  }

  uint32_t version_lock;
  uint32_t bitmap;          // allocation bitmap + pointer bitmap + counter
  uint8_t finger_array[18]; /*only use the first 14 bytes, can be accelerated by
                               SSE instruction,0-13 for finger, 14-17 for
                               overflowed*/
  uint8_t overflowBitmap;
  uint8_t overflowIndex;
  uint8_t overflowMember; /*overflowmember indicates membership of the overflow
                             fingerprint*/
  uint8_t overflowCount;
  uint8_t unused[2];

  _Pair<T> _[kNumPairPerBucket];
};

template <class T>
struct Table;

template <class T>
struct Directory {
  typedef Table<T> *table_p;
  uint32_t global_depth;
  uint32_t version;
  uint32_t depth_count;
  table_p _[0];

  Directory(size_t capacity, size_t _version) {
    version = _version;
    global_depth = static_cast<size_t>(log2(capacity));
    depth_count = 0;
  }

  static void New(void **dir, size_t capacity, size_t version) {
    Allocator::ZAllocate(dir, kCacheLineSize, sizeof(Directory<T>) + sizeof(table_p) * capacity);
    new (*dir) Directory(capacity, version);
#ifdef AU
    segment_num += (sizeof(Directory<T>) + sizeof(table_p) * capacity + 15)/16;
#endif
  }
};


/* the segment class*/
template <class T>
struct Table {
  static void New(void **tbl, size_t depth, void *pp) {
    Allocator::ZAllocate(tbl, kCacheLineSize, sizeof(Table<T>));
    auto table_ptr = reinterpret_cast<Table<T> *>(*tbl);
    table_ptr->local_depth = depth;
    table_ptr->next = pp;
    table_ptr->state = -3; /*NEW*/
    int sumBucket = kNumBucket + stashBucket;
    for (int i = 0; i < sumBucket; ++i) {
      auto curr_bucket = table_ptr->bucket + i;
      memset(curr_bucket, 0, 64);
    }
  };
  ~Table(void) {}

  bool Acquire_and_verify(size_t _pattern) {
    bucket->get_lock();
    if (pattern != _pattern) {
      bucket->release_lock();
      return false;
    } else {
      return true;
    }
  }

  void Acquire_remaining_locks() {
    for (int i = 1; i < kNumBucket; ++i) {
      auto curr_bucket = bucket + i;
      curr_bucket->get_lock();
    }
  }

  void Release_all_locks() {
    for (int i = 0; i < kNumBucket; ++i) {
      auto curr_bucket = bucket + i;
      curr_bucket->release_lock();
    }
  }

  int Insert(T key, Value_t value, size_t key_hash, uint8_t meta_hash,
             Directory<T> **);
  void Insert4split(T key, Value_t value, size_t key_hash, uint8_t meta_hash);
  void Insert4splitWithCheck(T key, Value_t value, size_t key_hash,
                             uint8_t meta_hash); /*with uniqueness check*/
  void Insert4merge(T key, Value_t value, size_t key_hash, uint8_t meta_hash,
                    bool flag = false);
  Table<T> *Split(size_t);
  void HelpSplit(Table<T> *);
  void Merge(Table<T> *, bool flag = false);
  int Delete(T key, size_t key_hash, uint8_t meta_hash, Directory<T> **_dir);

  int Next_displace(Bucket<T> *target, Bucket<T> *neighbor,
                    Bucket<T> *next_neighbor, T key, Value_t value,
                    uint8_t meta_hash) {
    int displace_index = neighbor->Find_org_displacement();
    if ((GET_COUNT(next_neighbor->bitmap) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      next_neighbor->Insert(neighbor->_[displace_index].key,
                            neighbor->_[displace_index].value,
                            neighbor->finger_array[displace_index], true);
      next_neighbor->release_lock();

      neighbor->unset_hash(displace_index);
      neighbor->Insert_displace(key, value, meta_hash, displace_index, true);
      neighbor->release_lock();

      target->release_lock();
#ifdef COUNTING
      __sync_fetch_and_add(&number, 1);
#endif
      return 0;
    }
    return -1;
  }

  int Prev_displace(Bucket<T> *target, Bucket<T> *prev_neighbor,
                    Bucket<T> *neighbor, T key, Value_t value,
                    uint8_t meta_hash) {
    int displace_index = target->Find_probe_displacement();
    if ((GET_COUNT(prev_neighbor->bitmap) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      prev_neighbor->Insert(target->_[displace_index].key,
                            target->_[displace_index].value,
                            target->finger_array[displace_index], false);
      prev_neighbor->release_lock();

      target->unset_hash(displace_index);
      target->Insert_displace(key, value, meta_hash, displace_index, false);
      target->release_lock();

      neighbor->release_lock();
#ifdef COUNTING
      __sync_fetch_and_add(&number, 1);
#endif
      return 0;
    }
    return -1;
  }

  int Stash_insert(Bucket<T> *target, Bucket<T> *neighbor, T key, Value_t value,
                   uint8_t meta_hash, int stash_pos) {
    for (int i = 0; i < stashBucket; ++i) {
      Bucket<T> *curr_bucket =
          bucket + kNumBucket + ((stash_pos + i) & stashMask);
      if (GET_COUNT(curr_bucket->bitmap) < kNumPairPerBucket) {
        curr_bucket->Insert(key, value, meta_hash, false);

        target->set_indicator(meta_hash, neighbor, (stash_pos + i) & stashMask);
#ifdef COUNTING
        __sync_fetch_and_add(&number, 1);
#endif
        return 0;
      }
    }
    return -1;
  }

  void recoverMetadata() {
    Bucket<T> *curr_bucket, *neighbor_bucket;
    /*reset the lock and overflow meta-data*/
    uint64_t knumber = 0;
    for (int i = 0; i < kNumBucket; ++i) {
      curr_bucket = bucket + i;
      curr_bucket->resetLock();
      curr_bucket->resetOverflowFP();
      neighbor_bucket = bucket + ((i + 1) & bucketMask);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        int mask = curr_bucket->get_current_mask();
        if (CHECK_BIT(mask, j) && (neighbor_bucket->check_and_get(
            curr_bucket->finger_array[j],
            curr_bucket->_[j].key, true) != NONE)) {
          curr_bucket->unset_hash(j);
        }
      }

#ifdef COUNTING
      knumber += __builtin_popcount(GET_BITMAP(curr_bucket->bitmap));
#endif
    }

    /*scan the stash buckets and re-insert the overflow FP to initial buckets*/
    for (int i = 0; i < stashBucket; ++i) {
      curr_bucket = bucket + kNumBucket + i;
      curr_bucket->resetLock();
#ifdef COUNTING
      knumber += __builtin_popcount(GET_BITMAP(curr_bucket->bitmap));
#endif
      uint64_t key_hash;
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }
          /*compute the initial bucket*/
          auto bucket_ix = BUCKET_INDEX(key_hash);
          auto meta_hash = ((uint8_t)(key_hash & kMask));  // the last 8 bits
          auto org_bucket = bucket + bucket_ix;
          auto neighbor_bucket = bucket + ((bucket_ix + 1) & bucketMask);
          org_bucket->set_indicator(meta_hash, neighbor_bucket, i);
        }
      }
    }
#ifdef COUNTING
    number = knumber;
#endif
    /* No need to flush these meta-data because persistent or not does not
     * influence the correctness*/
  }

  char dummy[48];
  Bucket<T> bucket[kNumBucket + stashBucket];
  size_t local_depth;
  size_t pattern;
  int number;
  void* next;
  int state; /*-1 means this bucket is merging, -2 means this bucket is
                splitting (SPLITTING), 0 meanning normal bucket, -3 means new
                bucket (NEW)*/
};

/* it needs to verify whether this bucket has been deleted...*/
template <class T>
int Table<T>::Insert(T key, Value_t value, size_t key_hash, uint8_t meta_hash,
                     Directory<T> **_dir) {
  RETRY:
  /*we need to first do the locking and then do the verify*/
  auto y = BUCKET_INDEX(key_hash);
  Bucket<T> *target = bucket + y;
  Bucket<T> *neighbor = bucket + ((y + 1) & bucketMask);
  target->get_lock();
  if (!neighbor->try_get_lock()) {
    target->release_lock();
    return -2;
  }

  auto old_sa = *_dir;
  auto x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
  if (reinterpret_cast<Table<T> *>(reinterpret_cast<uint64_t>(old_sa->_[x]) &
                                   tailMask) != this) {
    neighbor->release_lock();
    target->release_lock();
    return -2;
  }

  /*unique check, needs to check 2 hash table*/
  auto ret =
      target->unique_check(meta_hash, key, neighbor, bucket + kNumBucket);
  if (ret == -1) {
    neighbor->release_lock();
    target->release_lock();
    return -3; /* duplicate insert*/
  }

  if (((GET_COUNT(target->bitmap)) == kNumPairPerBucket) &&
      ((GET_COUNT(neighbor->bitmap)) == kNumPairPerBucket)) {
    Bucket<T> *next_neighbor = bucket + ((y + 2) & bucketMask);
    // Next displacement
    if (!next_neighbor->try_get_lock()) {
      neighbor->release_lock();
      target->release_lock();
      return -2;
    }
    auto ret =
        Next_displace(target, neighbor, next_neighbor, key, value, meta_hash);
    if (ret == 0) {
      return 0;
    }
    next_neighbor->release_lock();

    Bucket<T> *prev_neighbor;
    int prev_index;
    if (y == 0) {
      prev_neighbor = bucket + kNumBucket - 1;
      prev_index = kNumBucket - 1;
    } else {
      prev_neighbor = bucket + y - 1;
      prev_index = y - 1;
    }
    if (!prev_neighbor->try_get_lock()) {
      target->release_lock();
      neighbor->release_lock();
      return -2;
    }

    ret = Prev_displace(target, prev_neighbor, neighbor, key, value, meta_hash);
    if (ret == 0) {
      return 0;
    }

    Bucket<T> *stash = bucket + kNumBucket;
    if (!stash->try_get_lock()) {
      neighbor->release_lock();
      target->release_lock();
      prev_neighbor->release_lock();
      return -2;
    }
    ret = Stash_insert(target, neighbor, key, value, meta_hash, y & stashMask);

    stash->release_lock();
    neighbor->release_lock();
    target->release_lock();
    prev_neighbor->release_lock();
    return ret;
  }

  /* the fp+bitmap are persisted after releasing the lock of one bucket but
   * still guarantee the correctness of avoidance of "use-before-flush" since
   * the search operation could only proceed only if both target bucket and
   * probe bucket are released
   */
  if (GET_COUNT(target->bitmap) <= GET_COUNT(neighbor->bitmap)) {
    target->Insert(key, value, meta_hash, false);
    target->release_lock();

    neighbor->release_lock();
  } else {
    neighbor->Insert(key, value, meta_hash, true);
    neighbor->release_lock();

    target->release_lock();
  }
#ifdef COUNTING
  __sync_fetch_and_add(&number, 1);
#endif
  return 0;
}

template <class T>
void Table<T>::Insert4splitWithCheck(T key, Value_t value, size_t key_hash,
                                     uint8_t meta_hash) {
  auto y = BUCKET_INDEX(key_hash);
  Bucket<T> *target = bucket + y;
  Bucket<T> *neighbor = bucket + ((y + 1) & bucketMask);
  auto ret =
      target->unique_check(meta_hash, key, neighbor, bucket + kNumBucket);
  if (ret == -1) return;
  Bucket<T> *insert_target;
  bool probe = false;
  if (GET_COUNT(target->bitmap) <= GET_COUNT(neighbor->bitmap)) {
    insert_target = target;
  } else {
    insert_target = neighbor;
    probe = true;
  }

  /*some bucket may be overflowed?*/
  if (GET_COUNT(insert_target->bitmap) < kNumPairPerBucket) {
    insert_target->_[GET_COUNT(insert_target->bitmap)].key = key;
    insert_target->_[GET_COUNT(insert_target->bitmap)].value = value;
    insert_target->set_hash(GET_COUNT(insert_target->bitmap), meta_hash, probe);
#ifdef COUNTING
    ++number;
#endif
  } else {
    /*do the displacement or insertion in the stash*/
    Bucket<T> *next_neighbor = bucket + ((y + 2) & bucketMask);
    int displace_index;
    displace_index = neighbor->Find_org_displacement();
    if (((GET_COUNT(next_neighbor->bitmap)) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      next_neighbor->Insert_with_noflush(
          neighbor->_[displace_index].key, neighbor->_[displace_index].value,
          neighbor->finger_array[displace_index], true);
      neighbor->unset_hash(displace_index);
      neighbor->Insert_displace_with_noflush(key, value, meta_hash,
                                             displace_index, true);
#ifdef COUNTING
      ++number;
#endif
      return;
    }
    Bucket<T> *prev_neighbor;
    int prev_index;
    if (y == 0) {
      prev_neighbor = bucket + kNumBucket - 1;
      prev_index = kNumBucket - 1;
    } else {
      prev_neighbor = bucket + y - 1;
      prev_index = y - 1;
    }

    displace_index = target->Find_probe_displacement();
    if (((GET_COUNT(prev_neighbor->bitmap)) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      prev_neighbor->Insert_with_noflush(
          target->_[displace_index].key, target->_[displace_index].value,
          target->finger_array[displace_index], false);
      target->unset_hash(displace_index);
      target->Insert_displace_with_noflush(key, value, meta_hash,
                                           displace_index, false);
#ifdef COUNTING
      ++number;
#endif
      return;
    }

    Stash_insert(target, neighbor, key, value, meta_hash, y & stashMask);
  }
}

/*the insert needs to be perfectly balanced, not destory the power of balance*/
template <class T>
void Table<T>::Insert4split(T key, Value_t value, size_t key_hash,
                            uint8_t meta_hash) {
  auto y = BUCKET_INDEX(key_hash);
  Bucket<T> *target = bucket + y;
  Bucket<T> *neighbor = bucket + ((y + 1) & bucketMask);
  Bucket<T> *insert_target;
  bool probe = false;
  if (GET_COUNT(target->bitmap) <= GET_COUNT(neighbor->bitmap)) {
    insert_target = target;
  } else {
    insert_target = neighbor;
    probe = true;
  }

  /*some bucket may be overflowed?*/
  if (GET_COUNT(insert_target->bitmap) < kNumPairPerBucket) {
    insert_target->_[GET_COUNT(insert_target->bitmap)].key = key;
    insert_target->_[GET_COUNT(insert_target->bitmap)].value = value;
    insert_target->set_hash(GET_COUNT(insert_target->bitmap), meta_hash, probe);
#ifdef COUNTING
    ++number;
#endif
  } else {
    /*do the displacement or insertion in the stash*/
    Bucket<T> *next_neighbor = bucket + ((y + 2) & bucketMask);
    int displace_index;
    displace_index = neighbor->Find_org_displacement();
    if (((GET_COUNT(next_neighbor->bitmap)) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      next_neighbor->Insert_with_noflush(
          neighbor->_[displace_index].key, neighbor->_[displace_index].value,
          neighbor->finger_array[displace_index], true);
      neighbor->unset_hash(displace_index);
      neighbor->Insert_displace_with_noflush(key, value, meta_hash,
                                             displace_index, true);
#ifdef COUNTING
      ++number;
#endif
      return;
    }
    Bucket<T> *prev_neighbor;
    int prev_index;
    if (y == 0) {
      prev_neighbor = bucket + kNumBucket - 1;
      prev_index = kNumBucket - 1;
    } else {
      prev_neighbor = bucket + y - 1;
      prev_index = y - 1;
    }

    displace_index = target->Find_probe_displacement();
    if (((GET_COUNT(prev_neighbor->bitmap)) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      prev_neighbor->Insert_with_noflush(
          target->_[displace_index].key, target->_[displace_index].value,
          target->finger_array[displace_index], false);
      target->unset_hash(displace_index);
      target->Insert_displace_with_noflush(key, value, meta_hash,
                                           displace_index, false);
#ifdef COUNTING
      ++number;
#endif
      return;
    }

    Stash_insert(target, neighbor, key, value, meta_hash, y & stashMask);
  }
}

template <class T>
void Table<T>::Insert4merge(T key, Value_t value, size_t key_hash,
                            uint8_t meta_hash, bool unique_check_flag) {
  auto y = BUCKET_INDEX(key_hash);
  Bucket<T> *target = bucket + y;
  Bucket<T> *neighbor = bucket + ((y + 1) & bucketMask);

  if (unique_check_flag) {
    auto ret =
        target->unique_check(meta_hash, key, neighbor, bucket + kNumBucket);
    if (ret == -1) return;
  }

  Bucket<T> *insert_target;
  bool probe = false;
  if (GET_COUNT(target->bitmap) <= GET_COUNT(neighbor->bitmap)) {
    insert_target = target;
  } else {
    insert_target = neighbor;
    probe = true;
  }

  /*some bucket may be overflowed?*/
  if (GET_COUNT(insert_target->bitmap) < kNumPairPerBucket) {
    insert_target->Insert(key, value, meta_hash, probe);
#ifdef COUNTING
    ++number;
#endif
  } else {
    /*do the displacement or insertion in the stash*/
    Bucket<T> *next_neighbor = bucket + ((y + 2) & bucketMask);
    int displace_index;
    displace_index = neighbor->Find_org_displacement();
    if (((GET_COUNT(next_neighbor->bitmap)) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      next_neighbor->Insert_with_noflush(
          neighbor->_[displace_index].key, neighbor->_[displace_index].value,
          neighbor->finger_array[displace_index], true);
      neighbor->unset_hash(displace_index);
      neighbor->Insert_displace_with_noflush(key, value, meta_hash,
                                             displace_index, true);
#ifdef COUNTING
      ++number;
#endif
      return;
    }
    Bucket<T> *prev_neighbor;
    int prev_index;
    if (y == 0) {
      prev_neighbor = bucket + kNumBucket - 1;
      prev_index = kNumBucket - 1;
    } else {
      prev_neighbor = bucket + y - 1;
      prev_index = y - 1;
    }

    displace_index = target->Find_probe_displacement();
    if (((GET_COUNT(prev_neighbor->bitmap)) != kNumPairPerBucket) &&
        (displace_index != -1)) {
      prev_neighbor->Insert_with_noflush(
          target->_[displace_index].key, target->_[displace_index].value,
          target->finger_array[displace_index], false);
      target->unset_hash(displace_index);
      target->Insert_displace_with_noflush(key, value, meta_hash,
                                           displace_index, false);
#ifdef COUNTING
      ++number;
#endif
      return;
    }

    Stash_insert(target, neighbor, key, value, meta_hash, y & stashMask);
  }
}

template <class T>
void Table<T>::HelpSplit(Table<T> *next_table) {
  size_t new_pattern = (pattern << 1) + 1;
  size_t old_pattern = pattern << 1;

  size_t key_hash;
  uint32_t invalid_array[kNumBucket + stashBucket];
  for (int i = 0; i < kNumBucket; ++i) {
    auto *curr_bucket = bucket + i;
    auto mask = GET_BITMAP(curr_bucket->bitmap);
    uint32_t invalid_mask = 0;
    for (int j = 0; j < kNumPairPerBucket; ++j) {
      if (CHECK_BIT(mask, j)) {
        if constexpr (std::is_pointer_v<T>) {
          auto curr_key = curr_bucket->_[j].key;
          key_hash = h(curr_key->key, curr_key->length);
        } else {
          key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
        }

        if ((key_hash >> (64 - local_depth - 1)) == new_pattern) {
          invalid_mask = invalid_mask | (1 << j);
          next_table->Insert4splitWithCheck(curr_bucket->_[j].key,
                                            curr_bucket->_[j].value, key_hash,
                                            curr_bucket->finger_array[j]);
#ifdef COUNTING
          number--;
#endif
        }
      }
    }
    invalid_array[i] = invalid_mask;
  }

  for (int i = 0; i < stashBucket; ++i) {
    auto *curr_bucket = bucket + kNumBucket + i;
    auto mask = GET_BITMAP(curr_bucket->bitmap);
    uint32_t invalid_mask = 0;
    for (int j = 0; j < kNumPairPerBucket; ++j) {
      if (CHECK_BIT(mask, j)) {
        if constexpr (std::is_pointer_v<T>) {
          auto curr_key = curr_bucket->_[j].key;
          key_hash = h(curr_key->key, curr_key->length);
        } else {
          key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
        }
        if ((key_hash >> (64 - local_depth - 1)) == new_pattern) {
          invalid_mask = invalid_mask | (1 << j);
          next_table->Insert4splitWithCheck(curr_bucket->_[j].key,
                                            curr_bucket->_[j].value, key_hash,
                                            curr_bucket->finger_array[j]);
          auto bucket_ix = BUCKET_INDEX(key_hash);
          auto org_bucket = bucket + bucket_ix;
          auto neighbor_bucket = bucket + ((bucket_ix + 1) & bucketMask);
          org_bucket->unset_indicator(curr_bucket->finger_array[j],
                                      neighbor_bucket, curr_bucket->_[j].key,
                                      i);
#ifdef COUNTING
          number--;
#endif
        }
      }
    }
    invalid_array[kNumBucket + i] = invalid_mask;
  }
  next_table->pattern = new_pattern;
  pattern = old_pattern;

  size_t sumBucket = kNumBucket + stashBucket;
  for (int i = 0; i < sumBucket; ++i) {
    auto curr_bucket = bucket + i;
    curr_bucket->bitmap = curr_bucket->bitmap & (~(invalid_array[i] << 18)) &
                          (~(invalid_array[i] << 4));
    uint32_t count = __builtin_popcount(invalid_array[i]);
    curr_bucket->bitmap = curr_bucket->bitmap - count;
  }

}

template <class T>
Table<T> *Table<T>::Split(size_t _key_hash) {
  size_t new_pattern = (pattern << 1) + 1;
  size_t old_pattern = pattern << 1;

  for (int i = 1; i < kNumBucket; ++i) {
    (bucket + i)->get_lock();
  }
  state = -2; /*means the start of the split process*/
  Table<T>::New(&next, local_depth + 1, next);
  Table<T> *next_table = reinterpret_cast<Table<T> *>(next);

  next_table->state = -2;
  next_table->bucket
      ->get_lock(); /* get the first lock of the new bucket to avoid it
                 is operated(split or merge) by other threads*/
  size_t key_hash;
  uint32_t invalid_array[kNumBucket + stashBucket];
  for (int i = 0; i < kNumBucket; ++i) {
    auto *curr_bucket = bucket + i;
    auto mask = GET_BITMAP(curr_bucket->bitmap);
    uint32_t invalid_mask = 0;
    for (int j = 0; j < kNumPairPerBucket; ++j) {
      if (CHECK_BIT(mask, j)) {
        if constexpr (std::is_pointer_v<T>) {
          auto curr_key = curr_bucket->_[j].key;
          key_hash = h(curr_key->key, curr_key->length);
        } else {
          key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
        }

        if ((key_hash >> (64 - local_depth - 1)) == new_pattern) {
          invalid_mask = invalid_mask | (1 << j);
          next_table->Insert4split(
              curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
              curr_bucket->finger_array[j]); /*this shceme may destory the
                                                balanced segment*/
          // curr_bucket->unset_hash(j);
#ifdef COUNTING
          number--;
#endif
        }
      }
    }
    invalid_array[i] = invalid_mask;
  }

  for (int i = 0; i < stashBucket; ++i) {
    auto *curr_bucket = bucket + kNumBucket + i;
    auto mask = GET_BITMAP(curr_bucket->bitmap);
    uint32_t invalid_mask = 0;
    for (int j = 0; j < kNumPairPerBucket; ++j) {
      if (CHECK_BIT(mask, j)) {
        if constexpr (std::is_pointer_v<T>) {
          auto curr_key = curr_bucket->_[j].key;
          key_hash = h(curr_key->key, curr_key->length);
        } else {
          key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
        }
        if ((key_hash >> (64 - local_depth - 1)) == new_pattern) {
          invalid_mask = invalid_mask | (1 << j);
          next_table->Insert4split(
              curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
              curr_bucket->finger_array[j]); /*this shceme may destory the
                                                balanced segment*/
          auto bucket_ix = BUCKET_INDEX(key_hash);
          auto org_bucket = bucket + bucket_ix;
          auto neighbor_bucket = bucket + ((bucket_ix + 1) & bucketMask);
          org_bucket->unset_indicator(curr_bucket->finger_array[j],
                                      neighbor_bucket, curr_bucket->_[j].key,
                                      i);
#ifdef COUNTING
          number--;
#endif
        }
      }
    }
    invalid_array[kNumBucket + i] = invalid_mask;
  }
  next_table->pattern = new_pattern;
  pattern = old_pattern;

  size_t sumBucket = kNumBucket + stashBucket;
  for (int i = 0; i < sumBucket; ++i) {
    auto curr_bucket = bucket + i;
    curr_bucket->bitmap = curr_bucket->bitmap & (~(invalid_array[i] << 18)) &
                          (~(invalid_array[i] << 4));
    uint32_t count = __builtin_popcount(invalid_array[i]);
    curr_bucket->bitmap = curr_bucket->bitmap - count;
  }

  return next_table;
}

template <class T>
void Table<T>::Merge(Table<T> *neighbor, bool unique_check_flag) {
  /*Restore the split/merge procedure*/
  if (unique_check_flag) {
    size_t key_hash;
    for (int i = 0; i < kNumBucket; ++i) {
      auto *curr_bucket = neighbor->bucket + i;
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }

          Insert4merge(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j],
                       true); /*this shceme may destory
                           the balanced segment*/
        }
      }
    }

    for (int i = 0; i < stashBucket; ++i) {
      auto *curr_bucket = neighbor->bucket + kNumBucket + i;
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }
          Insert4merge(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j]); /*this shceme may destory
                                                         the balanced segment*/
        }
      }
    }

  } else {
    size_t key_hash;
    for (int i = 0; i < kNumBucket; ++i) {
      auto *curr_bucket = neighbor->bucket + i;
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }

          Insert4merge(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j]); /*this shceme may destory
                                                         the balanced segment*/
        }
      }
    }

    /*split the stash bucket, the stash must be full, right?*/
    for (int i = 0; i < stashBucket; ++i) {
      auto *curr_bucket = neighbor->bucket + kNumBucket + i;
      auto mask = GET_BITMAP(curr_bucket->bitmap);
      for (int j = 0; j < kNumPairPerBucket; ++j) {
        if (CHECK_BIT(mask, j)) {
          if constexpr (std::is_pointer_v<T>) {
            auto curr_key = curr_bucket->_[j].key;
            key_hash = h(curr_key->key, curr_key->length);
          } else {
            key_hash = h(&(curr_bucket->_[j].key), sizeof(Key_t));
          }
          Insert4merge(curr_bucket->_[j].key, curr_bucket->_[j].value, key_hash,
                       curr_bucket->finger_array[j]); /*this shceme may destory
                                                         the balanced segment*/
        }
      }
    }
  }
}

template <class T>
class Finger_EH : public Hash<T> {
 public:
  Finger_EH(void);
  Finger_EH(size_t);
  ~Finger_EH(void);
  inline int Insert(T key, Value_t value);
#ifdef AU
  int Insert0(T key, Value_t value);
#endif
  int Insert(T key, Value_t value, bool);
  inline bool Delete(T);
  bool Delete(T, bool);
  inline Value_t Get(T);
  Value_t Get(T key, bool is_in_epoch);
  int Get(T key, Value_t *ret_val, size_t *number){
    return 0;
  }
  void Directory_Doubling(int x, Table<T> *new_b, Table<T> *old_b);
  void Directory_Merge_Update(Directory<T> *_sa, uint64_t key_hash,
                              Table<T> *left_seg);
  void Directory_Update(Directory<T> *_sa, int x, Table<T> *new_b,
                        Table<T> *old_b);
  int FindAnyway(T key);
  void ShutDown() {
    clean = true;
  }
  uint64_t getNumber() {
    std::cout << "The size of the bucket is " << sizeof(struct Bucket<T>) << std::endl;
    size_t _count = 0;
    size_t seg_count = 0;
    Directory<T> *seg = dir;
    Table<T> **dir_entry = seg->_;
    Table<T> *ss;
    auto global_depth = seg->global_depth;
    size_t depth_diff;
    int capacity = pow(2, global_depth);
    for (int i = 0; i < capacity;) {
      ss = reinterpret_cast<Table<T> *>(
          reinterpret_cast<uint64_t>(dir_entry[i]) & tailMask);
      depth_diff = global_depth - ss->local_depth;
      for (int i = 0; i < kNumBucket + stashBucket; i++)
      {
        auto curr_bucket = ss->bucket + i;
        _count += GET_COUNT(curr_bucket->bitmap);
      }
      _count += ss->number;
      seg_count++;
      i += pow(2, depth_diff);
    }

    ss = reinterpret_cast<Table<T> *>(reinterpret_cast<uint64_t>(dir_entry[0]) &
                                      tailMask);
    uint64_t verify_seg_count = 1;
    while (ss->next) {
      verify_seg_count++;
      ss = reinterpret_cast<Table<T> *>((ss->next));
    }
    std::cout << "seg_count = " << seg_count << std::endl;
    std::cout << "verify_seg_count = " << verify_seg_count << std::endl;
#ifdef AU
    double average_utility = ((entry_num + 1.0) * entry_num ) / (2.0 * space_num);
    std::cout << "seg number = " << real_seg_num << std::endl;
    std::cout << "entry num = " << entry_num << std::endl;
    std::cout << "space unit = " << space_num << std::endl;
    std::cout << "average utility = " << average_utility << std::endl;
#endif
#ifdef COUNTING
    std::cout << "#items = " << _count << std::endl;
    std::cout << "load_factor = " <<
           (double)_count / (seg_count * kNumPairPerBucket * (kNumBucket + 2)) << std::endl;
    std::cout << "Raw_Space: ",
           (double)(_count * 16) / (seg_count * sizeof(Table<T>)) << std::endl;
#endif
    return 0;
  }

  void recoverTable(Table<T> **target_table, size_t, size_t, Directory<T> *);
  void Recovery();

  inline int Test_Directory_Lock_Set(void) {
    uint32_t v = __atomic_load_n(&lock, __ATOMIC_ACQUIRE);
    return v & lockSet;
  }

  inline bool try_get_directory_read_lock(){
    uint32_t v = __atomic_load_n(&lock, __ATOMIC_ACQUIRE);
    uint32_t old_value = v & lockMask;
    auto new_value = ((v & lockMask) + 1) & lockMask;
    return CAS(&lock, &old_value, new_value);
  }

  inline void release_directory_read_lock(){
    SUB(&lock, 1);
  }

  void Lock_Directory(){
    uint32_t v = __atomic_load_n(&lock, __ATOMIC_ACQUIRE);
    uint32_t old_value = v & lockMask;
    uint32_t new_value = old_value | lockSet;

    while (!CAS(&lock, &old_value, new_value)) {
      old_value = old_value & lockMask;
      new_value = old_value | lockSet;
    }

    //wait until the readers all exit the critical section
    v = __atomic_load_n(&lock, __ATOMIC_ACQUIRE);
    while(v & lockMask){
      v = __atomic_load_n(&lock, __ATOMIC_ACQUIRE);
    }
  }

  // just set the lock as 0
  void Unlock_Directory(){
    __atomic_store_n(&lock, 0, __ATOMIC_RELEASE);
  }

  Directory<T> *dir;
  uint32_t lock; // the MSB is the lock bit; remaining bits are used as the counter
  uint64_t
      crash_version; /*when the crash version equals to 0Xff => set the crash
                        version as 0, set the version of all entries as 1*/
  bool clean;
  /* directory allocation will write to here first,
   * in oder to perform safe directory allocation
   * */
  void* back_dir;
};

template <class T>
Finger_EH<T>::Finger_EH(size_t initCap) {
  Directory<T>::New(&back_dir, initCap, 0);
  dir = reinterpret_cast<Directory<T> *>(back_dir);
  back_dir = nullptr;
  lock = 0;
  crash_version = 0;
  clean = false;
  void* ptr;
  /*FIXME: make the process of initialization crash consistent*/
  Table<T>::New(&ptr, dir->global_depth, nullptr);
  dir->_[initCap - 1] = (Table<T> *)(ptr);
  dir->_[initCap - 1]->pattern = initCap - 1;
  dir->_[initCap - 1]->state = 0;
  /* Initilize the Directory*/
  for (int i = initCap - 2; i >= 0; --i) {
    Table<T>::New(&ptr, dir->global_depth, ptr);
    dir->_[i] = (Table<T> *)(ptr);
    dir->_[i]->pattern = i;
    dir->_[i]->state = 0;
  }
  dir->depth_count = initCap;
}

template <class T>
Finger_EH<T>::Finger_EH() {
  std::cout << "Reinitialize up" << std::endl;
}

template <class T>
Finger_EH<T>::~Finger_EH(void) {
  // TO-DO
}

template <class T>
void Finger_EH<T>::Directory_Doubling(int x, Table<T> *new_b, Table<T> *old_b) {
  Table<T> **d = dir->_;
  auto global_depth = dir->global_depth;
  // std::cout << "Directory_Doubling towards " << global_depth + 1 << std::endl;

  auto capacity = pow(2, global_depth);
  Directory<T>::New(&back_dir, 2 * capacity, dir->version + 1);
  Directory<T> *new_sa =
      reinterpret_cast<Directory<T> *>(back_dir);
  auto dd = new_sa->_;

  for (unsigned i = 0; i < capacity; ++i) {
    dd[2 * i] = d[i];
    dd[2 * i + 1] = d[i];
  }
  dd[2 * x + 1] = reinterpret_cast<Table<T> *>(
      reinterpret_cast<uint64_t>(new_b) | crash_version);
  new_sa->depth_count = 2;

  auto old_dir = dir;

  // Critical Path
  old_b->local_depth += 1;
  dir = new_sa;
}

template <class T>
void Finger_EH<T>::Directory_Update(Directory<T> *_sa, int x, Table<T> *new_b,
                                    Table<T> *old_b) {
  // LOG("Dir update");
  Table<T> **dir_entry = _sa->_;
  auto global_depth = _sa->global_depth;
  unsigned depth_diff = global_depth - new_b->local_depth;
  if (depth_diff == 0) {
    if (x % 2 == 0) {
              dir_entry[x + 1] = reinterpret_cast<Table<T> *>(
                  reinterpret_cast<uint64_t>(new_b) | crash_version);
              old_b->local_depth += 1;
    } else {
              dir_entry[x] = reinterpret_cast<Table<T> *>(
                  reinterpret_cast<uint64_t>(new_b) | crash_version);
              old_b->local_depth += 1;
    }
#ifdef COUNTING
    __sync_fetch_and_add(&_sa->depth_count, 2);
#endif
  } else {
    int chunk_size = pow(2, global_depth - (new_b->local_depth - 1));
    x = x - (x % chunk_size);
    int base = chunk_size / 2;
            for (int i = base - 1; i >= 0; --i) {
              dir_entry[x + base + i] = reinterpret_cast<Table<T> *>(
                  reinterpret_cast<uint64_t>(new_b) | crash_version);
            }
            old_b->local_depth += 1;
  }
  // printf("Done!directory update for %d\n", x);
}

template <class T>
void Finger_EH<T>::Directory_Merge_Update(Directory<T> *_sa, uint64_t key_hash,
                                          Table<T> *left_seg) {
  Table<T> **dir_entry = _sa->_;
  auto global_depth = _sa->global_depth;
  auto x = (key_hash >> (8 * sizeof(key_hash) - global_depth));
  uint64_t chunk_size = pow(2, global_depth - (left_seg->local_depth));
  auto left = x - (x % chunk_size);
  auto right = left + chunk_size / 2;

  for (int i = right; i < right + chunk_size / 2; ++i) {
    dir_entry[i] = left_seg;
  }

  if ((left_seg->local_depth + 1) == global_depth) {
    SUB(&_sa->depth_count, 2);
  }
}

template <class T>
void Finger_EH<T>::recoverTable(Table<T> **target_table, size_t key_hash,
                                size_t x, Directory<T> *old_sa) {
}

template <class T>
void Finger_EH<T>::Recovery() {
}

template <class T>
int Finger_EH<T>::Insert(T key, Value_t value, bool is_in_epoch) {
  return Insert(key, value);
}
template <class T>
int Finger_EH<T>::Insert(T key, Value_t value){
#ifdef AU
  int retval = Insert0(key, value);
  entry_num++;
  uint64_t old_space_num = space_num;
  space_num = old_space_num + segment_num;
  if (space_num < old_space_num){
    LOG("ERROR Length");
  }
  return retval;
}
template <class T>
int Finger_EH<T>::Insert0(T key, Value_t value) {
#endif
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }

  auto meta_hash = ((uint8_t)(key_hash & kMask));  // the last 8 bits
  RETRY:
  auto old_sa = dir;
  auto x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
  auto dir_entry = old_sa->_;
  Table<T> *target = reinterpret_cast<Table<T> *>(
      reinterpret_cast<uint64_t>(dir_entry[x]) & tailMask);

  if ((reinterpret_cast<uint64_t>(dir_entry[x]) & headerMask) !=
      crash_version) {
    recoverTable(&dir_entry[x], key_hash, x, old_sa);
    goto RETRY;
  }

  auto ret = target->Insert(key, value, key_hash, meta_hash, &dir);

  if(ret == -3){ /*duplicate insert, insertion failure*/
    return -1;
  }

  if (ret == -1) {
    if (!target->bucket->try_get_lock()) {
      goto RETRY;
    }

    /*verify procedure*/
    auto old_sa = dir;
    auto x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
    if (reinterpret_cast<Table<T> *>(reinterpret_cast<uint64_t>(old_sa->_[x]) &
                                     tailMask) != target) /* verify process*/
    {
      target->bucket->release_lock();
      goto RETRY;
    }

    auto new_b =
        target->Split(key_hash); /* also needs the verify..., and we use try
                                    lock for this rather than the spin lock*/
    /* update directory*/
    REINSERT:
    old_sa = dir;
    dir_entry = old_sa->_;
    x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
    if (target->local_depth < old_sa->global_depth) {
      if(!try_get_directory_read_lock()){
        goto REINSERT;
      }

      if (old_sa->version != dir->version) {
        // The directory has changed, thus need retry this update
        release_directory_read_lock();
        goto REINSERT;
      }

      Directory_Update(old_sa, x, new_b, target);
      release_directory_read_lock();
    } else {
      Lock_Directory();
      if (old_sa->version != dir->version) {
        Unlock_Directory();
        goto REINSERT;
      }
      Directory_Doubling(x, new_b, target);
      Unlock_Directory();
    }
    /*release the lock for the target bucket and the new bucket*/
    new_b->state = 0;
    target->state = 0;

    Bucket<T> *curr_bucket;
    for (int i = 0; i < kNumBucket; ++i) {
      curr_bucket = target->bucket + i;
      curr_bucket->release_lock();
    }
    curr_bucket = new_b->bucket;
    curr_bucket->release_lock();
    goto RETRY;
  } else if (ret == -2) {
    goto RETRY;
  }
  return 0;
}

template <class T>
Value_t Finger_EH<T>::Get(T key, bool is_in_epoch) {
  if (!is_in_epoch) {
    return Get(key);
  }
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }
  auto meta_hash = ((uint8_t)(key_hash & kMask));  // the last 8 bits
  RETRY:
  auto old_sa = dir;
  auto x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
  auto y = BUCKET_INDEX(key_hash);
  auto dir_entry = old_sa->_;
  auto old_entry = dir_entry[x];
  Table<T> *target = reinterpret_cast<Table<T> *>(
      reinterpret_cast<uint64_t>(old_entry) & tailMask);

  if ((reinterpret_cast<uint64_t>(old_entry) & headerMask) != crash_version) {
    recoverTable(&dir_entry[x], key_hash, x, old_sa);
    goto RETRY;
  }

  Bucket<T> *target_bucket = target->bucket + y;
  Bucket<T> *neighbor_bucket = target->bucket + ((y + 1) & bucketMask);

  uint32_t old_version =
      __atomic_load_n(&target_bucket->version_lock, __ATOMIC_ACQUIRE);
  uint32_t old_neighbor_version =
      __atomic_load_n(&neighbor_bucket->version_lock, __ATOMIC_ACQUIRE);

  if ((old_version & lockSet) || (old_neighbor_version & lockSet)) {
    goto RETRY;
  }

  /*verification procedure*/
  old_sa = dir;
  x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
  if (old_sa->_[x] != old_entry) {
    goto RETRY;
  }

  auto ret = target_bucket->check_and_get(meta_hash, key, false);
  if (target_bucket->test_lock_version_change(old_version)) {
    goto RETRY;
  }
  if (ret != NONE) {
    return ret;
  }

  /*no need for verification procedure, we use the version number of
   * target_bucket to test whether the bucket has ben spliteted*/
  ret = neighbor_bucket->check_and_get(meta_hash, key, true);
  if (neighbor_bucket->test_lock_version_change(old_neighbor_version)) {
    goto RETRY;
  }
  if (ret != NONE) {
    return ret;
  }

  if (target_bucket->test_stash_check()) {
    auto test_stash = false;
    if (target_bucket->test_overflow()) {
      /*this only occur when the bucket has more key-values than 10 that are
       * overfloed int he shared bucket area, therefore it needs to search in
       * the extra bucket*/
      test_stash = true;
    } else {
      /*search in the original bucket*/
      int mask = target_bucket->overflowBitmap & overflowBitmapMask;
      if (mask != 0) {
        for (int i = 0; i < 4; ++i) {
          if (CHECK_BIT(mask, i) &&
              (target_bucket->finger_array[14 + i] == meta_hash) &&
              (((1 << i) & target_bucket->overflowMember) == 0)) {
            Bucket<T> *stash =
                target->bucket + kNumBucket +
                ((target_bucket->overflowIndex >> (i * 2)) & stashMask);
            auto ret = stash->check_and_get(meta_hash, key, false);
            if (ret != NONE) {
              if (target_bucket->test_lock_version_change(old_version)) {
                goto RETRY;
              }
              return ret;
            }
          }
        }
      }

      mask = neighbor_bucket->overflowBitmap & overflowBitmapMask;
      if (mask != 0) {
        for (int i = 0; i < 4; ++i) {
          if (CHECK_BIT(mask, i) &&
              (neighbor_bucket->finger_array[14 + i] == meta_hash) &&
              (((1 << i) & neighbor_bucket->overflowMember) != 0)) {
            Bucket<T> *stash =
                target->bucket + kNumBucket +
                ((neighbor_bucket->overflowIndex >> (i * 2)) & stashMask);
            auto ret = stash->check_and_get(meta_hash, key, false);
            if (ret != NONE) {
              if (target_bucket->test_lock_version_change(old_version)) {
                goto RETRY;
              }
              return ret;
            }
          }
        }
      }
      goto FINAL;
    }
    TEST_STASH:
    if (test_stash == true) {
      for (int i = 0; i < stashBucket; ++i) {
        Bucket<T> *stash =
            target->bucket + kNumBucket + ((i + (y & stashMask)) & stashMask);
        auto ret = stash->check_and_get(meta_hash, key, false);
        if (ret != NONE) {
          if (target_bucket->test_lock_version_change(old_version)) {
            goto RETRY;
          }
          return ret;
        }
      }
    }
  }
  FINAL:
  return NONE;
}

template <class T>
Value_t Finger_EH<T>::Get(T key) {
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }
  auto meta_hash = ((uint8_t)(key_hash & kMask));  // the last 8 bits
  RETRY:
  auto old_sa = dir;
  auto x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
  auto y = BUCKET_INDEX(key_hash);
  auto dir_entry = old_sa->_;
  auto old_entry = dir_entry[x];
  Table<T> *target = reinterpret_cast<Table<T> *>(
      reinterpret_cast<uint64_t>(old_entry) & tailMask);

  if ((reinterpret_cast<uint64_t>(old_entry) & headerMask) != crash_version) {
    recoverTable(&dir_entry[x], key_hash, x, old_sa);
    goto RETRY;
  }

  Bucket<T> *target_bucket = target->bucket + y;
  Bucket<T> *neighbor_bucket = target->bucket + ((y + 1) & bucketMask);

  uint32_t old_version =
      __atomic_load_n(&target_bucket->version_lock, __ATOMIC_ACQUIRE);
  uint32_t old_neighbor_version =
      __atomic_load_n(&neighbor_bucket->version_lock, __ATOMIC_ACQUIRE);

  if ((old_version & lockSet) || (old_neighbor_version & lockSet)) {
    goto RETRY;
  }

  /*verification procedure*/
  old_sa = dir;
  x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
  if (old_sa->_[x] != old_entry) {
    goto RETRY;
  }

  auto ret = target_bucket->check_and_get(meta_hash, key, false);
  if (target_bucket->test_lock_version_change(old_version)) {
    goto RETRY;
  }
  if (ret != NONE) {
    return ret;
  }

  ret = neighbor_bucket->check_and_get(meta_hash, key, true);
  if (neighbor_bucket->test_lock_version_change(old_neighbor_version)) {
    goto RETRY;
  }
  if (ret != NONE) {
    return ret;
  }

  if (target_bucket->test_stash_check()) {
    auto test_stash = false;
    if (target_bucket->test_overflow()) {
      /*this only occur when the bucket has more key-values than 10 that are
       * overfloed int he shared bucket area, therefore it needs to search in
       * the extra bucket*/
      test_stash = true;
    } else {
      /*search in the original bucket*/
      int mask = target_bucket->overflowBitmap & overflowBitmapMask;
      if (mask != 0) {
        for (int i = 0; i < 4; ++i) {
          if (CHECK_BIT(mask, i) &&
              (target_bucket->finger_array[14 + i] == meta_hash) &&
              (((1 << i) & target_bucket->overflowMember) == 0)) {
            Bucket<T> *stash =
                target->bucket + kNumBucket +
                ((target_bucket->overflowIndex >> (i * 2)) & stashMask);
            auto ret = stash->check_and_get(meta_hash, key, false);
            if (ret != NONE) {
              if (target_bucket->test_lock_version_change(old_version)) {
                goto RETRY;
              }
              return ret;
            }
          }
        }
      }

      mask = neighbor_bucket->overflowBitmap & overflowBitmapMask;
      if (mask != 0) {
        for (int i = 0; i < 4; ++i) {
          if (CHECK_BIT(mask, i) &&
              (neighbor_bucket->finger_array[14 + i] == meta_hash) &&
              (((1 << i) & neighbor_bucket->overflowMember) != 0)) {
            Bucket<T> *stash =
                target->bucket + kNumBucket +
                ((neighbor_bucket->overflowIndex >> (i * 2)) & stashMask);
            auto ret = stash->check_and_get(meta_hash, key, false);
            if (ret != NONE) {
              if (target_bucket->test_lock_version_change(old_version)) {
                goto RETRY;
              }
              return ret;
            }
          }
        }
      }
      goto FINAL;
    }
    TEST_STASH:
    if (test_stash == true) {
      for (int i = 0; i < stashBucket; ++i) {
        Bucket<T> *stash =
            target->bucket + kNumBucket + ((i + (y & stashMask)) & stashMask);
        auto ret = stash->check_and_get(meta_hash, key, false);
        if (ret != NONE) {
          if (target_bucket->test_lock_version_change(old_version)) {
            goto RETRY;
          }
          return ret;
        }
      }
    }
  }
  FINAL:
  return NONE;
}

template <class T>
bool Finger_EH<T>::Delete(T key, bool is_in_epoch) {
  return Delete(key);
}

/*By default, the merge operation is disabled*/
template <class T>
bool Finger_EH<T>::Delete(T key) {
  /*Basic delete operation and merge operation*/
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }
  auto meta_hash = ((uint8_t)(key_hash & kMask));  // the last 8 bits
  RETRY:
  auto old_sa = dir;
  auto x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
  auto dir_entry = old_sa->_;
  Table<T> *target_table = reinterpret_cast<Table<T> *>(
      reinterpret_cast<uint64_t>(dir_entry[x]) & tailMask);

  if ((reinterpret_cast<uint64_t>(dir_entry[x]) & headerMask) !=
      crash_version) {
    recoverTable(&dir_entry[x], key_hash, x, old_sa);
    goto RETRY;
  }

  /*we need to first do the locking and then do the verify*/
  auto y = BUCKET_INDEX(key_hash);
  Bucket<T> *target = target_table->bucket + y;
  Bucket<T> *neighbor = target_table->bucket + ((y + 1) & bucketMask);
  target->get_lock();
  if (!neighbor->try_get_lock()) {
    target->release_lock();
    goto RETRY;
  }

  old_sa = dir;
  x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
  if (reinterpret_cast<Table<T> *>(reinterpret_cast<uint64_t>(old_sa->_[x]) &
                                   tailMask) != target_table) {
    target->release_lock();
    neighbor->release_lock();
    goto RETRY;
  }

  auto ret = target->Delete(key, meta_hash, false);
  if (ret == 0) {
#ifdef COUNTING
    auto num = SUB(&target_table->number, 1);
#endif
    target->release_lock();
    neighbor->release_lock();
    return true;
  }

  ret = neighbor->Delete(key, meta_hash, true);
  if (ret == 0) {
#ifdef COUNTING
    auto num = SUB(&target_table->number, 1);
#endif
    neighbor->release_lock();
    target->release_lock();
    return true;
  }

  if (target->test_stash_check()) {
    auto test_stash = false;
    if (target->test_overflow()) {
      /*this only occur when the bucket has more key-values than 10 that are
       * overfloed int he shared bucket area, therefore it needs to search in
       * the extra bucket*/
      test_stash = true;
    } else {
      /*search in the original bucket*/
      int mask = target->overflowBitmap & overflowBitmapMask;
      if (mask != 0) {
        for (int i = 0; i < 4; ++i) {
          if (CHECK_BIT(mask, i) &&
              (target->finger_array[14 + i] == meta_hash) &&
              (((1 << i) & target->overflowMember) == 0)) {
            test_stash = true;
            goto TEST_STASH;
          }
        }
      }

      mask = neighbor->overflowBitmap & overflowBitmapMask;
      if (mask != 0) {
        for (int i = 0; i < 4; ++i) {
          if (CHECK_BIT(mask, i) &&
              (neighbor->finger_array[14 + i] == meta_hash) &&
              (((1 << i) & neighbor->overflowMember) != 0)) {
            test_stash = true;
            break;
          }
        }
      }
    }

    TEST_STASH:
    if (test_stash == true) {
      Bucket<T> *stash = target_table->bucket + kNumBucket;
      stash->get_lock();
      for (int i = 0; i < stashBucket; ++i) {
        int index = ((i + (y & stashMask)) & stashMask);
        Bucket<T> *curr_stash = target_table->bucket + kNumBucket + index;
        auto ret = curr_stash->Delete(key, meta_hash, false);
        if (ret == 0) {
          /*need to unset indicator in original bucket*/
          stash->release_lock();
          auto bucket_ix = BUCKET_INDEX(key_hash);
          auto org_bucket = target_table->bucket + bucket_ix;
          assert(org_bucket == target);
          target->unset_indicator(meta_hash, neighbor, key, index);
#ifdef COUNTING
          auto num = SUB(&target_table->number, 1);
#endif
          neighbor->release_lock();
          target->release_lock();
          return true;
        }
      }
      stash->release_lock();
    }
  }
  neighbor->release_lock();
  target->release_lock();
  return false;
}

/*DEBUG FUNCTION: search the position of the key in this table and print
 * correspongdign informantion in this table, to test whether it is correct*/

template <class T>
int Finger_EH<T>::FindAnyway(T key) {
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    // key_hash = h(key, (reinterpret_cast<string_key *>(key))->length);
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }
  auto meta_hash = ((uint8_t)(key_hash & kMask));
  auto x = (key_hash >> (8 * sizeof(key_hash) - dir->global_depth));

  size_t _count = 0;
  size_t seg_count = 0;
  Directory<T> *seg = dir;
  Table<T> **dir_entry = seg->_;
  Table<T> *ss;
  auto global_depth = seg->global_depth;
  size_t depth_diff;
  int capacity = pow(2, global_depth);
  for (int i = 0; i < capacity;) {
    ss = dir_entry[i];
    Bucket<T> *curr_bucket;
    for (int j = 0; j < kNumBucket; ++j) {
      curr_bucket = ss->bucket + j;
      auto ret = curr_bucket->check_and_get(meta_hash, key, false);
      if (ret != NONE) {
        printf("successfully find in the normal bucket with false\n");
        printf(
            "the segment is %d, the bucket is %d, the local depth = %lld, the "
            "pattern is %lld\n",
            i, j, ss->local_depth, ss->pattern);
        return 0;
      }
      ret = curr_bucket->check_and_get(meta_hash, key, true);
      if (ret != NONE) {
        printf("successfully find in the normal bucket with true\n");
        printf(
            "the segment is %d, the bucket is %d, the local depth is %lld, the "
            "pattern is %lld\n",
            i, j, ss->local_depth, ss->pattern);
        return 0;
      }
    }

    for (int i = 0; i < stashBucket; ++i) {
      curr_bucket = ss->bucket + kNumBucket + i;
      auto ret = curr_bucket->check_and_get(meta_hash, key, false);
      if (ret != NONE) {
        printf("successfully find in the stash bucket\n");
        auto bucket_ix = BUCKET_INDEX(key_hash);
        auto org_bucket = ss->bucket + bucket_ix;
        auto neighbor_bucket = ss->bucket + ((bucket_ix + 1) & bucketMask);
        printf("the segment number is %d, the bucket_ix is %d\n", x, bucket_ix);

        printf("the image of org_bucket\n");
        int mask = org_bucket->overflowBitmap & overflowBitmapMask;
        for (int j = 0; j < 4; ++j) {
          printf(
              "the hash is %d, the pos bit is %d, the alloc bit is %d, the "
              "stash bucket info is %d, the real stash bucket info is %d\n",
              org_bucket->finger_array[14 + j],
              (org_bucket->overflowMember >> (j)) & 1,
              (org_bucket->overflowBitmap >> j) & 1,
              (org_bucket->overflowIndex >> (j * 2)) & stashMask, i);
        }

        printf("the image of the neighbor bucket\n");
        printf("the stash check is %d\n", neighbor_bucket->test_stash_check());
        mask = neighbor_bucket->overflowBitmap & overflowBitmapMask;
        for (int j = 0; j < 4; ++j) {
          printf(
              "the hash is %d, the pos bit is %d, the alloc bit is %d, the "
              "stash bucket info is %d, the real stash bucket info is %d\n",
              neighbor_bucket->finger_array[14 + j],
              (neighbor_bucket->overflowMember >> (j)) & 1,
              (neighbor_bucket->overflowBitmap >> j) & 1,
              (neighbor_bucket->overflowIndex >> (j * 2)) & stashMask, i);
        }

        if (org_bucket->test_overflow()) {
          printf("the org bucket has overflowed\n");
        }
        return 0;
      }
    }

    depth_diff = global_depth - ss->local_depth;
    _count += ss->number;
    seg_count++;
    i += pow(2, depth_diff);
  }
  return -1;
}

#undef BUCKET_INDEX
#undef GET_COUNT
#undef GET_BITMAP
#undef GET_MEMBER
#undef GET_INVERSE_MEMBER
}
