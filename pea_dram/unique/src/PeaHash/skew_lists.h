// Copyright (c) 2021-2023 Institute of Computing Technology, Chinese Academy of Sciences
// Pea Hash is licensed under Mulan PSL v2.

#pragma once

#include <libpmem.h>
#include <pair.h>

#include "space_manager.h"

namespace peas {

#ifdef AU
uint64_t type_cnt[32];
double type_avg[32];
uint64_t pointer_num = 0;
uint64_t space_num = 0;
uint64_t entry_num = 0;
uint64_t segment_num = 0;
uint64_t unique_entry_num = 0;
#endif

constexpr size_t kValSegment = 32 * 2048 - 1;

// in main bucket: segLeaf has 32 bit chunk id + 4 bit deg
struct SegLeaf {
  uint16_t number;
  uint16_t reserved;
  uint32_t next_chunk_id;
  Value_t value[0];

  static inline bool has_next(uint32_t next_c_id) { return (next_c_id != 0); }

  static uint64_t getNumber(SegLeaf *L) {
    uint64_t count = L->number;
    auto SM = SpaceManager::Get();
    uint32_t next_cid = L->next_chunk_id;
    while (next_cid) {
      auto nextL = reinterpret_cast<SegLeaf *>(SM->chunk_addr(next_cid));
      count += nextL->number;
      next_cid = nextL->next_chunk_id;
    }
    return count;
  }

  /**
   * @param ret_array the value array of all values
   * @param number input capacity, output number
   * @return 0: success; -1: not enough
   */
  static int inline lookup(Value_t *array, size_t *number, SegLeaf *L) {
    uint64_t count = L->number;
    auto SM = SpaceManager::Get();
    uint32_t next_cid = L->next_chunk_id;
    int left_cap = *number - count;
    if (left_cap < 0) {
      return -1;
    }
    memcpy(array, L->value, count * sizeof(Value_t));
    array += count;
    while (next_cid) {
      auto nextL = reinterpret_cast<SegLeaf *>(SM->chunk_addr(next_cid));
      BUCKET_PREF(nextL);
      left_cap -= kValSegment;
      if (left_cap < 0) {
        return -1;
      }
      memcpy(array, nextL->value, kValSegment * sizeof(Value_t));
      array += kValSegment;
      next_cid = nextL->next_chunk_id;
    }
    *number = *number - left_cap;
    return 0;
  }

  static int inline del_all(SegLeaf *L, uint32_t chunk_id, uint64_t *ptr) {
    auto SM = SpaceManager::Get();
    uint64_t *freeList = SM->getFreeList(degreeMax);
    SegLeaf *lastL = L;
    uint32_t next_cid = L->next_chunk_id;
    uint64_t last_cid = *ptr >> 32U;
    while (next_cid) {
      lastL = reinterpret_cast<SegLeaf *>(SM->chunk_addr(next_cid));
      last_cid = next_cid;
      next_cid = lastL->next_chunk_id;
    }
    uint64_t final_ptr = (last_cid << 32U) + degreeMax;
    SM->log_ptr_write(*ptr, final_ptr);
    // the last one's next is the free list
    lastL->next_chunk_id = *freeList;
    // free list hdr is L
    *freeList = *ptr;
    SM->log_ptr_clear();
    return 0;
  }

  /**
   *
   * @param v
   * @param L the segment block
   * @return 0: no expansion; > 0: 64bit pointer
   */
  static inline uint64_t insert(Value_t v, uint32_t chunk_id) {
    auto SM = SpaceManager::Get();
    auto L = reinterpret_cast<SegLeaf *>(SM->chunk_addr(chunk_id));
    if (L->number == kValSegment) {
      auto addr_ptr = SM->allocBlock(degreeMax);
      auto newL = reinterpret_cast<SegLeaf *>(addr_ptr.first);
      newL->value[0] = v;
      newL->next_chunk_id = chunk_id;
      newL->number = 1;
      pmem_persist(newL, 64);
#ifdef AU
      segment_num += 1 << degreeMax;
#endif
      return addr_ptr.second;
    } else {
      L->value[L->number] = v;
      pmem_persist(L->value + L->number, sizeof(Value_t));
      L->number++;
      pmem_persist(L, 64);
      return 0;
    }
  }
};

struct value_ptr {
  uint32_t id_num_deg;  // [31:21] 11 bit block id; 20: reserve; [19:4] 16 bit
  // number; [3:0] bit capacity in buckets
  uint32_t chunk_id;

  static inline uint32_t getNumber(value_ptr leaf_ptr) {
    return (leaf_ptr.id_num_deg & number_mask) >> 4;
  }

  static inline uint32_t getDeg(value_ptr leaf_ptr) {
    return leaf_ptr.id_num_deg & deg_mask;
  }

  static inline uint32_t getBlockId(value_ptr leaf_ptr) {
    return (leaf_ptr.id_num_deg & block_id_mask) >> 21;
  }

  static inline uint32_t set_num(uint32_t IdNumDeg, uint32_t new_num) {
    return (IdNumDeg & (~number_mask)) | (new_num << 4);
  }

  static inline uint64_t set_ptr_num(uint64_t ptr, uint64_t new_num) {
    return (ptr & number_64_mask) | (new_num << 4);
  }

  static uint64_t inline getValNum(value_ptr *leaf_ptr) {
    if (getDeg(*leaf_ptr) < degreeMax) {
      return getNumber(*leaf_ptr);
    } else {
      auto SM = SpaceManager::Get();
      auto segLeaf =
          reinterpret_cast<SegLeaf *>(SM->chunk_addr(leaf_ptr->chunk_id));
      return SegLeaf::getNumber(segLeaf);
    }
  }

  static int inline lookup(Value_t *array, size_t *number,
                           value_ptr *leaf_ptr) {
    uint32_t degree = getDeg(*leaf_ptr);
    auto SM = SpaceManager::Get();
    if (degree < degreeMax) {
      uint32_t block_id = getBlockId(*leaf_ptr);
      auto value_blocks = SM->getAddr(leaf_ptr->chunk_id, block_id, degree);
      uint64_t num = getNumber(*leaf_ptr);
      if (num > *number) {
        return -1;
      } else {
        memcpy(array, value_blocks, num * sizeof(Value_t));
        *number = num;
        return 0;
      }
    } else {
      auto segLeaf =
          reinterpret_cast<SegLeaf *>(SM->chunk_addr(leaf_ptr->chunk_id));
      return SegLeaf::lookup(array, number, segLeaf);
    }
  }

  static inline void New(value_ptr *ptr_addr) {
    auto SM = SpaceManager::Get();
    auto pair = SM->allocBlock(0);
    *ptr_addr = *reinterpret_cast<value_ptr *>(&pair.second);
#ifdef AU
    segment_num += 1;
#endif
  }

  /**
   * @param val
   * @return -1: success; >=0: current degree
   */
  static int inline insert(Value_t val, value_ptr *leaf_ptr) {
    value_ptr orig_ptr = *leaf_ptr;
    uint32_t degree = getDeg(orig_ptr);
    auto SM = SpaceManager::Get();
    if (degree < degreeMax) {
      uint32_t number = getNumber(orig_ptr);
      uint32_t cap = 32 << degree;
      if (number == cap) {
        // allocate new
        degree++;
        auto pair = SM->allocBlock(degree);
        uint32_t block_id = getBlockId(orig_ptr);
        auto value_array = reinterpret_cast<Value_t *>(
            SM->getAddr(orig_ptr.chunk_id, block_id, degree));
        if (degree < degreeMax) {
          // double to block
          auto new_val_arr = reinterpret_cast<Value_t *>(pair.first);
          memcpy(new_val_arr, value_array, number * sizeof(Value_t));
          new_val_arr[number] = val;
          pmem_persist(new_val_arr, (number + 1) * sizeof(Value_t));
          uint64_t ptr = set_ptr_num(pair.second, number + 1);
          *leaf_ptr = *reinterpret_cast<value_ptr *>(&ptr);
        } else {
          // double to chunk
          uint64_t *block_hdr = (uint64_t *)pair.first;
          auto new_val_arr = reinterpret_cast<Value_t *>(block_hdr + 1);
          memcpy(new_val_arr, value_array, number * sizeof(Value_t));
          new_val_arr[number] = val;
          number++;
          *block_hdr = number;
          pmem_persist(block_hdr, (number + 1) * sizeof(Value_t));
          uint64_t ptr = pair.second;
          *leaf_ptr = *reinterpret_cast<value_ptr *>(&ptr);
        }
        pmem_persist(leaf_ptr, 8);
        // clear allocating ptr
        SM->log_ptr_clear();
        // todo memory leakage, free val_array
#ifdef AU
        segment_num += 1 << (degree - 1);
#endif
        return 0;
      }
      // normal block insert
      uint32_t block_id = getBlockId(orig_ptr);
      auto value_array = reinterpret_cast<Value_t *>(
          SM->getAddr(orig_ptr.chunk_id, block_id, degree));
      value_array[number] = val;
      // persist
      pmem_persist(value_array + number, sizeof(Value_t));
      leaf_ptr->id_num_deg = set_num(orig_ptr.id_num_deg, number + 1);
      pmem_persist(leaf_ptr, 8);
      return 0;
    } else {
      uint64_t insert_flag = SegLeaf::insert(val, orig_ptr.chunk_id);
      if (insert_flag) {
        // incur another chunk link
        *leaf_ptr = *reinterpret_cast<value_ptr *>(&insert_flag);
        pmem_persist(leaf_ptr, 8);
        SM->log_ptr_clear();
      }
      // normal chunk insert
      return 0;
    }
  }

  static int inline insertBatch(Value_t val, value_ptr *leaf_ptr) {
    value_ptr orig_ptr = *leaf_ptr;
    uint32_t degree = getDeg(orig_ptr);
    auto SM = SpaceManager::Get();
    uint32_t number = getNumber(orig_ptr);
    uint32_t cap = 32 << degree;
    // normal block insert
    uint32_t block_id = getBlockId(orig_ptr);
    auto value_array = reinterpret_cast<Value_t *>(
        SM->getAddr(orig_ptr.chunk_id, block_id, degree));
    value_array[number] = val;
    // persist
    pmem_flush(value_array + number, sizeof(Value_t));
    leaf_ptr->id_num_deg = set_num(orig_ptr.id_num_deg, number + 1);
    return 0;
  }

  // the structure of free list or next ptr are the same as the main hash table
  static int inline del_all(value_ptr *leaf_ptr) {
    value_ptr orig_ptr = *leaf_ptr;
    uint32_t degree = getDeg(orig_ptr);
    auto SM = SpaceManager::Get();
    if (degree < degreeMax) {
      value_ptr *free_hdr =
          reinterpret_cast<value_ptr *>(SM->getFreeList(degree));
      // this.value[0] = free_list
      uint32_t block_id = getBlockId(orig_ptr);
      auto value_array_next = reinterpret_cast<value_ptr *>(
          SM->getAddr(orig_ptr.chunk_id, block_id, degree));
      SM->log_ptr_write(*reinterpret_cast<uint64_t *>(leaf_ptr),
                        *reinterpret_cast<uint64_t *>(leaf_ptr));
      *value_array_next = *free_hdr;
      // free_list = this.ptr
      *free_hdr = *leaf_ptr;
      SM->log_ptr_clear();
    } else {
      auto segLeaf =
          reinterpret_cast<SegLeaf *>(SM->chunk_addr(orig_ptr.chunk_id));
      SegLeaf::del_all(segLeaf, orig_ptr.chunk_id,
                       reinterpret_cast<uint64_t *>(leaf_ptr));
    }
    return 0;
  }
};

#define MAX_LINK_LEN 2000

bool in_map = false;

}  // namespace peas
