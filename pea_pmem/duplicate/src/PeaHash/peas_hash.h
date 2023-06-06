#pragma once

#include <hash.h>
#include <libpmem.h>
#include <pair.h>
#include <pref.h>

#include <queue>
#include <unordered_map>
#include <iomanip>
#include "skew_lists.h"
#include "space_manager.h"
#define PREF_INTERVAL 6

#ifdef ANA
uint64_t skew_meet = 0;
uint64_t skew_found = 0;
uint64_t unique_meet = 0;
uint8_t id_remain = 0;
#endif

namespace peas {

// 31:16 version; 15:8 ref-count; 7: lock bit; 6:0 thread_id;
uint32_t lock_map_table[4096];

// 31: pointer to Pea; 30: lock bit; 29:0 version
uint32_t dir_lock;

const uint32_t lockSet = ((uint32_t)1 << 7U);
const uint32_t dir_lock_set = (1U << 30U);
const uint32_t threadMask = ((uint32_t)1 << 7U) - 1;
constexpr size_t kNumBucket = 2048;
constexpr size_t segment_size = kNumBucket * 256;
constexpr uint16_t lock_table_mask = ((uint16_t)1 << 12U) - 1;
constexpr uint8_t unique_bit = 0;

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
struct _Pair {
  T key;
  Value_t value;
};

struct header_t {     // 16B
  uint16_t r_bitmap;  // [15:1]bitmap; 0 multi
  //  Header: 1 bit multi, 15 bit bitmap, 14B fingerprint
  // 1) unique = 1: all entries are <key,value>
  // 2) unique = 0: k multi -> k+1 bits: 10...0
  // 1 multi, bit_scan_forward, find how many pointers, normal are bitmaps
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
  _Pair<T> entry[16];  // <key, value> or <key, pointer>
  header_t header;
};

// 8B
struct seg_meta {
  uint8_t local_depth;
  uint8_t policy;
  uint16_t reserve;
  uint32_t chunk_id;
};

struct Directory {
  seg_meta _[0];
};

constexpr uint8_t SINGLE_HASH = 0;
constexpr uint8_t DOUBLE_HASH = 1;
constexpr uint8_t POLICY_COUNT = 2;

template <class T>
struct function {
  typedef Bucket<T> bucket_t;
  /**
   * @param hash the hash value of key
   * @param first_bucket the addr of the first bucket of the chunk
   * @param chunk_id the chunk_id in the SM
   * @param dir_lock_word the global dir lock
   * @return -2: concurrency failed; -1: full; 0: Success.
   */
  int (*insert)(T key, Value_t value, uint64_t hash, bucket_t *first_bucket,
                uint32_t chunk_id, uint32_t dir_lock_word);
  /*if success, then return 0, else return -1*/
  int (*del_all)(T key, uint64_t hash, bucket_t *first_bucket,
                 uint32_t chunk_id, uint32_t dir_lock_word);
  int (*del_one)(T key, Value_t value, uint64_t hash, bucket_t *first_bucket,
                 uint32_t chunk_id, uint32_t dir_lock_word);
  /**
   * @param ret_array the value array of all values
   * @param number input capacity, output number
   * @return 0: success; -1: not enough
   */
  int (*lookup)(T key, Value_t *ret_array, size_t *number, uint64_t hash,
                bucket_t *first_bucket, uint32_t chunk_id);
};

static void inline persist_bitmap(header_t *header, uint16_t bitmap) {
  header->r_bitmap = bitmap;
  pmem_flush(header, 64);
  pmem_drain();
}

// support variable length key
template <class T>
uint64_t hash(T key) {
  if constexpr (std::is_pointer_v<T>) {
    return h(key->key, key->length);
  } else {
    return h(&key, sizeof(key));
  }
}

template <class T>
int bucket_compact(Bucket<T> *bucket, T k, Value_t v, uint8_t fgprt);

template <class T>
int bucket_compact_unique(Bucket<T> *bucket, T k, Value_t v, uint8_t fgprt,
                          uint8_t pop_num);

template <class T>
class hashPolicy {
  typedef Bucket<T> bucket_t;
  typedef _Pair<T> pea_entry;

 public:
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
    uint32_t mask = ~((uint32_t)bitmap);
    return _bit_scan_forward((int)mask);
  }

  // from 15 to 1 is OK, if full: 0
  static uint8_t inline last_empty_slot_unique(uint16_t bitmap) {
    return _bit_scan_reverse((int)(((uint16_t)~bitmap) | 1U));
  }

  // from 15 to 2 is OK, if full: ret == pop_num; pop num > 0
  static uint8_t inline last_empty_slot_multi(uint16_t bitmap,
                                              uint8_t pop_num) {
    uint16_t new_bitmap = bitmap & ~(1U << pop_num);
    uint8_t last_empty = _bit_scan_reverse((int)((uint16_t)~new_bitmap));
    return last_empty;
  }

  // ------------------ Bucket Lock -------------------------
  // Spin wait until get lock
  static inline void get_spin_lock(uint16_t lock_id) {
    //    uint8_t thread_id = 1;
    uint32_t *lock = lock_map_table + lock_id;
    uint32_t new_value;
    uint32_t old_value = 0;
    do {
      while (true) {
        old_value = __atomic_load_n(lock, __ATOMIC_ACQUIRE);
        if (!(old_value & lockSet)) {
          new_value = old_value & 0xffff0000UL;
          new_value |= thread_id;
          new_value |= lockSet;
          break;
        }
        if ((old_value & threadMask) == thread_id) {
          new_value = old_value + 0x100UL;
          break;
        }
      }
    } while (!CAS(lock, &old_value, new_value));
  }

  static inline void acquire_lock(uint16_t id) {
    get_spin_lock(id & lock_table_mask);
  }

  static inline void acquire_2locks(uint16_t id0, uint16_t id1) {
    uint16_t lock_id0 = id0 & lock_table_mask;
    uint16_t lock_id1 = id1 & lock_table_mask;
    if (lock_id1 > lock_id0) {
      get_spin_lock(lock_id0);
      get_spin_lock(lock_id1);
    } else {
      get_spin_lock(lock_id1);
      get_spin_lock(lock_id0);
    }
  }

  static inline bool try_acquire_lock(uint16_t lock_id) {
    uint32_t *lock = lock_map_table + (lock_id & lock_table_mask);
    uint32_t new_value;
    uint32_t old_value = __atomic_load_n(lock, __ATOMIC_ACQUIRE);
    if (old_value & lockSet) {
      if ((old_value & threadMask) == thread_id) {
        new_value = old_value + 0x100UL;
        return CAS(lock, &old_value, new_value);
      }
      return false;
    }
    new_value = old_value & 0xffff0000UL;
    new_value |= thread_id;
    new_value |= lockSet;
    return CAS(lock, &old_value, new_value);
  }

  // increment version and release the lock
  static inline void release_lock(uint16_t lock_id) {
    uint32_t *lock = lock_map_table + (lock_id & lock_table_mask);
    uint32_t v = *lock + 0x10000UL;
    v = (v & 0xff00UL) ? v - 0x100UL : v - lockSet;
    __atomic_store_n(lock, v, __ATOMIC_RELEASE);
  }

  static inline bool belong_to_me(uint32_t lock_word) {
    return (lock_word & threadMask) == (uint32_t)thread_id;
  }

  static inline bool is_locked(uint32_t lock_word) {
    if (lock_word & lockSet) {
      if (belong_to_me(lock_word)) {
        return false;
      }
      return true;
    }
    return false;
  }

  static inline void release2locks(uint16_t lock1, uint16_t lock2) {
    release_lock(lock1);
    release_lock(lock2);
  }

  static inline uint32_t get_lock_word(uint16_t lock_id) {
    return __atomic_load_n(lock_map_table + (lock_id & lock_table_mask),
                           __ATOMIC_ACQUIRE);
  }

  static uint64_t inline expand_bits(uint64_t hash, uint8_t new_depth) {
    return hash << (30U + new_depth) >> 62U;
  }

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
    header_t *old_hdr = &old_bucket->header;
    uint16_t old_bitmap = old_hdr->r_bitmap;
    uint8_t k = 1;
    if (!get_bit_true(old_bitmap, unique_bit)) {
      uint8_t pop_num = _bit_scan_forward(old_bitmap);
      for (; k <= pop_num; k++) {
        // rehash multi to bucket
        old_entry = old_bucket->entry + k;
        old_key = old_entry->key;
        key_hash = hash(old_key);
        new_bucket = seg_ptrs[expand_bits(key_hash, new_depth)] + i;
        int ins_ret =
            list_rehash2DRAM(old_key, old_entry->value, key_hash, new_bucket);
        if (ins_ret == -1) LOG("What!");
      }
    }
    for (; k < 16; ++k) {
      if (get_bit_true(old_bitmap, k)) {
        old_entry = old_bucket->entry + k;
        old_key = old_entry->key;
        key_hash = hash(old_key);
        new_bucket = seg_ptrs[expand_bits(key_hash, new_depth)] + i;
        // rehash a single one entry to the bucket
        new_bucket->entry[k] = *old_entry;
        if (k != 15) {
          new_bucket->header.fingerprint[k - 1] = old_bucket->header.fingerprint[k - 1];
        }
        new_bucket->header.r_bitmap = bitmap_set1(new_bucket->header.r_bitmap, k);
      }
    }
  }

  static int list_rehash2DRAM(T key, Value_t value, uint64_t hash,
                              bucket_t *bucket) {
    uint8_t fgprt = hashPolicy<T>::FINGERPRINT(hash);
    header_t *header = &bucket->header;
    uint8_t *fgprt_arr = header->fingerprint;
    uint16_t bitmap = header->r_bitmap;
    uint8_t slot = first_empty_slot(bitmap);

    if (slot) {  // unique
      if (slot == 16) {
        return -1;
        //        return bucket_compactDRAM<T>(bucket, key, value, fgprt, true,
        //        0);
      }
      bucket->entry[slot] = {key, value};
      if (slot != 15) {
        fgprt_arr[slot - 1] = fgprt;
      }
      bitmap = bitmap_set0(bitmap, unique_bit);
      bitmap = bitmap_set1(bitmap, slot);
    } else {  // multi
      uint8_t pop_entry_num = _bit_scan_forward(bitmap);
      uint8_t check = last_empty_slot_multi(bitmap, pop_entry_num);
      if (check == pop_entry_num) {
        return -1;
        //        return bucket_compactDRAM(bucket, key, value, fgprt, true,
        //                                  pop_entry_num);
      }
      bucket->entry[pop_entry_num + 1] = {key, value};
      if (pop_entry_num != 14) {
        fgprt_arr[pop_entry_num] = fgprt;
      }
      bitmap = bitmap_set0(bitmap, pop_entry_num);
      bitmap = bitmap_set1(bitmap, pop_entry_num + 1);
    }
    header->r_bitmap = bitmap;
    return 0;
  }

  /* static int entry_rehash2DRAM(T key, Value_t value, uint64_t hash,
                                bucket_t *bucket) {
     uint8_t fgprt = hashPolicy<T>::FINGERPRINT(hash);
     header_t *header = &bucket->header;
     uint8_t *fgprt_arr = header->fingerprint;
     uint16_t bitmap = header->r_bitmap;
     uint8_t slot;

     if (get_bit_true(bitmap, unique_bit)) {
       slot = last_empty_slot_unique(bitmap);
       if (slot == 0) {
         //        return bucket_compactDRAM<T>(bucket, key, value, fgprt, false,
         //        0);
         return -1;
       }
     } else {  // if have popular entries
       uint8_t pop_entry_num = _bit_scan_forward(bitmap);
       // if unique is 0, search to find if exist
       pea_entry *this_entry;
       if (pop_entry_num == 15) {
         this_entry = bucket->entry + 15;
         if (key == this_entry->key) {
           value_ptr::insertBatch(
               value, reinterpret_cast<value_ptr *>(&(this_entry->value)));
           return 0;
         }
         pop_entry_num--;
       }
       for (int i = 1; i <= pop_entry_num; ++i) {
         if (fgprt_arr[i - 1] == fgprt) {
           this_entry = bucket->entry + i;
           if (this_entry->key == key) {
             value_ptr::insertBatch(
                 value, reinterpret_cast<value_ptr *>(&(this_entry->value)));
             return 0;
           }
         }
       }
       slot = hashPolicy<T>::last_empty_slot_multi(bitmap, pop_entry_num);
       if (slot == pop_entry_num) {
         //        return bucket_compactDRAM<T>(bucket, key, value, fgprt, false,
         //                                     pop_entry_num);
         return -1;
       }
     }

     bucket->entry[slot] = {.key = key, .value = value};
     if (slot != 15) {
       fgprt_arr[slot - 1] = fgprt;
     }
     bitmap = hashPolicy<T>::bitmap_set1(bitmap, slot);
     header->r_bitmap = bitmap;
     return 0;
   }
 */
  static int bucket_insert_NVM(T key, Value_t value, uint64_t hash,
                               bucket_t *bucket) {
    uint8_t fgprt = hashPolicy<T>::FINGERPRINT(hash);
    header_t *header = &bucket->header;

    uint16_t bitmap = header->r_bitmap;
    uint8_t slot = hashPolicy<T>::first_empty_slot(
        bitmap);  // 0 denote multi, > 0 is unique

    if (slot) {
      // unique
      if (slot == 16) {
        int ret = bucket_compact<T>(bucket, key, value, fgprt);
        return ret;
      }
      // enable entry moving
      uint8_t from;
      // write entry
      bucket->entry[slot] = {.key = key, .value = value};

      uint8_t lineid = slot / 4;
      if (lineid == 0) {
        header->fingerprint[slot - 1] = fgprt;
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, slot);
      } else {
        if (slot != 15) {
          header->fingerprint[slot - 1] = fgprt;
        }
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, slot);
        uint8_t *line_empty_slots =
            empty_slots_in_line[(uint8_t)(bitmap >> lineid * 4U) & 0xfU];
        for (from = 1; from <= line_empty_slots[0]; ++from) {
          uint8_t to = line_empty_slots[from];
          // copy value
          bucket->line[lineid].entry[to] = bucket->line[0].entry[from];
          uint8_t to_slotid = lineid * 4 + to;
          if (to_slotid != 15) {
            header->fingerprint[to_slotid - 1] = header->fingerprint[from - 1];
          }
          // validate and invalidate corresponding bitmap
          bitmap = hashPolicy<T>::bitmap_set1(bitmap, to_slotid);
          bitmap = hashPolicy<T>::bitmap_set0(bitmap, from);
        }
        // after entry moving, persist the new address
        pmem_flush(&bucket->line[lineid], 64);
        pmem_drain();
      }
      persist_bitmap(header, bitmap);
    } else {  // if have popular entries
      uint8_t pop_entry_num = _bit_scan_forward(bitmap);
      // if unique is 0, search to find if exist
      uint8_t *fgprt_arr = header->fingerprint;
      pea_entry *this_entry;
      if (pop_entry_num == 15) {
        this_entry = bucket->entry + 15;
        if (key == this_entry->key) {
          value_ptr::insert(value,
                            reinterpret_cast<value_ptr *>(&this_entry->value));
          return 0;
        }
        pop_entry_num--;
      }
      for (int i = 1; i <= pop_entry_num; ++i) {
        if (fgprt_arr[i - 1] == fgprt) {
          this_entry = bucket->entry + i;
          if (this_entry->key == key) {
            value_ptr::insert(
                value, reinterpret_cast<value_ptr *>(&this_entry->value));
            return 0;
          }
        }
      }

      // else insert WO entry moving
      slot = hashPolicy<T>::first_empty_slot(bitmap |
                                             ((1U << pop_entry_num) - 1U));
      if (slot == 16) {
        int ret =
            bucket_compact_unique<T>(bucket, key, value, fgprt, pop_entry_num);
        return ret;
      }
      bucket->entry[slot] = {.key = key, .value = value};
      if (slot != 15) {
        fgprt_arr[slot - 1] = fgprt;
      }
      bitmap = hashPolicy<T>::bitmap_set1(bitmap, slot);
      pmem_flush(bucket->line + slot / 4, 64);
      pmem_drain();
      persist_bitmap(header, bitmap);
    }
    return 0;
  }

  /**
   * @brief search all links
   * @return -1: space not enough; 0: multi search found and success; > 0: not
   * find, return pop entry num
   */
  static inline int bucket_lookup_multi(T key, Value_t *ret_array,
                                        size_t *number, uint8_t fgprt,
                                        bucket_t *bucket) {
    uint8_t j = 1;
    uint16_t bitmap = bucket->header.r_bitmap;
    uint8_t *fgprt_arr = bucket->header.fingerprint;
    pea_entry *this_entry;
    if (!get_bit_true(bitmap, unique_bit)) {
#ifdef ANA
    skew_meet++;
#endif
      // if multi
      uint8_t pop_entry_num = _bit_scan_forward(bitmap);
      if (pop_entry_num == 15) {
        this_entry = bucket->entry + 15;
        if (key == this_entry->key) {
          // if key is skew
          int ret = value_ptr::lookup(
              ret_array, number,
              reinterpret_cast<value_ptr *>(&this_entry->value));
#ifdef ANA
        skew_found++;
#endif
          return ret;
        }
        pop_entry_num--;
      }
      for (; j <= pop_entry_num; ++j) {
        if (fgprt_arr[j - 1] == fgprt) {
          this_entry = bucket->entry + j;
          if (key == this_entry->key) {
            int ret = value_ptr::lookup(
                ret_array, number,
                reinterpret_cast<value_ptr *>(&this_entry->value));
#ifdef ANA
        skew_found++;
#endif            
            return ret;
          }
        }
      }
    }
    return j;
  }

  /**
   * @brief search all entries
   * @return -1: space not enough; 0: search finished
   */
  static inline int bucket_lookup_unique(T key, Value_t *ret_array,
                                         size_t *number, uint8_t fgprt,
                                         bucket_t *bucket, uint8_t j) {
#ifdef ANA
    unique_meet++;
#endif
    uint16_t bitmap = bucket->header.r_bitmap;
    uint8_t *fgprt_arr = bucket->header.fingerprint;
    pea_entry *this_entry;
    size_t cap = *number;
    for (; j < 15; j++) {
      if (fgprt_arr[j - 1] == fgprt && get_bit_true(bitmap, j)) {
        // if they really equal, compare exact keys
        this_entry = bucket->entry + j;
        if (key == this_entry->key) {
          if (cap == 0) return -1;
          cap--;
          // ret_val add this_entry->value;
          *ret_array = this_entry->value;
          ret_array++;
        }
      }
    }
    if (get_bit_true(bitmap, 15)) {
      this_entry = bucket->entry + 15;
      if (key == this_entry->key) {
        if (cap == 0) return -1;
        cap--;
        // ret_val add this_entry->value;
        *ret_array = this_entry->value;
      }
    }
    *number = *number - cap;
    return 0;
  }

  /*if success, then return 0, if not found return -1, if multi delete return
   * -2*/
  //  todo simplify search
  static int bucket_del_all(T key, uint8_t fgprt, bucket_t *bucket,
                            uint32_t dir_lock_word) {
    uint8_t j = 1;
    header_t *header = &bucket->header;
    uint16_t bitmap = header->r_bitmap;
    uint8_t *fgprt_arr = header->fingerprint;
    pea_entry *this_entry;
    std::vector<uint8_t> exist_ids;
    uint8_t pop_num = 0;
    if (!get_bit_true(bitmap, unique_bit)) {
      // if multi
      uint8_t pop_entry_num = _bit_scan_forward(bitmap);
      pop_num = pop_entry_num;
      if (pop_entry_num == 15) {
        this_entry = bucket->entry + 15;
        if (key != this_entry->key) {
          // if key is skew
          exist_ids.push_back(15);
        } else {
          value_ptr::del_all(reinterpret_cast<value_ptr *>(&this_entry->value));
        }
        pop_entry_num--;
      }
      for (; j <= pop_entry_num; ++j) {
        if (fgprt_arr[j - 1] != fgprt) {
          exist_ids.push_back(j);
        } else {
          this_entry = bucket->entry + j;
          if (key != this_entry->key)
            exist_ids.push_back(j);
          else {
            value_ptr::del_all(
                reinterpret_cast<value_ptr *>(&this_entry->value));
          }
        }
      }
    }
    if (exist_ids.size() == pop_num) {
      // normal search
      bool found = false;
      for (; j < 15; j++) {
        if (fgprt_arr[j - 1] == fgprt && get_bit_true(bitmap, j)) {
          // if they really equal, compare exact keys
          this_entry = bucket->entry + j;
          if (key == this_entry->key) {
            bitmap = bitmap_set0(bitmap, j);
            found = true;
          }
        }
      }
      if (get_bit_true(bitmap, 15)) {
        this_entry = bucket->entry + 15;
        if (key == this_entry->key) {
          bitmap = bitmap_set0(bitmap, 15);
          found = true;
        }
      }
      if (found) {
        persist_bitmap(header, bitmap);
        return 0;
      } else {
        return -1;
      }
    } else {
      // construct a new bucket
      Bucket<T> new_bucket;
      uint16_t new_bitmap = 0;
      new_bitmap = bitmap_set1(new_bitmap, exist_ids.size());
      for (int i = 0; i < exist_ids.size(); ++i) {
        new_bucket.entry[1 + i] = bucket->entry[exist_ids[i]];
        uint8_t this_fgprt;
        if (exist_ids[i] == 15) {
          this_fgprt = FINGERPRINT(hash(new_bucket.entry[i + 1].key));
        } else {
          this_fgprt = bucket->header.fingerprint[exist_ids[i] - 1];
        }
        new_bucket.header.fingerprint[i] = this_fgprt;
      }
      for (; j < 15; j++) {
        if (get_bit_true(bitmap, j)) {
          this_entry = bucket->entry + j;
          if (fgprt_arr[j - 1] != fgprt) {
            new_bitmap = bitmap_set1(new_bitmap, j);
            new_bucket.entry[j] = *this_entry;
            new_bucket.header.fingerprint[j - 1] = fgprt_arr[j - 1];
          } else {
            if (key != this_entry->key) {
              new_bitmap = bitmap_set1(new_bitmap, j);
              new_bucket.entry[j] = *this_entry;
              new_bucket.header.fingerprint[j - 1] = fgprt_arr[j - 1];
            }
          }
        }
      }
      if (get_bit_true(bitmap, 15)) {
        this_entry = bucket->entry + 15;
        if (key != this_entry->key) {
          new_bitmap = bitmap_set1(new_bitmap, 15);
          new_bucket.entry[15] = *this_entry;
        }
      }
      new_bucket.header.r_bitmap = new_bitmap;
      auto SM = SpaceManager::Get();
      SM->log_buc_write(&new_bucket);
      *bucket = new_bucket;
      pmem_persist(bucket, 256);
      pmem_drain();
      SM->log_buc_clean();
      return 0;
    }
  }
};

template <class T>
class PeaHashing : public Hash<T> {
 public:
  PeaHashing();
  PeaHashing(size_t);
  ~PeaHashing();
  inline int Insert(T key, Value_t value);
  int Insert(T key, Value_t value, bool);
  inline bool Delete(T key);
  inline bool Delete(T key, Value_t value);
  bool Delete(T key, bool);
  inline int Get(T key, Value_t *ret_array, size_t *number);
  Value_t Get(T) { return nullptr; }
  Value_t Get(T key, bool is_in_epoch) { return nullptr; }
  void Recovery();
  uint64_t getNumber();
  void emphasizeDir();
  void mix_NVM_read(size_t size){
    auto SM = SpaceManager::Get();
    uint64_t seg_num = size;
    uint64_t sum = 0;
    for (size_t i = 0; i < seg_num; i += 2)
    {
      uint64_t * new_NVM_block = reinterpret_cast<uint64_t *>(SM->chunk_addr(i));
      for (size_t j = 0; j < 65536; j++)
      {
        sum += new_NVM_block[j];
      }
    }
    std::cout << "Mix NVM read: " << sum << std::endl;
  }

  void mix_NVM_write(size_t size){
    uint64_t seg_num = size / 524288;
    for (size_t i = 0; i < seg_num; i++)
    {
      uint64_t * new_NVM_block = reinterpret_cast<uint64_t *>(SM->chunk_addr(SM->allocChunk()));
      memset(new_NVM_block, 628, 524288);
    }
    LOG("Mix NVM write");
  }

  int Get_debug(T key, Value_t *ret_array, size_t *number);

  void Directory_Update(uint32_t, seg_meta *, uint8_t);

  void Directory_Doubling(uint32_t, seg_meta *);

 private:
  void init_function_list();

  // g bits in the bigger endian of [31: 32 - global_depth]hash
  // then the expand bit is in the little endian
  static uint32_t inline CHUNK_IDX(uint64_t hash, uint8_t mask_g) {
    return (hash << 32U) >> mask_g;
  }
  void free4chunk(seg_meta *segments) {
    for (int i = 0; i < 4; ++i) {
      SM->freeChunk(segments[i].chunk_id);
    }
  }

  uint8_t *global_depth;

  Directory *dir;

  SpaceManager *SM;

  function<T> class_list[POLICY_COUNT];

  static inline bool try_get_dir_lock() {
    uint32_t old_value = __atomic_load_n(&dir_lock, __ATOMIC_ACQUIRE);
    if (old_value & dir_lock_set) {
      return false;
    }
    uint32_t new_value = old_value | dir_lock_set;
    return CAS(&dir_lock, &old_value, new_value);
  }

  seg_meta *split4(Bucket<T> *this_bucket, uint8_t local_depth,
                   uint32_t chunk_id);
};

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
PeaHashing<T>::PeaHashing(size_t initCap) {
  // @warning: InitCap should be exponential of 2
  if (initCap < 2) {
    LOG_FATAL("[ERROR] The initCap must be larger than 1");
  }

  SM = SpaceManager::Get();

  global_depth = reinterpret_cast<uint8_t *>(SM->meta0);
  dir = reinterpret_cast<Directory *>(SM->dir0);
  uint8_t G = static_cast<size_t>(log2(initCap));
#ifdef AU
  segment_num += initCap * kNumBucket;
#endif
  /* Initialize the Directory*/
  for (int i = initCap - 1; i >= 0; --i) {
    uint32_t chunk_idx = SM->allocChunk();
    auto seg_begin = reinterpret_cast<Bucket<T> *>(SM->chunk_addr(chunk_idx));
    for (uint64_t j = 0; j < kNumBucket; ++j) {
      Bucket<T> *curr_bucket = seg_begin + j;
      curr_bucket->header.r_bitmap = 1U;
      pmem_flush(curr_bucket, 64);
    }
    dir->_[i].local_depth = G;
    dir->_[i].policy = SINGLE_HASH;
    dir->_[i].chunk_id = chunk_idx;
  }
  pmem_flush(dir, initCap * 8);
  pmem_drain();
  *global_depth = G;
  pmem_flush(global_depth, 64);
  pmem_drain();

  dir_lock = 0;
  memset(lock_map_table, 0, 16384);
  init_function_list();
}

template <class T>
PeaHashing<T>::~PeaHashing<T>() {
  std::cout << "Destructor." << std::endl;
}

/**
 * @return always 0: Success.
 */
template <class T>
int PeaHashing<T>::Insert(T key, Value_t value) {
  uint64_t key_hash = hash(key);
  uint32_t dir_lock_word;
  AFTER_EXPAND:
  dir_lock_word = get_dir_lock();
  if (dir_lock_word & dir_lock_set) goto AFTER_EXPAND;
  uint32_t idx = CHUNK_IDX(key_hash, 64 - *global_depth);
  seg_meta target_meta = dir->_[idx];
  uint32_t chunk_idx = target_meta.chunk_id;
  auto seg_begin = reinterpret_cast<Bucket<T> *>(SM->chunk_addr(chunk_idx));
  //  SEG_INSERT:
  // insert to segment
  uint8_t policy = target_meta.policy;
  int ret = class_list[policy].insert(key, value, key_hash, seg_begin,
                                      chunk_idx, dir_lock_word);
  if (ret == -2) goto AFTER_EXPAND;
  if (ret == -1) {
    if (policy == SINGLE_HASH) {
      dir->_[idx].policy = DOUBLE_HASH;
      goto AFTER_EXPAND;
    } else {
      // 1. allocate 2 new chunks, redistribute in it
      // crash: free 2 useless chunks and redo from 1
      uint8_t L = target_meta.local_depth;
      // lock this_seg, forbid all write to the seg
      seg_meta *new_segments = split4(seg_begin, L, chunk_idx);
      // WARNING what if insert in the old chunk
      if (new_segments->local_depth <= *global_depth) {
        // 2.1 Directory updating and release lock
        // crash: according to the insert key which lead to the split,
        //        continue the update process
        Directory_Update(idx, new_segments, L);
        release_dir_lock();
      } else {
        // 2.2. Directory doubling and change hash pointer to the new
        // directory(Atomic Write)
        // crash: redo from 2.2
        Directory_Doubling(idx, new_segments);
        release_flip_lock();
      }
      // 3. free old chunk, free new_segments pair
      // crash: redo from 3
      delete[] new_segments;
      SM->freeChunk(chunk_idx);
#ifdef AU
      segment_num += 3 * kNumBucket;
#endif
      goto AFTER_EXPAND;
    }
  }
#ifdef AU
  entry_num++;
  uint64_t old_space_num = space_num;
  space_num = old_space_num + segment_num * 16;
  if (space_num < old_space_num) LOG("ERROR, LENGTH");
#endif
  return ret;
}

template <class T>
seg_meta *PeaHashing<T>::split4(Bucket<T> *this_bucket, uint8_t local_depth,
                                uint32_t chunk_id) {
  // Allocate 4 segments
  uint32_t i;
  uint8_t j, k;
  uint8_t new_depth = local_depth + 2;

  seg_meta *ret_val = new seg_meta[4];
  Bucket<T> **seg_ptrs = new Bucket<T> *[4];

  for (j = 0; j < 4; j++) {
    uint32_t c_id = SM->allocChunk();
    ret_val[j].local_depth = new_depth;
    ret_val[j].policy = DOUBLE_HASH;
    ret_val[j].chunk_id = c_id;
    seg_ptrs[j] = (Bucket<T> *)malloc(segment_size);
    for (i = 0; i < kNumBucket; ++i) {
      Bucket<T> *cur_bucket = seg_ptrs[j] + i;
      cur_bucket->header.r_bitmap = 1;
    }
  }
  // double_rehash
  // default we dont link lists
  for (i = 0; i < kNumBucket; ++i) {
    Bucket<T> *old_bucket = this_bucket + i;
    hashPolicy<T>::bucket_rehash2segments_naive(old_bucket, new_depth, seg_ptrs,
                                                i);
  }

#ifdef AVX512F
  for (j = 0; j < 4; j++) {
    Bucket<T> *from_seg = seg_ptrs[j];
    auto to_seg =
        reinterpret_cast<Bucket<T> *>(SM->chunk_addr(ret_val[j].chunk_id));
    for (i = 0; i < kNumBucket; ++i) {
      _line<T> *from_line = from_seg[i].line;
      _line<T> *to_line = to_seg[i].line;
      for (k = 0; k < 4; ++k) {
        __m512i new_line0 = _mm512_loadu_si512(from_line + k);
        _mm512_stream_si512((__m512i *)(to_line + k), new_line0);
      }
    }
  }
#endif
  SM->log_ptr_clear();
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
  uint32_t new_chunk_size = old_chunk_size / 4;
  uint32_t other_ids = table00_start_idx + 1;
  for (uint32_t i = other_ids; i < table00_start_idx + new_chunk_size; ++i) {
    dir_entry[i] = new_s[0];
  }
  for (int i = 1; i < 4; ++i) {
    uint32_t start_idx = table00_start_idx + new_chunk_size * i;
    for (uint32_t j = start_idx; j < start_idx + new_chunk_size; ++j) {
      dir_entry[j] = new_s[i];
    }
  }
  pmem_flush(dir_entry + other_ids, (old_chunk_size - 1) * 8);
  pmem_drain();
  dir_entry[table00_start_idx] = new_s[0];
  pmem_flush(dir_entry + table00_start_idx, 64);
  pmem_drain();
}

template <class T>
void PeaHashing<T>::Directory_Doubling(uint32_t idx, seg_meta *new_s) {
  uint8_t old_G = *global_depth;
  // std::cout << "Directory_Doubling towards " << (uint32_t)old_G + 2
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
      new_d->_[new_i + j] = d[i];
    }
  }

  uint32_t idx0 = idx << 2U;
  for (int i = 0; i < 4; ++i) {
    new_d->_[idx0 + i] = new_s[i];
  }
  pmem_flush(new_d, new_capacity * 8);
  pmem_drain();
  *global_depth = old_G + 2;
  pmem_flush(global_depth, 64);
  pmem_drain();
  dir = new_d;
}

// delete all the entries of the key
template <class T>
bool PeaHashing<T>::Delete(T key) {
  uint64_t key_hash = hash(key);
  int ret_val;
  uint32_t dir_lock_word;
  DIR_CHECK:
  dir_lock_word = get_dir_lock();
  if (dir_lock_word & dir_lock_set) goto DIR_CHECK;
  uint32_t idx = CHUNK_IDX(key_hash, 64 - *global_depth);
  seg_meta target_meta = dir->_[idx];
  uint32_t chunk_idx = target_meta.chunk_id;
  auto seg_begin = reinterpret_cast<Bucket<T> *>(SM->chunk_addr(chunk_idx));
  ret_val = class_list[target_meta.policy].del_all(key, key_hash, seg_begin,
                                                   chunk_idx, dir_lock_word);
  // check whether the seg_lock is the same
  if (ret_val == -2) goto DIR_CHECK;
  return (ret_val == 0);
}

// delete the specific <key, value>
template <class T>
bool PeaHashing<T>::Delete(T key, Value_t value) {
  uint64_t key_hash = hash(key);
  int ret_val;
  uint32_t dir_lock_word;
  DIR_CHECK:
  dir_lock_word = get_dir_lock();
  if (dir_lock_word & dir_lock_set) goto DIR_CHECK;
  uint32_t idx = CHUNK_IDX(key_hash, 64 - *global_depth);
  seg_meta target_meta = dir->_[idx];
  uint32_t chunk_idx = target_meta.chunk_id;
  auto seg_begin = reinterpret_cast<Bucket<T> *>(SM->chunk_addr(chunk_idx));
  ret_val = class_list[target_meta.policy].del_one(
      key, value, key_hash, seg_begin, chunk_idx, dir_lock_word);
  // check whether the seg_lock is the same
  if (ret_val == -2) goto DIR_CHECK;
  return (ret_val == 0);
}

/**
 * @param key       one key
 * @param start     the index of value
 * @param ret_val   the empty array of values; return an array of values
 * @param number    the capacity of array; return the number of values in array
 * @return          the number of left values; -2: not found; -1:space not
 * enough
 */
template <class T>
int PeaHashing<T>::Get(T key, Value_t *ret_array, size_t *number) {
  uint64_t key_hash = hash(key);
  int ret_val;
  uint32_t dir_lock_word;
  size_t cap = *number;
  GET_CHECK:
  dir_lock_word = get_dir_lock();
  if (dir_lock_word & dir_lock_set) goto GET_CHECK;
  uint32_t idx = CHUNK_IDX(key_hash, 64 - *global_depth);
  seg_meta target_meta = dir->_[idx];
  uint32_t chunk_idx = target_meta.chunk_id;
  auto seg_begin = reinterpret_cast<Bucket<T> *>(SM->chunk_addr(chunk_idx));
  ret_val = class_list[target_meta.policy].lookup(
      key, ret_array, number, key_hash, seg_begin, chunk_idx);
  if (get_dir_lock() != dir_lock_word) goto GET_CHECK;
  return ret_val;
}

template <class T>
int PeaHashing<T>::Get_debug(T key, Value_t *ret_array, size_t *number){
  uint64_t key_hash = hash(key);
  int ret_val;
  uint32_t dir_lock_word;
  size_t cap = *number;
  GET_CHECK:
  dir_lock_word = get_dir_lock();
  if (dir_lock_word & dir_lock_set) goto GET_CHECK;
  uint32_t idx = CHUNK_IDX(key_hash, 64 - *global_depth);
  seg_meta target_meta = dir->_[idx];
  uint32_t chunk_idx = target_meta.chunk_id;
  auto seg_begin = reinterpret_cast<Bucket<T> *>(SM->chunk_addr(chunk_idx));
  if (target_meta.policy != DOUBLE_HASH){
    LOG("Single Hash occurred!!!");
  }
  ret_val = class_list[target_meta.policy].lookup(
      key, ret_array, number, key_hash, seg_begin, chunk_idx);
  if (get_dir_lock() != dir_lock_word) goto GET_CHECK;
  return ret_val;
}

using namespace std;

template <class T>
struct cmp {
  bool operator()(pair<uint8_t, _Pair<T>> a, pair<uint8_t, _Pair<T>> b) {
    return a.second.key <= b.second.key;
  }
};

/**
 * if original bucket is 0 multi
 * @return -1 overflow, 0 successful insert
 */
template <class T>
int bucket_compact(Bucket<T> *bucket, T k, Value_t v, uint8_t fgprt) {
  _Pair<T> new_pair = {.key = k, .value = v};
  priority_queue<pair<uint8_t, _Pair<T>>, vector<pair<uint8_t, _Pair<T>>>,
      cmp<T>>
      heap;
  _Pair<T> *old_entries = bucket->entry;
  for (int i = 1; i < 16; ++i) {
    heap.push(make_pair(i, old_entries[i]));
  }
  heap.push(make_pair(0, new_pair));

  Bucket<T> new_bucket;
  _Pair<T> *new_entries = new_bucket.entry;
  uint8_t *fgprt_arr = bucket->header.fingerprint;
  uint8_t *new_fgprts = new_bucket.header.fingerprint;
  auto check_key = heap.top();
  heap.pop();
  uint8_t multi_cnt = 0;
  value_ptr *buildingL = nullptr;
  uint8_t unique_ids[16];
  uint8_t unique_cnt = 0;
  uint16_t bitmap = 0;
  auto SM = SpaceManager::Get();
  do {
    auto i = heap.top();
    heap.pop();
    if (i.second.key == check_key.second.key) {
      if (!buildingL) {
        // allocate bucket for check_key
        uint8_t id = check_key.first;
        if (id == 15) id = i.first;
        if (multi_cnt < 14) {
          new_fgprts[multi_cnt] = (id == 0) ? fgprt : fgprt_arr[id - 1];
        }
        multi_cnt++;
        buildingL =
            reinterpret_cast<value_ptr *>(&(new_entries[multi_cnt].value));
        value_ptr::New(buildingL);
        new_entries[multi_cnt].key = check_key.second.key;
        value_ptr::insert(check_key.second.value, buildingL);
      }
      value_ptr::insert(i.second.value, buildingL);
    } else {
      if (buildingL) {
        buildingL = nullptr;
      } else {
        unique_ids[unique_cnt++] = check_key.first;
      }
      check_key = i;
    }
  } while (!heap.empty());

  if (!buildingL) unique_ids[unique_cnt++] = check_key.first;

  // construct bucket
  if (unique_cnt == 16) {
    return -1;
  }
  bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt);
  uint8_t last_id = 16;
  if (unique_cnt == 14) {
    last_id = unique_ids[13];
    if (last_id == 0) {
      new_entries[15] = {k, v};
    } else {
      new_entries[15] = old_entries[last_id];
    }
    unique_cnt--;
    bitmap = hashPolicy<T>::bitmap_set1(bitmap, 15);
  }
  for (int i = 0; i < unique_cnt; ++i) {
    uint8_t old_id = unique_ids[i];
    switch (old_id) {
      case 0:
        new_fgprts[multi_cnt + i] = fgprt;
        new_entries[multi_cnt + i + 1] = {k, v};
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        break;
      case 15:
        if (last_id == 0) {
          new_fgprts[multi_cnt + i] = fgprt;
          new_entries[multi_cnt + i + 1] = new_entries[15];
          bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        } else if (last_id == 16) {
        } else {
          new_fgprts[multi_cnt + i] = fgprt_arr[last_id - 1];
          new_entries[multi_cnt + i + 1] = old_entries[last_id];
          bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        }
        new_entries[15] = old_entries[15];
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, 15);
        break;
      default:
        new_fgprts[multi_cnt + i] = fgprt_arr[old_id - 1];
        new_entries[multi_cnt + i + 1] = old_entries[old_id];
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        break;
    }
  }
  new_bucket.header.r_bitmap = bitmap;
  SM->log_buc_write(&new_bucket);
  // write real NVM
  *bucket = new_bucket;
  pmem_flush(bucket, 256);
  pmem_drain();
  SM->log_buc_clean();
  SM->log_ptr_clear();
  return 0;
}

/**
 * if original bucket is 0 multi
 * @return -1 overflow, 0 successful insert
 */
template <class T>
int bucket_compact_double(Bucket<T> *bucket, Bucket<T> *first_bucket,
                          uint32_t bucket_id, T k, Value_t v, uint8_t fgprt) {
  _Pair<T> new_pair = {.key = k, .value = v};
  priority_queue<pair<uint8_t, _Pair<T>>, vector<pair<uint8_t, _Pair<T>>>,
      cmp<T>>
      heap;
  _Pair<T> *old_entries = bucket->entry;
  for (int i = 1; i < 16; ++i) {
    heap.push(make_pair(i, old_entries[i]));
  }
  heap.push(make_pair(0, new_pair));

  Bucket<T> new_bucket;
  _Pair<T> *new_entries = new_bucket.entry;
  uint8_t *fgprt_arr = bucket->header.fingerprint;
  uint8_t *new_fgprts = new_bucket.header.fingerprint;
  auto check_key = heap.top();
  heap.pop();
  uint8_t multi_cnt = 0;
  value_ptr *buildingL = nullptr;
  uint8_t unique_ids[16];
  uint8_t unique_cnt = 0;
  uint16_t bitmap = 0;
  auto SM = SpaceManager::Get();
  do {
    auto i = heap.top();
    heap.pop();
    if (i.second.key == check_key.second.key) {
      if (!buildingL) {
        // allocate bucket for check_key
        uint8_t id = check_key.first;
        // if id == 15, check key have no fingerprint, find alternative fgprt
        if (id == 15) id = i.first;
        uint8_t this_fgprt = (id == 0) ? fgprt : fgprt_arr[id - 1];
        if (multi_cnt < 14) {
          new_fgprts[multi_cnt] = this_fgprt;
        }
        multi_cnt++;
        buildingL =
            reinterpret_cast<value_ptr *>(&(new_entries[multi_cnt].value));
        value_ptr::New(buildingL);
        T this_key = check_key.second.key;
        new_entries[multi_cnt].key = this_key;
        value_ptr::insert(check_key.second.value, buildingL);
        // get alternative buckets to move entries
        uint64_t this_hash = hash(this_key);
        uint32_t bucket_id0 = hashPolicy<T>::bucket_idx0(this_hash);
        uint32_t bucket_id1 = hashPolicy<T>::bucket_idx1(this_hash);
        if (bucket_id0 == bucket_id1) goto NORMAL_INSERT;
        uint32_t alter_bucket_id =
            (bucket_id == bucket_id1) ? bucket_id0 : bucket_id1;
        // search in unique entries of the alter bucket
        Bucket<T> *alter_bucket = first_bucket + alter_bucket_id;
        SM->log_buc_write(alter_bucket);
        uint16_t alter_bitmap = alter_bucket->header.r_bitmap;
        uint8_t *alter_fgprt_arr = alter_bucket->header.fingerprint;
        uint8_t j = 1;
        if (!hashPolicy<T>::get_bit_true(alter_bitmap, unique_bit)) {
          uint8_t pop_entry_num = _bit_scan_forward(alter_bitmap);
          j = pop_entry_num;
        }
        _Pair<T> *alter_entry;
        uint16_t new_alter_bitmap = alter_bitmap;
        for (; j < 15; j++) {
          if (alter_fgprt_arr[j - 1] == this_fgprt &&
              hashPolicy<T>::get_bit_true(alter_bitmap, j)) {
            // if they really equal, compare exact keys
            alter_entry = alter_bucket->entry + j;
            if (this_key == alter_entry->key) {
              new_alter_bitmap =
                  hashPolicy<T>::bitmap_set0(new_alter_bitmap, j);
              value_ptr::insertBatch(alter_entry->value, buildingL);
            }
          }
        }
        if (hashPolicy<T>::get_bit_true(alter_bitmap, 15)) {
          alter_entry = alter_bucket->entry + 15;
          if (this_key == alter_entry->key) {
            new_alter_bitmap = hashPolicy<T>::bitmap_set0(new_alter_bitmap, 15);
            value_ptr::insertBatch(alter_entry->value, buildingL);
          }
        }
        alter_bucket->header.r_bitmap = new_alter_bitmap;
        pmem_persist(&alter_bucket->header, 64);
      }
      NORMAL_INSERT:
      value_ptr::insert(i.second.value, buildingL);
    } else {
      if (buildingL) {
        SM->log_buc_clean();
        buildingL = nullptr;
      } else {
        unique_ids[unique_cnt++] = check_key.first;
      }
      check_key = i;
    }
  } while (!heap.empty());

  if (!buildingL) unique_ids[unique_cnt++] = check_key.first;

  // construct bucket
  if (unique_cnt == 16) {
    return -1;
  }
  bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt);
  uint8_t last_id = 16;
  if (unique_cnt == 14) {
    last_id = unique_ids[13];
    if (last_id == 0) {
      new_entries[15] = {k, v};
    } else {
      new_entries[15] = old_entries[last_id];
    }
    unique_cnt--;
    bitmap = hashPolicy<T>::bitmap_set1(bitmap, 15);
  }
  for (int i = 0; i < unique_cnt; ++i) {
    uint8_t old_id = unique_ids[i];
    switch (old_id) {
      case 0:
        new_fgprts[multi_cnt + i] = fgprt;
        new_entries[multi_cnt + i + 1] = {k, v};
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        break;
      case 15:
        if (last_id == 0) {
          new_fgprts[multi_cnt + i] = fgprt;
          new_entries[multi_cnt + i + 1] = new_entries[15];
          bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        } else if (last_id == 16) {
        } else {
          new_fgprts[multi_cnt + i] = fgprt_arr[last_id - 1];
          new_entries[multi_cnt + i + 1] = old_entries[last_id];
          bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        }
        new_entries[15] = old_entries[15];
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, 15);
        break;
      default:
        new_fgprts[multi_cnt + i] = fgprt_arr[old_id - 1];
        new_entries[multi_cnt + i + 1] = old_entries[old_id];
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        break;
    }
  }
  new_bucket.header.r_bitmap = bitmap;
  SM->log_buc_write(&new_bucket);
  // write real NVM
  *bucket = new_bucket;
  pmem_flush(bucket, 256);
  pmem_drain();
  SM->log_buc_clean();
  SM->log_ptr_clear();
  return 0;
}

/**
 * if original bucket has some multi-keys
 * @return -1 overflow, 0 successful insert
 */
template <class T>
int bucket_compact_unique(Bucket<T> *bucket, T k, Value_t v, uint8_t fgprt,
                          uint8_t pop_num) {
  if (pop_num == 15){
    return -1;
  }
  _Pair<T> new_pair = {.key = k, .value = v};
  priority_queue<pair<uint8_t, _Pair<T>>, vector<pair<uint8_t, _Pair<T>>>,
      cmp<T>>
      heap;
  _Pair<T> *old_entries = bucket->entry;
  for (int i = pop_num + 1; i < 16; ++i) {
    heap.push(make_pair(i, old_entries[i]));
  }
  heap.push(make_pair(0, new_pair));

  Bucket<T> new_bucket;
  _Pair<T> *new_entries = new_bucket.entry;
  uint8_t *fgprt_arr = bucket->header.fingerprint;
  uint8_t *new_fgprts = new_bucket.header.fingerprint;
  auto check_key = heap.top();
  heap.pop();
  uint8_t multi_cnt = pop_num;
  value_ptr *buildingL = nullptr;
  uint8_t unique_ids[16];
  uint8_t unique_cnt = 0;
  uint16_t bitmap = 0;
  auto SM = SpaceManager::Get();
  do {
    auto i = heap.top();
    heap.pop();
    if (i.second.key == check_key.second.key) {
      if (!buildingL) {
        // allocate bucket for check_key
        uint8_t id = check_key.first;
        if (id == 15) id = i.first;
        if (multi_cnt < 14) {
          new_fgprts[multi_cnt] = (id == 0) ? fgprt : fgprt_arr[id - 1];
        }
        multi_cnt++;
        buildingL =
            reinterpret_cast<value_ptr *>(&(new_entries[multi_cnt].value));
        value_ptr::New(buildingL);
        new_entries[multi_cnt].key = check_key.second.key;
        value_ptr::insert(check_key.second.value, buildingL);
      }
      value_ptr::insert(i.second.value, buildingL);
    } else {
      if (buildingL) {
        buildingL = nullptr;
      } else {
        unique_ids[unique_cnt++] = check_key.first;
      }
      check_key = i;
    }
  } while (!heap.empty());

  if (!buildingL) unique_ids[unique_cnt++] = check_key.first;

  // construct bucket
  if (unique_cnt == 16 - pop_num) {
    return -1;
  }
  memcpy(new_entries + 1, old_entries + 1, pop_num * 16);
  memcpy(new_fgprts, fgprt_arr, pop_num);
  bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt);
  uint8_t last_id = 16;
  if (unique_cnt == 14 - pop_num && unique_cnt > 0) {
    last_id = unique_ids[unique_cnt - 1];
    if (last_id == 0) {
      new_entries[15] = {k, v};
    } else {
      new_entries[15] = old_entries[last_id];
    }
    unique_cnt--;
    bitmap = hashPolicy<T>::bitmap_set1(bitmap, 15);
  }
  for (int i = 0; i < unique_cnt; ++i) {
    uint8_t old_id = unique_ids[i];
    switch (old_id) {
      case 0:
        new_fgprts[multi_cnt + i] = fgprt;
        new_entries[multi_cnt + i + 1] = {k, v};
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        break;
      case 15:
        if (last_id == 0) {
          new_fgprts[multi_cnt + i] = fgprt;
          new_entries[multi_cnt + i + 1] = new_entries[15];
          bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        } else if (last_id == 16) {
        } else {
          new_fgprts[multi_cnt + i] = fgprt_arr[last_id - 1];
          new_entries[multi_cnt + i + 1] = old_entries[last_id];
          bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        }
        new_entries[15] = old_entries[15];
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, 15);
        break;
      default:
        new_fgprts[multi_cnt + i] = fgprt_arr[old_id - 1];
        new_entries[multi_cnt + i + 1] = old_entries[old_id];
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        break;
    }
  }
  new_bucket.header.r_bitmap = bitmap;
  SM->log_buc_write(&new_bucket);
  // write real NVM
  *bucket = new_bucket;
  pmem_flush(bucket, 256);
  pmem_drain();
  SM->log_buc_clean();
  SM->log_ptr_clear();
  return 0;
}

/**
 * if original bucket has some multi-keys
 * @return -1 overflow, 0 successful insert
 */
template <class T>
int bucket_compact_unique_double(Bucket<T> *bucket, Bucket<T> *first_bucket,
                                 uint32_t bucket_id, T k, Value_t v,
                                 uint8_t fgprt, uint8_t pop_num) {
  if (pop_num == 15){
    return -1;
  }
  _Pair<T> new_pair = {.key = k, .value = v};
  priority_queue<pair<uint8_t, _Pair<T>>, vector<pair<uint8_t, _Pair<T>>>,
      cmp<T>>
      heap;
  _Pair<T> *old_entries = bucket->entry;
  for (int i = pop_num + 1; i < 16; ++i) {
    heap.push(make_pair(i, old_entries[i]));
  }
  heap.push(make_pair(0, new_pair));

  Bucket<T> new_bucket;
  _Pair<T> *new_entries = new_bucket.entry;
  uint8_t *fgprt_arr = bucket->header.fingerprint;
  uint8_t *new_fgprts = new_bucket.header.fingerprint;
  auto check_key = heap.top();
  heap.pop();
  uint8_t multi_cnt = pop_num;
  value_ptr *buildingL = nullptr;
  uint8_t unique_ids[16];
  uint8_t unique_cnt = 0;
  uint16_t bitmap = 0;
  auto SM = SpaceManager::Get();
  do {
    auto i = heap.top();
    heap.pop();
    if (i.second.key == check_key.second.key) {
      if (!buildingL) {
        // allocate bucket for check_key
        uint8_t id = check_key.first;
        // if id == 15, check key have no fingerprint, find alternative fgprt
        if (id == 15) id = i.first;
        uint8_t this_fgprt = (id == 0) ? fgprt : fgprt_arr[id - 1];
        if (multi_cnt < 14) {
          new_fgprts[multi_cnt] = this_fgprt;
        }
        multi_cnt++;
        buildingL =
            reinterpret_cast<value_ptr *>(&(new_entries[multi_cnt].value));
        value_ptr::New(buildingL);
        T this_key = check_key.second.key;
        new_entries[multi_cnt].key = this_key;
        value_ptr::insert(check_key.second.value, buildingL);
        // get alternative buckets to move entries
        uint64_t this_hash = hash(this_key);
        uint32_t bucket_id0 = hashPolicy<T>::bucket_idx0(this_hash);
        uint32_t bucket_id1 = hashPolicy<T>::bucket_idx1(this_hash);
        if (bucket_id0 == bucket_id1) goto NORMAL_INSERT;
        uint32_t alter_bucket_id =
            (bucket_id == bucket_id1) ? bucket_id0 : bucket_id1;
        // search in unique entries of the alter bucket
        Bucket<T> *alter_bucket = first_bucket + alter_bucket_id;
        SM->log_buc_write(alter_bucket);
        uint16_t alter_bitmap = alter_bucket->header.r_bitmap;
        uint8_t *alter_fgprt_arr = alter_bucket->header.fingerprint;
        uint8_t j = 1;
        if (!hashPolicy<T>::get_bit_true(alter_bitmap, unique_bit)) {
          uint8_t pop_entry_num = _bit_scan_forward(alter_bitmap);
          j = pop_entry_num;
        }
        _Pair<T> *alter_entry;
        uint16_t new_alter_bitmap = alter_bitmap;
        for (; j < 15; j++) {
          if (alter_fgprt_arr[j - 1] == this_fgprt &&
              hashPolicy<T>::get_bit_true(alter_bitmap, j)) {
            // if they really equal, compare exact keys
            alter_entry = alter_bucket->entry + j;
            if (this_key == alter_entry->key) {
              new_alter_bitmap =
                  hashPolicy<T>::bitmap_set0(new_alter_bitmap, j);
              value_ptr::insertBatch(alter_entry->value, buildingL);
            }
          }
        }
        if (hashPolicy<T>::get_bit_true(alter_bitmap, 15)) {
          alter_entry = alter_bucket->entry + 15;
          if (this_key == alter_entry->key) {
            new_alter_bitmap = hashPolicy<T>::bitmap_set0(new_alter_bitmap, 15);
            value_ptr::insertBatch(alter_entry->value, buildingL);
          }
        }
        alter_bucket->header.r_bitmap = new_alter_bitmap;
        pmem_persist(&alter_bucket->header, 64);
      }
      NORMAL_INSERT:
      value_ptr::insert(i.second.value, buildingL);
    } else {
      if (buildingL) {
        SM->log_buc_clean();
        buildingL = nullptr;
      } else {
        unique_ids[unique_cnt++] = check_key.first;
      }
      check_key = i;
    }
  } while (!heap.empty());

  if (!buildingL) unique_ids[unique_cnt++] = check_key.first;

  // construct bucket
  if (unique_cnt == 16 - pop_num) {
    return -1;
  }
  memcpy(new_entries + 1, old_entries + 1, pop_num * 16);
  memcpy(new_fgprts, fgprt_arr, pop_num);
  bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt);
  uint8_t last_id = 16;
  if (unique_cnt == 14 - pop_num && unique_cnt > 0) {
    last_id = unique_ids[unique_cnt - 1];
    if (last_id == 0) {
      new_entries[15] = {k, v};
    } else {
      new_entries[15] = old_entries[last_id];
    }
    unique_cnt--;
    bitmap = hashPolicy<T>::bitmap_set1(bitmap, 15);
  }
  for (int i = 0; i < unique_cnt; ++i) {
    uint8_t old_id = unique_ids[i];
    switch (old_id) {
      case 0:
        new_fgprts[multi_cnt + i] = fgprt;
        new_entries[multi_cnt + i + 1] = {k, v};
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        break;
      case 15:
        if (last_id == 0) {
          new_fgprts[multi_cnt + i] = fgprt;
          new_entries[multi_cnt + i + 1] = new_entries[15];
          bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        } else if (last_id == 16) {
        } else {
          new_fgprts[multi_cnt + i] = fgprt_arr[last_id - 1];
          new_entries[multi_cnt + i + 1] = old_entries[last_id];
          bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        }
        new_entries[15] = old_entries[15];
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, 15);
        break;
      default:
        new_fgprts[multi_cnt + i] = fgprt_arr[old_id - 1];
        new_entries[multi_cnt + i + 1] = old_entries[old_id];
        bitmap = hashPolicy<T>::bitmap_set1(bitmap, multi_cnt + i + 1);
        break;
    }
  }
  new_bucket.header.r_bitmap = bitmap;
  SM->log_buc_write(&new_bucket);
  // write real NVM
  *bucket = new_bucket;
  pmem_flush(bucket, 256);
  pmem_drain();
  SM->log_buc_clean();
  SM->log_ptr_clear();
  return 0;
}

template <class T>
class singleHash : public hashPolicy<T> {
 public:
  typedef Bucket<T> bucket_t;
  typedef _Pair<T> pea_entry;
  singleHash() = default;
  ~singleHash() = default;

  static int lookup(T key, Value_t *ret_array, size_t *number, uint64_t hash,
                    bucket_t *first_bucket, uint32_t chunk_id) {
    uint32_t bucket_id = hashPolicy<T>::bucket_idx0(hash);
    bucket_t *bucket = &first_bucket[bucket_id];
    LINE_PREF(bucket);
    int pop_num;
    uint8_t fgprt = hashPolicy<T>::FINGERPRINT(hash);
    pop_num = hashPolicy<T>::bucket_lookup_multi(key, ret_array, number, fgprt,
                                                 bucket);
    if (pop_num <= 0) return pop_num;
    return hashPolicy<T>::bucket_lookup_unique(key, ret_array, number, fgprt,
                                               bucket, pop_num);
  }

  /**
   * @return -2: concurrency failed; -1: full; 0: Success.
   */
  static int insert(T key, Value_t value, uint64_t hash, bucket_t *first_bucket,
                    uint32_t chunk_id, uint32_t dir_lock_word) {
    bucket_t *bucket;
    uint32_t bucket_id = hashPolicy<T>::bucket_idx0(hash);
    bucket = &first_bucket[bucket_id];
    LINE_PREF(bucket);
    return hashPolicy<T>::bucket_insert_NVM(key, value, hash, bucket);
  }

  /*if success, then return 0, else return -1*/
  static int del_all(T key, uint64_t hash_val, bucket_t *first_bucket,
                     uint32_t chunk_id, uint32_t dir_lock_word) {
    uint32_t bucket_id = hashPolicy<T>::bucket_idx0(hash_val);
    bucket_t *bucket = &first_bucket[bucket_id];
    uint8_t fgprt = hashPolicy<T>::FINGERPRINT(hash_val);
    LINE_PREF(bucket);
    return hashPolicy<T>::bucket_del_all(key, fgprt, bucket, dir_lock_word);
  }

  static int del_one(T key, Value_t value, uint64_t hash,
                     bucket_t *first_bucket, uint32_t chunk_id,
                     uint32_t dir_lock_word) {
    return 0;
  }
};

template <class T>
class doubleHash : public hashPolicy<T> {
 public:
  typedef Bucket<T> bucket_t;
  typedef _Pair<T> pea_entry;

  doubleHash() = default;

  ~doubleHash() = default;

  static int lookup(T key, Value_t *ret_array, size_t *number, uint64_t hash,
                    bucket_t *first_bucket, uint32_t chunk_id) {
    // Calculate 2 indexes of buckets
    uint32_t bucket_id1 = hashPolicy<T>::bucket_idx1(hash);
    bucket_t *bucket1 = &first_bucket[bucket_id1];
    LINE_PREF(bucket1);
    uint32_t bucket_id0 = hashPolicy<T>::bucket_idx0(hash);
    int pop_num1, pop_num0;
    uint8_t fgprt = hashPolicy<T>::FINGERPRINT(hash);
    if (bucket_id0 == bucket_id1) {
      pop_num1 = hashPolicy<T>::bucket_lookup_multi(key, ret_array, number, fgprt,
                                                    bucket1);
      if (pop_num1 <= 0) return pop_num1;
      return hashPolicy<T>::bucket_lookup_unique(key, ret_array, number, fgprt,
                                                 bucket1, pop_num1);
    }
    bucket_t *bucket0 = &first_bucket[bucket_id0];
    LINE_PREF(bucket0);
    size_t cap = *number;
    size_t cnt = 0;
    pop_num1 = hashPolicy<T>::bucket_lookup_multi(key, ret_array, number, fgprt,
                                                  bucket1);
    if (pop_num1 <= 0) return pop_num1;
    pop_num0 = hashPolicy<T>::bucket_lookup_multi(key, ret_array, number, fgprt,
                                                  bucket0);
    if (pop_num0 <= 0) return pop_num0;
    if (hashPolicy<T>::bucket_lookup_unique(key, ret_array, &cap, fgprt,
                                            bucket1, pop_num1))
      return -1;
    cnt += cap;
    ret_array += cap;
    cap = *number - cap;
    if (hashPolicy<T>::bucket_lookup_unique(key, ret_array, &cap, fgprt,
                                            bucket0, pop_num0))
      return -1;
    *number = cnt + cap;
    return 0;
  }

  /**
   * @brief only insert into the multi entries
   */
  static inline int insert_in_multi(uint16_t bitmap, bucket_t *bucket, T key,
                                    Value_t value, uint8_t fgprt) {
    header_t *header = &bucket->header;
    uint8_t pop_entry_num = _bit_scan_forward(bitmap);
    // if unique is 0, search to find if exist
    uint8_t *fgprt_arr = header->fingerprint;
    pea_entry *this_entry;
    if (pop_entry_num == 15) {
      this_entry = bucket->entry + 15;
      if (key == this_entry->key) {
        value_ptr::insert(value,
                          reinterpret_cast<value_ptr *>(&this_entry->value));
        return 0;
      }
      pop_entry_num--;
    }
    for (int i = 1; i <= pop_entry_num; ++i) {
      if (fgprt_arr[i - 1] == fgprt) {
        this_entry = bucket->entry + i;
        if (this_entry->key == key) {
          value_ptr::insert(value,
                            reinterpret_cast<value_ptr *>(&this_entry->value));
          return 0;
        }
      }
    }
    return -1;
  }

  /**
   * @return -2: concurrency failed; -1: full; 0:
   * Success.
   */
  static int insert(T key, Value_t value, uint64_t hash, bucket_t *first_bucket,
                    uint32_t chunk_id, uint32_t dir_lock_word) {
    bucket_t *bucket0, *bucket1;
    // prefetch 8 cache line of 2 candidate buckets
    uint32_t bucket_id0 = hashPolicy<T>::bucket_idx0(hash);
    uint32_t bucket_id1 = hashPolicy<T>::bucket_idx1(hash);
    bucket0 = &first_bucket[bucket_id0];
    LINE_PREF(bucket0);
    bucket1 = &first_bucket[bucket_id1];
    LINE_PREF(bucket1);
    uint8_t fgprt = hashPolicy<T>::FINGERPRINT(hash);
    uint8_t from;

    uint16_t bitmap0 = bucket0->header.r_bitmap;
    uint16_t bitmap1 = bucket1->header.r_bitmap;
    uint8_t pseudo_slot0 = hashPolicy<T>::first_empty_slot(bitmap0);
    uint8_t pseudo_slot1;
    int ret = -1;
    if (!pseudo_slot0) {
      ret = insert_in_multi(bitmap0, bucket0, key, value, fgprt);
    }
    // if one of 2 bucket has multi, check if the key belongs to the multi key
    // if it is, insert into map
    if (ret == 0) {
      return ret;
    } else {
      pseudo_slot1 = hashPolicy<T>::first_empty_slot(bitmap1);
      if (!pseudo_slot1) {
        ret = insert_in_multi(bitmap1, bucket1, key, value, fgprt);
        if (ret == 0) return ret;
      }
    }
    uint8_t pop_num0 = _bit_scan_forward(bitmap0);
    uint8_t pop_num1 = _bit_scan_forward(bitmap1);
    uint8_t slot0 =
        (pseudo_slot0 == 0)
        ? hashPolicy<T>::first_empty_slot(bitmap0 | ((1U << pop_num0) - 1U))
        : pseudo_slot0;
    uint8_t slot1 =
        (pseudo_slot1 == 0)
        ? hashPolicy<T>::first_empty_slot(bitmap1 | ((1U << pop_num1) - 1U))
        : pseudo_slot1;
    uint16_t bitmap;
    bucket_t *bucket;
    uint8_t slot;
    // if not or no multi, choose the less occupied to normally insert
    if (slot0 < slot1) {
      bitmap = bitmap0;
      bucket = bucket0;
      slot = slot0;
    } else {
      bitmap = bitmap1;
      bucket = bucket1;
      slot = slot1;
    }
    if (slot == 16) {
      // if fail, check the bucket to compact
      if (pseudo_slot0 == 0)
        ret = bucket_compact_unique_double(bucket0, first_bucket, bucket_id0,
                                           key, value, fgprt, pop_num0);
      else
        ret = bucket_compact_double(bucket0, first_bucket, bucket_id0, key,
                                    value, fgprt);
      if (ret == 0) {
        return ret;
      }
      if (pseudo_slot1 == 0)
        ret = bucket_compact_unique_double(bucket1, first_bucket, bucket_id1,
                                           key, value, fgprt, pop_num1);
      else
        ret = bucket_compact_double(bucket1, first_bucket, bucket_id1, key,
                                    value, fgprt);
      return ret;
    }
    bucket->entry[slot] = {.key = key, .value = value};
    if (slot != 15) {
      bucket->header.fingerprint[slot - 1] = fgprt;
    }
    bitmap = hashPolicy<T>::bitmap_set1(bitmap, slot);
    pmem_flush(bucket->line + slot / 4, 64);
    pmem_drain();
    persist_bitmap(&bucket->header, bitmap);
    return 0;
  }

  /*if success, then return 0, else return -1*/
  static int del_all(T key, uint64_t hash_val, bucket_t *first_bucket,
                     uint32_t chunk_id, uint32_t dir_lock_word) {
    // Calculate 2 indexes of buckets
    uint8_t fgprt = hashPolicy<T>::FINGERPRINT(hash_val);
    uint32_t bucket_id1 = hashPolicy<T>::bucket_idx1(hash_val);
    bucket_t *bucket1 = &first_bucket[bucket_id1];
    LINE_PREF(bucket1);
    uint32_t bucket_id0 = hashPolicy<T>::bucket_idx0(hash_val);
    bucket_t *bucket0 = &first_bucket[bucket_id0];
    LINE_PREF(bucket0);

    hashPolicy<T>::bucket_del_all(key, fgprt, bucket1, dir_lock_word);
    hashPolicy<T>::bucket_del_all(key, fgprt, bucket0, dir_lock_word);
    return 0;
  }

  static int del_one(T key, Value_t value, uint64_t hash,
                     bucket_t *first_bucket, uint32_t chunk_id,
                     uint32_t dir_lock_word) {
    return 0;
  }
};

template <class T>
void PeaHashing<T>::init_function_list() {
  class_list[SINGLE_HASH].insert = singleHash<T>::insert;
  class_list[SINGLE_HASH].del_all = singleHash<T>::del_all;
  class_list[SINGLE_HASH].del_one = singleHash<T>::del_one;
  class_list[SINGLE_HASH].lookup = singleHash<T>::lookup;
  class_list[DOUBLE_HASH].insert = doubleHash<T>::insert;
  class_list[DOUBLE_HASH].del_all = doubleHash<T>::del_all;
  class_list[DOUBLE_HASH].del_one = doubleHash<T>::del_one;
  class_list[DOUBLE_HASH].lookup = doubleHash<T>::lookup;
}

// TODO
template <class T>
void PeaHashing<T>::Recovery() {
  SM = SpaceManager::Get();
  auto G0 = reinterpret_cast<uint8_t *>(SM->meta0);
  auto G1 = reinterpret_cast<uint8_t *>(SM->meta1);
  if (*G0 < *G1) {
    global_depth = G1;
    dir = reinterpret_cast<Directory *>(SM->dir1);
    dir_lock = 0x80000000;
  } else {
    global_depth = G0;
    dir = reinterpret_cast<Directory *>(SM->dir0);
    dir_lock = 0;
  }
  uint32_t g = *global_depth;
  seg_meta *dir_entry = dir->_;
  uint32_t capacity = 1U << g;
  uint8_t depth_diff;
  seg_meta s{};
  for (int i = 0; i < capacity;) {
    s = dir_entry[i];
    // modify bitmap
    SM->bitmap_set(s.chunk_id);
    depth_diff = g - s.local_depth;
    int j;
    for (j = i + 1; j < i + (1U << depth_diff); ++j) {
      dir_entry[j] = s;
    }
    i = j;
  }
  pmem_flush(dir, capacity * 8);
  pmem_drain();
  memset(lock_map_table, 0, 16384);
  init_function_list();
}

void check_all_lock_released() {
  for (unsigned int i : lock_map_table) {
    if (i & lockSet) LOG("The concurrency control is problem!");
  }
}

template <class T>
uint64_t PeaHashing<T>::getNumber() {
#ifdef AU
  for (int i = 0; i < 32; ++i) {
    type_cnt[i] = 0;
    type_avg[i] = 0;
  }
#endif
  assert(SpaceManager::Get()->log->ptr_log_cur == 0);
  seg_meta *dir_entry = dir->_;
  uint8_t g = *global_depth;
  uint8_t depth_diff;
  uint64_t count = 0;
  uint64_t seg_count = 0;
  seg_meta s{};
  int capacity = pow(2, g);
  for (int i = 0; i < capacity;) {
    s = dir_entry[i];
    depth_diff = g - s.local_depth;
    auto first_bucket =
        reinterpret_cast<Bucket<T> *>(SM->chunk_addr(s.chunk_id));
    for (uint64_t j = 0; j < kNumBucket; ++j) {
      Bucket<T> *check_bucket = first_bucket + j;
      header_t *check_header = &check_bucket->header;
      uint16_t bitmap = check_header->r_bitmap;
      // multi
      uint8_t k = 1;
      if (!hashPolicy<T>::get_bit_true(bitmap, unique_bit)) {
        uint8_t pop_entry_num = _bit_scan_forward(bitmap);
        for (k = 1; k <= pop_entry_num; ++k) {
          _Pair<T> *this_entry = check_bucket->entry + k;
#ifdef AU
          uint64_t arr_cnt = value_ptr::getValNum(
              reinterpret_cast<value_ptr *>(&(this_entry->value)));
          uint8_t index = _bit_scan_reverse(arr_cnt);
          index = (index == 0)? 0: _bit_scan_reverse(arr_cnt - 1);
          type_avg[index] = (arr_cnt + type_avg[index] * type_cnt[index]) /
                            (type_cnt[index] + 1.0);
          type_cnt[index]++;
          count += arr_cnt;
          pointer_num ++;
#else
          count += value_ptr::getValNum(
              reinterpret_cast<value_ptr *>(&(this_entry->value)));
#endif
        }
      }
      // if unique
      for (; k < 16; k++) {
        if (hashPolicy<T>::get_bit_true(bitmap, k)) {
#ifdef AU
          unique_entry_num++;
#endif
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
  std::cout << "pointer number = " << pointer_num << std::endl;
  std::cout << "unique entry num = " << unique_entry_num << std::endl;
  LOG("Value type count: ");
  for (int i = 0; i < 32; ++i) {
    std::cout << type_cnt[i] << ", ";
  }
  std::cout << std::endl;
  LOG("Value average count: ");
  for (int i = 0; i < 32; ++i) {
    std::cout << setprecision(14) << type_avg[i] << ", ";
  }
  std::cout << std::endl;
#endif
  return count;
}


template <class T>
void PeaHashing<T>::emphasizeDir() {
  LOG("emphasize dir!");
  seg_meta *dir_entry = dir->_;
  uint8_t g = *global_depth;
  int capacity = pow(2, g);
  uint64_t sum = 0;
  for (size_t i = 0; i < capacity; i++) {
    seg_meta *this_dir_entry = dir_entry + i;
    LINE_PREF(this_dir_entry);
    sum += *(uint64_t *)this_dir_entry;
  } 
  LOG(sum);
}

template <class T>
PeaHashing<T>::PeaHashing(void) {
  std::cout << "Reinitialize up" << std::endl;
}

template <class T>
int PeaHashing<T>::Insert(T key, Value_t value, bool is_in_epoch) {
  return Insert(key, value);
}

template <class T>
bool PeaHashing<T>::Delete(T key, bool flag) {
  return Delete(key);
}
}  // namespace peas
