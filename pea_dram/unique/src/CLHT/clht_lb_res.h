#pragma once

#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "clht_lb_base.h"
#include <cassert>

namespace clht_lb {

__thread size_t check_ht_status_steps = CLHT_STATUS_INVOK_IN;


inline int
is_power_of_two (unsigned int x)
{
  return ((x != 0) && !(x & (x - 1)));
}

static inline int
is_odd (int x)
{
  return x & 1;
}

/* Create a new bucket. */
bucket_t *
clht_bucket_create ()
{
  bucket_t *bucket = NULL;
  bucket = (bucket_t *) memalign (CACHE_LINE_SIZE, sizeof (bucket_t));
  if (bucket == NULL)
    {
      return NULL;
    }

  bucket->lock = 0;

  uint32_t j;
  for (j = 0; j < ENTRIES_PER_BUCKET; j++)
    {
      bucket->key[j] = 0;
    }
  bucket->next = NULL;

  return bucket;
}

bucket_t *
clht_bucket_create_stats (clht_hashtable_t * h, int *resize)
{
  bucket_t *b = clht_bucket_create ();
  if (IAF_U32 (&h->num_expands) == h->num_expands_threshold)
    {
      /* printf("      -- hit threshold (%u ~ %u)\n", h->num_expands,
       * h->num_expands_threshold); */
      *resize = 1;
    }
  return b;
}

clht_hashtable_t *clht_hashtable_create (uint64_t num_buckets);

clht_t *
clht_create (uint64_t num_buckets)
{
  clht_t *w = (clht_t *) memalign (CACHE_LINE_SIZE, sizeof (clht_t));
  if (w == NULL)
    {
      printf ("** malloc @ hatshtalbe\n");
      return NULL;
    }

  w->ht = clht_hashtable_create (num_buckets);
  if (w->ht == NULL)
    {
      free (w);
      return NULL;
    }
  w->resize_lock = LOCK_FREE;
  w->gc_lock = LOCK_FREE;
  w->status_lock = LOCK_FREE;
  w->version_list = NULL;
  w->version_min = 0;
  w->ht_oldest = w->ht;

  return w;
}

clht_hashtable_t *
clht_hashtable_create (uint64_t num_buckets)
{
  clht_hashtable_t *hashtable = NULL;

  if (num_buckets == 0)
    {
      return NULL;
    }

  /* Allocate the table itself. */
  hashtable =
    (clht_hashtable_t *) memalign (CACHE_LINE_SIZE,
				   sizeof (clht_hashtable_t));
  if (hashtable == NULL)
    {
      printf ("** malloc @ hatshtalbe\n");
      return NULL;
    }

  /* hashtable->table = calloc(num_buckets, (sizeof(bucket_t))); */
  hashtable->table =
    (bucket_t *) memalign (CACHE_LINE_SIZE,
			   num_buckets * (sizeof (bucket_t)));
  if (hashtable->table == NULL)
    {
      printf ("** alloc: hashtable->table\n");
      fflush (stdout);
      free (hashtable);
      return NULL;
    }

  memset (hashtable->table, 0, num_buckets * (sizeof (bucket_t)));

  uint64_t i;
  for (i = 0; i < num_buckets; i++)
    {
      hashtable->table[i].lock = LOCK_FREE;
      uint32_t j;
      for (j = 0; j < ENTRIES_PER_BUCKET; j++)
	{
	  hashtable->table[i].key[j] = 0;
	}
    }

  hashtable->num_buckets = num_buckets;
  hashtable->hash = num_buckets - 1;
  hashtable->version = 0;
  hashtable->table_tmp = NULL;
  hashtable->table_new = NULL;
  hashtable->table_prev = NULL;
  hashtable->num_expands = 0;
  hashtable->num_expands_threshold = (CLHT_PERC_EXPANSIONS * num_buckets);
  if (hashtable->num_expands_threshold == 0)
    {
      hashtable->num_expands_threshold = 1;
    }
  // LOG("num_expands_threshold: " << hashtable->num_expands_threshold);
  hashtable->is_helper = 1;
  hashtable->helper_done = 0;

  return hashtable;
}


uint64_t clht_hash(clht_hashtable_t* hashtable, clht_addr_t key) {
  return h(&key, sizeof(key)) & (hashtable->hash);
}

/* Retrieve a key-value entry from a hash table. */
clht_val_t clht_get(clht_hashtable_t* hashtable, clht_addr_t key) {
  size_t bin = clht_hash(hashtable, key);
  CLHT_GC_HT_VERSION_USED(hashtable);
  volatile bucket_t* bucket = hashtable->table + bin;

  uint32_t j;
  do {
    for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
      clht_val_t val = bucket->val[j];
#ifdef __tile__
      _mm_lfence();
#endif
      if (bucket->key[j] == key) {
        if (likely(bucket->val[j] == val)) {
          return val;
        } else {
          return 0;
        }
      }
    }

    bucket = bucket->next;
  } while (unlikely(bucket != NULL));
  return 0;
}

static inline int bucket_exists(volatile bucket_t* bucket, clht_addr_t key) {
  uint32_t j;
  do {
    for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
      if (bucket->key[j] == key) {
        return true;
      }
    }
    bucket = bucket->next;
  } while (unlikely(bucket != NULL));
  return false;
}

static uint32_t clht_put_seq(clht_hashtable_t* hashtable, clht_addr_t key,
                             clht_val_t val, uint64_t bin) {
  volatile bucket_t* bucket = hashtable->table + bin;
  uint32_t j;

  do {
    for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
      if (bucket->key[j] == 0) {
        bucket->val[j] = val;
        bucket->key[j] = key;
        return true;
      }
    }

    if (bucket->next == NULL) {
      DPP(put_num_failed_expand);
      int null;
      bucket->next = clht_bucket_create_stats(hashtable, &null);
      bucket->next->val[0] = val;
      bucket->next->key[0] = key;
      return true;
    }

    bucket = bucket->next;
  } while (true);
}

static int bucket_cpy(volatile bucket_t* bucket, clht_hashtable_t* ht_new) {
  if (!LOCK_ACQ_RES(&bucket->lock)) {
    return 0;
  }
  uint32_t j;
  do {
    for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
      clht_addr_t key = bucket->key[j];
      if (key != 0) {
        uint64_t bin = clht_hash(ht_new, key);
        clht_put_seq(ht_new, key, bucket->val[j], bin);
      }
    }
    bucket = bucket->next;
  } while (bucket != NULL);

  return 1;
}

size_t ht_status(clht_t* h, int resize_increase, int just_print);

int ht_resize_pes(clht_t* h, int is_increase, int by) {
  // LOG("by: " << by);
  check_ht_status_steps = CLHT_STATUS_INVOK;

  clht_hashtable_t* ht_old = h->ht;

  if (TRYLOCK_ACQ(&h->resize_lock)) {
    return 0;
  }

  size_t num_buckets_new;
  if (is_increase == true) {
    /* num_buckets_new = CLHT_RATIO_DOUBLE * ht_old->num_buckets; */
    num_buckets_new = by * ht_old->num_buckets;
  } else {
#if CLHT_HELP_RESIZE == 1
    ht_old->is_helper = 0;
#endif
    num_buckets_new = ht_old->num_buckets / CLHT_RATIO_HALVE;
  }

  printf("// resizing: from %8zu to %8zu buckets\n", ht_old->num_buckets,
  num_buckets_new);

  clht_hashtable_t* ht_new = clht_hashtable_create(num_buckets_new);
  ht_new->version = ht_old->version + 1;

#if CLHT_HELP_RESIZE == 1
  ht_old->table_tmp = ht_new;

  int32_t b;
  for (b = 0; b < ht_old->num_buckets; b++) {
    bucket_t* bu_cur = ht_old->table + b;
    if (!bucket_cpy(bu_cur,
                    ht_new)) /* reached a point where the helper is handling */
    {
      break;
    }
  }

  if (is_increase && ht_old->is_helper != 1) /* there exist a helper */
  {
    while (ht_old->helper_done != 1) {
      _mm_pause();
    }
  }

#else

  int32_t b;
  for (b = 0; b < ht_old->num_buckets; b++) {
    bucket_t* bu_cur = ht_old->table + b;
    bucket_cpy(bu_cur, ht_new);
  }
#endif

#if defined(DEBUG)
  /* if (clht_size(ht_old) != clht_size(ht_new)) */
  /*   { */
  /*     printf("**clht_size(ht_old) = %zu != clht_size(ht_new) = %zu\n",
   * clht_size(ht_old), clht_size(ht_new)); */
  /*   } */
#endif

  ht_new->table_prev = ht_old;

  int ht_resize_again = 0;
  if (ht_new->num_expands >= ht_new->num_expands_threshold) {
    /* printf("--problem: have already %u expands\n", ht_new->num_expands); */
    ht_resize_again = 1;
    /* ht_new->num_expands_threshold = ht_new->num_expands + 1; */
  }

  SWAP_U64((uint64_t*)h, (uint64_t)ht_new);
  ht_old->table_new = ht_new;
  TRYLOCK_RLS(h->resize_lock);

#if CLHT_DO_GC == 1
  clht_gc_collect(h);
#else
  // clht_gc_release(ht_old);
#endif

  if (ht_resize_again) {
    ht_status(h, 1, 0);
  }

  return 1;
}

size_t ht_status(clht_t* h, int resize_increase, int just_print) {
  if (TRYLOCK_ACQ(&h->status_lock) && !resize_increase) {
    return 0;
  }

  clht_hashtable_t* hashtable = h->ht;
  uint64_t num_buckets = hashtable->num_buckets;
  volatile bucket_t* bucket = NULL;
  size_t size = 0;
  int expands = 0;
  int expands_max = 0;

  uint64_t bin;
  for (bin = 0; bin < num_buckets; bin++) {
    bucket = hashtable->table + bin;

    int expands_cont = -1;
    expands--;
    uint32_t j;
    do {
      expands_cont++;
      expands++;
      for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
        if (bucket->key[j] > 0) {
          size++;
        }
      }

      bucket = bucket->next;
    } while (bucket != NULL);

    if (expands_cont > expands_max) {
      expands_max = expands_cont;
    }
  }

  double full_ratio =
      100.0 * size / ((hashtable->num_buckets) * ENTRIES_PER_BUCKET);

  // LOG("full ratio: " << full_ratio);

  if (just_print) {
    printf(
        "[STATUS-%02d] #bu: %7zu / #elems: %7zu / full%%: %8.4f%% / expands: "
        "%4d / max expands: %2d\n",
        99, hashtable->num_buckets, size, full_ratio, expands, expands_max);
  } else {
    if (full_ratio > 0 && full_ratio < CLHT_PERC_FULL_HALVE) {
      // printf(
      //     "[STATUS-%02d] #bu: %7zu / #elems: %7zu / full%%: %8.4f%% /
      //     expands: "
      //     "%4d / max expands: %2d\n",
      //     clht_gc_get_id(), hashtable->num_buckets, size, full_ratio,
      //     expands, expands_max);
      ht_resize_pes(h, 0, 33);
    } else if ((full_ratio > 0 && full_ratio > CLHT_PERC_FULL_DOUBLE) ||
               expands_max > CLHT_MAX_EXPANSIONS || resize_increase) {
      int inc_by = (full_ratio / CLHT_OCCUP_AFTER_RES);
      int inc_by_pow2 = pow2roundup(inc_by);

      // printf(
      //     "[STATUS-%02d] #bu: %7zu / #elems: %7zu / full%%: %8.4f%% /
      //     expands: "
      //     "%4d / max expands: %2d\n",
      //     clht_gc_get_id(), hashtable->num_buckets, size, full_ratio,
      //     expands, expands_max);
      if (inc_by_pow2 == 1) {
        inc_by_pow2 = 2;
      }
      ht_resize_pes(h, 1, inc_by_pow2);
    }
  }

  // if (!just_print) {
  //   clht_gc_collect(h);
  // }

  TRYLOCK_RLS(h->status_lock);
  return size;
}

/* Insert a key-value entry into a hash table. */
int clht_put(clht_t* h, clht_addr_t key, clht_val_t val) {
  clht_hashtable_t* hashtable = h->ht;
  size_t bin = clht_hash(hashtable, key);
  volatile bucket_t* bucket = hashtable->table + bin;

#if CLHT_READ_ONLY_FAIL == 1
  if (bucket_exists(bucket, key)) {
    return false;
  }
#endif

  clht_lock_t* lock = &bucket->lock;
  while (!LOCK_ACQ(lock, hashtable)) {
    hashtable = h->ht;
    size_t bin = clht_hash(hashtable, key);

    bucket = hashtable->table + bin;
    lock = &bucket->lock;
  }

  CLHT_GC_HT_VERSION_USED(hashtable);
  CLHT_CHECK_STATUS(h);
  clht_addr_t* empty = NULL;
  clht_val_t* empty_v = NULL;

  uint32_t j;
  do {
    for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
      if (bucket->key[j] == key) {
        LOCK_RLS(lock);
        return false;
      } else if (empty == NULL && bucket->key[j] == 0) {
        empty = (clht_addr_t*)&bucket->key[j];
        empty_v = &bucket->val[j];
      }
    }

    int resize = 0;
    if (likely(bucket->next == NULL)) {
      if (unlikely(empty == NULL)) {
        DPP(put_num_failed_expand);

        bucket_t* b = clht_bucket_create_stats(hashtable, &resize);
        b->val[0] = val;
#ifdef __tile__
        /* keep the writes in order */
        _mm_sfence();
#endif
        b->key[0] = key;
#ifdef __tile__
        /* make sure they are visible */
        _mm_sfence();
#endif
        bucket->next = b;
      } else {
        *empty_v = val;
#ifdef __tile__
        /* keep the writes in order */
        _mm_sfence();
#endif
        *empty = key;
      }

      LOCK_RLS(lock);
      if (unlikely(resize)) {
        /* ht_resize_pes(h, 1); */
        ht_status(h, 1, 0);
      }
      return true;
    }
    bucket = bucket->next;
  } while (true);
}

/* Remove a key-value entry from a hash table. */
clht_val_t clht_remove(clht_t* h, clht_addr_t key) {
  clht_hashtable_t* hashtable = h->ht;
  size_t bin = clht_hash(hashtable, key);
  volatile bucket_t* bucket = hashtable->table + bin;

#if CLHT_READ_ONLY_FAIL == 1
  if (!bucket_exists(bucket, key)) {
    return false;
  }
#endif

  clht_lock_t* lock = &bucket->lock;
  while (!LOCK_ACQ(lock, hashtable)) {
    hashtable = h->ht;
    size_t bin = clht_hash(hashtable, key);

    bucket = hashtable->table + bin;
    lock = &bucket->lock;
  }

  CLHT_GC_HT_VERSION_USED(hashtable);
  CLHT_CHECK_STATUS(h);

  uint32_t j;
  do {
    for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
      if (bucket->key[j] == key) {
        clht_val_t val = bucket->val[j];
        bucket->key[j] = 0;
        LOCK_RLS(lock);
        return val;
      }
    }
    bucket = bucket->next;
  } while (unlikely(bucket != NULL));
  LOCK_RLS(lock);
  return false;
}

void ht_resize_help(clht_hashtable_t* h) {
  if ((int32_t)FAD_U32((volatile uint32_t*)&h->is_helper) <= 0) {
    return;
  }

  int32_t b;
  /* hash = num_buckets - 1 */
  for (b = h->hash; b >= 0; b--) {
    bucket_t* bu_cur = h->table + b;
    if (!bucket_cpy(
            bu_cur,
            h->table_tmp)) { /* reached a point where the resizer is handling */
      /* printf("[GC-%02d] helped  #buckets: %10zu = %5.1f%%\n",  */
      /* 	 clht_gc_get_id(), h->num_buckets - b, 100.0 * (h->num_buckets
       * - b) / h->num_buckets); */
      break;
    }
  }

  h->helper_done = 1;
}

size_t clht_size(clht_hashtable_t* hashtable) {
  uint64_t num_buckets = hashtable->num_buckets;
  volatile bucket_t* bucket = NULL;
  size_t size = 0;

  uint64_t bin;
  for (bin = 0; bin < num_buckets; bin++) {
    bucket = hashtable->table + bin;

    uint32_t j;
    do {
      for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
        if (bucket->key[j] > 0) {
          size++;
        }
      }

      bucket = bucket->next;
    } while (bucket != NULL);
  }
  return size;
}

size_t clht_size_mem(clht_hashtable_t* h) /* in bytes */
{
  if (h == NULL) {
    return 0;
  }

  size_t size_tot = sizeof(clht_hashtable_t**);
  size_tot += (h->num_buckets + h->num_expands) * sizeof(bucket_t);
  return size_tot;
}

size_t clht_size_mem_garbage(clht_hashtable_t* h) /* in bytes */
{
  if (h == NULL) {
    return 0;
  }

  size_t size_tot = 0;
  clht_hashtable_t* cur = h->table_prev;
  while (cur != NULL) {
    size_tot += clht_size_mem(cur);
    cur = cur->table_prev;
  }

  return size_tot;
}

void clht_print(clht_hashtable_t* hashtable) {
  uint64_t num_buckets = hashtable->num_buckets;
  volatile bucket_t* bucket;

  printf("Number of buckets: %zu\n", num_buckets);
  uint64_t count = 0;
  uint64_t bin;
  for (bin = 0; bin < num_buckets; bin++) {
    bucket = hashtable->table + bin;

    // printf("[[%05zu]] ", bin);

    uint32_t j;
    do {
      count++;
      for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
        if (bucket->key[j]) {
          // printf("(%-5llu/%p)-> ", (long long unsigned int)bucket->key[j],
          //        (void*)bucket->val[j]);
        }
      }

      bucket = bucket->next;
      // printf(" ** -> ");
    } while (bucket != NULL);
    // printf("\n");
  }
  // LOG("bucket_count: " << count);
  fflush(stdout);
}


template <class T>
class CLHT : public Hash<T> {
  clht_t* hashtable;

 public:
  CLHT(void){
    std::cout << "Reintialize Up for CLHT" << std::endl;
  }
  CLHT(uint64_t num_buckets) {
    hashtable = clht_create(num_buckets);
    assert(hashtable != NULL);
    printf("#BEFORE size: %zu\n", clht_size(hashtable->ht));
    /* clht_print(hashtable, num_buckets); */
    // LOG("Size of key: " << sizeof(clht_addr_t));
    // LOG("Size of value: " << sizeof(clht_val_t));
    // LOG("Size of bucket: " << sizeof(bucket_t));
    assert(sizeof(bucket_t) % 64 == 0);
  }
  ~CLHT(void){}
  int Insert(T key, Value_t value) {
    return clht_put(hashtable, reinterpret_cast<clht_addr_t>(key), reinterpret_cast<clht_val_t>(value));
  }

  int Insert(T key, Value_t value, bool epoch) {
    return Insert(key, value);
  }
  bool Delete(T key) {
    int removed = clht_remove(hashtable, reinterpret_cast<clht_addr_t>(key));
    return (removed == 0)? false: true;
  }
  bool Delete(T key, bool epoch){
    return Delete(key);
  }
  Value_t Get(T key) {
    return reinterpret_cast<Value_t>(clht_get(hashtable->ht, reinterpret_cast<clht_addr_t>(key)));
  }
  Value_t Get(T key, bool is_in_epoch) {
    return Get(key);
  }

  void Recovery(void){
    return;
  }

  uint64_t getNumber() {
    clht_print(hashtable->ht);
    return 0;
  }
  int Get(T key, Value_t *ret_val, size_t *number){
    return 0;
  }
};



} // namespace