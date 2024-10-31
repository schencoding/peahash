
// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong. All rights reserved.
// Licensed under the MIT license.
#ifndef HASH_INTERFACE_H_
#define HASH_INTERFACE_H_

#include "../util/pair.h"
#ifdef PMEM
#include <libpmemobj.h>
#endif

/*
* Parent function of all hash indexes
* Used to define the interface of the hash indexes
*/

template <class T>
class Hash {
 public:
  Hash(void) = default;
  ~Hash(void) = default;
  virtual int Insert(T, Value_t) {
    return -1;
  }
  virtual int Insert(T, Value_t, bool) {
    return -1;
  }
 
  virtual void bootRestore(){

  };
  virtual void reportRestore(){

  };
  virtual bool Delete(T) {
    return false;
  }
  virtual bool Delete(T, bool) {
    return false;
  }
  virtual Value_t Get(T) {
    return NONE;
  }
  virtual Value_t Get(T key, bool is_in_epoch) {
    return NONE;
  }
  virtual int Get(T key, Value_t *ret_val, size_t *number) {
    return 0;
  }
  virtual void Recovery() {}
  virtual uint64_t getNumber() {
    return 0;
  }
};

#endif  // _HASH_INTERFACE_H_
