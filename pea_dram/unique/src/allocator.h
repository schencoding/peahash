
// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong. All rights reserved.
// Licensed under the MIT license.
#pragma once
#include <sys/mman.h>
#include "../util/utils.h"
#include "x86intrin.h"
#include <cstring>

struct Allocator {
 public:
  static Allocator* instance_;
  static Allocator* Get() { return instance_; }

  static void Allocate(void** ptr, uint32_t alignment, size_t size) {
    posix_memalign(ptr, alignment, size);
  }

  static void ZAllocate(void** ptr, uint32_t alignment, size_t size) {
    int ret = posix_memalign(ptr, alignment, size);
    if (ret) {
      fprintf (stderr, "posix_memalign: %s\n",
             strerror (ret));
    }
    memset(*ptr, 0, size);
  }
};

Allocator* Allocator::instance_ = nullptr;