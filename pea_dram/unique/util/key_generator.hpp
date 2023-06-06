
// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong.
// All rights reserved. Licensed under the MIT license.

#pragma once

#include "random.h"
#include "zipfian_int_distribution.hpp"
#include <optional>

/*
 * Used to generate different key distribution used in the mini benchmark.
 * Including: (1) uniform random distribution (2) zipfian distribution (3) range
 * distribution (generate the interger one by one in a range specified by the
 * user)
 */

class key_generator_t {
 public:
  key_generator_t() {}

  virtual uint64_t next_uint64() = 0;
};

class uniform_key_generator_t : public key_generator_t {
 public:
  uniform_key_generator_t() {
    unsigned long long init[4] = {0x12345ULL, 0x23456ULL, 0x34567ULL,
                                  0x45678ULL},
                       length = 4;
    init_by_array64(init, length);
  }

  uint64_t next_uint64() { return genrand64_int64(); }
};

class range_key_generator_t : public key_generator_t {
 public:
  range_key_generator_t(size_t start) { start_number = start; }

  uint64_t next_uint64() { return start_number++; }

 private:
  size_t start_number = 0;
};

class zipfian_key_generator_t : public key_generator_t {
 public:
  // control the seed to test different skew
  // zipfian_key_generator_t(size_t start, size_t end, float skew) :
  // dist_(start, end, skew), generator_(time(0)){

  // }
  zipfian_key_generator_t(size_t start, size_t end, float skew)
      : dist_(start, end, skew),
        u(0.0, 1.0),
        generator_(1234),
        start_(start),
        end_(end),
        skew_alp(skew) {}

  uint64_t next_uint64() { return dist_(generator_); }

  uint64_t next_uint64_other() {
    return start_ + zipf(skew_alp, (int)(end_ - start_));
  }

  uint64_t next_uint64_skew(std::optional<uint64_t> seed = std::nullopt) {
    if (seed == std::nullopt) {
      seed = (uint16_t)time(0) % (end_ - start_);
    }
    uint64_t origin = dist_(generator_);
    uint64_t retval = origin + seed.value();
    if (retval > end_) {
      retval = retval - end_ - 1 + start_;
    }
    return retval;
  }

  int zipf(double alpha, int n) {
    static int first = true;   // Static first time flag
    static double c = 0;       // Normalization constant
    static double *sum_probs;  // Pre-calculated sum of probabilities
    double z;                  // Uniform random number (0 < z < 1)
    int zipf_value;            // Computed exponential value to be returned
    int i;                     // Loop counter
    int low, high, mid;        // Binary-search bounds

    // Compute normalization constant on first call only
    if (first == true) {
      for (i = 1; i <= n; i++) c = c + (1.0 / pow((double)i, alpha));
      c = 1.0 / c;

      sum_probs = static_cast<double *>(malloc((n + 1) * sizeof(*sum_probs)));
      sum_probs[0] = 0;
      for (i = 1; i <= n; i++) {
        sum_probs[i] = sum_probs[i - 1] + c / pow((double)i, alpha);
      }
      first = false;
    }

    // Pull a uniform random number (0 < z < 1)
    do {
      z = u(generator_);
    } while ((z == 0) || (z == 1));

    // Map z to the value
    low = 1, high = n, mid;
    do {
      mid = floor((low + high) / 2);
      if (sum_probs[mid] >= z && sum_probs[mid - 1] < z) {
        zipf_value = mid;
        break;
      } else if (sum_probs[mid] >= z) {
        high = mid - 1;
      } else {
        low = mid + 1;
      }
    } while (low <= high);

    // Assert that zipf_value is between 1 and N
    assert((zipf_value >= 1) && (zipf_value <= n));

    return (zipf_value);
  }

 private:
  std::default_random_engine generator_;
  std::uniform_real_distribution<double> u;
  double skew_alp;
  zipfian_int_distribution<uint64_t> dist_;
  size_t start_;
  size_t end_;
};
