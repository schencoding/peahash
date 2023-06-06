//
// Created by liu on 2020-11-08.
//

#ifndef PEHASH_PREF_H
#define PEHASH_PREF_H

/* --------------------------------------------------------------------- */
/*                      Prefetch Instructions                            */
/* --------------------------------------------------------------------- */
/*
    T0 (temporal data):
        prefetch data into all levels of the cache hierarchy.
        Pentium III processor-1st- or 2nd-level cache.
        Pentium 4 and Intel Xeon processors-2nd-level cache.
    NTA (non-temporal data with respect to all cache levels):
        prefetch data into non-temporal cache structure and
        into a location close to the processor, minimizing cache pollution.
 */

#ifndef NO_PREFETCH

#define prefetcht0(mem_var)     \
        __asm__ __volatile__ ("prefetcht0 %0": :"m"(mem_var))
#define prefetchnta(mem_var)    \
        __asm__ __volatile__ ("prefetchnta %0": :"m"(mem_var))

#define pref(mem_var)      prefetcht0(mem_var)
#define prefnta(mem_var)   prefetchnta(mem_var)

#else // NO_PREFETCH

#define pref(mem_var)
#define prefnta(mem_var)

#endif

/* ---------------------------------------------------------------------- */

#define BUCKET_LINE_NUM 4
#define CACHE_LINE_SIZE 64

static void inline LINE_PREF(/*register*/ void *bp){        //from base pointer
    pref (* ((char *)bp));
}

static void inline DOUBLE_LINE_PREF(/*register*/ void *bp){ //from base pointer
    pref (* ((char *)bp));
    pref (* ((char *)bp + CACHE_LINE_SIZE));
}

// for bucket
static void inline BUCKET_PREF(/*register*/ void *bp)       //from base pointer
{
    pref (* ((char *)bp));
#    if BUCKET_LINE_NUM >= 2
    pref (* ((char *)bp + CACHE_LINE_SIZE));
#    endif
#    if BUCKET_LINE_NUM >= 3
    pref (* ((char *)bp + CACHE_LINE_SIZE*2));
#    endif
#    if BUCKET_LINE_NUM >= 4
    pref (* ((char *)bp + CACHE_LINE_SIZE*3));
#    endif
#    if BUCKET_LINE_NUM >= 5
    pref (* ((char *)bp + CACHE_LINE_SIZE*4));
#    endif
#    if BUCKET_LINE_NUM >= 6
    pref (* ((char *)bp + CACHE_LINE_SIZE*5));
#    endif
#    if BUCKET_LINE_NUM >= 7
    pref (* ((char *)bp + CACHE_LINE_SIZE*6));
#    endif
#    if BUCKET_LINE_NUM >= 8
    pref (* ((char *)bp + CACHE_LINE_SIZE*7));
#    endif
#    if BUCKET_LINE_NUM >= 9
    pref (* ((char *)bp + CACHE_LINE_SIZE*8));
#    endif
#    if BUCKET_LINE_NUM >= 10
    pref (* ((char *)bp + CACHE_LINE_SIZE*9));
#    endif
#    if BUCKET_LINE_NUM >= 11
    pref (* ((char *)bp + CACHE_LINE_SIZE*10));
#    endif
#    if BUCKET_LINE_NUM >= 12
    pref (* ((char *)bp + CACHE_LINE_SIZE*11));
#    endif
#    if BUCKET_LINE_NUM >= 13
    pref (* ((char *)bp + CACHE_LINE_SIZE*12));
#    endif
#    if BUCKET_LINE_NUM >= 14
    pref (* ((char *)bp + CACHE_LINE_SIZE*13));
#    endif
#    if BUCKET_LINE_NUM >= 15
    pref (* ((char *)bp + CACHE_LINE_SIZE*14));
#    endif
#    if BUCKET_LINE_NUM >= 16
    pref (* ((char *)bp + CACHE_LINE_SIZE*15));
#    endif
#    if BUCKET_LINE_NUM >= 17
    pref (* ((char *)bp + CACHE_LINE_SIZE*16));
#    endif
#    if BUCKET_LINE_NUM >= 18
    pref (* ((char *)bp + CACHE_LINE_SIZE*17));
#    endif
#    if BUCKET_LINE_NUM >= 19
    pref (* ((char *)bp + CACHE_LINE_SIZE*18));
#    endif
#    if BUCKET_LINE_NUM >= 20
    pref (* ((char *)bp + CACHE_LINE_SIZE*19));
#    endif
#    if BUCKET_LINE_NUM >= 21
    pref (* ((char *)bp + CACHE_LINE_SIZE*20));
#    endif
#    if BUCKET_LINE_NUM >= 22
    pref (* ((char *)bp + CACHE_LINE_SIZE*21));
#    endif
#    if BUCKET_LINE_NUM >= 23
    pref (* ((char *)bp + CACHE_LINE_SIZE*22));
#    endif
#    if BUCKET_LINE_NUM >= 24
    pref (* ((char *)bp + CACHE_LINE_SIZE*23));
#    endif
#    if BUCKET_LINE_NUM >= 25
    pref (* ((char *)bp + CACHE_LINE_SIZE*24));
#    endif
#    if BUCKET_LINE_NUM >= 26
    pref (* ((char *)bp + CACHE_LINE_SIZE*25));
#    endif
#    if BUCKET_LINE_NUM >= 27
    pref (* ((char *)bp + CACHE_LINE_SIZE*26));
#    endif
#    if BUCKET_LINE_NUM >= 28
    pref (* ((char *)bp + CACHE_LINE_SIZE*27));
#    endif
#    if BUCKET_LINE_NUM >= 29
    pref (* ((char *)bp + CACHE_LINE_SIZE*28));
#    endif
#    if BUCKET_LINE_NUM >= 30
    pref (* ((char *)bp + CACHE_LINE_SIZE*29));
#    endif
#    if BUCKET_LINE_NUM >= 31
    pref (* ((char *)bp + CACHE_LINE_SIZE*30));
#    endif
#    if BUCKET_LINE_NUM >= 32
    pref (* ((char *)bp + CACHE_LINE_SIZE*31));
#    endif
#    if BUCKET_LINE_NUM >= 33
#    error "BUCKET_LINE_NUM must be <= 32!"
#    endif
}

#endif //PEHASH_PREF_H
