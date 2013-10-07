/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef RANDOM_INCLUDED
#define RANDOM_INCLUDED

#ifndef PLATFORM_DARWIN
#include <x86intrin.h>
#endif

#include <pmmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include "sse.hpp"

#define SFMT_MEXP 19937
#include "SFMT.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

#define USE_SSE2
#include "sse_mathfun.h"

// RNG wraps Boost Random Number Generation
// RNG_SFMT wraps SFMT by Saito and Matsumoto
//  Provides block generation and buffering for performance, plus commonly-used accessor functions
//  eg. unit vectors, different types, etc
//  Produces its floating-point types by direct bit manipulation of the random number instead of division (faster)


typedef union SSEReg_t {
    __m128  m128_f;
    __m128i m128_si; 
    __m128d m128_d;
    float f[4];
    double d[2];
    uint8_t i8[16];
    uint32_t i32[4];
    uint64_t i64[2];
} SSEReg;

class RNG {
    boost::random::mt19937 rng;
    boost::random::uniform_01<double> uni01d;
    boost::random::uniform_01<float>  uni01f;
    public:
    // requirements for Boost RNG concept
    typedef boost::random::mt19937::result_type result_type;
    result_type min() const { return rng.min(); }
    result_type max() const { return rng.max(); }
    result_type operator()(){ return rng(); } 


    inline double draw_double_u01();
    inline float  draw_float_u01();
    inline __m128 draw_m128f3_u01();
    inline __m128 draw_m128f4_u01();
    inline __m128 draw_m128f3_uvect();
};


// The size of array (32b words) generated must be at least 624 and a multiple of 4
// Nbuf must be at least 156

class RNG_SFMT {
    SSEReg * const randBuf,* const lastRand;
    SSEReg *nextRand;
    unsigned Nbuf;
    sfmt_t sfmt;

    float       __attribute__((aligned(16))) f_log[4];
    float       __attribute__((aligned(16))) f[4];
    double      __attribute__((aligned(16))) d[2];
    SSEReg i32;

    unsigned char f_count,d_count,i32_count,f_log_count;

    void refill();
    inline SSEReg  draw();
    inline SSEReg* drawPtr();
    
    public:

    static const unsigned int  exp_float24 = 0x40000000;
    static const unsigned int  exp_float =0x3f800000;
    static const uint32_t exp_double_h = 0x3ff00000;
    static const uint32_t exp_double_l = 0x00000000;

    RNG_SFMT(unsigned int Nbuf_=1024,unsigned seed_=1) : randBuf(new SSEReg[Nbuf_]),
        lastRand(randBuf+Nbuf_),
        nextRand(randBuf),
        Nbuf(Nbuf_),
        f_count(4),
        d_count(2),
        i32_count(4),
        f_log_count(4)
        {
//		std::cout << "Initialized SFMT RNG with seed of " << seed_ << std::endl;
            sfmt_init_gen_rand(&sfmt,seed_);
            refill();
        }

    // requirements for Boost RNG concept
    typedef uint32_t result_type;
    result_type min() const { return std::numeric_limits<uint32_t>::min(); }
    result_type max() const { return std::numeric_limits<uint32_t>::max(); }
    result_type operator()(){ return draw_uint32(); } 

    void refill_log_u01();

    inline uint32_t draw_uint32();
    inline double draw_double_u01();
    inline float  draw_float_u01();
    inline __m128 draw_m128f3_u01();
    inline __m128 draw_m128f4_u01();
    inline __m128 draw_m128f1_u01();

    inline __m128 draw_m128f1_log_u01();

    inline __m128d draw_m128d1_u01();
    inline __m128d draw_m128d2_u01();

    inline __m128 draw_m128f4_pm1();

    inline __m128 draw_m128f2_uvect();
    inline __m128 draw_m128f3_uvect();

    inline const uint64_t* draw_u64_2();
    inline const uint32_t* draw_u32_4();
};

inline SSEReg RNG_SFMT::draw()
{
    if (nextRand == lastRand)
        refill();
    return *(nextRand++);
}

inline SSEReg* RNG_SFMT::drawPtr()
{
    if (nextRand == lastRand)
        refill();
    return nextRand++;
}

inline __m128d RNG_SFMT::draw_m128d2_u01()
{
    __m128d rnd=_mm_castps_pd(draw().m128_f);
    __m128d  exp = _mm_castsi128_pd(_mm_set_epi32(exp_double_h,exp_double_l,exp_double_h,exp_double_l));
    __m128d offs = _mm_set1_pd(1.0);
    rnd = _mm_castsi128_pd(_mm_srli_epi64(_mm_castpd_si128(rnd),12));
    rnd = _mm_or_pd(rnd,exp);
    __m128d res = _mm_sub_pd(rnd,offs);
    return res;
}

// pm1 -> plus/minus 1

inline __m128 RNG_SFMT::draw_m128f4_pm1()
{
    __m128 rnd=draw().m128_f;
    __m128  exp = __m128(_mm_set_epi32(exp_float24,exp_float24,exp_float24,exp_float24));
    __m128 offs = _mm_set1_ps(3.0);
    rnd = _mm_castsi128_ps(_mm_srli_epi32(_mm_castps_si128(rnd),9));
    rnd = _mm_or_ps(rnd,exp);
    rnd = _mm_sub_ps(rnd,offs);
    return rnd;
}

inline __m128 RNG_SFMT::draw_m128f4_u01()
{
    __m128 rnd = draw().m128_f;
    __m128  exp = __m128(_mm_set1_epi32(exp_float));
    __m128 offs = _mm_set1_ps(1.0);
    rnd = __m128(_mm_srli_epi32(_mm_castps_si128(rnd),9));
    rnd = _mm_or_ps(rnd,exp);
    rnd = _mm_sub_ps(rnd,offs);
    return rnd;
}

inline __m128 RNG_SFMT::draw_m128f3_u01()
{
    __m128 rnd = draw().m128_f;
    __m128  exp = __m128(_mm_set_epi32(0,exp_float,exp_float,exp_float));
    __m128 offs = _mm_set_ps(0,1.0,1.0,1.0);
    rnd = _mm_castsi128_ps(_mm_srli_epi32(_mm_castps_si128(rnd),9));
    rnd = _mm_or_ps(rnd,exp);
    rnd = _mm_sub_ps(rnd,offs);
    return rnd;
}
/*
const uint64_t* RNG_SFMT::draw_u64_2()
{
    if(nextRand == lastRand)
        refill();
    return (nextRand++)->i64;
}

const uint32_t* RNG_SFMT::draw_u32_4()
{
    if (nextRand == lastRand)
        refill();
    return (nextRand++)->i32;
}*/

__m128 RNG_SFMT::draw_m128f1_log_u01()
{
    __m128 r,l;
    if (f_log_count == 4)
    {
        r = _mm_sub_ps(_mm_set1_ps(1.0),draw_m128f4_u01());     // 1 - [0,1) => (0,1] to avoid -Inf
        l = log_ps(r);
        _mm_store_ps(f_log,l);
        f_log_count=0;
    }
    return _mm_load1_ps(f_log+(f_log_count++));
}

// Return a single float
inline float RNG_SFMT::draw_float_u01()
{
    if (f_count == 4)
    {
        _mm_store_ps(f,draw_m128f4_u01());
        f_count=0;
    }
    return f[f_count++];
}

inline double RNG_SFMT::draw_double_u01()
{
    if (d_count == 2)
    {
        _mm_store_pd(d,draw_m128d2_u01());
        d_count=0;
    }
    return d[d_count++];
}

inline uint32_t RNG_SFMT::draw_uint32()
{
    if (i32_count == 4)
    {
        _mm_store_ps(i32.f,draw().m128_f);
        i32_count=0;
    }
    return i32.i32[i32_count++];
}

// this version returning a scalar is a bit faster than returning a float
inline __m128 RNG_SFMT::draw_m128f1_u01()
{
    if (f_count == 4)
    {
        _mm_store_ps(f,draw_m128f4_u01());
        f_count = 0;
    }
    return _mm_load_ss(f+(f_count++));
}

inline __m128d RNG_SFMT::draw_m128d1_u01()
{
    if (d_count == 2)
    {
        _mm_store_pd(d,draw_m128d2_u01());
        d_count = 0;
    }
    return _mm_load_sd(d+(d_count++));
}

inline __m128 RNG_SFMT::draw_m128f2_uvect()
{
    __m128 ones = _mm_set1_ps(1.0);
    __m128 rnd,norm2,rnd2;
    unsigned mask;

    do {
        rnd = draw_m128f4_pm1();
        rnd2 = _mm_mul_ps(rnd,rnd);
        norm2 = _mm_hadd_ps(rnd2,rnd2);
    }
    while (!(mask=_mm_movemask_ps(_mm_cmpgt_ps(ones,norm2))&3));

    if (mask & 2)
        rnd = _mm_shuffle_ps(rnd,norm2,_MM_SHUFFLE(3,3,3,2));
    else
        rnd = _mm_shuffle_ps(rnd,norm2,_MM_SHUFFLE(2,2,1,0));

    __m128 res = _mm_div_ps(rnd,_mm_sqrt_ps(_mm_movehl_ps(rnd,rnd)));

    return _mm_and_ps(__m128(_mm_set_epi32(0,0,-1,-1)),res);
}

inline __m128 RNG_SFMT::draw_m128f3_uvect()
{
    __m128 ones = _mm_set1_ps(1.0);
    __m128 rnd,norm2,rnd2,twice_rnd;
    unsigned mask;

    do {
        rnd = draw_m128f4_pm1();
        rnd2 = _mm_mul_ps(rnd,rnd);
        norm2 = _mm_hadd_ps(rnd2,rnd2);
    }
    while(!(mask=_mm_movemask_ps(_mm_cmpgt_ps(ones,norm2))&3));

    if (mask & 2)
        rnd = _mm_shuffle_ps(rnd,norm2,_MM_SHUFFLE(3,3,3,2));
    else
        rnd = _mm_shuffle_ps(rnd,norm2,_MM_SHUFFLE(2,2,1,0));

    // now rnd  = x1  x2  (x1^2 + x2^2)  (x1^2 + x2^2)

    twice_rnd = _mm_add_ps(rnd,rnd);
    // and 2rnd = 2(x1)  2(x2)  2(x1^2 + x2^2)  2(x1^2 + x2^2)

    __m128 one_minus = _mm_movehl_ps(_mm_setzero_ps(),_mm_sqrt_ps(_mm_sub_ps(ones,rnd)));
    // one_minus = sqrt(1-x1^2-x2^2) sqrt(1-x1^2-x2^2)  0  0

    __m128 xy = _mm_mul_ps(one_minus,twice_rnd);

    // xy = 2(x1)sqrt(1-x1^2-x2^2) 2(x2)sqrt(1-x1^2-x2^2) 0 0

    return _mm_shuffle_ps(xy,_mm_sub_ss(_mm_set_ss(1.0),_mm_movehl_ps(twice_rnd,twice_rnd)),_MM_SHUFFLE(3,0,1,0));
}

#endif
