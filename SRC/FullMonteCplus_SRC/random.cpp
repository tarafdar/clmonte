/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "random.hpp"
#include "sse.hpp"

// Class RNG uses a Boost random number generator
// Class RNG_SFMT is a wrapper for Saito & Matsumoto SFMT

double RNG::draw_double_u01(){ return uni01d(rng); }
float  RNG::draw_float_u01(){ return uni01f(rng);  }
__m128 RNG::draw_m128f3_u01(){ float f[4] = { uni01f(rng),uni01f(rng),uni01f(rng),0 }; return _mm_load_ps(f); }
__m128 RNG::draw_m128f4_u01(){ float f[4] = { uni01f(rng),uni01f(rng),uni01f(rng),uni01f(rng) }; return _mm_load_ps(f); }
__m128 RNG::draw_m128f3_uvect(){ return normalize(_mm_sub_ps(draw_m128f3_u01(),_mm_set_ps(0,0.5,0.5,0.5))); }

void RNG_SFMT::refill()
{
    nextRand = randBuf;
    sfmt_fill_array32(&sfmt,(uint32_t*)randBuf,Nbuf*4);
}
