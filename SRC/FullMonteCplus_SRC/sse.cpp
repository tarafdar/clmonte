/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "sse.hpp"
#include <iomanip>

std::ostream& operator<<(std::ostream& os,__m128 x)
{
    float __attribute__((aligned(16))) f[4];
    _mm_store_ps(f,x);
    os << std::fixed << std::setprecision(6) << '<' << f[0] << ',' << f[1] << ',' << f[2] << ',' << f[3] << '>';
    return os;
}

std::ostream& operator<<(std::ostream& os,__m128d x)
{
	double __attribute__((aligned(16))) d[2];
	_mm_store_pd(d,x);
	return os << std::fixed << std::setprecision(3) << '<' << d[0] << ',' << d[1] << '>';
}

std::ostream& operator<<(std::ostream& os,__m128i x)
{
	unsigned char __attribute__((aligned(16))) cx[16];
	_mm_store_si128((__m128i*)cx,x);
	for(int i=15;i >= 0;--i)
		os << std::setw(2) << std::hex << std::setfill('0') << (unsigned)cx[i];
	return os;
}
