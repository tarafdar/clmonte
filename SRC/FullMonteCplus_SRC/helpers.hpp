/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef HELPERS_INCLUDED
#define HELPERS_INCLUDED
#include <cmath>  
#include <vector>
#include <cstdlib>
#include <cassert>
#include <limits>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <emmintrin.h>

using namespace std;

bool hasSuffix(string str,string sfx);
string getExtension(string str);

string bigIntSuffix(unsigned long long i_);

class Stats {
    long N,Nz;
    double sum_x,sum_xx,min_x,max_x;
    public:

    Stats() : N(0),Nz(0),sum_x(0.0),sum_xx(0.0),min_x(std::numeric_limits<double>::infinity()),max_x(-std::numeric_limits<double>::infinity()){}

    void operator()(double x) {
        ++N;
        Nz += (x==0.0);
        sum_x  += x;
        sum_xx += x*x;
        min_x = min(x,min_x);
        max_x = max(x,max_x);
    }

    long getN()      const { return N; }
    long getNz()     const { return Nz; }
    long getNnz()    const { return N-Nz; }
    double getMean() const { return sum_x/(double)N; }
    double getVar()  const { double mu=getMean(); return sum_xx/(double)N - mu*mu; }
    double getStd()  const { return sqrt(getVar()); }
    double getMin()  const { return min_x; }
    double getMax()  const { return max_x; }
    double getCV()   const { return getStd()/getMean(); }
};

ostream& operator<<(ostream&,const Stats&);

double rand01();
__m128 rand_m128f3();

template<class T> inline T clipRange(T l,T u,T x)
{
    return (x > u ? u : (x < l ? l : x));
}

template<class T>class Tolerance;
template<class T>bool operator==(T a,Tolerance<T> b);

template<class T>class Tolerance {
    T x,eps;
    public:

    Tolerance(T x_,T eps_) : x(x_),eps(eps_){}

    bool operator==(T a) const;

    bool operator()(T a,T b) const;

    // I don't think the line below is quite correct, but it's the only way it seems to work
    template<class U>friend bool operator==(U a,Tolerance<U> b);
};

template<class T>bool operator==(T a,Tolerance<T> b){ return b==a; }

#endif
