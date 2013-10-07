/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include <utility>
#include <iostream>
#include <iomanip>

using namespace std;

// HitCount<C,V>
//  Accumulates a value-type V, while keeping track of how many times it has been added to using counter-type C
//  += pair<C,V>
//  += <V> (assume C=1)
//  += HitCount<C,V>
//
// Accessors getHits(), getValue()

template<class CounterType,class ValueType>class HitCount;
template<class CounterType,class ValueType>CounterType getHits(const HitCount<CounterType,ValueType>& hc);
template<class CounterType,class ValueType>ValueType getValue(const HitCount<CounterType,ValueType>& hc);

template<class CounterType,class ValueType>class HitCount {
    CounterType hits;
    ValueType   value;

    public:

    HitCount(CounterType hits_=CounterType(),const ValueType& value_=ValueType()) : hits(hits_),value(value_){}

    HitCount& operator+=(const ValueType& rhs){ hits++; value += rhs; return *this; }
    HitCount& operator+=(const HitCount& rhs){ hits += rhs.hits; value += rhs.value; return *this; }
    template<class C,class V>HitCount& operator+=(const pair<C,V>& rhs)
        { hits += rhs.first; value += rhs.second; return *this; }

	HitCount& operator=(const ValueType& rhs){ hits=1; value = rhs; return *this; }

    // returns a std::pair of (hits,value) for any types implicitly convertible from CounterType,ValueType
    template<class C,class V> operator pair<C,V>() const { return std::pair<C,V>(hits,value); }

    // Helper functions (useful for iterators for instance)
    friend CounterType getHits<>(const HitCount&);
    friend ValueType   getValue<>(const HitCount&);

    CounterType getHits() const { return hits; }
    ValueType getValue()  const { return value; }

    template<class C,class V>friend ostream& operator<<(ostream&,const HitCount<C,V>&);
};

template<class CounterType,class ValueType>CounterType getHits(const HitCount<CounterType,ValueType>& hc)
{ return hc.hits; }
template<class T>T getHits(const T&){ return 0; }

template<class CounterType,class ValueType>ValueType getValue(const HitCount<CounterType,ValueType>& hc)
{ return hc.value; }
template<class T>T getValue(const T& t){ return t; }

template<class CounterType,class ValueType>ostream& operator<<(ostream& os,const HitCount<CounterType,ValueType>& hc)
{
    return os << hc.hits << ' ' << hc.value << endl;
}
