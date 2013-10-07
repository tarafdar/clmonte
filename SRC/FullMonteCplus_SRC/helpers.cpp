/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include <cmath>
#include <string>
#include <algorithm>

#include "sse.hpp"
#include "helpers.hpp"

using namespace std;

template<>bool Tolerance<double>::operator==(double a) const { return fabs(a-x)<eps; }
template<>bool Tolerance<float>::operator==(float a) const { return fabs(a-x)<eps; }
template<>bool Tolerance<long double>::operator==(long double a) const { return fabs(a-x)<eps; }

template<>bool Tolerance<double>::operator()(double a,double b) const { return fabs(a-b)<eps; }

ostream& operator<<(ostream& os,const Stats& s)
{
    return os << "Range [" << s.getMin() << ',' << s.getMax() << "] mean " << s.getMean() << " std " << s.getStd() << " Nz=" << s.getNz() << endl;
}

string getExtension(string str)
{
    size_t pos=str.find_last_of('.');
    if (pos==string::npos)
        return string("");
    return str.substr(pos+1);
}

bool hasSuffix(string str,string sfx)
{
    for(string::const_reverse_iterator a=str.rbegin(),b=sfx.rbegin(); b != str.rend(); ++a,++b)
    {
        if (a==str.rend() || *a != *b)
            return false;
    }
    return true;
}

string bigIntSuffix(unsigned long long i_)
{
    unsigned j=0;
    double i=i_;
    char suffixes[] = {" kMGTP"};

    for(j=0; i >= 1000; ++j)
        i /= 1000;

    stringstream ss;
    if (j != 0)
        ss << fixed << setw(6) << setprecision(2) << i << suffixes[j];
    else
        ss << i_;
    return ss.str();
}
