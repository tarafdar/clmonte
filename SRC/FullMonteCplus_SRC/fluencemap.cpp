/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include <iostream>
#include "linefile.hpp"
#include "newgeom.hpp"
#include "fluencemap.hpp"

using namespace std;

// define data types for the database (1=surface, 2=volume)
template<>const int FluenceMap<FaceByPointID>::dtype=1;
template<>const int FluenceMap<TetraByPointID>::dtype=2;

template<>const int PointIDLookup<FaceByPointID>::dtype=1;
template<>const int PointIDLookup<TetraByPointID>::dtype=2;

// fluence for an area element is just energy divided by area
template<>double GetFluence<FaceByPointID>::operator()(const type& p) const
{
    return p.second/mesh.getFaceArea(p.first);
}

// fluence for a volume element is energy (J) / volume(cm3) / mu_a (1/cm) = J/cm2
template<>double GetFluence<TetraByPointID>::operator()(const type& p) const
{
    unsigned IDt=p.first, IDm=mesh.getMaterial(IDt);
    return p.second/mesh.getTetraVolume(IDt)/mats[IDm].getMuA();
}

// copies from a vector, squeezing out the zeros
FluenceMapBase::FluenceMapBase(const vector<double>& d)
{
    unsigned i=0;
    for(vector<double>::const_iterator it=d.begin(); it != d.end(); ++it,++i)
        if(*it != 0.0)
            F.insert(make_pair(i,*it));
}

template<>double AreaMult<TetraByPointID>::getArea(unsigned IDt) const { return mesh.getTetraVolume(IDt); }
template<>double  AreaMult<FaceByPointID>::getArea(unsigned IDf) const { return mesh.getFaceArea(IDf); }

template<>FaceByPointID  PointIDLookup<FaceByPointID>::operator() (unsigned IDf) const { return mesh.getFacePointIDs(IDf);  }
template<>TetraByPointID PointIDLookup<TetraByPointID>::operator()(unsigned IDt) const { return mesh.getTetraPointIDs(IDt); }

template<>unsigned PointIDLookup<FaceByPointID>::operator()(const FaceByPointID& IDps) const
{
    int IDf=mesh.getFaceID(IDps.getSort());
    assert(IDf);
    if (IDf<0)
        return -IDf;
    return IDf;
}

template<>unsigned PointIDLookup<TetraByPointID>::operator()(const TetraByPointID& IDps) const
{
    unsigned IDt=mesh.getTetraID(IDps.getSort());
    assert(IDt);
    return IDt;
}

FluenceMapBase& FluenceMapBase::operator*=(double k)
{
    for(iterator it=F.begin(); it != F.end(); ++it)
        it->second *= k;
    return *this;
}

FluenceMapBase& FluenceMapBase::operator-=(const FluenceMapBase& m)
{
    for(const_iterator it=m.begin(); it != m.end(); ++it)
        operator[](it->first) -= it->second;
    return *this;
}

FluenceMapBase& FluenceMapBase::operator+=(const FluenceMapBase& m)
{
    for(const_iterator it=m.begin(); it != m.end(); ++it)
        operator[](it->first) += it->second;
    return *this;
}

// computes a/b (this / RHS)
FluenceMapBase& FluenceMapBase::operator/=(const FluenceMapBase& m)
{
    const_iterator b=m.begin();

    for (iterator a=F.begin(); a != F.end(); ++a)
    {
        if (b == m.end() || b->first > a->first)
            a->second = std::numeric_limits<double>::infinity();
        else
        {
            for(; b->first < a->first; ++b)
                F.insert(a,make_pair(b->first,0));
            if (a->first == b->first)
            {
                a->second /= b->second;
                ++b;
            }
        }
    }
    return *this;
}

void FluenceMapBase::absdiff()
{
    for(iterator it=F.begin(); it != F.end(); ++it)
        it->second = abs(it->second-1);
}

// serialize the fluence map in binary
Blob FluenceMapBase::toBinary() const
{
    uint8_t* p = new uint8_t[12*F.size()];
    uint8_t* t=p;

    const uint8_t* tmp;

    for(map<unsigned,double>::const_iterator it=F.begin(); it != F.end(); ++it,t += 12)
    {
        *(uint32_t*)t   = it->first;
        tmp = (const uint8_t*)&it->second;
        copy (tmp,tmp+8,t+4);
    }
    return Blob(12*F.size(),p);
}

// convert back from binary format
bool FluenceMapBase::fromBinary(const Blob& blob) 
{
    unsigned N=blob.getSize();
    const uint8_t *p=blob.getPtr();

    double tmp;

    if (N % 12 != 0)
        throw InvalidBlobSize();

    for(unsigned i=0; i<N/12; ++i, p += 12)
    {
        copy(p+4,p+12,(uint8_t*)&tmp);
        F.insert(make_pair(*(uint32_t*)p,tmp));
    }

    return true;
}

// Writes output in TIM-OS format, except that it excludes faces with no fluence
template<class T>void FluenceMap<T>::writeTIMOS(const TetraMesh& mesh,string fn)
{
    ofstream os(fn.c_str());

    os << 1 << ' ' << F.size() << ' ' << 1 << endl;

    for(const_IDp_fluence_iterator it=fluenceByIDpBegin(); it != fluenceByIDpEnd(); ++it)
    {
        pair<T,double> p=*it;
        os << plainwhite << p.first << mesh.getFaceArea(p.first) << ' ' << p.second << endl;
    }
}

Blob HitMap::toBinary() const
{
	Blob b(size()*12);	// 12= uint32 (4) + uint64 (8)
	uint8_t* p=b.getWritePtr();
	for(const_iterator it=begin(); it != end(); ++it)
	{
		*(uint32_t*)p 		= it->first;
		*(uint64_t*)(p+4)	= it->second;
		p+= 12;
	}
	return b;
}

void HitMap::fromBinary(const Blob& b)
{
	iterator it;
	clear();
	if (b.getSize() == 0)
		cerr << "Error, zero size for hit map!" << endl;
	else if (b.getSize() % 12 != 0)
		cerr << "Error, size not a multiple of 12 for hit map!" << endl;
	else
		for(const uint8_t* p=b.getPtr(); p<b.getEndPtr(); p+=12)
			it = insert(it,make_pair(*(const uint32_t*)p,*(const uint64_t*)(p+4)));
}

// converts a fluence map of (ID,value) into a vector where v[ID]=value
vector<double> FluenceMapBase::toVector(unsigned N) const
{
    vector<double> v(N+1,0);

    for(const_iterator it=begin(); it != end(); ++it)
        if (it->first > 0 && it->first <= N)
            v[it->first] = it->second;
        else
            cerr << "Error in map -> vector conversion; index " << it->first << " exceeds max (" << N << ")" << endl;
    return v;
}
