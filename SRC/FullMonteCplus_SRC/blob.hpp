/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef BLOB_INCLUDED
#define BLOB_INCLUDED
#include <boost/shared_array.hpp>
#include <inttypes.h>
#include <string>
#include <cstdio>

// Blob = Binary Large OBject
//  Just a bunch of bytes as a Boost shared array, with a method provided to compute SHA1 hash

class Blob 
{
    unsigned                        Nb;
    boost::shared_array<uint8_t>    p;

    public:

    Blob() : Nb(0){}
    Blob(unsigned Nb_,bool append_null=false) : Nb(Nb_),p(new uint8_t[Nb_]){}
    Blob(unsigned Nb_,uint8_t* p_) : Nb(Nb_),p(p_){}
    Blob(const Blob& b) : Nb(b.Nb),p(b.p){}
    Blob(std::string fn,bool append_null=false);

    void sha1_160(uint8_t* md) const;
	std::string sha1_160() const;

    void writeFile(std::string fn) const;

    void release();

    unsigned        getSize() const { return Nb; };
	uint8_t*		getWritePtr() { return p.get(); }
    const uint8_t*  getPtr()  const { return p.get(); };
    const uint8_t*  getEndPtr() const { return p.get()+Nb; }
};

#endif
