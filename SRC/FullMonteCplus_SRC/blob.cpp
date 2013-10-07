/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "blob.hpp"
#include <fstream>
#include <sys/stat.h> // for stat/creat
#include <fcntl.h>
#include <sys/types.h>
#include <openssl/sha.h>

using namespace std;

Blob::Blob(string fn,bool append_null)
{
    struct stat fileStats;
    if(stat(fn.c_str(),&fileStats))
        throw string("Failed to stat() "+fn);
    Nb = fileStats.st_size + (append_null ? 1 : 0);
    p = boost::shared_array<uint8_t>(new uint8_t[Nb]);
    ifstream is(fn.c_str());
    if(!is.good())
        throw string("Failed to open "+fn);
    is.read((char *)p.get(),fileStats.st_size);
    if(is.fail())
        throw string("Failed to read "+fn);
    if(append_null)
        p[fileStats.st_size]=0;
}

void Blob::writeFile(std::string fn) const
{
    ofstream os(fn.c_str());
    if(!os.good())
        throw std::string("Failed to open file for write");

    os.write((const char*)p.get(),Nb);

    if(!os.good())
        throw string("Failed to write");
}

void Blob::sha1_160(uint8_t* md) const
{
	SHA1(p.get(),Nb,md);
}

std::string Blob::sha1_160() const
{
    unsigned char md[20];
    char str[41];
    SHA1(p.get(),Nb,md);
    for(unsigned i=0;i<20;++i)
        sprintf(str+(i<<1),"%02x",(unsigned)md[i]);
    return std::string(str);
}

void Blob::release()
{
	Nb=0;
	p=boost::shared_array<uint8_t>(NULL);
}

