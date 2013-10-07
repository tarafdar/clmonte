/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef IO_TIMOS_INCLUDED
#define IO_TIMOS_INCLUDED
#include "source.hpp"
#include "graph.hpp"

#include "linefile.hpp"
#include "fluencemap.hpp"

// io_timos.hpp; Contains function prototypes for reading from TIM-OS format files
//  see https://sites.google.com/a/imaging.sbes.vt.edu/tim-os/usage for TIM-OS details

// load a surface with time-varying fluence values; sum across all time steps
template<class T>double readTIMOSOutputFlatten(const TetraMesh& mesh,FluenceMap<T>& F,LineFile& is,unsigned Nd,unsigned Ns);

void readTIMOSOutput(string fn,const TetraMesh&,SurfaceFluenceMap& surf,VolumeFluenceMap& vol);

vector<Source*> readTIMOSSource(string fn);
vector<Material> readTIMOSMaterials(string fn);
bool writeTIMOSMaterials(string fn,const vector<Material>& mat);

bool writeTIMOSSource(string fn,const vector<Source*>& s,long long packetcount=0);

#endif
