/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "LoggerSurface.hpp"

void writeFileBin(string,const vector<Point<3,double> >&,const map<FaceByPointID,double>&);

void LoggerSurface::fluenceMap(SurfaceFluenceMap& F)
{
    F.clear();

    cout << "Creating fluence map; counts.size() = " << counts.size() << endl;

    // copy from vector to map, eliminating zero entries and dividing by face area
    for(TetraMesh::boundary_f_const_iterator it=mesh.boundaryFaceBegin(); it != mesh.boundaryFaceEnd(); ++it)
    {
        unsigned IDf = it->first;
        if (std::isnan(getValue(counts[IDf])))
            cerr << "NaN reported on surface for element " << IDf << " (max val " << mesh.getNf() << "; surface size " << counts.size() << ")" << endl;
        if (getValue(counts[IDf]) != 0)
            F[IDf] = getValue(counts[IDf]);
    }
}

void LoggerSurface::hitMap(map<unsigned,unsigned long long>& m)
{
    m.clear();
    map<unsigned,unsigned long long>::iterator m_it=m.begin();
    unsigned i=1;
    if(counts.size()<2)
        return;
    for(vector<FluenceCountType>::const_iterator it=counts.begin()+1; it != counts.end(); ++it,++i)
        if (getHits(*it) != 0)
            m_it = m.insert(m_it,make_pair(i,getHits(*it)));
}

/*void LoggerSurface::resultMap(map<FaceByPointID,double>& m,bool per_area)
{
    m.clear();

    // copy from vector to map
    for(TetraMesh::boundary_f_const_iterator it=mesh.boundaryFaceBegin(); it != mesh.boundaryFaceEnd(); ++it)
    {
        if (getValue(counts[it->first]) != 0)
            m.insert(make_pair(mesh.getFacePointIDs(it->first),getValue(counts[it->first])));
    }

    for(map<FaceByPointID,double>::iterator it=m.begin(); per_area && it != m.end(); ++it)
    {
        const Point<3,double> &A = mesh.getPoint(it->first[0]),
            &B = mesh.getPoint(it->first[1]),
            &C = mesh.getPoint(it->first[2]);
        Vector<3,double> AB(A,B), AC(A,C);
        double area=cross(AB,AC).norm_l2()/2;
        assert(area != 0.0);
        it->second /= area;
    }
}*/

// Writes out a VTK datafile (v3.0)
/*void LoggerSurface::writeFileVTK(string fn)
{
    map<FaceByPointID,double> m;

    for(TetraMesh::boundary_f_const_iterator it=mesh.boundaryFaceBegin(); it != mesh.boundaryFaceEnd(); ++it)
    {
        FaceByPointID f = *it;

        const Point<3,double> &A = mesh.P[mesh.F_p[it->first][0]], &B = mesh.P[mesh.F_p[it->first][1]], &C = mesh.P[mesh.F_p[it->first][2]];
        Vector<3,double> AB(A,B), AC(A,C);
        double area=cross(AB,AC).norm_l2();
        assert(area != 0.0);
        if(area != 0.0)
            m.insert(make_pair(mesh.F_p[it->first],counts[it->first].second/area));
        else
            m.insert(make_pair(*it,0));
    }

    ::writeFileVTK(fn,mesh.P,m);
}*/
