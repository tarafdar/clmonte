/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "source.hpp"
#include "newgeom.hpp"
#include "graph.hpp"
#include <limits>
#include <sstream>
#include "random.hpp"

string PencilBeamSource::timos_str(unsigned long long Np) const
{
    stringstream ss;
    ss << "11 " << IDt << ' ' << plainwhite << pointFrom(pkt.p) << ' ' << plainwhite << uvectFrom(pkt.d) << ' ' << (Np == 0 ? (unsigned long long)getPower() : Np);
    return ss.str();
}

string IsotropicPointSource::timos_str(unsigned long long Np) const
{
    stringstream ss;
    ss << "1 " << plainwhite << getOrigin() << ' ' << (Np == 0 ? (unsigned long long)getPower() : Np);
    return ss.str();
}

string VolumeSource::timos_str(unsigned long long Np) const
{
    stringstream ss;
    ss << "2 " << IDt << ' ' << (Np == 0 ? (unsigned long long) getPower() : Np);
    return ss.str();
}

string FaceSource::timos_str(unsigned long long Np) const
{
    stringstream ss;
    ss << "12 " << f[0] << ' ' << f[1] << ' ' << f[2] << ' ' << (Np == 0 ? (unsigned long long) getPower() : Np);
    return ss.str();
}

bool PencilBeamSource::prepare(const TetraMesh& m)
{
    Point<3,double> Q;
    double t,t_min=numeric_limits<double>::infinity();

    unsigned i=1;

    // Iterate through all points
    for(TetraMesh::face_id_const_iterator it=m.faceIDBegin(); it != m.faceIDEnd(); ++it,++i)
    {
        Point<3,double> tmp[3] = { m.getPoint((*it)[0]), m.getPoint((*it)[1]), m.getPoint((*it)[2]) };

        if (PointInTriangle(pointFrom(pkt.p),uvectFrom(pkt.d),tmp,Q,t)){
            if (t < t_min)
            {
                t_min=t;
                IDf=i;
                if (dot(m.getFace(IDf).getNormal(),uvectFrom(pkt.d)) < 0)
                    IDf *= -1;
                IDt=m.getTetraIDFromFaceID(IDf);
            }
        }
    }

    return true;
}

bool VolumeSource::prepare(const TetraMesh& m)
{
    Point<3,double> A=m.getTetraPoint(IDt,0);
    Point<3,double> B=m.getTetraPoint(IDt,1);
    Point<3,double> C=m.getTetraPoint(IDt,2);
    Point<3,double> D=m.getTetraPoint(IDt,3);

    // copy first point
    P0[0]=A[0];
    P0[1]=A[1];
    P0[2]=A[2];

    // copy AB, AC, AD
    M[0][0]=B[0]-P0[0];
    M[1][0]=B[1]-P0[1];
    M[2][0]=B[2]-P0[2];

    M[0][1]=C[0]-P0[0];
    M[1][1]=C[1]-P0[1];
    M[2][1]=C[2]-P0[2];

    M[0][2]=D[0]-P0[0];
    M[1][2]=D[1]-P0[1];
    M[2][2]=D[2]-P0[2];
    return true;
}

bool FaceSource::prepare(const TetraMesh& m,bool force_boundary)
{
    // find face in TetraMesh
    IDf = m.getFaceID(f);

    if(force_boundary)
    {
        unsigned tetraDown = m.getTetraFromFace(-IDf), tetraUp = m.getTetraFromFace(IDf);
        if (tetraUp){
            IDt=tetraUp;
            if (tetraDown)
                cerr << "Surprise: source is not on a boundary" << endl;
            else
                IDf = -IDf;
        }
        else if (tetraDown)
            IDt=tetraDown;
        else
            cerr << "Surprise: face borders no tetras!" << endl;
    }
    else 
        IDt = m.getTetraFromFace(IDf);

    // get normal
    n = -m.getFace(IDf).getNormal();

    // get anchor point (A)
    P0 = m.getPoint(f[0]);

    const Point<3,double> &B=m.getPoint(f[1]), &C=m.getPoint(f[2]);

    // create shear matrix
    M[0][0] = B[0]-P0[0];
    M[0][1] = C[0]-P0[0];
    M[1][0] = B[1]-P0[1];
    M[1][1] = C[1]-P0[1];
    M[2][0] = B[2]-P0[2];
    M[2][1] = C[2]-P0[2];

    pkt.setDirection(n);
    return true;
}

pair<Point<3,double>,unsigned > FaceSource::getOrigin(RNG_Type& rng) const
{
    double s=rng.draw_float_u01(),t=rng.draw_float_u01();
    if (s+t > 1)
        s=1-s,t=1-t;

    double P[4] = {
        P0[0] + M[0][0]*s + M[0][1]*t,
        P0[1] + M[1][0]*s + M[1][1]*t,
        P0[2] + M[2][0]*s + M[2][1]*t,
        0
    };
    return make_pair(Point<3,double>(P),IDt);
}

pair<Packet,unsigned> FaceSource::emit(RNG_Type& rng) const
{
    Packet p=pkt;
    p.p = getOrigin(rng).first;
    return make_pair(p,IDt);
}

ostream& IsotropicPointSource::print(ostream& os)
{
	return os << "Isotropic point source located " << getOrigin() << " with weight " << getPower();
}

ostream& PencilBeamSource::print(ostream& os)
{
	return os << "Pencil beam source located " << pointFrom(pkt.p) << " aimed " << uvectFrom(pkt.d) << " entering at IDf=" << IDf << " IDt=" << IDt << " with weight " << getPower();
}

ostream& VolumeSource::print(ostream& os)
{
    return os << "Volume source located in tetrahedron IDt=" << IDt << " with weight " << getPower();
}

ostream& operator<<(ostream& os,Source& src)
{
	return src.print(os);
}

pair<Packet,unsigned> IsotropicPointSource::emit(RNG_Type& rng) const
{
    pair<Packet,unsigned> tmp;
    Packet& pkt=tmp.first;
    tmp.second=getTetraID();

    pkt.p = getOrigin();
    pkt.setDirection(rng.draw_m128f3_uvect());

    return tmp;
}

UnitVectorType IsotropicSource::getDirection(RNG_Type& rng) const
{
    return rng.draw_m128f3_uvect();
}

bool PointSource::prepare(const TetraMesh& m)
{
    origin.second=m.findEnclosingTetra(origin.first);
    return (origin.second != 0);
}

pair<Point<3,double>,unsigned> VolumeSource::getOrigin(RNG_Type& rng) const
{
    double p[3],src[3];
    double rnd[3] = { rng.draw_double_u01(), rng.draw_double_u01(), rng.draw_double_u01() };
    if (rnd[0]+rnd[1] > 1.0)
    {
        rnd[0] = 1.0-rnd[0];
        rnd[1] = 1.0-rnd[1];
    }

    if (rnd[0]+rnd[1]+rnd[2] > 1.0)
    {
        if(rnd[1]+rnd[2] < 1.0)
        {
            p[0] = 1-rnd[1]-rnd[2];
            p[1] = rnd[1];
            p[2] = rnd[0]+rnd[1]+rnd[2] - 1;
        }
        else {
            p[0] = rnd[0];
            p[1] = 1-rnd[2];
            p[2] = 1-rnd[0]-rnd[1];
        }
    }
    else {
        p[0]=rnd[0];
        p[1]=rnd[1];
        p[2]=rnd[2];
    }

    assert(p[0] >= 0.0 && p[1] >= 0.0 && p[2] >= 0.0 && p[0]+p[1]+p[2] <= 1.0);

    src[0] = P0[0] + M[0][0]*p[0] + M[0][1]*p[1] + M[0][2]*p[2];
    src[1] = P0[1] + M[1][0]*p[0] + M[1][1]*p[1] + M[1][2]*p[2];
    src[2] = P0[2] + M[2][0]*p[0] + M[2][1]*p[1] + M[2][2]*p[2];

    return make_pair(Point<3,double>(src),IDt);
}

ostream& FaceSource::print(ostream& os)
{
    return os << "Face source, IDf=" << IDf << " Points " << f << " Direction " << n;
}


bool SourceMulti::prepare(const TetraMesh& m)
{
    bool result=true;
    for(vector<Source*>::const_iterator it=sources.begin(); it != sources.end(); ++it)
        result &= (*it)->prepare(m);
    return result;
}

string SourceMulti::timos_str(unsigned long long Npacket) const
{
    string str;
    for(vector<Source*>::const_iterator it=sources.begin(); it != sources.end(); ++it)
        str += (*it)->timos_str() + "\n";
    return str;
}

// select a source from the distribution and call its emit() member
pair<Packet,unsigned> SourceMulti::emit(RNG_Type& rng) const
{
    return sources[source_dist(rng)]->emit(rng);
}

pair<Packet,unsigned> PencilBeamSource::emit(RNG_Type& rng) const
{
    return make_pair(pkt,IDt);
}

ostream& SourceMulti::print(ostream& os)
{
    for(vector<Source*>::const_iterator it=sources.begin(); it != sources.end(); ++it)
        (*it)->print(os) << endl;
    return os;
}
