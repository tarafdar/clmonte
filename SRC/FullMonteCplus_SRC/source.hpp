/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef SOURCE_INCLUDED
#define SOURCE_INCLUDED
#include "graph.hpp"
#include "newgeom.hpp"

#include <boost/iterator/transform_iterator.hpp>
#include <boost/random/discrete_distribution.hpp>

#include "sse.hpp"

// Source classes (isotropic point, directed face, isotropic volume, directed point (pencil beam)
// A collection of sources can also be a Source (SourceMulti)
// Framework should be relatively simple to introduce new types

typedef __m128 UnitVectorType;

typedef RNG_SFMT RNG_Type;

class Source {
    double w;

    protected:

    public:
    Source(double w_=1.0) : w(w_){}
    double getPower() const { return w; }

    virtual bool prepare(const TetraMesh& m){ return true; }

    virtual pair<Packet,unsigned>            emit(RNG_Type&) const=0;

    virtual string operator()() const=0;
    virtual string timos_str(unsigned long long=0) const=0;
	virtual ostream& print(ostream& os)=0;
};

class PointSource : virtual public Source {
    pair<Point<3,double>, unsigned> origin;

    protected:
    virtual pair<Point<3,double>,unsigned> getOrigin(RNG_Type&) const { return origin; }

    public:
    PointSource(const Point<3,double>& p_,double w_=1.0) : Source(w_),origin(make_pair(p_,0)){}

    virtual bool prepare(const TetraMesh& m);

    Point<3,double> getOrigin() const { return origin.first; }
    unsigned getTetraID() const { return origin.second; }
};

class IsotropicSource : virtual public Source {
    protected:
    virtual UnitVectorType getDirection(RNG_Type&) const;

    public:
    IsotropicSource(double w_=1.0) : Source(w_){}
};

class IsotropicPointSource : public IsotropicSource, public PointSource {
    public:
    IsotropicPointSource(Point<3,double>& p_,double w_=1.0) : PointSource(p_,w_){};

    virtual bool prepare(const TetraMesh& m){ return IsotropicSource::prepare(m) & PointSource::prepare(m); }
    virtual string operator()() const { return "Isotropic point source"; }
    virtual string timos_str(unsigned long long=0) const;

    virtual pair<Packet,unsigned> emit(RNG_Type&) const;
	virtual ostream& print(ostream&);
};

class PencilBeamSource : public PointSource {
    Packet pkt;
    int IDt,IDf;
    public:

    PencilBeamSource(Ray<3,double> r_,double w=1.0) : PointSource(r_.getOrigin(),w),pkt(r_){ }
    PencilBeamSource(Point<3,double> p_,UnitVectorType d_,double w=1.0,int IDt_=0) : PointSource(p_,w),
        pkt(Ray<3,double>(p_,uvectFrom(d_))),IDt(IDt_){ }

    virtual string operator()() const { return "Pencil beam source"; }

    virtual pair<Packet,unsigned> emit(RNG_Type&) const;

    virtual bool prepare(const TetraMesh& m);

    virtual string timos_str(unsigned long long=0) const;

	virtual ostream& print(ostream&);

    unsigned getIDt() const { return IDt; }
};

ostream& operator<<(ostream& os,Source& src);

class VolumeSource : public IsotropicSource {
    // creates a random tetrahedral source by shearing the unit cube
    double M[3][3];
    double P0[3];
    unsigned IDt;

    protected:
    virtual pair<Point<3,double>,unsigned> getOrigin(RNG_Type&) const;

    public:
    VolumeSource(unsigned IDt_=0,double w=1.0) : Source(w),IDt(IDt_){}
    virtual string timos_str(unsigned long long=0) const;
    string operator()() const { return "Volume source"; }

    bool prepare(const TetraMesh&);

    virtual pair<Packet,unsigned> emit(RNG_Type& rng) const
        {
            Packet pkt;
            pkt.setDirection(rng.draw_m128f3_uvect());
            pkt.p = getOrigin(rng).first;
            return make_pair(pkt,IDt);
        }

    virtual ostream& print(ostream&);
    unsigned getIDt() const { return IDt; }
};

class FaceSource : public Source {
    Packet pkt;
    FaceByPointID f;
    unsigned IDt;
    int IDf;
    UnitVectorType n;
    Point<3,double> P0;
    double M[3][2];

    protected:
    virtual pair<Point<3,double>,unsigned> getOrigin(RNG_Type&) const;

    public:

    // if force_boundary is set, requires the face to be pointing in from an object boundary
    virtual bool prepare(const TetraMesh& m) { return prepare(m,true); }
    bool prepare(const TetraMesh&,bool force_boundary=true);

    double w;
    public:
    FaceSource(FaceByPointID f_,double w_=1.0) : Source(w_),f(f_),IDt(0),IDf(0){}
    virtual string operator()() const { return "Face Source"; }

    virtual pair<Packet,unsigned> emit(RNG_Type&) const;
	virtual ostream& print(ostream& os);
    virtual string timos_str(unsigned long long=0) const;

    unsigned getIDt() const { return IDt; }
    FaceByPointID getIDps() const { return f; }
};

class SourceMulti : public Source {
    vector<Source*> sources;
    boost::random::discrete_distribution<unsigned> source_dist;
    double w_total;

    static double _f_getPower(Source* s){ return s->getPower(); }

    public:
    // Create from a pair of iterators that dereference to a Source*
    template<class ConstIterator> SourceMulti(ConstIterator begin,ConstIterator end) :
        sources(begin,end),
        source_dist(boost::make_transform_iterator(begin,_f_getPower),
            boost::make_transform_iterator(end,_f_getPower)),
        w_total(0.0)
        { for(; begin != end; ++begin) w_total += (*begin)->getPower(); }

    virtual string operator()() const { return "Multiple sources"; }
    virtual string timos_str(unsigned long long=0) const;

    virtual bool prepare(const TetraMesh& m);

    virtual pair<Packet,unsigned> emit(RNG_Type&) const;

	virtual ostream& print(ostream& os);
};

#endif
