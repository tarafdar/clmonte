/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef FLUENCEMAP_INCLUDED
#define FLUENCEMAP_INCLUDED
#include <inttypes.h>
#include <map>

#include "newgeom.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "graph.hpp"

#include "blob.hpp"

using namespace std;

template<class T>class FluenceMap;

template<class T>class PointIDLookup {
    protected:
    const TetraMesh& mesh;
    public:
    static const int dtype;

    typedef pair<unsigned,double> type;
    typedef pair<T,double> result_type;

    PointIDLookup(const PointIDLookup& lu_) : mesh(lu_.mesh){}
    PointIDLookup(const TetraMesh* mesh_)   : mesh(*mesh_){ if (!mesh_) throw typename FluenceMap<T>::InvalidMesh(); }

    // looks up the IDps (FaceByPointID/TetraByPointID) type, returning an unsigned ID
    unsigned operator()(const T&) const;

    pair<T,double> operator()(const pair<unsigned,double>& p) const { return make_pair(operator()(p.first),p.second); }
    T operator()(unsigned) const;
};

template<class T>class GetFluence : public PointIDLookup<T> {
    using PointIDLookup<T>::mesh;
    const vector<Material>& mats;
    public:
    typedef double result_type;
    typedef pair<unsigned,double> type;

    GetFluence(const GetFluence& gf_,const vector<Material>& mats_)  : PointIDLookup<T>(&gf_.mesh),mats(mats_){}
    GetFluence(const TetraMesh* mesh_,const vector<Material>& mats_) : PointIDLookup<T>(mesh_),mats(mats_){}

    result_type operator()(const type& p) const;
};

template<>double GetFluence<FaceByPointID>::operator()(const type& p) const;
template<>double GetFluence<TetraByPointID>::operator()(const type& p) const;


template<class T>class AreaMult : public PointIDLookup<T> {
    using PointIDLookup<T>::mesh;
    double getArea(unsigned) const;
    public:
    typedef pair<unsigned,double> result_type;
    typedef pair<unsigned,double> type;
    AreaMult(const AreaMult& am_)    : PointIDLookup<T>(&am_.mesh){}
    AreaMult(const TetraMesh* mesh_) : PointIDLookup<T>(mesh_){}
    double operator()(unsigned ID) const { return getArea(ID); }
    pair<unsigned,double> operator()(const pair<unsigned,double>& p) const
        { return make_pair(p.first,p.second*getArea(p.first)); }
};


template<class A,class B>const B& getFirst(const pair<A,B>&  p){ return p.first;  }
template<class A,class B>const B& getSecond(const pair<A,B>& p){ return p.second; }

class FluenceMapBase {
    protected:
    // Maps from the unsigned Face ID (IDf) to a double-precision fluence value
    map<unsigned,double> F;

    // provides a mapping from FaceByPointID (sorted) to an unsigned faceID
    const TetraMesh* mesh;

    public:
    class InvalidMesh {};
    class InvalidBlobSize {};

    FluenceMapBase(const TetraMesh* mesh_=NULL) : mesh(mesh_){}
    FluenceMapBase(const vector<double>& d);
    FluenceMapBase(const map<unsigned,double>& F_,const TetraMesh* mesh_) : F(F_),mesh(mesh_){}

    void clear(){ F.clear(); }
    unsigned size(){ return F.size(); }
    void setMesh(const TetraMesh* mesh_){ mesh=mesh_; }

    // dereferences to a pair<unsigned,double>
    typedef map<unsigned,double>::const_iterator    const_iterator;
    typedef map<unsigned,double>::iterator          iterator;

    iterator insert(iterator it,pair<unsigned,double> p){ return F.insert(it,p); }

    iterator begin() { return F.begin(); }
    iterator   end() { return F.end();   }

    const_iterator begin() const { return F.begin(); }
    const_iterator end()   const { return F.end(); }

    double& operator[](unsigned IDf)    { return F[IDf]; }
    double& operator[](int IDf)         { return F[unsigned(abs(IDf))]; }

    // iterator returning only the fluence values, dereferences to a double
    typedef boost::transform_iterator<const double&(*)(const pair<unsigned,double>&),const_iterator> const_value_iterator;
    const_value_iterator valuesBegin() const { return boost::make_transform_iterator(F.begin(),&getSecond<unsigned,double>); }
    const_value_iterator valuesEnd()   const { return boost::make_transform_iterator(F.end(),  &getSecond<unsigned,double>); }

    // Add/Subtract another fluence map
    FluenceMapBase& operator-=(const FluenceMapBase& m);
    FluenceMapBase& operator+=(const FluenceMapBase& m);
    FluenceMapBase& operator/=(const FluenceMapBase& m);

    FluenceMapBase& operator*=(double);

    void absdiff();

    void writeASCII(string fn);

    // convert to vector
    vector<double> toVector(unsigned N) const;

    // Serialization to/from binary format
    Blob    toBinary()              const;
    bool    fromBinary(const Blob&);
};

template<class T>class FluenceMap : public FluenceMapBase {
    public:
    static const int dtype;

    using FluenceMapBase::operator[];
	FluenceMap(const TetraMesh& mesh_,const Blob& b_) : FluenceMapBase(&mesh_){ fromBinary(b_); }

    FluenceMap(const TetraMesh* mesh_) : FluenceMapBase(mesh_){}

    typedef T PointIDType;
    typedef PointIDLookup<T> PointIDLookupType;

//    typedef boost::transform_iterator<AreaMult<T>,const_iterator> const_energy_iterator;
    typedef map<unsigned,double>::const_iterator const_energy_iterator;

    // returns a pair<PointIDType,double> with the fluence values
    typedef boost::transform_iterator<PointIDLookupType,const_iterator> const_IDp_fluence_iterator;

    const_IDp_fluence_iterator fluenceByIDpBegin() const
        { return boost::make_transform_iterator(F.begin(),PointIDLookupType(mesh)); }
    const_IDp_fluence_iterator   fluenceByIDpEnd() const
        { return boost::make_transform_iterator(F.end(),PointIDLookupType(mesh)); }
    
    // returns an iterator for the energy (fluence*area or fluence*volume)
    const_energy_iterator energyBegin() const
//        { return boost::make_transform_iterator(F.begin(),AreaMult<PointIDType>(mesh)); }
        { return F.begin(); }
    const_energy_iterator energyEnd()   const
//        { return boost::make_transform_iterator(F.end(),AreaMult<PointIDType>(mesh)); }
        { return F.end(); }

    // read/write in TIMOS-compatible form
    void    writeTIMOS(const TetraMesh&,string fn);
    void    writeText(string fn,const vector<Material>&);

    // finds the total energy
    double getTotalEnergy() const {
        double sum=0.0;
        for(const_energy_iterator it=energyBegin(); it != energyEnd(); ++it)
            sum += it->second;
        return sum;
    };

    // Element access
    double& operator[](PointIDType f)   { return F[PointIDLookupType(mesh)(f)]; }
};

class HitMap : public map<unsigned,unsigned long long> {
	public:

	Blob toBinary() const;
	void fromBinary(const Blob&);
};

typedef FluenceMap<FaceByPointID>  SurfaceFluenceMap;
typedef FluenceMap<TetraByPointID> VolumeFluenceMap;

// Writes out in a plaintext format
//  <point IDs: 3 (surface) or 4 (volume)> <element area/volume> <fluence J/cm2> <total energy J>
//
// eg for a the volume element enclosed by points 1,2,5,10 volume 0.375 
//
// For volume elements, the fluence is found by dividing the absorption coeff (in 1/cm)

template<class ElementType>void FluenceMap<ElementType>::writeText(string fn,const vector<Material>& mats)
{
    ofstream os(fn);
    PointIDLookup<ElementType> pID(mesh);
    AreaMult<ElementType> elSize(mesh);
    GetFluence<ElementType> fluence(mesh,mats);

    ElementType IDps;

    for(map<unsigned,double>::const_iterator it=begin(); it != end(); ++it)
    {
        IDps = pID(it->first);
        os << plainwhite << IDps << ' ' << elSize(it->first) << ' ' << fluence(*it) << ' '
            << it->second << endl;
    }
    os.close();
}

#endif
