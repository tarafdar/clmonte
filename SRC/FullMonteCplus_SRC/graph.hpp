/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef GRAPH_HPP_INCLUDE
#define GRAPH_HPP_INCLUDE
#include <emmintrin.h>
#include <xmmintrin.h>
#include <mmintrin.h>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <inttypes.h>

#include <boost/shared_array.hpp>

#include "optics.hpp"

#include "newgeom.hpp"
#include "blob.hpp"
#include "random.hpp"

using namespace std;

class Face {
	// plane definition
	UnitVector<3,double> normal;        // normal points into the tetrahedron
	double               C;

	public:

	Face() : C(std::numeric_limits<double>::quiet_NaN()) {};
    Face(const Point<3,double>&,const Point<3,double>&,const Point<3,double>&);
	Face(const Point<3,double>&,const Point<3,double>&,const Point<3,double>&,const Point<3,double>&);

	// flip the face orientation
	void flip(){ C=-C; normal=-normal; }

    // get the normal
    UnitVector<3,double> getNormal() const { return normal; }
    double getConstant() const { return C; }

	// query points
	bool   pointAbove(const Point<3,double>&)  const;
	double pointHeight(const Point<3,double>&) const;

	// Test if a ray intersects the face
	pair<bool,double>           rayIntersect(const Ray<3,double>&,double=numeric_limits<double>::infinity(),bool=false) const;

    // Gets the point where a ray intersects the face (returns false if no intersection)
	pair<bool,Point<3,double> > rayIntersectPoint(const Ray<3,double>&,bool=false) const;

    Face operator-() const { Face f; f.normal = -normal; f.C=-C; return f; }

	// get projections of a vector onto the face
	Vector<3,double> project(const Vector<3,double>&) const;
	Vector<3,double> normalComponent(const Vector<3,double>&) const;
	pair<UnitVector<2,double>,UnitVector<3,double> > reflectionBasis(const UnitVector<3,double>&,bool=false) const;

    double nx() const { return normal[0]; }
    double ny() const { return normal[1]; }
    double nz() const { return normal[2]; }
    double nC() const { return C; }

	friend ostream& operator<<(ostream&,const Face&);
};

class Material {
    static const float const_c0;    // units of mm/ns (3e8 m/s = 3e11 mm/s = 3e2 mm/ns)
    float mu_s, mu_a, mu_p, mu_t, g, n, one_minus_gg, one_plus_gg, recip_2g, albedo, absfrac; // 11x4 = 40B
    bool matchedboundary,isscattering;          // 2x1B = 2B
    public:

    __m128 s_prop,s_init;

    Material(double mu_a_=0,double mu_s_=0,double g_=0,double n_=1.0,double mu_p_=0.0,bool matchedboundary_=false) :
        mu_s(mu_s_),mu_a(mu_a_),mu_p(mu_p_),mu_t(mu_s_+mu_a_+mu_p_),
        g(g_),n(n_),one_minus_gg(1.0-g_*g_),one_plus_gg(1.0+g_*g_),recip_2g(1.0/2.0/g_),
        albedo((mu_s_+mu_p_)/(mu_a_+mu_s_+mu_p_)),
        absfrac(mu_a_/(mu_a_+mu_s_+mu_p_)),
        matchedboundary(matchedboundary_),
        isscattering(!(g_ == 1.0 || mu_s_ == 0)),
        s_prop(_mm_set_ps(0,n_/const_c0,-mu_t,-1)),
        s_init(_mm_set_ps(0,0,-1,-1/mu_t))
        {}

    bool isScattering() const { return isscattering; }

    float getParam_g() const { return g; }
    float getn()       const { return n; }
    float getMuA()     const { return mu_a; }
    float getMuS()     const { return mu_s; }
    float getMuT()     const { return mu_t; }
    float getg()       const { return g; }

    float getAlbedo()              const { return albedo; }
    float getAbsorbedFraction()    const { return absfrac; }

	Packet Scatter(Packet pkt,float rnd0,__m128 cosphi_sinphi) const;

    bool isMatched() const { return matchedboundary; }
    void setMatched(bool m=true){ matchedboundary=m; }

    friend ostream& operator<<(ostream& os,const Material& mat);
};

inline Packet Material::Scatter(Packet pkt,float rnd0,__m128 cosphi_sinphi) const
{
    float costheta,P=2*rnd0-1;

    float t=one_minus_gg/(1+g*P);

    // choose angles: HG function for component along d0, uniform circle for normal components
    costheta = g==0 ? P : recip_2g * (one_plus_gg-t*t);
    assert (costheta <= 1.0 && costheta >= -1.0);

    return matspin(pkt,costheta,cosphi_sinphi);
}

    typedef struct { 
        __m128 Pe;          // 4x16B = 64 B
        __m128 distance;
        int IDfe;           // 4B
        unsigned IDte;      // 4B
        int idx;            // 4B
    } StepResult;

struct Tetra {
    __m128 nx,ny,nz,C;      // 4 x 16B = 64B
//    __m128 Fn[4];           // 4x16B = 64B         Don't keep these here; only need occasionally -> transpose as needed
    TetraByFaceID   IDfs;   // 4 x 4B = 16 B
    unsigned adjTetras[4];  // 4 x 4B = 16 B
    unsigned matID;         // 4 B

    bool pointWithin(__m128);

    StepResult getIntersection(const Ray<3,double>& r,double s,unsigned IDfe=0) const;
    StepResult getIntersection(__m128,__m128,__m128 s) const;
} __attribute__ ((aligned(64)));

class TetraMesh {
    vector<unsigned>            T_m;        // tetra -> material mapping
    vector<unsigned>            T_r;        // tetra -> region mapping
	vector<Point<3,double> >    P;          // point vector
	vector<TetraByFaceID>	    T_f;        // tetra -> 4 face IDs
	vector<TetraByPointID>      T_p;        // tetra -> 4 point IDs
    vector<FaceByPointID>       F_p;        // face ID -> 3 point IDs
	vector<Face>			    F;          // faces (with normals and constants)
	vector<pair<int,int> >      vecFaceID_Tetra;
    vector<Tetra>               tetras;     // new SSE-friendly data structure

    // boundary data structures; map with key=volume-set ID, value=surface-set ID
    map<unsigned,unsigned> P_boundary_ID;
	map<unsigned,unsigned> F_boundary_ID;

    map<TetraByPointID,unsigned> tetraMap;

	// file input/output
	bool readFileMatlabTP(string);

	int  tetrasToFaces(vector<Face>&,vector<TetraByPointID>&,const vector<Point<3,double> >&,vector<TetraByFaceID>&);

    map<FaceByPointID,unsigned> faceMap;

	public:
    bool writeFileMatlabTP(string) const;

    // iterators
    typedef vector<Point<3,double> >::const_iterator point_const_iterator;
    point_const_iterator pointBegin() const { return P.begin()+1; }
    point_const_iterator pointEnd()   const { return P.end(); }

    typedef vector<Face>::const_iterator face_const_iterator;
    face_const_iterator faceBegin() const { return F.begin()+1; }
    face_const_iterator faceEnd()   const { return F.end(); }

    typedef vector<FaceByPointID>::const_iterator face_id_const_iterator;
    face_id_const_iterator faceIDBegin() const { return F_p.begin()+1; }
    face_id_const_iterator faceIDEnd()   const { return F_p.end(); }

    typedef vector<TetraByPointID>::const_iterator tetra_const_iterator;
    tetra_const_iterator tetraIDBegin() const { return T_p.begin()+1; }
    tetra_const_iterator tetraIDEnd()   const { return T_p.end();     }

    class boundary_f_const_iterator : public map<unsigned,unsigned>::const_iterator {
        const TetraMesh& __m;

        public:
        using map<unsigned,unsigned>::const_iterator::operator++;
        using map<unsigned,unsigned>::const_iterator::operator==;
        using map<unsigned,unsigned>::const_iterator::operator!=;
        boundary_f_const_iterator(const TetraMesh& m_,map<unsigned,unsigned>::const_iterator it_) :
            map<unsigned,unsigned>::const_iterator(it_),__m(m_){ };

        // returns the current face in terms of its surface point IDs (indexes into the boundary point set)
        FaceByPointID operator*() const {
            assert((*this)->first != 0);
            FaceByPointID f_original(__m.F_p[(*this)->first]);

            map<unsigned,unsigned>::const_iterator it=__m.P_boundary_ID.find((unsigned)abs((int)f_original[0]));
            assert(it != __m.P_boundary_ID.end());
            f_original[0] = it->second;

            it = __m.P_boundary_ID.find((unsigned)abs((int)f_original[1]));
            assert(it != __m.P_boundary_ID.end());
            f_original[1] = it->second;

            it = __m.P_boundary_ID.find((unsigned)abs((int)f_original[2]));
            assert(it != __m.P_boundary_ID.end());
            f_original[2] = it->second;

            return f_original;
            };
    };

    boundary_f_const_iterator boundaryFaceBegin() const { return boundary_f_const_iterator(*this,F_boundary_ID.begin()); }
    boundary_f_const_iterator boundaryFaceEnd()   const { return boundary_f_const_iterator(*this,F_boundary_ID.end());   }

	bool checkValid() const;
	bool writeFileMatlabF(string) const; 

	enum TetraFileType { MatlabTP };

	TetraMesh(){};
    TetraMesh(unsigned Np_,unsigned Nt_,unsigned Nf_) : T_m(Nt_+1),P(Np_+1),T_f(Nf_+1),T_p(Nt_+1),F_p(Nf_+1),F(Nf_+1),
        vecFaceID_Tetra(Nf_+1),tetras(Nt_+1){}
	TetraMesh(string,TetraFileType);
	TetraMesh(const vector<Point<3,double> >& P_,const vector<TetraByPointID>& T_p_,const vector<unsigned>& T_m_)
		: T_m(T_m_),T_r(T_m),P(P_),T_p(T_p_) { tetrasToFaces(F,T_p,P,T_f); }
    TetraMesh(const double*,unsigned Np,const unsigned*,unsigned Nt);
    ~TetraMesh();

    void fromBinary(const Blob& pts,const Blob& tetras,const Blob& faces=Blob());

	// query size of mesh
	unsigned getNp() const { return P.size()-1; };
	unsigned getNt() const { return T_f.size()-1; };
	unsigned getNf() const { return F.size()-1; }

    unsigned getNf_boundary() const { return F_boundary_ID.size(); }
    unsigned getNp_boundary() const { return P_boundary_ID.size(); }

	// Accessors for various point/face constructs
	const Point<3,double>&  getPoint(unsigned id)           const { return P[id]; }
    Face                    getFace(int id)                 const { Face f = F[abs(id)]; if(id<0){ f.flip(); } return f; }
    int                     getFaceID(FaceByPointID)        const;
    const FaceByPointID&    getFacePointIDs(unsigned id)    const { assert(id < F_p.size()); return F_p[id]; }
	const TetraByFaceID&    getTetraByFaceID(unsigned id)   const { return T_f[id]; }
	const TetraByPointID&   getTetraPointIDs(unsigned id)  const { return T_p[id]; }
    unsigned                getTetraFromFace(int IDf)       const;
    Point<3,double>         getTetraPoint(unsigned IDt,unsigned i) const { return P[T_p[IDt][i]]; }

    unsigned                getTetraID(TetraByPointID IDt) const {
        map<TetraByPointID,unsigned>::const_iterator it = tetraMap.find(IDt);
        return (it == tetraMap.end() ? 0 : it->second);
    }

    Tetra                   getTetra(unsigned IDt) const { return tetras[IDt]; }
    unsigned                getMaterial(unsigned IDt) const { return T_m[IDt]; }
    unsigned                getRegion(unsigned IDt) const { return T_r[IDt]; }

    // returns the tetra that the given face points into
    unsigned                getTetraIDFromFaceID(int IDf) const
        { return IDf > 0 ? vecFaceID_Tetra[IDf].first : vecFaceID_Tetra[-IDf].second; }

    double                  getFaceArea(const FaceByPointID&)   const;
    double                  getFaceArea(int IDf)                const { return getFaceArea(F_p[abs(IDf)]); }
    double                  getFaceArea(unsigned IDf)           const { return getFaceArea(F_p[IDf]); }

    double                  getTetraVolume(TetraByPointID IDps) const
        { return abs(scalartriple(P[IDps[0]],P[IDps[1]],P[IDps[2]],P[IDps[3]])/6); }
    double                  getTetraVolume(unsigned IDt) const { return getTetraVolume(T_p[IDt]); }

    const vector<Point<3,double> >& getPoints() const { return P; }

	// walk the graph (test code)
	void walk(Point<3,double>,const UnitVector<3,double>&,double) const;

	// find nearest point or enclosing tetra
	unsigned findEnclosingTetra(const Point<3,double>&) const;
	unsigned findNearestPoint  (const Point<3,double>&) const;
	pair<Point<3,double>,double> findNearestDelaunay(const Point<3,double>&) const;

    // checks if a point is within a given tetra by expressing as a linear combination of the corner points
    bool isWithinByPoints(int,const Point<3,double>&) const;

    // check if a ray intersects a face of a tetra
    StepResult getIntersection(unsigned,const Ray<3,double>&,double=std::numeric_limits<double>::infinity(),int=0) const;

    // find the surface element hit by an incoming ray
    pair<pair<unsigned,int>,Point<3,double> > getSurfaceElement(const Ray<3,double>&) const;

	// checks if faces are oriented correctly
	bool checkFaces() const;

    // does a number of data-structure integrity checks
    bool checkIntegrity(bool printResults=true) const;

    // functions for saving tetramesh representations
    pair<unsigned,boost::shared_array<const uint8_t> > tetrasAsBinary() const;
    pair<unsigned,boost::shared_array<const uint8_t> > pointsAsBinary() const;
};


template<class T>int signum(T a)
{
    return (a>T(0)) - (a<T(0));
}

#endif
