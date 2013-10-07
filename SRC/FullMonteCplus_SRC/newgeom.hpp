/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef INCLUDE_NEWGEOM
#define INCLUDE_NEWGEOM

#ifdef SSE
#include "sse.hpp"
#endif 

#include <math.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <utility>
#include <cassert>

#include "helpers.hpp"

// newgeom.hpp      Provides all the geometry-related classes (Point, arrays of IDs, etc)

using namespace std;

class Packet;

template<int D,class T>class Vector;

template<int D,class T>class FixedArray
{
	protected:
	T p[D];

	public:
	FixedArray()                    {};
	FixedArray(T t)                 { fill(p,p+D,t);     };
	FixedArray(const FixedArray& P) { copy(P.p,P.p+D,p); };
	FixedArray(const T* x_)         { if (x_!=NULL){ copy(x_,x_+D,p);   } };

    static bool lexCompare(FixedArray& a,FixedArray& b) { for(unsigned i=0;i<D;++i){ if(a[i] < b[i]) return true; else if (a[i]>b[i]) return false; } }

	const FixedArray& operator= (const FixedArray& P)       { copy(P.p,P.p+D,p); return *this; };
	bool              operator==(const FixedArray& P) const { for(int i=0; i<D; ++i){ if(P[i]!=p[i]){ return(false); } } return true; };

	T&       operator[](int i)       { return p[i]; };
	const T& operator[](int i) const { return p[i]; };

	template<int D_,class T_>friend std::ostream& operator<<(std::ostream&,const FixedArray<D_,T_>&);
	template<int D_,class T_>friend std::istream& operator>>(std::istream&,FixedArray<D_,T_>&);
};

template<int D,class T>class FixedArrayID : public FixedArray<D,T> {
	protected:
	using FixedArray<D,T>::p;

	public:
	using FixedArray<D,T>::operator[];
	FixedArrayID()                       : FixedArray<D,T>()  {};
	FixedArrayID(const T* p_)            : FixedArray<D,T>(p_){};
	FixedArrayID(const FixedArrayID& P_) : FixedArray<D,T>(P_){};

	// find min/max elements
	unsigned findMin() const;
	unsigned findMax() const;

    unsigned getMin() const { return p[findMin()]; }
    unsigned getMax() const { return p[findMax()]; }

	// rotate so min/max comes first
	FixedArrayID getRotateMin() const { FixedArrayID<D,T> tmp(*this); tmp.rotateMin(); return tmp; }
	void rotateMin();
	FixedArrayID getRotateMax() const { FixedArrayID<D,T> tmp(*this); tmp.rotateMax(); return tmp; }
	void rotateMax();
	FixedArrayID getSort()      const { FixedArrayID<D,T> tmp(*this); sort(tmp.p,tmp.p+D); return tmp; }

    // orderCount counts the number of times that a[i+1] < a[i], wrapping back to include a[N-1] < a[0]
    //   for a set of 3, this indicates face ordering
    //
    // ABC / BCA / CAB all count 1
    // ACB / CBA / BAC all count 2 (opposite ordering to set above)

    unsigned orderCount() const { unsigned n=0; for(unsigned i=0;i<D;++i){ n += p[(i+1)%D]>p[i]; } return n; }

    bool isSorted() const { for(unsigned i=1;i<D;++i){ if(p[i-1] > p[i]) return false; } return true; }

	// sort ascending
	typedef T* iterator;
	typedef const T* const_iterator;

	iterator begin() { return p;   }
	iterator   end() { return p+D; }
	const_iterator begin() const { return p;   }
	const_iterator   end() const { return p+D; }


	// lexicographical compare
	bool operator<(const FixedArrayID<D,T>& P) const { return lexicographical_compare(p,   p+D,   P.p, P.p+D); }
	bool operator>(const FixedArrayID<D,T>& P) const { return lexicographical_compare(P.p, P.p+D, p,   p+D);   }
};

class FaceByPointID;

class TetraByPointID : public FixedArrayID<4,unsigned> {
	using FixedArrayID<4,unsigned>::p;
	public:
	using FixedArrayID<4,unsigned>::operator=;
	TetraByPointID()                   : FixedArrayID<4,unsigned>()  {};
	TetraByPointID(const unsigned* p_) : FixedArrayID<4,unsigned>(p_){};
    TetraByPointID(const FixedArrayID<4,unsigned>& p_) : FixedArrayID<4,unsigned>(p_){}

	FaceByPointID getFace         (unsigned faceNum);
	unsigned      getOppositePoint(unsigned faceNum) const;
};



class TetraByFaceID : public FixedArrayID<4,int> {
	using FixedArrayID<4,int>::p;
	public:
	using FixedArrayID<4,int>::operator=;
	TetraByFaceID()              : FixedArrayID<4,int>()  {};
	TetraByFaceID(const int* p_) : FixedArrayID<4,int>(p_){};

	void flipFace(unsigned faceNum){ p[faceNum]=-p[faceNum]; }

	pair<bool,unsigned> getFace(unsigned faceNum){ bool inv=p[faceNum]<0; unsigned id=inv?-p[faceNum] : p[faceNum]; return make_pair(inv,id); }
};

class FaceByPointID : public FixedArrayID<3,unsigned> {
	public:
	using FixedArrayID<3,unsigned>::operator[];
	using FixedArrayID<3,unsigned>::operator=;
	FaceByPointID()                   : FixedArrayID<3,unsigned>()  {};
	FaceByPointID(const unsigned* p_) : FixedArrayID<3,unsigned>(p_){};
	FaceByPointID(const FixedArrayID<3,unsigned>& t_) : FixedArrayID<3,unsigned>(t_){};
	FaceByPointID(unsigned a,unsigned b,unsigned c){ p[0]=a; p[1]=b; p[2]=c; }

    // flips the orientation of the face [A,B,C] -> [A,C,B]
    FaceByPointID flip() const { return FaceByPointID(p[0],p[2],p[1]); }
};

template<int D,class T>unsigned FixedArrayID<D,T>::findMin() const
{
	unsigned m=std::numeric_limits<T>::max(),j=0;
	for(unsigned i=0;i<D;++i)
	{
		if (p[i]<m){
			m=p[i];
			j=i;
		}
	}
	return j;
}

template<int D,class T>unsigned FixedArrayID<D,T>::findMax() const
{
	unsigned m=std::numeric_limits<T>::min(),j=0;
	for(unsigned i=0;i<D;++i)
	{
		if (p[i]>m){
			m=p[i];
			j=i;
		}
	}
	return j;
}

template<int D,class T>void FixedArrayID<D,T>::rotateMin()
{
	T tmp[D];
	unsigned j=findMin();
	copy(p,p+D,tmp);
	for(unsigned i=0;i<D;++i)
		p[i]=tmp[(i+j)%D];
}

template<int D,class T>void FixedArrayID<D,T>::rotateMax()
{
	unsigned j=findMax();
	T tmp[D];
	copy(p,p+D,tmp);
	for(unsigned i=0;i<D;++i)
		p[i]=tmp[(i+j)%D];
}


// A Point is a FixedArray which represents a Cartesian position that can be translated by a vector
template<int D,class T>class Point : public FixedArray<D,T>
{
	protected:
	using FixedArray<D,T>::p;

	public:
	using FixedArray<D,T>::operator==;
	using FixedArray<D,T>::operator=;

	Point()                : FixedArray<D,T>(){};
	Point(const Point& P_) : FixedArray<D,T>(P_){};
	Point(const T* p_)     : FixedArray<D,T>(p_){};

#ifdef SSE
    operator __m128() const { return _mm_set_ps(0.0,p[2],p[1],p[0]); }
    void set(__m128 r){ float f[4]; _mm_store_ps(f,r); for(unsigned i=0;i<3;++i){ p[i]=f[i]; } }
#endif

	Point operator+(const Vector<D,T>& v) const { Point t; for(unsigned i=0;i<D;++i){ t[i]=p[i]+v[i]; } return t; }
	Point operator-(const Vector<D,T>& v) const { Point t; for(unsigned i=0;i<D;++i){ t[i]=p[i]-v[i]; } return t; }
};

template<int D,class T>std::ostream& operator<<(std::ostream& os,const FixedArray<D,T>& P)
{
	os << '(';
	for(int i=0; i<D; ++i){ os << P[i] << (i<D-1? ',' : ')'); }
	return os;
}

template<int D,class T>std::istream& operator>>(std::istream& is,FixedArray<D,T>& P)
{
	is >> std::skipws;
	bool paren=false;

	if(is.peek()=='('){ paren=true; is.ignore(1); }
	for(int i=0; i<D; ++i){ is >> P.p[i]; if(i < D-1 && is.peek()==','){ is.ignore(1); } }
	if (paren){ is.ignore(1); }
	return is;
}

// A Vector extends the Point class with a norm, dot product, cross product, add/sub and scalar multiply/divide
//    vector can be defined as going between two points, or implicitly as the origin (0,0,0) to a point
template<int D,class T>class Vector : public Point<D,T>
{
	protected:
	using Point<D,T>::p;

	public:
	using Point<D,T>::operator=;
	using Point<D,T>::operator==;
	using Point<D,T>::operator[];
	using Point<D,T>::operator-;
	using Point<D,T>::operator+;

	Vector()                                        : Point<D,T>()  {};
	Vector(const T* x_)                             : Point<D,T>(x_){};
	Vector(const Point<D,T>& P)                     : Point<D,T>(P) {};
	Vector(const Point<D,T>& A,const Point<D,T>& B){ for(int i=0;i<D;++i){ p[i]=B[i]-A[i]; } }

	// norms and dots
	T norm_l2()                 const { T s=0; for(int i=0;i<D;++i){ s += p[i]*p[i]; } return(sqrt(s)); }
	T norm_l1()                 const { T s=0; for(int i=0;i<D;++i){ s += abs(p[i]); } return abs(s);   }
	T norm2_l2()                const { T s=0; for(int i=0;i<D;++i){ s += p[i]*p[i]; } return s;        }
	T dot(const Vector<D,T>& a) const { T s=0; for(int i=0;i<D;++i){ s += a[i]*p[i]; } return s;        }

	// unary negate
	Vector operator-()           { Vector v; for(int i=0; i<D; ++i){ v[i]=-p[i]; } return v; }

	// vector += / -= operations
	const Vector& operator+=(const Vector& k) { for(int i=0; i<D; ++i){ p[i]+=k[i]; } return *this; }
	const Vector& operator-=(const Vector& k) { for(int i=0; i<D; ++i){ p[i]-=k[i]; } return *this; }

	// scalar operations
	const Vector& operator*=(const T& k) { for(int i=0; i<D; ++i){ p[i]*=k; } return *this; }
	const Vector& operator/=(const T& k) { for(int i=0; i<D; ++i){ p[i]/=k; } return *this; }
	Vector operator* (T k)        const { Vector<D,T> t(*this); return t *= k; }
	Vector operator/ (const T& k) const { Vector<D,T> t(*this); return t /= k; }

	// vector cross product
	Vector cross(const Vector<D,T>&) const;
	template<int D_,class U>friend Vector<D_,U> cross(const Vector<D_,U>&,const Vector<D_,U>&);
};

// helpers
template<int D,class T>T dot(const Vector<D,T>& a,const Vector<D,T>& b) { return a.dot(b); }
template<int D,class T>T norm_l2(const Vector<D,T>& a) { return a.norm_l2(); }
template<int D,class T>T norm_l1(const Vector<D,T>& a) { return a.norm_l1(); }
template<int D,class T>T norm2_l2(const Vector<D,T>& a){ return a.norm2_l2(); }

template<int D,class T>T norm2_l2(const Point<D,T>& a,const Point<D,T>& b){ T s; for(unsigned i=0;i<D;++i){ s += (a[i]-b[i])*(a[i]-b[i]); } return s; }


// A UnitVector is a Vector that is guaranteed to always have L2 norm 1
template<unsigned D,class T>class UnitVector : public Vector<D,T>
{
	protected:
	using Vector<D,T>::p;

	public:
	using Vector<D,T>::dot;
	using Vector<D,T>::norm_l2;
	using Vector<D,T>::norm2_l2;
    using Vector<D,T>::operator-;
	using Vector<D,T>::operator*;
	using Vector<D,T>::operator/;
	using Vector<D,T>::operator==;

    const static Tolerance<T> eps;

	UnitVector()                    { p[0]=1;          for(unsigned i=1;i<D;++i){ p[i]=0;      } };
    UnitVector(const T* v_,bool alreadyUnit=false){
        for (unsigned i=0;i<D;++i){ p[i]=v_[i]; }
        if (!alreadyUnit)
        {
            T L = norm_l2();
            for (unsigned i=0;i<D;++i){ p[i] /= L; }
        }
//        assert(norm_l2()==eps);
    }
	UnitVector(const Vector<D,T>& v,bool alreadyUnit=false){
        if (alreadyUnit)
            for(unsigned i=0;i<D;++i){ p[i]=v[i]; }
        else
        {
            T L=v.norm_l2(); for(unsigned i=0;i<D;++i){ p[i]=v[i]/L; }
        }
        assert(norm_l2()==eps);
    };
    UnitVector(const UnitVector<D,T>& v){
        for(unsigned i=0;i<D;++i){ p[i]=v.p[i]; }
    };

    UnitVector operator-() const { UnitVector t(*this); for(unsigned i=0;i<D;++i) t.p[i] = -t.p[i]; return t; }
};

UnitVector<3,double> uvect3FromPolar(double phi,double lambda);

// set tolerance for checking unit vectors
//template<>
//Tolerance<double> UnitVector<3,double>::eps<>(1.0,1e-9);

template<int D,class T> Vector<D,T> cross(const Vector<D,T>& a,const Vector<D,T>& b)
{
	return a.cross(b);
}

template<int D,class T>Vector<D,T> Vector<D,T>::cross(const Vector<D,T>& x) const
{
	T cp[3]= { p[1]*x[2]-p[2]*x[1], p[2]*x[0]-x[2]*p[0], p[0]*x[1]-p[1]*x[0] };
	return Vector<D,T>(cp);
}

// Computes the scalar triple product for points A,B,C,D = DA dot (DB cross DC)
template<class T>T scalartriple(const Point<3,T>& A,const Point<3,T>& B,const Point<3,T>& C,const Point<3,T>& D)
{
    return scalartriple(Vector<3,T>(D,A),Vector<3,T>(D,B),Vector<3,T>(D,C));
}

// Computes the scalar triple product for vectors a,b,c = a dot (b cross c)
template<class T>T scalartriple(const Vector<3,T>& a,const Vector<3,T>& b,const Vector<3,T>& c)
{
    return dot(a,cross(b,c));
}


// A Ray represents the geometric idea of a ray, a semi-infinite line starting from a point and extending along a unit vector
template<int D,class T>class Ray
{
	Point<D,T>      P;
	UnitVector<D,T> d;

	public:
	// construct from an origin point and a direction unit vector
    Ray(){};
	Ray(const Point<D,T>& P_,const UnitVector<D,T>& d_) : P(P_),d(d_){};

	// returns the point that is T units along the ray
	Point<D,T> operator()(T t) const { return Point<D,T>(P+d*t); };
//	const Ray& operator=()(const Ray& r){ P=r.P; d=r.d; return *this; }

	const Point<D,T>&      getOrigin()    const { return P; }
	const UnitVector<D,T>& getDirection() const { return d; }

    operator Packet() const;

    void setOrigin(Point<3,double> p_){ P=p_; }

	void print(ostream& os) const { os << "Ray: " << P << " " << d << std::endl; }
};
Point<3,double> pointFrom(__m128 p);
Ray<3,double> rayFrom(__m128 p,__m128 d);


UnitVector<3,double> uvectFrom(__m128 v);

// returns true if P+td falls within triangle defined by points T for some t
bool PointInTriangle(Point<3,double> p,UnitVector<3,double> d,Point<3,double> T[3],Point<3,double>& Q,double& t);


// IO manipulator for printing points (parentheses or not, commas or not, etc)

class GeomManip {
    char delimchar;
    bool parens;
    char parenchar[2];
    char uvparenchar[2];
    char idparenchar[2];
    ostream& os;
    public:

    explicit GeomManip(ostream& os_=cout) : delimchar(','),parens(false),os(os_){
        parenchar[0]='(';     parenchar[1]=')';
        uvparenchar[0] = '<'; uvparenchar[1]='>';
        idparenchar[0] = '['; idparenchar[1]=']';}

    static GeomManip plainwhite(){
        GeomManip gm;
        gm.parens=false;
        gm.delimchar=' ';
        return gm;
    }

    friend GeomManip operator<<(ostream&,const GeomManip&);

    template<int D,class T>friend ostream& operator<<(const GeomManip& gm,const Point<D,T>& p);
    template<int D,class T>friend ostream& operator<<(const GeomManip& gm,const FixedArrayID<D,T>& p);
};

template<int D,class T>ostream& operator<<(const GeomManip& gm,const UnitVector<D,T>& u)
{
    if(gm.parens)
        gm.os << gm.uvparenchar[0];
    for(unsigned i=0;i<(unsigned)D; ++i)
    {
        gm.os << u[i];
        if (i<D-1)
            gm.os << gm.delimchar;
    }
    if(gm.parens)
        gm.os << gm.uvparenchar[1];
    return gm.os;
}

template<int D,class T>ostream& operator<<(const GeomManip& gm,const Point<D,T>& p)
{
    if(gm.parens)
        gm.os << gm.parenchar[0];
    for(unsigned i=0;i<(unsigned)D; ++i)
    {
        gm.os << p[i];
        if (i<D-1)
            gm.os << gm.delimchar;
    }
    if(gm.parens)
        gm.os << gm.parenchar[1];
    return gm.os;
}

template<int D,class T>ostream& operator<<(const GeomManip& gm,const FixedArrayID<D,T>& p)
{
    if(gm.parens)
        gm.os << gm.parenchar[0];
    for(unsigned i=0;i<(unsigned)D; ++i)
    {
        gm.os << p[i];
        if (i<D-1)
            gm.os << gm.delimchar;
    }
    if(gm.parens)
        gm.os << gm.parenchar[1];
    return gm.os;
}

extern GeomManip plainwhite;

class Packet {
    public:
    __m128 d,a,b,p;     // 4x16B = 64B
    __m128 s;           // 16B
    double w;           // 8B

    Packet() : s(_mm_setzero_ps()),w(1.0){}
    Packet(const Ray<3,double>& r) : s(_mm_setzero_ps()),w(1.0){ setRay(r); }

    void setRay(const Ray<3,double>& r)
        { p = r.getOrigin(); setDirection(r.getDirection()); }

    void setDirection(__m128 d_)
        { d=d_; a=getNormalTo(d_); b=cross(d,a); }
};

// r0,r1 are the random numbers
inline Packet matspin(Packet pkt,float costheta,__m128 cosphi_sinphi)
{
    Packet res=pkt;
    // colums of matrix M (appearance in cout output below is transposed)
    __m128 M0,M1,M2;
    const __m128 d0=pkt.d, a0=pkt.a, b0=pkt.b;

    // rows of matrix M
//    M0 = _mm_setr_ps(costheta,sintheta,0,0);
//    M1 = _mm_setr_ps(-sintheta*cosphi,costheta*cosphi,sinphi,0);
//    M2 = _mm_setr_ps(sinphi*sintheta,-sinphi*costheta,cosphi,0);

	__m128 costheta_000 = _mm_set_ss(costheta);

	// calculation from inputs
	__m128 cost_sint = _mm_sqrt_ss(
		_mm_sub_ps(
			_mm_unpacklo_ps(_mm_set_ss(1.0),costheta_000),
			_mm_mul_ss(costheta_000,costheta_000)));

	__m128 trig = _mm_shuffle_ps(cost_sint,cosphi_sinphi,_MM_SHUFFLE(1,0,0,1));
	__m128 zero = _mm_setzero_ps();
	__m128 strig = _mm_addsub_ps(zero,trig);	// (-sin phi) (cos phi) (-sin theta) (cos theta)

	__m128 prods = _mm_mul_ps(strig,_mm_shuffle_ps(strig,strig,_MM_SHUFFLE(1,0,2,3)));

	__m128 cp_0_sp_0 = _mm_unpackhi_ps(trig,zero);	// (cos phi) 0 (sin phi) 0

	M0 = _mm_movelh_ps(trig,zero);
	M1 = _mm_shuffle_ps(prods,cp_0_sp_0,_MM_SHUFFLE(3,2,2,1));
	M2 = _mm_shuffle_ps(prods,cp_0_sp_0,_MM_SHUFFLE(3,0,0,3));

    res.d = _mm_mul_ps(d0,_mm_shuffle_ps(M0,M0,_MM_SHUFFLE(0,0,0,0)));
    res.d = _mm_add_ps(res.d,_mm_mul_ps(a0,_mm_shuffle_ps(M1,M1,_MM_SHUFFLE(0,0,0,0))));
    res.d = _mm_add_ps(res.d,_mm_mul_ps(b0,_mm_shuffle_ps(M2,M2,_MM_SHUFFLE(0,0,0,0))));

    res.a = _mm_mul_ps(d0,_mm_shuffle_ps(M0,M0,_MM_SHUFFLE(1,1,1,1)));
    res.a = _mm_add_ps(res.a,_mm_mul_ps(a0,_mm_shuffle_ps(M1,M1,_MM_SHUFFLE(1,1,1,1))));
    res.a = _mm_add_ps(res.a,_mm_mul_ps(b0,_mm_shuffle_ps(M2,M2,_MM_SHUFFLE(1,1,1,1))));

    res.b = _mm_mul_ps(d0,_mm_shuffle_ps(M0,M0,_MM_SHUFFLE(2,2,2,2)));
    res.b = _mm_add_ps(res.b,_mm_mul_ps(a0,_mm_shuffle_ps(M1,M1,_MM_SHUFFLE(2,2,2,2))));
    res.b = _mm_add_ps(res.b,_mm_mul_ps(b0,_mm_shuffle_ps(M2,M2,_MM_SHUFFLE(2,2,2,2))));

    return res;
}


#endif
