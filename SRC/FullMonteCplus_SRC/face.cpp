/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "graph.hpp"
#include <limits>
#include <map>
#include <cassert>

using namespace std;

// TODO: Add a routine to calculate face area using cross product
// TODO: Map from uniform square coords to face points (for face source type)

double Face::pointHeight(const Point<3,double>& P) const
{
	return dot((Vector<3,double>)P,normal)-C;
}

// queries whether a point is above the face or not
bool Face::pointAbove(const Point<3,double>& p) const
{
	return pointHeight(p)>=0;
}

// creates a Face with normal and constant by using cross-product
// face coordinates should be given in clockwise order, such that the remaining point is above AB x AC
Face::Face(const Point<3,double>& Pa,const Point<3,double>& Pb,const Point<3,double>& Pc) :
    normal(cross(Vector<3,double>(Pa,Pb),Vector<3,double>(Pa,Pc))),
    C (dot(normal,Vector<3,double>(Pa)))
    { }

// creates a Face from four points, ensuring that AD (dot) n > 0, ie A is above ABC
Face::Face(const Point<3,double>& Pa,const Point<3,double>& Pb,const Point<3,double>& Pc,const Point<3,double>& Pd) :
	normal(cross(Vector<3,double>(Pa,Pb),Vector<3,double>(Pa,Pc))),
	C (dot(normal,Vector<3,double>(Pa)))
{
	if (!pointAbove(Pd))
		flip();
}

// Intersects a ray with the face plane
// NOTE: does not guarantee that the point is actually inside the face points
pair<bool,Point<3,double> > Face::rayIntersectPoint(const Ray<3,double>& r,bool invert) const
{
	double costheta=-dot(r.getDirection(),normal);
	double h=pointHeight(r.getOrigin());
	if (signum(costheta)!=signum(h))
		return make_pair(false,r.getOrigin());
	return make_pair(true,r(h/costheta*(invert ? -1 : 1)));
}

// checks if a ray intersects the face plane within t units of travel
// NOTE: does not guarantee that the point is actually inside the face points
pair<bool,double> Face::rayIntersect(const Ray<3,double>& r,double t,bool invert) const
{
	double costheta=-dot(r.getDirection(),normal);
	double h=pointHeight(r.getOrigin());
	if(signum(costheta)!=signum(h) || (invert ^ (h >= t*costheta)))
		return make_pair(false,t);
	else
		return make_pair(true,h/costheta);
}

// pair<UnitVector<2,double>,UnitVector<3,double> > Face::reflectionBasis(v,flip)
//
// v        Incoming direction unit vector
// flip     Boolean flag, indicating if face must be flipped

// sign conventions:
//  normal is positive going from m1 to m2 if flip is false

// result.first     2D unit vector giving components normal (cos theta_i) and along plane (sin theta_i)
// result.second    3D unit vector indicating direction of along-plane component
//
// NOTES
//  in-plane component (sin theta) is always positive
//  normal component (cos theta) may be positive or negative depending on face orientation

/*pair<UnitVector<2,double>,UnitVector<3,double> > Face::reflectionBasis(const UnitVector<3,double>& v,bool flip) const
{
    double d=dot(v,normal);
    assert(abs(d)<=1.0);

    // calculate in-plane component
    Vector<3,double> v_p;
    if (abs(d) > 0.999999)
        v_p = UnitVector<3,double>();
    else
	    v_p=v - normal * d;
    double tmp[2] = { d,sqrt(1-d*d) };
    return make_pair (UnitVector<2,double>(tmp,true), UnitVector<3,double>(v_p,false));
}

Vector<3,double> Face::project(const Vector<3,double>& v) const
{
	return v-v*dot(v,normal);
}

Vector<3,double> Face::normalComponent(const Vector<3,double>& v) const
{
	return v*dot(v,normal);
}
*/
ostream& operator<<(ostream& os,const Face& f)
{
	return os << "n=" << f.normal << " C=" << f.C;
}

