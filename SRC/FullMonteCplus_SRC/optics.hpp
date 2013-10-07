/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef INCLUDE_OPTICS
#define INCLUDE_OPTICS
#include "newgeom.hpp"
#include <boost/math/constants/constants.hpp>

template<class T>UnitVector<2,T> Refract(double n1,double n2,const UnitVector<2,T>& v_i);
template<class T>UnitVector<2,T> Reflect(const UnitVector<2,T>& v_i);
template<class T>double FresnelReflect(double n1,double n2,const UnitVector<2,T>& v_i);

template<class T>UnitVector<3,T> Scatter(double g,const UnitVector<3,T>& d0,double r0,double r1);
template<class T>double FresnelReflect(double n1,double n2,const UnitVector<2,T>& v_i,const UnitVector<2,T>& v_t);
template<class T>double FresnelReflect(double n1,double n2,const UnitVector<2,T>& v_i);

/* THIS IS THE **OLD** scattering routine; new one does away with branches and sqrt (contained in UnitVector constructor) */

template<class T>UnitVector<3,T> Scatter(const UnitVector<3,T>& d0,double costheta,double sintheta,double cosphi,double sinphi)
{
    const double &dx=d0[0],&dy=d0[1],&dz=d0[2];
    double d[3] = { d0[0]*costheta,d0[1]*costheta,d0[2]*costheta };

    // find largest Cartesian component to use for cross product
    unsigned min_i=0;
    double min_u=abs(d0[0]),tmp;
    if ((tmp=abs(d0[1])) < min_u)
        min_i=1,min_u=tmp;
    if ((tmp=abs(d0[2])) < min_u)
        min_i=2,min_u=tmp;

    double k=1/sqrt(1-min_u*min_u);

    switch(min_i){
        case 0: // use X component, n1 = (d) cross (i); n2 = (d) cross (n1)
        d[0]+=k*sintheta*(0 - sinphi*(dy*dy+dz*dz));
        d[1]+=k*sintheta*(cosphi*dz + sinphi*dx*dy);
        d[2]+=k*sintheta*(-cosphi*dy + sinphi*dx*dz); 
        break;
        case 1: // use Y component
        d[0]+=k*sintheta*(-cosphi*dz + sinphi*dx*dy);
        d[1]+=k*sintheta*(0 - sinphi*(dx*dx+dz*dz));
        d[2]+=k*sintheta*(cosphi*dx + sinphi*dy*dz);
        break;
        case 2: // use Z component
        d[0]+=k*sintheta*(cosphi*dy + sinphi*dx*dz);
        d[1]+=k*sintheta*(-cosphi*dx + sinphi*dy*dz);
        d[2]+=k*sintheta*(0 - sinphi*(dx*dx+dy*dy)); 
        break;
        default:
        assert(0);
        break;
    }

    // includes assertion check that it actually is unit ('true' argument below does not scale by length)
    UnitVector<3,double> d_v(d,false);          // TODO: See if must be false

    return d_v;
}

// For the optics routines below, the following unit vector names are used:
//    v_i    Incidence
//    v_r    Reflection
//    v_t    Transmission
//
// Vector v has components (cos (theta), cos (phi))
//    v[0]   theta   Angle to normal / component along normal
//    v[1]   phi     Angle to surface / component perpendicular to normal (cos phi = sin theta)
//
// In the surface basis defined above, v_i[0] and v_i[1] are always positive facing into the interface

// Refract(n1,n2,d)
//    Calculates refraction at a boundary
template<class T>UnitVector<2,T> Refract(double n1,double n2,const UnitVector<2,T>& v_i)
{
	// Snell's law n_i * sin(theta_i) = n_t * sin(theta_t)
	T cosines[2];
	cosines[1] = n1/n2*v_i[1];

    assert(!isnan(v_i[0]) && !isnan(v_i[1]));
    assert(abs(cosines[1])<=1.0);

	cosines[0] = sqrt(1 - cosines[1]*cosines[1]);
	return UnitVector<2,T>(cosines,true);
}

// Reflection takes a unit vector in the direction of incidence v_i = (n,p)
//   where n is the component along the normal, and p is the component projected on the surface

template<class T>UnitVector<2,T> Reflect(const UnitVector<2,T>& v_i)
{
	UnitVector<2,T> v_r;
	v_r[0]=-v_i[0];
	v_r[1]=v_i[1];
	return v_r;
}

// double FresnelReflect(n1,n2,v_i,v_t)
//
//  n1      Refractive indices
//  n2      
//  v_i     Incident 2D unit vector (incidence-plane coords)
//  v_t     Transmitted 2D unit vector (incidence-plane coords)

template<class T>double FresnelReflect(double n1,double n2,const UnitVector<2,T>& v_i,const UnitVector<2,T>& v_t)
{
	// Rs = [(n1 cos(theta_i) - n2 cos(theta_t))/((n1 cos(theta_i) + n2 cos(theta_t)))] ^ 2
	// Rp = [(n1 cos(theta_t) - n2 cos(theta_i))/((n1 cos(theta_t) + n2 cos(theta_i)))] ^ 2

	double Rs= (n1*v_i[0] - n2*v_t[0]) / (n1*v_i[0] + n2*v_t[0]);
	Rs *= Rs;
	double Rp= (n1*v_t[0] - n2*v_i[0]) / (n1*v_t[0] + n2*v_i[0]);
	Rp *= Rp;

	return (Rs+Rp)/2;
}

// Version of FresnelReflect where the refraction vectors are not given; calculates them first

template<class T>double FresnelReflect(double n1,double n2,const UnitVector<2,T>& v_i)
{
	UnitVector<2,T> v_t=Refract(n1,n2,v_i);
	return FresnelReflect(n1,n2,v_i,v_t);
}

template<class T>UnitVector<3,T> Scatter(double g,const UnitVector<3,T>& d0,double rnd0,double rnd1)
{
    double costheta,sintheta,phi;

    double t=(1-g*g)/(1+g*(2*rnd0-1));

    // choose angles: HG function for component along d0, uniform circle for normal components
    costheta = (g==Tolerance<double>(0.0,1e-6)) ? 2*rnd0-1 : 1.0/2.0/g*(1+g*g-t*t);
    assert (costheta <= 1.0 && costheta >= -1.0);
    sintheta = sqrt(1-costheta*costheta);

    phi = 2*boost::math::constants::pi<double>()*rnd1;

    return Scatter(d0,costheta,sintheta,cos(phi),sin(phi));
}

#endif
