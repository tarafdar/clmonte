/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "newgeom.hpp"

// Points are ABC clockwise, with D above ABC
// Faces are 0=ABC, 1=ACD, 2=ADB, 3=BDC

template<> const Tolerance<double> UnitVector<2,double>::eps(1.0,1e-4);
template<> const Tolerance<double> UnitVector<3,double>::eps(1.0,1e-4);
template<> const Tolerance<float>  UnitVector<2,float>::eps(1.0,1e-3);

// for a triangle with points T[0..2], determines if P+td is in the triangle for some t

// tries to solve A + c0(B-A) + c1(C-A) = P - c2*d
// equivalently (1-c0-c1)A + c0(B) + c1(C) + c2(d) = P
bool PointInTriangle(Point<3,double> P,UnitVector<3,double> d,Point<3,double> T[3],Point<3,double>& Q,double& t)
{
    double a[3] = { T[1][0]-T[0][0], T[1][1]-T[0][1], T[1][2]-T[0][2] };    // vector from point T[0] to T[1]
    double b[3] = { T[2][0]-T[0][0], T[2][1]-T[0][1], T[2][2]-T[0][2] };    // vector from T[0] to T[2]
    double M[3][3] = { { a[0],b[0],-d[0] }, {a[1],b[1],-d[1]}, {a[2],b[2],-d[2]} };

    // make determinant
    double det = M[0][0]*(M[2][2]*M[1][1]-M[2][1]*M[1][2])
        -M[1][0]*(M[2][2]*M[0][1]-M[2][1]*M[0][2])
        +M[2][0]*(M[1][2]*M[0][1]-M[1][1]*M[0][2]);

    // make inverse (needs to multiplied by 1/det)
    double inv[3][3] = { { M[2][2]*M[1][1]-M[2][1]*M[1][2], -(M[2][2]*M[0][1]-M[2][1]*M[0][2]), M[1][2]*M[0][1]-M[1][1]*M[0][2] },
        { -(M[2][2]*M[1][0]-M[2][0]*M[1][2]), M[2][2]*M[0][0]-M[2][0]*M[0][2], -(M[1][2]*M[0][0]-M[1][0]*M[0][2]) },
        {   M[2][1]*M[1][0]-M[2][0]*M[1][1],-(M[2][1]*M[0][0]-M[2][0]*M[0][1]),  M[1][1]*M[0][0]-M[1][0]*M[0][1] } };

    double A[3] = { P[0]-T[0][0], P[1]-T[0][1],P[2]-T[0][2] };

    double c[3] = { (inv[0][0]*A[0]+inv[0][1]*A[1]+inv[0][2]*A[2])/det,
        (inv[1][0]*A[0]+inv[1][1]*A[1]+inv[1][2]*A[2])/det,
        (inv[2][0]*A[0]+inv[2][1]*A[1]+inv[2][2]*A[2])/det };

    // c[0], c[1] = coefficients for AB, AC respectively
    // c[2] = distance from P to intersection along d

    double tmp2[3] = {
        (1-c[0]-c[1])*T[0][0] + c[0]*T[1][0] + c[1]*T[2][0],
        (1-c[0]-c[1])*T[0][1] + c[0]*T[1][1] + c[1]*T[2][1],
        (1-c[0]-c[1])*T[0][2] + c[0]*T[1][2] + c[1]*T[2][2]};
    Point<3,double> q_test(tmp2);


    // test that P = (1-c0-c1)A + c0 B + c1 C - c2 (d)
    //   (ie. that travelling c2 units from d gets you to the face

    double tmp[3] = {
        (1-c[0]-c[1])*T[0][0] + c[0]*T[1][0] + c[1]*T[2][0] - c[2]*d[0],
        (1-c[0]-c[1])*T[0][1] + c[0]*T[1][1] + c[1]*T[2][1] - c[2]*d[1],
        (1-c[0]-c[1])*T[0][2] + c[0]*T[1][2] + c[1]*T[2][2] - c[2]*d[2]};

    Point<3,double> p_test(tmp);

    bool is_inside= (c[0] <= 1.0) & (c[0] >= 0.0) & (c[1] <= 1.0) & (c[1] >= 0.0) && (c[0]+c[1] <= 1.0);

    if(is_inside)
    {
        t=c[2];
        Q=q_test;
    }

/*    cout << "P=" << P << endl;
    cout << " Coeffs " << c[0] << ' ' << c[1] << ' ' << c[2] << (is_inside ? " IN " : " OUT")  << endl;
    cout << " Q_test= " << q_test << endl << " P_test=" << p_test << endl;*/

    return is_inside;
}

FaceByPointID TetraByPointID::getFace(unsigned faceNum)
{
	unsigned tmp[3] = { p[0],p[1],p[2] };
	switch(faceNum){
		case 0: break;
		case 1: tmp[1]=p[2]; tmp[2]=p[3]; break;
		case 2: tmp[1]=p[3]; tmp[2]=p[1]; break;
		case 3: tmp[0]=p[1]; tmp[1]=p[3]; break;
		default: assert(0);
	}
	return FaceByPointID(tmp);
}

unsigned TetraByPointID::getOppositePoint(unsigned faceNum) const
{
	switch(faceNum){
		case 0: return p[3];
		case 1: return p[1];
		case 2: return p[2];
		case 3: return p[0]; 
		default: assert(0); 
	}
    return -1;
}

UnitVector<3,double> uvect3FromPolar(double phi,double lambda)
{
    double p[3] = { cos(phi)*sin(lambda), cos(phi)*cos(lambda), sin(phi) };
    return UnitVector<3,double>(p);
}

GeomManip plainwhite=GeomManip::plainwhite();

GeomManip operator<<(ostream& os,const GeomManip& gm_)
{
    GeomManip gm(os);
    gm.delimchar = gm_.delimchar;
    gm.parens = gm_.parens;
    gm.parenchar[0] = gm_.parenchar[0];
    gm.parenchar[1] = gm_.parenchar[1];
    gm.uvparenchar[0] = gm_.uvparenchar[0];
    gm.uvparenchar[1] = gm_.uvparenchar[1];
    gm.idparenchar[0] = gm_.idparenchar[0];
    gm.idparenchar[1] = gm_.idparenchar[1];
    return gm;
}

Point<3,double> pointFrom(__m128 p)
{
    Point<3,double> P;
    float tmp_f[4];
    _mm_store_ps(tmp_f,p);

    for(unsigned i=0;i<3;++i)
        P[i]=tmp_f[i];
    return P;
}

Ray<3,double> rayFrom(__m128 p,__m128 d)
{
    Point<3,double> pv;
    float tmp_f[4];
    double tmp_d[4];
    _mm_store_ps(tmp_f,d);
    for(unsigned i=0;i<3;++i)
        tmp_d[i]=tmp_f[i];
    UnitVector<3,double> dv(tmp_d,true);

    _mm_store_ps(tmp_f,p);
    for(unsigned i=0;i<3;++i)
        pv[i] = tmp_f[i];
    return Ray<3,double>(pv,dv);
}

UnitVector<3,double> uvectFrom(__m128 v)
{
    float f[4];
    _mm_store_ps(f,v);
    double d[3];
    for(unsigned i=0;i<3;++i)
        d[i]=f[i];
    return UnitVector<3,double>(d);
}

// r0,r1 are the random numbers
/*Packet matspin(Packet pkt,double costheta,double sintheta,double cosphi,double sinphi)
{
    Packet res=pkt;
    // colums of matrix M (appearance in cout output below is transposed)
    __m128 M0,M1,M2;
    const __m128 d0=pkt.d, a0=pkt.a, b0=pkt.b;

    // rows of matrix M
    M0 = _mm_setr_ps(costheta,sintheta,0,0);
    M1 = _mm_setr_ps(-sintheta*cosphi,costheta*cosphi,sinphi,0);
    M2 = _mm_setr_ps(sinphi*sintheta,-sinphi*costheta,cosphi,0);

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
}*/

// r0,r1 are the random numbers
/*inline Packet matspin(Packet pkt,float costheta,__m128 cosphi_sinphi)
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

    // rows of matrix M
	// "golden" copy
//    M0 = _mm_setr_ps(costheta,sintheta,0,0);
//    M1 = _mm_setr_ps(-sintheta*cosphi,costheta*cosphi,sinphi,0);	// 0 sinphi (costheta * cosphi)  (-sintheta * cosphi)
//    M2 = _mm_setr_ps(sinphi*sintheta,-sinphi*costheta,cosphi,0);	// 0 cosphi (-sinphi * costheta) (sinphi * sintheta)

	__m128 zero = _mm_setzero_ps();

	__m128 strig = _mm_addsub_ps(zero,trig);	// (-sin phi) (cos phi) (-sin theta) (cos theta)

	__m128 prods = _mm_mul_ps(strig,_mm_shuffle_ps(strig,strig,_MM_SHUFFLE(1,0,2,3)));
		// prods = (sintheta*sinphi) (costheta*cosphi) (-sintheta*cosphi) (-costheta*sinphi)

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
}*/

template<>Ray<3,double>::operator Packet() const
{
    Packet p;

    float f[4];
    f[0] = P[0];
    f[1] = P[1];
    f[2] = P[2];
    f[3] = 0;

    p.p=_mm_load_ps(f);

    f[0] = d[0];
    f[1] = d[1];
    f[2] = d[2];
    p.setDirection(_mm_load_ps(f));
    return p;
}
