/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "graph.hpp"
#include <limits>
#include <map>
#include <set>
#include <cassert>
#include <signal.h>

#include "blob.hpp"
#include "sse.hpp"

using namespace std;

const float Material::const_c0=299.792458;

// compares two FixedArrays lexicographically
    template<unsigned D>class FACompare {
        public:
        bool operator()(const FixedArray<D,double>&a,const FixedArray<D,double>&b){
            for(unsigned i=0;i<D;++i){
                if(a[i] < b[i])
                    return true;
                else if (a[i] > b[i])
                    return false;
            }
            return false;
        }
    };

// constructor: loads a file of specified type & converts to internal representation
TetraMesh::TetraMesh(string fn,TetraFileType type)
{
	switch(type){
		case MatlabTP:
		readFileMatlabTP(fn);
		
		default:
		break;
	}
	tetrasToFaces(F,T_p,P,T_f);
}

TetraMesh::TetraMesh(const double* p,unsigned Np_,const unsigned* t,unsigned Nt_)
{
    P.clear();
    P.resize(Np_+1);
    T_m.clear();
    T_m.resize(Nt_+1);
    T_r.clear();
    T_r.resize(Nt_+1);
    T_p.clear();
    T_p.resize(Nt_+1);

    unsigned zeros[4]={0,0,0,0};

    P[0] = Point<3,double>();
    T_m[0] = 0;
    T_p[0] = TetraByPointID(zeros);
    T_r[0] = 0;

    for(unsigned i=1;i<=Np_;++i,p+=3)
        P[i] = Point<3,double>(p);

    for(unsigned i=1;i<=Nt_;++i,t+=5)
    {
        T_p[i]=TetraByPointID(t);        // direct copy from unsigned*
        T_m[i]=t[4];
        T_r[i]=T_m[i];
    }

	tetrasToFaces(F,T_p,P,T_f);
}

TetraMesh::~TetraMesh()
{
}

void TetraMesh::fromBinary(const Blob& pts,const Blob& tetras,const Blob& faces)
{
    P.clear();
    P.resize(pts.getSize()/24+1);
    T_m.clear();
    T_m.resize(tetras.getSize()/20+1);
    T_p.clear();
    T_p.resize(tetras.getSize()/20+1);

    unsigned zeros[4]={0,0,0,0};
    P[0] = Point<3,double>();
    T_m[0] = 0;
    T_p[0] = TetraByPointID(zeros);

    unsigned i=1;

    // copy points over
    for(const double* p=(const double*)pts.getPtr(); p < (const double*)pts.getEndPtr(); p += 3,++i)
    {
        P[i][0] = p[0];
        P[i][1] = p[1];
        P[i][2] = p[2];
    }

    i=1;
    for(const uint32_t* p=(const uint32_t*)tetras.getPtr(); p < (const uint32_t*)tetras.getEndPtr(); p += 5,++i)
    {
        T_p[i][0] = p[0];
        T_p[i][1] = p[1];
        T_p[i][2] = p[2];
        T_p[i][3] = p[3];
        T_m[i] = p[4];
    }
	tetrasToFaces(F,T_p,P,T_f);
}

vector<Material> loadMatFile(string fn)
{
    vector<Material> mat;
    float n,g,mu_s,mu_a;
    unsigned IDm;

    ifstream is(fn.c_str());

	while(!is.eof()){
		// ignore comments
		while(is.good() && is.peek()=='#')
		{
			char buf[100];
			is.get(buf,100,'\n');
			cout << "ignoring line || " << buf << endl;
			is.ignore(100,'\n');
		}

		// check for default statement
		if(is.good())
		{
			if(is.peek()=='?')
			{
				is.ignore(100,'\t');
				if(is.eof())
					break;
				Material m_default;
				is >> n >> g >> mu_s >> mu_a;
                m_default = Material(mu_a,mu_s,g,n,0,false);
				fill(mat.begin(),mat.end(),m_default);
			}
			else {
				is >> IDm;
				if (is.good())
                {
					is >> n >> g >> mu_s >> mu_a;
                    mat[IDm] = Material(mu_a,mu_s,g,n,0,false);
                }
				else
					break;
			}
		}
	}
    return mat;
}

/*
	readFileMatlabTP(fn,ids,coords)

Args
	fn			File name to open

Returns
	ids			Vector of TetraByPointID, specifies the tetrahedra in terms of their points
	coords		Vector of Point<3,double> giving the actual point coordinates

	pair<int,int> contains (number of faces, number of points)


	NOTE: Arrays are 1-based so they will agree with Matlab/Octave
			0 is used as a special value (for not-found or exited)

*/

bool TetraMesh::readFileMatlabTP(string fn)
{
	ifstream is(fn.c_str(),ios_base::in);

    if(!is.good()){
        cerr << "Failed to open " << fn << " for reading; abort" << endl;
        return false;
    }
	int Nt,Np;

    // read sizes
    is >> Np;
	is >> Nt;

	// read point coordinates -- uses 1-based addressing
	P.resize(Np+1);
	P[0]=Point<3,double>();
	for (vector<Point<3,double> >::iterator it = P.begin()+1; it != P.end(); ++it)
		is >> *it;

	T_p.resize(Nt+1);
    T_m.resize(Nt+1);
    tetraMap.clear();
	unsigned t[4]={0,0,0,0},i=1,max_m=0;
	T_p[0]=TetraByPointID(t);
    T_m[0]=0;
    TetraByPointID IDps;
	for (vector<TetraByPointID>::iterator it=T_p.begin()+1; it != T_p.end(); ++it,++i)
	{
		is >> IDps >> T_m[i];
        *it=IDps.getSort();
        tetraMap.insert(make_pair(*it,i));
		max_m = max(max_m,T_m[i]);
	}

	return true;
}

bool TetraMesh::writeFileMatlabTP(string fn) const
{
    ofstream os(fn.c_str());

    if(!os.good())
        return false;

    os << P.size()-1 << endl << T_p.size()-1 << endl;
    
    for(vector<Point<3,double> >::const_iterator it=P.begin()+1; it != P.end(); ++it)
        os << (*it)[0] << ' ' << (*it)[1] << ' ' << (*it)[2] << endl;
    
    vector<TetraByPointID>::const_iterator T_it;
    vector<unsigned>::const_iterator M_it;
    for(M_it=T_m.begin()+1, T_it=T_p.begin()+1; T_it != T_p.end(); ++M_it,++T_it)
        os << (*T_it)[0] << ' ' << (*T_it)[1] << ' ' << (*T_it)[2] << ' ' << (*T_it)[3] << ' ' << *M_it << endl;

    return !os.fail();
}

// writes the face mapping to a Matlab-compatible file
bool TetraMesh::writeFileMatlabF(string fn) const
{
	ofstream of(fn.c_str(),ios_base::out);

	if (!of.good())
		return false;

	assert(F.size()==vecFaceID_Tetra.size());

	of << vecFaceID_Tetra.size()-1 << endl;
	for(vector<pair<int,int> >::const_iterator it=vecFaceID_Tetra.begin()+1; it != vecFaceID_Tetra.end(); ++it)
		of << it->first << " " << it->second << endl;

	return true;
}



/* 	pair<unsigned id,bool flip> getFaceID(m,f)

Args
	m		Map of FaceByPointID to unsigned
	f		FaceByPointID to find

Returns
	id		Face ID
	flip	True if orientation is opposite
*/

int TetraMesh::getFaceID(FaceByPointID f) const
{
	bool neg;

	FaceByPointID tmpR = f.getRotateMin();

	map<FaceByPointID,unsigned>::const_iterator it;

	if ((neg = (tmpR[2] < tmpR[1])))
		swap(tmpR[2],tmpR[1]);

	it = faceMap.find(tmpR);
	assert(it != faceMap.end());

    return it == faceMap.end() ? 0 : (neg ? it->second : -it->second);
}


/* TetrasToFaces(v,vecTetra_PointID,P)

	v					Vector of Face objects (ID# is index)
	vecTetra_PointID	Vector of TetraByPointID (4-element lists of point IDs corresponding to v)
	P					Vector of Point<3,double> coordinates referenced by values above

*/

int TetraMesh::tetrasToFaces(vector<Face>& v,vector<TetraByPointID>& vecTetra_PointID,const vector<Point<3,double> >& P,vector<TetraByFaceID>& vecTetra_FaceID)
{
	typedef pair<FaceByPointID,int> value_type;
	typedef FaceByPointID mm_key_type;
	typedef pair<mm_key_type,unsigned> mm_pair_type;
	multimap <mm_key_type,unsigned> m;

	unsigned i=1;
	// loop through list of tetrahedra defined by point IDs, creating the tetraGraph
	for(vector<TetraByPointID>::const_iterator it=vecTetra_PointID.begin()+1; it != vecTetra_PointID.end(); ++it,++i)
	{

        assert((*it)[0] > 0 && (*it)[1] > 0 && (*it)[2] > 0 && (*it)[3] > 0);
        assert((*it)[0] <= P.size() && (*it)[1] <= P.size() && (*it)[2] <= P.size() && (*it)[3] <= P.size());

		// Create list of faces from tetrahedra (always sorting ascending)
		TetraByPointID tmp = it->getSort();

		// multimap m will map a sorted FaceByPointID to a Tetra vertex descriptor & the corresponding opposite point ID
		m.insert(mm_pair_type(FaceByPointID(tmp[0],tmp[1],tmp[2]),tmp[3]));
		m.insert(mm_pair_type(FaceByPointID(tmp[0],tmp[2],tmp[3]),tmp[1]));
		m.insert(mm_pair_type(FaceByPointID(tmp[0],tmp[1],tmp[3]),tmp[2]));
		m.insert(mm_pair_type(FaceByPointID(tmp[1],tmp[2],tmp[3]),tmp[0]));

        tetraMap.insert(make_pair(tmp,i));

		faceMap.insert(value_type(FaceByPointID(tmp[0],tmp[1],tmp[2]),0));
		faceMap.insert(value_type(FaceByPointID(tmp[0],tmp[2],tmp[3]),0));
		faceMap.insert(value_type(FaceByPointID(tmp[0],tmp[1],tmp[3]),0));
		faceMap.insert(value_type(FaceByPointID(tmp[1],tmp[2],tmp[3]),0));
	}

    // number faces from 1
    i=1;
    for(map<FaceByPointID,unsigned>::iterator it=faceMap.begin(); it != faceMap.end(); ++it,++i)
        it->second=i;

	// fill in edges of tetraGraph
	mm_pair_type a,b;
	unsigned Nf_regular=0;

    // P_boundary_ID maps a point ID to a new 1-based index of the boundary set
	for(multimap<FaceByPointID,unsigned>::const_iterator it=m.begin(); it != m.end();)
	{
		a=*(it++);
		if (it != m.end() && (*it).first==a.first)
		{
			++it;
			++Nf_regular;

			// check that there is not a face connecting more than two tetrahedra
			assert(!((*it).first == a.first));
		}
		else
        {
			// we do not have two tetrahedra adjoining this face: it is a boundary
            P_boundary_ID.insert(make_pair(a.first[0],0));
            P_boundary_ID.insert(make_pair(a.first[1],0));
            P_boundary_ID.insert(make_pair(a.first[2],0));

            map<FaceByPointID,unsigned>::const_iterator dIt = faceMap.find(a.first);
            assert(dIt != faceMap.end());

			F_boundary_ID.insert(make_pair(dIt->second,0));
        }
	}

    // renumber points from one
    i=1;
    for(map<unsigned,unsigned>::iterator it=P_boundary_ID.begin(); it != P_boundary_ID.end(); ++it,++i)
        it->second=i;

    // renumber faces from one
    i=1;
    for(map<unsigned,unsigned>::iterator it=F_boundary_ID.begin(); it != F_boundary_ID.end(); ++it,++i)
        it->second=i;

    // make sure point ID 0 is not included
    assert(P_boundary_ID.begin()->first != 0);

    // check that all faces are accounted for
	assert(2*Nf_regular+F_boundary_ID.size()==m.size());

	cout << "Of " << faceMap.size() << " faces, there are " << Nf_regular << " interior and " << F_boundary_ID.size() << " on the boundary " << endl;
    cout << "Of " << P.size() << " points, there are " << P_boundary_ID.size() << " on the boundary surface" << endl;

	//================================================================================
	// Create Face vector

	v.clear();
	v.resize(faceMap.size()+1); 
	vector<Face>::iterator vIt=v.begin()+1;

    F_p.clear();
	F_p.resize(faceMap.size()+1);
	vector<FaceByPointID>::iterator fIt=F_p.begin()+1;

	// Assign FaceByPointIDs to distinct faces, create vector of Face objects
	i=1;
	for(map<FaceByPointID,unsigned>::iterator it=faceMap.begin(); it != faceMap.end(); ++it)
	{
        assert(it->first[0] > 0 && it->first[0] <= P.size());
        assert(it->first[1] > 0 && it->first[1] <= P.size());
        assert(it->first[2] > 0 && it->first[2] <= P.size());
		it->second = i++;
		*(vIt++)=Face(P[it->first[0]],P[it->first[1]],P[it->first[2]]);
		*(fIt++)=it->first;
	}

	// Create vector of Tetra objects holding faceID tuples
	// and a vector mapping from faceID to Tetra
	// vecFaceID_Tetra[faceID].first is the up-face, and .second is the down-face
	vecTetra_FaceID.clear();
	vecTetra_FaceID.resize(vecTetra_PointID.size());
	vecFaceID_Tetra.clear();
	vecFaceID_Tetra.resize(v.size());

	i=1;
	for(vector<TetraByPointID>::const_iterator it=vecTetra_PointID.begin()+1; it != vecTetra_PointID.end(); ++i,++it){
		unsigned tmp[4];
		copy(&(*it)[0],&(*it)[4],tmp);

		// find the IDs of the 4 faces enclosing this tetra
		int id[4]= {
			getFaceID(FaceByPointID(tmp[0],tmp[1],tmp[2])),
			getFaceID(FaceByPointID(tmp[0],tmp[2],tmp[3])),
			getFaceID(FaceByPointID(tmp[0],tmp[3],tmp[1])),
			getFaceID(FaceByPointID(tmp[1],tmp[3],tmp[2])) };

		vecTetra_FaceID[i] = TetraByFaceID(id);

		// flip faces so that height of the point opposite each face is positive
		for(int j=0;j<4;++j){
			unsigned pID=it->getOppositePoint(j);
			assert(i>0);
			bool flip=vecTetra_FaceID[i][j]<0;
			if ((v[abs(vecTetra_FaceID[i][j])].pointHeight(P[pID])<0) ^ flip)
			{
				vecTetra_FaceID[i][j] = -vecTetra_FaceID[i][j];
				id[j]=-id[j];
			}
		}

		assert((unsigned)abs(id[0]) <= v.size() && (unsigned)abs(id[1]) <= v.size() &&
            (unsigned)abs(id[2]) <= v.size() && (unsigned)abs(id[3]) <= v.size()); 

		for(unsigned j=0;j<4;++j){
			if (id[j]>0)
				vecFaceID_Tetra[id[j]].first=i;
			else
				vecFaceID_Tetra[-id[j]].second=i;
		}

	}

	/* Validate by calculating height of point over opposite face */
	cout << setw(100) << setfill('=') << "" << setfill(' ') << endl;
	cout << "Post-check" << endl;
	checkFaces();

    cout << "Loading tetra structure" << endl;

    tetras.resize(getNt()+1);

    for(unsigned i=1;i<getNt()+1; ++i)
    {
        tetras[i].IDfs = getTetraByFaceID(i);
        tetras[i].matID = getMaterial(i);

        UnitVector<3,double> n[4];
        double C[4];
        Face F_tmp;
        for(unsigned j=0;j<4;++j)
        {
            F_tmp = getFace(tetras[i].IDfs[j]);
            tetras[i].adjTetras[j] = getTetraIDFromFaceID(-tetras[i].IDfs[j]);
            n[j] = F_tmp.getNormal();
            C[j] = F_tmp.getConstant();
        }

        tetras[i].nx=_mm_setr_ps(n[0][0],n[1][0],n[2][0],n[3][0]);
        tetras[i].ny=_mm_setr_ps(n[0][1],n[1][1],n[2][1],n[3][1]);
        tetras[i].nz=_mm_setr_ps(n[0][2],n[1][2],n[2][2],n[3][2]);
        tetras[i].C =_mm_setr_ps(C[0],C[1],C[2],C[3]);
    }


	return v.size();
}

// checks validity of TetraMesh construct
bool TetraMesh::checkValid() const
{
	bool valid=true;
	return valid;
}

bool TetraMesh::checkFaces() const
{
	bool status_ok=true;
	unsigned i=1;
	for(vector<TetraByPointID>::const_iterator it=T_p.begin()+1; it != T_p.end(); ++it,++i){
		for(int j=0;j<4;++j)
		{
			unsigned pID=it->getOppositePoint(j);
			Point<3,double> pt=P[pID];
			int f=T_f[i][j];
			if ((f<0?-1:1)*F[abs(f)].pointHeight(pt) < 0){
				status_ok=false;
				cout << "Error: height of opposite point to face " << f << " on tetrahedron " << i << " is negative" << endl;
			}
		}
	}
	return status_ok;
}


// Very naive search to find closest point
unsigned TetraMesh::findNearestPoint(const Point<3,double>& p) const
{
	double d2=numeric_limits<double>::infinity(),t;
	unsigned id=0,c=1;

	for(vector<Point<3,double> >::const_iterator it=P.begin()+1; it != P.end(); ++it,++c)
	{
		if ((t=norm2_l2(Vector<3,double>(*it,p))) < d2)
		{
			d2 = t;
			id = c;
		}
	}
	return id;
}

unsigned TetraMesh::findEnclosingTetra(const Point<3,double>& p) const
{
	unsigned id=0,c=1;
	for(vector<TetraByFaceID>::const_iterator it=T_f.begin()+1; it != T_f.end(); ++it,++c)
	{
		int i;
		for(i=0;i<4 && (F[abs((*it)[i])].pointAbove(p) ^ (((*it)[i])<0)); ++i){}
		if (i==4)
			id=c;
	}
	return id;
}

// verifies that a point is within the specified tetrahedron
bool Tetra::pointWithin(__m128 p)
{
    // compute p (dot) n_i minus C_i for i in [0,3]
    __m128 dot =         _mm_mul_ps(nx,_mm_shuffle_ps(p,p,_MM_SHUFFLE(0,0,0,0)));
    dot = _mm_add_ps(dot,_mm_mul_ps(ny,_mm_shuffle_ps(p,p,_MM_SHUFFLE(1,1,1,1))));
    dot = _mm_add_ps(dot,_mm_mul_ps(nz,_mm_shuffle_ps(p,p,_MM_SHUFFLE(2,2,2,2))));
    dot = _mm_sub_ps(dot,C);

    return _mm_movemask_ps(dot) == 0;
}

StepResult Tetra::getIntersection(const Ray<3,double>& r,double s,unsigned) const
{
    // convert to SSE vectors
    __m128 p = r.getOrigin();
    __m128 d = r.getDirection();

    return getIntersection(p,d,_mm_set_ss(s));
}

StepResult Tetra::getIntersection(__m128 p,__m128 d,__m128 s) const
{
    StepResult result;

    result.idx=-1;

    s = _mm_shuffle_ps(s,s,_MM_SHUFFLE(0,0,0,0));

    // calculate dot = n (dot) d, height = n (dot) p - C
    __m128 dot    =             _mm_mul_ps(nx,_mm_shuffle_ps(d,d,_MM_SHUFFLE(0,0,0,0)));
    __m128 h1 =             _mm_mul_ps(nx,_mm_shuffle_ps(p,p,_MM_SHUFFLE(0,0,0,0)));

    dot    = _mm_add_ps(dot,    _mm_mul_ps(ny,_mm_shuffle_ps(d,d,_MM_SHUFFLE(1,1,1,1))));
    h1 = _mm_add_ps(h1, _mm_mul_ps(ny,_mm_shuffle_ps(p,p,_MM_SHUFFLE(1,1,1,1))));

    dot    = _mm_add_ps(dot,    _mm_mul_ps(nz,_mm_shuffle_ps(d,d,_MM_SHUFFLE(2,2,2,2))));
    h1 = _mm_add_ps(h1, _mm_mul_ps(nz,_mm_shuffle_ps(p,p,_MM_SHUFFLE(2,2,2,2))));

    // height (=C - d dot n) should be negative if inside tetra, may occasionally be (small) positive due to numerical error
    // dot negative means facing outwards
    h1 = _mm_sub_ps(C,h1);

    // dist = height/dot
    __m128 dist = _mm_div_ps(h1,dot);

//  selects dist where dist>0 and dot<0 (facing outwards), s otherwise
    // very, very rarely ( < 1e-8? ) gives an error where no intersection is found
    dist = _mm_blendv_ps(s,dist,_mm_and_ps(_mm_cmpgt_ps(dist,_mm_setzero_ps()),dot));

    // at most three of the dot products should be negative
    // ideally none of the heights should be negative (assuming we're in the tetra)

    //      height  dot     h/dot   meaning
    //      +       +       +       OK: inside, facing away (no possible intersection)
    //      +       -       -       OK: inside, facing towards (intersection possible)
    //      -       +       -       OK: outside, facing in (this is the entry face with roundoff error, no possible intersection)
    //      -       -       +       ERROR: outside, facing out (problem!! this must be the entry face, but exiting!)

    // require p dot n - C > 0 (above face) and d dot n < 0

    pair<unsigned,__m128> min_idx_val = getMinIndex4p(dist);

    if (_mm_ucomilt_ss(min_idx_val.second,s))
    {
        result.IDfe = IDfs[min_idx_val.first&3];
        result.IDte = adjTetras[min_idx_val.first&3];
        result.idx = min_idx_val.first;					// will be 4 if no min found
    }
    result.distance=min_idx_val.second;
    result.Pe = _mm_add_ps(p,_mm_mul_ps(d,result.distance));

    return result;
}

// check linear combination of points to verify we're within tetra
//   used only for testing, very slow
bool TetraMesh::isWithinByPoints(int tID,const Point<3,double>& p) const
{
    float M[3][4];
    const Point<3,double> &A=P[T_p[tID][0]], &B=P[T_p[tID][1]], &C=P[T_p[tID][2]], &D=P[T_p[tID][3]];
    Vector<3,double> e[3];

    // calculate edge vectors
    e[0]=B-A;
    e[1]=C-A;
    e[2]=D-A;

    // build basis matrix for tetra
    M[0][0] = e[0][0];
    M[1][0] = e[0][1];
    M[2][0] = e[0][2];
    M[0][3] = p[0]-A[0];

    M[0][1] = e[1][0];
    M[1][1] = e[1][1];
    M[2][1] = e[1][2];
    M[1][3] = p[1]-A[1];

    M[0][2] = e[2][0];
    M[1][2] = e[2][1];
    M[2][2] = e[2][2];
    M[2][3] = p[2]-A[2];

    double c;

    // eliminate
    for(unsigned i=0;i<3;++i)
    {
        // subtract rows above
        for(unsigned j=0;j<i;++j)
        {
            c=M[i][j];
            for(unsigned k=0;k<4;++k)
                M[i][k] -= c*M[j][k];
        }

        // normalize the row
        c=M[i][i];
        for(unsigned j=i;j<4;++j)
            M[i][j] /= c;
    }

    // backsub
    for(int i=1;i>=0;--i)
        for(unsigned j=i+1;j<3;++j)
        {
            c=M[i][j];
            for(unsigned k=0;k<4;++k)
                M[i][k] -= c*M[j][k];
        }

    float coeff[4];
    coeff[0]=1.0;

    for(unsigned i=1;i<4;++i)
    {
        coeff[i]=M[i-1][3];
        coeff[0] -= coeff[i];
    }

    printf("Coeffs are: %9.5f %9.5f %9.5f %9.5f ",coeff[0],coeff[1],coeff[2],coeff[3]);

    bool within=true,onedge=false;

    for(unsigned i=0;i<4;++i)
    {
        within &= (coeff[i] > -1e-4) & (coeff[i] < 1+1e-4);
        onedge |= (abs(coeff[i]) < 1e-4);
    }

    cout << " - " << (within ? (onedge ? "on edge" : "within") : "OUTSIDE") << endl;

    return within;
}

// get the surface element hit by an incoming ray

// alg needs improvement; currently scans all boundary faces
// returns a pair of <tetraID,faceID>
pair<pair<unsigned,int>,Point<3,double> > TetraMesh::getSurfaceElement(const Ray<3,double>& r) const
{
    const UnitVector<3,double> &d=r.getDirection();
    const Point<3,double> &p=r.getOrigin();
    Point<3,double> Q;
    int IDf=0;
    unsigned IDt=0;
    double t,t_min = std::numeric_limits<double>::infinity();

    for(map<unsigned,unsigned>::const_iterator it=F_boundary_ID.begin(); it != F_boundary_ID.end(); ++it)
    {
        Point<3,double> A(P[F_p[it->first][0]]),B(P[F_p[it->first][1]]),C(P[F_p[it->first][2]]);
        Vector<3,double> OA(p,A),OB(p,B),OC(p,C);
        Vector<3,double> AB(A,B),AC(A,C);
        if (F[it->first].pointHeight(p) < 0 &&            // make sure we're on correct side
            dot(F[it->first].getNormal(),d) > 0        // and face is directed the correct way
            )
            {
                double c1,c2;
                pair<bool,Point<3,double> > tmp(F[it->first].rayIntersectPoint(r,t_min));
                Vector<3,double> PQ(p,tmp.second);
                if ((t=PQ.norm_l2()) < t_min){
                    Vector<3,double> AQ(A,tmp.second);
                    c1 = dot(AB,AQ)/AB.norm2_l2();
                    c2 = dot(Vector<3,double>(PQ-AB*c1),AC)/AC.norm2_l2();
                    // check that intersection point is within the face

                    if (0.0 <= c1 && c1 <= 1.0 && 0.0 <= c2 && c2 <= 1.0 && c1+c2 <= 1.0)
                    {
                        t_min = t;
                        IDf=it->first;
                        Q=tmp.second;
                    }
                }
            }
    }

    assert(IDf>0);
    IDt = vecFaceID_Tetra[IDf].first;

    return make_pair(make_pair(IDt,IDf),Q);
}

ostream& operator<<(ostream& os,const Material& mat)
{
    return os << "mu_a=" << mat.mu_a << " mu_s=" << mat.mu_s << " g=" << mat.g << " n=" << mat.n;
}

bool sameOrientation(FaceByPointID f0,FaceByPointID f1)
{
    return f0.orderCount() == f1.orderCount();
}

unsigned TetraMesh::getTetraFromFace(int IDf) const
{
    return IDf > 0 ? vecFaceID_Tetra[IDf].first : vecFaceID_Tetra[-IDf].second;
}

double TetraMesh::getFaceArea(const FaceByPointID& f) const
{
    Vector<3,double> AB(P[f[0]],P[f[1]]);
    Vector<3,double> AC(P[f[0]],P[f[2]]);

    return cross(AB,AC).norm_l2()/2;
}

pair<unsigned,boost::shared_array<const uint8_t> > TetraMesh::pointsAsBinary() const
{
    // create large object for points
    unsigned Nb = 3*sizeof(double)*getNp();
    double *p=new double[3*getNp()],*q=p;

    for(point_const_iterator it=pointBegin(); it != pointEnd(); ++it)
    {
        *(q++) = (*it)[0];
        *(q++) = (*it)[1];
        *(q++) = (*it)[2];
    }
    return make_pair(Nb,boost::shared_array<const uint8_t>((const uint8_t*)p));
}

pair<unsigned,boost::shared_array<const uint8_t> > TetraMesh::tetrasAsBinary() const
{
    unsigned Nb,Nt = getNt();
    unsigned *t = new unsigned[5*Nt],*u = t;
    Nb=5*sizeof(unsigned)*Nt;
    vector<unsigned>::const_iterator Mit=T_m.begin()+1;

    for(tetra_const_iterator it=tetraIDBegin(); it != tetraIDEnd(); ++it)
    {
        *(u++) = (*it)[0];
        *(u++) = (*it)[1];
        *(u++) = (*it)[2];
        *(u++) = (*it)[3];
        *(u++) = *(Mit++);
    }
    return make_pair(Nb,boost::shared_array<const uint8_t>((const uint8_t*)t));
}

// TODO: Add more integrity checks
bool TetraMesh::checkIntegrity(bool printResults) const
{
    bool status_ok=true,status_fp=true;

    if (printResults)
        cout << "Checking face ordering in F_p - " << flush;
    for(vector<FaceByPointID>::const_iterator it=F_p.begin()+1; it != F_p.end(); ++it)
        status_fp &= ((*it)[0] < (*it)[1]) & ((*it)[1] < (*it)[2]);

    if (printResults)
        cout << (status_fp ? "OK" : "Error") << endl;
    status_ok &= status_fp;


    if (printResults)
        cout << "Integrity check complete - status " << (status_ok ? "OK" : "Error") << endl;

    return status_ok;
}

