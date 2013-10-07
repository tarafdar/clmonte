/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <math.h>
#include "linefile.hpp"
#include "graph.hpp"
#include "source.hpp"
#include "newgeom.hpp"
#include "fluencemap.hpp"

#include "io_timos.hpp"

using namespace std;

// Reads a TIM-OS .source file
vector<Source*> readTIMOSSource(string fn)
{
    unsigned Ns,i,stype,Np,IDt;
    Point<3,double> pos;
    UnitVector<3,double> dir;
    FaceByPointID f;
    vector<Source*> sources;

    LineFile is(fn,'%',cerr);

    // First line: number of sources
    is >> Ns >> LineFile::LF_EOL;
    sources.reserve(Ns);

    // Loop over sources, one per line
    for(i=0;i<Ns;++i)
    {
        is >> stype;
        switch(stype){
            case 1:
            is >> pos >> Np;
            sources.push_back(new IsotropicPointSource(pos,Np));
			break;

            case 2:
            is >> IDt >> Np;
            sources.push_back(new VolumeSource(IDt,Np));
            break;

            case 11:
            is >> IDt >> pos >> dir >> Np;
            // pos is the position where it hits the geometry, not source location!
            sources.push_back(new PencilBeamSource(pos,dir,Np,IDt));
            break;

            case 12:
            is >> f >> Np;
            sources.push_back(new FaceSource(f,Np));
            break;

            default:
            is.error() << "Unrecognized source type " << stype << endl;
            break;
        }
        if (i < Ns-1)
            is >> LineFile::LF_EOL;         // Check for end of line
        else
            is >> LineFile::LF_EOLS;        // May be followed by zero or more newlines
    }
    is >> LineFile::LF_EOF;

    double sum_Np=0.0;
    for(vector<Source*>::const_iterator it=sources.begin(); it != sources.end(); ++it)
        sum_Np += (*it)->getPower();

    return sources;
}

bool writeTIMOSSource(string fn,const vector<Source*>& s,long long packetcount)
{
    double w_sum=0.0;
    unsigned long Np;
    long long Np_total=0;
    ofstream os(fn.c_str());

    if (!os.good())
        return false;

    // First line: number of sources
    os << s.size() << endl;

    if (packetcount != 0){
        for(vector<Source*>::const_iterator it=s.begin(); it != s.end(); ++it)
            w_sum += (*it)->getPower();
    
        for(vector<Source*>::const_iterator it=s.begin(); it != s.end()-1; ++it)
        {
            Np = rint(((*it)->getPower()/w_sum) * (double)packetcount);
            Np_total += Np;
            os << (*it)->timos_str(Np) << endl;
        }
        os << s.back()->timos_str(packetcount-Np_total);
    }
    else {
        for(vector<Source*>::const_iterator it=s.begin(); it != s.end(); ++it)
            os << (*it)->timos_str() << endl;
    }

    return os.good();
}

// Reads a TIM-OS Material (.opt) file
// Material 0 is surrounding environment
vector<Material> readTIMOSMaterials(string fn)
{
    unsigned Nm,i,ostyle,etype;
    double mu_a,mu_s,g,n,n_e;
    LineFile is(fn,'%',cerr);

    // First line: material parameter style (must be 1=per-region)
    is >> ostyle >> LineFile::LF_EOL;
    if (ostyle != 1)
    {
        is.error() << "Unsupported material parameter style " << ostyle << ": must be 1 (per-region)" << endl << "2 (per-tetra) is not supported" << endl;
        return vector<Material>();
    }

    // Second line: number of materials
    is >> Nm >> LineFile::LF_EOL;

    vector<Material> mat(Nm+1);

    for(i=1;i<=Nm;++i)
    {
        is >> mu_a >> mu_s >> g >> n >> LineFile::LF_EOL;
        mat[i] = Material(mu_a,mu_s,g,n);
    }

    // After materials: 1 line for environment type
    is >> etype >> LineFile::LF_EOL;

    if (etype == 1)
        is >> n_e >> LineFile::LF_EOL;
    else if (etype == 2)
    {
        is.error() << "Unvalidated environment type (matched boundary)" << endl;
        mat[0].setMatched(true);
    }
    else
    {
        is.error() << "Invalid environment type; must be 1 (uniform) or 2 (matched)" << endl;
        return vector<Material>();
    }

	mat[0] = Material(0,0,0,n_e);

    return mat;
}

// writes out TIM-OS compatible materials defs
bool writeTIMOSMaterials(string fn,const vector<Material>& mat)
{
    ofstream os(fn.c_str(),ios_base::out);

    if (!os.good())
        return false;

    os << 1 << endl << mat.size()-1 << endl;

    for(vector<Material>::const_iterator it=mat.begin()+1; it != mat.end(); ++it)
    {
        if (it->getg() == 1.0)
            os << "0.001 1 1 " << it->getn() << endl;
        else
            os << it->getMuA() << ' ' << it->getMuS() << ' ' << it->getg() << ' ' << it->getn() << endl;
    }

    if(mat[0].isMatched())
        os << 2 << endl;
    else
        os << 1 << endl << mat[0].getn() << endl;

    return os.good();
}

// Reads a TIM-OS results file (result.dat) - expects surface results first, then volume results
void readTIMOSOutput(string fn,const TetraMesh& mesh,SurfaceFluenceMap& surf,VolumeFluenceMap& vol)
{
    unsigned dtype,Ns,Nd;
    LineFile lf(fn,'%',cerr);

    surf.clear();
    vol.clear();

    while(!lf.eof())
    {
        lf >> dtype >> Nd >> Ns >> LineFile::LF_EOL;
        cout << "Reading data segment of type " << dtype << " with " << Nd << " points and " << Ns << " time steps" << endl;
        if (dtype == 1)
            readTIMOSOutputFlatten(mesh,surf,lf,Nd,Ns);
        else if (dtype == 2)
            readTIMOSOutputFlatten(mesh,vol,lf,Nd,Ns);
        else
            cerr << "Unexpected output format" << endl;
    }
}

// Read and flatten (sum across time steps for each element)
// Template parameter T specifies volume or surface element (either FaceByPointID or TetraByPointID)
template<class T>double readTIMOSOutputFlatten(const TetraMesh& mesh,FluenceMap<T>& F,LineFile& is,unsigned Nd,unsigned Ns)
{
    T IDps;
    F.clear();

    // read data lines and sum energy across all timesteps
    // IDp IDp IDp A f[Ns]
    double f, f_t, A, E_sum=0.0;
    for(unsigned i=0;i<Nd;++i)
    {
        is >> IDps >> A;
        if (!IDps.isSorted())
            cerr << "ERROR: Unsorted IDs in readTIMOSOutputFlatten" << endl;
        assert(IDps.isSorted());

        f=0;
        for(unsigned t=0;t<Ns;++t)
        {
            is >> f_t;
            f += f_t;
        }
        is >> LineFile::LF_EOL;
        E_sum += f*A;

        if (f != 0)
            F[IDps] = f;
    }
    return E_sum;
}
