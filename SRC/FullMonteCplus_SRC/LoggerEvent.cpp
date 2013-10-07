/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "LoggerEvent.hpp"

const LoggerEvent& LoggerEvent::operator+=(const LoggerEvent& rhs)
{
    Nlaunch     += rhs.Nlaunch;
    Nabsorb     += rhs.Nabsorb;
    Nscatter    += rhs.Nscatter;
    Nbound      += rhs.Nbound;
    Ntir        += rhs.Ntir;
    Nfresnel    += rhs.Nfresnel;
    Nrefr       += rhs.Nrefr;
    Ninterface  += rhs.Ninterface;
    Nexit       += rhs.Nexit;
    Ndie        += rhs.Ndie;
    Nwin        += rhs.Nwin;
    return *this;
}

ostream& operator<<(ostream& os,const LoggerEvent& le)
{
    os << "Launched: " << le.Nlaunch << endl;

    os << "Boundary (same):      " << le.Nbound << endl;
    os << "Boundary (different): " << le.Ninterface << endl;
    os << "  TIR:     " << le.Ntir << endl;
    os << "  Fresnel: " << le.Nfresnel << endl;
    os << "  Refract: " << le.Nrefr << endl;
    os << "  Balance (bound - [TIR + fresnel + refract]): " << le.Ninterface-le.Ntir-le.Nfresnel-le.Nrefr << endl;

    os << "Absorption: " << le.Nabsorb << endl;
    os << "Scatter:    " << le.Nscatter << endl;

    os << "Roulette results" << endl;
    os << "  Win:  " << le.Nwin << endl;
    os << "  Lose: " << le.Ndie << endl;

    os << "End results" << endl;
    os << "Died:   " << le.Ndie << endl;
    os << "Exited: " << le.Nexit << endl;
    os << "Balance ([launch] - [die + exit]): " << le.Nlaunch-le.Ndie-le.Nexit << endl;

    return os;
}
