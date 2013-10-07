/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "LoggerMemTrace.hpp"

LoggerMemTrace::~LoggerMemTrace()
{
    os_tetra.close();
    os_exit.close();
}

// On launch, write out previous tetra record unless it's the same tetra
void LoggerMemTrace::eventLaunch(const Ray3,unsigned IDt,double)
{
    if (Nabs == (unsigned)(-1))
        Nabs=0;
    else if (IDt_current != IDt)
    {
        logTetraHits(IDt_current,Nabs);
        Nabs=0;
    }

    IDt_current=IDt;
}

// At a boundary (no refractive index difference or matched boundary), write out previous tetra record
void LoggerMemTrace::eventBoundary(const Point3,int,unsigned,unsigned IDte)
{
    logTetraHits(IDt_current,Nabs);
    IDt_current=IDte;
    Nabs=0;
}

// At interface, mark the next IDt but don't score it yet (we may reflect back into current tetra)
void LoggerMemTrace::eventInterface(const Ray3,int,unsigned IDt)
{
    IDt_next=IDt;
}

// Record an absorption - IDt_current is always correct
void LoggerMemTrace::eventAbsorb(const Point3 p,unsigned IDt,double,double)
{
    if (IDt != IDt_current)
        cerr << "LoggerMemTrace::eventAbsorb - Oops! IDt_current is not correct!" << endl;
    ++Nabs;
}

// In case of refraction, we have entered the new tetra recorded at the interface
void LoggerMemTrace::eventRefract(const Point3,UVect3)
{
    logTetraHits(IDt_current,Nabs);
    Nabs=0;
    IDt_current=IDt_next;
}

// Record exit
void LoggerMemTrace::eventExit(const Ray3,int IDf,double)
{
    if (IDf==IDf_last_exit)
        ++Nexit;
    else {
        if (IDf_last_exit != 0)
            logExit(IDf_last_exit,Nexit);
        Nexit=1;
        IDf_last_exit=IDf;
    }
}

void LoggerMemTrace::logTetraHits(unsigned IDt,unsigned Nabs)
{
    TetraRecord tmp = { IDt,Nabs} ;
    os_tetra.write((const char*)&tmp,sizeof(TetraRecord));
    Nabs=0;
}

void LoggerMemTrace::logExit(int IDf,unsigned Nexit)
{
    ExitRecord tmp = {abs(IDf),Nexit};
    os_exit.write((const char*)&tmp,sizeof(ExitRecord));
    Nexit=0;
}
