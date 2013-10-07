/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include <iostream>
#include <fstream>
#include "logger.hpp"

using namespace std;

typedef struct {
    unsigned IDt;
    unsigned Nabs;
} TetraRecord;

typedef struct {
    int IDf;
    unsigned Nexit;
} ExitRecord;

class LoggerMemTrace : public LoggerNull {
    string fnTetra,fnExit;

    ofstream os_tetra,os_exit;

    unsigned IDt_current,IDt_next;
    int IDf_last_exit;
    unsigned Nabs,Nexit;

    // helper functions
    void logTetraHits(unsigned,unsigned);
    void logExit(int,unsigned);

    public:
    LoggerMemTrace(string fnTetra_,string fnExit_) : fnTetra(fnTetra_),fnExit(fnExit_),os_tetra(fnTetra_),os_exit(fnExit_),
        IDf_last_exit(0),Nabs(-1),Nexit(-1){
            if (!os_tetra.good())
                cerr << "Failed to open " << fnTetra_ << " for writing" << endl;
            if (!os_exit.good())
                cerr << "Failed to open " << fnExit_ << " for writing" << endl; }
    LoggerMemTrace(const LoggerMemTrace& l_) : fnTetra(l_.fnTetra),fnExit(l_.fnExit),os_tetra(l_.fnTetra),os_exit(l_.fnExit),
        IDf_last_exit(0),Nabs(-1),Nexit(-1){
            if (!os_tetra.good())
                cerr << "Failed to open " << l_.fnTetra << " for writing" << endl;
            if (!os_exit.good())
                cerr << "Failed to open " << l_.fnExit << " for writing" << endl; }
    ~LoggerMemTrace();

    void eventLaunch(const Ray3 r,unsigned IDt,double w);
    void eventAbsorb(const Point3 p,unsigned IDt,double w0,double dw);
    void eventBoundary(const Point3 p,int,unsigned,unsigned);
    void eventInterface(const Ray3,int,unsigned);
    void eventRefract(const Point3,UVect3);

    // termination events
    void eventExit(const Ray3,int,double);
};

class LoggerMemTraceMT {
    public:

    // ThreadWorker is just a single-thread instance but with a different file name
    class ThreadWorker : public LoggerMemTrace {
        public:
        ThreadWorker(string a,string b) : LoggerMemTrace(a,b){}
    };

    // thread ## writes to files tetra.trace.##.bin and exit.trace.##.bin
    ThreadWorker getThreadWorkerInstance(unsigned i){
        stringstream ss_tetra,ss_exit;
        ss_tetra << "tetra.trace." << i << ".bin";
        ss_exit  << "exit.trace."  << i << ".bin";
        return ThreadWorker(ss_tetra.str(),ss_exit.str());
    }
};
