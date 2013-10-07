/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "logger.hpp"

class LoggerSurface {
    const TetraMesh& mesh;

//    protected:
    public:
    vector<FluenceCountType> counts;

    public:
    LoggerSurface(const TetraMesh& mesh_) : mesh(mesh_),counts(mesh_.getNf()+1,FluenceCountType()){ };
    
    // get results as a map; optional arg per_area specifies to return fluence (per-area) or total energy per patch
//    void resultMap(map<FaceByPointID,double>& m,bool per_area=true);
    void fluenceMap(SurfaceFluenceMap&);

    void hitMap(map<unsigned,unsigned long long>& m);
};


// LoggerSurfaceST (single-thread version) overrides only the eventExit method to record packets exiting the volume
class LoggerSurfaceST : public LoggerSurface,public LoggerNull {
    public:
    LoggerSurfaceST(const TetraMesh& mesh_) : LoggerSurface(mesh_){}

    inline void eventExit(const Ray3,int IDf,double w){
        IDf=abs(IDf);
		counts[IDf] += w;
    }
};

// LoggerSurfaceMT (multi-thread) keeps just one record with a mutex to protect it
//  This isn't the best solution since locking/unlocking is potentially expensive
//  However we exit at most once per packet launched vs. ~100-400 absorption events per pkt

class LoggerSurfaceMT : public LoggerSurfaceST,private boost::mutex {
    public:
    LoggerSurfaceMT(const TetraMesh& mesh_) : LoggerSurfaceST(mesh_){ }

    class ThreadWorker : public LoggerNull {
        LoggerSurfaceMT& parent;

        public:
        ThreadWorker(LoggerSurfaceMT& parent_) : parent(parent_){ }
        ThreadWorker(const ThreadWorker& tw_) : parent(tw_.parent){}
        
        inline void commit(){}
        inline void eventExit(const Ray3 r,int IDf,double w){
            parent.lock();
            parent.eventExit(r,IDf,w);
            parent.unlock();
        }
    };

    // Since we protect each access with a mutex, thread worker is just a reference to the LoggerSurfaceMT object
    ThreadWorker getThreadWorkerInstance(unsigned){ return ThreadWorker(*this); }
};
