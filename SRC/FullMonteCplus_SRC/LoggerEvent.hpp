/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "logger.hpp"

class LoggerEvent {
    public:
    unsigned long long Nlaunch,Nabsorb,Nscatter,Nbound,Ntir,Nfresnel,Nrefr,Ninterface,Nexit,Ndie,Nwin;

    void clear(){ *this = LoggerEvent(); }

    LoggerEvent() : Nlaunch(0),Nabsorb(0),Nscatter(0),Nbound(0),Ntir(0),Nfresnel(0),Nrefr(0),Ninterface(0),
        Nexit(0),Ndie(0),Nwin(0){};

    const LoggerEvent& operator+=(const LoggerEvent& rhs);
    friend ostream& operator<<(ostream&,const LoggerEvent&);
};

class LoggerEventST : public LoggerEvent,public LoggerNull {
    public:
    inline void eventLaunch(const Ray3 r,unsigned IDt,double w){ ++Nlaunch; };   // launch new packet

    inline void eventAbsorb(const Point3 p,unsigned IDt,double w0,double dw){ ++Nabsorb; };     // packet absorbed
    inline void eventScatter(const UVect3 d0,const UVect3 d,double g){ ++Nscatter; };

    inline void eventBoundary(const Point3 p,int,int,int){ ++Nbound; };   // boundary (same material)
    inline void eventInterface(const Ray3,int,unsigned){ ++Ninterface; }; // found a material interface; possible results are:
    inline void eventRefract(const Point3,UVect3){ ++Nrefr; };                //      refracted
    inline void eventReflectInternal(const Point3,const UVect3){ ++Ntir; };        //      internal reflection
    inline void eventReflectFresnel(const Point3,UVect3){ ++Nfresnel; };         //      fresnel reflection

    // termination events
    inline void eventExit(const Ray3,int,double){ ++Nexit; };        // exited geometry
    inline void eventDie(double){ ++Ndie; };                                    // lost Russian roulette
    inline void eventRouletteWin(double,double){ ++Nwin; };                     // won roulette
};

class LoggerEventMT : public LoggerEvent, private boost::mutex {
    public:

    // each worker has a copy of the event counter, a reference to the parent, and a commit() method
    class ThreadWorker : public LoggerEventST {
        LoggerEventMT& parent;
        public:
        ThreadWorker(LoggerEventMT& parent_) : parent(parent_){}
        ~ThreadWorker() { commit(); }

        void commit()
        {
            parent.lock();
            parent += *this;
            clear();
            parent.unlock();
        }
    };

    ThreadWorker getThreadWorkerInstance(unsigned){ return ThreadWorker(*this); }
};
