/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef LOGGER_INCLUDED
#define LOGGER_INCLUDED
#include <iostream>
#include <fstream>
#include <vector>
#include "graph.hpp"
#include "newgeom.hpp"
#include "fluencemap.hpp"
#include "hitcounter.hpp"
#include <boost/thread.hpp>
#include <boost/tuple/tuple.hpp>

// abstract Logger class, with member functions for all relevant events
// Logger is not required or assumed to be reentrant (SINGLE-THREAD ONLY)

typedef __m128 Point3;
typedef __m128 UVect3;
typedef pair<__m128,__m128> Ray3;

// If NOHIT is defined, just count fluence (faster); if not defined, will also count memory accesses

#ifdef NOHIT
typedef double FluenceCountType;
#else
typedef HitCount<unsigned,double> FluenceCountType;
#endif

class LoggerEvent;
class LoggerVolume;
ostream& operator<<(ostream&,const LoggerEvent&);
ostream& operator<<(ostream&,const LoggerVolume&);

// A single-thread Logger must have:
//  move constructor
//  result type
//  += operator (for accumulating multi-thread results together; possibly no-op if log methods already do this)
//  methods overriding LoggerNull to handle events

// A multithread Logger must have:
//  a type ThreadWorker which meets requirements of Logger
//  a member getThreadWorkerInstance(unsigned threadNum) to return a ThreadWorker

template<class H,class T>class LoggerCons;

// LoggerParentCons gives the ability to group a set of loggers and get thread workers for the whole set
template<class H,class T>class LoggerParentCons {
    H& head;
    T& tail;
    public:
    typedef H head_type;
    typedef T tail_type;
    typedef LoggerCons<typename H::ThreadWorker,typename T::ThreadWorker> ThreadWorker;
    LoggerParentCons(H& head_,T& tail_) : head(head_),tail(tail_){}

    ThreadWorker getThreadWorkerInstance(unsigned i)
        { return ThreadWorker(head.getThreadWorkerInstance(i),
            tail.getThreadWorkerInstance(i)); }
};

template<class H,class T>class LoggerCons {
    H head;
    T tail;

    public:
    typedef H head_type;
    typedef T tail_type;
    typedef LoggerCons<typename H::ThreadWorker,typename tail_type::ThreadWorker> ThreadWorker;

    LoggerCons(LoggerCons&& l_) : head(std::move(l_.head)),tail(std::move(l_.tail)){}
    LoggerCons(H&& head_,T&& tail_)  : head(std::move(head_)),tail(std::move(tail_)){}

    // event handlers
    inline void eventLaunch(const Ray3 r,unsigned IDt,double w)
        { head.eventLaunch(r,IDt,w); tail.eventLaunch(r,IDt,w); };

    inline void eventAbsorb(const Point3 p,unsigned IDt,double w0,double dw)
        { head.eventAbsorb(p,IDt,w0,dw); tail.eventAbsorb(p,IDt,w0,dw); };

    inline void eventScatter(const UVect3 d0,const UVect3 d,double g)
        { head.eventScatter(d0,d,g); tail.eventScatter(d0,d,g); };

    inline void eventBoundary(const Point3 p,int IDf,int IDts,int IDte)        // boundary (same material)
        { head.eventBoundary(p,IDf,IDts,IDte); tail.eventBoundary(p,IDf,IDts,IDte); }

    inline void eventInterface(const Ray3 r,int IDf,unsigned a)          // found a material interface; possible results are:
        { head.eventInterface(r,IDf,a); tail.eventInterface(r,IDf,a); }

    inline void eventRefract(const Point3 p,UVect3 d)                //      refracted
        { head.eventRefract(p,d); tail.eventRefract(p,d); }

    inline void eventReflectInternal(const Point3 p,const UVect3 d)  //      internal reflection
        { head.eventReflectInternal(p,d); tail.eventReflectInternal(p,d); }

    inline void eventReflectFresnel(const Point3 p,UVect3 d)         //      fresnel reflection
        { head.eventReflectFresnel(p,d); tail.eventReflectFresnel(p,d); }

    // termination events
    inline void eventExit(const Ray3 r,int IDf,double w)            // exited geometry
        { head.eventExit(r,IDf,w); tail.eventExit(r,IDf,w); }

    inline void eventDie(double w)                                     // lost Russian roulette
        { head.eventDie(w); tail.eventDie(w); }

    inline void eventRouletteWin(double w0,double w)                      // won roulette
        { head.eventRouletteWin(w0,w); tail.eventRouletteWin(w0,w); }

    const LoggerCons& operator+=(const LoggerCons& rhs)
        { head += rhs.head; tail += rhs.tail; return *this; }

    template<class H0,class T0>friend ostream& operator<<(ostream&,const LoggerCons<H0,T0>& l);
};

// make_logger functions
/*template<class L0>L0 make_logger(L0&& l0){ return l0; }
template<class L0,class L1>LoggerParentCons<L0,L1> make_logger(L0&& l0,L1&& l1){ return LoggerParentCons<L0,L1>(l0,l1); }
template<class L0,class L1,class L2>LoggerParentCons<L0,LoggerParentCons<L1,L2>> make_logger(L0&& l0,L1&& l1,L2&& l2)
{
    return LoggerParentCons<L0,LoggerParentCons<L1,L2>>(std::move(l0),make_logger<L1,L2>(l1,l2));
}
template<class L0,class L1,class L2,class L3>LoggerParentCons<L0,LoggerParentCons<L1,LoggerParentCons<L2,L3>>> make_logger(L0&& l0,L1&& l1,L2&& l2,L3& l3)
{
    return LoggerParentCons<L0,LoggerParentCons<L1,LoggerParentCons<L2,L3>>>(std::move(l0),make_logger<L1,L2,L3>(l1,l2,l3));
}*/

template<class H,class T> ostream& operator<<(ostream& os,const LoggerCons<H,T>& l)
{
    return os << l.head << endl << l.tail << endl;
}

class LoggerNull {
    public:
    virtual ~LoggerNull(){}

    inline void eventLaunch(const Ray3 r,unsigned IDt,double w){};

    inline void eventAbsorb(const Point3 p,unsigned IDt,double w0,double dw){};
    inline void eventScatter(const UVect3 d0,const UVect3 d,double g){};

    inline void eventBoundary(const Point3 p,int,int,int){};        // boundary (same material)

    inline void eventInterface(const Ray3,int,unsigned){};          // found a material interface; possible results are:
    inline void eventRefract(const Point3,UVect3){};                //      refracted
    inline void eventReflectInternal(const Point3,const UVect3){};  //      internal reflection
    inline void eventReflectFresnel(const Point3,UVect3){};         //      fresnel reflection

    // termination events
    inline void eventExit(const Ray3,int,double){};            // exited geometry
    inline void eventDie(double){};                                     // lost Russian roulette
    inline void eventRouletteWin(double,double){};                      // won roulette

    typedef LoggerNull ThreadWorker;

    ThreadWorker getThreadWorkerInstance(unsigned) const { return ThreadWorker(); }

    const LoggerNull& operator+=(const LoggerNull&){ return *this; };
};



// Buffer is a support class

template<class T>class Buffer {
    protected:
    T*          first;
    T*          last;
    T*          current;    // always points to a valid place [first.get() .. last-1]

    unsigned N;

    virtual void atBufferEnd(const T*,const T*)=0;

    public:

    Buffer(unsigned N_,bool doInit_) :
        first(new T[N_]),
        last(first+N_),
        current(doInit_ ? last-1 : first),
        N(N_){}

    Buffer(Buffer&& b_) : first(b_.first),last(b_.last),current(b_.current),N(b_.N){
        b_.first = b_.current = NULL;
        b_.N = 0;
    }

    virtual ~Buffer() { delete[] first; }

	void flush() { atBufferEnd(first,current); current=first; }

    // get next available slot; if at last slot, flush buffer and go back to start
    T* getNext(){
        if (current == last-1)
        {
            atBufferEnd(first,last);
            current=first;
        }
		else
			++current;
        return current;
    }
};


#endif
