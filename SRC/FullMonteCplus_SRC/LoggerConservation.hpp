/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "logger.hpp"

// LoggerConservation - checks that energy is conserved by logging total energy:
//  launched
//  absorbed
//  expired through roulette (w_die)
//  exiting volume
//  added by roulette wins (w_roulette)

class LoggerConservation : public LoggerNull {
    double w_launch,w_absorb,w_die,w_exit,w_roulette;

    public:
    inline void eventLaunch(Ray3 r,unsigned IDt,double w) { w_launch += w; };
    inline void eventAbsorb(Point3 p,unsigned IDt,double w0,double dw) { w_absorb += dw; };
    inline void eventExit(Ray3 r,int IDf,double w) { w_exit += w; };
    inline void eventDie(double w){ w_die += w; };
    inline void eventRouletteWin(double w0,double w){ w_roulette += w-w0; };

    LoggerConservation() : w_launch(0.0),w_absorb(0.0),w_die(0.0),w_exit(0.0),w_roulette(0.0){}
    LoggerConservation(string fn) { LoggerConservation(); };

    LoggerConservation& operator+=(const LoggerConservation& lc){
        w_launch += lc.w_launch;
        w_absorb += lc.w_absorb;
        w_die    += lc.w_die;
        w_exit   += lc.w_exit;
        w_roulette += lc.w_roulette;
        return *this;
    }

    void clear(){ *this=LoggerConservation(); }

    friend ostream& operator<<(ostream&,const LoggerConservation&);
};


// LoggerConservationMT: Multithreaded instance of conservation logger
//  just keeps one set of counters for each thread, then merges through the += operator

class LoggerConservationMT : public LoggerConservation,private boost::mutex {
    public:
    class ThreadWorker : public LoggerConservation {
        LoggerConservationMT& parent;
        public:
        ThreadWorker(LoggerConservationMT& parent_) : parent(parent_){}
        ~ThreadWorker(){ commit(); }

        void commit()
        {
            parent.lock();
            parent += *this;
            parent.unlock();
            clear();
        }
    };
    ThreadWorker getThreadWorkerInstance(unsigned) { return ThreadWorker(*this); }
};
