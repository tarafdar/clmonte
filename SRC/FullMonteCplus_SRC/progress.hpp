/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef PROGRESS_INCLUDED
#define PROGRESS_INCLUDED
#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <signal.h>
#include <sys/types.h>

#ifdef PLATFORM_DARWIN
#undef POSIX_TIMER
#else
#define POSIX_TIMER
#endif

using namespace std;

// Provides a progress counter using the POSIX timer facilities

class Progress {
    struct timeval t_start;
    public:
    Progress(){ gettimeofday(&t_start,NULL); }

    double getElapsedTime() const {
        struct timeval t_last;
        gettimeofday(&t_last,NULL);
        return (t_last.tv_sec-t_start.tv_sec)+(t_last.tv_usec-t_start.tv_usec)*1e-6;
    }
};

class ProgressPrinter : public Progress {
    ostream& os;
    const unsigned* counter;
    int w;
    public:

    ProgressPrinter(ostream& os_,const unsigned* counter_,int w_=-1) : os(os_),counter(counter_),w(w_){}

    void operator()() const;
};

class ProgressCountdown : public Progress {
    unsigned long long start,*remain;
    public:
    ProgressCountdown(unsigned long long* remain_) : start(*remain_),remain(remain_){}

    // inspect remaining count
    unsigned long long getRemainCount() const { return *remain; }
    unsigned long long getStartCount() const { return start; }
    unsigned long long getCompletedCount() const { return start-*remain; }

    // calculate as percentage
    double getCompletedPercent() const { return 100.0*(start-*remain)/start; }
    double getRemainingPercent() const { return 100.0*(*remain)/start; }

    // get remaining effort as a ratio to the amount done
    double getRemainingRatio()      const { return 1.0*(*remain)/(start-*remain); }

    friend ostream& operator<<(ostream&,const ProgressCountdown&);
};

class POSIX_Timer {
    double interval;
    unsigned long long Ntick;
#ifdef POSIX_TIMER
    struct itimerspec timerval;
    timer_t timerid;
#endif

    // the actual timer function (updates time and calls tick() )
    static void timer_function(union sigval sigev_value);

    static struct itimerspec tval_zero;

    public:

    POSIX_Timer(double interval_,bool auto_start=false);
    ~POSIX_Timer();

    void setInterval(double interval_);

    void start();
    void start(double interval_){ setInterval(interval_); start(); }
    void pause();
    void stop();

    virtual void tick(){ cout << "ping "; }
};

template<class Callable>class NewTimer : public POSIX_Timer {
    Callable f;
    public:
    NewTimer(double interval_,Callable f_,bool auto_start=false) :
        POSIX_Timer(interval_,auto_start),
        f(f_){}
    virtual void tick() { f(); }
};

#endif
