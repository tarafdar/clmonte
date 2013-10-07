/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include "progress.hpp"
#include <signal.h>
#include <cmath>

// Macro POSIX_TIMER must be defined if POSIX timer is present (Mac OS X does not have)
// If not, all code creating timers is commented out

using namespace std;

// prints as "       1/1000000 (0.01%)"
ostream& operator<<(ostream& os,const ProgressCountdown& pc)
{
    return os << setfill(' ') << setw(8) << (pc.start-*pc.remain) << "/" << setw(8) << pc.start
        << " (" << setw(5) << fixed << setprecision(1) << pc.getCompletedPercent() << "%) est "
        << setw(8) << setprecision(1) << pc.getElapsedTime()*pc.getRemainingRatio() << " sec left";
}

#ifdef POSIX_TIMER
struct itimerspec POSIX_Timer::tval_zero = {{0,0},{0,0}};
#endif

POSIX_Timer::POSIX_Timer(double interval_,bool auto_start)
{
#ifdef POSIX_TIMER
    struct sigevent sigdef;
    sigdef.sigev_notify=SIGEV_THREAD;
    sigdef.sigev_signo=0; // not relevant?
    sigdef.sigev_value.sival_ptr=this;
    sigdef.sigev_notify_function=timer_function;
    sigdef.sigev_notify_attributes=NULL;

    if(timer_create(CLOCK_REALTIME,&sigdef,&timerid))
        cerr << "Failed to create timer" << endl;
    else {
        setInterval(interval_);
        if(auto_start)
            start();
    }
#endif
}

void POSIX_Timer::setInterval(double interval_)
{
    interval=interval_;
#ifdef POSIX_TIMER
    timerval.it_interval.tv_sec  = timerval.it_value.tv_sec  = floor(interval_);
    timerval.it_interval.tv_nsec = timerval.it_value.tv_nsec = interval_-floor(interval_);
#endif
}

void POSIX_Timer::start()
{
#ifdef POSIX_TIMER
    if (timer_settime(timerid,0,&timerval,NULL))
        cerr << "Failed to set timer" << endl;
#endif
    tick();
}

void POSIX_Timer::pause()
{
#ifdef POSIX_TIMER
    timer_settime(timerid,0,&timerval,NULL);
#endif
}

void POSIX_Timer::timer_function(union sigval sigev_value)
{
    class POSIX_Timer* prog=(POSIX_Timer*)(sigev_value.sival_ptr);
    ++prog->Ntick;
    prog->tick();
}

void POSIX_Timer::stop()
{
#ifdef POSIX_TIMER
    if(timer_settime(timerid,0,&tval_zero,NULL))
        cerr << "Failed to set timer value in stop()" << endl;
#endif
}

POSIX_Timer::~POSIX_Timer()
{
#ifdef POSIX_TIMER
    if(timer_delete(timerid))
        cerr << "Failed to delete timer" << endl;
#endif
}

// prints a message:
//  ######## completed after #.##s; throughput is #######/sec

void ProgressPrinter::operator()() const
{
    double t=getElapsedTime();
    unsigned c=*counter;
    if (w != -1)
        os << setw(w);
    os << c << " completed after " << fixed << setprecision(2) << t << "s; throughput is "
        << scientific << (double(c)/t) << "/sec" << endl;
}

