/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef RUNRESULTS_INCLUDED
#define RUNRESULTS_INCLUDED

#include "blob.hpp"

// Struct to hold the results of a run

class RunResults {
    public:
    RunResults() : t_wall(0),t_user(0),t_system(0),E_launch(0),E_exit(0),E_absorb(0),E_die(0),runid(0),exitcode(0),Np(0),
        Nintersection(0),Nscatter(0),Nabsorb(0),Ntir(0),Nfresnel(0),Nexit(0),Nwin(0),Nrefr(0),Ndie(0){}
    double t_wall,t_user,t_system;
    double E_launch,E_exit,E_absorb,E_die;
    unsigned runid,exitcode;

    unsigned long long Np,Nintersection,Nhop,Nscatter,Nabsorb,Ntir,Nfresnel,Nexit,Nwin,Nrefr,Ndie;

    Blob log_stdout;
    Blob log_stderr;
};

#endif
