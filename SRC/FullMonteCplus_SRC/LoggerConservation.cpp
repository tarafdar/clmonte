/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "LoggerConservation.hpp"

ostream& operator<<(ostream& os,const LoggerConservation& log)
{
    double w_dispose = log.w_absorb+log.w_exit;
    double w_diff = w_dispose-log.w_launch;

    os << "Energy conservation" << endl;
    os << "  Launched: " << log.w_launch << endl;
    os << "  Disposed: " << w_dispose << endl;
    os << "    Absorbed " << log.w_absorb << endl;
    os << "    Exited   " << log.w_exit << endl;
    os << "  Difference: " << w_diff << " (" << 100.0*w_diff/log.w_launch << "%)" << endl;
    os << endl;

    os << "Roulette difference " << log.w_roulette-log.w_die << endl;
    os << "  Died " << log.w_die << endl;
    os << "  Added " << log.w_roulette << endl;
    return os;
}


