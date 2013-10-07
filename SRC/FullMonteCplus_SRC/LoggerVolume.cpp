/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "LoggerVolume.hpp"

ostream& operator<<(ostream& os,const LoggerVolume&){ cout << "Volume logger report" << endl; return os; }

void LoggerVolume::fluenceMap(VolumeFluenceMap& F,const vector<Material>& mat,bool per_volume)
{
    F.clear();
    map<unsigned,double>::const_iterator m_it=F.begin();
    unsigned i=1;

    for(vector<FluenceCountType>::const_iterator v_it=counts.begin()+1; v_it != counts.end(); ++v_it,++i)
    {
        if (getValue(*v_it) != 0)
            F[i] = getValue(*v_it);
    }
}

void LoggerVolume::hitMap(map<unsigned,unsigned long long>& m)
{
    m.clear();
    map<unsigned,unsigned long long>::iterator m_it=m.begin();
    unsigned IDt=1;
    if (counts.size()<2)
        return;
    for(vector<FluenceCountType>::const_iterator it=counts.begin()+1; it != counts.end(); ++it,++IDt)
        if (getHits(*it) != 0)
            m_it = m.insert(m_it,make_pair(IDt,getHits(*it)));
}
