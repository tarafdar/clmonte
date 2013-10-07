/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "logger.hpp"

// LoggerVolume tracks the energy absorbed in each tetrahedral mesh element

class LoggerVolume {
    const TetraMesh& mesh;

    protected:
    vector<FluenceCountType> counts;

    public:
    LoggerVolume(LoggerVolume&& lv_) : mesh(lv_.mesh),counts(std::move(lv_.counts)){};
    LoggerVolume(const TetraMesh& mesh_) : mesh(mesh_),counts(mesh_.getNt()+1){};
    
    // get results as a map; optional arg per_area specifies to return fluence (per-area) or total energy per patch
    void fluenceMap(VolumeFluenceMap&,const vector<Material>&,bool=true);
    void hitMap(map<unsigned,unsigned long long>& m);

    friend ostream& operator<<(ostream&,LoggerVolume&);
};


// LoggerVolumeST (single-thread) overrides eventAbsorb
class LoggerVolumeST : public LoggerVolume,public LoggerNull {
    public:
    LoggerVolumeST(const TetraMesh& mesh_) : LoggerVolume(mesh_){}

    // log only absorption events
    inline void eventAbsorb(const Point3,unsigned IDt,double w0,double dw){ counts[IDt] += dw; }
};

// LoggerVolumeMT (multi-thread) provides access to elements of type ThreadWorker
//  Each threadworker has a buffer of 1M (1024*1024) absorption events
//  When the buffer is full, it locks a mutex on the absorption map, updates, and unlocks

class LoggerVolumeMT : public LoggerVolume,private boost::mutex {
    typedef pair<unsigned,FluenceCountType> BufElType;
    unsigned bufsize;

    public:
    LoggerVolumeMT(const TetraMesh& mesh_,unsigned bufsize_=1024*1024) : LoggerVolume(mesh_),bufsize(bufsize_){}

    class ThreadWorker : private Buffer<BufElType>,public LoggerNull {
	    using Buffer<BufElType>::atBufferEnd;
        LoggerVolumeMT& parent;
        unsigned ID_last;
        virtual void atBufferEnd(const BufElType* begin,const BufElType* end){ parent.addValues(begin,end); ID_last=-1; }

        public:

        ThreadWorker(LoggerVolumeMT& parent_,unsigned N_) : Buffer<BufElType>(N_,false),parent(parent_),ID_last(-1){};
        ThreadWorker(ThreadWorker&& tw_) : Buffer<BufElType>(std::move(tw_)),parent(tw_.parent),ID_last(tw_.ID_last){}
        ~ThreadWorker() { flush(); }

        inline void eventAbsorb(Point3 p,unsigned IDt,double w0,double dw)
        {
            if (ID_last == IDt)
				current->second += dw;
            else
            {
                BufElType *tmp = ((ID_last == (unsigned)(-1)) ? current : getNext());
                tmp->first=ID_last=IDt;
                tmp->second=dw;
            }
        };

        void flush()
        {
            atBufferEnd(first,current);
            current=first;
        }

        void commit()
        {
            flush();
        }
    };

    // create a new thread worker with its own buffer
    ThreadWorker getThreadWorkerInstance(unsigned)
    {
        return ThreadWorker(*this,bufsize);
    }

    // adds values to the log
    void addValues(const pair<unsigned,FluenceCountType>* p,const pair<unsigned,FluenceCountType>* p_last)
    {
        lock();
        for(; p != p_last; ++p)
	        counts[p->first] += p->second;
        unlock();
    }
};

ostream& operator<<(ostream&,LoggerVolume&);
