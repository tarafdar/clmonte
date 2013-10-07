/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "graph.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/timer/timer.hpp>

#include <pthread.h>

#include "runresults.hpp"

#include "progress.hpp"

namespace globalopts {
    extern double Npkt;         // number of packets
    extern double prwin;        // Probability of winning roulette
    extern double wmin;         // Minimum weight for roulette
    extern long Nk;             // number of packets as long int
    extern unsigned Nthread;
    extern unsigned randseed;
    extern string logFN;        // log filename
}

// I am really sorry the worker/manager thread system here is SO ugly
// meant to be much more flexible and elegant
// needs to be fixed!!

template<class RNG,class Logger>boost::timer::cpu_times MonteCarloLoop(unsigned long long Np,Logger& logger,const TetraMesh& mesh,const vector<Material>&,const Source&);

class RunConfig;
template<class Logger,class RNG>class Manager;

template<class LoggerType,class RNG>int doOnePacket(const RunConfig& cfg,Packet pkt,
    LoggerType& logger,unsigned IDt,RNG& rng);


class RunConfig {
    public:
    const TetraMesh&        mesh;
    const vector<Material>& mat;
    const Source&           source;

    unsigned seed;
    unsigned Nthread;

    double wmin,pr_win;


    // default copy constructor is fine

    // prepare sources before creating run config
    RunConfig(const TetraMesh& mesh_,const vector<Material>& mat_,const Source& src_,unsigned seed_=1,unsigned Nthread_=1,
        double wmin_=globalopts::wmin,double pr_win_=globalopts::prwin) :
        mesh(mesh_),
        mat(mat_),
        source(src_),
        seed(seed_),
        Nthread(Nthread_),
        wmin(wmin_),
        pr_win(pr_win_)
        { }

    template<class LoggerType,class RNG>friend int doOnePacket(const RunConfig& cfg,Packet pkt,
        LoggerType& logger,unsigned IDt,RNG& rng);
};

template<class Logger,class RNG>class Worker;

template<class Logger,class RNG>class Manager {
    boost::timer::cpu_timer t;

    protected:
    const RunConfig cfg;
    unsigned long long Npacket;

    public:

    Manager(const RunConfig& cfg_,unsigned long long Npacket_) :
        cfg(cfg_),
        Npacket(Npacket_)
        {}

    virtual void run(Logger& logger)=0;
//    boost::timer::cpu_times start(Logger& logger);

    virtual unsigned long long getProgressCount() const { return 0; }

    virtual string getDetails() const { return ""; }

    double getProgressPercent() const  { return 100.0*double(getProgressCount())/double(Npacket); }
    double getRemainingPercent() const { return 100.0-getProgressPercent(); }
};

template<class Logger,class RNG>class Manager_MT : public Manager<Logger,RNG> {
    using Manager<Logger,RNG>::cfg;
    static void* threadStartFcn(void*);

    typedef Worker<Logger,RNG> WorkerType;
    typedef typename Logger::ThreadWorker LoggerType;

    pair<pthread_t,WorkerType*>* workers;

    public:
    virtual void run(Logger&);

    Manager_MT(const RunConfig& cfg_,unsigned long long Npacket_) :
        Manager<Logger,RNG>(cfg_,Npacket_),
        workers(new pair<pthread_t,WorkerType*>[cfg_.Nthread])
        {}

    unsigned long long getProgressCount() const {
        unsigned long long sum=0;
        for(unsigned i=0;i<cfg.Nthread && workers[i].second;++i)
            sum += workers[i].second->getProgressCount();
        return sum;
    }
    boost::timer::cpu_times start(Logger& logger);

    string getDetails() const {
        stringstream ss;
        for(unsigned i=0;i<cfg.Nthread && workers[i].second;++i)
            ss << workers[i].second->getProgressCount() << " ";
        return ss.str();
    }
};

template<class Logger,class RNG>void Manager_MT<Logger,RNG>::run(Logger& logger)
{
    typedef typename Logger::ThreadWorker ThreadWorker;
    ThreadWorker* l[cfg.Nthread];

    cout << "Running with " << cfg.Nthread << " threads" << endl;
    // create the workers
    for(unsigned i=0;i<cfg.Nthread;++i)
    {
        l[i] = new ThreadWorker(logger.getThreadWorkerInstance(i));
        workers[i].second = new WorkerType(Manager<Logger,RNG>::cfg,
            *(l[i]),
            Manager<Logger,RNG>::Npacket/cfg.Nthread,
            Manager<Logger,RNG>::cfg.seed+100*i,
            this);

        pthread_create(&workers[i].first,NULL,threadStartFcn,workers[i].second);
    }

    for(unsigned i=0;i<cfg.Nthread;++i)
    {
        pthread_join(workers[i].first,NULL);
        delete workers[i].second;
        delete l[i];
    }
    delete[] workers;
}

template<class Logger,class RNG>class Worker {
    public:
    const RunConfig cfg;
    unsigned long long i,Npacket;
    typedef typename Logger::ThreadWorker ThreadWorker;
    ThreadWorker& logger;
    RNG rng;
    Manager_MT<Logger,RNG>* manager;

    public:

    Worker(const RunConfig& cfg_,ThreadWorker& logger_,unsigned long long Npacket_,unsigned seed_,
        Manager_MT<Logger,RNG>* manager_) :
        cfg(cfg_),
        i(0),
        Npacket(Npacket_),
        logger(logger_),
        rng(1024,seed_),
        manager(manager_)
        {}

    unsigned long long getProgressCount() const { return i; }

    void start();
};

template<class Logger,class RNG>class MgrProgressUpdate {
    const Manager<Logger,RNG>* p;
    public:
    MgrProgressUpdate(const Manager<Logger,RNG>* p_) : p(p_){}

    void operator()() const
        { cout << "\rProgress " << setw(5) << fixed << setprecision(2) << p->getProgressPercent() << "%" << flush; };
};

template<class Logger,class RNG>boost::timer::cpu_times Manager_MT<Logger,RNG>::start(Logger& logger)
{
    boost::timer::cpu_timer runTimer;

    NewTimer<MgrProgressUpdate<Logger,RNG> > t(1.0,MgrProgressUpdate<Logger,RNG>(this),false);

    t.start();

    runTimer.start();
    run(logger);

    runTimer.stop();
    t.stop();
    cout << runTimer.format() << endl;
    return runTimer.elapsed();
}
/*
template<class Logger,class RNG>void Manager<Logger,RNG>::run(Logger& logger)
{
    Worker<Logger,RNG> w(cfg,logger,Npacket,cfg.seed,this);
    w.start();
}*/

template<class Logger,class RNG>void Worker<Logger,RNG>::start()
{
    pair<Packet,unsigned> tmp;
    Packet& pkt=tmp.first;
    unsigned& IDt=tmp.second;

    for(i=0;i<Npacket;++i)
    {
        tmp = cfg.source.emit(rng);
        doOnePacket(cfg,pkt,logger,IDt,rng);
    }
}

// thread launch function
template<class Logger,class RNG>void* Manager_MT<Logger,RNG>::threadStartFcn(void* arg)
{
    Worker<Logger,RNG>* w = (Worker<Logger,RNG>*) arg;
    pthread_t tid=pthread_self();
    w->start();
    return NULL;
}

template<class RNG,class Logger>boost::timer::cpu_times MonteCarloLoop(unsigned long long Np,Logger& logger,const TetraMesh& mesh,const vector<Material>& mat,const Source& src)
{
    RunConfig cfg(mesh,mat,src,globalopts::randseed,globalopts::Nthread,globalopts::wmin,globalopts::prwin);

    Manager_MT<Logger,RNG> mgr(cfg,Np);

    return mgr.start(logger);
}
