/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <iomanip>
#include "progress.hpp"

// NOTE: These switches are but partially implemented. You'll need to mess around with lines involving LoggerParentCons
// in runSimulation if you actually want to disable features; sorry, to be fixed!!
#define LOG_EVENT
#define LOG_VOLUME
#define LOG_SURFACE
//#define LOG_MEMTRACE
#define LOG_CONSERVATION

#include "logger.hpp"
#include "LoggerConservation.hpp"
#include "LoggerVolume.hpp"
#include "LoggerSurface.hpp"
#include "LoggerEvent.hpp"

#include "graph.hpp"
#include "source.hpp"
#include "random.hpp"
#include "io_timos.hpp"
#include <signal.h>
#include <boost/program_options.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/timer/timer.hpp>

#include "LoggerMemTrace.cpp"

#include "mainloop.cpp"
#include <map>

RunResults runSimulation(const TetraMesh& mesh,const vector<Material>& materials,const Source*,unsigned long long Nk);

namespace po=boost::program_options;

namespace globalopts {
    double Npkt=1e6;           // number of packets
    long Nk=0;                 // number of packets as long int
    unsigned Nthread=1;
    unsigned randseed=1;
    string outfn("output");     // output file name base
    double wmin=0.00001;
    double prwin=0.1;
}

using namespace std;

void banner()
{
    cout << "FullMonte v0.1" << endl;
    cout << "$Id: fullmonte.cpp 302 2013-03-03 05:40:12Z jcassidy $ (tree " << SVNVERSION << ')' << endl;
    cout << "(c) Jeffrey Cassidy, 2013" << endl;
    cout << endl;
}

int main(int argc,char **argv)
{
    // if you want to launch from a terminal and then close while it's still running
    // signal(SIGHUP,SIG_IGN);

    string fn_materials,fn_sources,fn_mesh;
    vector<Source*> sources;
    vector<Material> materials;

    banner();

    // define command-line options
    po::options_description cmdline("Command-line options");
    po::positional_options_description pos;
    pos.add("input",1).add("materials",1).add("sourcefile",1);

    cmdline.add_options()
        ("help,h","Display option help")
        ("input,i",po::value<string>(&fn_mesh),"Input file")
        ("N,N",po::value<double>(&globalopts::Npkt),"Number of packets")
        ("sourcefile,s",po::value<string>(&fn_sources),"Source location file (TIM-OS .source type)")
        ("materials,m",po::value<string>(&fn_materials),"Materials file (TIM-OS .opt type)")
        ("rngseed,r",po::value<unsigned>(&globalopts::randseed),"RNG seed (int)")
        ("threads,t",po::value<unsigned>(&globalopts::Nthread),"Thread count")
        ("output,o",po::value<string>(&globalopts::outfn),"Output file name base")
        ("wmin,w",po::value<double>(&globalopts::wmin),"Minimum weight for roulette")
        ("prwin,p",po::value<double>(&globalopts::prwin),"Probability of winning roulette")
        ("ignore-hup","Ignore SIGHUP (eg. if launching from terminal & then disconnecting)")
        ;

    // parse options (if an option is already set, subsequent store() calls will not overwrite
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc,argv).options(cmdline).positional(pos).run(),vm);
        po::notify(vm);
    }
    catch (po::error& e){
        cerr << "Caught an exception in options processing" << endl;
        cout << cmdline << endl;
        return -1;
    }

    // apply options
    if (vm.count("help"))
    {
        cout << cmdline << endl;
        return -1;
    }
	cout << endl;

    globalopts::Nk=globalopts::Npkt;

    TetraMesh* M;
    Source* src;

    if (vm.count("input"))
    {
        cout << "Loading mesh from " << fn_mesh << endl;
        M = new TetraMesh(fn_mesh,TetraMesh::MatlabTP);
    }

    if (vm.count("sourcefile"))
    {
        cout << "Loading sources from " << fn_sources << endl;
        vector<Source*>  sources = readTIMOSSource(fn_sources);         // Possible/probable memory leak here
        src = (sources.size() > 1 ? new SourceMulti(sources.begin(),sources.end()) : sources.front());
        src->prepare(*M);
    }

    if (vm.count("materials"))
    {
        cout << "Loading materials from " << fn_materials << endl;
        materials = readTIMOSMaterials(fn_materials);
    }

    cout << "Packets: " << bigIntSuffix(globalopts::Nk) << endl;

    if (!(materials.size() > 0 && src && M))
    {
        cerr << "Error: requires material definitions, sources, and mesh" << endl;
        exit(-1);
    }

    RunResults res = runSimulation(*M,materials,src,globalopts::Nk);

    return 0;
}

RunResults runSimulation(const TetraMesh& mesh,const vector<Material>& materials,const Source* source,unsigned long long Nk)
{

// DEBUG: Print materials and sources
//    cout << "Materials: " << endl;
//    unsigned i=0;
//    for(vector<Material>::const_iterator it=materials.begin(); it != materials.end(); ++it,++i)
//        cout << setw(2) << i << ": " << *it << endl;
//
//    cout << "Sources: " << endl;
//    cout << *source << endl;

    // Create logger types

#ifdef LOG_SURFACE
    LoggerSurfaceMT ls(mesh);
#endif
#ifdef LOG_VOLUME
    LoggerVolumeMT  lv(mesh);
#endif
#ifdef LOG_EVENT
    LoggerEventMT   le;
#endif
#ifdef LOG_MEMTRACE
    LoggerMemTraceMT lmt;
#endif
#ifdef LOG_CONSERVATION
    LoggerConservationMT lc;
#endif

    // Sorry this is REALLY ugly at the moment...
    //  intending to make it much more elegant
    LoggerParentCons<LoggerSurfaceMT,LoggerVolumeMT> l1(ls,lv);
    auto logger2=LoggerParentCons<LoggerParentCons<LoggerSurfaceMT,LoggerVolumeMT>,LoggerEventMT>(l1,le);
    LoggerParentCons<LoggerConservationMT,LoggerParentCons<LoggerParentCons<LoggerSurfaceMT,LoggerVolumeMT>,LoggerEventMT>> logger(lc,logger2);

    // Run it
    boost::timer::cpu_times t = MonteCarloLoop<RNG_SFMT>(Nk,logger,mesh,materials,*source);

    // Gather results
    RunResults res;
    res.Np = Nk;
    res.runid = 0;

    res.exitcode=0;
    res.t_wall=t.wall/1e9;
    res.t_user=t.user/1e9;
    res.t_system=t.system/1e9;

    // handle event logger
#ifdef LOG_CONSERVATION
    cout << "Conservation check: " << endl;
    cout << lc << endl << endl;
#endif

#ifdef LOG_EVENT
    cout << le << endl;
    res.Nintersection=le.Nbound;
    res.Nabsorb=le.Nabsorb;
    res.Nscatter=le.Nscatter;
    res.Ntir=le.Ntir;
    res.Nfresnel=le.Nfresnel;
    res.Nexit=le.Nexit;
    res.Nwin=le.Nwin;
    res.Nrefr=le.Nrefr;
    res.Ndie=le.Ndie;
#endif

#ifdef LOG_EVENT
    cout << "Total launched: " << le.Nlaunch << endl;
#else
    cout << "Total launched: " << Nk << endl;
#endif

#ifdef LOG_SURFACE
    SurfaceFluenceMap surf(&mesh);
	ls.fluenceMap(surf);
    surf.writeText(globalopts::outfn+".surf.txt",materials);
    cout << "Total exited:   " << surf.getTotalEnergy() << endl;
#endif

#ifdef LOG_VOLUME
    VolumeFluenceMap vol(&mesh);
    lv.fluenceMap(vol,materials);
    vol.writeText(globalopts::outfn+".vol.txt",materials);
    cout << "Total absorbed: " << vol.getTotalEnergy() << endl;
#endif

    return res;
}
