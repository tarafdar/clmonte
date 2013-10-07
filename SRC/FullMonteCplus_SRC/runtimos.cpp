/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#define DB_DEF_SIMULATOR 1
#include <iostream>
#include <string>
//#include <pthread.h>

#include <algorithm>
#include <signal.h>

#include <boost/scoped_array.hpp>

#include <boost/timer/timer.hpp>

#include <spawn.h>

#include "graph.hpp"
#include "io_timos.hpp"

#include <unistd.h>
#include <sys/wait.h>

#include "source.hpp"

#include "fm-postgres/fm-postgres.hpp"

#include "fmdb.hpp"

#include <sys/stat.h> // for stat/creat
#include <fcntl.h>
#include <sys/types.h>

#ifdef PLATFORM_DARWIN
extern char** environ;
#endif

const char timos_bin[] = {"/home/jcassidy/src/refcode/TIMOS_Package/source/timos"};

TetraMesh* exportMesh(PGConnection&,unsigned);
int exportSources(PGConnection&,unsigned,vector<Source*>&,long long Npacket=0);
int exportMaterials(PGConnection&,unsigned,vector<Material>&);

using namespace std;
int runtimos(string path,string meshfn,string sourcefn,string matfn,string outputfn,string logfn,string errfn,int Nthread,unsigned suiteid,unsigned caseorder,unsigned flightid,unsigned long long randseed=1);
int runcase(PGConnection&,unsigned,int,unsigned,unsigned IDflight,unsigned IDsg,unsigned long long Np=0,unsigned Nthread=8,unsigned long long randseed=1);

boost::shared_ptr<PGConnection> dbconn;

int main(int argc,char **argv)
{
    signal(SIGHUP,SIG_IGN);
    try {
        dbconn = PGConnect();
    }
    catch(PGConnection::PGConnectionException e)
    {
        cerr << e.msg << endl;
        return -1;
    }

    unsigned IDs_request=2;

    if (argc >= 2)
        IDs_request=atoi(argv[1]);

    PGConnection::ResultType res = dbconn->execParams("SELECT suiteid,name,comment FROM suites WHERE suiteid=$1;",boost::tuples::make_tuple(IDs_request));

    unsigned IDs;
    string name,comment;

    unpackSinglePGRow(res,boost::tuples::tie(IDs,name,comment));

    cout << "Contents of suite \"" << name << "\" (ID " << IDs << ")" << endl;

    unsigned IDflight;

    try {
        IDflight = db_startFlight(dbconn.get());
    }
    catch(PGConnection::PGConnectionException e)
    {
        cerr << "Database failure: " << e.msg << endl;
    }


    // Get cases within the suite

    res = dbconn->execParams("SELECT suites_map.caseorder,cases.caseid,cases.meshid,cases.sourcegroupid,cases.materialsetid,"\
        "casename,packets,meshes.name,threads,seed FROM suites_map "\
        "JOIN cases ON cases.caseid=suites_map.caseid "\
        "JOIN materialsets ON materialsets.materialsetid=cases.materialsetid "\
        "JOIN sourcegroups ON sourcegroups.sourcegroupid=cases.sourcegroupid "\
        "JOIN meshes ON meshes.meshid=cases.meshid "\
        "WHERE suites_map.suiteid=$1 ORDER BY caseorder;",boost::tuples::make_tuple(IDs_request));


    int Nr = PQntuples(res.get());
    unsigned caseorder,IDmesh,IDsg,IDmatset,IDc,Nthread;
    unsigned long long randseed;
    unsigned long long Np;
    string casename,meshfn,meshname,sourcefn,matfn,matname;

    if(Nr <= 0)
    {
        cerr << "Incorrect number of rows found: expecting one or more" << endl;
        return -1;
    }

    for(int i=0;i<Nr;++i)
    {
        unpackPGRow(res,boost::tuples::tie(
            caseorder,
            IDc,
            IDmesh,
            IDsg,
            IDmatset,
            casename,
            Np,
            meshname,
            Nthread,
            randseed
            ), i);
        cout << "Case order (" << caseorder << ") with " << bigIntSuffix(Np) << " packets" << endl;
        cout << "  Mesh \"" << meshname << "\" (" << IDmesh << ") from " <<  meshfn << endl;
        cout << "  Source set (" << IDsg << ") from " << sourcefn << endl;
        cout << "  Material set \"" << matname << "\" (" << IDmatset << ") from " << matfn << endl;

        try {
            runcase(*dbconn,IDs_request,caseorder,IDc,IDflight,IDsg,Np,Nthread,randseed);
        }
        catch (string s)
        {
            cerr << "Caught exception string: " << s << endl;
        }
        catch (...){
            cerr << "Caught an unknown exception!" << endl;
        }
    }
}

int runcase(PGConnection& dbconn, unsigned IDs,int caseorder,unsigned IDc,unsigned IDflight,unsigned IDsg,unsigned long long Np,unsigned Nthread,unsigned long long randseed)
{
    vector<Material> materials;
    vector<Source*> sources;
//    int Nthread=8;
    char tmpdirname[] = "timos_XXXXXX";
    TetraMesh *m=NULL;
    try {
        cout << "Loading mesh - " << flush;
        m = exportMesh(dbconn,IDc);
        cout << "Done" << endl << "Loading sources - " << flush;
        exportSources(dbconn,IDsg,sources);
        cout << "Done" << endl << "Loading materials - " << flush;
        exportMaterials(dbconn,IDc,materials);
        cout << "Done" << endl;
    }
    catch(PGConnection::PGConnectionException e)
    {
        cerr << e.msg << endl;
    }

    if(!mkdtemp(tmpdirname))
    {
        cerr << "Failed to create temp dir" << endl;
        return -1;
    }

    string dirname(tmpdirname);

    string meshfn = dirname+"/input.mesh";
    string sourcefn = dirname+"/input.source";
    string matfn = dirname+"/input.opt";
    string outfn = dirname+"/result.dat";
    string logfn = dirname+"/log.out";
    string errfn = dirname+"/err.out";

    writeTIMOSMaterials(matfn,materials);
    writeTIMOSSource(sourcefn,sources,Np);
    if(!m->writeFileMatlabTP(meshfn))
    {
        cerr << "Failed to write mesh file" << endl;
        return -1;
    }

    cout << "Launching TIM-OS into folder " << dirname << endl;

    int runid = runtimos(dirname,meshfn,sourcefn,matfn,outfn,logfn,errfn,Nthread,IDs,caseorder,IDflight,randseed);

    SurfaceFluenceMap fm_surf(m);
    VolumeFluenceMap  fm_vol(m);

    // load fluence map
//    fm.readTIMOS(outfn,m);
    readTIMOSOutput(outfn,*m,fm_surf,fm_vol);

    // write fluence results to database
    db_writeResult(&dbconn,runid,fm_surf);
    db_writeResult(&dbconn,runid,fm_vol);

    return 0;
}

// spawns a TIM-OS process and pipes the output
int runtimos(string path,string meshfn,string sourcefn,string matfn,string outputfn,string logfn,string errfn,int Nthread,unsigned suiteid,unsigned caseorder,unsigned flightid,unsigned long long randseed)
{
    char Nthread_str[5],randseed_str[10];
    pid_t pid;
    sprintf(Nthread_str,"%d",Nthread);
    sprintf(randseed_str,"%lld",randseed);

    const char* const cmdline[] = {
        timos_bin,
        "-p", matfn.c_str(),
        "-f", meshfn.c_str(),
        "-s", sourcefn.c_str(),
        "-t", Nthread_str,
        "-o", outputfn.c_str(),
        "-m","si",
        "-r", randseed_str,
        NULL
    };

    // set up pipes

    posix_spawn_file_actions_t fd_map;
    posix_spawn_file_actions_init(&fd_map);

    int logfd = creat(logfn.c_str(), S_IRWXU);
    int errfd = creat(errfn.c_str(), S_IRWXU);

    if (logfd<0 || errfd<0)
        throw string("Failed to open log files");

    posix_spawn_file_actions_adddup2(&fd_map,logfd,1);   // std output
    posix_spawn_file_actions_adddup2(&fd_map,errfd,2);   // std error

    boost::timer::cpu_timer timer;

    string argstr(cmdline[0]);
    for(unsigned i=1;cmdline[i] != NULL;++i)
    {
        argstr += " ";
        argstr += cmdline[i];
    }

    timer.start();
    if(posix_spawn(&pid,timos_bin,&fd_map,NULL,(char*const*)cmdline,(char*const*)environ)) throw "Failed to spawn() process";

    close(logfd);
    close(errfd);

    unsigned runid = db_startRun(dbconn.get(),path,argstr,suiteid,caseorder,pid,flightid);
 
    cout << "Waiting for program termination" << endl;

    int exitstatus;

    waitpid(pid,&exitstatus,0);
    timer.stop();

    // get exit code
    unsigned exitcode=WEXITSTATUS(exitstatus);

    // read times
    boost::timer::cpu_times t=timer.elapsed();

    // read log files and parse out number of intersections, number of steps
    Blob logText(logfn,true);
    Blob errText(errfn,true);

    const char intsearch[]  = "Num of Intersection: ";
    const char stepsearch[] = "Num of Step: ";

    unsigned long long Nintersections=0,Nsteps=0;

    const uint8_t* p = search(errText.getPtr(),errText.getEndPtr(),intsearch,intsearch+sizeof(intsearch)-1);
    if (p != errText.getEndPtr())
        sscanf((const char*)(p+sizeof(intsearch)),"%Lu",&Nintersections);

    p = search(p,errText.getEndPtr(),stepsearch,stepsearch+sizeof(stepsearch)-1);
    if (p != errText.getEndPtr())
        sscanf((const char*)(p+sizeof(stepsearch)),"%Lu",&Nsteps);

	RunResults results;

    results.t_wall    = (double)t.wall/1e9;
    results.t_user    = (double)t.user/1e9;
    results.t_system  = (double)t.system/1e9;
	results.Nhop = Nsteps;
	results.Nintersection = Nintersections;

	results.exitcode=exitcode;
	results.log_stdout = logText;
	results.log_stderr = errText;

	db_finishRun(dbconn.get(),runid,results);

    // print to stdout
    cout << "INFO: Total execution time " << setprecision(2) << fixed << setw(8) << results.t_wall << " seconds" << endl;
    cout << "INFO: Program terminated with status code " << exitcode << endl;

    return runid;
}

