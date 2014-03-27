
//#define INPUT_PARAMETERS

#define MAX_SOURCE_SIZE (0x100000)
#define NUM_THREADS_PER_BLOCK 560 //Keep above 192 to eliminate global memory access overhead
#define NUM_BLOCKS 48 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
#define NUM_THREADS 26880

#define NUMSTEPS_GPU 500000
//#define NUMSTEPS_GPU 270
//#define NUMSTEPS_GPU 25000
#define NUMSTEPS_CPU 500000

#define PI 3.14159265f

#define TEMP 201 //ceil(TMAX/DT), precalculated to avoid dynamic memory allocation (fulhack)

#ifndef INPUT_PARAMETERS
    #define G 0.9f	
    #define MUS_MAX 90.0f	//[1/cm]
    #define V 0.0214f		//[cm/ps] (c=0.03 [cm/ps] v=c/n) here n=1.4
    #define COS_CRIT 0.6999f	//the critical angle for total internal reflection at the border cos_crit=sqrt(1-(nt/ni)^2)
    #define N 1.4f
#endif

#define EVENT_LOGGING
//#define SPATIAL_HISTOGRAM  //SPATIAL HISTOGEAM IS DR*DT 2D histogram
//#define RUN_HOST
//#define RETURN_ON_GPU_TIME 
#define LINUX

#define DR 0.00125f
#define MAXR 0.25f

#ifdef SPATIAL_HISTOGRAM
    #define TMAX 125.0f //[ps] Maximum time of flight
    #define DT 0.625f //[ps] Time binning resolution
#else
    #define TMAX 2000.0f //[ps] Maximum time of flight
    #define DT 10.0f //[ps] Time binning resolution
#endif

#define RBUCKETS 200
#define RBUCKETSXTEMP 40200
