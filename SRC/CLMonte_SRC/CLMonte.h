//#define NUM_THREADS_PER_BLOCK 320 //Keep above 192 to eliminate global memory access overhead
//#define NUM_THREADS_PER_BLOCK 448 //Keep above 192 to eliminate global memory access overhead
//#define NUM_BLOCKS 84 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
//#define NUM_BLOCKS 30 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
//#define NUM_THREADS 80640
//#define NUM_THREADS NUM_THREADS_PER_BLOCK
//#define NUMSTEPS_GPU 100
//#define NUMSTEPS_CPU 100
#include "defines.h"
struct float3{
	float x;
	float y;
	float z;

};
// forward declaration of the host code
void MCh(unsigned int*,unsigned int*,unsigned int*,unsigned int*,unsigned int*, clock_t GPU_time);
float rand_MWC_och(unsigned long long*,unsigned int*);
float rand_MWC_coh(unsigned long long*,unsigned int*);
void LaunchPhotonh(float3*, float3*, float*);
void Spinh(float3*,float*,unsigned long long*,unsigned int*);
unsigned int Reflecth(float3*, float3*, float*, float*, float*, float*, unsigned long long*,unsigned int*,unsigned int*);
