/////////////////////////////////////////////////////////////
//
//		CUDA-based Monte Carlo simulation of photon migration in semi infinite media.
//	
//			This is the version of the code used in the letter submitted to JBO-letters 2008
//			Currently the code is in an experimental state, i.e. the code is not always pretty 
//			or efficient and some ways of implementing certain aspects of the code are far 
//			from desirable. Still it should provide a good starting point for anyone interested 
//			in CUDA-based Monte Carlo simulations of photon migration. 
//
//			For the JBO-letters article the code was run on a NVIDIA 8800GT and the number of 
//			threads are hence optimized for this particular card.
//
//			We apologize for the lack of comment in the current code. We will soon re-relese 
//			this code with detailed explanations of the implementation as well as proper commenting.
//
//			To compile and run this code, please visit www.nvidia.com and download the necessary 
//			CUDA Toolkit and SKD. I also highly recommend the Visual Studio wizard 
//			(available at:http://forums.nvidia.com/index.php?showtopic=69183) 
//			if you use Visual Studio 2005 
//			(The express edition is available for free at: http://www.microsoft.com/express/2005/). 
//
//			This code is distributed under the terms of the GNU General Public Licence (see
//			below). If you use this code for academic purposes, we would greatly appreciate a 
//			citation of our letter describing GPU-based Monte Carlo simulations of photon migration. 
//
//
///////////////////////////////////////////////////////////////

/*	This file is part of CUDAMC.

    CUDAMC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CUDAMC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CUDAMC.  If not, see <http://www.gnu.org/licenses/>.*/

#include <CL/cl.h>
#include <cmath>
//#define NUM_THREADS_PER_BLOCK 320 //Keep above 192 to eliminate global memory access overhead
//#define NUM_THREADS_PER_BLOCK 448 //Keep above 192 to eliminate global memory access overhead
#define NUM_THREADS_PER_BLOCK 560 //Keep above 192 to eliminate global memory access overhead
//#define NUM_THREADS_PER_BLOCK 128 //Keep above 192 to eliminate global memory access overhead
//#define NUM_BLOCKS 84 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
//#define NUM_BLOCKS 30 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
#define NUM_BLOCKS 48 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
#define NUM_THREADS 26880
//#define NUM_THREADS 1024
#define NUM_THREADS_TIMES_THREE 80640
//#define NUM_THREADS_TIMES_THREE 3072
//#define NUM_THREADS 80640
//#define NUM_THREADS NUM_THREADS_PER_BLOCK
#define NUMSTEPS_GPU 500000
#define NUMSTEPS_CPU 500000
//#define NUMSTEPS_GPU 100
//#define NUMSTEPS_CPU 100

#define MUS_MAX 90.0f	//[1/cm]
#define V 0.0214f		//[cm/ps] (c=0.03 [cm/ps] v=c/n) here n=1.4
#define COS_CRIT 0.6999f	//the critical angle for total internal reflection at the border cos_crit=sqrt(1-(nt/ni)^2)
#define G 0.9f	
#define N 1.4f
#define PI 3.14159265f

#define TMAX 2000.0f //[ps] Maximum time of flight
#define DT 10.0f //[ps] Time binning resolution
#define TEMP 201 //ceil(TMAX/DT), precalculated to avoid dynamic memory allocation (fulhack)

#define DIRX -0.7f
#define DIRY -0.7f
#define DIRZ -0.7f

#define LINUX
//#define JORDAN
//#define NAIF
unsigned int xtest[NUM_THREADS];
unsigned int ctest[NUM_THREADS];
unsigned int atest[NUM_THREADS];

float posxtest[NUM_THREADS];
float dirxtest[NUM_THREADS];

float posytest[NUM_THREADS];
float dirytest[NUM_THREADS];

float posztest[NUM_THREADS];
float dirztest[NUM_THREADS];
struct float3{
	float x;
	float y;
	float z;

};

// forward declaration of the device code


#include <math.h>
#include <limits.h>
// forward declaration of the host code
void MCh(unsigned int*,unsigned int*,unsigned int*,unsigned int*,unsigned int*);
float rand_MWC_och(unsigned long long*,unsigned int*);
float rand_MWC_coh(unsigned long long*,unsigned int*);
void LaunchPhotonh(float3*, float3*, float*);
void Spinh(float3*,float*,unsigned long long*,unsigned int*);
unsigned int Reflecth(float3*, float3*, float*, float*, float*, float*, unsigned long long*,unsigned int*,unsigned int*);

//#include "math_constants.h"

#include <stdio.h>
//#include <cutil.h>
#include "CLMonte_goldstandard.c"
#include <time.h>
#include <stdlib.h>
#define MAX_SOURCE_SIZE (0x100000)
// wrapper for device code
int MC(unsigned int* x,unsigned int* c,unsigned int* a)
{
	unsigned int num[NUM_THREADS];
	unsigned long long num_tot, hist_tot;
	unsigned int i;
	unsigned int hist[TEMP];
	float posx[NUM_THREADS];
        float dirx[NUM_THREADS];
	float posy[NUM_THREADS];
        float diry[NUM_THREADS];
	float posz[NUM_THREADS];
        float dirz[NUM_THREADS];

	unsigned int numh[NUM_THREADS];
	unsigned int histh[TEMP];

	FILE *fp;
    char *source_str;
    size_t source_size;



#ifdef JORDAN
    fp = fopen("C:\\Users\\Jordan\\Documents\\GitHub\\ece496\\CLMonteTransport.cl", "r");
#endif
#ifdef LINUX
    fp = fopen("SRC/CLMonte_SRC/UnitTests/CLMonteBenchmark_reflect.cl", "r");

#endif
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source_str = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
    fclose( fp );


	FILE *file; 
	FILE *print_file;
	FILE *build_file;
	
	file = fopen("outp.txt", "w");
	print_file = fopen("print_out.txt", "w");
	build_file = fopen("build.txt", "w");

	clock_t time1,time2,GPUtime,CPUtime;
    int size;
    size=NUM_THREADS*sizeof(unsigned int);
	
    
	// Get platform and device information
    cl_platform_id platform_id = NULL;
    cl_device_id device_id = NULL;   
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_GPU, 1, 
            &device_id, &ret_num_devices);


	// Create an OpenCL context
    cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);

    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

    cl_mem posx_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, 
            NUM_THREADS*sizeof(float), NULL, &ret);
	if(ret!=CL_SUCCESS){
		printf("create pos buff fail xd\n");
		return -1;
	}
    
   cl_mem dirx_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, 
            NUM_THREADS*sizeof(float), NULL, &ret);
	if(ret!=CL_SUCCESS){
		printf("create dir buff fail xd\n");
		return -1;
	}
    cl_mem posy_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, 
            NUM_THREADS*sizeof(float), NULL, &ret);
	if(ret!=CL_SUCCESS){
		printf("create pos buff fail xd\n");
		return -1;
	}
    
   cl_mem diry_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, 
            NUM_THREADS*sizeof(float), NULL, &ret);
	if(ret!=CL_SUCCESS){
		printf("create dir buff fail xd\n");
		return -1;
	}
    
   cl_mem posz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, 
            NUM_THREADS*sizeof(float), NULL, &ret);
	if(ret!=CL_SUCCESS){
		printf("create pos buff fail xd\n");
		return -1;
	}
    
   cl_mem dirz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, 
            NUM_THREADS*sizeof(float), NULL, &ret);
	if(ret!=CL_SUCCESS){
		printf("create dir buff fail xd\n");
		return -1;
	}
	// Create memory buffers on the device for each vector 
    cl_mem xd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, 
            size, NULL, &ret);

    
    cl_mem cd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, 
            size, NULL, &ret);

    cl_mem ad_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, 
            size, NULL, &ret);

    cl_mem histd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, 
            TEMP*sizeof(unsigned int), NULL, &ret);
	
    ret = clEnqueueWriteBuffer(command_queue, xd_mem_obj, CL_TRUE, 0,
            size, xtest, 0, NULL, NULL);

   cl_mem retd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, 
            size, NULL, &ret);
	if(ret!=CL_SUCCESS){
		printf("Write Buff fail xd\n");
		return -1;
	}

	ret = clEnqueueWriteBuffer(command_queue, posx_mem_obj, CL_TRUE, 0,
            NUM_THREADS*sizeof(float), posxtest, 0, NULL, NULL);

	if(ret!=CL_SUCCESS){
		printf("Write Buff fail ad\n");
		return -1;
	}
	
	ret = clEnqueueWriteBuffer(command_queue, posy_mem_obj, CL_TRUE, 0,
            NUM_THREADS*sizeof(float), posytest, 0, NULL, NULL);

	if(ret!=CL_SUCCESS){
		printf("Write Buff fail ad\n");
		return -1;
	}
	ret = clEnqueueWriteBuffer(command_queue, posz_mem_obj, CL_TRUE, 0,
            NUM_THREADS*sizeof(float), posztest, 0, NULL, NULL);

	if(ret!=CL_SUCCESS){
		printf("Write Buff fail ad\n");
		return -1;
	}
	ret = clEnqueueWriteBuffer(command_queue, dirx_mem_obj, CL_TRUE, 0,
            NUM_THREADS*sizeof(float), dirxtest, 0, NULL, NULL);

	if(ret!=CL_SUCCESS){
		printf("Write Buff fail ad\n");
		return -1;
	}
	ret = clEnqueueWriteBuffer(command_queue, diry_mem_obj, CL_TRUE, 0,
            NUM_THREADS*sizeof(float), dirytest, 0, NULL, NULL);

	if(ret!=CL_SUCCESS){
		printf("Write Buff fail ad\n");
		return -1;
	}
	ret = clEnqueueWriteBuffer(command_queue, dirz_mem_obj, CL_TRUE, 0,
            NUM_THREADS*sizeof(float), dirztest, 0, NULL, NULL);

	if(ret!=CL_SUCCESS){
		printf("Write Buff fail ad\n");
		return -1;
	}

	ret = clEnqueueWriteBuffer(command_queue, ad_mem_obj, CL_TRUE, 0,
            size, atest, 0, NULL, NULL);

	if(ret!=CL_SUCCESS){
		printf("Write Buff fail ad\n");
		return -1;
	}


	ret = clEnqueueWriteBuffer(command_queue, cd_mem_obj, CL_TRUE, 0,
            size, ctest, 0, NULL, NULL);

	if(ret!=CL_SUCCESS){
		printf("Write Buff fail cd\n");
		return -1;
	}

	for(i=0;i<TEMP;i++)hist[i]=0;
	
	ret = clEnqueueWriteBuffer(command_queue, histd_mem_obj, CL_TRUE, 0,
            TEMP*sizeof(unsigned int), hist, 0, NULL, NULL);

	if(ret!=CL_SUCCESS){
		printf("Write Buff fail histd\n");
		return -1;
	}
    //load initial values
    
    

	//cudaThreadSynchronize(); //probably not necessary

	// Create a program from the kernel source
    cl_program program = clCreateProgramWithSource(context, 1, 
            (const char **)&source_str, (const size_t *)&source_size, &ret);

    // Build the program
//    ret = clBuildProgram(program, 1, &device_id, "-g -s \"C:\\Users\\Naif\\Documents\\Visual Studio 2010\\Projects\\CudaMC\\CudaMC\\CUDAMC\\CUDAMCtransport.cl\"", NULL, NULL);

	ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
	//ret = clBuildProgram(program, 1, &device_id, "-cl-unsafe-math-optimizations -cl-finite-math-only", NULL, NULL);


	if(ret!=CL_SUCCESS){
		char *build_log;
		size_t ret_val_size;
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size);
		build_log = new char[ret_val_size+1];
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL);
		build_log[ret_val_size] = '\0';
		fputs (build_log, build_file);
		delete build_log;
		fclose(build_file);
		printf("Build Program fail ad\n");
		return -1;
	}

    // Create the OpenCL kernel
    cl_kernel kernel = clCreateKernel(program, "Reflect_test1", &ret);

    // Set the arguments of the kernel
    
ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&dirx_mem_obj);
    
	if(ret!=CL_SUCCESS){
		printf("set arg 0 fail ad\n");
		return -1;
	}

ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&diry_mem_obj);
    
	if(ret!=CL_SUCCESS){
		printf("set arg 1 fail ad\n");
		return -1;
	}

ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&dirz_mem_obj);
    
	if(ret!=CL_SUCCESS){
		printf("set arg 2 fail ad\n");
		return -1;
	}

ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&posx_mem_obj);
    
	if(ret!=CL_SUCCESS){
		printf("set arg 3 fail ad\n");
		return -1;
	}

ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&posy_mem_obj);
    
	if(ret!=CL_SUCCESS){
		printf("set arg 4 fail ad\n");
		return -1;
	}
ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&posz_mem_obj);
    
	if(ret!=CL_SUCCESS){
		printf("set arg 5 fail ad\n");
		return -1;
	}
ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&xd_mem_obj);
    
	if(ret!=CL_SUCCESS){
		printf("set arg 6 fail ad\n");
		return -1;
	}

	ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&cd_mem_obj);
    
	if(ret!=CL_SUCCESS){
		printf("set arg 7 fail ad\n");
		return -1;
	}

	ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&ad_mem_obj);
    
	
	if(ret!=CL_SUCCESS){
		printf("set arg 2 fail ad\n");
		return -1;
	}

	ret = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&retd_mem_obj);
	
	if(ret!=CL_SUCCESS){
		printf("set arg 3 fail ad\n");
		return -1;
	}

	// Execute the OpenCL kernel on the list
	cl_event kern_event;
    size_t global_item_size = NUM_THREADS; // Process the entire lists
    size_t local_item_size = NUM_THREADS_PER_BLOCK; // Process in groups of NUM_THREADS_PER_BLOCk
    printf("Running code on Accelerator in OpenCL\n");
	fprintf(print_file, "Running code on Accelerator in OpenCL\n");
	ret = clFlush(command_queue);
    ret = clFinish(command_queue);
	time1=clock();
	ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, 
            &global_item_size, &local_item_size, 0, NULL, &kern_event);
	
	if(ret!=CL_SUCCESS){
		printf("enquendrange  fail ad\n");
		return -1;
	}
	clWaitForEvents(1,&kern_event);
	
	time2=clock();


	ret = clEnqueueReadBuffer(command_queue, posx_mem_obj, CL_TRUE, 0, 
            NUM_THREADS*sizeof(float), posx, 1, &kern_event, NULL);
	if(ret!=CL_SUCCESS){
		printf("read buff pos fail ad\n");
		return -1;
	}
	
	ret = clEnqueueReadBuffer(command_queue, dirx_mem_obj, CL_TRUE, 0, 
            NUM_THREADS*sizeof(float), dirx, 1, &kern_event, NULL);
	if(ret!=CL_SUCCESS){
		printf("read buff dir fail ad\n");
		return -1;
	}
	
	ret = clEnqueueReadBuffer(command_queue, posy_mem_obj, CL_TRUE, 0, 
            NUM_THREADS*sizeof(float), posy, 1, &kern_event, NULL);
	if(ret!=CL_SUCCESS){
		printf("read buff pos fail ad\n");
		return -1;
	}
	
	ret = clEnqueueReadBuffer(command_queue, diry_mem_obj, CL_TRUE, 0, 
            NUM_THREADS*sizeof(float), diry, 1, &kern_event, NULL);
	if(ret!=CL_SUCCESS){
		printf("read buff dir fail ad\n");
		return -1;
	}
	
	ret = clEnqueueReadBuffer(command_queue, posz_mem_obj, CL_TRUE, 0, 
            NUM_THREADS*sizeof(float), posz, 1, &kern_event, NULL);
	if(ret!=CL_SUCCESS){
		printf("read buff pos fail ad\n");
		return -1;
	}
	
	ret = clEnqueueReadBuffer(command_queue, dirz_mem_obj, CL_TRUE, 0, 
            NUM_THREADS*sizeof(float), dirz, 1, &kern_event, NULL);
	if(ret!=CL_SUCCESS){
		printf("read buff dir fail ad\n");
		return -1;
	}
	
	ret = clFlush(command_queue);
    ret = clFinish(command_queue);
	//cudaThreadSynchronize(); //probably not necessary
	
	
	
	hist_tot=0;
	double tot_reflected = 0.0f;
	double measured_r;
        int j;
	for(i=0; i < NUM_THREADS; i++){
		//tot_posz+=dirz[i];
		if(dirz[i] == -DIRZ)
		  tot_reflected+=1.0f; 
		//printf("dir %f, pos %f\n", dirz[i], posz[i]);

	}

//	double n1 = 1.0f;
//	double n2 = 1.4f;	
//	double avg_posz = tot_posz/NUM_THREADS;
//	double old_posz=DIRZ;
//	double expected_posz = cos(asin(n1*sin(acos(old_posz))/n2));

	float sinangle_i = sqrtf(1.0f-DIRZ*DIRZ);
	float sinangle_t = N*sinangle_i;
	float cosangle_t = sqrtf(1.0f-sinangle_t*sinangle_t);
	
	
	float cossumangle = (-DIRZ*cosangle_t) - sinangle_i*sinangle_t;
	float cosdiffangle = (-DIRZ*cosangle_t) + sinangle_i*sinangle_t;
	float sinsumangle = sinangle_i*cosangle_t + (-DIRZ*sinangle_t);
	float sindiffangle = sinangle_i*cosangle_t - (-DIRZ*sinangle_t); 

	double r = 0.5*sindiffangle*sindiffangle* ((cosdiffangle*cosdiffangle+cossumangle*cossumangle)/(sinsumangle*sinsumangle*cosdiffangle*cosdiffangle));

	measured_r = tot_reflected/NUM_THREADS;


	printf("r %f,measured_r %f\n", r, measured_r); 	
	ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(xd_mem_obj);
	ret = clReleaseMemObject(cd_mem_obj);
	ret = clReleaseMemObject(ad_mem_obj);
	ret = clReleaseMemObject(histd_mem_obj);
	ret = clReleaseCommandQueue(command_queue);

}




void initialize(void)//Straight from Steven Gratton's code
{
    FILE *fp;
    unsigned int begin=0u;
    unsigned long long int xinit=1ull;
    unsigned int cinit=0u;
    unsigned int fora,tmp1,tmp2;
#ifdef JORDAN
    fp=fopen("C:\\Users\\Jordan\\Documents\\GitHub\\ece496\\safeprimes_base32.txt","r");//use an expanded list containing 50000 safeprimes instead of Steven's shorter list
#endif
#ifdef LINUX
	fp=fopen("SRC/HelperFiles/safeprimes_base32.txt","r");//use an expanded list containing 50000 safeprimes instead of Steven's shorter list
#endif

// use begin as a multiplier to generate the initial x's for 
// the other generators...
	if(fp!=NULL)
		fscanf(fp,"%u %u %u",&begin,&tmp1,&tmp2);
	else{
		printf("bro null\n");
		exit(1);
	}


int j;

    for (unsigned int i=0;i<NUM_THREADS;i++)
    {

	posxtest[i]=0.0f;
	dirxtest[i]=DIRX;
	posytest[i]=0.0f;
	dirytest[i]=DIRY;
	posztest[i]=0.0f;
	dirztest[i]=DIRZ;

	xinit=xinit*begin+cinit;
	cinit=xinit>>32;
	xinit=xinit&0xffffffffull;
	xtest[i]=(unsigned int) xinit;
	fscanf(fp,"%u %u %u",&fora,&tmp1,&tmp2);
	atest[i]=fora;

	xinit=xinit*begin+cinit;
	cinit=xinit>>32;
	xinit=xinit&0xffffffffull;
	ctest[i]=(unsigned int) ((((double)xinit)/UINT_MAX)*fora);

    }
    fclose(fp);
}



int main(int argc,char* argv[])
{
	//do all the initialization for the RNG's (one MWC per thread)
    initialize();
    MC(xtest,ctest,atest);
	return 0;
    
}
