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
#include <cstring>


//#define NUM_THREADS_PER_BLOCK 320 //Keep above 192 to eliminate global memory access overhead
//#define NUM_THREADS_PER_BLOCK 448 //Keep above 192 to eliminate global memory access overhead
#define NUM_THREADS_PER_BLOCK 100 //Keep above 192 to eliminate global memory access overhead
//#define NUM_THREADS_PER_BLOCK 20 //Keep above 192 to eliminate global memory access overhead
//#define NUM_THREADS_PER_BLOCK 128 //Keep above 192 to eliminate global memory access overhead
//#define NUM_BLOCKS 84 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
//#define NUM_BLOCKS 30 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
#define NUM_BLOCKS 48 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
#define NUM_THREADS 100
//#define NUM_THREADS 100
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
//#define G 0.9f	
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

typedef struct{
	float x;
	float y;
	float z;

}float3;

unsigned int xtest[NUM_THREADS];
unsigned int ctest[NUM_THREADS];
unsigned int atest[NUM_THREADS];

float posxtest[NUM_THREADS];
float dirxtest[NUM_THREADS];

float posytest[NUM_THREADS];
float dirytest[NUM_THREADS];

float posztest[NUM_THREADS];
float dirztest[NUM_THREADS];

float axtest[NUM_THREADS];
float aytest[NUM_THREADS];
float aztest[NUM_THREADS];

float bxtest[NUM_THREADS];
float bytest[NUM_THREADS];
float bztest[NUM_THREADS];

#ifdef LINUX
#include "CLMonte_host_func.cpp"
#endif
#ifdef NAIF
#include "C:\\Users\\Naif\\Documents\\ece496\\SRC\\CLMonte_SRC\\UnitTests\\CLMonte_host_func.cpp"
#endif
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
//#include "CLMonte_goldstandard.c"
#include <time.h>
#include <stdlib.h>
#define MAX_SOURCE_SIZE (0x100000)
// wrapper for device code

// takes 2 floats x1,y1,z1 and x2,y2,z2 returns result of cross product in float3 struct
void cross_product (float x1, float y1, float z1, float x2, float y2, float z2, float* ret_x, float* ret_y, float* ret_z) {
    *ret_x = y1*z2 - z1*y2;
    *ret_y = z1*x2 - x1*z2;
    *ret_z = x1*y2 - y1*x2;
}

float dot_product (float x1, float y1, float z1, float x2, float y2, float z2) {
    float ret = 0;
    ret += x1*x2;
    ret += y1*y2;
    ret += z1*z2;
    return ret;
}    

void initialize_vectors(float posx, float posy, float posz, float dirx, float diry, float dirz, float ax, float ay, float az, float bx, float by, float bz){
	int i;
    for(i=0; i<NUM_THREADS; i++){
        posxtest[i]=posx;
	    dirxtest[i]=dirx;
	    posytest[i]=posy;
	    dirytest[i]=diry;
	    posztest[i]=posz;
	    dirztest[i]=dirz;
        axtest[i]=ax;
        aytest[i]=ay;
        aztest[i]=az;
        bxtest[i]=bx;
        bytest[i]=by;
        bztest[i]=bz;
    }
}

void initialize_vectors_reflect () {
    int i;
    float temp;
    for (i=0; i<NUM_THREADS; i++) {
        dirxtest[i] = ((float) rand() / (float)(RAND_MAX));
        dirytest[i] = ((float) rand() / (float)(RAND_MAX));
        dirztest[i] = ((float) rand() / (float)(RAND_MAX));
        temp = rsqrtf(dirxtest[i]*dirxtest[i] + dirytest[i]*dirytest[i] + dirztest[i]*dirztest[i]);
        dirxtest[i] = dirxtest[i]*temp;
        dirytest[i] = dirytest[i]*temp;
        dirztest[i] = dirztest[i]*temp;
        axtest[i]=0.0f;
        aytest[i]=0.0f;
        aztest[i]=0.0f;
        bxtest[i]=0.0f;
        bytest[i]=0.0f;
        bztest[i]=0.0f;
    }
}

int compare_to_expected (FILE* fp, const char* var, float3 calculated, float3 expected) {
    int ret = 0;
    if (fabs(calculated.x - expected.x) >= 0.000001f){
        printf("%s x-component not as expected\n\texpected %f, got %f\n",var,expected.x,calculated.x);
        fprintf(fp,"%s x-component not as expected\n\texpected %f, got %f\n",var,expected.x,calculated.x);
        ret = -1;
    }        
    if (fabs(calculated.y - expected.y) >= 0.000001f){
        printf("%s y-component not as expected\n\texpected %f, got %f\n",var,expected.y,calculated.y);
        fprintf(fp,"%s y-component not as expected\n\texpected %f, got %f\n",var,expected.y,calculated.y);
        ret = -1;
    }        
    if (fabs(calculated.z - expected.z) >= 0.000001f){
        printf("%s z-component not as expected\n\texpected %f, got %f\n",var,expected.z,calculated.z);
        fprintf(fp,"%s z-component not as expected\n\texpected %f, got %f\n",var,expected.z,calculated.z);
        ret = -1;
    }
    return ret;
}

int compare_to_expected_float (FILE* fp, const char* var, float calculated, float expected) {
    int ret = 0;
    if (fabs(calculated - expected) >= 0.000001f){
        printf("%s not as expected\n\texpected %f, got %f\n",var,expected,calculated);
        fprintf(fp,"%s not as expected\n\texpected %f, got %f\n",var,expected,calculated);
        ret = -1;
    }        
    return ret;
}

void check_return(cl_int ret, const char* message){
    if(ret!=CL_SUCCESS){
        printf("%d %s\n",(int)ret, message);
        exit(-1);
    }
}



int reflect(float posx, float posy, float posz, float dirx, float diry, float dirz, unsigned int* x,unsigned int* c,unsigned int* a)
{
	unsigned int num[NUM_THREADS];
	unsigned long long num_tot, hist_tot;
	unsigned int i;
	unsigned int hist[TEMP];
	float posx_array[NUM_THREADS];
    float dirx_array[NUM_THREADS];
	float posy_array[NUM_THREADS];
    float diry_array[NUM_THREADS];
	float posz_array[NUM_THREADS];
    float dirz_array[NUM_THREADS];
    

	unsigned int numh[NUM_THREADS];
	unsigned int histh[TEMP];

	FILE *fp;
    char *source_str;
    size_t source_size;

    initialize_vectors(posx, posy, posz, dirx, diry, dirz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

#ifdef LINUX
    fp = fopen("SRC/CLMonte_SRC/UnitTests/CLMonteUT.cl", "r");
#endif
#ifdef NAIF
	fp = fopen("C:\\Users\\Naif\\Documents\\ece496\\SRC\\CLMonte_SRC\\UnitTests\\CLMonte_host_func.cpp", "r");
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
	
#ifdef LINUX
    file = fopen("outp.txt", "w");
    print_file = fopen("print_out.txt", "w");
    build_file = fopen("build.txt", "w");
#endif

#ifdef NAIF
	file =fopen("C:\\Users\\Naif\\Documents\\ece496\\build.txt", "w");
	print_file = fopen("C:\\Users\\Naif\\Documents\\ece496\\print_out.txt", "w");
    build_file = fopen("C:\\Users\\Naif\\Documents\\ece496\\build.txt", "w");
#endif
//    clock_t time1,time2,GPUtime,CPUtime;
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

    cl_mem posx_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
	check_return(ret, "create posx buff fail\n");
    
    cl_mem dirx_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create dirx buff fail xd\n");
   
    cl_mem posy_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
	check_return(ret,"create posy buff fail xd\n");
    
    cl_mem diry_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret,"create diry buff fail xd\n");
    
    cl_mem posz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret,"create posz buff fail xd\n");
    
    cl_mem dirz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create dirz buff fail xd\n");
	// Create memory buffers on the device for each vector 
    cl_mem xd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create xd buff fail\n");

    
    cl_mem cd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create cd buff fail\n");

    cl_mem ad_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create ad buff fail\n");

    cl_mem histd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, TEMP*sizeof(unsigned int), NULL, &ret);
    check_return(ret, "create histd buff fail\n");
	
    ret = clEnqueueWriteBuffer(command_queue, xd_mem_obj, CL_TRUE, 0,size, xtest, 0, NULL, NULL);
    check_return(ret, "write buff xd fail\n");

    cl_mem retd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
    check_return(ret,"write buff retd fail \n");

	ret = clEnqueueWriteBuffer(command_queue, posx_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), posxtest, 0, NULL, NULL);
    check_return(ret,"write buff posx fail \n");
	
	ret = clEnqueueWriteBuffer(command_queue, posy_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), posytest, 0, NULL, NULL);
    check_return(ret,"write buff posy fail \n");
	
    ret = clEnqueueWriteBuffer(command_queue, posz_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), posztest, 0, NULL, NULL);
    check_return(ret,"write buff posz fail \n");

	ret = clEnqueueWriteBuffer(command_queue, dirx_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), dirxtest, 0, NULL, NULL);
    check_return(ret,"write buff dirx fail \n");

	ret = clEnqueueWriteBuffer(command_queue, diry_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirytest, 0, NULL, NULL);
    check_return(ret,"write buff diry fail \n");

	ret = clEnqueueWriteBuffer(command_queue, dirz_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), dirztest, 0, NULL, NULL);
    check_return(ret,"write buff dirz fail \n");

	ret = clEnqueueWriteBuffer(command_queue, ad_mem_obj, CL_TRUE, 0, size, atest, 0, NULL, NULL);
    check_return(ret,"write buff ad fail \n");

	ret = clEnqueueWriteBuffer(command_queue, cd_mem_obj, CL_TRUE, 0, size, ctest, 0, NULL, NULL);
    check_return(ret,"write buff cd fail \n");


	for(i=0;i<TEMP;i++)hist[i]=0;
	
	ret = clEnqueueWriteBuffer(command_queue, histd_mem_obj, CL_TRUE, 0, TEMP*sizeof(unsigned int), hist, 0, NULL, NULL);
    check_return(ret,"write buff histd fail \n");

	// Create a program from the kernel source
    cl_program program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
#ifdef LINUX
	ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
#endif
#ifdef NAIF
    ret = clBuildProgram(program, 1, &device_id, "-g -s \"C:\\Users\\Naif\\Documents\\Visual Studio 2010\\Projects\\CudaMC\\CudaMC\\CUDAMC\\CUDAMCtransport.cl\"", NULL, NULL);
#endif
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
    check_return(ret, "create kernel Reflect_test1 fail"); 

    // Set the arguments of the kernel
    
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&dirx_mem_obj);
    check_return(ret, "set arg 0 fail\n"); 

    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&diry_mem_obj);
    check_return(ret, "set arg 1 fail\n"); 
    
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&dirz_mem_obj);
    check_return(ret, "set arg 2 fail\n"); 
    
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&posx_mem_obj);
    check_return(ret, "set arg 3 fail\n"); 

    ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&posy_mem_obj);
    check_return(ret, "set arg 4 fail\n"); 
    
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&posz_mem_obj);
    check_return(ret, "set arg 5 fail\n"); 
    
    ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&xd_mem_obj);
    check_return(ret, "set arg 6 fail\n"); 
    
	ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&cd_mem_obj);
    check_return(ret, "set arg 7 fail\n"); 
    
	ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&ad_mem_obj);
    check_return(ret, "set arg 8 fail\n"); 

	ret = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&retd_mem_obj);
    check_return(ret, "set arg 9 fail\n"); 
	
	// Execute the OpenCL kernel on the list
	cl_event kern_event;
    size_t global_item_size = NUM_THREADS; // Process the entire lists
    size_t local_item_size = NUM_THREADS_PER_BLOCK; // Process in groups of NUM_THREADS_PER_BLOCk
    printf("Running code on Accelerator in OpenCL\n");
	fprintf(print_file, "Running code on Accelerator in OpenCL\n");
	ret = clFlush(command_queue);
    ret = clFinish(command_queue);
	ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, &kern_event);
    check_return(ret, "EnqueuendRangeFail"); 
	
	clWaitForEvents(1,&kern_event);
	

	ret = clEnqueueReadBuffer(command_queue, posx_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), posx_array, 1, &kern_event, NULL);
    check_return(ret, "read buff posx fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, dirx_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirx_array, 1, &kern_event, NULL);
    check_return(ret, "read buff dirx fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, posy_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), posy_array, 1, &kern_event, NULL);
    check_return(ret, "read buff posy fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, diry_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), diry_array, 1, &kern_event, NULL);
    check_return(ret, "read buff diry fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, posz_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), posz_array, 1, &kern_event, NULL);
    check_return(ret, "read buff posz fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, dirz_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirz_array, 1, &kern_event, NULL);
    check_return(ret, "read buff dirz fail\n"); 
	
	ret = clFlush(command_queue);
    ret = clFinish(command_queue);
	
	
	
	hist_tot=0;
	double tot_reflected = 0.0f;
	double measured_r;
   int j;
	for(i=0; i < NUM_THREADS; i++){
		if(dirz_array[i] == -dirz)
		  tot_reflected+=1.0f; 
		//printf("dir %f, pos %f\n", dirz[i], posz[i]);

	}


	float sinangle_i = sqrtf(1.0f-dirz*dirz);
	float sinangle_t = N*sinangle_i;
	float cosangle_t = sqrtf(1.0f-sinangle_t*sinangle_t);
	
	
	float cossumangle = (-dirz*cosangle_t) - sinangle_i*sinangle_t;
	float cosdiffangle = (-dirz*cosangle_t) + sinangle_i*sinangle_t;
	float sinsumangle = sinangle_i*cosangle_t + (-dirz*sinangle_t);
	float sindiffangle = sinangle_i*cosangle_t - (-dirz*sinangle_t); 

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
int spinCPU(float dirx, float diry, float dirz, unsigned int* x,unsigned int* c,unsigned int* a, float G){
    float3 dir[NUM_THREADS];
    
    float cost_array[NUM_THREADS];
    float sint_array[NUM_THREADS];
    float cosp_array[NUM_THREADS];
    float sinp_array[NUM_THREADS];
    int i;
    for(i=0; i<NUM_THREADS; i++){
	    dir[i].x=dirx;
	    dir[i].y=diry;
	    dir[i].z=dirz;
    }


    Spin_testh(dir, cost_array, sint_array, cosp_array, sinp_array, G);

	float total= 0.0f;
	for(i=0; i < NUM_THREADS; i++){
          if(cost_array[i]<0)
              total-=cost_array[i];
          else
              total+=cost_array[i];
          //printf("cost_array[%d] %f\n", cost_array[i], i);
	}
    
    printf("Spin Test1 CPU - Testing average cost vs expected value G: average %f, expected %f \n", total/NUM_THREADS, G);
    
	return 0;
}


int spin(float dirx, float diry, float dirz, unsigned int* x,unsigned int* c,unsigned int* a, float G)
{
	
    unsigned int num[NUM_THREADS];
	unsigned long long num_tot, hist_tot;
	unsigned int i;
	unsigned int hist[TEMP];
    float dirx_array[NUM_THREADS];
    float diry_array[NUM_THREADS];
    float dirz_array[NUM_THREADS];
    
    float cost_array[NUM_THREADS];
    float sint_array[NUM_THREADS];
    float cosp_array[NUM_THREADS];
    float sinp_array[NUM_THREADS];
//    float3 a_array[NUM_THREADS];
//    float3 b_array[NUM_THREADS];

	unsigned int numh[NUM_THREADS];
	unsigned int histh[TEMP];

	FILE *fp;
    char *source_str;
    size_t source_size;

    initialize_vectors(0.0f, 0.0f, 0.0f, dirx, diry, dirz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    fp = fopen("C:\\Users\\Naif\\Documents\\ece496\\SRC\\CLMonte_SRC\\UnitTests\\CLMonteUT.cl", "r");

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
	

    file = fopen("C:\\Users\\Naif\\Documents\\ece496\\outp.txt", "w");
    print_file = fopen("C:\\Users\\Naif\\Documents\\ece496\\print_out.txt", "w");
    build_file = fopen("C:\\Users\\Naif\\Documents\\ece496\\build.txt", "w");

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

    
    cl_mem dirx_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create dirx buff fail xd\n");
   
    
    cl_mem diry_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret,"create diry buff fail xd\n");
    
    
    cl_mem dirz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create dirz buff fail xd\n");
	// Create memory buffers on the device for each vector 
    cl_mem xd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create xd buff fail\n");

    
    cl_mem cd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create cd buff fail\n");

    cl_mem ad_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create ad buff fail\n");

    cl_mem sintd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
    check_return(ret, "create sintd buff fail\n");
    
    cl_mem costd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
    check_return(ret, "create costd buff fail\n");
    
    cl_mem cospd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
    check_return(ret, "create cospd buff fail\n");
	
    cl_mem sinpd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
    check_return(ret, "create sinpd buff fail\n");
    
    
    ret = clEnqueueWriteBuffer(command_queue, xd_mem_obj, CL_TRUE, 0,size, xtest, 0, NULL, NULL);
    check_return(ret, "write buff xd fail\n");

	ret = clEnqueueWriteBuffer(command_queue, dirx_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), dirxtest, 0, NULL, NULL);
    check_return(ret,"write buff dirx fail \n");

	ret = clEnqueueWriteBuffer(command_queue, diry_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirytest, 0, NULL, NULL);
    check_return(ret,"write buff diry fail \n");

	ret = clEnqueueWriteBuffer(command_queue, dirz_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), dirztest, 0, NULL, NULL);
    check_return(ret,"write buff dirz fail \n");

	ret = clEnqueueWriteBuffer(command_queue, ad_mem_obj, CL_TRUE, 0, size, atest, 0, NULL, NULL);
    check_return(ret,"write buff ad fail \n");

	ret = clEnqueueWriteBuffer(command_queue, cd_mem_obj, CL_TRUE, 0, size, ctest, 0, NULL, NULL);
    check_return(ret,"write buff cd fail \n");


	for(i=0;i<TEMP;i++)hist[i]=0;
	

	// Create a program from the kernel source
    cl_program program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);

    ret = clBuildProgram(program, 1, &device_id, "-g -s \"C:\\Users\\Naif\\Documents\\ece496\\SRC\\CLMonte_SRC\\UnitTests\\CLMonteUT.cl", NULL, NULL);

	//ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
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
    cl_kernel kernel = clCreateKernel(program, "Spin_test1", &ret);
    check_return(ret, "create kernel Spin_test1 fail"); 

    // Set the arguments of the kernel
    
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&dirx_mem_obj);
    check_return(ret, "set arg 0 fail\n"); 

    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&diry_mem_obj);
    check_return(ret, "set arg 1 fail\n"); 
    
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&dirz_mem_obj);
    check_return(ret, "set arg 2 fail\n"); 
    
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&xd_mem_obj);
    check_return(ret, "set arg 3 fail\n"); 
    
	ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&ad_mem_obj);
    check_return(ret, "set arg 4 fail\n"); 
	
	ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&cd_mem_obj);
    check_return(ret, "set arg 5 fail\n");
     
	ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&costd_mem_obj);
    check_return(ret, "set arg 6 fail\n");
	
	ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&sintd_mem_obj);
    check_return(ret, "set arg 7 fail\n");
	
    ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&cospd_mem_obj);
    check_return(ret, "set arg 8 fail\n");
    
    ret = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&sinpd_mem_obj);
    check_return(ret, "set arg 9 fail\n");
    
    ret = clSetKernelArg(kernel, 10, sizeof(float), (void *)&G);
    check_return(ret, "set arg 10 fail\n");
    // Execute the OpenCL kernel on the list
	cl_event kern_event;
    size_t global_item_size = NUM_THREADS; // Process the entire lists
    size_t local_item_size = NUM_THREADS_PER_BLOCK; // Process in groups of NUM_THREADS_PER_BLOCk
    printf("Running code on Accelerator in OpenCL\n");
	fprintf(print_file, "Running code on Accelerator in OpenCL\n");
	ret = clFlush(command_queue);
    ret = clFinish(command_queue);
	ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, &kern_event);
    check_return(ret, "EnqueuendRangeFail"); 
	
	clWaitForEvents(1,&kern_event);
	
	ret = clEnqueueReadBuffer(command_queue, dirx_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirx_array, 1, &kern_event, NULL);
    check_return(ret, "read buff dirx fail\n"); 
	
	
	ret = clEnqueueReadBuffer(command_queue, diry_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), diry_array, 1, &kern_event, NULL);
    check_return(ret, "read buff diry fail\n"); 
	
	
	ret = clEnqueueReadBuffer(command_queue, dirz_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirz_array, 1, &kern_event, NULL);
    check_return(ret, "read buff dirz fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, costd_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), cost_array, 1, &kern_event, NULL);
    check_return(ret, "read buff cost fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, sintd_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), sint_array, 1, &kern_event, NULL);
    check_return(ret, "read buff sint fail\n"); 
    
	ret = clEnqueueReadBuffer(command_queue, cospd_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), cosp_array, 1, &kern_event, NULL);
    check_return(ret, "read buff cosp fail\n"); 
    
	ret = clEnqueueReadBuffer(command_queue, sinpd_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), sinp_array, 1, &kern_event, NULL);
    check_return(ret, "read buff sinp fail\n"); 
    
    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
	
	
	float total= 0.0f;
	for(i=0; i < NUM_THREADS; i++){
          total+=cost_array[i];
          //printf("cost_array[%d] %f\n", cost_array[i], i);
	}
    
    printf("Spin Test1 - Testing average cost vs expected value G: average %f, expected %f \n", total/NUM_THREADS, G);
    
    FILE *results;
    float3 a_v;
    float3 b_v;
    float3 a_normalized;
    float3 b_normalized;
    float d_dot_a;
    float d_dot_b;
    int test_case2_fail = 0;

    cross_product(dirx, diry, dirz, 0, 0, 1, &a_v.x, &a_v.y, &a_v.z);
    cross_product(dirx, diry, dirz, a_v.x, a_v.y, a_v.z, &b_v.x, &b_v.y, &b_v.z);
   
    
    
    a_normalized.x = a_v.x/(sqrt(dot_product(a_v.x, a_v.y, a_v.z, a_v.x, a_v.y, a_v.z)));
    a_normalized.y = a_v.y/(sqrt(dot_product(a_v.x, a_v.y, a_v.z, a_v.x, a_v.y, a_v.z)));
    a_normalized.z = a_v.z/(sqrt(dot_product(a_v.x, a_v.y, a_v.z, a_v.x, a_v.y, a_v.z)));

    b_normalized.x = b_v.x/(sqrt(dot_product(b_v.x, b_v.y, b_v.z, b_v.x, b_v.y, b_v.z)));
    b_normalized.y = b_v.y/(sqrt(dot_product(b_v.x, b_v.y, b_v.z, b_v.x, b_v.y, b_v.z)));
    b_normalized.z = b_v.z/(sqrt(dot_product(b_v.x, b_v.y, b_v.z, b_v.x, b_v.y, b_v.z)));
    
    printf("a dot product with itself = %f, b dot product with itself %f.\n", dot_product(a_v.x, a_v.y, a_v.z, a_v.x, a_v.y, a_v.z), dot_product(b_v.x, b_v.y, b_v.z, b_v.x, b_v.y, b_v.z));

    results=fopen("spintest.txt","w");
    //create a and b arrays
    for (i=0; i < NUM_THREADS; i++) {
        //a_array[i] = cross_product(dirx_array[i], diry_array[i], dirz_array[i], 0, 0, 1);
        //b_array[i] = cross_product(dirx_array[i], diry_array[i], dirz_array[i], a_array[i].x, a_array[i].y, a_array[i].z);
//        printf("a[%d].x = %f, a[%d].y = %f, a[%d].z = %f; b[%d].x = %f, b[%d].y = %f, b[%d].z = %f\n", i,a_array[i].x,i,a_array[i].y,i,a_array[i].z,i,b_array[i].x,i,b_array[i].y,i,b_array[i].z);  
        //d_dot_a = dot_product(dirx_array[i], diry_array[i], dirz_array[i], a_array[i].x, a_array[i].y, a_array[i].z);
        //d_dot_b = dot_product(dirx_array[i], diry_array[i], dirz_array[i], b_array[i].x, b_array[i].y, b_array[i].z);
        d_dot_a = dot_product(dirx_array[i], diry_array[i], dirz_array[i], a_v.x, a_v.y, a_v.z);
        d_dot_b = dot_product(dirx_array[i], diry_array[i], dirz_array[i], b_v.x, b_v.y, b_v.z);
        //printf("d_dot_a = %f, d_dot_b = %f\n",d_dot_a, d_dot_b);
        if ((d_dot_a*d_dot_a + d_dot_b*d_dot_b) - (sint_array[i]*sint_array[i]) >= 0.000001f){ 
            printf("Spin Test2 failed on i = %d, left side (d_dot_a)^2+(d_dot_b)^2 = %f, right side sin^2(theta) = %f diff %f\n",
            i, (d_dot_a*d_dot_a + d_dot_b*d_dot_b), (sint_array[i]*sint_array[i]), 
            (d_dot_a*d_dot_a + d_dot_b*d_dot_b) - (sint_array[i]*sint_array[i]));
            int test_case2_fail = 1;
        }     
        d_dot_a = dot_product(dirx_array[i], diry_array[i], dirz_array[i], a_normalized.x, a_normalized.y, a_normalized.z);
        d_dot_b = dot_product(dirx_array[i], diry_array[i], dirz_array[i], b_normalized.x, b_normalized.y, b_normalized.z);
        fprintf (results, "%f %f\n", d_dot_a, d_dot_b);
        
        
    }
    if (test_case2_fail) 
        printf("Spin Test2 failed\n");
    else 
        printf("Spin Test2 successful, left side (d_dot_a)^2+(d_dot_b)^2 = right side sin^2(theta) for all threads\n");    
    
    fclose(results);
     
    
	ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(xd_mem_obj);
	ret = clReleaseMemObject(cd_mem_obj);
	ret = clReleaseMemObject(ad_mem_obj);
	ret = clReleaseCommandQueue(command_queue);
	return 0;
}


int spin_new(float dirx, float diry, float dirz, float ax, float ay, float az, float bx, float by, float bz, unsigned int* x,unsigned int* c,unsigned int* a, float G)
{
	
    unsigned int num[NUM_THREADS];
	unsigned long long num_tot, hist_tot;
	unsigned int i;
	unsigned int hist[TEMP];
    float dirx_array[NUM_THREADS];
    float diry_array[NUM_THREADS];
    float dirz_array[NUM_THREADS];
    
    float ax_array[NUM_THREADS];
    float ay_array[NUM_THREADS];
    float az_array[NUM_THREADS];
    
    float bx_array[NUM_THREADS];
    float by_array[NUM_THREADS];
    float bz_array[NUM_THREADS];
    
    float cost_array[NUM_THREADS];
    float sint_array[NUM_THREADS];
    float cosp_array[NUM_THREADS];
    float sinp_array[NUM_THREADS];

	unsigned int numh[NUM_THREADS];
	unsigned int histh[TEMP];

	FILE *fp;
    char *source_str;
    size_t source_size;

    initialize_vectors(0.0f, 0.0f, 0.0f, dirx, diry, dirz, ax, ay, az, bx, by, bz);

#ifdef LINUX
    fp = fopen("SRC/CLMonte_SRC/UnitTests/CLMonteUT.cl", "r");
#endif
#ifdef NAIF
    fp = fopen("C:\\Users\\Naif\\Documents\\ece496\\SRC\\CLMonte_SRC\\UnitTests\\CLMonteUT.cl", "r");
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

#ifdef LINUX    
    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &ret_num_devices);
#endif
#ifdef NAIF    
    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_CPU, 1, &device_id, &ret_num_devices);
#endif

	// Create an OpenCL context
    cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);

    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

    
    cl_mem dirx_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create dirx buff fail xd\n");
	if(ret!=CL_SUCCESS){
		printf("create dirx buff fail \n");
		exit(-1);
	}
   
    
    cl_mem diry_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret,"create diry buff fail xd\n");
    
    
    cl_mem dirz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create dirz buff fail xd\n");
 
 	
    cl_mem ax_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create ax buff fail xd\n");
  
    
    cl_mem ay_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret,"create ay buff fail xd\n");
    
    
    cl_mem az_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create az buff fail xd\n");
   
    
    cl_mem bx_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create bx buff fail xd\n");
  
    
    cl_mem by_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret,"create by buff fail xd\n");
    
    
    cl_mem bz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create bz buff fail xd\n");
    
    // Create memory buffers on the device for each vector 
    cl_mem xd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create xd buff fail\n");

    
    cl_mem cd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create cd buff fail\n");

    cl_mem ad_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create ad buff fail\n");

    cl_mem sintd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
    check_return(ret, "create sintd buff fail\n");
    
    cl_mem costd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
    check_return(ret, "create costd buff fail\n");
    
    cl_mem cospd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
    check_return(ret, "create cospd buff fail\n");
	
    cl_mem sinpd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
    check_return(ret, "create sinpd buff fail\n");
    
    
    ret = clEnqueueWriteBuffer(command_queue, xd_mem_obj, CL_TRUE, 0,size, xtest, 0, NULL, NULL);
    check_return(ret, "write buff xd fail\n");

	ret = clEnqueueWriteBuffer(command_queue, dirx_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), dirxtest, 0, NULL, NULL);
    check_return(ret,"write buff dirx fail \n");

	ret = clEnqueueWriteBuffer(command_queue, diry_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirytest, 0, NULL, NULL);
    check_return(ret,"write buff diry fail \n");

	ret = clEnqueueWriteBuffer(command_queue, dirz_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), dirztest, 0, NULL, NULL);
    check_return(ret,"write buff dirz fail \n");
	
    ret = clEnqueueWriteBuffer(command_queue, ax_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), axtest, 0, NULL, NULL);
    check_return(ret,"write buff ax fail \n");

	ret = clEnqueueWriteBuffer(command_queue, ay_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), aytest, 0, NULL, NULL);
    check_return(ret,"write buff ay fail \n");

	ret = clEnqueueWriteBuffer(command_queue, az_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), aztest, 0, NULL, NULL);
    check_return(ret,"write buff az fail \n");

    ret = clEnqueueWriteBuffer(command_queue, bx_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), bxtest, 0, NULL, NULL);
    check_return(ret,"write buff bx fail \n");

	ret = clEnqueueWriteBuffer(command_queue, by_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), bytest, 0, NULL, NULL);
    check_return(ret,"write buff by fail \n");

	ret = clEnqueueWriteBuffer(command_queue, bz_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), bztest, 0, NULL, NULL);
    check_return(ret,"write buff bz fail \n");
	
    ret = clEnqueueWriteBuffer(command_queue, ad_mem_obj, CL_TRUE, 0, size, atest, 0, NULL, NULL);
    check_return(ret,"write buff ad fail \n");

	ret = clEnqueueWriteBuffer(command_queue, cd_mem_obj, CL_TRUE, 0, size, ctest, 0, NULL, NULL);
    check_return(ret,"write buff cd fail \n");


	for(i=0;i<TEMP;i++)hist[i]=0;
	

	// Create a program from the kernel source
    cl_program program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);

#ifdef LINUX
    ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
#endif
#ifdef NAIF
    ret = clBuildProgram(program, 1, &device_id, "-g -s \"C:\\Users\\Naif\\Documents\\ece496\\SRC\\CLMonte_SRC\\UnitTests\\CLMonteUT.cl\"", NULL, NULL);
#endif
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
    cl_kernel kernel = clCreateKernel(program, "newSpin_test1", &ret);
    check_return(ret, "create kernel newSpin_test1 fail"); 

    // Set the arguments of the kernel
    
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&dirx_mem_obj);
    check_return(ret, "set arg 0 fail\n"); 

    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&diry_mem_obj);
    check_return(ret, "set arg 1 fail\n"); 
    
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&dirz_mem_obj);
    check_return(ret, "set arg 2 fail\n"); 
    
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&ax_mem_obj);
    check_return(ret, "set arg 3 fail\n"); 

    ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&ay_mem_obj);
    check_return(ret, "set arg 4 fail\n"); 
    
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&az_mem_obj);
    check_return(ret, "set arg 5 fail\n"); 
    
    ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&bx_mem_obj);
    check_return(ret, "set arg 6 fail\n"); 

    ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&by_mem_obj);
    check_return(ret, "set arg 7 fail\n"); 
    
    ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&bz_mem_obj);
    check_return(ret, "set arg 8 fail\n"); 
    
    ret = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&xd_mem_obj);
    check_return(ret, "set arg 9 fail\n"); 
    
	ret = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&ad_mem_obj);
    check_return(ret, "set arg 10 fail\n"); 
	
	ret = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *)&cd_mem_obj);
    check_return(ret, "set arg 11 fail\n");
     
	ret = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *)&costd_mem_obj);
    check_return(ret, "set arg 12 fail\n");
	
	ret = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *)&sintd_mem_obj);
    check_return(ret, "set arg 13 fail\n");
	
    ret = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *)&cospd_mem_obj);
    check_return(ret, "set arg 14 fail\n");
    
    ret = clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *)&sinpd_mem_obj);
    check_return(ret, "set arg 15 fail\n");
    
    ret = clSetKernelArg(kernel, 16, sizeof(float), (void *)&G);
    check_return(ret, "set arg 16 fail\n");
    // Execute the OpenCL kernel on the list
	cl_event kern_event;
    size_t global_item_size = NUM_THREADS; // Process the entire lists
    size_t local_item_size = NUM_THREADS_PER_BLOCK; // Process in groups of NUM_THREADS_PER_BLOCk
    printf("Running code on Accelerator in OpenCL\n");
	fprintf(print_file, "Running code on Accelerator in OpenCL\n");
	ret = clFlush(command_queue);
    ret = clFinish(command_queue);
	ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, &kern_event);
  //  check_return(ret, "EnqueuendRangeFail\n"); 
	
	ret = clFlush(command_queue);
    if(ret!=CL_SUCCESS){
        printf("%d bro 1\n", ret);
        exit(-1);
    }
    ret = clFinish(command_queue);
//    check_return(ret, "EnqueuendRangeFail\n"); 
//    check_return(ret, "read buff dirx fail\n"); 	
	
    if(ret!=CL_SUCCESS){
        printf("%d bro 2\n", ret);
        exit(-1);
    }
//exit(0);	
	ret = clEnqueueReadBuffer(command_queue, dirx_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirx_array, 1, &kern_event, NULL);
    check_return(ret, "read buff dirx fail\n"); 	
//exit(0);	

	ret = clEnqueueReadBuffer(command_queue, diry_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), diry_array, 1, &kern_event, NULL);
    check_return(ret, "read buff diry fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, dirz_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirz_array, 1, &kern_event, NULL);
    check_return(ret, "read buff dirz fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, ax_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), ax_array, 1, &kern_event, NULL);
    check_return(ret, "read buff ax fail\n"); 	
	
	ret = clEnqueueReadBuffer(command_queue, ay_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), ay_array, 1, &kern_event, NULL);
    check_return(ret, "read buff ay fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, az_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), az_array, 1, &kern_event, NULL);
    check_return(ret, "read buff az fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, bx_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), bx_array, 1, &kern_event, NULL);
    check_return(ret, "read buff bx fail\n"); 	
	
	ret = clEnqueueReadBuffer(command_queue, by_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), by_array, 1, &kern_event, NULL);
    check_return(ret, "read buff by fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, bz_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), bz_array, 1, &kern_event, NULL);
    check_return(ret, "read buff bz fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, costd_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), cost_array, 1, &kern_event, NULL);
    check_return(ret, "read buff cost fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, sintd_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), sint_array, 1, &kern_event, NULL);
    check_return(ret, "read buff sint fail\n"); 
    
	ret = clEnqueueReadBuffer(command_queue, cospd_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), cosp_array, 1, &kern_event, NULL);
    check_return(ret, "read buff cosp fail\n"); 
    
	ret = clEnqueueReadBuffer(command_queue, sinpd_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), sinp_array, 1, &kern_event, NULL);
    check_return(ret, "read buff sinp fail\n"); 
    
    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
	
	
	float total= 0.0f;
	for(i=0; i < NUM_THREADS; i++){
          total+=cost_array[i];
          //printf("cost_array[%d] %f\n", cost_array[i], i);
	}
    
    printf("Spin Test1 - Testing average cost vs expected value G: average %f, expected %f \n", total/NUM_THREADS, G);
    
    FILE *results;
    results = fopen("spintest.txt","w");
    float3 a_exp;
    float3 b_exp;
    float3 d_exp;
    float3 d_calc;
    float3 a_calc;
    float3 b_calc;
    float temp;
    bool testcase2fail = false;
    for (i=0; i< NUM_THREADS; i++) {

        d_exp.x = dot_product(cost_array[i], -sint_array[i]*cosp_array[i], sint_array[i]*sinp_array[i], dirx, ax, bx) ;
	    d_exp.y = dot_product(cost_array[i], -sint_array[i]*cosp_array[i], sint_array[i]*sinp_array[i], diry, ay, by) ;
	    d_exp.z = dot_product(cost_array[i], -sint_array[i]*cosp_array[i], sint_array[i]*sinp_array[i], dirz, az, bz) ;
	
        a_exp.x = dot_product(sint_array[i], cost_array[i]*cosp_array[i], -cost_array[i]*sinp_array[i], dirx, ax, bx) ;
        a_exp.y = dot_product(sint_array[i], cost_array[i]*cosp_array[i], -cost_array[i]*sinp_array[i], diry, ay, by) ;
        a_exp.z = dot_product(sint_array[i], cost_array[i]*cosp_array[i], -cost_array[i]*sinp_array[i], dirz, az, bz) ;
	
        b_exp.x = dot_product(0, sint_array[i], cosp_array[i], dirx, ax, bx) ;
        b_exp.y = dot_product(0, sint_array[i], cosp_array[i], diry, ay, by) ;
        b_exp.z = dot_product(0, sint_array[i], cosp_array[i], dirz, az, bz) ;
    
    
	    temp= rsqrtf(d_exp.x*d_exp.x + d_exp.y*d_exp.y + d_exp.z*d_exp.z);
	    d_exp.x=d_exp.x*temp;
	    d_exp.y=d_exp.y*temp;
	    d_exp.z=d_exp.z*temp;
	
	    temp= rsqrtf(a_exp.x*a_exp.x + a_exp.y*a_exp.y + a_exp.z*a_exp.z);
	    a_exp.x=a_exp.x*temp;
	    a_exp.y=a_exp.y*temp;
	    a_exp.z=a_exp.z*temp;
	
        temp= rsqrtf(b_exp.x*b_exp.x + b_exp.y*b_exp.y + b_exp.z*b_exp.z);
	    b_exp.x=b_exp.x*temp;
	    b_exp.y=b_exp.y*temp;
	    b_exp.z=b_exp.z*temp;
        
        d_calc.x = dirx_array[i];
        d_calc.y = diry_array[i];
        d_calc.z = dirz_array[i];

        a_calc.x = ax_array[i];
        a_calc.y = ay_array[i];
        a_calc.z = az_array[i];
        
        b_calc.x = bx_array[i];
        b_calc.y = by_array[i];
        b_calc.z = bz_array[i];
        
        if (compare_to_expected(results,"dir", d_calc, d_exp) != 0) testcase2fail = true;
        if (compare_to_expected(results,"a", a_calc, a_exp) != 0) testcase2fail = true;
        if (compare_to_expected(results,"b", b_calc, b_exp) != 0) testcase2fail = true;	
        //printf("expected d ( %f, %f, %f ) , expected a ( %f, %f, %f ), expected b ( %f, %f, %f )\n", 
    
        //printf("final d ( %f, %f, %f) , final a ( %f, %f, %f ), final b ( %f, %f, %f )\n", dirx_array[i], diry_array[i], dirz_array[i],
        //ax_array[i], ay_array[i], az_array[i], bx_array[i], by_array[i], bz_array[i]);
	    //printf("sint %f , sinp %f, cost %f, cosp %f\n", sint_array[i], sinp_array[i], cost_array[i], cosp_array[i]);
    }
    if (testcase2fail) 
        printf ("Spin Test 2 failed - tested expected value of direction, a, and b unit vectors, difference between vectors found\n");
    else 
        printf ("Spin Test 2 passed - tested expected value of direction, a, and b unit vectors, no difference between vectors found\n");    
	
    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(xd_mem_obj);
	ret = clReleaseMemObject(cd_mem_obj);
	ret = clReleaseMemObject(ad_mem_obj);
	ret = clReleaseCommandQueue(command_queue);
	return 0;
}

int reflect_new(unsigned int* x,unsigned int* c,unsigned int* a)
{
    unsigned int num[NUM_THREADS];
	unsigned long long num_tot, hist_tot;
	unsigned int hist[TEMP];
	unsigned int i;
    float dirx_array[NUM_THREADS];
    float diry_array[NUM_THREADS];
    float dirz_array[NUM_THREADS];
    
    float ax_array[NUM_THREADS];
    float ay_array[NUM_THREADS];
    float az_array[NUM_THREADS];
    
    float bx_array[NUM_THREADS];
    float by_array[NUM_THREADS];
    float bz_array[NUM_THREADS];

	unsigned int numh[NUM_THREADS];
	unsigned int histh[TEMP];
	
    FILE *fp;
    char *source_str;
    size_t source_size;

    initialize_vectors_reflect();

#ifdef LINUX
    fp = fopen("SRC/CLMonte_SRC/UnitTests/CLMonteUT.cl", "r");
#endif
#ifdef NAIF
    fp = fopen("C:\\Users\\Naif\\Documents\\ece496\\SRC\\CLMonte_SRC\\UnitTests\\CLMonteUT.cl", "r");
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

#ifdef LINUX    
    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &ret_num_devices);
#endif
#ifdef NAIF    
    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_CPU, 1, &device_id, &ret_num_devices);
#endif

	// Create an OpenCL context
    cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);

    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

    
    cl_mem dirx_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create dirx buff fail xd\n");
   
    
    cl_mem diry_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret,"create diry buff fail xd\n");
    
    
    cl_mem dirz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create dirz buff fail xd\n");
 
 	
    cl_mem ax_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create ax buff fail xd\n");
  
    
    cl_mem ay_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret,"create ay buff fail xd\n");
    
    
    cl_mem az_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create az buff fail xd\n");
   
    
    cl_mem bx_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create bx buff fail xd\n");
  
    
    cl_mem by_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret,"create by buff fail xd\n");
    
    
    cl_mem bz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_THREADS*sizeof(float), NULL, &ret);
    check_return(ret, "create bz buff fail xd\n");
    
    // Create memory buffers on the device for each vector 
    cl_mem xd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create xd buff fail\n");

    
    cl_mem cd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create cd buff fail\n");

    cl_mem ad_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create ad buff fail\n");

    
    ret = clEnqueueWriteBuffer(command_queue, xd_mem_obj, CL_TRUE, 0,size, xtest, 0, NULL, NULL);
    check_return(ret, "write buff xd fail\n");

	ret = clEnqueueWriteBuffer(command_queue, dirx_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), dirxtest, 0, NULL, NULL);
    check_return(ret,"write buff dirx fail \n");

	ret = clEnqueueWriteBuffer(command_queue, diry_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirytest, 0, NULL, NULL);
    check_return(ret,"write buff diry fail \n");

	ret = clEnqueueWriteBuffer(command_queue, dirz_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), dirztest, 0, NULL, NULL);
    check_return(ret,"write buff dirz fail \n");
	
    ret = clEnqueueWriteBuffer(command_queue, ax_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), axtest, 0, NULL, NULL);
    check_return(ret,"write buff ax fail \n");

	ret = clEnqueueWriteBuffer(command_queue, ay_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), aytest, 0, NULL, NULL);
    check_return(ret,"write buff ay fail \n");

	ret = clEnqueueWriteBuffer(command_queue, az_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), aztest, 0, NULL, NULL);
    check_return(ret,"write buff az fail \n");

    ret = clEnqueueWriteBuffer(command_queue, bx_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), bxtest, 0, NULL, NULL);
    check_return(ret,"write buff bx fail \n");

	ret = clEnqueueWriteBuffer(command_queue, by_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), bytest, 0, NULL, NULL);
    check_return(ret,"write buff by fail \n");

	ret = clEnqueueWriteBuffer(command_queue, bz_mem_obj, CL_TRUE, 0,NUM_THREADS*sizeof(float), bztest, 0, NULL, NULL);
    check_return(ret,"write buff bz fail \n");
	
    ret = clEnqueueWriteBuffer(command_queue, ad_mem_obj, CL_TRUE, 0, size, atest, 0, NULL, NULL);
    check_return(ret,"write buff ad fail \n");

	ret = clEnqueueWriteBuffer(command_queue, cd_mem_obj, CL_TRUE, 0, size, ctest, 0, NULL, NULL);
    check_return(ret,"write buff cd fail \n");


	

	// Create a program from the kernel source
    cl_program program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);

#ifdef LINUX
    ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
#endif
#ifdef NAIF
    ret = clBuildProgram(program, 1, &device_id, "-g -s \"C:\\Users\\Naif\\Documents\\ece496\\SRC\\CLMonte_SRC\\UnitTests\\CLMonteUT.cl\"", NULL, NULL);
#endif
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
		printf("Build Program fail\n");
		return -1;
	}
	
    // Create the OpenCL kernel
    cl_kernel kernel = clCreateKernel(program, "Reflect_new", &ret);
    check_return(ret, "create kernel Reflect_new fail"); 

    // Set the arguments of the kernel
    
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&dirx_mem_obj);
    check_return(ret, "set arg 0 fail\n"); 

    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&diry_mem_obj);
    check_return(ret, "set arg 1 fail\n"); 
    
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&dirz_mem_obj);
    check_return(ret, "set arg 2 fail\n"); 
    
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&ax_mem_obj);
    check_return(ret, "set arg 3 fail\n"); 

    ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&ay_mem_obj);
    check_return(ret, "set arg 4 fail\n"); 
    
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&az_mem_obj);
    check_return(ret, "set arg 5 fail\n"); 
    
    ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&bx_mem_obj);
    check_return(ret, "set arg 6 fail\n"); 

    ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&by_mem_obj);
    check_return(ret, "set arg 7 fail\n"); 
    
    ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&bz_mem_obj);
    check_return(ret, "set arg 8 fail\n"); 
    
    ret = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&xd_mem_obj);
    check_return(ret, "set arg 9 fail\n"); 
    
	ret = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&ad_mem_obj);
    check_return(ret, "set arg 10 fail\n"); 
	
	ret = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *)&cd_mem_obj);
    check_return(ret, "set arg 11 fail\n");
     
    // Execute the OpenCL kernel on the list
	cl_event kern_event;
    size_t global_item_size = NUM_THREADS; // Process the entire lists
    size_t local_item_size = NUM_THREADS_PER_BLOCK; // Process in groups of NUM_THREADS_PER_BLOCk
    printf("Running code on Accelerator in OpenCL\n");
	fprintf(print_file, "Running code on Accelerator in OpenCL\n");
	ret = clFlush(command_queue);
    ret = clFinish(command_queue);
	ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, &kern_event);
    check_return(ret, "EnqueuendRangeFail\n"); 
	
	ret = clFlush(command_queue);
    if(ret!=CL_SUCCESS){
        printf("%d bro 1\n", ret);
        exit(-1);
    }
    ret = clFinish(command_queue);
//    check_return(ret, "EnqueuendRangeFail\n"); 
//    check_return(ret, "read buff dirx fail\n"); 	
	
    if(ret!=CL_SUCCESS){
        printf("%d bro 2\n", ret);
        exit(-1);
    }
//exit(0);	
	ret = clEnqueueReadBuffer(command_queue, dirx_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirx_array, 1, &kern_event, NULL);
    check_return(ret, "read buff dirx fail\n"); 	
//exit(0);	

	ret = clEnqueueReadBuffer(command_queue, diry_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), diry_array, 1, &kern_event, NULL);
    check_return(ret, "read buff diry fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, dirz_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), dirz_array, 1, &kern_event, NULL);
    check_return(ret, "read buff dirz fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, ax_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), ax_array, 1, &kern_event, NULL);
    check_return(ret, "read buff ax fail\n"); 	
	
	ret = clEnqueueReadBuffer(command_queue, ay_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), ay_array, 1, &kern_event, NULL);
    check_return(ret, "read buff ay fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, az_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), az_array, 1, &kern_event, NULL);
    check_return(ret, "read buff az fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, bx_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), bx_array, 1, &kern_event, NULL);
    check_return(ret, "read buff bx fail\n"); 	
	
	ret = clEnqueueReadBuffer(command_queue, by_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), by_array, 1, &kern_event, NULL);
    check_return(ret, "read buff by fail\n"); 
	
	ret = clEnqueueReadBuffer(command_queue, bz_mem_obj, CL_TRUE, 0, NUM_THREADS*sizeof(float), bz_array, 1, &kern_event, NULL);
    check_return(ret, "read buff bz fail\n"); 
	
    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
	/*
	double tot_reflected = 0.0f;
	double measured_r;
    int j;
	for(i=0; i < NUM_THREADS; i++){
		if(dirz_array[i] == -dirztest[i])
		  tot_reflected+=1.0f; 
	}


	float sinangle_i = sqrtf(1.0f-dirz*dirz);
	float sinangle_t = N*sinangle_i;
	float cosangle_t = sqrtf(1.0f-sinangle_t*sinangle_t);
	
	
	float cossumangle = (-dirz*cosangle_t) - sinangle_i*sinangle_t;
	float cosdiffangle = (-dirz*cosangle_t) + sinangle_i*sinangle_t;
	float sinsumangle = sinangle_i*cosangle_t + (-dirz*sinangle_t);
	float sindiffangle = sinangle_i*cosangle_t - (-dirz*sinangle_t); 

	double r = 0.5*sindiffangle*sindiffangle* ((cosdiffangle*cosdiffangle+cossumangle*cossumangle)/(sinsumangle*sinsumangle*cosdiffangle*cosdiffangle));

	measured_r = tot_reflected/NUM_THREADS;


	printf("Reflect Test 1 - r %f,measured_r %f\n", r, measured_r); 	
	*/
    
    FILE *results;
    results = fopen("reflecttest.txt","w");
    float3 a_exp;
    float3 b_exp;
    float3 d_exp;
    float3 d_calc;
    float3 a_calc;
    float3 b_calc;

    float temp;
    float3 temp3;
    bool testcasefail = false;
    for (i=0; i< NUM_THREADS; i++) {
        if(dirz_array[i] == -dirztest[i]) {
            temp = sqrtf(dot_product(dirx_array[i], diry_array[i], dirz_array[i], dirx_array[i], diry_array[i], dirz_array[i]));
            if (compare_to_expected_float(results,"dir unit vector", temp, 1.0f) != 0) testcasefail = true;
        
            temp = sqrtf(dot_product(ax_array[i], ay_array[i], az_array[i], ax_array[i], ay_array[i], az_array[i]));
            if (compare_to_expected_float(results,"a unit vector", temp, 1.0f) != 0) testcasefail = true;
        
            temp = sqrtf(dot_product(bx_array[i], by_array[i], bz_array[i], bx_array[i], by_array[i], bz_array[i]));
            if (compare_to_expected_float(results,"b unit vector", temp, 1.0f) != 0) testcasefail = true;
        
            //test if d is orthog to a
            cross_product(dirx_array[i], diry_array[i], dirz_array[i], ax_array[i], ay_array[i], az_array[i], &temp3.x, &temp3.y, &temp3.z);
            temp = sqrtf(dot_product(temp3.x, temp3.y, temp3.z, temp3.x, temp3.y, temp3.z));
            if (compare_to_expected_float(results,"dir cross a", temp, 1.0f) != 0) testcasefail = true;

            //test if d is orthog to b
            cross_product(dirx_array[i], diry_array[i], dirz_array[i], bx_array[i], by_array[i], bz_array[i], &temp3.x, &temp3.y, &temp3.z);
            temp = sqrtf(dot_product(temp3.x, temp3.y, temp3.z, temp3.x, temp3.y, temp3.z));
            if (compare_to_expected_float(results,"dir cross b", temp, 1.0f) != 0) testcasefail = true;
	
            //test if a is orthog to b
            cross_product(ax_array[i], ay_array[i], az_array[i], bx_array[i], by_array[i], bz_array[i], &temp3.x, &temp3.y, &temp3.z);
            temp = sqrtf(dot_product(temp3.x, temp3.y, temp3.z, temp3.x, temp3.y, temp3.z));
            if (compare_to_expected_float(results,"a cross b", temp, 1.0f) != 0) testcasefail = true;
        }
    }
    if (testcasefail) 
        printf ("Reflect Test 2 failed - tested that direction, a, and b are unit vectors, and orthogonal. Inconsistency found\n");
    else 
        printf ("Reflect Test 2 passed - tested that direction, a, and b are unit vectors, and orthogonal. No inconsistency found\n");
	
    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(xd_mem_obj);
	ret = clReleaseMemObject(cd_mem_obj);
	ret = clReleaseMemObject(ad_mem_obj);
	ret = clReleaseCommandQueue(command_queue);
	return 0;
}

int MC(unsigned int* x,unsigned int* c,unsigned int* a){
    //reflect(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, DIRZ, x, c, a);
    //spin(sqrt(0.5f), 0.5f, 0.5f, x, c, a, 0.5f);
    //spin_new(0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, x, c ,a, 1.0f);
    reflect_new(x,c,a);
    //spinCPU(sqrt(0.5f), 0.5f, 0.5f, x, c, a, 0.5f);
	return 0;

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
#ifdef NAIF
	fp=fopen("C:\\Users\\Naif\\Documents\\ece496\\SRC\\HelperFiles\\safeprimes_base32.txt", "r");
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
