/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   GPU-based Monte Carlo simulation of photon migration in multi-layered media (GPU-MCML)
//   Copyright (C) 2010
//	
//   || DEVELOPMENT TEAM: 
//   --------------------------------------------------------------------------------------------------
//   Erik Alerstam, David Han, and William C. Y. Lo
//   
//   This code is the result of the collaborative efforts between 
//   Lund University and the University of Toronto.  
//
//   || DOCUMENTATION AND USER MANUAL: 
//   --------------------------------------------------------------------------------------------------
//	 Detailed "Wiki" style documentation is being developed for GPU-MCML 
//   and will be available on our webpage soon:
//   http://code.google.com/p/gpumcml 
// 
//   || NEW FEATURES: 
//   --------------------------------------------------------------------------------------------------
//    - Supports the Fermi GPU architecture 
//    - Backward compatible with pre-Fermi graphics cards
//    - Supports linux and Windows environment (Visual Studio)
//   
//   || PREVIOUS WORK: 
//   --------------------------------------------------------------------------------------------------
//	 This code is the fusion of our earlier, preliminary implementations and combines the best features 
//   from each implementation.  
//
//   W. C. Y. Lo, T. D. Han, J. Rose, and L. Lilge, "GPU-accelerated Monte Carlo simulation for photodynamic
//   therapy treatment planning," in Proc. of SPIE-OSA Biomedical Optics, vol. 7373.
//   
//   and 
//
//   http://www.atomic.physics.lu.se/biophotonics/our_research/monte_carlo_simulations/gpu_monte_carlo/
//	 E. Alerstam, T. Svensson and S. Andersson-Engels, "Parallel computing with graphics processing
//	 units for high-speed Monte Carlo simulations of photon migration", Journal of Biomedical Optics
//	 Letters, 13(6) 060504 (2008).
//
//   || CITATION: 
//   --------------------------------------------------------------------------------------------------
//	 We encourage the use, and modification of this code, and hope it will help 
//	 users/programmers to utilize the power of GPGPU for their simulation needs. While we
//	 don't have a scientific publication describing this code yet, we would very much appreciate it
//	 if you cite our original papers above if you use this code or derivations 
//   thereof for your own scientific work
//
//	 To compile and run this code, please visit www.nvidia.com and download the necessary 
//	 CUDA Toolkit, SDK, and Developer Drivers 
//
//	 If you use Visual Studio, the express edition is available for free at 
//   http://www.microsoft.com/express/Downloads/). 
//  	
//   This code is distributed under the terms of the GNU General Public Licence (see below). 
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/*	 
*   This file is part of GPUMCML.
* 
*   GPUMCML is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   GPUMCML is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with GPUMCML.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <float.h> //for FLT_MAX 
#include <stdio.h>
#include <time.h> 
#include <string.h>
#include <math.h>

//#include <cuda_runtime.h>

//#ifdef _WIN32 
//#include "gpumcml_io.c"
//#include "cutil-win32/cutil.h"
//#else 
//#include <cutil.h>
//#endif

#include <CL/cl.h>
#include "kernel.h"
#include "gpumcml.h"
#define MAX_SOURCE_SIZE 0x100000

cl_context context;
cl_command_queue command_queue;
cl_mem simparam_mem_obj;
cl_mem layerspecs_mem_obj;
cl_mem num_photons_left_mem_obj;
cl_mem a_mem_obj;
cl_mem x_mem_obj;
cl_mem A_rz_mem_obj;
cl_mem Rd_ra_mem_obj;
cl_mem Tt_ra_mem_obj;
cl_mem photon_x_mem_obj;
cl_mem photon_y_mem_obj;
cl_mem photon_z_mem_obj;
cl_mem photon_ux_mem_obj;
cl_mem photon_uy_mem_obj;
cl_mem photon_uz_mem_obj;
cl_mem photon_w_mem_obj;
cl_mem photon_sleft_mem_obj;
cl_mem photon_layer_mem_obj;
cl_mem is_active_mem_obj;
cl_mem tetra_mesh_mem_obj;
cl_mem materials_mem_obj;
cl_kernel initkernel;
cl_kernel kernel;
cl_program program;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//   Initialize random number generator 
//////////////////////////////////////////////////////////////////////////////
int init_RNG(UINT64 *x, UINT32 *a, 
             const UINT32 n_rng, const char *safeprimes_file, UINT64 xinit)
{
  FILE *fp;
  UINT32 begin=0u;
  UINT32 fora,tmp1,tmp2;

  if (strlen(safeprimes_file) == 0)
  {
    // Try to find it in the local directory
    safeprimes_file = "safeprimes_base32.txt";
  }

  fp = fopen(safeprimes_file, "r");

  if(fp == NULL)
  {
    printf("Could not find the file of safeprimes (%s)! Terminating!\n", safeprimes_file);
    exit(-1);
    return 1;
  }

  fscanf(fp,"%u %u %u",&begin,&tmp1,&tmp2);

  // Here we set up a loop, using the first multiplier in the file to generate x's and c's
  // There are some restictions to these two numbers:
  // 0<=c<a and 0<=x<b, where a is the multiplier and b is the base (2^32)
  // also [x,c]=[0,0] and [b-1,a-1] are not allowed.

  //Make sure xinit is a valid seed (using the above mentioned restrictions)
  if((xinit == 0ull) | (((UINT32)(xinit>>32))>=(begin-1)) | (((cl_uint)xinit)>=0xfffffffful))
  {
    //xinit (probably) not a valid seed! (we have excluded a few unlikely exceptions)
    printf("%llu not a valid seed! Terminating!\n",xinit);
    return 1;
  }

  for (UINT32 i=0;i < n_rng;i++)
  {
    fscanf(fp,"%u %u %u",&fora,&tmp1,&tmp2);
    a[i]=fora;
    x[i]=0;
    while( (x[i]==0) | (((UINT32)(x[i]>>32))>=(fora-1)) | (((cl_uint)x[i])>=0xfffffffful))
    {
      //generate a random number
      xinit=(xinit&0xffffffffull)*(begin)+(xinit>>32);

      //calculate c and store in the upper 32 bits of x[i]
      x[i]=(UINT32) floor((((double)((cl_uint)xinit))/(double)0x100000000)*fora);//Make sure 0<=c<a
      x[i]=x[i]<<32;

      //generate a random number and store in the lower 32 bits of x[i] (as the initial x of the generator)
      xinit=(xinit&0xffffffffull)*(begin)+(xinit>>32);//x will be 0<=x<b, where b is the base 2^32
      x[i]+=(UINT32) xinit;
    }
    //if(i<10)printf("%llu\n",x[i]);
  }
  fclose(fp);

  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//   Supports 1 GPU only
//   Calls RunGPU with HostThreadState parameters
//////////////////////////////////////////////////////////////////////////////
int RunGPUi(HostThreadState *hstate)
{
  SimState *HostMem = &(hstate->host_sim_state);

  printf("HostMem num photons %d\n", *HostMem->n_photons_left); 
  FILE *fp;
  FILE *build_file;
  build_file = fopen("build.txt", "w");
  char *source_str;
  size_t source_size = MAX_SOURCE_SIZE;

  fp = fopen("kernel.cl", "r");
  if(!fp){
    fprintf(stderr, "Failed to read kernel file.\n");
    exit(-1);

  }
  source_str = (char*)malloc(MAX_SOURCE_SIZE);
  source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);

  cl_platform_id platform_id = NULL;
  cl_device_id device_id = NULL;
  UINT32 ret_num_devices;
  UINT32 ret_num_platforms;
  cl_int ret= clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
  if(ret != CL_SUCCESS){
    printf("get platform id fail, exiting\n");
    exit(-1);
  }
  
  ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &ret_num_devices);
  if(ret != CL_SUCCESS){
    printf("get device id fail, exiting\n");
    exit(-1);
  }
  

  context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);
  if(ret != CL_SUCCESS){
    printf("Create context fail, exiting\n");
    exit(-1);
  }
  
  command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
  if(ret != CL_SUCCESS){
    printf("Create Command Queue fail, exiting\n");
    exit(-1);
  }
  
  // Init the remaining states.
  InitSimStates(HostMem, hstate->sim, context, command_queue, &num_photons_left_mem_obj, 
          &a_mem_obj, &x_mem_obj, &A_rz_mem_obj, &Rd_ra_mem_obj, &Tt_ra_mem_obj, &photon_x_mem_obj, &photon_y_mem_obj,
          &photon_z_mem_obj, &photon_ux_mem_obj, &photon_uy_mem_obj, &photon_uz_mem_obj, &photon_w_mem_obj, 
          &photon_sleft_mem_obj, &photon_layer_mem_obj, &is_active_mem_obj);

  InitDCMem(hstate->sim, context, command_queue, &simparam_mem_obj, &layerspecs_mem_obj, &tetra_mesh_mem_obj, &materials_mem_obj);

  program = clCreateProgramWithSource(context, 1, (const char**)&source_str, (const size_t *)&source_size, &ret);
  if(ret != CL_SUCCESS){
    printf("create program fail, exiting\n");
    exit(-1);
  }
  
  ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
  
  if(ret!=CL_SUCCESS){
    printf("Build Program Fail, exiting\n");
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
    exit(-1);
  }

  kernel = clCreateKernel(program, "MCMLKernel", &ret);
  if(ret != CL_SUCCESS){
    printf("create mcml kernel fail, exiting\n");
    exit(-1);
  }
  
 // initkernel = clCreateKernel(program, "InitThreadState", &ret);
 // if(ret != CL_SUCCESS){
 //   printf("create initthreadstate kernel fail, exiting\n");
 //   exit(-1);
 // }
  
 // ret = clSetKernelArg(initkernel, 0, sizeof(cl_mem),(void *)&photon_x_mem_obj);
 // if(ret != CL_SUCCESS){
 //   printf("Error setting photon x init kernel argument, exiting\n");
 //   exit(-1);
 // }

 // ret = clSetKernelArg(initkernel, 1, sizeof(cl_mem),(void *)&photon_y_mem_obj);
 // if(ret != CL_SUCCESS){
 //   printf("Error setting photon y init kernel argument, exiting\n");
 //   exit(-1);
 // }

 // ret = clSetKernelArg(initkernel, 2, sizeof(cl_mem),(void *)&photon_z_mem_obj);
 // if(ret != CL_SUCCESS){
 //   printf("Error setting photon z init kernel argument, exiting\n");
 //   exit(-1);
 // }

 // ret = clSetKernelArg(initkernel, 3, sizeof(cl_mem),(void *)&photon_ux_mem_obj);
 // if(ret != CL_SUCCESS){
 //   printf("Error setting photon ux init kernel argument, exiting\n");
 //   exit(-1);
 // }

 // ret = clSetKernelArg(initkernel, 4, sizeof(cl_mem),(void *)&photon_uy_mem_obj);
 // if(ret != CL_SUCCESS){
 //   printf("Error setting photon uy init kernel argument, exiting\n");
 //   exit(-1);
 // }

 // ret = clSetKernelArg(initkernel, 5, sizeof(cl_mem),(void *)&photon_uz_mem_obj);
 // if(ret != CL_SUCCESS){
 //   printf("Error setting photon uz init kernel argument, exiting\n");
 //   exit(-1);
 // }

 // ret = clSetKernelArg(initkernel, 6, sizeof(cl_mem),(void *)&photon_w_mem_obj);
 // if(ret != CL_SUCCESS){
 //   printf("Error setting photon w init kernel argument, exiting\n");
 //   exit(-1);
 // }

 // ret = clSetKernelArg(initkernel, 7, sizeof(cl_mem),(void *)&photon_sleft_mem_obj);
 // if(ret != CL_SUCCESS){
 //   printf("Error setting photon sleft init kernel argument, exiting\n");
 //   exit(-1);
 // }

 // ret = clSetKernelArg(initkernel, 8, sizeof(cl_mem),(void *)&photon_layer_mem_obj);
 // if(ret != CL_SUCCESS){
 //   printf("Error setting photon layer init kernel argument, exiting\n");
 //   exit(-1);
 // }

 // ret = clSetKernelArg(initkernel, 9, sizeof(cl_mem),(void *)&is_active_mem_obj);
 // if(ret != CL_SUCCESS){
 //   printf("Error setting is active init kernel argument, exiting\n");
 //   exit(-1);
 // }

 // ret = clSetKernelArg(initkernel, 10, sizeof(cl_mem),(void *)&simparam_mem_obj);
 // if(ret != CL_SUCCESS){
 //   printf("Error setting simparam init kernel argument, exiting\n");
 //   exit(-1);
 // }
  int argnum = 0;
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&simparam_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting simparam kernel argument, exiting\n");
    exit(-1);
  }
  
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&layerspecs_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting layerspecs kernel argument, exiting\n");
    exit(-1);
  }
  
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&(num_photons_left_mem_obj));
  if(ret != CL_SUCCESS){
    printf("Error setting photons left kernel argument, exiting\n");
    exit(-1);
  }
  
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&x_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting x kernel argument, exiting\n");
    exit(-1);
  }
  
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&a_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting a kernel argument, exiting\n");
    exit(-1);
  }
  
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&Rd_ra_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting Rd_ra kernel argument, exiting\n");
    exit(-1);
  }
  
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&A_rz_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting A_rz kernel argument, exiting\n");
    exit(-1);
  }
  
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&Tt_ra_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting Tt_ra kernel argument, exiting\n");
    exit(-1);
  }
  
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem), (void *)&tetra_mesh_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting tetra mesh kernel argument, exiting\n");
    exit(-1);
  }
  
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem), (void *)&materials_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting materials kernel argument, exiting\n");
    exit(-1);
  }
  //ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&photon_x_mem_obj);
  //if(ret != CL_SUCCESS){
  //  printf("Error setting photon_x kernel argument, exiting\n");
  //  exit(-1);
  //}
  //
  //ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&photon_y_mem_obj);
  //if(ret != CL_SUCCESS){
  //  printf("Error setting photon_y kernel argument, exiting\n");
  //  exit(-1);
  //}
  //
  //ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&photon_z_mem_obj);
  //if(ret != CL_SUCCESS){
  //  printf("Error setting photon_z kernel argument, exiting\n");
  //  exit(-1);
  //}
  //
  //ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&photon_ux_mem_obj);
  //if(ret != CL_SUCCESS){
  //  printf("Error setting photon_ux kernel argument, exiting\n");
  //  exit(-1);
  //}
  //
  //ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&photon_uy_mem_obj);
  //if(ret != CL_SUCCESS){
  //  printf("Error setting photon_uy kernel argument, exiting\n");
  //  exit(-1);
  //}
  //
  //ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&photon_uz_mem_obj);
  //if(ret != CL_SUCCESS){
  //  printf("Error setting photon_uz kernel argument, exiting\n");
  //  exit(-1);
  //}
  //
  //ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&photon_w_mem_obj);
  //if(ret != CL_SUCCESS){
  //  printf("Error setting photon_w kernel argument, exiting\n");
  //  exit(-1);
  //}
  //
  //ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&photon_sleft_mem_obj);
  //if(ret != CL_SUCCESS){
  //  printf("Error setting photon_sleft kernel argument, exiting\n");
  //  exit(-1);
  //}
  //
  //ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&photon_layer_mem_obj);
  //if(ret != CL_SUCCESS){
  //  printf("Error setting photon_layer kernel argument, exiting\n");
  //  exit(-1);
  //}
  //
  //ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&is_active_mem_obj);
  //if(ret != CL_SUCCESS){
  //  printf("Error setting is_active kernel argument, exiting\n");
  //  exit(-1);
  //}

  size_t global_size = NUM_BLOCKS*NUM_THREADS_PER_BLOCK;
  size_t local_size = NUM_THREADS_PER_BLOCK;

//  ret = clEnqueueNDRangeKernel(command_queue, initkernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
//  if(ret != CL_SUCCESS){
//    printf("Error enqueundrange of initkernel, exiting\n");
//    exit(-1);
//  }
//
//   clFlush(command_queue);
//   clFinish(command_queue);
//
//    ret = clEnqueueReadBuffer(command_queue, num_photons_left_mem_obj, CL_TRUE, 0, sizeof(UINT32), HostMem->n_photons_left, 0, NULL, NULL);
//    if(ret != CL_SUCCESS){
//        printf("test read failed %d\n", ret);
//        exit(-1);
//
//    }
//    
//
//    printf("num photons left %d after initthreadstate\n", *HostMem->n_photons_left);
    //exit(0);
  // Initialize the remaining thread states.

  int i=0;
  for (i=0; *HostMem->n_photons_left > 0; ++i)
  {
    // Run the kernel.
    ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
    if(ret != CL_SUCCESS){
      printf("Error enqueundrange of kernel in loop iteration %d, exiting\n", i);
      exit(-1);
    }
    clFinish(command_queue); //good practice for synchronization purposes in OpenCL

    // Copy the number of photons left from device to host.
    ret = clEnqueueReadBuffer(command_queue, num_photons_left_mem_obj, CL_TRUE, 0, sizeof(UINT32), HostMem->n_photons_left, 0, NULL, NULL);
    if(ret != CL_SUCCESS){
      printf("Error reading number of photons left of kernel in loop iteration %d, exiting %d\n", i, ret);
      exit(-1);
    }

    ret = clEnqueueReadBuffer(command_queue, x_mem_obj, CL_TRUE, 0, NUM_THREADS * sizeof(UINT64), HostMem->x, 0, NULL, NULL);
    if(ret != CL_SUCCESS){
      printf("Error reading x buffer, exiting\n");
      exit(-1);
    }
    printf("x[1] = %E\n", (double)HostMem->x[1]);
  
    ret = clEnqueueReadBuffer(command_queue, a_mem_obj, CL_TRUE, 0, NUM_THREADS * sizeof(UINT32), HostMem->a, 0, NULL, NULL);
    if(ret != CL_SUCCESS){
      printf("Error reading a buffer, exiting\n");
      exit(-1);
    }
    printf("a[1] = %E\n", (double)HostMem->a[1]);
    ////////////////////////////////////////////////////////////

    printf("[GPU] batch %d, number of photons left %d\n",i, *HostMem->n_photons_left);
  }
  
  printf("[GPU] simulation done!\n");

  CopyDeviceToHostMem(HostMem, hstate->sim, command_queue, A_rz_mem_obj, Rd_ra_mem_obj, Tt_ra_mem_obj, x_mem_obj);
  FreeDeviceSimStates(context, command_queue, initkernel,kernel, program, simparam_mem_obj, layerspecs_mem_obj,num_photons_left_mem_obj, a_mem_obj, x_mem_obj, A_rz_mem_obj, Rd_ra_mem_obj, Tt_ra_mem_obj, photon_x_mem_obj, photon_y_mem_obj, photon_z_mem_obj, photon_ux_mem_obj, 
photon_uy_mem_obj, photon_uz_mem_obj, photon_w_mem_obj, photon_sleft_mem_obj, photon_layer_mem_obj, is_active_mem_obj, tetra_mesh_mem_obj, materials_mem_obj);
  // We still need the host-side structure.
  return i;
}

//////////////////////////////////////////////////////////////////////////////
//   Perform MCML simulation for one run out of N runs (in the input file)
//////////////////////////////////////////////////////////////////////////////
static void DoOneSimulation(int sim_id, SimulationStruct* simulation,
                            unsigned long long *x, unsigned int *a)
{
  printf("\n------------------------------------------------------------\n");
  printf("        Simulation #%d\n", sim_id);
  printf("        - number_of_photons = %u\n", simulation->number_of_photons);
  printf("------------------------------------------------------------\n\n");

  // Start simulation kernel exec timer
  clock_t start, end;
  double time_elapsed;
  start = clock();
  // For each GPU, init the host-side structure.
  HostThreadState* hstates;
  hstates = (HostThreadState*)malloc(sizeof(HostThreadState));

  hstates->sim = simulation;

  SimState *hss = &(hstates->host_sim_state);

  // number of photons responsible 
  hss->n_photons_left = (int*)malloc(sizeof(int)); //change to int
  *(hss->n_photons_left) = simulation->number_of_photons; 

  // random number seeds
  hss->x = &x[0]; hss->a = &a[0];

  printf("simulation number of photons %d\n", *(hss->n_photons_left));
  // Launch simulation
  int number_of_iterations = RunGPUi (hstates);

  // End the timer.
  //CUT_SAFE_CALL( cutStopTimer(execTimer) );
  end = clock();
  //float elapsedTime = cutGetTimerValue(execTimer);
  float elapsedTime = ((float) end - start)/CLOCKS_PER_SEC;
  float total_steps = (float)number_of_iterations*(float)NUM_STEPS*(float)NUM_THREADS;
  printf( "\n\n>>>>>>Simulation time: %f (s)\n", elapsedTime);
  printf("total num of iterations = %d\n", number_of_iterations);
  printf("NUM_STEPS = %d\n", NUM_STEPS);
  printf("NUM_THREADS = %d\n", NUM_THREADS);
  printf( ">>>>>>Simulation Speed: %e photon events per second\n", total_steps/elapsedTime);
  
  Write_Simulation_Results(hss, simulation, elapsedTime);

  //CUT_SAFE_CALL( cutDeleteTimer(execTimer) );

  // Free SimState structs.
  FreeHostSimState(hss);
  free(hstates);
}

void banner()
{
    printf("FullMonte OpenCL v0.0\n");
    printf("(c) Yu Wu, Emil Salavat, Li Chen, 2015\n\n");
}

//////////////////////////////////////////////////////////////////////////////
//   Perform MCML simulation for one run out of N runs (in the input file)
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  banner();

  Point *points;
  Tetra *tetra_mesh;
  unsigned long Np, Nt;	//number of points, number of tetrahedra
  //PopulateTetraFromMeshFile("", *points, tetra_mesh, &Np, &Nt);

  char* filename = NULL;
  unsigned long long seed = (unsigned long long) time(NULL);
  int ignoreAdetection = 0;
  
  SimulationStruct* simulations;
  int n_simulations;

  int i;

  // Parse command-line arguments.
  if (interpret_arg(argc, argv, &filename,&seed, &ignoreAdetection))
  {
    usage(argv[0]);
    return 1;
  }

  // Output the execution configuration.
  printf("\n====================================\n");
  printf("EXECUTION MODE:\n");
  printf("  ignore A-detection:      %s\n", ignoreAdetection ? "YES" : "NO");
  printf("  seed:                    %llu\n", seed);
  printf("====================================\n\n");

  // Read the simulation inputs.
  n_simulations = read_simulation_data(filename, &simulations, ignoreAdetection);
  if(n_simulations == 0)
  {
    printf("Something wrong with read_simulation_data!\n");
    return 1;
  }
  printf("Read %d simulations\n",n_simulations);

  // Allocate and initialize RNG seeds.
  unsigned int len = NUM_THREADS;

  unsigned long long *x = (unsigned long long*)malloc(len * sizeof(unsigned long long));
  unsigned int *a = (unsigned int*)malloc(len * sizeof(unsigned int));

  if (init_RNG(x, a, len, "executable/safeprimes_base32.txt", seed)) return 1;
  
  printf("Using the MWC random number generator ...\n");

  //perform all the simulations
  for(i=0;i<n_simulations;i++)
  {
    // Run a simulation
    DoOneSimulation(i, &simulations[i], x, a);
  }

  // Free the random number seed arrays.
  free(x); free(a);
  FreeSimulationStruct(simulations, n_simulations);

  return 0; 
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

