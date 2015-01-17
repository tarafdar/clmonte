#include <float.h> //for FLT_MAX 
#include <stdio.h>
#include <time.h> 
#include <string.h>
#include <math.h>
#include <CL/cl.h>
#include "kernel.h"
#include "gpumcml.h"
#define MAX_SOURCE_SIZE 0x100000

cl_context context;
cl_command_queue command_queue;
float *debug;
cl_mem simparam_mem_obj;
cl_mem num_photons_left_mem_obj;
cl_mem a_mem_obj;
cl_mem x_mem_obj;
cl_mem debug_mem_obj;
cl_mem tetra_mesh_mem_obj;
cl_mem materials_mem_obj;
cl_mem scaled_w_mem_obj;
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
int RunGPUi(HostThreadState *hstate, Tetra *tetra_mesh, Material *materialspec)
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
          &a_mem_obj, &x_mem_obj, &scaled_w_mem_obj, &debug_mem_obj);

  InitDCMem(hstate->sim, tetra_mesh, materialspec, context, command_queue, &simparam_mem_obj, &tetra_mesh_mem_obj, &materials_mem_obj);

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
  
  int argnum = 0;
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem),(void *)&simparam_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting simparam kernel argument, exiting\n");
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
  
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem), (void *)&scaled_w_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting scaled weight kernel argument, exiting\n");
    exit(-1);
  }
  
  ret = clSetKernelArg(kernel, argnum++, sizeof(cl_mem), (void *)&debug_mem_obj);
  if(ret != CL_SUCCESS){
    printf("Error setting debug kernel argument, exiting\n");
    exit(-1);
  }  

  size_t global_size = NUM_BLOCKS*NUM_THREADS_PER_BLOCK;
  size_t local_size = NUM_THREADS_PER_BLOCK;

  int i=0;
  for (i=0; *HostMem->n_photons_left > 0 && i<1; ++i)
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

  CopyDeviceToHostMem(HostMem, hstate->sim, command_queue, x_mem_obj, scaled_w_mem_obj, debug_mem_obj, tetra_mesh);
  FreeDeviceSimStates(context, command_queue, initkernel,kernel, program, simparam_mem_obj, num_photons_left_mem_obj, a_mem_obj, x_mem_obj, tetra_mesh_mem_obj, materials_mem_obj, scaled_w_mem_obj, debug_mem_obj);
  // We still need the host-side structure.
  return i;
}

//////////////////////////////////////////////////////////////////////////////
//   Perform MCML simulation for one run out of N runs (in the input file)
//////////////////////////////////////////////////////////////////////////////
static void DoOneSimulation(int sim_id, SimulationStruct* simulation,
                            unsigned long long *x, unsigned int *a, Tetra *tetra_mesh, Material *materialspec)
{
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
  int number_of_iterations = RunGPUi (hstates, tetra_mesh, materialspec);

  // End the timer.
  end = clock();
  float elapsedTime = ((float) end - start)/CLOCKS_PER_SEC;
  float total_steps = (float)number_of_iterations*(float)NUM_STEPS*(float)NUM_THREADS;
  printf( "\n\n>>>>>>Simulation time: %f (s)\n", elapsedTime);
  printf("total num of iterations = %d\n", number_of_iterations);
  printf("NUM_STEPS = %d\n", NUM_STEPS);
  printf("NUM_THREADS = %d\n", NUM_THREADS);
  printf( ">>>>>>Simulation Speed: %e photon events per second\n", total_steps/elapsedTime);
  
  Write_Simulation_Results(hss, simulation, elapsedTime);

  // Free SimState structs.
  FreeHostSimState(hss);
  free(hstates);
}

void banner()
{
    printf("FullMonte OpenCL v0.0\n");
    printf("(c) Yu Wu, Emil Salavat, Li Chen, 2015\n\n");
}

void OutputTetraMesh(Tetra *tetra_mesh, int Nt)	//debug code
{
  int i;
  for(i=1; i<=Nt; i++)
  {
    Tetra &t = tetra_mesh[i];
    printf("Tetra%d: material: %d\n", i, t.matID);
    printf("  %f %f %f %f %d\n", t.face[0][0], t.face[0][1], t.face[0][2], t.face[0][3], t.adjTetras[0]);
    printf("  %f %f %f %f %d\n", t.face[1][0], t.face[1][1], t.face[1][2], t.face[1][3], t.adjTetras[1]);
    printf("  %f %f %f %f %d\n", t.face[2][0], t.face[2][1], t.face[2][2], t.face[2][3], t.adjTetras[2]);
    printf("  %f %f %f %f %d\n", t.face[3][0], t.face[3][1], t.face[3][2], t.face[3][3], t.adjTetras[3]);
  }
}

//////////////////////////////////////////////////////////////////////////////
//   Perform MCML simulation for one run out of N runs (in the input file)
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  banner();

  Tetra *tetra_mesh;
  int Np, Nt, Nm;	//number of points, number of tetrahedra, number of materials
  PopulateTetraFromMeshFile("one_layer_18_18_1_2.mesh", &tetra_mesh, &Np, &Nt);
  Material *materialspec;
  OutputTetraMesh(tetra_mesh, Nt);	//debug code

  Source *sourcepoint;
  //ParseSource("FourLayer.source", &sourcepoint);

  Material *mat;
  //ParseMaterial("FourLayer.opt", &mat);

  PopulateMaterialFromInput("", &materialspec, &Nm);
  printf("mat n: %f\n", materialspec[1].n);
  char* filename = NULL;
  unsigned long long seed = (unsigned long long) time(NULL);
  
  SimulationStruct simulation;

  int i;

  // Parse command-line arguments.
  if (interpret_arg(argc, argv, &simulation))
  {
    usage(argv[0]);
    return 1;
  }
  // Output the execution configuration.
  printf("\n====================================\n");
  printf("EXECUTION MODE:\n");
  printf("  seed:                    %llu\n", seed);
  printf("====================================\n\n");

  // Read the simulation inputs.
  printf("Number of photons: %d\n", simulation.number_of_photons);

  // Allocate and initialize RNG seeds.
  unsigned int len = NUM_THREADS;

  unsigned long long *x = (unsigned long long*)malloc(len * sizeof(unsigned long long));
  unsigned int *a = (unsigned int*)malloc(len * sizeof(unsigned int));

  if (init_RNG(x, a, len, "executable/safeprimes_base32.txt", seed)) return 1;
  
  printf("Using the MWC random number generator ...\n");

  //perform the simulation
  simulation.nTetras = Nt;
  simulation.nMaterials = Nm;
  // Run a simulation
  DoOneSimulation(i, &simulation, x, a, tetra_mesh, materialspec);

  // Free the random number seed arrays.
  free(materialspec);
  free(tetra_mesh);
  free(x); free(a);

  return 0; 
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

