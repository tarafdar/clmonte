
#include <stdio.h>
#include <math.h>
#include <CL/cl.h>
#include "kernel.h"
#include "gpumcml.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//   Initialize Device Constant Memory with read-only data
//////////////////////////////////////////////////////////////////////////////
extern float *debug;
int InitDCMem(SimulationStruct *sim, Tetra *tetra_mesh, Material *materialspec, cl_context context, cl_command_queue command_queue, cl_mem *simparam_mem_obj, cl_mem *tetra_mesh_mem_obj, cl_mem *materials_mem_obj)
{
  SimParamGPU h_simparam;

  h_simparam.dz = sim->det.dz;
  h_simparam.dr = sim->det.dr;
  
  
  h_simparam.na = sim->det.na;
  h_simparam.nz = sim->det.nz;
  h_simparam.nr = sim->det.nr;

  h_simparam.originX = 0;
  h_simparam.originY = 0;
  h_simparam.originZ = 0;
  h_simparam.init_tetraID = 1;
  h_simparam.terminationThresh = 0.5;
  h_simparam.proulettewin = 0.5;
  h_simparam.weight_scale = 10000000;

  cl_int ret;
  *simparam_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(SimParamGPU), NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating simparam buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueWriteBuffer(command_queue, *simparam_mem_obj, CL_TRUE, 0, sizeof(SimParamGPU), &h_simparam, 0, NULL, NULL); 
  ////CUDA_SAFE_CALL( cudaMemcpyToSymbol(d_simparam,
  ////  &h_simparam, sizeof(SimParamGPU)) );
  //if(ret!= CL_SUCCESS){
  //  printf("Error writing to simparam buffer, exiting\n");
  //  exit(-1);
  //}

  *tetra_mesh_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Tetra)*(sim->nTetras+1), NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error create tetra mesh buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueWriteBuffer(command_queue, *tetra_mesh_mem_obj, CL_TRUE, 0, sizeof(Tetra)*(sim->nTetras+1), tetra_mesh, 0, NULL, NULL);
  // Copy tetra mesh data to constant device memory
  if(ret!= CL_SUCCESS){
    printf("Error writing to tetraspecs buffer, exiting\n");
    exit(-1);
  }

  *materials_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Material)*(sim->nMaterials+1), NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error create materials buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueWriteBuffer(command_queue, *materials_mem_obj, CL_TRUE, 0, sizeof(Material)*(sim->nMaterials+1), materialspec, 0, NULL, NULL);
  // Copy material data to constant device memory
  if(ret!= CL_SUCCESS){
    printf("Error writing to materialspecs buffer, exiting\n");
    exit(-1);
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//   Initialize Device Memory (global) for read/write data
//////////////////////////////////////////////////////////////////////////////
int InitSimStates(SimState* HostMem, SimulationStruct* sim, cl_context context, cl_command_queue command_queue, 
        cl_mem *num_photons_left_mem_obj, cl_mem *a_mem_obj, cl_mem *x_mem_obj, cl_mem *A_rz_mem_obj, cl_mem *Rd_ra_mem_obj, cl_mem *Tt_ra_mem_obj, 
        cl_mem *photon_x_mem_obj ,cl_mem *photon_y_mem_obj, cl_mem *photon_z_mem_obj, cl_mem *photon_ux_mem_obj, 
        cl_mem *photon_uy_mem_obj, cl_mem *photon_uz_mem_obj, cl_mem *photon_w_mem_obj, cl_mem *photon_sleft_mem_obj, 
        cl_mem *photon_layer_mem_obj, cl_mem *is_active_mem_obj, cl_mem *scaled_w_mem_obj, cl_mem *debug_mem_obj
        )
{
  int rz_size = sim->det.nr * sim->det.nz;
  int ra_size = sim->det.nr * sim->det.na;
  unsigned int size;
  cl_int ret;
  size = sizeof(UINT32); 
  // Allocate n_photons_left (on device only)
  *num_photons_left_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photons left buffer, exiting\n");
    exit(-1);
  }

  ret = clEnqueueWriteBuffer(command_queue, *num_photons_left_mem_obj, CL_TRUE, 0, size, HostMem->n_photons_left, 0, NULL, NULL); 
  if(ret!= CL_SUCCESS){
    printf("Error writing to a mem buffer, exiting\n");
    exit(-1);
  }

  // random number generation (on device only)
  size = NUM_THREADS * sizeof(UINT32);

  *a_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating a buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueWriteBuffer(command_queue, *a_mem_obj, CL_TRUE, 0, size, HostMem->a, 0, NULL, NULL); 
  if(ret!= CL_SUCCESS){
    printf("Error writing to a mem buffer, exiting\n");
    exit(-1);
  }

  size = NUM_THREADS * sizeof(UINT64);
  
  *x_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating x buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueWriteBuffer(command_queue, *x_mem_obj, CL_TRUE, 0, size, HostMem->x, 0, NULL, NULL); 
  if(ret!= CL_SUCCESS){
    printf("Error writing to x mem buffer, exiting\n");
    exit(-1);
  }

  size = (sim->nTetras+1) * sizeof(UINT64);
  *scaled_w_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!=CL_SUCCESS){
    printf("Error creating scaled weight buffer, exiting\n");
    exit(-1);
  }
  HostMem->scaled_w = (UINT64*)malloc(size);
  ret = clEnqueueWriteBuffer(command_queue, *scaled_w_mem_obj, CL_TRUE, 0, size, HostMem->scaled_w, 0, NULL, NULL);
  if(ret!=CL_SUCCESS){
    printf("Error writing to scaled weight mem buffer, exiting\n");
    exit(-1);
  }

  // Allocate A_rz on host and device
  size = rz_size * sizeof(UINT64);
  
  *A_rz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating A_rz buffer, exiting\n");
    exit(-1);
  }
  
  HostMem->A_rz = (UINT64*)malloc(size);
  if (HostMem->A_rz == NULL)
  {
    fprintf(stderr, "Error allocating HostMem->A_rz");
    exit(1);
  }
  for(int i=0; i<rz_size; i++){
    HostMem->A_rz[i] = 0;
  }
  ret = clEnqueueWriteBuffer(command_queue, *A_rz_mem_obj, CL_TRUE, 0, size, HostMem->A_rz, 0, NULL, NULL); 
  
  // Allocate debug on host and device
  size = sizeof(float)*40;
  *debug_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating debug buffer, exiting\n");
    exit(-1);
  }
  debug = (float*)malloc(size);
  if (debug == NULL)
  {
    fprintf(stderr, "Error allocating debug");
    exit(1);
  }
  ret = clEnqueueWriteBuffer(command_queue, *debug_mem_obj, CL_TRUE, 0, size, debug, 0, NULL, NULL);   
  
  // On the device, we allocate multiple copies for less access contention.

  // Allocate Rd_ra on host and device
  
  size = ra_size * sizeof(UINT64);
  HostMem->Rd_ra = (UINT64*)malloc(size);
  if(HostMem->Rd_ra==NULL){printf("Error allocating HostMem->Rd_ra"); exit (1);}
  for(int i=0; i<ra_size; i++){
    HostMem->Rd_ra[i] = 0;
  }
  *Rd_ra_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating Rd_ra buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueWriteBuffer(command_queue, *Rd_ra_mem_obj, CL_TRUE, 0, size, HostMem->Rd_ra, 0, NULL, NULL); 
  

  // Allocate Tt_ra on host and device
  size = ra_size * sizeof(UINT64);
  HostMem->Tt_ra = (UINT64*)malloc(size);
  for(int i=0; i<ra_size; i++){
    HostMem->Tt_ra[i] = 0;
  }
  if(HostMem->Tt_ra==NULL){printf("Error allocating HostMem->Tt_ra"); exit (1);}
  
  *Tt_ra_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating Tt_ra buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueWriteBuffer(command_queue, *Tt_ra_mem_obj, CL_TRUE, 0, size, HostMem->Tt_ra, 0, NULL, NULL); 
  

  /* Allocate and initialize GPU thread states on the device.
  *
  * We only initialize rnd_a and rnd_x here. For all other fields, whose
  * initial value is a known constant, we use a kernel to do the
  * initialization.
  */

  // photon structure
  size = NUM_THREADS * sizeof(float);
  
  *photon_x_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_x buffer, exiting\n");
    exit(-1);
  }
  
  *photon_y_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_y buffer, exiting\n");
    exit(-1);
  }
  *photon_z_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_z buffer, exiting\n");
    exit(-1);
  }
  *photon_ux_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_ux buffer, exiting\n");
    exit(-1);
  }
  *photon_uy_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_uy buffer, exiting\n");
    exit(-1);
  }
  *photon_uz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_uz buffer, exiting\n");
    exit(-1);
  }
  *photon_w_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_w buffer, exiting\n");
    exit(-1);
  }
  *photon_sleft_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_sleft buffer, exiting\n");
    exit(-1);
  }
  size = NUM_THREADS * sizeof(UINT32);
  *photon_layer_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_layer buffer, exiting\n");
    exit(-1);
  }

  // thread active
  *is_active_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating is_active buffer, exiting\n");
    exit(-1);
  }

  return 1;
}

//////////////////////////////////////////////////////////////////////////////
//   Transfer data from Device to Host memory after simulation
//////////////////////////////////////////////////////////////////////////////
int CopyDeviceToHostMem(SimState* HostMem,SimulationStruct* sim, cl_command_queue command_queue, cl_mem A_rz_mem_obj, cl_mem Rd_ra_mem_obj, cl_mem Tt_ra_mem_obj, cl_mem x_mem_obj, cl_mem scaled_w_mem_obj, cl_mem debug_mem_obj, Tetra *tetra_mesh)
{
  int rz_size = sim->det.nr*sim->det.nz;
  int ra_size = sim->det.nr*sim->det.na;


  cl_int ret;
  ret = clEnqueueReadBuffer(command_queue, scaled_w_mem_obj, CL_TRUE, 0, (sim->nTetras+1)*sizeof(UINT64), HostMem->scaled_w, 0, NULL, NULL);
  if(ret != CL_SUCCESS){
    printf("Error reading scaled weight buffer, exiting\n");
    exit(-1);
  }
  // Copy A_rz, Rd_ra and Tt_ra
  ret = clEnqueueReadBuffer(command_queue, A_rz_mem_obj, CL_TRUE, 0, rz_size*sizeof(UINT64), HostMem->A_rz, 0, NULL, NULL);
  if(ret != CL_SUCCESS){
    printf("Error reading A_rz buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueReadBuffer(command_queue, debug_mem_obj, CL_TRUE, 0, 40*sizeof(float), debug, 0, NULL, NULL);
  if(ret != CL_SUCCESS){
    printf("Error reading debug buffer, exiting\n");
    exit(-1);
  }
  
  for(int i = 0; i < 4; i++)
  {
  printf("dx: %f\n", debug[9*i+0]);
  printf("dy: %f\n", debug[9*i+1]);
  printf("dz: %f\n", debug[9*i+2]);
  printf("ax: %f\n", debug[9*i+3]);
  printf("ay: %f\n", debug[9*i+4]);
  printf("az: %f\n", debug[9*i+5]);
  printf("bx: %f\n", debug[9*i+6]);
  printf("by: %f\n", debug[9*i+7]);
  printf("bz: %f\n", debug[9*i+8]);
  printf("|a|: %f\n", debug[9*i+3]*debug[9*i+3]+debug[9*i+4]*debug[9*i+4]+debug[9*i+5]*debug[9*i+5]);
  printf("|b|: %f\n", debug[9*i+6]*debug[9*i+6]+debug[9*i+7]*debug[9*i+7]+debug[9*i+8]*debug[9*i+8]);
  printf("a dot d: %f\n", debug[9*i+3]*debug[9*i+0]+debug[9*i+4]*debug[9*i+1]+debug[9*i+5]*debug[9*i+2]);
  printf("a dot b: %f\n", debug[9*i+3]*debug[9*i+6]+debug[9*i+4]*debug[9*i+7]+debug[9*i+5]*debug[9*i+8]);
  printf("b dot d: %f\n", debug[9*i+6]*debug[9*i+0]+debug[9*i+7]*debug[9*i+1]+debug[9*i+8]*debug[9*i+2]);
  printf("\n");
  }
  ret = clEnqueueReadBuffer(command_queue, Rd_ra_mem_obj, CL_TRUE, 0, ra_size*sizeof(UINT64), HostMem->Rd_ra, 0, NULL, NULL);
  if(ret != CL_SUCCESS){
    printf("Error reading Rd_ra buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueReadBuffer(command_queue, Tt_ra_mem_obj, CL_TRUE, 0, ra_size*sizeof(UINT64), HostMem->Tt_ra, 0, NULL, NULL);
  if(ret != CL_SUCCESS){
    printf("Error reading Tt_ra buffer, exiting\n");
    exit(-1);
  }

  //Also copy the state of the RNG's
  ret = clEnqueueReadBuffer(command_queue, x_mem_obj, CL_TRUE, 0, NUM_THREADS * sizeof(UINT64), HostMem->x, 0, NULL, NULL);
  if(ret != CL_SUCCESS){
    printf("Error reading x buffer, exiting\n");
    exit(-1);
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//   Free Host Memory
//////////////////////////////////////////////////////////////////////////////
void FreeHostSimState(SimState *hstate)
{
  if (hstate->n_photons_left != NULL)
  {
    free(hstate->n_photons_left); hstate->n_photons_left = NULL;
  }

  // DO NOT FREE RANDOM NUMBER SEEDS HERE.

  if (hstate->A_rz != NULL)
  {
    free(hstate->A_rz); hstate->A_rz = NULL;
  }
  if (hstate->Rd_ra != NULL)
  {
    free(hstate->Rd_ra); hstate->Rd_ra = NULL;
  }
  if (hstate->Tt_ra != NULL)
  {
    free(hstate->Tt_ra); hstate->Tt_ra = NULL;
  }
  if (hstate->scaled_w != NULL)
  {
    free(hstate->scaled_w); hstate->scaled_w = NULL;
  }
}

//////////////////////////////////////////////////////////////////////////////
//   Free GPU Memory
//////////////////////////////////////////////////////////////////////////////
void FreeDeviceSimStates(cl_context context, cl_command_queue command_queue,cl_kernel initkernel, cl_kernel kernel, 
        cl_program program, cl_mem simparam_mem_obj, cl_mem num_photons_left_mem_obj, 
        cl_mem a_mem_obj, cl_mem x_mem_obj, cl_mem A_rz_mem_obj, cl_mem Rd_ra_mem_obj, cl_mem Tt_ra_mem_obj, 
        cl_mem photon_x_mem_obj, cl_mem photon_y_mem_obj,cl_mem photon_z_mem_obj, cl_mem photon_ux_mem_obj, 
        cl_mem photon_uy_mem_obj, cl_mem photon_uz_mem_obj, cl_mem photon_w_mem_obj, cl_mem photon_sleft_mem_obj,
        cl_mem photon_layer_mem_obj, cl_mem is_active_mem_obj, cl_mem tetra_mesh_mem_obj, cl_mem materials_mem_obj,
        cl_mem scaled_w_mem_obj, cl_mem debug_mem_obj
     )
{
 cl_int ret;
 ret = clFlush(command_queue);
 if(ret!= CL_SUCCESS){
    printf("Error flushing queue in FreeDeviceSimStates, exiting\n");
    exit(-1);
 }
 ret = clFinish(command_queue);
 if(ret!= CL_SUCCESS){
    printf("Error finishing queue in FreeDeviceSimStates, exiting\n");
    exit(-1);
 }
// ret = clReleaseKernel(initkernel);
// if(ret!= CL_SUCCESS){
//    printf("Error releasing kernel, exiting\n");
//    exit(-1);
// }
 ret = clReleaseKernel(kernel);
 if(ret!= CL_SUCCESS){
    printf("Error releasing init kernel, exiting\n");
    exit(-1);
 }
 ret = clReleaseProgram(program);
 if(ret!= CL_SUCCESS){
    printf("Error releasing program, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(simparam_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing simparam mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(num_photons_left_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing num_photons_left mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(tetra_mesh_mem_obj);
 if(ret!=CL_SUCCESS){
   printf("Error releasing tetra_mesh mem obj, exiting\n");
   exit(-1);
 }
 ret = clReleaseMemObject(materials_mem_obj);
 if(ret!=CL_SUCCESS){
   printf("Error releasing materials mem obj, exiting\n");
   exit(-1);
 }
 ret = clReleaseMemObject(a_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing a mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(x_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing x mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(scaled_w_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing scaled weight mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(A_rz_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing A_rz mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(debug_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing debug mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(Rd_ra_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing Rd_ra mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(Tt_ra_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing Tt_ra mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_x_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_x mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_y_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_y mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_z_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_z mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_ux_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_ux mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_uy_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_uy mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_uz_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_uz mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_w_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_w mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_sleft_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_sleft mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_layer_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_layer mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(is_active_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing is_active mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseCommandQueue(command_queue);
 if(ret!= CL_SUCCESS){
    printf("Error releasing command_queue mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseContext(context);
 if(ret!= CL_SUCCESS){
    printf("Error releasing context mem obj, exiting\n");
    exit(-1);
 }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

