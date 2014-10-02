
#include <stdio.h>
#include <string>
#include <math.h>
#include <CL/cl.h>
#include "gpumcml.h"
#include "kernel.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//   Initialize Device Constant Memory with read-only data
//////////////////////////////////////////////////////////////////////////////
int InitDCMem(SimulationStruct *sim, cl_context context, cl_command_queue command_queue, cl_mem *simparam_mem_obj, cl_mem *layerspecs_mem_obj )
{
  // Make sure that the number of layers is within the limit.
  UINT32 n_layers = sim->n_layers + 2;
  if (n_layers > MAX_LAYERS) return 1;

  SimParamGPU h_simparam;


  h_simparam.num_layers = sim->n_layers;  // not plus 2 here
  h_simparam.init_photon_w = sim->start_weight;
  h_simparam.dz = sim->det.dz;
  h_simparam.dr = sim->det.dr;
  
  
  h_simparam.na = sim->det.na;
  h_simparam.nz = sim->det.nz;
  h_simparam.nr = sim->det.nr;

  printf("number of simulation layers %d", sim->n_layers);
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

  LayerStructGPU h_layerspecs[MAX_LAYERS];
  float h_layerspecs_floats[MAX_LAYERS*9];
  for (UINT32 i = 0; i < n_layers; ++i)
  {
    h_layerspecs[i].z0 = sim->layers[i].z_min;
    h_layerspecs[i].z1 = sim->layers[i].z_max;
    float n1 = sim->layers[i].n;
    h_layerspecs[i].n = n1;

    // TODO: sim->layer should not do any pre-computation.
    float rmuas = sim->layers[i].mutr;
    h_layerspecs[i].muas = FP_ONE / rmuas;
    h_layerspecs[i].rmuas = rmuas;
    h_layerspecs[i].mua_muas = sim->layers[i].mua * rmuas;

    h_layerspecs[i].g = sim->layers[i].g;

    if (i == 0 || i == n_layers-1)
    {
      h_layerspecs[i].cos_crit0 = MCML_FP_ZERO;
      h_layerspecs[i].cos_crit1 = MCML_FP_ZERO;
    }
    else
    {
      float n2 = sim->layers[i-1].n;
      h_layerspecs[i].cos_crit0 = (n1 > n2) ?
        sqrt(FP_ONE - n2*n2/(n1*n1)) : MCML_FP_ZERO;
      n2 = sim->layers[i+1].n;
      h_layerspecs[i].cos_crit1 = (n1 > n2) ?
        sqrtf(FP_ONE - n2*n2/(n1*n1)) : MCML_FP_ZERO;
    }
  }

  *layerspecs_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(LayerStructGPU)*n_layers, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating layerspecs buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueWriteBuffer(command_queue, *layerspecs_mem_obj, CL_TRUE, 0, sizeof(LayerStructGPU) *n_layers, h_layerspecs, 0, NULL, NULL); 
  // Copy layer data to constant device memory
  if(ret!= CL_SUCCESS){
    printf("Error writing to layerspecs buffer, exiting\n");
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
        cl_mem *photon_layer_mem_obj, cl_mem *is_active_mem_obj
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
    printf("Error creating layerspecs buffer, exiting\n");
    exit(-1);
  }

  ret = clEnqueueWriteBuffer(command_queue, *num_photons_left_mem_obj, CL_TRUE, 0, size, HostMem->n_photons_left, 0, NULL, NULL); 
  if(ret!= CL_SUCCESS){
    printf("Error writing to a mem buffer, exiting\n");
    exit(-1);
  }

  // random number generation (on device only)
  size = NUM_THREADS * sizeof(UINT32);

  
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
  for(int i=0; i<ra_size; i++){
    HostMem->A_rz[i] = 0;
  }
  ret = clEnqueueWriteBuffer(command_queue, *A_rz_mem_obj, CL_TRUE, 0, size, HostMem->A_rz, 0, NULL, NULL); 
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
int CopyDeviceToHostMem(SimState* HostMem,SimulationStruct* sim, cl_command_queue command_queue, cl_mem A_rz_mem_obj, cl_mem Rd_ra_mem_obj, cl_mem Tt_ra_mem_obj, cl_mem x_mem_obj)
{
  int rz_size = sim->det.nr*sim->det.nz;
  int ra_size = sim->det.nr*sim->det.na;

  // Copy A_rz, Rd_ra and Tt_ra
  cl_int ret;
  ret = clEnqueueReadBuffer(command_queue, A_rz_mem_obj, CL_TRUE, 0, rz_size*sizeof(UINT32), HostMem->A_rz, 0, NULL, NULL);
  if(ret != CL_SUCCESS){
    printf("Error reading A_rz buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueReadBuffer(command_queue, Rd_ra_mem_obj, CL_TRUE, 0, ra_size*sizeof(UINT32), HostMem->Rd_ra, 0, NULL, NULL);
  if(ret != CL_SUCCESS){
    printf("Error reading Rd_ra buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueReadBuffer(command_queue, Tt_ra_mem_obj, CL_TRUE, 0, ra_size*sizeof(UINT32), HostMem->Tt_ra, 0, NULL, NULL);
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
}

//////////////////////////////////////////////////////////////////////////////
//   Free GPU Memory
//////////////////////////////////////////////////////////////////////////////
void FreeDeviceSimStates(cl_context context, cl_command_queue command_queue,cl_kernel initkernel, cl_kernel kernel, 
        cl_program program, cl_mem simparam_mem_obj,  cl_mem layerspecs_mem_obj, cl_mem num_photons_left_mem_obj, 
        cl_mem a_mem_obj, cl_mem x_mem_obj, cl_mem A_rz_mem_obj, cl_mem Rd_ra_mem_obj, cl_mem Tt_ra_mem_obj, 
        cl_mem photon_x_mem_obj, cl_mem photon_y_mem_obj,cl_mem photon_z_mem_obj, cl_mem photon_ux_mem_obj, 
        cl_mem photon_uy_mem_obj, cl_mem photon_uz_mem_obj, cl_mem photon_w_mem_obj, cl_mem photon_sleft_mem_obj,
        cl_mem photon_layer_mem_obj, cl_mem is_active_mem_obj
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
 ret = clReleaseMemObject(layerspecs_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing layerspecs mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(num_photons_left_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing num_photons_left mem obj, exiting\n");
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
 ret = clReleaseMemObject(A_rz_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing A_rz mem obj, exiting\n");
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

