
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
int InitDCMem(SimulationStruct *sim, Source *p_src, Tetra *tetra_mesh, Material *materialspec, cl_context context, cl_command_queue command_queue, cl_mem *simparam_mem_obj, cl_mem *tetra_mesh_mem_obj, cl_mem *materials_mem_obj)
{
  SimParamGPU h_simparam;

  h_simparam.originX = p_src->x;
  h_simparam.originY = p_src->y;
  h_simparam.originZ = p_src->z;
  h_simparam.init_tetraID = p_src->IDt;

  cl_int ret;
  *simparam_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(SimParamGPU), NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating simparam buffer, exiting\n");
    exit(-1);
  }
  ret = clEnqueueWriteBuffer(command_queue, *simparam_mem_obj, CL_TRUE, 0, sizeof(SimParamGPU), &h_simparam, 0, NULL, NULL); 

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
        cl_mem *num_photons_left_mem_obj, cl_mem *a_mem_obj, cl_mem *x_mem_obj, cl_mem *absorption_mem_obj, cl_mem *transmittance_mem_obj, cl_mem *debug_mem_obj, cl_mem *photon_x_mem_obj ,cl_mem *photon_y_mem_obj, cl_mem *photon_z_mem_obj, cl_mem *photon_dx_mem_obj, 
        cl_mem *photon_dy_mem_obj, cl_mem *photon_dz_mem_obj, cl_mem *photon_w_mem_obj,
        cl_mem *photon_tetra_id_mem_obj, cl_mem *photon_mat_id_mem_obj, cl_mem *is_active_mem_obj)
{  
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

  // output arrays absorption and transmittance
  size = (sim->nTetras+1) * sizeof(UINT64);
  *absorption_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!=CL_SUCCESS){
    printf("Error creating absorption mem buffer, exiting\n");
    exit(-1);
  }
  HostMem->absorption = (UINT64*)malloc(size);
  for(int i = 1; i < sim->nTetras+1; i++)
  {
    HostMem->absorption[i] = 0;
  }
  ret = clEnqueueWriteBuffer(command_queue, *absorption_mem_obj, CL_TRUE, 0, size, HostMem->absorption, 0, NULL, NULL);
  if(ret!=CL_SUCCESS){
    printf("Error writing to absorption mem buffer, exiting\n");
    exit(-1);
  }

  size = (sim->nTetras)* 4 * sizeof(UINT64);
  *transmittance_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!=CL_SUCCESS){
    printf("Error creating transmittance mem buffer, exiting\n");
    exit(-1);
  }
  HostMem->transmittance = (UINT64*)malloc(size);
  for(int i = 0; i < (sim->nTetras)* 4; i++)
  {
    HostMem->transmittance[i] = 0;
  }
  ret = clEnqueueWriteBuffer(command_queue, *transmittance_mem_obj, CL_TRUE, 0, size, HostMem->transmittance, 0, NULL, NULL);
  if(ret!=CL_SUCCESS){
    printf("Error writing to transmittance mem buffer, exiting\n");
    exit(-1);
  }

  // Allocate debug on host and device
  size = sizeof(float)*80;
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
  for(int i = 0; i < 80; i++)
  {
    debug[i] = 0;
  }
  ret = clEnqueueWriteBuffer(command_queue, *debug_mem_obj, CL_TRUE, 0, size, debug, 0, NULL, NULL);   
  
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
  *photon_dx_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_dx buffer, exiting\n");
    exit(-1);
  }
  *photon_dy_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_dy buffer, exiting\n");
    exit(-1);
  }
  *photon_dz_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_dz buffer, exiting\n");
    exit(-1);
  }
  *photon_w_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating photon_w buffer, exiting\n");
    exit(-1);
  }
  size = NUM_THREADS * sizeof(UINT32);
  *photon_tetra_id_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating tetra_id buffer, exiting\n");
    exit(-1);
  }

  *photon_mat_id_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if(ret!= CL_SUCCESS){
    printf("Error creating mat_id buffer, exiting\n");
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
int CopyDeviceToHostMem(SimState* HostMem, SimulationStruct* sim, cl_command_queue command_queue, cl_mem x_mem_obj, cl_mem absorption_mem_obj, cl_mem transmittance_mem_obj, cl_mem debug_mem_obj)
{
  cl_int ret;
  ret = clEnqueueReadBuffer(command_queue, absorption_mem_obj, CL_TRUE, 0, (sim->nTetras+1)*sizeof(UINT64), HostMem->absorption, 0, NULL, NULL);
  if(ret != CL_SUCCESS){
    printf("Error reading absorption buffer, exiting\n");
    exit(-1);
  }

  ret = clEnqueueReadBuffer(command_queue, transmittance_mem_obj, CL_TRUE, 0, (sim->nTetras)*4*sizeof(UINT64), HostMem->transmittance, 0, NULL, NULL);
  if(ret != CL_SUCCESS){
    printf("Error reading transmittance buffer, exiting\n");
    exit(-1);
  }

  ret = clEnqueueReadBuffer(command_queue, debug_mem_obj, CL_TRUE, 0, 80*sizeof(float), debug, 0, NULL, NULL);
  if(ret != CL_SUCCESS){
    printf("Error reading debug buffer, exiting\n");
    exit(-1);
  }
  printf("bad d: %d\n", (int)debug[0]);
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
  if (hstate->absorption != NULL)
  {
    free(hstate->absorption); hstate->absorption = NULL;
  }

  if (hstate->transmittance != NULL)
  {
    free(hstate->transmittance); hstate->transmittance = NULL;
  }
}

//////////////////////////////////////////////////////////////////////////////
//   Free GPU Memory
//////////////////////////////////////////////////////////////////////////////
void FreeDeviceSimStates(cl_context context, cl_command_queue command_queue, cl_kernel initkernel, cl_kernel kernel, 
        cl_program program, cl_mem simparam_mem_obj, cl_mem num_photons_left_mem_obj, 
        cl_mem a_mem_obj, cl_mem x_mem_obj, cl_mem tetra_mesh_mem_obj, cl_mem materials_mem_obj,
        cl_mem absorption_mem_obj, cl_mem transmittance_mem_obj, cl_mem debug_mem_obj,
        cl_mem photon_x_mem_obj, cl_mem photon_y_mem_obj,cl_mem photon_z_mem_obj, cl_mem photon_dx_mem_obj, 
        cl_mem photon_dy_mem_obj, cl_mem photon_dz_mem_obj, cl_mem photon_w_mem_obj,
        cl_mem photon_tetra_id_mem_obj, cl_mem photon_mat_id_mem_obj, cl_mem is_active_mem_obj
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
 ret = clReleaseMemObject(absorption_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing absorption mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(transmittance_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing transmittance mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(debug_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing debug mem obj, exiting\n");
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
 ret = clReleaseMemObject(photon_dx_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_dx mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_dy_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_dy mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_dz_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_dz mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_w_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing photon_w mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_tetra_id_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing tetra_id mem obj, exiting\n");
    exit(-1);
 }
 ret = clReleaseMemObject(photon_mat_id_mem_obj);
 if(ret!= CL_SUCCESS){
    printf("Error releasing mat_id mem obj, exiting\n");
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

