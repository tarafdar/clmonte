
#ifndef COMM_H
#define COMM_H

#include <CL/cl.h>
#include "host/h/gpumcml.h"
#include "h/misc.h"

typedef struct DevDataSeg
{
	UINT32 num_photons;
	//UINT32 num_tetras;
	//UINT32 num_mats;
	cl_mem source;
	cl_mem materials;
	cl_mem mesh;
	cl_mem absorption;
	//cl_mem transmittance;
	UINT32 launcher_rng_init;
	UINT32 stepper_rng_init;
	UINT32 dropper_rng_init;
	UINT32 spinner_rng_init;
	UINT32 reflactor_rng_init;
	
	// Used for testing absorbance cache
	cl_mem absorption_cacheless;
} DevDataSeg;

void read_from_device(cl_command_queue command_queue, size_t payload_size, void* payload, cl_mem* mem_buf);
void write_to_device(cl_context context, cl_command_queue command_queue, cl_mem_flags flags, size_t payload_size, void* payload, cl_mem* mem_buf);

void init_device_data_seg(cl_context context, cl_command_queue command_queue, UINT32 num_photons, FPGA::Source* source, Tetra* mesh, Material* mats, UINT32* absorption, UINT32* transmittance, DevDataSeg* dev_dseg, UINT32 num_tetras, UINT32 num_mats, UINT32 launcher_rng_init, UINT32 stepper_rng_init, UINT32 dropper_rng_init, UINT32 spinner_rng_init, UINT32 reflactor_rng_init);
void dtor_device_data_seg(DevDataSeg* dev_dseg);
void set_kernel_args(cl_kernel* kernel, DevDataSeg* dev_dseg, UINT32 num_tetras, UINT32 num_mats);
void read_kernel_results(cl_command_queue command_queue, DevDataSeg* dev_dseg, UINT32* absorption, UINT32* transmittance, UINT32 num_tetras);

#endif // COMM_H





