#include <stdio.h>
#include <math.h>
#include <CL/cl.h>
#include "host/h/gpumcml.h"
#include "host/h/comm.h"
#include <assert.h>

void read_from_device(cl_command_queue command_queue, size_t payload_size, void* payload, cl_mem* mem_buf)
{
	cl_int ret;
	ret = clEnqueueReadBuffer(command_queue, *mem_buf, CL_TRUE, 0, payload_size, payload, 0, NULL, NULL);
	assert(ret == CL_SUCCESS);
}

void write_to_device(cl_context context, cl_command_queue command_queue, cl_mem_flags flags, size_t payload_size, void* payload, cl_mem* mem_buf)
{
	cl_int ret;
	*mem_buf = clCreateBuffer(context, flags, payload_size, NULL, &ret);
	assert(ret == CL_SUCCESS);

	ret = clEnqueueWriteBuffer(command_queue, *mem_buf, CL_TRUE, 0, payload_size, payload, 0, NULL, NULL);
	assert(ret == CL_SUCCESS);
}

void set_kernel_arg(cl_kernel* kernel, unsigned int argnum, size_t payload_size, void* payload)
{
	printf("Setting kern arg %u\n", argnum);
	cl_int ret = clSetKernelArg(*kernel, argnum, payload_size, payload);
	switch (ret)
	{
	case CL_INVALID_KERNEL:
		printf("CL_INVALID_KERNEL\n");
		break;
	case CL_INVALID_ARG_INDEX:
		printf("CL_INVALID_ARG_INDEX\n");
		break;
	case CL_INVALID_ARG_VALUE:
		printf("CL_INVALID_ARG_VALUE\n");
		break;
	case CL_INVALID_MEM_OBJECT:
		printf("CL_INVALID_MEM_OBJECT\n");
		break;
	case CL_INVALID_SAMPLER:
		printf("CL_INVALID_SAMPLER\n");
		break;
	case CL_INVALID_ARG_SIZE:
		printf("CL_INVALID_ARG_SIZE\n");
		break; 
	default:
		break;
	}
	
	assert(ret == CL_SUCCESS);
}

void init_device_data_seg(cl_context context, cl_command_queue command_queue, UINT32 num_photons, FPGA::Source* source, Tetra* mesh, Material* mats, UINT32* absorption, UINT32* transmittance, DevDataSeg* dev_dseg, unsigned int num_tetras, unsigned int num_mats)
{
	write_to_device(context, command_queue, CL_MEM_READ_ONLY, sizeof(unsigned int), &num_tetras, &(dev_dseg->num_tetras));
	write_to_device(context, command_queue, CL_MEM_READ_ONLY, sizeof(unsigned int), &num_mats, &(dev_dseg->num_mats));
	write_to_device(context, command_queue, CL_MEM_READ_ONLY, sizeof(Source), source, &(dev_dseg->source));
	write_to_device(context, command_queue, CL_MEM_READ_WRITE, sizeof(UINT32)*(num_tetras + 1), absorption, &(dev_dseg->absorption));
	write_to_device(context, command_queue, CL_MEM_READ_WRITE, sizeof(UINT32)*num_tetras*4, transmittance, &(dev_dseg->transmittance));
	write_to_device(context, command_queue, CL_MEM_READ_ONLY, sizeof(Material)*num_mats, mats, &(dev_dseg->materials));
	write_to_device(context, command_queue, CL_MEM_READ_ONLY, sizeof(Tetra)*num_tetras, mesh, &(dev_dseg->mesh));
	write_to_device(context, command_queue, CL_MEM_READ_ONLY, sizeof(UINT32), &num_photons, &(dev_dseg->num_photons));
}

void set_kernel_args(cl_kernel* kernel, DevDataSeg* dev_dseg)
{
	unsigned int argnum = 0;
	set_kernel_arg(kernel, argnum++, sizeof(cl_mem), &dev_dseg->num_photons);
	set_kernel_arg(kernel, argnum++, sizeof(cl_mem), &dev_dseg->source);
	set_kernel_arg(kernel, argnum++, sizeof(cl_mem), &dev_dseg->materials);
	set_kernel_arg(kernel, argnum++, sizeof(cl_mem), &dev_dseg->mesh);
	set_kernel_arg(kernel, argnum++, sizeof(cl_mem), &dev_dseg->num_tetras);
	set_kernel_arg(kernel, argnum++, sizeof(cl_mem), &dev_dseg->num_mats);
	set_kernel_arg(kernel, argnum++, sizeof(cl_mem), &dev_dseg->absorption);
	set_kernel_arg(kernel, argnum++, sizeof(cl_mem), &dev_dseg->transmittance);
}

void read_kernel_results(cl_command_queue command_queue, DevDataSeg* dev_dseg, UINT32* absorption, UINT32* transmittance, UINT32 num_tetras)
{
	read_from_device(command_queue, sizeof(UINT32)*(num_tetras + 1), absorption, &dev_dseg->absorption);
	read_from_device(command_queue, sizeof(UINT32)*num_tetras*4, transmittance, &dev_dseg->transmittance);
}

