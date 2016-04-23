#include <float.h> //for FLT_MAX 
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <CL/cl.h>
#include "host/h/gpumcml.h"
#include "host/h/comm.h"
#include "h/misc.h"
#include <iostream>

#include "CL/opencl.h"
#include "AOCLUtils/aocl_utils.h"

using namespace std;
using namespace aocl_utils;
#define MAX_SOURCE_SIZE 0x100000
#define STRING_BUFFER_LEN 1024

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
cl_mem absorption_mem_obj;
cl_mem transmittance_mem_obj;

cl_mem photon_x_mem_obj;
cl_mem photon_y_mem_obj;
cl_mem photon_z_mem_obj;
cl_mem photon_dx_mem_obj;
cl_mem photon_dy_mem_obj;
cl_mem photon_dz_mem_obj;
cl_mem photon_w_mem_obj;
cl_mem photon_tetra_id_mem_obj;
cl_mem photon_mat_id_mem_obj;
cl_mem is_active_mem_obj;

cl_kernel kernel;
cl_program program;

static cl_platform_id platform = NULL;
static cl_device_id device = NULL;

void cleanup();
static void device_info_ulong( cl_device_id device, cl_device_info param, const char* name);
static void device_info_uint( cl_device_id device, cl_device_info param, const char* name);
static void device_info_bool( cl_device_id device, cl_device_info param, const char* name);
static void device_info_string( cl_device_id device, cl_device_info param, const char* name);
static void display_device_info( cl_device_id device );
//////////////////////////////////////////////////////////////////////////////
// time measuring helper function
//////////////////////////////////////////////////////////////////////////////

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        printf("Error getting time of day!\n");
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

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
int RunGPUi(HostThreadState *hstate, Tetra *tetra_mesh, Material *materialspec, Source *p_src)
{
	unsigned int num_tetras = hstate->sim->nTetras;
	unsigned int num_mats = hstate->sim->nMaterials;

	FPGA::Source source;
	source.pos.x = p_src->x;
	source.pos.y = p_src->y;
	source.pos.z = p_src->z;
	source.tid = p_src->IDt;

	cl_int status;
	cl_int ret;
	SimState *HostMem = &(hstate->host_sim_state);

	// Get the OpenCL platform.
	platform = findPlatform("Altera");
	if(platform == NULL) {
		printf("ERROR: Unable to find Altera OpenCL platform.\n");
		return false;
	}

  // User-visible output - Platform information
  {
    char char_buffer[STRING_BUFFER_LEN]; 
    printf("Querying platform for info:\n");
    printf("==========================\n");
    clGetPlatformInfo(platform, CL_PLATFORM_NAME, STRING_BUFFER_LEN, char_buffer, NULL);
    printf("%-40s = %s\n", "CL_PLATFORM_NAME", char_buffer);
    clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, STRING_BUFFER_LEN, char_buffer, NULL);
    printf("%-40s = %s\n", "CL_PLATFORM_VENDOR ", char_buffer);
    clGetPlatformInfo(platform, CL_PLATFORM_VERSION, STRING_BUFFER_LEN, char_buffer, NULL);
    printf("%-40s = %s\n\n", "CL_PLATFORM_VERSION ", char_buffer);
  }
  
	// Query the available OpenCL devices.
	scoped_array<cl_device_id> devices;
	cl_uint num_devices;
	devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &num_devices));

	// We'll just use the first device.
	device = devices[0];

	// Display some device information.
	display_device_info(device);

	// Create the context.
	context = clCreateContext(NULL, 1, &device, &oclContextCallback, NULL, &status);
	checkError(status, "Failed to create context");  
  
	// Create the command queue.
	command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
	checkError(status, "Failed to create command queue");
  
	HostMem->absorption = (UINT32*)malloc(sizeof(UINT32)*(num_tetras + 1));
	memset(HostMem->absorption, 0, sizeof(UINT32)*(num_tetras + 1));

	HostMem->transmittance = (UINT32*)malloc(sizeof(UINT32)*num_tetras*4);
	memset(HostMem->transmittance, 0, sizeof(UINT32)*num_tetras*4);

	DevDataSeg dev_dseg;
	init_device_data_seg(context, command_queue, *HostMem->n_photons_left, &source, tetra_mesh, materialspec, HostMem->absorption, HostMem->transmittance, &dev_dseg, num_tetras, num_mats);

	// Create the program.
	const char* monte_aocx = getenv("MONTE_AOCX");
	if(monte_aocx == NULL) {
		printf("MONTE_AOCX environment variable not set!");
		return 1;
	}
	std::string binary_file = getBoardBinaryFile(monte_aocx, device);
	printf("Using AOCX: %s\n", binary_file.c_str());
	program = createProgramFromBinary(context, binary_file.c_str(), &device, 1);

	// Build the program that was just created.
	status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
	checkError(status, "Failed to build program");

	// Create the kernel - name passed in here must match kernel name in the
	// original CL file, that was compiled into an AOCX file using the AOC tool
	const char *kernel_name = "entry";  // Kernel name, as defined in the CL file
	kernel = clCreateKernel(program, kernel_name, &status);
	checkError(status, "Failed to create kernel");

	set_kernel_args(&kernel, &dev_dseg);

	size_t global_size = 1;
	size_t local_size = 1;
	ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
	if(ret != CL_SUCCESS){
		printf("Error enqueundrange of kernel, exiting\n");
		exit(-1);
	}

	clFlush(command_queue);
	clFinish(command_queue);

	read_kernel_results(command_queue, &dev_dseg, HostMem->absorption, HostMem->transmittance, num_tetras);
	clFlush(command_queue);
	clFinish(command_queue);

	// TODO: cleanup cl_mem objects and kernels/programs etc...
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
//   Perform MCML simulation for one run out of N runs (in the input file)
//////////////////////////////////////////////////////////////////////////////
static void DoOneSimulation(int sim_id, SimulationStruct* simulation,
                            unsigned long long *x, unsigned int *a, Tetra *tetra_mesh, Material *materialspec,
                            TriNode *trinodes, TetraNode *tetranodes, Source *p_src)
{
	// Start simulation kernel exec timer
	double wall0 = get_wall_time();
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
  int number_of_iterations = RunGPUi (hstates, tetra_mesh, materialspec, p_src);

  // End the timer.
  double wall1 = get_wall_time();
  double elapsed = wall1 - wall0;
  float total_steps = (float)number_of_iterations*NUM_STEPS*NUM_THREADS;
  cout << "Wall Time = " << elapsed << " seconds\n";
  printf("total num of iterations = %d\n", number_of_iterations);
  printf("NUM_STEPS = %d\n", NUM_STEPS);
  printf("NUM_THREADS = %d\n", NUM_THREADS);
  printf( ">>>>>>Simulation Speed: %e photon events per second\n", total_steps/elapsed);
  
  Write_Simulation_Results(hss, simulation, elapsed, trinodes, tetranodes, materialspec, tetra_mesh, simulation->outp_filename);

  Conservation_Of_Energy(hss, simulation, trinodes);
/*
  // Free SimState structs.
  FreeHostSimState(hss);
  free(hstates);
*/
}

void banner()
{
    printf("FullMonte OpenCL v0.0\n");
    printf("(c) Christopher Monardo, Robert Matyjewicz 2015-2016\n\n");
}

void OutputTetraMeshHelper(Tetra* tetra_mesh, int currentID, bool* visited)
{
	Tetra& current_tetra = tetra_mesh[currentID];
	visited[currentID] = true;
	printf("Tetra %d: Adjacencies: [%d, %d, %d, %d]\n", currentID, current_tetra.adjTetras[0], current_tetra.adjTetras[1], current_tetra.adjTetras[2], current_tetra.adjTetras[3]);

	for(int i = 0;i < 4;i++) {
		if(visited[current_tetra.adjTetras[i]])
			continue;

		OutputTetraMeshHelper(tetra_mesh, current_tetra.adjTetras[i], visited);
	}
}

void OutputTetraMesh(Tetra *tetra_mesh, int Nt)	//debug code
{
	bool* visited = (bool*)malloc(sizeof(bool)*(Nt+1));
	memset(visited, 0, sizeof(bool)*(Nt + 1));

	OutputTetraMeshHelper(tetra_mesh, 1, visited);
	return;
  int i;
  for(i=48; i<=48; i++)
  {
    Tetra &t = tetra_mesh[i];
    printf("Tetra%d: material: %d\n", i, t.matID);
    printf("  %f %f %f %f %d\n", t.face[0][0], t.face[0][1], t.face[0][2], t.face[0][3], t.adjTetras[0]);
    printf("  %f %f %f %f %d\n", t.face[1][0], t.face[1][1], t.face[1][2], t.face[1][3], t.adjTetras[1]);
    printf("  %f %f %f %f %d\n", t.face[2][0], t.face[2][1], t.face[2][2], t.face[2][3], t.adjTetras[2]);
    printf("  %f %f %f %f %d\n", t.face[3][0], t.face[3][1], t.face[3][2], t.face[3][3], t.adjTetras[3]);
  }
  printf("\n");
}

void OutputMaterial(Material *mats, int Nm)
{
  int i;
  for(i=1; i<=Nm; i++)
  {
    Material &mat = mats[i];
    printf("Material %d:  mua %f  mus %f  muas %f  rmuas %f\n",i, mat.mu_a, mat.mu_s, mat.mu_as, mat.rmu_as);
    //printf("  n %f\n  g %f\n  HGCoeff1 %f\n  HGCoeff2 %f\n  absfrac %f\n", mat.n, mat.g, mat.HGCoeff1, mat.HGCoeff2, mat.absfrac);
  }
  printf("\n");
}

void OutputSource(Source *p_src)
{
  printf("Source\n");
  printf("  initial weight: %d\n", p_src->Np);
  printf("  x: %f\n  y: %f\n  z: %f\n", p_src->x, p_src->y, p_src->z);
  printf("  dx: %f\n  dy: %f\n  dz: %f\n", p_src->dx, p_src->dy, p_src->dz);
  printf("  initial tetra ID: %d\n", p_src->IDt);
  printf("\n");
}

//////////////////////////////////////////////////////////////////////////////
//   Perform MCML simulation for one run out of N runs (in the input file)
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  banner();

  SimulationStruct simulation;
  // Parse command-line arguments.
  if (interpret_arg(argc, argv, &simulation))
  {
    usage(argv[0]);
    return 1;
  }

	const char* monte_aocx = getenv("MONTE_AOCX");
	if(monte_aocx == NULL) {
		printf("MONTE_AOCX environment variable not set!\n");
		return 1;
	}

  Tetra *tetra_mesh;
  Point *points;
  TriNode *trinodes;
  TetraNode *tetranodes;
  int Np, Nt, Nm;	//number of points, number of tetrahedra, number of materials
  char filename[256];
  strncpy(filename, simulation.inp_filename, 250);
  strcat(filename, ".mesh");
  PopulateTetraFromMeshFile(filename, &tetra_mesh, &points, &trinodes, &tetranodes, &Np, &Nt);
  Material *materialspec;
  //OutputTetraMesh(tetra_mesh, Nt);	//debug code

  strncpy(filename, simulation.inp_filename, 248);
  strcat(filename, ".source");
  Source *sourcepoint;
  ParseSource(filename, &sourcepoint, tetra_mesh, Nt, points, tetranodes);
  //OutputSource(sourcepoint);	//debug code

  //Material *mat;
  strncpy(filename, simulation.inp_filename, 251);
  strcat(filename, ".opt");
  ParseMaterial(filename, &materialspec, &Nm);
  //OutputMaterial(materialspec, Nm);	//debug code
  
  unsigned long long seed = 5;//(unsigned long long) time(NULL);

  int i;

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

  //if (init_RNG(x, a, len, "executable/safeprimes_base32.txt", seed)) return 1;
  
  printf("Using the MWC random number generator ...\n");

  //perform the simulation
  simulation.nTetras = Nt;
  simulation.nMaterials = Nm;
  // Run a simulation
  DoOneSimulation(i, &simulation, x, a, tetra_mesh, materialspec, trinodes,
  tetranodes, sourcepoint);

  // Free the random number seed arrays.
  free(materialspec);
  free(tetra_mesh);
  free(trinodes);
  free(tetranodes);
  free(x); free(a);
  free(points);
  return 0; 
}

//////////////////////////////////////////////////////////////////////////////

void cleanup() {
  if(kernel) {
    clReleaseKernel(kernel);  
  }
  if(program) {
    clReleaseProgram(program);
  }
  if(command_queue) {
    clReleaseCommandQueue(command_queue);
  }
  if(context) {
    clReleaseContext(context);
  }
}

// Helper functions to display parameters returned by OpenCL queries
static void device_info_ulong( cl_device_id device, cl_device_info param, const char* name) {
   cl_ulong a;
   clGetDeviceInfo(device, param, sizeof(cl_ulong), &a, NULL);
   printf("%-40s = %lu\n", name, a);
}
static void device_info_uint( cl_device_id device, cl_device_info param, const char* name) {
   cl_uint a;
   clGetDeviceInfo(device, param, sizeof(cl_uint), &a, NULL);
   printf("%-40s = %u\n", name, a);
}
static void device_info_bool( cl_device_id device, cl_device_info param, const char* name) {
   cl_bool a;
   clGetDeviceInfo(device, param, sizeof(cl_bool), &a, NULL);
   printf("%-40s = %s\n", name, (a?"true":"false"));
}
static void device_info_string( cl_device_id device, cl_device_info param, const char* name) {
   char a[STRING_BUFFER_LEN]; 
   clGetDeviceInfo(device, param, STRING_BUFFER_LEN, &a, NULL);
   printf("%-40s = %s\n", name, a);
}

// Query and display OpenCL information on device and runtime environment
static void display_device_info( cl_device_id device ) {

   printf("Querying device for info:\n");
   printf("========================\n");
   device_info_string(device, CL_DEVICE_NAME, "CL_DEVICE_NAME");
   device_info_string(device, CL_DEVICE_VENDOR, "CL_DEVICE_VENDOR");
   device_info_uint(device, CL_DEVICE_VENDOR_ID, "CL_DEVICE_VENDOR_ID");
   device_info_string(device, CL_DEVICE_VERSION, "CL_DEVICE_VERSION");
   device_info_string(device, CL_DRIVER_VERSION, "CL_DRIVER_VERSION");
   device_info_uint(device, CL_DEVICE_ADDRESS_BITS, "CL_DEVICE_ADDRESS_BITS");
   device_info_bool(device, CL_DEVICE_AVAILABLE, "CL_DEVICE_AVAILABLE");
   device_info_bool(device, CL_DEVICE_ENDIAN_LITTLE, "CL_DEVICE_ENDIAN_LITTLE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHE_SIZE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_SIZE, "CL_DEVICE_GLOBAL_MEM_SIZE");
   device_info_bool(device, CL_DEVICE_IMAGE_SUPPORT, "CL_DEVICE_IMAGE_SUPPORT");
   device_info_ulong(device, CL_DEVICE_LOCAL_MEM_SIZE, "CL_DEVICE_LOCAL_MEM_SIZE");
   device_info_ulong(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, "CL_DEVICE_MAX_CLOCK_FREQUENCY");
   device_info_ulong(device, CL_DEVICE_MAX_COMPUTE_UNITS, "CL_DEVICE_MAX_COMPUTE_UNITS");
   device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_ARGS, "CL_DEVICE_MAX_CONSTANT_ARGS");
   device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, "CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE");
   device_info_uint(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS");
   device_info_uint(device, CL_DEVICE_MEM_BASE_ADDR_ALIGN, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS");
   device_info_uint(device, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, "CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE");

   {
      cl_command_queue_properties ccp;
      clGetDeviceInfo(device, CL_DEVICE_QUEUE_PROPERTIES, sizeof(cl_command_queue_properties), &ccp, NULL);
      printf("%-40s = %s\n", "Command queue out of order? ", ((ccp & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE)?"true":"false"));
      printf("%-40s = %s\n", "Command queue profiling enabled? ", ((ccp & CL_QUEUE_PROFILING_ENABLE)?"true":"false"));
   }
}


//////////////////////////////////////////////////////////////////////////////

