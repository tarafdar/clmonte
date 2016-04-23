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

#define MAX_SOURCE_SIZE 0x100000
#define STRING_BUFFER_LEN 1024

using namespace std;
using namespace aocl_utils;

bool debug_mode = false;
bool use_cache = false;
bool cache_test = false;

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

time_t get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        printf("Error getting time of day!\n");
        return 0;
    }
    return time.tv_sec*1000000 + time.tv_usec;
}

// Supports 1 FPGA only
int RunFPGA(HostThreadState *hstate, Tetra *tetra_mesh, Material *materialspec, Source *p_src)
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
	/*char char_buffer[STRING_BUFFER_LEN]; 
	printf("Querying platform for info:\n");
	printf("==========================\n");
	clGetPlatformInfo(platform, CL_PLATFORM_NAME, STRING_BUFFER_LEN, char_buffer, NULL);
	printf("%-40s = %s\n", "CL_PLATFORM_NAME", char_buffer);
	clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, STRING_BUFFER_LEN, char_buffer, NULL);
	printf("%-40s = %s\n", "CL_PLATFORM_VENDOR ", char_buffer);
	clGetPlatformInfo(platform, CL_PLATFORM_VERSION, STRING_BUFFER_LEN, char_buffer, NULL);
	printf("%-40s = %s\n\n", "CL_PLATFORM_VERSION ", char_buffer);*/

	// Query the available OpenCL devices.
	scoped_array<cl_device_id> devices;
	cl_uint num_devices;
	devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &num_devices));

	// We'll just use the first device.
	device = devices[0];

	// Display some device information.
	//display_device_info(device);

	// Create the context.
	context = clCreateContext(NULL, 1, &device, &oclContextCallback, NULL, &status);
	checkError(status, "Failed to create context");  
  
	// Create the command queue.
	command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
	checkError(status, "Failed to create command queue");
  
	HostMem->absorption = (UINT32*)alignedMalloc(sizeof(UINT32)*(num_tetras + 1));
	memset(HostMem->absorption, 0, sizeof(UINT32)*(num_tetras + 1));

	HostMem->transmittance = (UINT32*)alignedMalloc(sizeof(UINT32)*num_tetras*4);
	memset(HostMem->transmittance, 0, sizeof(UINT32)*num_tetras*4);


	UINT32 launcher_rng_init;
	UINT32 stepper_rng_init;
	UINT32 dropper_rng_init;
	UINT32 spinner_rng_init;
	UINT32 reflactor_rng_init;
	if(getenv("MONTE_DETERMINISTIC")) {
		launcher_rng_init = 79;
		stepper_rng_init = 80;
		dropper_rng_init = 81;
		spinner_rng_init = 82;
		reflactor_rng_init = 83;
	}
	else {
		srand((UINT32)time(NULL));
		launcher_rng_init = (UINT32)rand();
		stepper_rng_init = (UINT32)rand();
		dropper_rng_init = (UINT32)rand();
		spinner_rng_init = (UINT32)rand();
		reflactor_rng_init = (UINT32)rand();
	}
	
	printf("Using random seeds (%u, %u, %u, %u, %u)\n", launcher_rng_init, stepper_rng_init, dropper_rng_init, spinner_rng_init, reflactor_rng_init);

	DevDataSeg dev_dseg;
	init_device_data_seg(context, command_queue, hstate->sim->number_of_photons, &source, tetra_mesh, materialspec, HostMem->absorption, HostMem->transmittance, &dev_dseg, num_tetras, num_mats, launcher_rng_init, stepper_rng_init, dropper_rng_init, spinner_rng_init, reflactor_rng_init);

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

	set_kernel_args(&kernel, &dev_dseg, num_tetras, num_mats);

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

	// Cleanup
	clReleaseKernel(kernel);
	clReleaseProgram(program);
	dtor_device_data_seg(&dev_dseg);
	clReleaseCommandQueue(command_queue);
	clReleaseContext(context);
	return 0;
}

static void DoOneSimulation(int sim_id, SimulationStruct* simulation, Tetra *tetra_mesh, Material *materialspec,
                            TriNode *trinodes, TetraNode *tetranodes, Source *p_src)
{
	// Start simulation kernel exec timer
	time_t wall0 = get_wall_time();
	// For each GPU, init the host-side structure.
	HostThreadState* hstates;
	hstates = (HostThreadState*)malloc(sizeof(HostThreadState));

	hstates->sim = simulation;

	SimState *hss = &(hstates->host_sim_state);

	// Launch simulation
	RunFPGA (hstates, tetra_mesh, materialspec, p_src);

  // End the timer.
  time_t wall1 = get_wall_time();
  double elapsed = (wall1 - wall0) / 1000000.0;
  cout << "Wall Time = " << elapsed << " seconds\n";
  
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

	// If emulation enable debug mode
	debug_mode = (getenv("MONTE_DEBUG") != NULL);
	cache_test = (getenv("MONTE_CACHE_TEST") != NULL);
	use_cache = (cache_test || (getenv("MONTE_CACHE") != NULL));

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
  
	// Post process add adj material n so algorithm can avoid tetra lookup and directly get 
	// neighbouring material interface index coefficient in the case of reflection
	UINT32 z = 1;
	for(z = 1; z <= Nt;z++) {
		tetra_mesh[z].adjN[0] = materialspec[tetra_mesh[tetra_mesh[z].adjTetras[0]].matID].n;
		tetra_mesh[z].adjN[1] = materialspec[tetra_mesh[tetra_mesh[z].adjTetras[1]].matID].n;
		tetra_mesh[z].adjN[2] = materialspec[tetra_mesh[tetra_mesh[z].adjTetras[2]].matID].n;
		tetra_mesh[z].adjN[3] = materialspec[tetra_mesh[tetra_mesh[z].adjTetras[3]].matID].n;
	}

  unsigned long long seed = 5;//(unsigned long long) time(NULL);

  int i;

  // Read the simulation inputs.
  printf("Number of photons: %d\n", simulation.number_of_photons);

  //perform the simulation
  simulation.nTetras = Nt;
  simulation.nMaterials = Nm;
  // Run a simulation
  DoOneSimulation(i, &simulation, tetra_mesh, materialspec, trinodes,
  tetranodes, sourcepoint);

  alignedFree(materialspec);
  alignedFree(tetra_mesh);
  free(trinodes);
  free(tetranodes);
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

