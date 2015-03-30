/*****************************************************************************
*
*   Header file for common data structures and constants (CPU and GPU)
*
****************************************************************************/
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


#ifndef _GPUMCML_H_
#define _GPUMCML_H_

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



// Critical weight for roulette
#define WEIGHT 1E-4F        

// scaling factor for photon weight, which is then converted to integer
#define WEIGHT_SCALE 32768//12000000

#define PI_const 3.1415926F
#define RPI 0.318309886F

//NOTE: Single Precision
#define COSNINETYDEG 1.0E-6F
#define COSZERO (1.0F - 1.0E-6F)   
#define CHANCE 0.1F

#define MCML_FP_ZERO 0.0F
#define FP_ONE  1.0F
#define FP_TWO  2.0F

#define STR_LEN 200


#define NUM_STEPS 500000  //Use 5000 for faster response time

#define NUM_BLOCKS 30
#define NUM_THREADS_PER_BLOCK 512
#define NUM_THREADS (NUM_BLOCKS * NUM_THREADS_PER_BLOCK)

typedef unsigned long long UINT64;
typedef unsigned int UINT32;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// Simulation input parameters 
typedef struct
{
	char outp_filename[STR_LEN];
	char inp_filename[STR_LEN];

	UINT32 number_of_photons;

	int nTetras;	//total number of tetrahedra
	int nMaterials;	//total number of materials
} SimulationStruct;

// Per-GPU simulation states
// One instance of this struct exists in the host memory, while the other
// in the global memory.
typedef struct
{
	// points to a scalar that stores the number of photons that are not
	// completed (i.e. either on the fly or not yet started)
	int *n_photons_left;

	// per-thread seeds for random number generation
	// arrays of length NUM_THREADS
	// We put these arrays here as opposed to in GPUThreadStates because
	// they live across different simulation runs and must be copied back
	// to the host.
	UINT64 *x;
	UINT32 *a;

	// output data
	UINT64* absorption; // Array of nTetras + 1
	UINT64* transmittance; // Array of nTetras * 4
} SimState;

// Everything a host thread needs to know in order to run simulation on
// one GPU (host-side only)
typedef struct
{
	// those states that will be updated
	SimState host_sim_state;

	// simulation input parameters
	SimulationStruct *sim;

} HostThreadState;

// Used in io.cpp PopulateTetraFromMeshFile
class TwoPointIDsToTetraID
{
public:
	int lowerPointID;
	int higherPointID;
	int TetraID;
	int otherPointID;
	//Overload operator in order to enable list.remove
	bool operator==(const TwoPointIDsToTetraID &value)
	{
		return (this->lowerPointID == value.lowerPointID)
			&& (this->higherPointID == value.higherPointID)
			&& (this->TetraID == value.TetraID)
			&& (this->otherPointID == value.otherPointID);
	}
};

typedef struct
{
	float x;
	float y;
	float z;
} Point;

class TriNode
{
public:
	//Convention here is N0<N1<N2 must be guaranteed.
	int N0;
	int N1;
	int N2;
	float area;
	float fluence;
	bool operator<(const TriNode &value)
	{
		return (this->N0<value.N0)
			|| (this->N0 == value.N0 && this->N1<value.N1)
			|| (this->N0 == value.N0 && this->N1 == value.N1 && this->N2<value.N2);
	}
	int TetraID;
};

typedef struct
{
	int N0;
	int N1;
	int N2;
	int N3;
	float volume;
	float allSurfaceFluence;
} TetraNode;

typedef struct
{
	// The face i's plane is defined by equation (face[i][0]) * x + (face[i][1]) * y + (face[i][2]) * z = face[i][3]
	// The face normal vectors (face[i][0],face[i][1],face[i][2]) always point into the tetrahedron. They are always unit vector.
	float face[4][4];

	UINT32 adjTetras[4];
	UINT32 matID;
} Tetra;

typedef struct
{
	//float mu_a;	//absorption coefficient
	//float mu_s;	//scatter coefficient
	float mu_as;	//attenuation coefficient: result of mu_a + mu_s
	float rmu_as;	//reciprocal of mu_as, store this to get rid of slow division arithmetic 
	float n;	//index of refraction
	float g;	//anisotropy constant
	float HGCoeff1;	// HGCoeff1 = (1+g^2)/(2g)
	float HGCoeff2;	// HGCoeff2 = (1-g^2)^2 /(2g). So cos(theta) = HGCoeff1 - HGCoeff2 / (1-g * rand(-1,1))
	float absfrac;	//absorb fraction = 1- albedo = 1 - mus / (mus+mua)
} Material;

typedef struct
{

	float originX;
	float originY;
	float originZ;	//coordinate of the source emitter location

	float UX, UY, UZ;     //initial direction for pencil source

	int stype;            //source type

	UINT32 init_tetraID;
} SimParamGPU;

typedef struct {
	int stype;
	int Np;	//means weight (in terms of number of photons in this packet)

	// Source Position -> IFF stype = 1 or 11
	float x;
	float y;
	float z;

	// Source Direction -> IFF stype = 11
	float dx;
	float dy;
	float dz;

	// Source Location -> IFF stype = 2
	int IDt;
} Source;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


extern void usage(const char *prog_name);

// Parse the command-line arguments.
// Return 0 if successfull or a +ive error code.
extern int interpret_arg(int argc, char* argv[], SimulationStruct *p_simulation, int *global_size, int *local_size);

extern int Write_Simulation_Results(SimState* HostMem, SimulationStruct* sim, double simulation_time, TriNode* TriNodeList, TetraNode* TetraNodeList, Material* material_spec, Tetra* tetra_mesh, char* output_filename);

int Conservation_Of_Energy(SimState* HostMem, SimulationStruct* sim, TriNode* TriNodeList);

int InitDCMem(SimulationStruct *sim, Source *p_src, Tetra *tetra_mesh, Material *materialspec, cl_context context, cl_command_queue command_queue, cl_mem *simparam_mem_obj, cl_mem *tetra_mesh_mem_obj, cl_mem *materials_mem_obj);

int InitSimStates(SimState* HostMem, SimulationStruct* sim, cl_context context, cl_command_queue command_queue,
	cl_mem *num_photons_simulated_mem_obj, cl_mem *a_mem_obj, cl_mem *x_mem_obj, cl_mem *absorption_mem_obj, cl_mem *transmittance_mem_obj,
	cl_mem *debug_mem_obj, cl_mem *photon_x_mem_obj, cl_mem *photon_y_mem_obj, cl_mem *photon_z_mem_obj, cl_mem *photon_dx_mem_obj,
	cl_mem *photon_dy_mem_obj, cl_mem *photon_dz_mem_obj, cl_mem *photon_w_mem_obj,
	cl_mem *photon_tetra_id_mem_obj, cl_mem *is_active_mem_obj);

int CopyDeviceToHostMem(SimState* HostMem, SimulationStruct* sim, cl_command_queue command_queue, cl_mem x_mem_obj, cl_mem absorption_mem_obj, cl_mem transmittance_mem_obj, cl_mem debug_mem_obj, UINT32 *absorption_buf, UINT32 *transmittance_buf);

void FreeHostSimState(SimState *hstate);

void FreeDeviceSimStates(cl_context context, cl_command_queue command_queue, cl_kernel initkernel, cl_kernel kernel,
	cl_program program, cl_mem simparam_mem_obj, cl_mem num_photons_simulated_mem_obj,
	cl_mem a_mem_obj, cl_mem x_mem_obj, cl_mem tetra_mesh_mem_obj, cl_mem materials_mem_obj,
	cl_mem absorption_mem_obj, cl_mem transmittance_mem_obj, cl_mem debug_mem_obj,
	cl_mem photon_x_mem_obj, cl_mem photon_y_mem_obj, cl_mem photon_z_mem_obj, cl_mem photon_dx_mem_obj,
	cl_mem photon_dy_mem_obj, cl_mem photon_dz_mem_obj, cl_mem photon_w_mem_obj,
	cl_mem photon_tetra_id_mem_obj, cl_mem is_active_mem_obj
	);
void PopulateTetraFromMeshFile(const char* filename, Tetra **p_tetra_mesh, Point **p_points, TriNode **p_trinodes, TetraNode **p_tetranodes, int *p_Np, int *p_Nt);
void PopulateMaterialFromInput(const char*, Material **p_material_spec, int *Nm);
void ParseSource(const char* filename, Source** sourcepoint, const Tetra *tetra_mesh, const int Nt, const Point *points, const TetraNode *tetranodes);
void ParseMaterial(const char *fileName, Material **mat, int *p_Nm);
#endif  // _GPUMCML_H_
