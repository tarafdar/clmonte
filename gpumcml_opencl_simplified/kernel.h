
#ifndef _GPUMCML_KERNEL_H_
#define _GPUMCML_KERNEL_H_

//#include "gpumcml.h"

#define NUM_STEPS 500000  //Use 5000 for faster response time

#define NUM_BLOCKS 30
#define NUM_THREADS_PER_BLOCK 512

//#define NUM_BLOCKS 30
//#define NUM_THREADS_PER_BLOCK 528

#define NUM_THREADS (NUM_BLOCKS * NUM_THREADS_PER_BLOCK)

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

typedef unsigned long long UINT64;
typedef unsigned int UINT32;

typedef struct
{
  float dz;                 // z grid separation.[cm] 
  float init_photon_w;      // initial photon weight
  float dr;                 // r grid separation.[cm] 

  UINT32 na;                // array range 0..na-1. 
  UINT32 nz;                // array range 0..nz-1. 
  UINT32 nr;                // array range 0..nr-1. 

  UINT32 num_layers;        // number of layers.

  float originX;
  float originY;
  float originZ;	//coordinate of the source emitter location
  UINT32 init_tetraID; 
  float terminationThresh;	//weight threshold to start roulette termination
  float proulettewin;	//probability to win the roulette
} SimParamGPU;

typedef struct 
{
  float z0, z1;             // z coordinates of a layer. [cm] 
  float n;                  // refractive index of a layer. 

  float muas;               // mua + mus 
  float rmuas;              // 1/(mua+mus) 
  float mua_muas;           // mua/(mua+mus)

  float g;                  // anisotropy.

  float cos_crit0, cos_crit1;
} LayerStructGPU;

// The max number of layers supported (MAX_LAYERS including 2 ambient layers)
#define MAX_LAYERS 100

//__constant SimParamGPU d_simparam;
//__constant LayerStructGPU d_layerspecs[MAX_LAYERS];

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// Thread-private states that live across batches of kernel invocations
// Each field is an array of length NUM_THREADS.
//
// We use a struct of arrays as opposed to an array of structs to enable
// global memory coalescing.
//
typedef struct
{
  // cartesian coordinates of the photon [cm]
  float *photon_x;
  float *photon_y;
  float *photon_z;

  // directional cosines of the photon
  float *photon_ux;
  float *photon_uy;
  float *photon_uz;

  float *photon_w;            // photon weight
  float *photon_sleft;        // leftover step size [cm]

  // index to layer where the photon resides
  UINT32 *photon_layer;

  UINT32 *is_active;          // is this thread active?
} GPUThreadStates;


typedef struct
{
  // cartesian coordinates of the packet [cm]
  float x;
  float y;
  float z;

  // azimuthal auxiliary vectors a and b
  float ax;
  float ay;
  float az;
  
  float bx;
  float by;
  float bz;

  // direction of the packet
  float dx;
  float dy;
  float dz;

  float w;            // packet weight

  float s;            // step size [cm] to be consumed in next Hop
  float sleft;        // leftover step size [cm]

  // index to layer where the photon resides
  UINT32 layer;

  // flag to indicate if packet hits a boundary
  UINT32 hit;
  // id of the tetrahedron where the packet is currently in
  UINT32 tetraID;
  // index (0,1,2,3) of the current tetra's face to hit
  UINT32 faceIndexToHit;
} Packet;

typedef struct
{
  // The face i's plane is defined by equation (face[i][0]) * x + (face[i][1]) * y + (face[i][2]) * z = face[i][3]
  // The face normal vectors (face[i][0],face[i][1],face[i][2]) always point into the tetrahedron. They are always unit vector.
  float face[4][4];
  
  unsigned adjTetras[4];
  unsigned matID;    
} Tetra ;

typedef struct
{
  float mu_as;	//result of mu_a + mu_s
  float rmu_as;	//reciprocal of mu_as, store this to get rid of slow division arithmetic 
  float n;	//index of refraction
  float g;	//anisotropy constant
  float HGCoeff1;	// HGCoeff1 = (1+g^2)/(2g)
  float HGCoeff2;	// HGCoeff2 = (1-g^2)^2 /(2g). So cos(theta) = HGCoeff1 - HGCoeff2 / (1-g * rand(-1,1))
  float absfrac;	//absorb fraction = 1- albedo = 1 - mus / (mus+mua)
} Material;

#endif // _GPUMCML_KERNEL_H_

