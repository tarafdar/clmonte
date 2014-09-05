
#ifndef _GPUMCML_KERNEL_H_
#define _GPUMCML_KERNEL_H_

//#include "gpumcml.h"

#define NUM_STEPS 50000  //Use 5000 for faster response time

#define NUM_BLOCKS 150
#define NUM_THREADS_PER_BLOCK 512

//#define NUM_BLOCKS 30
//#define NUM_THREADS_PER_BLOCK 528

#define NUM_THREADS (NUM_BLOCKS * NUM_THREADS_PER_BLOCK)

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

typedef struct
{
  float init_photon_w;      // initial photon weight 

  float dz;                 // z grid separation.[cm] 
  float dr;                 // r grid separation.[cm] 

  UINT32 na;                // array range 0..na-1. 
  UINT32 nz;                // array range 0..nz-1. 
  UINT32 nr;                // array range 0..nr-1. 

  UINT32 num_layers;        // number of layers. 
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
  // cartesian coordinates of the photon [cm]
  float x;
  float y;
  float z;

  float ax;
  float ay;
  float az;

  float bx;
  float by;
  float bz;

  // directional cosines of the photon
  float ux;
  float uy;
  float uz;

  float w;            // photon weight

  float s;            // step size [cm]
  float sleft;        // leftover step size [cm]

  // index to layer where the photon resides
  UINT32 layer;

  // flag to indicate if photon hits a boundary
  UINT32 hit;
} PhotonStructGPU;

#endif // _GPUMCML_KERNEL_H_

