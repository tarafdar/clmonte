typedef ulong UINT64;
typedef uint UINT32;
// Critical weight for roulette
#define WEIGHT 1E-4F        

// scaling factor for photon weight, which is then converted to integer
#define WEIGHT_SCALE 12000000

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


#include "kernel.h"
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

float dot_product (float x1, float y1, float z1, float x2, float y2, float z2) {
    return  x1*x2 + y1*y2 + z1*z2;
}    

//////////////////////////////////////////////////////////////////////////////
//   Generates a random number between 0 and 1 [0,1) 
//////////////////////////////////////////////////////////////////////////////
float rand_MWC_co(UINT64* x,UINT32* a)
{
  *x=(*x&0xfffffffful)*(*a)+(*x>>32);
  return native_divide(convert_float_rtz((UINT32)(*x)),(float)0x100000000);
  
  //*x=((*x)&0xfffffffful)*(*a)+((*x)>>32);
  //return((float)((unsigned int)((*x)&0xfffffffful))/(UINT_MAX));
  
  // The typecast will truncate the x so that it is 0<=x<(2^32-1),
  // __uint2FLOAT_rz ensures a round towards zero since 32-bit FLOATing point 
  // cannot represent all integers that large. 
  // Dividing by 2^32 will hence yield [0,1)
} 

//////////////////////////////////////////////////////////////////////////////
//   Generates a random number between 0 and 1 (0,1]
//////////////////////////////////////////////////////////////////////////////
float rand_MWC_oc(UINT64* x,UINT32* a)
{
  return 1.0f-rand_MWC_co(x,a);
} 

//////////////////////////////////////////////////////////////////////////////
//   AtomicAdd for Unsigned Long Long (ULL) data type
//   Note: 64-bit atomicAdd to global memory are only supported 
//   in graphics card with compute capability 1.2 or above
//////////////////////////////////////////////////////////////////////////////
//void AtomicAddULL(UINT32* address, UINT32 add)
//{
//  if (atomic_add((UINT32*)address,add) +add < add)
//  {
//    atomic_add(((UINT32*)address)+1, 1U);
//  }
//  //atomic_add(address, (int)add);
//}

//////////////////////////////////////////////////////////////////////////////
//   Initialize photon position (x, y, z), direction (ux, uy, uz), weight (w), 
//   step size remainder (sleft), and current layer (layer) 
//   Note: Infinitely narrow beam (pointing in the +z direction = downwards)
//////////////////////////////////////////////////////////////////////////////
void LaunchPhoton(PhotonStructGPU *photon, SimParamGPU d_simparam)
{
  photon->x = MCML_FP_ZERO;
  photon->y = MCML_FP_ZERO;
  photon->z = MCML_FP_ZERO;
  photon->ux = MCML_FP_ZERO; 
  photon->uy = MCML_FP_ZERO;
  photon->uz = FP_ONE;
  photon->w = d_simparam.init_photon_w;
  photon->sleft = MCML_FP_ZERO;
  photon->layer = 1;

  photon->ax = MCML_FP_ZERO;
  photon->az = MCML_FP_ZERO;
  photon->ay = FP_ONE;


  photon->by = MCML_FP_ZERO;
  photon->bz = MCML_FP_ZERO;
  photon->bx = FP_ONE;
}

//////////////////////////////////////////////////////////////////////////////
//   Compute the step size for a photon packet when it is in tissue
//   If sleft is 0, calculate new step size: -log(rnd)/(mua+mus).
//   Otherwise, pick up the leftover in sleft.
//////////////////////////////////////////////////////////////////////////////
void ComputeStepSize(PhotonStructGPU *photon,
                                UINT64 *rnd_x, UINT32 *rnd_a, __global const LayerStructGPU *d_layerspecs)
{
  // Make a new step if no leftover.
  if (photon->sleft == MCML_FP_ZERO)
  {
    //TEMP CHANGE: CHANGE BACK!!! CHANGED BACK :)
    float rand = rand_MWC_oc(rnd_x, rnd_a);
    //float rand = 0;
    photon->s = -log(rand) * d_layerspecs[photon->layer].rmuas;
  }
  else {
    photon->s = photon->sleft * d_layerspecs[photon->layer].rmuas;
    photon->sleft = MCML_FP_ZERO;
  }
}

//////////////////////////////////////////////////////////////////////////////
//   Check if the step size calculated above will cause the photon to hit the 
//   boundary between 2 layers.
//   Return 1 for a hit, 0 otherwise.
//   If the projected step hits the boundary, the photon steps to the boundary
//   and the remainder of the step size is stored in sleft for the next iteration
//////////////////////////////////////////////////////////////////////////////
int HitBoundary(PhotonStructGPU *photon, __global const LayerStructGPU *d_layerspecs)
{
  /* step size to boundary. */
  float dl_b; 

  /* Distance to the boundary. */
  float z_bound = (photon->uz > MCML_FP_ZERO) ?
    d_layerspecs[photon->layer].z1 : d_layerspecs[photon->layer].z0;
  dl_b = native_divide(z_bound - photon->z, photon->uz);     // dl_b > 0

  UINT32 hit_boundary = (photon->uz != MCML_FP_ZERO) && (photon->s > dl_b);
  if (hit_boundary)
  {
    // No need to multiply by (mua + mus), as it is later
    // divided by (mua + mus) anyways (in the original version).
    photon->sleft = (photon->s - dl_b) * d_layerspecs[photon->layer].muas;
    photon->s = dl_b;
  }

  return hit_boundary;
}

//////////////////////////////////////////////////////////////////////////////
//   Move the photon by step size (s) along direction (ux,uy,uz) 
//////////////////////////////////////////////////////////////////////////////
void Hop(PhotonStructGPU *photon)
{
  photon->x += photon->s * photon->ux;
  photon->y += photon->s * photon->uy;
  photon->z += photon->s * photon->uz;
}

//////////////////////////////////////////////////////////////////////////////
//   Drop a part of the weight of the photon to simulate absorption
//////////////////////////////////////////////////////////////////////////////
void Drop(PhotonStructGPU *photon, __global UINT64 *g_A_rz, __global const LayerStructGPU *d_layerspecs, SimParamGPU d_simparam)
{
   
  float dwa = photon->w * d_layerspecs[photon->layer].mua_muas;
  photon->w -= dwa;

  // If scoring of absorption is specified (no -A flag in command line)
    UINT32 iz = native_divide(photon->z, d_simparam.dz);
    UINT32 ir = native_divide(sqrt(photon->x * photon->x + photon->y * photon->y),d_simparam.dr);

    // Only record if photon is not at the edge!!
    // This will be ignored anyways.
    if (iz < d_simparam.nz && ir < d_simparam.nr)
    {
      UINT32 addr = ir * d_simparam.nz + iz;

      // Write to the global memory.
      atomic_add(&g_A_rz[addr], (UINT64)(dwa * WEIGHT_SCALE));
    }

}

//////////////////////////////////////////////////////////////////////////////
//   UltraFast version (featuring reduced divergence compared to CPU-MCML)
//   If a photon hits a boundary, determine whether the photon is transmitted
//   into the next layer or reflected back by computing the internal reflectance
//////////////////////////////////////////////////////////////////////////////
void FastReflectTransmit(SimParamGPU d_simparam, __global const LayerStructGPU *d_layerspecs, 
                                     PhotonStructGPU *photon, __global UINT64 *d_state_Rd_ra, __global UINT64 *d_state_Tt_ra,
                                    UINT64 *rnd_x, UINT32 *rnd_a)
{
  /* Collect all info that depend on the sign of "uz". */
  float cos_crit;
  UINT32 new_layer;
  if (photon->uz > MCML_FP_ZERO)
  {
    cos_crit = d_layerspecs[photon->layer].cos_crit1;
    new_layer = photon->layer+1;
  }
  else
  {
    cos_crit = d_layerspecs[photon->layer].cos_crit0;
    new_layer = photon->layer-1;
  }

  // cosine of the incident angle (0 to 90 deg)
  float ca1 = fabs(photon->uz);

  // The default move is to reflect.
  photon->uz = -photon->uz;

  // Moving this check down to "RFresnel = MCML_FP_ZERO" slows down the
  // application, possibly because every thread is forced to do
  // too much.
  if (ca1 > cos_crit)
  {
    /* Compute the Fresnel reflectance. */

    // incident and transmit refractive index
    float ni = d_layerspecs[photon->layer].n;
    float nt = d_layerspecs[new_layer].n;
    float ni_nt = native_divide(ni, nt);   // reused later

    float sa1 = sqrt(FP_ONE-ca1*ca1);
    if (ca1 > COSZERO) sa1 = MCML_FP_ZERO;
    float sa2 = min(ni_nt * sa1, FP_ONE);
    float uz1 = sqrt(FP_ONE-sa2*sa2);    // uz1 = ca2

    float ca1ca2 = ca1 * uz1;
    float sa1sa2 = sa1 * sa2;
    float sa1ca2 = sa1 * uz1;
    float ca1sa2 = ca1 * sa2;

    // normal incidence: [(1-ni_nt)/(1+ni_nt)]^2
    // We ensure that ca1ca2 = 1, sa1sa2 = 0, sa1ca2 = 1, ca1sa2 = ni_nt
    if (ca1 > COSZERO)
    {
      sa1ca2 = FP_ONE;
      ca1sa2 = ni_nt;
    }

    float cam = ca1ca2 + sa1sa2; /* c- = cc + ss. */
    float sap = sa1ca2 + ca1sa2; /* s+ = sc + cs. */
    float sam = sa1ca2 - ca1sa2; /* s- = sc - cs. */

    float rFresnel = native_divide(sam, sap*cam);
    rFresnel *= rFresnel;
    rFresnel *= (ca1ca2*ca1ca2 + sa1sa2*sa1sa2);

    // In this case, we do not care if "uz1" is exactly 0.
    if (ca1 < COSNINETYDEG || sa2 == FP_ONE) rFresnel = FP_ONE;

    float rand = rand_MWC_co(rnd_x, rnd_a);

    if (rFresnel < rand)
    {
      // The move is to transmit.
      photon->layer = new_layer;

      // Let's do these even if the photon is dead.
      photon->ux *= ni_nt;
      photon->uy *= ni_nt;
      photon->uz = -copysign(uz1, photon->uz);

      if (photon->layer == 0 || photon->layer > d_simparam.num_layers)
      {
        // transmitted
        float uz2 = photon->uz;
        __global UINT64 *ra_arr = d_state_Tt_ra;
        if (photon->layer == 0)
        {
          // diffuse reflectance
          uz2 = -uz2;
          ra_arr = d_state_Rd_ra;
        }

        UINT32 ia = acos(uz2) * FP_TWO * RPI * d_simparam.na;
        UINT32 ir = native_divide(sqrt(photon->x*photon->x+photon->y*photon->y), d_simparam.dr);
        if (ir >= d_simparam.nr) ir = d_simparam.nr - 1;

        //AtomicAddULL(&ra_arr[ia * d_simparam.nr + ir],
        //  (UINT32)(photon->w * WEIGHT_SCALE));
        
        atomic_add(&ra_arr[ia * d_simparam.nr + ir], (UINT32)(photon->w * WEIGHT_SCALE));

        // Kill the photon.
        photon->w = MCML_FP_ZERO;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
//   Computing the scattering angle and new direction by 
//	 sampling the polar deflection angle theta and the
// 	 azimuthal angle psi.
//////////////////////////////////////////////////////////////////////////////
void Spin(float g, PhotonStructGPU *photon,
                     UINT64 *rnd_x, UINT32 *rnd_a)
{
  float cost, sint; // cosine and sine of the polar deflection angle theta
  float cosp, sinp; // cosine and sine of the azimuthal angle psi
  float psi;
  float SIGN;
  float temp;
  float last_ux, last_uy, last_uz;
  float rand;

  /***********************************************************
  *	>>>>>>> SpinTheta
  *  Choose (sample) a new theta angle for photon propagation
  *	according to the anisotropy.
  *
  *	If anisotropy g is 0, then
  *		cos(theta) = 2*rand-1.
  *	otherwise
  *		sample according to the Henyey-Greenstein function.
  *
  *	Returns the cosine of the polar deflection angle theta.
  ****/

  rand = rand_MWC_co(rnd_x, rnd_a);

  cost = FP_TWO * rand - FP_ONE;

  if (g != MCML_FP_ZERO)
  {
    temp = native_divide((FP_ONE - g * g), FP_ONE + g*cost);
    cost = native_divide(FP_ONE + g * g - temp*temp, FP_TWO * g);
    cost = max(cost, -FP_ONE);
    cost = min(cost, FP_ONE);
  }
  sint = sqrt(FP_ONE - cost * cost);

  /* spin psi 0-2pi. */
  rand = rand_MWC_co(rnd_x, rnd_a);

  psi = FP_TWO * PI_const * rand;
  sinp = sincos(psi, &cosp);

  float stcp = sint * cosp;
  float stsp = sint * sinp;

  last_ux = photon->ux;
  last_uy = photon->uy;
  last_uz = photon->uz;

  if (fabs(last_uz) > COSZERO) 
    /* normal incident. */
  {
    photon->ux = stcp;
    photon->uy = stsp;
    SIGN = ((last_uz) >= MCML_FP_ZERO ? FP_ONE : -FP_ONE);
    photon->uz = cost * SIGN;
  }
  else 
    /* regular incident. */
  {
    temp = rsqrt(FP_ONE - last_uz * last_uz);
    photon->ux = (stcp * last_ux * last_uz - stsp * last_uy) * temp
      + last_ux * cost;
    photon->uy = (stcp * last_uy * last_uz + stsp * last_ux) * temp
      + last_uy * cost;
    photon->uz = native_divide(-stcp, temp) + last_uz * cost;
  }
}


void NewSpin(float g, PhotonStructGPU *photon,
                     UINT64 *rnd_x, UINT32 *rnd_a)
{
  float cost, sint; // cosine and sine of the polar deflection angle theta
  float cosp, sinp; // cosine and sine of the azimuthal angle psi
  float psi;
  float SIGN;
  float temp;
  float last_ux, last_uy, last_uz;
  float last_ax, last_ay, last_az;
  float last_bx, last_by, last_bz;
  float rand;

  /***********************************************************
  *	>>>>>>> SpinTheta
  *  Choose (sample) a new theta angle for photon propagation
  *	according to the anisotropy.
  *
  *	If anisotropy g is 0, then
  *		cos(theta) = 2*rand-1.
  *	otherwise
  *		sample according to the Henyey-Greenstein function.
  *
  *	Returns the cosine of the polar deflection angle theta.
  ****/

  rand = rand_MWC_co(rnd_x, rnd_a);

  cost = FP_TWO * rand - FP_ONE;

  if (g != MCML_FP_ZERO)
  {
    temp = native_divide((FP_ONE - g * g), FP_ONE + g*cost);
    cost = native_divide(FP_ONE + g * g - temp*temp, FP_TWO * g);
    cost = max(cost, -FP_ONE);
    cost = min(cost, FP_ONE);
  }
  sint = sqrt(FP_ONE - cost * cost);

  /* spin psi 0-2pi. */
  rand = rand_MWC_co(rnd_x, rnd_a);

  psi = FP_TWO * PI_const * rand;
  sinp = sincos(psi, &cosp);

  float stcp = sint * cosp;
  float stsp = sint * sinp;

  last_ux = photon->ux;
  last_uy = photon->uy;
  last_uz = photon->uz;

  last_ax = photon->ax;
  last_ay = photon->ay;
  last_az = photon->az;
  
  if (fabs(last_uz) > COSZERO) 
    /* normal incident. */
  {
    photon->ux = stcp;
    photon->uy = stsp;
    SIGN = ((last_uz) >= MCML_FP_ZERO ? FP_ONE : -FP_ONE);
    photon->uz = cost * SIGN;
  }
  else{
     photon->ux = dot_product(cost, -sint*cosp, sint*sinp, last_ux, last_ax, last_bx) ;
     photon->uy = dot_product(cost, -sint*cosp, sint*sinp, last_uy, last_ay, last_by) ;
     photon->uz = dot_product(cost, -sint*cosp, sint*sinp, last_uz,last_az, last_bz) ;
       
     photon->ax = dot_product(sint, cost*cosp, -cost*sinp, last_ux, last_ax, last_bx) ;
     photon->ay = dot_product(sint, cost*cosp, -cost*sinp, last_uy, last_ay, last_by) ;
     photon->az = dot_product(sint, cost*cosp, -cost*sinp, last_uz, last_az, last_bz) ;
       
     photon->bx = dot_product(0, sinp, cosp, last_ux, last_ax, last_bx) ;
     photon->by = dot_product(0, sinp, cosp, last_uy, last_ay, last_by) ;
     photon->bz = dot_product(0, sinp, cosp, last_uz, last_az, last_bz) ;
  }
  float temp1 = photon->ux * photon->ux + photon->uy*photon->uy + photon->uz*photon->uz;
  temp= rsqrt(temp1);
  photon->ux=photon->ux*temp;
  photon->uy=photon->uy*temp;
  photon->uz=photon->uz*temp;
 
  temp1 = photon->ax * photon->ax + photon->ay*photon->ay + photon->az*photon->az;
  temp= rsqrt(temp1);
  photon->ax=photon->ax*temp;
  photon->ay=photon->ay*temp;
  photon->az=photon->az*temp;
 // 
  temp1 = photon->bx * photon->bx + photon->by*photon->by + photon->bz*photon->bz;
  temp= rsqrt(temp1);
  photon->bx=photon->bx*temp;
  photon->by=photon->by*temp;
  photon->bz=photon->bz*temp;
}
//////////////////////////////////////////////////////////////////////////////
//   Initialize thread states (tstates), created to allow a large 
//   simulation to be broken up into batches 
//   (avoiding display driver time-out errors)
//////////////////////////////////////////////////////////////////////////////
__kernel void InitThreadState(__global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                 __global float *tstates_photon_ux, __global float *tstates_photon_uy, __global float *tstates_photon_uz,
                                 __global float *tstates_photon_w, __global float *tstates_photon_sleft,
                                 __global UINT32 *tstates_photon_layer, __global UINT32 *tstates_is_active, __global SimParamGPU *d_simparam_addr)  
{
  PhotonStructGPU photon_temp; 
  SimParamGPU d_simparam = d_simparam_addr[0];
  // Initialize the photon and copy into photon_<parameter x>
  LaunchPhoton(&photon_temp, d_simparam);

  // This is the unique ID for each thread (or thread ID = tid)
  UINT32 tid = get_global_id(0);

  tstates_photon_x[tid] = photon_temp.x;
  tstates_photon_y[tid] = photon_temp.y;
  tstates_photon_z[tid] = photon_temp.z;
  tstates_photon_ux[tid] = photon_temp.ux;
  tstates_photon_uy[tid] = photon_temp.uy;
  tstates_photon_uz[tid] = photon_temp.uz;
  tstates_photon_w[tid] = photon_temp.w;
  tstates_photon_sleft[tid] = photon_temp.sleft;
  tstates_photon_layer[tid] = photon_temp.layer;
//
  tstates_is_active[tid] = 1;
}

//////////////////////////////////////////////////////////////////////////////
//   Save thread states (tstates), by copying the current photon 
//   data from registers into global memory
//////////////////////////////////////////////////////////////////////////////
 void SaveThreadState(__global UINT64 *d_state_x, __global UINT32 *d_state_a,
                                 __global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                 __global float *tstates_photon_ux, __global float *tstates_photon_uy, __global float *tstates_photon_uz,
                                 __global float *tstates_photon_w, __global float *tstates_photon_sleft,
                                 __global UINT32 *tstates_photon_layer, __global UINT32 *tstates_is_active,  
                                   PhotonStructGPU *photon,
                                   UINT64 rnd_x, UINT32 rnd_a,
                                   UINT32 is_active)
{
  UINT32 tid = get_global_id(0);

  d_state_x[tid] = rnd_x;
  d_state_a[tid] = rnd_a;

  tstates_photon_x[tid] = photon->x;
  tstates_photon_y[tid] = photon->y;
  tstates_photon_z[tid] = photon->z;
  tstates_photon_ux[tid] = photon->ux;
  tstates_photon_uy[tid] = photon->uy;
  tstates_photon_uz[tid] = photon->uz;
  tstates_photon_w[tid] = photon->w;
  tstates_photon_sleft[tid] = photon->sleft;
  tstates_photon_layer[tid] = photon->layer;

  tstates_is_active[tid] = is_active;
}

//////////////////////////////////////////////////////////////////////////////
//   Restore thread states (tstates), by copying the latest photon 
//   data from global memory back into the registers
//////////////////////////////////////////////////////////////////////////////
 void RestoreThreadState(__global UINT64 *d_state_x, __global UINT32 *d_state_a,
                                 __global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                 __global float *tstates_photon_ux, __global float *tstates_photon_uy, __global float *tstates_photon_uz,
                                 __global float *tstates_photon_w, __global float *tstates_photon_sleft,
                                 __global UINT32 *tstates_photon_layer, __global UINT32 *tstates_is_active,  
                                   PhotonStructGPU *photon,
                                   UINT64 *rnd_x, UINT32 *rnd_a,
                                   UINT32 *is_active)
{
  UINT32 tid = get_global_id(0);

  *rnd_x = d_state_x[tid];
  *rnd_a = d_state_a[tid];

  photon->x = tstates_photon_x[tid];
  photon->y = tstates_photon_y[tid];
  photon->z = tstates_photon_z[tid];
  photon->ux = tstates_photon_ux[tid];
  photon->uy = tstates_photon_uy[tid];
  photon->uz = tstates_photon_uz[tid];
  photon->w = tstates_photon_w[tid];
  photon->sleft = tstates_photon_sleft[tid];
  photon->layer = tstates_photon_layer[tid];

  *is_active = tstates_is_active[tid];
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//   Main Kernel for MCML (Calls the above inline device functions)
//////////////////////////////////////////////////////////////////////////////

__kernel void MCMLKernel(__global const SimParamGPU *d_simparam_addr,__global const LayerStructGPU *d_layerspecs,
                                  //__global SimState d_state, 
                                  __global UINT32 *d_state_n_photons_left_addr, __global UINT64 *d_state_x, 
                                  __global UINT32 *d_state_a,__global UINT64 *d_state_Rd_ra, 
                                  __global UINT64 *d_state_A_rz, __global UINT64 *d_state_Tt_ra,
                                  //__global GPUThreadStates tstates
                                 __global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                 __global float *tstates_photon_ux, __global float *tstates_photon_uy, __global float *tstates_photon_uz,
                                 __global float *tstates_photon_w, __global float *tstates_photon_sleft,
                                 __global UINT32 *tstates_photon_layer, __global UINT32 *tstates_is_active  
                                  )
{
  // photon structure stored in registers
  PhotonStructGPU photon; 
  SimParamGPU d_simparam = d_simparam_addr[0];
  //UINT32 d_state_n_photons_left = d_state_n_photons_left_addr[0];
  // random number seeds
  UINT64 rnd_x;
  UINT32 rnd_a;

  // Flag to indicate if this thread is active
  UINT32 is_active;

  // Restore the thread state from global memory.
  RestoreThreadState(d_state_x, d_state_a, 
                                             tstates_photon_x, tstates_photon_y, tstates_photon_z,
                                             tstates_photon_ux, tstates_photon_uy, tstates_photon_uz, 
                                             tstates_photon_w, tstates_photon_sleft, 
                                             tstates_photon_layer, tstates_is_active,
                                             &photon, &rnd_x, &rnd_a, &is_active);


  for (int iIndex = 0; iIndex < NUM_STEPS; ++iIndex)
  {
    // Only process photon if the thread is active.
    if (is_active)
    {
      //>>>>>>>>> StepSizeInTissue() in MCML
      ComputeStepSize(&photon,&rnd_x, &rnd_a, d_layerspecs);

      //>>>>>>>>> HitBoundary() in MCML

      photon.hit = HitBoundary(&photon, d_layerspecs);

      Hop(&photon);

      if (photon.hit){
        FastReflectTransmit(d_simparam, d_layerspecs,
                                        &photon, d_state_Rd_ra, d_state_Tt_ra, &rnd_x, &rnd_a);
      
        }
      else
      {
        Drop (&photon, d_state_A_rz, d_layerspecs, d_simparam); 


        Spin(d_layerspecs[photon.layer].g, &photon, &rnd_x, &rnd_a);
        //NewSpin(d_layerspecs[photon.layer].g, &photon, &rnd_x, &rnd_a);
      }

      /***********************************************************
      *  >>>>>>>>> Roulette()
      *  If the photon weight is small, the photon packet tries
      *  to survive a roulette.
      ****/
      if (photon.w < WEIGHT)
      {
        float rand = rand_MWC_co(&rnd_x, &rnd_a);

        // This photon survives the roulette.
        if (photon.w != MCML_FP_ZERO && rand < CHANCE)
          photon.w *= (FP_ONE / CHANCE);
        // This photon is terminated.
        else if (atomic_sub(d_state_n_photons_left_addr, 1) > NUM_THREADS){
            
          LaunchPhoton(&photon, d_simparam); // Launch a new photon.
        
        }
        // No need to process any more photons.
        else{
            is_active = 0;

        }
      }
    }

    //////////////////////////////////////////////////////////////////////
  } // end of the main loop

  barrier(CLK_GLOBAL_MEM_FENCE);
//  d_state_n_photons_left_addr[0]=0;
//////////////////////////////////////////////////////////////////////////

//  d_state_n_photons_left_addr[0] = 0;
// Save the thread state to the global memory.
  SaveThreadState(d_state_x, d_state_a, 
                                             tstates_photon_x, tstates_photon_y, tstates_photon_z,
                                             tstates_photon_ux, tstates_photon_uy, tstates_photon_uz, 
                                             tstates_photon_w, tstates_photon_sleft, 
                                             tstates_photon_layer, tstates_is_active,
                                             &photon, rnd_x, rnd_a, is_active);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

