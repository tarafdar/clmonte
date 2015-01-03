typedef ulong UINT64;
typedef uint UINT32;
// Critical weight for roulette
#define WEIGHT 1E-4F        

// scaling factor for packet weight, which is then converted to integer
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
float rand_MWC_co(__global UINT64* x,__global UINT32* a)
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
float rand_MWC_oc(__global UINT64* x,__global UINT32* a)
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
//   Initialize packet position (x, y, z), direction (dx, dy, dz), weight (w), 
//   step size remainder (sleft), current layer (layer), and auxiliary vectors a,b
//   Note: Infinitely narrow beam (pointing in the +z direction = downwards)
//////////////////////////////////////////////////////////////////////////////
void LaunchPacket(Packet *pkt, SimParamGPU d_simparam, __global UINT64 *rnd_x, __global UINT32 *rnd_a)
{
  pkt->x = d_simparam.originX;
  pkt->y = d_simparam.originY;
  pkt->z = d_simparam.originZ;
  pkt->tetraID = d_simparam.init_tetraID;
  
  float rand, theta, phi;  
  rand = rand_MWC_co(rnd_x, rnd_a);
  theta = PI_const * rand;
  phi = FP_TWO * PI_const * rand;
  
  pkt->dx = native_sin(phi) * native_cos(theta); 
  pkt->dy = native_sin(phi) * native_sin(theta);
  pkt->dz = native_cos(phi);
  pkt->w = d_simparam.init_photon_w;
  pkt->sleft = MCML_FP_ZERO;
  pkt->layer = 1;

  // vector a = (vector d) cross (positive z axis unit vector) and normalize it 
  float4 d = (float4)(pkt->dx, pkt->dy, pkt->dz, 0);
  float4 crossProduct = cross(d, (float4)(0,0,1,0));
  float4 a = normalize(crossProduct);
  pkt->ax = a.x;
  pkt->ay = a.y;
  pkt->az = a.z;

  // vector b = (vector d) cross (vector a)
  float4 b = cross(d,a);
  pkt->bx = b.x;
  pkt->by = b.y;
  pkt->bz = b.z;
}

//////////////////////////////////////////////////////////////////////////////
//   Compute the step size for a packet when it is in tissue
//   If sleft is 0, calculate new step size: -log(rnd)/(mua+mus).
//   Otherwise, pick up the leftover in sleft.
//////////////////////////////////////////////////////////////////////////////
void ComputeStepSize(Packet *pkt,
                                __global UINT64 *rnd_x, __global UINT32 *rnd_a, __global const LayerStructGPU *d_layerspecs,
                                __global const Tetra *d_tetra_mesh, __global const Material *d_materials)
{
  // Make a new step if no leftover.
  if (pkt->sleft == MCML_FP_ZERO)
  {
    //TEMP CHANGE: CHANGE BACK!!! CHANGED BACK :)
    float rand = rand_MWC_oc(rnd_x, rnd_a);
    //float rand = 0;
    pkt->s = -log(rand) * d_layerspecs[pkt->layer].rmuas;
    /* for full monte
    UINT32 materialID = d_tetra_mesh[pkt->tetraID].matID;
    pkt->s = -log(rand) * d_materials[materialID].rmuas;
    */    
  }
  else {
    pkt->s = pkt->sleft * d_layerspecs[pkt->layer].rmuas;
    /* for full monte
    UINT32 materialID = d_tetra_mesh[pkt->tetraID].matID;
    pkt->s = pkt->sleft * d_materials[materialID].rmuas;
    */
    pkt->sleft = MCML_FP_ZERO;
  }
}

//////////////////////////////////////////////////////////////////////////////
//   Check if the step size calculated above will cause the packet to hit the 
//   boundary between 2 layers.
//   Return 1 for a hit, 0 otherwise.
//   If the projected step hits the boundary, the packet steps to the boundary
//   and the remainder of the step size is stored in sleft for the next iteration
//////////////////////////////////////////////////////////////////////////////
int HitBoundary(Packet *pkt, __global const LayerStructGPU *d_layerspecs)
{
  /* step size to boundary. */
  float dl_b; 

  /* Distance to the boundary. */
  float z_bound = (pkt->dz > MCML_FP_ZERO) ?
    d_layerspecs[pkt->layer].z1 : d_layerspecs[pkt->layer].z0;
  dl_b = native_divide(z_bound - pkt->z, pkt->dz);     // dl_b > 0

  UINT32 hit_boundary = (pkt->dz != MCML_FP_ZERO) && (pkt->s > dl_b);
  if (hit_boundary)
  {
    // No need to multiply by (mua + mus), as it is later
    // divided by (mua + mus) anyways (in the original version).
    pkt->sleft = (pkt->s - dl_b) * d_layerspecs[pkt->layer].muas;
    pkt->s = dl_b;
  }

  return hit_boundary;
}

//////////////////////////////////////////////////////////////////////////////
//   Move the packet by step size (s) along direction (dx,dy,dz) 
//////////////////////////////////////////////////////////////////////////////
void Hop(Packet *pkt)
{
  pkt->x += pkt->s * pkt->dx;
  pkt->y += pkt->s * pkt->dy;
  pkt->z += pkt->s * pkt->dz;
}

//////////////////////////////////////////////////////////////////////////////
//   Drop a part of the weight of the packet to simulate absorption
//////////////////////////////////////////////////////////////////////////////
void Drop(Packet *pkt, __global UINT64 *g_A_rz, __global const LayerStructGPU *d_layerspecs, SimParamGPU d_simparam)
{
   
  float dwa = pkt->w * d_layerspecs[pkt->layer].mua_muas;
  pkt->w -= dwa;

  // If scoring of absorption is specified (no -A flag in command line)
    UINT32 iz = native_divide(pkt->z, d_simparam.dz);
    UINT32 ir = native_divide(sqrt(pkt->x * pkt->x + pkt->y * pkt->y),d_simparam.dr);

    // Only record if packet is not at the edge!!
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
//   If a packet hits a boundary, determine whether the packet is transmitted
//   into the next layer or reflected back by computing the internal reflectance
//////////////////////////////////////////////////////////////////////////////
void FastReflectTransmit(SimParamGPU d_simparam, __global const LayerStructGPU *d_layerspecs, 
                                     Packet *pkt, __global UINT64 *d_state_Rd_ra, __global UINT64 *d_state_Tt_ra,
                                    __global UINT64 *rnd_x, __global UINT32 *rnd_a)
{
  /* Collect all info that depend on the sign of "dz". */
  float cos_crit;
  UINT32 new_layer;
  if (pkt->dz > MCML_FP_ZERO)
  {
    cos_crit = d_layerspecs[pkt->layer].cos_crit1;
    new_layer = pkt->layer+1;
  }
  else
  {
    cos_crit = d_layerspecs[pkt->layer].cos_crit0;
    new_layer = pkt->layer-1;
  }

  // cosine of the incident angle (0 to 90 deg)
  float ca1 = fabs(pkt->dz);

  // The default move is to reflect.
  pkt->dz = -pkt->dz;

  // Moving this check down to "RFresnel = MCML_FP_ZERO" slows down the
  // application, possibly because every thread is forced to do
  // too much.
  if (ca1 > cos_crit)
  {
    /* Compute the Fresnel reflectance. */

    // incident and transmit refractive index
    float ni = d_layerspecs[pkt->layer].n;
    float nt = d_layerspecs[new_layer].n;
    float ni_nt = native_divide(ni, nt);   // reused later

    float sa1 = sqrt(FP_ONE-ca1*ca1);
    if (ca1 > COSZERO) sa1 = MCML_FP_ZERO;
    float sa2 = min(ni_nt * sa1, FP_ONE);
    float dz1 = sqrt(FP_ONE-sa2*sa2);    // dz1 = ca2

    float ca1ca2 = ca1 * dz1;
    float sa1sa2 = sa1 * sa2;
    float sa1ca2 = sa1 * dz1;
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

    // In this case, we do not care if "dz1" is exactly 0.
    if (ca1 < COSNINETYDEG || sa2 == FP_ONE) rFresnel = FP_ONE;

    float rand = rand_MWC_co(rnd_x, rnd_a);

    if (rFresnel < rand)
    {
      // The move is to transmit.
      pkt->layer = new_layer;

      // Let's do these even if the packet is dead.
      pkt->dx *= ni_nt;
      pkt->dy *= ni_nt;
      pkt->dz = -copysign(dz1, pkt->dz);

      if (pkt->layer == 0 || pkt->layer > d_simparam.num_layers)
      {
        // transmitted
        float dz2 = pkt->dz;
        __global UINT64 *ra_arr = d_state_Tt_ra;
        if (pkt->layer == 0)
        {
          // diffuse reflectance
          dz2 = -dz2;
          ra_arr = d_state_Rd_ra;
        }

        UINT32 ia = acos(dz2) * FP_TWO * RPI * d_simparam.na;
        UINT32 ir = native_divide(sqrt(pkt->x*pkt->x+pkt->y*pkt->y), d_simparam.dr);
        if (ir >= d_simparam.nr) ir = d_simparam.nr - 1;

        //AtomicAddULL(&ra_arr[ia * d_simparam.nr + ir],
        //  (UINT32)(pkt->w * WEIGHT_SCALE));
        
        atomic_add(&ra_arr[ia * d_simparam.nr + ir], (UINT64)(pkt->w * WEIGHT_SCALE));

        // Kill the packet.
        pkt->w = MCML_FP_ZERO;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
//   Computing the scattering angle and new direction by 
//	 sampling the polar deflection angle theta and the
// 	 azimuthal angle psi.
//////////////////////////////////////////////////////////////////////////////
void Spin(float g, Packet *pkt,
                     __global UINT64 *rnd_x, __global UINT32 *rnd_a)
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
  *  Choose (sample) a new theta angle for packet propagation
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

  last_ux = pkt->dx;
  last_uy = pkt->dy;
  last_uz = pkt->dz;

  if (fabs(last_uz) > COSZERO) 
    /* normal incident. */
  {
    pkt->dx = stcp;
    pkt->dy = stsp;
    SIGN = ((last_uz) >= MCML_FP_ZERO ? FP_ONE : -FP_ONE);
    pkt->dz = cost * SIGN;
  }
  else 
    /* regular incident. */
  {
    temp = rsqrt(FP_ONE - last_uz * last_uz);
    pkt->dx = (stcp * last_ux * last_uz - stsp * last_uy) * temp
      + last_ux * cost;
    pkt->dy = (stcp * last_uy * last_uz + stsp * last_ux) * temp
      + last_uy * cost;
    pkt->dz = native_divide(-stcp, temp) + last_uz * cost;
  }
}


void NewSpin(float g, Packet *pkt,
                     __global UINT64 *rnd_x, __global UINT32 *rnd_a)
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
  *  Choose (sample) a new theta angle for packet propagation
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

  last_ux = pkt->dx;
  last_uy = pkt->dy;
  last_uz = pkt->dz;

  last_ax = pkt->ax;
  last_ay = pkt->ay;
  last_az = pkt->az;
  
  if (fabs(last_uz) > COSZERO) 
    /* normal incident. */
  {
    pkt->dx = stcp;
    pkt->dy = stsp;
    SIGN = ((last_uz) >= MCML_FP_ZERO ? FP_ONE : -FP_ONE);
    pkt->dz = cost * SIGN;
  }
  else{
     pkt->dx = dot_product(cost, -sint*cosp, sint*sinp, last_ux, last_ax, last_bx) ;
     pkt->dy = dot_product(cost, -sint*cosp, sint*sinp, last_uy, last_ay, last_by) ;
     pkt->dz = dot_product(cost, -sint*cosp, sint*sinp, last_uz,last_az, last_bz) ;
       
     pkt->ax = dot_product(sint, cost*cosp, -cost*sinp, last_ux, last_ax, last_bx) ;
     pkt->ay = dot_product(sint, cost*cosp, -cost*sinp, last_uy, last_ay, last_by) ;
     pkt->az = dot_product(sint, cost*cosp, -cost*sinp, last_uz, last_az, last_bz) ;
       
     pkt->bx = dot_product(0, sinp, cosp, last_ux, last_ax, last_bx) ;
     pkt->by = dot_product(0, sinp, cosp, last_uy, last_ay, last_by) ;
     pkt->bz = dot_product(0, sinp, cosp, last_uz, last_az, last_bz) ;
  }
  float temp1 = pkt->dx * pkt->dx + pkt->dy*pkt->dy + pkt->dz*pkt->dz;
  temp= rsqrt(temp1);
  pkt->dx=pkt->dx*temp;
  pkt->dy=pkt->dy*temp;
  pkt->dz=pkt->dz*temp;
 
  temp1 = pkt->ax * pkt->ax + pkt->ay*pkt->ay + pkt->az*pkt->az;
  temp= rsqrt(temp1);
  pkt->ax=pkt->ax*temp;
  pkt->ay=pkt->ay*temp;
  pkt->az=pkt->az*temp;
 // 
  temp1 = pkt->bx * pkt->bx + pkt->by*pkt->by + pkt->bz*pkt->bz;
  temp= rsqrt(temp1);
  pkt->bx=pkt->bx*temp;
  pkt->by=pkt->by*temp;
  pkt->bz=pkt->bz*temp;
}
//////////////////////////////////////////////////////////////////////////////
//   Initialize thread states (tstates), created to allow a large 
//   simulation to be broken up into batches 
//   (avoiding display driver time-out errors)
//////////////////////////////////////////////////////////////////////////////
/*__kernel void InitThreadState(__global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                 __global float *tstates_photon_ux, __global float *tstates_photon_uy, __global float *tstates_photon_uz,
                                 __global float *tstates_photon_w, __global float *tstates_photon_sleft,
                                 __global UINT32 *tstates_photon_layer, __global UINT32 *tstates_is_active, __global SimParamGPU *d_simparam_addr)  
{
  Packet photon_temp; 
  SimParamGPU d_simparam = d_simparam_addr[0];
  // Initialize the photon and copy into photon_<parameter x>
  LaunchPacket(&photon_temp, d_simparam);

  // This is the unique ID for each thread (or thread ID = tid)
  UINT32 tid = get_global_id(0);

  tstates_photon_x[tid] = photon_temp.x;
  tstates_photon_y[tid] = photon_temp.y;
  tstates_photon_z[tid] = photon_temp.z;
  tstates_photon_ux[tid] = photon_temp.dx;
  tstates_photon_uy[tid] = photon_temp.dy;
  tstates_photon_uz[tid] = photon_temp.dz;
  tstates_photon_w[tid] = photon_temp.w;
  tstates_photon_sleft[tid] = photon_temp.sleft;
  tstates_photon_layer[tid] = photon_temp.layer;
//
  tstates_is_active[tid] = 1;
}*/

//////////////////////////////////////////////////////////////////////////////
//   Save thread states (tstates), by copying the current photon 
//   data from registers into global memory
//////////////////////////////////////////////////////////////////////////////
 void SaveThreadState(__global UINT64 *d_state_x, __global UINT32 *d_state_a,
                                 __global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                 __global float *tstates_photon_ux, __global float *tstates_photon_uy, __global float *tstates_photon_uz,
                                 __global float *tstates_photon_w, __global float *tstates_photon_sleft,
                                 __global UINT32 *tstates_photon_layer, __global UINT32 *tstates_is_active,  
                                   Packet *pkt,
                                   UINT64 rnd_x, UINT32 rnd_a,
                                   UINT32 is_active)
{
  UINT32 tid = get_global_id(0);

  d_state_x[tid] = rnd_x;
  d_state_a[tid] = rnd_a;

  tstates_photon_x[tid] = pkt->x;
  tstates_photon_y[tid] = pkt->y;
  tstates_photon_z[tid] = pkt->z;
  tstates_photon_ux[tid] = pkt->dx;
  tstates_photon_uy[tid] = pkt->dy;
  tstates_photon_uz[tid] = pkt->dz;
  tstates_photon_w[tid] = pkt->w;
  tstates_photon_sleft[tid] = pkt->sleft;
  tstates_photon_layer[tid] = pkt->layer;

  tstates_is_active[tid] = is_active;
}

//////////////////////////////////////////////////////////////////////////////
//   Restore thread states (tstates), by copying the latest packet 
//   data from global memory back into the registers
//////////////////////////////////////////////////////////////////////////////
 void RestoreThreadState(__global UINT64 *d_state_x, __global UINT32 *d_state_a,
                                 __global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                 __global float *tstates_photon_ux, __global float *tstates_photon_uy, __global float *tstates_photon_uz,
                                 __global float *tstates_photon_w, __global float *tstates_photon_sleft,
                                 __global UINT32 *tstates_photon_layer, __global UINT32 *tstates_is_active,  
                                   Packet *pkt,
                                   UINT64 *rnd_x, UINT32 *rnd_a,
                                   UINT32 *is_active)
{
  UINT32 tid = get_global_id(0);

  *rnd_x = d_state_x[tid];
  *rnd_a = d_state_a[tid];

  pkt->x = tstates_photon_x[tid];
  pkt->y = tstates_photon_y[tid];
  pkt->z = tstates_photon_z[tid];
  pkt->dx = tstates_photon_ux[tid];
  pkt->dy = tstates_photon_uy[tid];
  pkt->dz = tstates_photon_uz[tid];
  pkt->w = tstates_photon_w[tid];
  pkt->sleft = tstates_photon_sleft[tid];
  pkt->layer = tstates_photon_layer[tid];

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
                                  __global UINT64 *d_state_A_rz, __global UINT64 *d_state_Tt_ra, __global const Tetra *d_tetra_mesh,
                                  __global const Material *d_materialspecs
                                  //__global GPUThreadStates tstates
                                // __global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                // __global float *tstates_photon_ux, __global float *tstates_photon_uy, __global float *tstates_photon_uz,
                                // __global float *tstates_photon_w, __global float *tstates_photon_sleft,
                                // __global UINT32 *tstates_photon_layer, __global UINT32 *tstates_is_active  
                                  )
{
  // packet structure stored in registers
  Packet pkt; 
  SimParamGPU d_simparam = d_simparam_addr[0];
  //UINT32 d_state_n_photons_left = d_state_n_photons_left_addr[0];
  // random number seeds
  //UINT64 rnd_x;
  //UINT32 rnd_a;

  UINT32 tid = get_global_id(0);

  //rnd_x = d_state_x[tid];
  //rnd_a = d_state_a[tid];
  LaunchPacket(&pkt, d_simparam, &d_state_x[tid], &d_state_a[tid]); // Launch a new packet.
  // Flag to indicate if this thread is active
  UINT32 is_active ;
  is_active = 1;
  // Restore the thread state from global memory.
//  RestoreThreadState(d_state_x, d_state_a, 
//                                             tstates_photon_x, tstates_photon_y, tstates_photon_z,
//                                             tstates_photon_ux, tstates_photon_uy, tstates_photon_uz, 
//                                             tstates_photon_w, tstates_photon_sleft, 
//                                             tstates_photon_layer, tstates_is_active,
//                                             &pkt, &rnd_x, &rnd_a, &is_active);
//

  for (int iIndex = 0; iIndex < NUM_STEPS; ++iIndex)
  {
    // Only process packet if the thread is active.
    if (is_active)
    {
      //>>>>>>>>> StepSizeInTissue() in MCML
      ComputeStepSize(&pkt,&d_state_x[tid], &d_state_a[tid], d_layerspecs, d_tetra_mesh, d_materialspecs);

      
      //>>>>>>>>> HitBoundary() in MCML

      pkt.hit = HitBoundary(&pkt, d_layerspecs);

      Hop(&pkt);

      if (pkt.hit){
        FastReflectTransmit(d_simparam, d_layerspecs,
                                        &pkt, d_state_Rd_ra, d_state_Tt_ra, &d_state_x[tid], &d_state_a[tid]);
      
        }
      else
      {
        Drop (&pkt, d_state_A_rz, d_layerspecs, d_simparam); 


        Spin(d_layerspecs[pkt.layer].g, &pkt, &d_state_x[tid], &d_state_a[tid]);
        //NewSpin(d_layerspecs[pkt.layer].g, &pkt, &rnd_x, &rnd_a);
      }

      /***********************************************************
      *  >>>>>>>>> Roulette()
      *  If the pkt weight is small, the packet tries
      *  to survive a roulette.
      ****/
      if (pkt.w < WEIGHT)
      {
        float rand = rand_MWC_co(&d_state_x[tid], &d_state_a[tid]);

        // This pkt survives the roulette.
        if (pkt.w != MCML_FP_ZERO && rand < CHANCE)
          pkt.w *= (FP_ONE / CHANCE);
        // This pkt is terminated.
        else if (atomic_sub(d_state_n_photons_left_addr, 1) > NUM_THREADS){
            
          LaunchPacket(&pkt, d_simparam, &d_state_x[tid], &d_state_a[tid]); // Launch a new packet.
        
        }
        // No need to process any more packets.
        else{
            is_active = 0;

        }
      }
    }

    //////////////////////////////////////////////////////////////////////
  } // end of the main loop

  //barrier(CLK_GLOBAL_MEM_FENCE);
//  d_state_n_photons_left_addr[0]=0;
//////////////////////////////////////////////////////////////////////////

//  d_state_n_photons_left_addr[0] = 0;
// Save the thread state to the global memory.
//  SaveThreadState(d_state_x, d_state_a, 
//                                             tstates_photon_x, tstates_photon_y, tstates_photon_z,
//                                             tstates_photon_ux, tstates_photon_uy, tstates_photon_uz, 
//                                             tstates_photon_w, tstates_photon_sleft, 
//                                             tstates_photon_layer, tstates_is_active,
//                                             &pkt, rnd_x, rnd_a, is_active);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

