typedef ulong UINT64CL;
typedef uint UINT32CL;
// Critical weight for roulette
#define WEIGHT 1E-4F        
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
float rand_MWC_co(__global UINT64CL* x,__global UINT32CL* a)
{
  *x=(*x&0xfffffffful)*(*a)+(*x>>32);
  return native_divide(convert_float_rtz((UINT32CL)(*x)),(float)0x100000000);
  
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
float rand_MWC_oc(__global UINT64CL* x,__global UINT32CL* a)
{
  return 1.0f-rand_MWC_co(x,a);
} 

//////////////////////////////////////////////////////////////////////////////
//   AtomicAdd for Unsigned Long Long (ULL) data type
//   Note: 64-bit atomicAdd to global memory are only supported 
//   in graphics card with compute capability 1.2 or above
//////////////////////////////////////////////////////////////////////////////
//void AtomicAddULL(UINT32CL* address, UINT32CL add)
//{
//  if (atomic_add((UINT32CL*)address,add) +add < add)
//  {
//    atomic_add(((UINT32CL*)address)+1, 1U);
//  }
//  //atomic_add(address, (int)add);
//}

//////////////////////////////////////////////////////////////////////////////
//   Initialize packet position (x, y, z), direction (dx, dy, dz), weight (w), 
//   step size remainder (sleft), current layer (layer), and auxiliary vectors a,b
//   Note: Infinitely narrow beam (pointing in the +z direction = downwards)
//////////////////////////////////////////////////////////////////////////////
void LaunchPacket(Packet *pkt, SimParamGPU d_simparam, __global UINT64CL *rnd_x, __global UINT32CL *rnd_a)
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
  pkt->w = FP_ONE;
  pkt->sleft = MCML_FP_ZERO;

  // vector a = (vector d) cross (positive z axis unit vector) and normalize it 
  float4 d = (float4)(pkt->dx, pkt->dy, pkt->dz, 0);
  float4 crossProduct;
  if(d.x==0 && d.y==0)
    crossProduct = cross(d, (float4)(1,0,0,0));
  else
    crossProduct = cross(d, (float4)(0,0,1,0));
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
                                __global UINT64CL *rnd_x, __global UINT32CL *rnd_a,
                                __global const Tetra *d_tetra_mesh, __global const Material *d_materials,
                                __global float *debug, int index)
{
  // Make a new step if no leftover.
  if (pkt->sleft == MCML_FP_ZERO)
  {
    //TEMP CHANGE: CHANGE BACK!!! CHANGED BACK :)
    float rand = rand_MWC_oc(rnd_x, rnd_a);
    //float rand = 0;
    UINT32CL materialID = d_tetra_mesh[pkt->tetraID].matID;
    pkt->s = -log(rand) * d_materials[materialID].rmu_as;
  }
  else {
    UINT32CL materialID = d_tetra_mesh[pkt->tetraID].matID;
    pkt->s = pkt->sleft * d_materials[materialID].rmu_as;
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
int HitBoundary(Packet *pkt, __global const Tetra *d_tetra_mesh, __global float *debug, int i)
{
  /* step size to boundary. */
  //cos of the angle between direction and normal vector.
  //assuming the packet inside the tetrahedron, cosdn<0 when moving towards that face; cosdn>0 when moving away.
  float cosdn[4];
  
  //orthogonal distance from the packet's position to the face
  float orth_dis[4];

  //distance from the packet's position to the face along its direction
  float move_dis[4];

  float4 d = (float4)(pkt->dx, pkt->dy, pkt->dz, 0);
  float4 p = (float4)(pkt->x, pkt->y, pkt->z, 0);
  Tetra tetra = d_tetra_mesh[pkt->tetraID];
  float4 n[4];
  n[0] = (float4)(tetra.face[0][0],tetra.face[0][1],tetra.face[0][2],tetra.face[0][3]);
  n[1] = (float4)(tetra.face[1][0],tetra.face[1][1],tetra.face[1][2],tetra.face[1][3]);
  n[2] = (float4)(tetra.face[2][0],tetra.face[2][1],tetra.face[2][2],tetra.face[2][3]);
  n[3] = (float4)(tetra.face[3][0],tetra.face[3][1],tetra.face[3][2],tetra.face[3][3]);
  cosdn[0] = dot(d,n[0]);
  cosdn[1] = dot(d,n[1]);
  cosdn[2] = dot(d,n[2]);
  cosdn[3] = dot(d,n[3]);

  orth_dis[0] = dot(p,n[0])-n[0].w;
  orth_dis[1] = dot(p,n[1])-n[1].w;
  orth_dis[2] = dot(p,n[2])-n[2].w;
  orth_dis[3] = dot(p,n[3])-n[3].w;

  //if the packet is moving away from a face, its distance to that face's intersection is infinity.
  move_dis[0] = cosdn[0]>=0 ? FLT_MAX : -native_divide(orth_dis[0], cosdn[0]);
  move_dis[1] = cosdn[1]>=0 ? FLT_MAX : -native_divide(orth_dis[1], cosdn[1]);
  move_dis[2] = cosdn[2]>=0 ? FLT_MAX : -native_divide(orth_dis[2], cosdn[2]);
  move_dis[3] = cosdn[3]>=0 ? FLT_MAX : -native_divide(orth_dis[3], cosdn[3]);

  UINT32CL localMinIndex1, localMinIndex2, minIndex;
  localMinIndex1 = move_dis[0]<move_dis[1] ? 0 : 1;
  localMinIndex2 = move_dis[2]<move_dis[3] ? 2 : 3;
  minIndex = move_dis[localMinIndex1]<move_dis[localMinIndex2] ? localMinIndex1 : localMinIndex2;
  
  //TODO in tuning: Should we also store move_dis[minIndex] into Packet to avoid recomputing?
  pkt->faceIndexToHit = minIndex;
  pkt->nextTetraID = tetra.adjTetras[minIndex];
  pkt->matID = tetra.matID;
  if(move_dis[minIndex] > pkt->s)
  {
    pkt->s = pkt->sleft;
    pkt->sleft = 0;
    return 0;	//means won't hit in this step
  }
  else
  {
    pkt->sleft = pkt->s - move_dis[minIndex];
    pkt->s = move_dis[minIndex];
    return 1;	//means will hit the face specified by minIndex
  }
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
void absorb(Packet *pkt, __global const Material *d_materialspecs, __global UINT64CL *absorption, __global float *debug, int iIndex) {
    
    float w0 = pkt->w;
    float dw = w0*(d_materialspecs[pkt->matID].absfrac);

    pkt->w = w0 - dw;

    //float dw = dw*WEIGHT_SCALE;
    // Store the absorbed weight in the material -> score in the absorption array
    // UINT64 used to maintain compatibility with the "atomic_add" method
    atomic_add(&(absorption[pkt->tetraID]), (UINT64CL)(dw*WEIGHT_SCALE));
}

float GetCosCrit(float ni, float nt)
{
  float sin_crit = native_divide(nt,ni);
  float cos_crit = sqrt(FP_ONE - sin_crit*sin_crit);
  return cos_crit;
}


///<return>ID of the tetrahedron that the packet is entering (same as current ID if reflect, otherwise the correct adjacent tetra ID</return>
void ReflectTransmit(SimParamGPU d_simparam, Packet *pkt, __global UINT64CL *transmittance,
                     __global UINT64CL *rnd_x, __global UINT32CL *rnd_a, __global const Tetra *d_tetra_mesh,
                     __global const Material *d_materialspecs, __global float *debug, int iIndex)

{
  //TODO: Should we store the cos_crit in each Tetra for each face?
  Tetra tetra = d_tetra_mesh[pkt->tetraID];
  Tetra nextTetra;
  int aboutToLeaveSurface = 0;
  float ni, nt; //refractive indices
  ni = d_materialspecs[tetra.matID].n;
  
  if (tetra.adjTetras[pkt->faceIndexToHit] == 0 ) {
  	//about to leave the surface
        aboutToLeaveSurface = 1;
  } 
  nextTetra = d_tetra_mesh[tetra.adjTetras[pkt->faceIndexToHit]];
  nt = d_materialspecs[nextTetra.matID].n;  

  float crit_cos=0;	//initialize as cos=0, means total internal reflection never occurs.
  if (nt<ni)	//total internal reflection may occur
  {
    crit_cos = GetCosCrit(ni, nt);
  }
  
  float *normal = tetra.face[pkt->faceIndexToHit];
  float costheta = -dot((float4)(pkt->dx,pkt->dy,pkt->dz,0), (float4)(normal[0],normal[1],normal[2],0));

  float ni_nt = native_divide(ni, nt);

  float ca1 = costheta; //ca1 is cos(theta incidence) 
  float sa1 = sqrt(FP_ONE-ca1*ca1); //sa1 is sin(theta incidence)
  if (ca1 > COSZERO) sa1 = MCML_FP_ZERO;
  float sa2 = min(ni_nt * sa1, FP_ONE); //sa2 is sin(theta transmit)
  float ca2 = sqrt(FP_ONE-sa2*sa2); //ca2 is cos(theta transmit)
  
  if (costheta <= crit_cos)	//total internal reflection occurs
  {
    //reflected direction = original_direction - 2(original_direction dot normal)*normal
    pkt->dx = pkt->dx + 2*ca1*normal[0]; //here the costheta must be the same as the one Li defined during the TIR step, which is negative d dot n, so the plus sign in the middle is consistent with the minus sign in the formula that's in the comment above
    pkt->dy = pkt->dy + 2*ca1*normal[1];
    pkt->dz = pkt->dz + 2*ca1*normal[2];
    // auxilary updating function
    // vector a = (vector d) cross (positive z axis unit vector) and normalize it 
    float4 d = (float4)(pkt->dx, pkt->dy, pkt->dz, 0);
    float4 crossProduct;
    if(d.x==0 && d.y==0)
      crossProduct = cross(d, (float4)(1,0,0,0));
    else
      crossProduct = cross(d, (float4)(0,0,1,0));
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
  else
  {    
    //Jeff's way of calculating Fresnel
    // Rs = [(n1 cos(theta_i) - n2 cos(theta_t))/((n1 cos(theta_i) + n2 cos(theta_t)))] ^ 2
    // Rp = [(n1 cos(theta_t) - n2 cos(theta_i))/((n1 cos(theta_t) + n2 cos(theta_i)))] ^ 2
    float Rs = (ni*ca1 - nt*ca2)/(ni*ca1 + nt*ca2);
    Rs *= Rs;
    float Rp = (ni*ca2 - nt*ca1)/(ni*ca2 + nt*ca1);
    Rp *= Rp;
    float rFresnel = (Rs+Rp)/2;

    if (ca1 < COSNINETYDEG || sa2 == FP_ONE) rFresnel = FP_ONE;

    float rand = rand_MWC_co(rnd_x, rnd_a);

    if (rFresnel < rand) //refract
    {
      // refracted direction = n1/n2*original_direction - [n1/n2*cos(theta incidence) + sqrt(1-sin^2(theta transmitt))]*normal
      pkt->dx = ni_nt*(pkt->dx) - (ni_nt*(-ca1)+ca2)*normal[0]; 
      pkt->dy = ni_nt*(pkt->dy) - (ni_nt*(-ca1)+ca2)*normal[1]; 
      pkt->dz = ni_nt*(pkt->dz) - (ni_nt*(-ca1)+ca2)*normal[2]; 
      
      // auxilary updating function
      // vector a = (vector d) cross (positive z axis unit vector) and normalize it 
      float4 d = (float4)(pkt->dx, pkt->dy, pkt->dz, 0);
      float4 crossProduct;
      if(d.x==0 && d.y==0)
        crossProduct = cross(d, (float4)(1,0,0,0));
      else
        crossProduct = cross(d, (float4)(0,0,1,0));
      float4 a = normalize(crossProduct);
      pkt->ax = a.x;
      pkt->ay = a.y;
      pkt->az = a.z;

      // vector b = (vector d) cross (vector a)
      float4 b = cross(d,a);
      pkt->bx = b.x;
      pkt->by = b.y;
      pkt->bz = b.z;
      
      if (aboutToLeaveSurface == 1) {
      
      	//store surface fluence TODO: find more efficient way of doing this!
      	atomic_add(&(transmittance[(pkt->tetraID - 1) * 4 + pkt->faceIndexToHit]), (UINT64CL)(pkt->w * WEIGHT_SCALE));
      
      
      	// Kill the packet.
      	pkt->w = MCML_FP_ZERO;
      }
      pkt->tetraID = pkt->nextTetraID;
    }
    else //reflect
    {
      //reflected direction = original_direction - 2(original_direction dot normal)*normal
      pkt->dx = pkt->dx + 2*ca1*normal[0]; //here the costheta must be the same as the one Li defined during the TIR step, which is negative d dot n, so the plus sign in the middle is consistent with the minus sign in the formula that's in the comment above
      pkt->dy = pkt->dy + 2*ca1*normal[1];
      pkt->dz = pkt->dz + 2*ca1*normal[2];
      
      // auxilary updating function
      // vector a = (vector d) cross (positive z axis unit vector) and normalize it 
      float4 d = (float4)(pkt->dx, pkt->dy, pkt->dz, 0);
      float4 crossProduct;
      if(d.x==0 && d.y==0)
        crossProduct = cross(d, (float4)(1,0,0,0));
      else
        crossProduct = cross(d, (float4)(0,0,1,0));
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

  }
}

//////////////////////////////////////////////////////////////////////////////
//   Computing the scattering angle and new direction by 
//	 sampling the polar deflection angle theta and the
// 	 azimuthal angle psi.
//////////////////////////////////////////////////////////////////////////////
void Spin(Packet *pkt, __global UINT64CL *rnd_x, __global UINT32CL *rnd_a,
                     __global const Tetra *d_tetra_mesh, __global const Material *d_materialspec)
{
  float cost, sint; // cosine and sine of the polar deflection angle theta
  float cosp, sinp; // cosine and sine of the azimuthal angle psi
  float psi;
  float SIGN;
  float temp;
  float last_dx, last_dy, last_dz, last_ax, last_ay, last_az, last_bx, last_by, last_bz;
  float rand;

  rand = FP_TWO * rand_MWC_co(rnd_x, rnd_a) - FP_ONE;	//rand is sampled from -1 to 1

  Material mat = d_materialspec[pkt->matID];
  cost = mat.HGCoeff1 - native_divide(mat.HGCoeff2, (1-mat.g * rand));
  
  sint = sqrt(FP_ONE - cost * cost);

  /* spin psi 0-2pi. */
  rand = rand_MWC_co(rnd_x, rnd_a);

  psi = FP_TWO * PI_const * rand;
  sinp = sincos(psi, &cosp);

  float stcp = sint * cosp;
  float stsp = sint * sinp;

  last_dx = pkt->dx;
  last_dy = pkt->dy;
  last_dz = pkt->dz;
  last_ax = pkt->ax;
  last_ay = pkt->ay;
  last_az = pkt->az;
  last_bx = pkt->bx;
  last_by = pkt->by;
  last_bz = pkt->bz;
  pkt->dx = cost*last_dx - sint*cosp*last_ax + sint*sinp*last_bx;
  pkt->dy = cost*last_dy - sint*cosp*last_ay + sint*sinp*last_by;
  pkt->dz = cost*last_dz - sint*cosp*last_az + sint*sinp*last_bz;
  pkt->ax = sint*last_dx +cost*cosp*last_ax - cost*sinp*last_bx;
  pkt->ay = sint*last_dy +cost*cosp*last_ay - cost*sinp*last_by;
  pkt->az = sint*last_dz +cost*cosp*last_az - cost*sinp*last_bz;
  pkt->bx = sinp*last_ax + cosp*last_bx;
  pkt->by = sinp*last_ay + cosp*last_by;
  pkt->bz = sinp*last_az + cosp*last_bz;
}

//////////////////////////////////////////////////////////////////////////////
//   Main Kernel for MCML (Calls the above inline device functions)
//////////////////////////////////////////////////////////////////////////////

__kernel void MCMLKernel(__global const SimParamGPU *d_simparam_addr,
                                  __global int *d_state_n_photons_left_addr, __global UINT64CL *d_state_x, 
                                  __global UINT32CL *d_state_a, __global const Tetra *d_tetra_mesh,
                                  __global const Material *d_materialspecs, __global UINT64CL *absorption,__global UINT64CL *transmittance, __global float *debug)
{
  // packet structure stored in registers
  Packet pkt; 
  SimParamGPU d_simparam = *d_simparam_addr;
  UINT32CL tid = get_global_id(0);

  LaunchPacket(&pkt, d_simparam, &d_state_x[tid], &d_state_a[tid]); // Launch a new packet.

  // Flag to indicate if this thread is active
  UINT32CL is_active ;
  is_active = 1;

  for (int iIndex = 0; iIndex < 50000; ++iIndex)
  {
    // Only process packet if the thread is active.
    if (is_active)
    {
      ComputeStepSize(&pkt,&d_state_x[tid], &d_state_a[tid], d_tetra_mesh, d_materialspecs, debug, iIndex);
      
      pkt.hit = HitBoundary(&pkt, d_tetra_mesh, debug, iIndex);

      Hop(&pkt);

      if (pkt.hit)
      {
        ReflectTransmit(d_simparam, &pkt, transmittance, &d_state_x[tid], &d_state_a[tid], d_tetra_mesh, d_materialspecs, debug, iIndex);
        if (pkt.tetraID == 0)	//exited the tissue
        {
          if (atomic_sub(d_state_n_photons_left_addr, 1) > NUM_THREADS)
          {
            LaunchPacket(&pkt, d_simparam, &d_state_x[tid], &d_state_a[tid]); // Launch a new packet.
          }
          // No need to process any more packets.
          else
          {
            is_active = 0;
            break;
          }
        }
      }
      else
      {
        absorb (&pkt, d_materialspecs, absorption, debug, iIndex);
        Spin(&pkt, &d_state_x[tid], &d_state_a[tid], d_tetra_mesh, d_materialspecs);
      }

      // >>>>>>>>> Roulette()
      // If the pkt weight is small, the packet tries
      // to survive a roulette.
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
            break;
        }
      }
    }
  } 
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

