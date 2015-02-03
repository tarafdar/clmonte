typedef ulong UINT64CL;
typedef uint UINT32CL;
// Critical weight for roulette
#define WEIGHT 1E-4F        
#define WEIGHT_SCALE 8388608//12000000

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
float rand_MWC_co(UINT64CL* x,UINT32CL* a)
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
float rand_MWC_oc(UINT64CL* x, UINT32CL* a)
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

float GetCosCrit(float ni, float nt)
{
  float sin_crit = native_divide(nt,ni);
  float cos_crit = sqrt(FP_ONE - sin_crit*sin_crit);
  return cos_crit;
}

__kernel void InitThreadState(__global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                 __global float *tstates_photon_dx, __global float *tstates_photon_dy, __global float *tstates_photon_dz,
                                 __global float *tstates_photon_w,
                                 __global UINT32CL *tstates_photon_tetra_id, __global UINT32CL *tstates_photon_mat_id, 
                                 __global UINT32CL *tstates_is_active, __global const SimParamGPU *d_simparam_addr,
                                 __global UINT64CL *d_state_x, __global UINT32CL *d_state_a)
{
  Packet pkt; 
  SimParamGPU d_simparam = *d_simparam_addr;
  // This is the unique ID for each thread (or thread ID = tid)
  UINT32CL tid = get_global_id(0);
  UINT64CL rnd_x=d_state_x[tid];
  UINT32CL rnd_a=d_state_a[tid];
  // Initialize the photon and copy into photon_<parameter x>
  
  float rand, theta, phi, sinp, cosp, sint, cost;  
  rand = rand_MWC_co(&rnd_x, &rnd_a);
  theta = PI_const * rand;
  rand = rand_MWC_co(&rnd_x, &rnd_a);
  phi = FP_TWO * PI_const * rand;

  sint = sincos(theta, &cost);
  sinp = sincos(phi, &cosp);

  tstates_photon_x[tid] = d_simparam.originX;
  tstates_photon_y[tid] = d_simparam.originY;
  tstates_photon_z[tid] = d_simparam.originZ;
  tstates_photon_dx[tid] = sinp*cost;
  tstates_photon_dy[tid] = sinp*sint;
  tstates_photon_dz[tid] = cosp;
  tstates_photon_w[tid] = FP_ONE;
  tstates_photon_tetra_id[tid] = d_simparam.init_tetraID;
//
  tstates_is_active[tid] = 1;
}

//////////////////////////////////////////////////////////////////////////////
//   Save thread states (tstates), by copying the current photon 
//   data from registers into global memory
//////////////////////////////////////////////////////////////////////////////
void SaveThreadState(__global UINT64CL *d_state_x, __global UINT32CL *d_state_a,
                                 __global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                 __global float *tstates_photon_dx, __global float *tstates_photon_dy, __global float *tstates_photon_dz,
                                 __global float *tstates_photon_w,
                                 __global UINT32CL *tstates_photon_tetra_id, __global UINT32CL *tstates_photon_mat_id,
                                 __global UINT32CL *tstates_is_active,  
                                 Packet *pkt, UINT64CL rnd_x, UINT32CL rnd_a, UINT32CL is_active)
{
  UINT32CL tid = get_global_id(0);

  d_state_x[tid] = rnd_x;
  d_state_a[tid] = rnd_a;
  
  tstates_photon_x[tid] = pkt->x;
  tstates_photon_y[tid] = pkt->y;
  tstates_photon_z[tid] = pkt->z;
  tstates_photon_dx[tid] = pkt->dx;
  tstates_photon_dy[tid] = pkt->dy;
  tstates_photon_dz[tid] = pkt->dz;
  tstates_photon_w[tid] = pkt->w;
  tstates_photon_tetra_id[tid] = pkt->tetraID;

  tstates_is_active[tid] = is_active;
}

//////////////////////////////////////////////////////////////////////////////
//   Restore thread states (tstates), by copying the latest photon 
//   data from global memory back into the registers
//////////////////////////////////////////////////////////////////////////////
void RestoreThreadState(__global UINT64CL *d_state_x, __global UINT32CL *d_state_a,
                                 __global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                 __global float *tstates_photon_dx, __global float *tstates_photon_dy, __global float *tstates_photon_dz,
                                 __global float *tstates_photon_w,
                                 __global UINT32CL *tstates_photon_tetra_id, __global UINT32CL *tstates_photon_mat_id,
                                 __global UINT32CL *tstates_is_active,  
                                 Packet *pkt, UINT64CL *rnd_x, UINT32CL *rnd_a, UINT32CL *is_active)
{
  UINT32CL tid = get_global_id(0);
  
  *rnd_x = d_state_x[tid];
  *rnd_a = d_state_a[tid];
  
  pkt->x = tstates_photon_x[tid];
  pkt->y = tstates_photon_y[tid];
  pkt->z = tstates_photon_z[tid];
  pkt->dx = tstates_photon_dx[tid];
  pkt->dy = tstates_photon_dy[tid];
  pkt->dz = tstates_photon_dz[tid];
  pkt->w = tstates_photon_w[tid];
  pkt->tetraID = tstates_photon_tetra_id[tid];

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

  *is_active = tstates_is_active[tid];
}

//////////////////////////////////////////////////////////////////////////////
//   Main Kernel for MCML (Calls the above inline device functions)
//////////////////////////////////////////////////////////////////////////////

__kernel void MCMLKernel(__global const SimParamGPU *d_simparam_addr,
                                  __global int *d_state_n_photons_left_addr, __global UINT64CL *d_state_x, 
                                  __global UINT32CL *d_state_a, __global const Tetra *d_tetra_mesh,
                                  __global const Material *d_materialspecs, __global UINT64CL *absorption,__global UINT64CL *transmittance, __global float *debug,
                                   //__global GPUThreadStates tstates
                                 __global float *tstates_photon_x, __global float *tstates_photon_y, __global float *tstates_photon_z,
                                 __global float *tstates_photon_dx, __global float *tstates_photon_dy, __global float *tstates_photon_dz,
                                 __global float *tstates_photon_w,
                                 __global UINT32CL *tstates_photon_tetra_id, __global UINT32CL *tstates_photon_mat_id, 
                                 __global UINT32CL *tstates_is_active  )
{
  // packet structure stored in registers
  Packet pkt; 
  SimParamGPU d_simparam = *d_simparam_addr;
  UINT32CL tid = get_global_id(0);

  UINT64CL rnd_x;
  UINT32CL rnd_a;
  
  // Flag to indicate if this thread is active
  UINT32CL is_active ;
  RestoreThreadState(d_state_x, d_state_a, 
     tstates_photon_x, tstates_photon_y, tstates_photon_z,
     tstates_photon_dx, tstates_photon_dy, tstates_photon_dz, 
     tstates_photon_w, 
     tstates_photon_tetra_id, tstates_photon_mat_id, tstates_is_active,
     &pkt, &rnd_x, &rnd_a, &is_active);
  float rand;
  
  for (int iIndex = 0; iIndex < 32768; ++iIndex)
  {
    // Only process packet if the thread is active.
    if (is_active)
    {
      //ComputeStepSize(&pkt,&d_state_x[tid], &d_state_a[tid], d_tetra_mesh, d_materialspecs);
      if (pkt.s == MCML_FP_ZERO)
      {
        rand = rand_MWC_oc(&rnd_x, &rnd_a);
        pkt.s = -log(rand);
      }
      
      float cosdn[4], dis[4];
      Tetra tetra = d_tetra_mesh[pkt.tetraID];
      Material mat = d_materialspecs[tetra.matID];
      
      cosdn[0] = dot_product(pkt.dx, pkt.dy, pkt.dz, tetra.face[0][0], tetra.face[0][1], tetra.face[0][2]);
      cosdn[1] = dot_product(pkt.dx, pkt.dy, pkt.dz, tetra.face[1][0], tetra.face[1][1], tetra.face[1][2]);
      cosdn[2] = dot_product(pkt.dx, pkt.dy, pkt.dz, tetra.face[2][0], tetra.face[2][1], tetra.face[2][2]);
      cosdn[3] = dot_product(pkt.dx, pkt.dy, pkt.dz, tetra.face[3][0], tetra.face[3][1], tetra.face[3][2]);

      dis[0] = dot_product(pkt.x,pkt.y,pkt.z,tetra.face[0][0],tetra.face[0][1],tetra.face[0][2])-tetra.face[0][3];
      dis[1] = dot_product(pkt.x,pkt.y,pkt.z,tetra.face[1][0],tetra.face[1][1],tetra.face[1][2])-tetra.face[1][3];
      dis[2] = dot_product(pkt.x,pkt.y,pkt.z,tetra.face[2][0],tetra.face[2][1],tetra.face[2][2])-tetra.face[2][3];
      dis[3] = dot_product(pkt.x,pkt.y,pkt.z,tetra.face[3][0],tetra.face[3][1],tetra.face[3][2])-tetra.face[3][3];

      //if the packet is moving away from a face, its distance to that face's intersection is infinity.
      dis[0] = cosdn[0]>=0 ? FLT_MAX : -native_divide(dis[0], cosdn[0]);
      dis[1] = cosdn[1]>=0 ? FLT_MAX : -native_divide(dis[1], cosdn[1]);
      dis[2] = cosdn[2]>=0 ? FLT_MAX : -native_divide(dis[2], cosdn[2]);
      dis[3] = cosdn[3]>=0 ? FLT_MAX : -native_divide(dis[3], cosdn[3]);

      UINT32CL localMinIndex1, localMinIndex2, minIndex;
      localMinIndex1 = dis[0]<dis[1] ? 0 : 1;
      localMinIndex2 = dis[2]<dis[3] ? 2 : 3;
      minIndex = dis[localMinIndex1]<dis[localMinIndex2] ? localMinIndex1 : localMinIndex2;
  
      pkt.faceIndexToHit = minIndex;
      UINT32CL nextTetraID = tetra.adjTetras[minIndex];
  
      float canmove = pkt.s * mat.rmu_as;
      if(dis[minIndex] > canmove)
      {
        pkt.s = MCML_FP_ZERO;
        pkt.x += canmove *pkt.dx;
        pkt.y += canmove *pkt.dy;
        pkt.z += canmove *pkt.dz;
        
        float dw = pkt.w*mat.absfrac;
        pkt.w = pkt.w - dw;
        atomic_add(&(absorption[pkt.tetraID]), (UINT64CL)(dw*WEIGHT_SCALE));
        
        if (pkt.w < WEIGHT)
        {
          rand = rand_MWC_co(&rnd_x, &rnd_a);

          // This pkt survives the roulette.
          if (pkt.w != MCML_FP_ZERO && rand < CHANCE)
            pkt.w *= (FP_ONE / CHANCE);
          // This pkt is terminated.
          else if (atomic_sub(d_state_n_photons_left_addr, 1) > 0)
          {
          
          
            pkt.x = d_simparam.originX;
            pkt.y = d_simparam.originY;
            pkt.z = d_simparam.originZ;
            pkt.tetraID = d_simparam.init_tetraID;
  
            float rand, theta, phi, sinp, cosp, sint, cost;  
            rand = rand_MWC_co(&rnd_x, &rnd_a);
            theta = PI_const * rand;
            rand = rand_MWC_co(&rnd_x, &rnd_a);
            phi = FP_TWO * PI_const * rand;

            sint = sincos(theta, &cost);
            sinp = sincos(phi, &cosp);

            pkt.dx = sinp * cost; 
            pkt.dy = sinp * sint;
            pkt.dz = cosp;
            pkt.w = FP_ONE;

            // vector a = (vector d) cross (positive z axis unit vector) and normalize it 
            float4 d = (float4)(pkt.dx, pkt.dy, pkt.dz, 0);
            float4 crossProduct;
            if(d.x==0 && d.y==0)
              crossProduct = cross(d, (float4)(1,0,0,0));
            else
              crossProduct = cross(d, (float4)(0,0,1,0));
            float4 a = normalize(crossProduct);
            pkt.ax = a.x;
            pkt.ay = a.y;
            pkt.az = a.z;

            // vector b = (vector d) cross (vector a)
            float4 b = cross(d,a);
            pkt.bx = b.x;
            pkt.by = b.y;
            pkt.bz = b.z;
          
          
          }
          else
            is_active = 0;
        }
        else
        {
          float cost, sint, cosp, sinp;
          float last_dx, last_dy, last_dz, last_ax, last_ay, last_az, last_bx, last_by, last_bz;
          rand = FP_TWO * rand_MWC_co(&rnd_x, &rnd_a) - FP_ONE;
          cost = mat.HGCoeff1 - native_divide(mat.HGCoeff2, (1-mat.g*rand)*(1-mat.g*rand));
          sint = sqrt(FP_ONE - cost * cost);

          /* spin psi 0-2pi. */
          rand = rand_MWC_co(&rnd_x, &rnd_a);
          sinp = sincos(FP_TWO * PI_const * rand, &cosp);

          float stcp = sint * cosp;
          float stsp = sint * sinp;
          float ctcp = cost * cosp;
          float ctsp = cost * sinp;

          last_dx = pkt.dx;
          last_dy = pkt.dy;
          last_dz = pkt.dz;
          last_ax = pkt.ax;
          last_ay = pkt.ay;
          last_az = pkt.az;
          last_bx = pkt.bx;
          last_by = pkt.by;
          last_bz = pkt.bz;
          pkt.dx = cost*last_dx - stcp*last_ax + stsp*last_bx;
          pkt.dy = cost*last_dy - stcp*last_ay + stsp*last_by;
          pkt.dz = cost*last_dz - stcp*last_az + stsp*last_bz;
          pkt.ax = sint*last_dx +ctcp*last_ax - ctsp*last_bx;
          pkt.ay = sint*last_dy +ctcp*last_ay - ctsp*last_by;
          pkt.az = sint*last_dz +ctcp*last_az - ctsp*last_bz;
          pkt.bx = sinp*last_ax + cosp*last_bx;
          pkt.by = sinp*last_ay + cosp*last_by;
          pkt.bz = sinp*last_az + cosp*last_bz;
        }
      }
      else
      {
        pkt.s = pkt.s - dis[minIndex] * mat.mu_as;
        pkt.x += dis[minIndex] * pkt.dx;
        pkt.y += dis[minIndex] * pkt.dy;
        pkt.z += dis[minIndex] * pkt.dz;
        
        
        Tetra nextTetra = d_tetra_mesh[nextTetraID];
        float ni, nt; //refractive indices
        ni = mat.n;  
        nt = d_materialspecs[nextTetra.matID].n;  
  
        if (ni==nt)
        {
          pkt.tetraID = nextTetraID;
          if (nextTetraID == 0)
          {
            atomic_add(&(transmittance[(pkt.tetraID - 1) * 4 + pkt.faceIndexToHit]), (UINT64CL)(pkt.w * WEIGHT_SCALE));
            pkt.w = MCML_FP_ZERO;
          }
        }
        else
        {
        float crit_cos=nt<ni ? 0 : GetCosCrit(ni,nt);
  
        float *normal = tetra.face[pkt.faceIndexToHit];
        float ca1 = -dot_product(pkt.dx,pkt.dy,pkt.dz,normal[0],normal[1],normal[2]); //ca1 is cos(theta incidence) 
        float sa1 = sqrt(FP_ONE-ca1*ca1); //sa1 is sin(theta incidence)
        if (ca1 > COSZERO) sa1 = MCML_FP_ZERO; 
  
        if (ca1 <= crit_cos)	//total internal reflection occurs
        {
          //reflected direction = original_direction - 2(original_direction dot normal)*normal
          pkt.dx = pkt.dx + 2*ca1*normal[0]; 
          pkt.dy = pkt.dy + 2*ca1*normal[1];
          pkt.dz = pkt.dz + 2*ca1*normal[2];
          // auxilary updating function
          // vector a = (vector d) cross (positive z axis unit vector) and normalize it 
          float4 d = (float4)(pkt.dx, pkt.dy, pkt.dz, 0);
          float4 crossProduct;
          if(d.x==0 && d.y==0)
            crossProduct = cross(d, (float4)(1,0,0,0));
          else
            crossProduct = cross(d, (float4)(0,0,1,0));
          float4 a = normalize(crossProduct);
          pkt.ax = a.x;
          pkt.ay = a.y;
          pkt.az = a.z;

          // vector b = (vector d) cross (vector a)
          float4 b = cross(d,a);
          pkt.bx = b.x;
          pkt.by = b.y;
          pkt.bz = b.z;
        }
        else
        {    
          // Rs = [(n1 cos(theta_i) - n2 cos(theta_t))/((n1 cos(theta_i) + n2 cos(theta_t)))] ^ 2
          // Rp = [(n1 cos(theta_t) - n2 cos(theta_i))/((n1 cos(theta_t) + n2 cos(theta_i)))] ^ 2
          float ni_nt = native_divide(ni, nt);
          float sa2 = min(ni_nt * sa1, FP_ONE); //sa2 is sin(theta transmit)
          float ca2 = sqrt(FP_ONE-sa2*sa2); //ca2 is cos(theta transmit)
          float Rs = native_divide(ni*ca1 - nt*ca2, ni*ca1 + nt*ca2);
          Rs *= Rs;
          float Rp = native_divide(ni*ca2 - nt*ca1, ni*ca2 + nt*ca1);
          Rp *= Rp;
          float rFresnel = (Rs+Rp)/2;

          if (ca1 < COSNINETYDEG || sa2 == FP_ONE) rFresnel = FP_ONE;

          rand = rand_MWC_co(&rnd_x, &rnd_a);

          if (rFresnel < rand) //refract
          {
            // refracted direction = n1/n2*original_direction - [n1/n2*cos(theta incidence) + sqrt(1-sin^2(theta transmitt))]*normal
            pkt.dx = ni_nt*(pkt.dx) - (ni_nt*(-ca1)+ca2)*normal[0]; 
            pkt.dy = ni_nt*(pkt.dy) - (ni_nt*(-ca1)+ca2)*normal[1]; 
            pkt.dz = ni_nt*(pkt.dz) - (ni_nt*(-ca1)+ca2)*normal[2]; 
            
            // auxilary updating function
            // vector a = (vector d) cross (positive z axis unit vector) and normalize it 
            float4 d = (float4)(pkt.dx, pkt.dy, pkt.dz, 0);
            float4 crossProduct;
            if(d.x==0 && d.y==0)
              crossProduct = cross(d, (float4)(1,0,0,0));
            else
              crossProduct = cross(d, (float4)(0,0,1,0));
            float4 a = normalize(crossProduct);
            pkt.ax = a.x;
            pkt.ay = a.y;
            pkt.az = a.z;

	        // vector b = (vector d) cross (vector a)
	        float4 b = cross(d,a);
	        pkt.bx = b.x;
	        pkt.by = b.y;
	        pkt.bz = b.z;
      
            if (nextTetraID == 0) {
      
          	  //store surface fluence TODO: find more efficient way of doing this!
          	  atomic_add(&(transmittance[(pkt.tetraID - 1) * 4 + pkt.faceIndexToHit]), (UINT64CL)(pkt.w * WEIGHT_SCALE));
              // Kill the packet.
          	  pkt.w = MCML_FP_ZERO;
            }
            pkt.tetraID = nextTetraID;
          }
          else //reflect
          {
            //reflected direction = original_direction - 2(original_direction dot normal)*normal
            pkt.dx = pkt.dx + 2*ca1*normal[0]; 
            pkt.dy = pkt.dy + 2*ca1*normal[1];
            pkt.dz = pkt.dz + 2*ca1*normal[2];
      
            // auxilary updating function
            // vector a = (vector d) cross (positive z axis unit vector) and normalize it 
            float4 d = (float4)(pkt.dx, pkt.dy, pkt.dz, 0);
            float4 crossProduct;
            if(d.x==0 && d.y==0)
              crossProduct = cross(d, (float4)(1,0,0,0));
            else
              crossProduct = cross(d, (float4)(0,0,1,0));
            float4 a = normalize(crossProduct);
            pkt.ax = a.x;
            pkt.ay = a.y;
            pkt.az = a.z;

            // vector b = (vector d) cross (vector a)
            float4 b = cross(d,a);
            pkt.bx = b.x;
            pkt.by = b.y;
            pkt.bz = b.z; 
          }
        }
        }
      }
    }
  }
  
  //barrier(CLK_GLOBAL_MEM_FENCE); 
  SaveThreadState(d_state_x, d_state_a, 
     tstates_photon_x, tstates_photon_y, tstates_photon_z,
     tstates_photon_dx, tstates_photon_dy, tstates_photon_dz, 
     tstates_photon_w, 
     tstates_photon_tetra_id, tstates_photon_mat_id, tstates_is_active,
     &pkt, rnd_x, rnd_a, is_active);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

