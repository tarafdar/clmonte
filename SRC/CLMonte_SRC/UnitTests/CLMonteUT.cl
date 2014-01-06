
/*	This file is part of CUDAMC.

    CUDAMC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CUDAMC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CUDAMC.  If not, see <http://www.gnu.org/licenses/>.*/


#define NUM_THREADS_PER_BLOCK 320 //Keep above 192 to eliminate global memory access overhead
#define NUM_BLOCKS 84 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
#define NUM_THREADS 1000
#define NUM_THREADS_TIMES_THREE 80640
#define NUMSTEPS_GPU 500000
#define NUMSTEPS_CPU 500000
#define PI 3.14159265f

#define TMAX 2000.0f //[ps] Maximum time of flight
#define DT 10.0f //[ps] Time binning resolution
#define TEMP 201 //ceil(TMAX/DT), precalculated to avoid dynamic memory allocation (fulhack)


#define MUS_MAX 90.0f	//[1/cm]
#define V 0.0214f		//[cm/ps] (c=0.03 [cm/ps] v=c/n) here n=1.4
#define COS_CRIT 0.6999f	//the critical angle for total internal reflection at the border cos_crit=sqrt(1-(nt/ni)^2)
#define N 1.4f


float dot_product (float x1, float y1, float z1, float x2, float y2, float z2) {
	float x = x1*x2;
	float y = y1*y2;
	float z = z1*z2;
    //return  x1*x2 + y1*y2 + z1*z2;
	return  x + y + z;
}    

float divide(float x, float y)
{
	//return (float)x/y;
	return native_divide(x,y);
}

float logf(float x)
{
	//return (float)log(x);
	return native_log(x);

}


float sqrtf(float x)
{

	return (float)sqrt(x);

}

float rsqrtf(float x)
{
	return (float)rsqrt(x);

}
float rand_MWC_co(unsigned long* x,//unsigned int* c,
		       unsigned int* a)
{
		//Generate a random number [0,1)
		//this implementation seems to be faster
		*x=((*x)&0xfffffffful)*(*a)+((*x)>>32);
		return((float)((unsigned int)((*x)&0xfffffffful))/(UINT_MAX));

}//end __device__ rand_MWC_co

float rand_MWC_oc(unsigned long* x,//unsigned int* c,
		       unsigned int* a)
{
		//Generate a random number (0,1]
		*x=((*x)&0xfffffffful)*(*a)+((*x)>>32);
		return(1.0f-(float)((unsigned int)((*x)&0xfffffffful))/(UINT_MAX));
}//end __device__ rand_MWC_oc


void LaunchPhoton(float3* pos, float3* dir, float* t)
{
	pos->x=0.0f;
	pos->y=0.0f;
	pos->z=0.0f;

	dir->x=0.0f;
	dir->y=0.0f;
	dir->z=1.0f;

	*t=0.0f;
}





__kernel void Reflect_test1(__global float* dirx_array,__global float* diry_array, __global float* dirz_array, __global float* posx_array,__global float* posy_array, __global float* posz_array,  __global unsigned int* xd, __global unsigned int* cd, __global unsigned int* ad, __global unsigned int* retval) {

	int global_id = get_global_id(0);
        
	float3 dirvar;
	float3 * dir;
	dirvar.x =dirx_array[global_id]; 
	dirvar.y =diry_array[global_id]; 
	dirvar.z =dirz_array[global_id];
	dir = &dirvar;
 
	float3 posvar;
	float3 * pos;
	posvar.x =posx_array[global_id]; 
	posvar.y =posy_array[global_id]; 
	posvar.z =posz_array[global_id];
	pos = &posvar;

	
	unsigned long int x = cd[global_id];
	x=(x<<32) + xd[global_id];


        unsigned int a=ad[global_id];
        int ret;
	float r;
	float fibre_separtion=1.0f;//[cm]
	float fibre_diameter=0.05f;//[cm]
	float abs_return;
	if(-dir->z<=COS_CRIT)
		r=1.0f; //total internal reflection
	else
	{
		if(-dir->z==1.0f)//normal incident
		{		
			r = divide((1.0f-N),(1+N));
			r *= r;//square
		}
		else
		{
			//long and boring calculations of r
			float sinangle_i = sqrtf(1.0f-dir->z*dir->z);
			float sinangle_t = N*sinangle_i;
			float cosangle_t = sqrtf(1.0f-sinangle_t*sinangle_t);
			
			float cossumangle = (-dir->z*cosangle_t) - sinangle_i*sinangle_t;
			float cosdiffangle = (-dir->z*cosangle_t) + sinangle_i*sinangle_t;
			float sinsumangle = sinangle_i*cosangle_t + (-dir->z*sinangle_t);
			float sindiffangle = sinangle_i*cosangle_t - (-dir->z*sinangle_t); 
			
			r = 0.5*sindiffangle*sindiffangle* divide((cosdiffangle*cosdiffangle+cossumangle*cossumangle),(sinsumangle*sinsumangle*cosdiffangle*cosdiffangle));
		
		}
	}
	if(r<1.0f)
	{
		if(rand_MWC_co(&x,&a)<=r)//reflect
			r=1.0f;
		else//transmitt
		{
			//calculate x and y where the photon escapes the medium
			
			r= divide(pos->z,-dir->z);//dir->z must be finite since we have a boundary cross!
			pos->x+=dir->x*r;
			pos->y+=dir->y*r;
			//*t+= divide(r,V); //calculate the time when the photon exits

			r=sqrtf(pos->x*pos->x+pos->y*pos->y);
			
			//check for detection here
			if((fabs((r-fibre_separtion)))<=fibre_diameter)
			{
				//photon detected!
			//	atomic_add( histd + (unsigned int)floor(divide((*t), DT)), 1);
				ret=1;
			}
			else
			{
				ret=2;
			}	
		}
	}
	if(r==1.0f)//reflect (mirror z and dz in reflection plane)
	{
		pos->z *= -1;//mirror the z-coordinate in the z=0 plane, equal to a reflection.
		dir->z *= -1;// do the same to the z direction vector
	}
	ret=0;
	retval[global_id]=ret;

	
	dirx_array[global_id]=dirvar.x; 
	diry_array[global_id]=dirvar.y; 
	dirz_array[global_id]=dirvar.z;
	posx_array[global_id]=posvar.x; 
	posy_array[global_id]=posvar.y; 
	posz_array[global_id]=posvar.z;
}


//EXPECTED RESULTS


//z= costheta
//snells law sintheta1/sintheta2=n2/n1
//theta2 = arcsin(sintheta1*n1/n2)
//       = arcsin(sin(arccos(z1))*n1/n2)
//z2 = cos(theta2)




__kernel void Spin_test1(__global float* dirx_array,__global float* diry_array, __global float* dirz_array, __global unsigned int* xd,__global unsigned int* ad, __global unsigned int* cd, __global float* cost_array, __global float* sint_array, __global float* cosp_array, __global float* sinp_array, const float g ) {

	int global_id = get_global_id(0);
    float3 dirvar;
    float3*dir;
    dirvar.x = dirx_array[global_id];
    dirvar.y =diry_array[global_id]; 
	dirvar.z =dirz_array[global_id];
	dir = &dirvar;
    float P;

	

	unsigned long int x = cd[global_id];
	x=(x<<32) + xd[global_id];
    unsigned int a=ad[global_id];

//START HG

    float cost, sint;	// cosine and sine of the 
						// polar deflection angle theta. 
	float cosp, sinp;	// cosine and sine of the 
						// azimuthal angle psi. 
	float temp;



	//This is more efficient for g!=0 but of course less efficient for g==0

    P = 2.0f* rand_MWC_co(&x,&a) - 1.0f;

	if((g)==0.0f)
		cost = P;
    else{
	   // temp = divide((1.0f-(g)*(g)),(1.0f-(g)+2.0f*(g)*rand_MWC_co(&x,&a)));//Should be close close????!!!!!
        temp = divide((1.0f-(g)*(g)),(1.0f+(g*P)));//Should be close close????!!!!!
	    cost = divide((1.0f+(g)*(g) - temp*temp),(2.0f*(g)));
//        temp = ((1.0f-(g)*(g))/(1.0f+(g*P)));//Should be close close????!!!!!
//	    cost = ((1.0f+(g)*(g) - temp*temp)/(2.0f*(g)));

    }

	sint = sqrtf(1.0f - cost*cost);

	cosp= sincos(2.0f*PI*rand_MWC_co(&x,&a),&sinp);
	
    cost_array[global_id] = cost;
    sint_array[global_id] = sint;
    cosp_array[global_id] = cosp;
    sinp_array[global_id] = sinp;

//END HG

	float tempdir=dir->x;
	temp = sqrtf(1.0f - dir->z*dir->z);

	if(temp==0.0f)// normal incident.
	{
		dir->x = sint*cosp;
		dir->y = sint*sinp;
		dir->z = copysign(cost,dir->z*cost);
	}
	else // regular incident.
	{
		dir->x = divide(sint*(dir->x*dir->z*cosp - dir->y*sinp),temp) + dir->x*cost;
		dir->y = divide(sint*(dir->y*dir->z*cosp + tempdir*sinp),temp) + dir->y*cost;
		dir->z = (-1)*sint*cosp*temp + dir->z*cost;
	}

	//normalisation seems to be required as we are using floats! Otherwise the small numerical error will accumulate
	temp= rsqrtf(dir->x*dir->x+dir->y*dir->y+dir->z*dir->z);
	dir->x=dir->x*temp;
	dir->y=dir->y*temp;
	dir->z=dir->z*temp;
	
	dirx_array[global_id]=dirvar.x; 
	diry_array[global_id]=dirvar.y; 
	dirz_array[global_id]=dirvar.z;
}


//__kernel void newSpintest1(float3* dir, unsigned long* x, unsigned int* a, float3* dir_a, float3* dir_b)
__kernel void newSpin_test1(
__global float* dirx_array,__global float* diry_array, __global float* dirz_array,
__global float* ax_array,__global float* ay_array, __global float* az_array, 
__global float* bx_array,__global float* by_array, __global float* bz_array,  
__global unsigned int* xd,__global unsigned int* ad, __global unsigned int* cd,
__global float* cost_array, __global float* sint_array, __global float* cosp_array, __global float* sinp_array, const float g ) {


	int global_id = get_global_id(0);
    
    float3 dirvar; 
    float3 avar;
    float3 bvar;
    float3*dir;
    float3* dir_a;
    float3* dir_b;
    dirvar.x = dirx_array[global_id];
    dirvar.y =diry_array[global_id]; 
	dirvar.z =dirz_array[global_id];
	dir = &dirvar;
    avar.x = ax_array[global_id];
    avar.y = ay_array[global_id]; 
	avar.z = az_array[global_id];
	dir_a = &avar;
    bvar.x = bx_array[global_id];
    bvar.y = by_array[global_id]; 
	bvar.z = bz_array[global_id];
	dir_b = &bvar;
	
    float cost, sint;	// cosine and sine of the 
						// polar deflection angle theta. 
	float cosp, sinp;	// cosine and sine of the 
						// azimuthal angle psi. 
	float temp;

	float tempdir=dir->x;

	float3 dir_old;
	float3 a_old;
	float3 b_old;
	
	dir_old.x = dir->x;
	dir_old.y = dir->y;
	dir_old.z = dir->z;
	
	a_old.x = dir_a->x;
	a_old.y = dir_a->y;
	a_old.z = dir_a->z;
	
	b_old.x = dir_b->x;
	b_old.y = dir_b->y;
	b_old.z = dir_b->z;

	unsigned long int x = cd[global_id];
	x=(x<<32) + xd[global_id];
    unsigned int a=ad[global_id];

	//This is more efficient for g!=0 but of course less efficient for g==0
	temp = divide((1.0f-(g)*(g)),(1.0f-(g)+2.0f*(g)*rand_MWC_co(&x,&a)));//Should be close close????!!!!!
	cost = divide((1.0f+(g)*(g) - temp*temp),(2.0f*(g)));
	//cost = (1.0f+(G)*(G) - temp*temp)*ONE_OVER_2G;

	if((g)==0.0f)
		cost = 2.0f*rand_MWC_co(&x,&a) -1.0f;


	sint = sqrtf(1.0f - cost*cost);

	cosp= sincos(2.0f*PI*rand_MWC_co(&x,&a),&sinp);

    cost_array[global_id] = cost;
    sint_array[global_id] = sint;
    cosp_array[global_id] = cosp;
    sinp_array[global_id] = sinp;


/*    	
    float P = 2.0f* rand_MWC_co(x,a) - 1.0f;
	if((g)==0.0f)
		cost = P;
    else{
        temp = divide((1.0f-(g)*(g)),(1.0f+(g*P)));//Should be close close????!!!!!
	    cost = divide((1.0f+(g)*(g) - temp*temp),(2.0f*(g)));


    }

	sint = sqrtf(1.0f - cost*cost);

	cosp= sincos(2.0f*PI*rand_MWC_co(x,a),&sinp);
*/
	//temp = sqrtf(1.0f - dir->z*dir->z);
//	dir->x = cost*dir_old.x + (-sint*cosp*a_old.x) + (sint*cosp*b_old.x);
//	dir->y = cost*dir_old.y + (-sint*cosp*a_old.y) + (sint*cosp*b_old.y);
//	dir->z = cost*dir_old.z + (-sint*cosp*a_old.z) + (sint*cosp*b_old.z);

    dir->x = dot_product(cost, -sint*cosp, sint*sinp, dir_old.x, a_old.x, b_old.x) ;
	dir->y = dot_product(cost, -sint*cosp, sint*sinp, dir_old.y, a_old.y, b_old.y) ;
	dir->z = dot_product(cost, -sint*cosp, sint*sinp, dir_old.z, a_old.z, b_old.z) ;
	
//  dir_a->x = sint*dir_old.x + (cost*cosp*a_old.x) + (-cost*sinp*b_old.x);
//	dir_a->y = sint*dir_old.y + (cost*cosp*a_old.y) + (-cost*sinp*b_old.y);
//	dir_a->z = sint*dir_old.z + (cost*cosp*a_old.z) + (-cost*sinp*b_old.z);
	
    dir_a->x = dot_product(sint, cost*cosp, -cost*sinp, dir_old.x, a_old.x, b_old.x) ;
    dir_a->y = dot_product(sint, cost*cosp, -cost*sinp, dir_old.y, a_old.y, b_old.y) ;
    dir_a->z = dot_product(sint, cost*cosp, -cost*sinp, dir_old.z, a_old.z, b_old.z) ;
	
//  dir_b->x = sint*a_old.x + cosp*b_old.x;
//	dir_b->y = sint*a_old.y + cosp*b_old.y;
//	dir_b->z = sint*a_old.z + cosp*b_old.z;

    dir_b->x = dot_product(0, sint, cosp, dir_old.x, a_old.x, b_old.x) ;
    dir_b->y = dot_product(0, sint, cosp, dir_old.y, a_old.y, b_old.y) ;
    dir_b->z = dot_product(0, sint, cosp, dir_old.z, a_old.z, b_old.z) ;



//	if(temp==0.0f)// normal incident.
//	{
//		dir->x = sint*cosp;
//		dir->y = sint*sinp;
//		dir->z = copysign(cost,dir->z*cost);
//	}
//	else // regular incident.
//	{
//		dir->x = divide(sint*(dir->x*dir->z*cosp - dir->y*sinp),temp) + dir->x*cost;
//		dir->y = divide(sint*(dir->y*dir->z*cosp + tempdir*sinp),temp) + dir->y*cost;
//		dir->z = (-1)*sint*cosp*temp + dir->z*cost;
//	}
//
	//normalisation seems to be required as we are using floats! Otherwise the small numerical error will accumulate
	/*
	temp= rsqrtf(dir->x*dir->x+dir->y*dir->y+dir->z*dir->z);
	dir->x=dir->x*temp;
	dir->y=dir->y*temp;
	dir->z=dir->z*temp;
	
	temp= rsqrtf(dir_a->x*dir_a->x+dir_a->y*dir_a->y+dir_a->z*dir_a->z);
	dir_a->x=dir_a->x*temp;
	dir_a->y=dir_a->y*temp;
	dir_a->z=dir_a->z*temp;
	*/
    temp= rsqrtf(dir_b->x*dir_b->x+dir_b->y*dir_b->y+dir_b->z*dir_b->z);
	dir_b->x=dir_b->x*temp;
	dir_b->y=dir_b->y*temp;
	dir_b->z=dir_b->z*temp;
	

	dirx_array[global_id]=dirvar.x; 
	diry_array[global_id]=dirvar.y; 
	dirz_array[global_id]=dirvar.z;
	
    ax_array[global_id]=avar.x; 
	ay_array[global_id]=avar.y; 
	az_array[global_id]=avar.z;
	
    bx_array[global_id]=bvar.x; 
	by_array[global_id]=bvar.y; 
	bz_array[global_id]=bvar.z;
	
}

__kernel void newSpin_test2(
__global float* dirx_array,__global float* diry_array, __global float* dirz_array,
__global float* ax_array,__global float* ay_array, __global float* az_array, 
__global float* bx_array,__global float* by_array, __global float* bz_array,  
__global unsigned int* xd,__global unsigned int* ad, __global unsigned int* cd,
__global float* cost_array, __global float* sint_array, __global float* cosp_array, __global float* sinp_array, const float g ) {


	int global_id = get_global_id(0);
    
    float3 dirvar; 
    float3 avar;
    float3 bvar;
    float3*dir;
    float3* dir_a;
    float3* dir_b;
    dirvar.x = dirx_array[global_id];
    dirvar.y =diry_array[global_id]; 
	dirvar.z =dirz_array[global_id];
	dir = &dirvar;
    avar.x = ax_array[global_id];
    avar.y = ay_array[global_id]; 
	avar.z = az_array[global_id];
	dir_a = &avar;
    bvar.x = bx_array[global_id];
    bvar.y = by_array[global_id]; 
	bvar.z = bz_array[global_id];
	dir_b = &bvar;
	
    float cost, sint;	// cosine and sine of the 
						// polar deflection angle theta. 
	float cosp, sinp;	// cosine and sine of the 
						// azimuthal angle psi. 
	float temp;

	float tempdir=dir->x;

	float3 dir_old;
	float3 a_old;
	float3 b_old;
	
	dir_old.x = dir->x;
	dir_old.y = dir->y;
	dir_old.z = dir->z;
	
	a_old.x = dir_a->x;
	a_old.y = dir_a->y;
	a_old.z = dir_a->z;
	
	b_old.x = dir_b->x;
	b_old.y = dir_b->y;
	b_old.z = dir_b->z;

	/*
	unsigned long int x = cd[global_id];
	x=(x<<32) + xd[global_id];
    unsigned int a=ad[global_id];

	//This is more efficient for g!=0 but of course less efficient for g==0
	temp = divide((1.0f-(g)*(g)),(1.0f-(g)+2.0f*(g)*rand_MWC_co(&x,&a)));//Should be close close????!!!!!
	cost = divide((1.0f+(g)*(g) - temp*temp),(2.0f*(g)));
	//cost = (1.0f+(G)*(G) - temp*temp)*ONE_OVER_2G;

	if((g)==0.0f)
		cost = 2.0f*rand_MWC_co(&x,&a) -1.0f;


	sint = sqrtf(1.0f - cost*cost);

	cosp= sincos(2.0f*PI*rand_MWC_co(&x,&a),&sinp);
	*/

    cost = 0.0f;
    sint = 1.0f;
    cosp = 0.0f;
    sinp = 1.0f;

	cost_array[global_id] = cost;
    sint_array[global_id] = sint;
    cosp_array[global_id] = cosp;
    sinp_array[global_id] = sinp;

    dir->x = dot_product(cost, -sint*cosp, sint*sinp, dir_old.x, a_old.x, b_old.x) ;
	dir->y = dot_product(cost, -sint*cosp, sint*sinp, dir_old.y, a_old.y, b_old.y) ;
	dir->z = dot_product(cost, -sint*cosp, sint*sinp, dir_old.z, a_old.z, b_old.z) ;
	
	
    dir_a->x = dot_product(sint, cost*cosp, -cost*sinp, dir_old.x, a_old.x, b_old.x) ;
    dir_a->y = dot_product(sint, cost*cosp, -cost*sinp, dir_old.y, a_old.y, b_old.y) ;
    dir_a->z = dot_product(sint, cost*cosp, -cost*sinp, dir_old.z, a_old.z, b_old.z) ;
	
    dir_b->x = dot_product(0, sint, cosp, dir_old.x, a_old.x, b_old.x) ;
    dir_b->y = dot_product(0, sint, cosp, dir_old.y, a_old.y, b_old.y) ;
    dir_b->z = dot_product(0, sint, cosp, dir_old.z, a_old.z, b_old.z) ;


	
	//normalisation seems to be required as we are using floats! Otherwise the small numerical error will accumulate
	temp= rsqrtf(dir->x*dir->x+dir->y*dir->y+dir->z*dir->z);
	dir->x=dir->x*temp;
	dir->y=dir->y*temp;
	dir->z=dir->z*temp;
	
	temp= rsqrtf(dir_a->x*dir_a->x+dir_a->y*dir_a->y+dir_a->z*dir_a->z);
	dir_a->x=dir_a->x*temp;
	dir_a->y=dir_a->y*temp;
	dir_a->z=dir_a->z*temp;
	
    temp= rsqrtf(dir_b->x*dir_b->x+dir_b->y*dir_b->y+dir_b->z*dir_b->z);
	dir_b->x=dir_b->x*temp;
	dir_b->y=dir_b->y*temp;
	dir_b->z=dir_b->z*temp;
	

	dirx_array[global_id]=dirvar.x; 
	diry_array[global_id]=dirvar.y; 
	dirz_array[global_id]=dirvar.z;
	
    ax_array[global_id]=avar.x; 
	ay_array[global_id]=avar.y; 
	az_array[global_id]=avar.z;
	
    bx_array[global_id]=bvar.x; 
	by_array[global_id]=bvar.y; 
	bz_array[global_id]=bvar.z;

}


