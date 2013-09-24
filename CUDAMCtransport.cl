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
#define NUM_THREADS 26880
#define NUMSTEPS_GPU 500000
#define NUMSTEPS_CPU 500000
#define PI 3.14159265f

#define TMAX 2000.0f //[ps] Maximum time of flight
#define DT 10.0f //[ps] Time binning resolution
#define TEMP 201 //ceil(TMAX/DT), precalculated to avoid dynamic memory allocation (fulhack)

#define MUS_MAX 90.0f	//[1/cm]
#define V 0.0214f		//[cm/ps] (c=0.03 [cm/ps] v=c/n) here n=1.4
#define COS_CRIT 0.6999f	//the critical angle for total internal reflection at the border cos_crit=sqrt(1-(nt/ni)^2)
#define G 0.9f	
#define N 1.4f

#define ONE_OVER_2G				( 1.0f / (2.0f * G) )
#define ONE_MINUS_N_OVER_1_PLUS_N ( (1.0f - N) / (1.0f + N) )
#define ONE_OVER_V				( 1.0f / V )
#define ONE_OVER_DT				( 1.0f / DT )
#define ONE_OVER_MUS_MAX			( 1.0f / MUS_MAX )

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



void Spin(float3* dir, unsigned long* x,//unsigned int* c,
		       unsigned int* a)
{
	float cost, sint;	// cosine and sine of the 
						// polar deflection angle theta. 
	float cosp, sinp;	// cosine and sine of the 
						// azimuthal angle psi. 
	float temp;

	float tempdir=dir->x;


	//This is more efficient for g!=0 but of course less efficient for g==0
	temp = divide((1.0f-(G)*(G)),(1.0f-(G)+2.0f*(G)*rand_MWC_co(x,a)));//Should be close close????!!!!!
	//cost = divide((1.0f+(G)*(G) - temp*temp),(2.0f*(G)));
	cost = (1.0f+(G)*(G) - temp*temp)*ONE_OVER_2G;

	//temp = (1.0f-(*g)*(*g))/(1.0f-(*g)+2.0f*(*g)*rand_MWC_co(x,a));//Should be close close????!!!!!
	//cost = (1.0f+(*g)*(*g) - temp*temp)(2.0f*(*g)));


	if((G)==0.0f)
		cost = 2.0f*rand_MWC_co(x,a) -1.0f;


	sint = sqrtf(1.0f - cost*cost);

	cosp= sincos(2.0f*PI*rand_MWC_co(x,a),&sinp);
	

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
	
}


unsigned int Reflect(float3* dir, float3* pos, float* t,  unsigned long* x, unsigned int* a, __global unsigned int* histd)
{
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
			//r = divide((1.0f-N),(1+N));
			r = ONE_MINUS_N_OVER_1_PLUS_N;
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
		if(rand_MWC_co(x/*,c*/,a)<=r)//reflect
			r=1.0f;
		else//transmitt
		{
			//calculate x and y where the photon escapes the medium
			
			r= divide(pos->z,-dir->z);//dir->z must be finite since we have a boundary cross!
			pos->x+=dir->x*r;
			pos->y+=dir->y*r;
			//*t+= divide(r,V); //calculate the time when the photon exits
			*t+= r*ONE_OVER_V;

			r=sqrtf(pos->x*pos->x+pos->y*pos->y);
			
			//check for detection here
			if((fabs((r-fibre_separtion)))<=fibre_diameter)
			{
				//photon detected!
				//atomic_add( histd + (unsigned int)(floor((divide((*t),DT)) , 1)));//&histd[(unsigned int)floorf(native_divide((t*),DT))],(unsigned int)1);
				//unsigned int offset;
				//offset= (unsigned int)(floor(native_divide((*t),DT)))
				//atomic_add( histd + (unsigned int)floor(divide((*t), DT)), 1);
				atomic_add( histd + (unsigned int)floor(*t*ONE_OVER_DT), 1);
				return 1;
			}
			else
			{
				return 2;
			}	
		}
	}
	if(r==1.0f)//reflect (mirror z and dz in reflection plane)
	{
		pos->z *= -1;//mirror the z-coordinate in the z=0 plane, equal to a reflection.
		dir->z *= -1;// do the same to the z direction vector
	}
	return 0;
}





__kernel void MCd(__global unsigned int* xd,__global unsigned int* cd,__global unsigned int* ad,__global unsigned int* numd,__global unsigned int* histd)
{
	int global_id= get_global_id(0);
    //for loops
    unsigned int ii=0;

    //First element processed by the block
    //int begin=NUM_THREADS_PER_BLOCK*bx;

    unsigned long int x=cd[global_id];
	
	x=(x<<32)+xd[global_id];


    unsigned int a=ad[global_id];


	float3 pos; //float triplet to store the position
	float3 dir; //float triplet to store the direction
	float t;	//float to store the time of flight
	float s;	//step length
	
	unsigned int num_det_photons=0;
	unsigned int flag=0;

	LaunchPhoton(&pos, &dir, &t);//Launch the photon
	
	for(ii=0;ii<NUMSTEPS_GPU;ii++) //this is the main while loop
	{
		//num_det_photons++;
		//s = divide(- logf(rand_MWC_oc(&x,&a)), MUS_MAX);//sample step length 
		s = (- logf(rand_MWC_oc(&x,&a))*ONE_OVER_MUS_MAX);

		//Perform boundary crossing check here
		if((pos.z+dir.z*s)<=0)//photon crosses boundary within the next step
		{
			flag=Reflect(&dir,&pos,&t,&x,&a,histd);
		}
		
		//Move (we can move the photons that have been terminated above since it improves our performance and does not affect our results)
		pos.x += s*dir.x;
		pos.y += s*dir.y;
		pos.z += s*dir.z;
		t += divide(s,V); 

		Spin(&dir,&x,&a);

		if(t >= TMAX || flag>=1)//Kill photon and launch a new one
		{
			num_det_photons++;
			flag=0;
			LaunchPhoton(&pos, &dir, &t);//Launch the photon
		}
		

	}
	
	//end main for loop!
	

	//barrier(CLK_GLOBAL_MEM_FENCE);//necessary?

	numd[global_id]/*[begin+tx]*/=num_det_photons; 

}//end MCd

