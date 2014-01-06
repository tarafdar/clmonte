
#include <CL/cl.h>
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>



float copysign(float x, float y){

	if (y >= 0)
		return abs(x);
	else
		return -abs(x);
}


float rand_MWC_coh(unsigned long long* x, unsigned int* a)
{
		//this implementation seems to be faster.
		*x=(*x&0xffffffffull)*(*a)+(*x>>32);
		return((float)((unsigned int)(*x&0xffffffffull))/(UINT_MAX));

}//end __device__ rand_MWC_coh


float rand_MWC_och(unsigned long long* x, unsigned int* a)
{
		*x=(*x&0xffffffffull)*(*a)+(*x>>32);
		return(1.0f-(float)((unsigned int)(*x&0xffffffffull))/(UINT_MAX));//return(1.0f-(float)((unsigned int)(*x&0xffffffffull))/(UINT_MAX));
}//end __device__ rand_MWC_och

float rsqrtf(float x){
	return (float)1.0/(sqrtf(x));

}

void Spinh(float3* dir, float* g, unsigned long long* x, unsigned int* a, float* cost_ret, float* sint_ret, float* cosp_ret, float* sinp_ret)
{
	float cost, sint;	// cosine and sine of the 
						// polar deflection angle theta. 
	float cosp, sinp;	// cosine and sine of the 
						// azimuthal angle psi. 
	float temp;

	float tempdir=dir->x;

	if((*g)==0.0f)
		cost = 2.0f*rand_MWC_coh(x,a) -1.0f;//Should be close close??!!!!!
	else
	{
		temp = (1.0f-(*g)*(*g))/(1.0f-(*g)+2.0f*(*g)*rand_MWC_coh(x,a));//Should be close close????!!!!!
		cost = (1.0f+(*g)*(*g) - temp*temp)/(2.0f*(*g));
	}
	sint = sqrtf(1.0f - cost*cost);

	temp = 2.0f*PI*rand_MWC_coh(x,a); // spin psi [0-2*PI)
	cosp = cosf(temp);
	sinp = sinf(temp);
	
	temp = sqrtf(1.0f - dir->z*dir->z);
	
    *(cost_ret) = cost;
    *(sint_ret) = sint;
    *(cosp_ret) = cosp;
    *(sinp_ret) = sinp;
    
	if(temp==0.0f)// normal incident.
	{
		dir->x = sint*cosp;
		dir->y = sint*sinp;
		dir->z = copysign(cost,dir->z*cost);
		
	}
	else // regular incident.
	{
		dir->x = sint*(dir->x*dir->z*cosp - dir->y*sinp)/temp + dir->x*cost;
		dir->y = sint*(dir->y*dir->z*cosp + tempdir*sinp)/temp + dir->y*cost;
		dir->z = -sint*cosp*temp + dir->z*cost;
	}

	//normalisation seems to be required as we are using floats! Otherwise the small numerical error will accumulate
	temp=rsqrtf(dir->x*dir->x+dir->y*dir->y+dir->z*dir->z);
	dir->x=dir->x*temp;
	dir->y=dir->y*temp;
	dir->z=dir->z*temp;
	

}

void Spinh2(float3* dir, float* g, unsigned long long* x, unsigned int* a, float* cost_ret, float* sint_ret, float* cosp_ret, float* sinp_ret)
{
	float cost, sint;	// cosine and sine of the 
						// polar deflection angle theta. 
	float cosp, sinp;	// cosine and sine of the 
						// azimuthal angle psi. 
	float temp;

	float tempdir=dir->x;

    float P = 2.0f* rand_MWC_coh(x,a) - 1.0f;
	if((*g)==0.0f)
		cost = P;//Should be close close??!!!!!
	else
	{
        temp = ((1.0f-(*g)*(*g))/(1.0f+(*g*P)));//Should be close close????!!!!!
	    cost = ((1.0f+(*g)*(*g) - temp*temp)/(2.0f*(*g)));
	}
	sint = sqrtf(1.0f - cost*cost);

	temp = 2.0f*PI*rand_MWC_coh(x,a); // spin psi [0-2*PI)
	cosp = cosf(temp);
	sinp = sinf(temp);
	
	temp = sqrtf(1.0f - dir->z*dir->z);
	
    *(cost_ret) = cost;
    *(sint_ret) = sint;
    *(cosp_ret) = cosp;
    *(sinp_ret) = sinp;
    
	if(temp==0.0f)// normal incident.
	{
		dir->x = sint*cosp;
		dir->y = sint*sinp;
		dir->z = copysign(cost,dir->z*cost);
		
	}
	else // regular incident.
	{
		dir->x = sint*(dir->x*dir->z*cosp - dir->y*sinp)/temp + dir->x*cost;
		dir->y = sint*(dir->y*dir->z*cosp + tempdir*sinp)/temp + dir->y*cost;
		dir->z = -sint*cosp*temp + dir->z*cost;
	}

	//normalisation seems to be required as we are using floats! Otherwise the small numerical error will accumulate
	temp=rsqrtf(dir->x*dir->x+dir->y*dir->y+dir->z*dir->z);
	dir->x=dir->x*temp;
	dir->y=dir->y*temp;
	dir->z=dir->z*temp;
	

}


void Spin_testh(float3* dir, float* cost_array, float* sint_array, float* cosp_array, float* sinp_array, float G)
{
    int i;	

    float g = G;
	unsigned long long int x;
	unsigned int a;


	for(i=0;i<NUM_THREADS;i++)
	{
		
		x=ctest[i];
		x=(x<<32)+xtest[i];

		a=atest[i];



		//Spinh(&dir[i],&g,&x,&a, &cost_array[i], &sint_array[i], &cosp_array[i], &sinp_array[i]);
		Spinh2(&dir[i],&g,&x,&a, &cost_array[i], &sint_array[i], &cosp_array[i], &sinp_array[i]);


	}//end for block
	
}//end MCh

