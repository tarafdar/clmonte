#include "defines.h"

float dot_product (float x1, float y1, float z1, float x2, float y2, float z2) {
    return  x1*x2 + y1*y2 + z1*z2;
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


void LaunchPhoton(float3* pos, float3* dir, float* t, float3* dir_a, float3* dir_b)
{
	pos->x=0.0f;
	pos->y=0.0f;
	pos->z=0.0f;

	dir->x=0.0f;
	dir->y=0.0f;
	dir->z=1.0f;

	dir_a->x=1.0f;
	dir_a->y=0.0f;
	dir_a->z=0.0f;

	dir_b->x=0.0f;
	dir_b->y=1.0f;
	dir_b->z=0.0f;

	*t=0.0f;
}



void Spin(float3* dir, unsigned long* x, unsigned int* a, float3* dir_a, float3* dir_b)
{
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


	//This is more efficient for g!=0 but of course less efficient for g==0
	//cost = (1.0f+(G)*(G) - temp*temp)*ONE_OVER_2G;

	if((G)==0.0f)
		cost = 2.0f*rand_MWC_co(x,a) -1.0f;
	else
    {    
	    temp = divide((1.0f-(G)*(G)),(1.0f-(G)+2.0f*(G)*rand_MWC_co(x,a)));//Should be close close????!!!!!
        cost = divide((1.0f+(G)*(G) - temp*temp),(2.0f*(G)));
    }

	sint = sqrtf(1.0f - cost*cost);

	cosp= sincos(2.0f*PI*rand_MWC_co(x,a),&sinp);

    dir->x = dot_product(cost, -sint*cosp, sint*sinp, dir_old.x, a_old.x, b_old.x) ;
	dir->y = dot_product(cost, -sint*cosp, sint*sinp, dir_old.y, a_old.y, b_old.y) ;
	dir->z = dot_product(cost, -sint*cosp, sint*sinp, dir_old.z, a_old.z, b_old.z) ;
	
    dir_a->x = dot_product(sint, cost*cosp, -cost*sinp, dir_old.x, a_old.x, b_old.x) ;
    dir_a->y = dot_product(sint, cost*cosp, -cost*sinp, dir_old.y, a_old.y, b_old.y) ;
    dir_a->z = dot_product(sint, cost*cosp, -cost*sinp, dir_old.z, a_old.z, b_old.z) ;
	
    dir_b->x = dot_product(0, sinp, cosp, dir_old.x, a_old.x, b_old.x) ;
    dir_b->y = dot_product(0, sinp, cosp, dir_old.y, a_old.y, b_old.y) ;
    dir_b->z = dot_product(0, sinp, cosp, dir_old.z, a_old.z, b_old.z) ;

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
}

void cross_product(float3* a, float3* b, float3* c){
    c->x = a->y*b->z - a->z*b->y;
    c->y = a->z*b->x - a->x*b->z;
    c->z = a->x*b->y - a->y*b->x;
}

void cross_product_unit(float3* a, int z, float3* c){
    if(z == 1) { //0,0,1
        c->x = a->y;
        c->y = 0 - a->x;
        c->z = 0;
    }
    else{ //cross with 0,1,0
        c->x = 0 - a->z;
        c->y = 0;
        c->z = a->x;
    }

}

unsigned int Reflect(float3* dir, float3* pos, float* t,  unsigned long* x, unsigned int* a, __global unsigned int* histd, float3* dir_a, float3* dir_b)
{
	float r;
	float fibre_separtion=1.0f;//[cm]
	float fibre_diameter=0.05f;//[cm]
	float abs_return;
    
    float temp;
    float3 z_unit;
    z_unit.x = 0;
    z_unit.y = 0;
    z_unit.z = 1;
	
    
    if(-dir->z<=COS_CRIT)
		r=1.0f; //total internal reflection
	else
	{
		if(-dir->z==1.0f)//normal incident
		{		
			r = divide((1.0f-N),(1+N));
			//r = ONE_MINUS_N_OVER_1_PLUS_N_SQUARED;
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
			
            *t+= divide(r,V); //calculate the time when the photon exits
			//*t+= r*ONE_OVER_V;

			r=sqrtf(pos->x*pos->x+pos->y*pos->y);
	        	
             
        
#ifdef SPATIAL_HISTOGRAM
        
        
            if(r <= MAXR){
                atomic_add(histd + (unsigned int)floor(r/DR)*TEMP + (unsigned int)floor(divide((*t), DT)), 1);
                return 1;   
            }
            else
                return 2;
        

#else	
			//check for detection here
			if((fabs((r-fibre_separtion)))<=fibre_diameter)
			{
				//photon detected!
				atomic_add( histd + (unsigned int)floor(divide((*t), DT)), 1);
				return 1;
			}
			else
			{
				return 2;
			}
#endif	
		}
	}
	if(r==1.0f)//reflect (mirror z and dz in reflection plane)
	{
		pos->z *= -1;//mirror the z-coordinate in the z=0 plane, equal to a reflection.
		dir->z *= -1;// do the same to the z direction vector
        
        if(dir->x != 0 && dir->y != 0 && dir->z != 1) 
            cross_product_unit(dir, 1, dir_a);
        else
            cross_product_unit(dir, 0, dir_a); 
	    
        temp= rsqrtf(dir_a->x*dir_a->x+dir_a->y*dir_a->y+dir_a->z*dir_a->z);
	    dir_a->x=dir_a->x*temp;
	    dir_a->y=dir_a->y*temp;
	    dir_a->z=dir_a->z*temp;
        
        cross_product(dir, dir_a, dir_b); 
        temp= rsqrtf(dir_b->x*dir_b->x+dir_b->y*dir_b->y+dir_b->z*dir_b->z);
	    dir_b->x=dir_b->x*temp;
	    dir_b->y=dir_b->y*temp;
	    dir_b->z=dir_b->z*temp;
	}
	return 0;
}





__kernel void MCd(__global unsigned int* xd,__global unsigned int* cd,__global unsigned int* ad,
__global unsigned int* histd
#ifdef EVENT_LOGGING
,__global unsigned int* numd,
 __global unsigned int* scatteringEvents
#endif
)
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
	float3 dir_a; //float triplet to store the position
	float3 dir_b; //float triplet to store the position
	float t;	//float to store the time of flight
	float s;	//step length
#ifdef EVENT_LOGGING	
	unsigned int num_det_photons=0;
    unsigned int num_scatters=0;
#endif
	unsigned int flag=0;
	LaunchPhoton(&pos, &dir, &t, &dir_a, &dir_b);//Launch the photon
	
	for(ii=0;ii<NUMSTEPS_GPU;ii++) //this is the main while loop
	{
		//num_det_photons++;
		s = divide(- logf(rand_MWC_oc(&x,&a)), MUS_MAX);//sample step length 
		//s = (- logf(rand_MWC_oc(&x,&a))*ONE_OVER_MUS_MAX);

		//Perform boundary crossing check here
		if((pos.z+dir.z*s)<=0)//photon crosses boundary within the next step
		{
			flag=Reflect(&dir,&pos,&t,&x,&a,histd, &dir_a, &dir_b);
		}
		
		//Move (we can move the photons that have been terminated above since it improves our performance and does not affect our results)
		pos.x += s*dir.x;
		pos.y += s*dir.y;
		pos.z += s*dir.z;
		t += divide(s,V); 

		Spin(&dir,&x,&a, &dir_a, &dir_b);
#ifdef EVENT_LOGGING
        num_scatters++;
#endif		
        if(t >= TMAX || flag>=1)//Kill photon and launch a new one
		{
#ifdef EVENT_LOGGING
			num_det_photons++;
#endif
			flag=0;
			LaunchPhoton(&pos, &dir, &t, &dir_a, &dir_b);//Launch the photon
		}
		

	}
	
	//end main for loop!
	

	//barrier(CLK_GLOBAL_MEM_FENCE);//necessary?
#ifdef EVENT_LOGGING    
	numd[global_id]/*[begin+tx]*/=num_det_photons; 
    scatteringEvents[global_id]=num_scatters;
#endif
}//end MCd

