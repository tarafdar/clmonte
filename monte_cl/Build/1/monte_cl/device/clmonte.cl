
#include <misc.h>

#define PI 3.1415926F
#define CHANCE 0.1F
#define WEIGHT 1E-4F
#define WEIGHT_SCALE 65536//12000000

// TODO: Implementation specific to GPU for taylor series expansion approximation?
// TODO: Can FPGA give exact values?
#define COSNINETYDEG 1.0E-6F
#define COSZERO (1.0F - 1.0E-6F)

#define RAND_MAX 2147483647.0F

// From glibc
int rand_r (UINT32 *seed)
{
	UINT32 next = *seed;
	int result;

	next *= 1103515245;
	next += 12345;
	result = (UINT32) (next / 65536) % 2048;

	next *= 1103515245;
	next += 12345;
	result <<= 10;
	result ^= (UINT32) (next / 65536) % 1024;

	next *= 1103515245;
	next += 12345;
	result <<= 10;
	result ^= (UINT32) (next / 65536) % 1024;

	*seed = next;

	return result;
}

typedef struct Packet
{
	Point pos;
	Point dir;

	float step;
	float weight;

	Point a;
	Point b;

	TetraID current_tid;
} Packet;

float dot_product(Point a, Point b)
{
	// TODO: dot Altera OpenCL function more optimal?
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

Point cross_product(Point a, Point b)
{
	Point cp;
	cp.x = a.y*b.z - a.z*b.y;
	cp.y = a.z*b.x - a.x*b.z;
	cp.z = a.x*b.y - a.y*b.x;
	return cp;
}

Point normalize_(Point a)
{
	float mag = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
	a.x /= mag;
	a.y /= mag;
	a.z /= mag;
	return a;
}

float unit_rand(UINT32* seed)
{
	float result = rand_r(seed)/RAND_MAX;
	return result;
}

float draw_step(UINT32* seed)
{
	float r = unit_rand(seed);
	return -log(r);
}

Packet update_spin_helpers(Packet pkt)
{
	// Spin Vectors
	// vector a = (vector d) cross (positive z axis unit vector) and normalize it 
	Point crossProduct;
	if(pkt.dir.x == 0 && pkt.dir.y == 0) {
		crossProduct.x = 1;
		crossProduct.y = 0;
		crossProduct.z = 0;
	}
	else {
		crossProduct.x = 0;
		crossProduct.y = 0;
		crossProduct.z = 1;
	}

	crossProduct = cross_product(pkt.dir, crossProduct);
	pkt.a = normalize_(crossProduct);

	// vector b = (vector d) cross (vector a)
	pkt.b = cross_product(pkt.dir,pkt.a);

	return pkt;
}

Packet launch(Source source, UINT32* seed)
{
	Packet pkt;
	pkt.step = draw_step(seed);
	pkt.weight = 1.0F;
	pkt.current_tid = source.tid;

	// Position
	pkt.pos.x = source.pos.x;
	pkt.pos.y = source.pos.y;
	pkt.pos.z = source.pos.z;

	// Direction	
	float sint, sinp, cost, cosp;
	float theta = PI*unit_rand(seed);
	float phi = 2.0F*PI*unit_rand(seed);
	sint = sincos(theta, &cost);
	sinp = sincos(phi, &cosp);

	pkt.dir.x = sinp*cost;
	pkt.dir.y = sinp*sint;
	pkt.dir.z = cosp;
printf("LAUNCHING\n");
	pkt = update_spin_helpers(pkt);
	return pkt;
}

// Project a onto b where |b| = 1
float project_on_unit(Point a, float b[])
{
	Point temp;
	temp.x = b[0];
	temp.y = b[1];
	temp.z = b[2];
	return dot_product(a, temp);
}

float distance_to_plane(Point p, float plane[])
{
	return project_on_unit(p, plane) - plane[3];
}

Packet reflect_packet(Packet pkt, float ca1, float normal[])
{
	// Reflected direction = original_direction - 2(original_direction dot normal)*normal
	// Normal is oriented inwards hence addition
	pkt.dir.x += 2*ca1*normal[0];
	pkt.dir.y += 2*ca1*normal[1];
	pkt.dir.z += 2*ca1*normal[2];
	return pkt;
}

Packet move_packet(Packet pkt, float distance)
{
	pkt.pos.x += distance*pkt.dir.x;
	pkt.pos.y += distance*pkt.dir.y;
	pkt.pos.z += distance*pkt.dir.z;
	return pkt;
}

Point scalar_multiply(Point point, float scalar)
{
	point.x *= scalar;
	point.y *= scalar;
	point.z *= scalar;
	return point;
}

Point add_points(Point a, Point b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}
Point sub_points(Point a, Point b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}

float GetCosCrit(float ni, float nt)
{
	float sin_crit = native_divide(nt,ni);
	float cos_crit = sqrt(1.0F - sin_crit*sin_crit);
	return cos_crit;
}

__kernel void entry(__global UINT32* photon_pool, __global const Source* source, __global const Material* materials, __global const Tetra* mesh, __global const UINT32* num_tetras, __global const UINT32* num_mats, __global UINT32* absorption, __global UINT32* transmittance) {

	UINT32 seed_t = 80;
	UINT32* seed = &seed_t;
  unsigned thread_id = get_global_id(0);

    printf("Thread #%u: Hello from Altera's OpenCL Compiler!\n", thread_id);
	printf("Args: nT %u, nM %u, init_tid %u, pool_size %u\n", *num_tetras, *num_mats, source->tid, *photon_pool);

	Source launcher = *source;

	Packet pkt = launch(launcher, seed);
	
	bool dead = false;
	while(!dead) {
		if(pkt.current_tid == 0 || pkt.weight == 0) {
			if(*photon_pool == 0) {
				printf("PHOTON DIED 1\n");
				dead = true;
				break;
			}
			else if (atomic_sub(photon_pool, 1) > 0) {
				printf("LAUNCH 1: %u\n", *photon_pool);
				// This pkt is terminated, launch another
				pkt = launch(launcher, seed);
			}
			
		}

		//printf("Location ID: %u\n", pkt.current_tid);
		if(pkt.step == 0) {
			pkt.step = draw_step(seed);
		}

		float cosdn[4], dis[4];
		Tetra tetra = mesh[pkt.current_tid];
		Material mat = materials[tetra.matID];
      
		// = |dir|cos0
		cosdn[0] = project_on_unit(pkt.dir, tetra.face[0]);
		cosdn[1] = project_on_unit(pkt.dir, tetra.face[1]);
		cosdn[2] = project_on_unit(pkt.dir, tetra.face[2]);
		cosdn[3] = project_on_unit(pkt.dir, tetra.face[3]);

		// (p*n) - d(plane constant) = |p-p0|cos0
		dis[0] = distance_to_plane(pkt.pos, tetra.face[0]);
		dis[1] = distance_to_plane(pkt.pos, tetra.face[1]);
		dis[2] = distance_to_plane(pkt.pos, tetra.face[2]);
		dis[3] = distance_to_plane(pkt.pos, tetra.face[3]);
	
		
		// if the packet is moving away from a face, its distance to that face's intersection is infinity.
		// Distance divided by speed (direction) = time to intersect.
		float time[4];
		time[0] = cosdn[0]>=0 ? FLT_MAX : -native_divide(dis[0], cosdn[0]);
		time[1] = cosdn[1]>=0 ? FLT_MAX : -native_divide(dis[1], cosdn[1]);
		time[2] = cosdn[2]>=0 ? FLT_MAX : -native_divide(dis[2], cosdn[2]);
		time[3] = cosdn[3]>=0 ? FLT_MAX : -native_divide(dis[3], cosdn[3]);

		// Select the minimum time to intersection
		// Cascaded comparison tree
		UINT32 localMinIndex1, localMinIndex2, minIndex;
		localMinIndex1 = time[0] < time[1] ? 0 : 1;
		localMinIndex2 = time[2] < time[3] ? 2 : 3;
		minIndex = time[localMinIndex1] < time[localMinIndex2] ? localMinIndex1 : localMinIndex2;
  
		TetraID nextTetraID = tetra.adjTetras[minIndex];

		float travel_distance = pkt.step * mat.rmu_as;
		float r = unit_rand(seed);
		if(travel_distance < dis[minIndex])
		{
			// Remain in same tetra
			pkt.step = 0.0F;
			pkt = move_packet(pkt, travel_distance);
        
			float dw = pkt.weight*mat.absfrac;
			pkt.weight -= dw;
			atomic_add(&(absorption[pkt.current_tid]), (UINT32)(dw*WEIGHT_SCALE));
	
			if (pkt.weight < WEIGHT)
			{
				// Cons. of Energy guaranteed?
				// This pkt survives the roulette.
				if (pkt.weight != 0.0F && r < CHANCE) {
					pkt.weight *= (1.0F / CHANCE);
				}
				else if(*photon_pool == 0) {
					printf("PHOTON DIED 2\n");
					dead = true;
					break;
				}
				else if (atomic_sub(photon_pool, 1) > 0) {
					printf("LAUNCH 2: %u\n", *photon_pool);
					// This pkt is terminated, launch another
					pkt = launch(launcher, seed);
				}
			}
			else
			{
				// Spin
				float cost, sint, cosp, sinp;
				Point last_dir, last_a, last_b;
				//rand = FP_TWO * rand_MWC_co(&rnd_x, &rnd_a) - FP_ONE;
				r = 2.0F * r - 1.0F;
				cost = mat.HGCoeff1 - native_divide(mat.HGCoeff2, (1-mat.g*r)*(1-mat.g*r));
				sint = sqrt(1.0F - cost * cost);

				// spin psi 0-2pi.
				r = unit_rand(seed);
				sinp = sincos(2.0F * PI * r, &cosp);

				float stcp = sint * cosp;
				float stsp = sint * sinp;
				float ctcp = cost * cosp;
				float ctsp = cost * sinp;

				last_dir = pkt.dir;
				last_a = pkt.a;
				last_b = pkt.b;
				//pkt.d = cost*last_d - stcp*last_a + stsp*last_b;
				pkt.dir = add_points(sub_points(scalar_multiply(last_dir, cost), scalar_multiply(last_a, stcp)), scalar_multiply(last_b,stsp));
				//pkt.a = sint*last_d + ctcp*last_a - ctsp*last_b;
				pkt.a = sub_points(add_points(scalar_multiply(last_dir,sint), scalar_multiply(last_a, ctcp)), scalar_multiply(last_b,ctsp));
				//pkt.b = sinp*last_a + cosp*last_b;
				pkt.b = add_points(scalar_multiply(last_a,sinp), scalar_multiply(last_b,cosp));
			}
		}
		else
		{
			// Intersection
			pkt.step -= dis[minIndex] * mat.mu_as;
			pkt = move_packet(pkt, dis[minIndex]);

			Tetra nextTetra = mesh[nextTetraID];
			float ni, nt; //refractive indices
			ni = mat.n;
			nt = materials[nextTetra.matID].n;
			if (ni==nt) {
				// Direction doesn't change
				pkt.current_tid = nextTetraID;
				if (nextTetraID == 0) {
					//UINT32 trans_index = (pkt.current_tid - 1) * 4 + minIndex;
					//float add_2 = (UINT32)(pkt.weight * WEIGHT_SCALE);
					//atomic_add(&(transmittance[trans_index]), add_2);
					pkt.weight = 0.0F;
				}
			}
			else {
				float crit_cos = nt < ni ? 0.0F : GetCosCrit(ni,nt);

				float ca1 = -project_on_unit(pkt.dir, tetra.face[minIndex]); //ca1 is cos(theta incidence)

				// TODO: any special reason for doing in this order?
				// TODO: Waste of compute to do sqrt before check
				float sa1 = sqrt(1.0F-ca1*ca1); //sa1 is sin(theta incidence)
				if (ca1 > COSZERO) sa1 = 0.0F; 

				if (ca1 <= crit_cos) {
					//total internal reflection occurs
					//pkt = reflect_packet(pkt, ca1, normal);
					pkt = reflect_packet(pkt, ca1, tetra.face[minIndex]);
        			}
			        else
			        {
					// Calculate Fresnel
					// Rs = [(n1 cos(theta_i) - n2 cos(theta_t))/((n1 cos(theta_i) + n2 cos(theta_t)))] ^ 2
					// Rp = [(n1 cos(theta_t) - n2 cos(theta_i))/((n1 cos(theta_t) + n2 cos(theta_i)))] ^ 2

					// TODO: Already done by get_crit_cos?
					float ni_nt = native_divide(ni, nt);
					float sa2 = min(ni_nt * sa1, 1.0F); //sa2 is sin(theta transmit)
					float ca2 = sqrt(1.0F - sa2*sa2); //ca2 is cos(theta transmit)
					float Rs = native_divide(ni*ca1 - nt*ca2, ni*ca1 + nt*ca2);
					Rs *= Rs;
					float Rp = native_divide(ni*ca2 - nt*ca1, ni*ca2 + nt*ca1);
					Rp *= Rp;

					// TODO: Why in this order? Wasted above computation					
					float rFresnel = (Rs + Rp)/2;
					if (ca1 < COSNINETYDEG || sa2 == 1.0F) rFresnel = 1.0F;

					// Computation simplification, reflect or refract entirely, no partial components
					// Monte Carlo approximation will correct for this
					if (rFresnel < r) //refract
					{
						// refracted direction = n1/n2*original_direction - [n1/n2*cos(theta incidence) + sqrt(1-sin^2(theta transmitt))]*normal
						//pkt.dir = ni_nt*(pkt.dir) - fma(ni_nt,-ca1,ca2)*(float4)(normal[0],normal[1],normal[2],0);
						float4 temp = fma(ni_nt,-ca1,ca2)*(float4)(tetra.face[minIndex][0],tetra.face[minIndex][1],tetra.face[minIndex][2],0);
						Point temp2;
						temp2.x = temp.x;
						temp2.y = temp.y;
						temp2.z = temp.z;
						pkt.dir = sub_points(scalar_multiply(pkt.dir, ni_nt), temp2);

						if (nextTetraID == 0) {
      
							//store surface fluence TODO: find more efficient way of doing this!
							atomic_add(&(transmittance[(pkt.current_tid - 1) * 4 + minIndex]), (UINT32)(pkt.weight * WEIGHT_SCALE));
							// Kill the packet.
							pkt.weight = 0.0F;
						}
						pkt.current_tid = nextTetraID;
					}
					else //reflect
					{
						//reflected direction = original_direction - 2(original_direction dot normal)*normal
						pkt = reflect_packet(pkt, ca1, tetra.face[minIndex]);
					}
				}

				// Spin Stuff
				pkt = update_spin_helpers(pkt);
			}
		}
	}

	printf("Simulation Complete!\n");
}

