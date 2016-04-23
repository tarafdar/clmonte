
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

// Atomic is not required with only one pipeline
// #define ATOMIC_REQ
#ifdef ATOMIC_REQ
#define ATOMIC_ADD(a,b) atomic_add(&(a), (b))
#else
#define ATOMIC_ADD(a,b) (a = ((a) + (b)))
#endif

typedef struct Spinner
{
	Point a;
	Point b;
	Point dir;
} Spinner;

typedef struct RNG_R
{
	UINT32 seed;
	UINT32 rand;
} RNG_R;

typedef struct RNG 
{
	UINT32 seed;
	float rand;
} RNG;

// From glibc
RNG_R rand_r (UINT32 seed)
{
	UINT32 next = seed;
	UINT32 result;

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

	RNG_R new_rng;
	new_rng.seed = next;
	new_rng.rand = result;
	return new_rng;
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

RNG unit_rand(UINT32 seed)
{
	RNG_R rng_r = rand_r(seed);

	RNG rng;
	rng.seed = rng_r.seed;
	rng.rand = rng_r.rand / RAND_MAX;
	return rng;
}

Spinner update_spin_helpers(Point dir)
{
	// Spin Vectors
	// vector a = (vector d) cross (positive z axis unit vector) and normalize it 
	Point crossProduct;
	if(dir.x == 0 && dir.y == 0) {
		crossProduct.x = 1;
		crossProduct.y = 0;
		crossProduct.z = 0;
	}
	else {
		crossProduct.x = 0;
		crossProduct.y = 0;
		crossProduct.z = 1;
	}

	crossProduct = cross_product(dir, crossProduct);

	Spinner spinner;
	spinner.a = normalize_(crossProduct);

	// vector b = (vector d) cross (vector a)
	spinner.b = cross_product(dir, spinner.a);

	return spinner;
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

Point reflect_packet(Packet pkt, float ca1, float normal[])
{
	// Reflected direction = original_direction - 2(original_direction dot normal)*normal
	// Normal is oriented inwards hence addition
	Point p;
	p.x += 2*ca1*normal[0];
	p.y += 2*ca1*normal[1];
	p.z += 2*ca1*normal[2];
	return p;
}

Point move_packet(Point pos, Point dir, float distance)
{
	pos.x += distance*dir.x;
	pos.y += distance*dir.y;
	pos.z += distance*dir.z;
	return pos;
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

typedef struct Interfacer
{
	UINT32 index;
	float distance;
} Interfacer;

Interfacer get_face_intercept_index(Tetra tetra, Point pos, Point dir)
{
	float dis[4];
	float cosdn[4];
	// = |dir|cos0
	cosdn[0] = project_on_unit(dir, tetra.face[0]);
	cosdn[1] = project_on_unit(dir, tetra.face[1]);
	cosdn[2] = project_on_unit(dir, tetra.face[2]);
	cosdn[3] = project_on_unit(dir, tetra.face[3]);

	// (p*n) - d(plane constant) = |p-p0|cos0
	dis[0] = distance_to_plane(pos, tetra.face[0]);
	dis[1] = distance_to_plane(pos, tetra.face[1]);
	dis[2] = distance_to_plane(pos, tetra.face[2]);
	dis[3] = distance_to_plane(pos, tetra.face[3]);

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

	Interfacer ifacer;
	ifacer.index = minIndex;
	ifacer.distance = dis[minIndex];
	return ifacer;
}

Spinner spin(Packet pkt, Material mat, float rand)
{
	// Spin
	float cost, sint, cosp, sinp;
	Point last_dir, last_a, last_b;
	//rand = FP_TWO * rand_MWC_co(&rnd_x, &rnd_a) - FP_ONE;
	rand = 2.0F * rand - 1.0F;
	cost = mat.HGCoeff1 - native_divide(mat.HGCoeff2, (1-mat.g*rand)*(1-mat.g*rand));
	sint = sqrt(1.0F - cost * cost);

	// spin psi 0-2pi.
	sinp = sincos(2.0F * PI * rand, &cosp);

	float stcp = sint * cosp;
	float stsp = sint * sinp;
	float ctcp = cost * cosp;
	float ctsp = cost * sinp;

	Spinner spinner;
	last_dir = pkt.dir;
	last_a = pkt.a;
	last_b = pkt.b;
	//pkt.d = cost*last_d - stcp*last_a + stsp*last_b;
	spinner.dir = add_points(sub_points(scalar_multiply(last_dir, cost), scalar_multiply(last_a, stcp)), scalar_multiply(last_b,stsp));
	//pkt.a = sint*last_d + ctcp*last_a - ctsp*last_b;
	spinner.a = sub_points(add_points(scalar_multiply(last_dir,sint), scalar_multiply(last_a, ctcp)), scalar_multiply(last_b,ctsp));
	//pkt.b = sinp*last_a + cosp*last_b;
	spinner.b = add_points(scalar_multiply(last_a,sinp), scalar_multiply(last_b,cosp));

	return spinner;
}

// TODO: Optimize, make unique struct for each stage to carry only necessary info forward
typedef struct Q_Payload
{
	bool bubble;
	bool jump;
	UINT32 ifaceID;
	float step;
	Packet pkt;
	Tetra tetra;
	Material mat;

	// TODO: this is wasteful pass only indices if possible
	Material mat2;
} Q_Payload;

typedef struct Q
{
	Q_Payload payload;
	bool full;
	UINT32 size;
} Q;

void q_push(Q* q, Q_Payload payload)
{
	// TODO: Check not full
	q->payload = payload;
	q->full = true;
	q->size++;
}

Q_Payload q_pop(Q* q)
{
	q->size--;
	q->full = false;
	return q->payload;
}

__kernel void entry(__global const UINT32* restrict photon_pool, __global const Source* restrict source, __global const Material* restrict materials, __global const Tetra* restrict mesh, __global const UINT32* restrict num_tetras, __global const UINT32* restrict num_mats, __global UINT32* restrict absorption, __global UINT32* restrict transmittance) {

	// PIPELINE
	//
	// Launch
	// |
	// |<-------------------^
	// v			|
	// Step			|
	// |			|
	// |<-----------^	|
	// |	  Look	|	|
	// v	  Ref	|	|
	// Hop -> Mat--->	|
	// |			|
	// v			|
	// Drop			|
	// Spin ---------------->

	Q launch2step_q;
	Q step2hop_q;
	Q hop2mat_q;
	Q mat2ref_q;
	Q ref2look_q;
	Q look2hop_q;
	Q hop2drop_q;
	Q drop2spin_q;
	Q spin2step_q;

	Q_Payload launcher;
	Q_Payload stepper;
	Q_Payload hopper;
	Q_Payload matter;
	Q_Payload reflactor;
	Q_Payload looker;
	Q_Payload dropper;
	Q_Payload spinner;

	RNG launcher_rng;
	RNG stepper_rng;
	RNG dropper_rng;
	RNG spinner_rng;
	RNG reflactor_rng;

	// RNG initialization
	// TODO: properly initialize seeds from host based on time
	launcher_rng.seed = 79;
	stepper_rng.seed = 80;
	dropper_rng.seed = 81;
	spinner_rng.seed = 82;
	reflactor_rng.seed = 83;

	// Launcher initialization
	launcher.tetra = mesh[source->tid];
	launcher.mat = materials[launcher.tetra.matID];

	launcher.pkt.weight = 1.0F;
	launcher.pkt.current_tid = source->tid;

	launcher.pkt.pos.x = source->pos.x;
	launcher.pkt.pos.y = source->pos.y;
	launcher.pkt.pos.z = source->pos.z;

	UINT32 launcher_photon_pool = *photon_pool;


	// Don't want to copy the queue's, pass pointer
	// TODO: Due to pointer usage, will this end up in memory (bad)? or just somewhere local (good)? Experiment and evaluate performance/area.
	// TODO: conditional '?' more efficient by making mux obvious? Experiment or lookup RTL or google

	bool done = false;
	while(!done) {

	// Launcher Pipelining
	bool launcher_stall = !launcher.bubble && launch2step_q.full;
	bool launcher_push = !launcher.bubble && !launch2step_q.full;
	if(launcher_push) q_push(&launch2step_q, launcher);

	// Stepper Pipelining
	bool stepper_stall = !stepper.bubble && step2hop_q.full;
	bool stepper_push = !stepper.bubble && !step2hop_q.full;
	if(stepper_push) q_push(&step2hop_q, stepper);
	if(!stepper_stall) {
		stepper.bubble = (spin2step_q.size == 0) && (launch2step_q.size == 0);
		if(!stepper.bubble) {
			stepper = (spin2step_q.size > 0) ? q_pop(&spin2step_q) : q_pop(&launch2step_q);
		}
	}

	// Hopper Pipelining
	bool hopper_stall = !hopper.bubble && ((hopper.hop && hop2mat_q.full) || (!hopper.hop && hop2drop_q.full));
	bool hopper_push = !hopper.bubble && ((hopper.hop && !hop2mat_q.full) || (!hopper.hop && !hop2drop_q.full));
	Q* hopper2q = hopper.hop ? &hop2mat_q : &hop2drop_q;
	if(hopper_push) q_push(hopper2q, hopper);
	if(!hopper_stall) {
		hopper.bubble = (step2hop_q.size == 0) && (look2hop_q.size == 0);
		if(!hopper.bubble) {
			// Prioritize
			hopper = (look2hop_q.size > 0) ? q_pop(&look2hop_q) : q_pop(&step2hop_q);
		}
	}

	// Dropper Pipelining
	bool dropper_stall = !dropper.bubble && drop2spin_q.full;
	bool dropper_push = !dropper.bubble && !drop2spin_q.full;
	if(dropper_push) q_push(&drop2spin_q, dropper);
	if(!dropper_stall) {
		dropper.bubble = (hop2drop_q.size == 0);
		if(!dropper.bubble) {
			dropper = q_pop(&hop2drop_q);
		}
	}

	// Spinner Pipelining
	bool spinner_stall = !spinner.bubble && spin2step_q.full;
	bool spinner_push = !spinner.bubble && !spin2step_q.full;
	if(spinner_push) q_push(&spin2step_q, spinner);
	if(!spinner_stall) {
		spinner.bubble = (drop2spin_q.size == 0);
		if(!spinner.bubble) {
			spinner = q_pop(&drop2spin_q);
		}
	}

	// Matter Pipelining
	bool matter_stall = !matter.bubble && mat2ref_q.full;
	bool matter_push = !spinner.bubble && !mat2ref_q.full;
	if(matter_push) q_push(&mat2ref_q, matter);
	if(!matter_stall) {
		matter.bubble = (hop2mat_q.size == 0);
		if(!matter.bubble) {
			matter = q_pop(&hop2mat_q);
		}
	}

	// Reflactor Pipelining
	bool reflactor_stall = !reflactor.bubble && ref2look_q.full;
	bool reflactor_push = !reflactor.bubble && !ref2look_q.full;
	if(reflactor_push) q_push(&ref2look_q, reflactor);
	if(!reflactor_stall) {
		reflactor.bubble = (mat2ref_q.size == 0);
		if(!reflactor.bubble) {
			reflactor = q_pop(&mat2ref_q);
		}
	}

	// Looker Pipelining
	bool looker_stall = !looker.bubble && look2hop_q.full;
	bool looker_push = !looker.bubble && !look2hop_q.full;
	if(looker_push) q_push(&look2hop_q, looker);
	if(!looker_stall) {
		looker.bubble = (ref2look_q.size == 0);
		if(!looker.bubble) {
			looker = q_pop(&ref2look_q);
		}
	}

	// Launcher work
	if(!launcher_stall) {
		launcher.bubble = *photon_pool > 0;
		launcher_photon_pool--;

		launcher_rng = unit_rand(launcher_rng.seed);
		UINT32 rand1 = launcher_rng.rand;
		launcher_rng = unit_rand(launcher_rng.seed);
		UINT32 rand2 = launcher_rng.rand;
		
		// Direction	
		float sint, sinp, cost, cosp;
		float theta = PI*rand1;
		float phi = 2.0F*PI*rand2;
		sint = sincos(theta, &cost);
		sinp = sincos(phi, &cosp);
	
		launcher.pkt.dir.x = sinp*cost;
		launcher.pkt.dir.y = sinp*sint;
		launcher.pkt.dir.z = cosp;
		Spinner spin_helper = update_spin_helpers(launcher.pkt.dir);
		launcher.pkt.a = spin_helper.a;
		launcher.pkt.b = spin_helper.b;
	}
	
	// Stepper work
	if(!stepper_stall) {
		stepper_rng = unit_rand(stepper_rng.seed);
		stepper.step = -log(stepper_rng.rand);
	}
	
	// Hopper work
	if(!hopper_stall) {
		// Get the minimum distances from the pkt to each face and the index of which will hit first
		Interfacer ifacer  = get_face_intercept_index(hopper.tetra, hopper.pkt.pos, hopper.pkt.dir);
		hopper.ifaceID = ifacer.index;
		float iface_distance = ifacer.distance;

		float travel_distance = hopper.step * hopper.mat.rmu_as;
		bool hop = travel_distance >= iface_distance;
		float move_distance = hop ? iface_distance : travel_distance;
		hopper.pkt.pos = move_packet(hopper.pkt.pos, hopper.pkt.dir, move_distance);

		// Pipeline simplification
		// Always subtract from step. This is correct if it hops.
		// If it doesn't hop a new step will simply be drawn
		hopper.step -= iface_distance * hopper.mat.mu_as;
	}
	
	// Dropper work
	if(!dropper_stall) {
		float dw = dropper.pkt.weight * dropper.mat.absfrac;

		// This is a global memory access
		ATOMIC_ADD((absorption[dropper.pkt.current_tid]), (UINT32)(dw*WEIGHT_SCALE));

		dropper.pkt.weight -= dw;

		// Cons. of Energy guaranteed? No, just monte carlo approximated.
		// This pkt survives the roulette.
		dropper_rng = unit_rand(dropper_rng.seed);
		bool roulette = (dropper.pkt.weight < WEIGHT && dropper.pkt.weight != 0.0F);
		bool survived = roulette && (dropper_rng.rand < CHANCE);
		dropper.bubble = (dropper.pkt.weight == 0.0F) || (roulette && !survived);
		if(roulette && survived) {
			// TODO: Precompute 1/chance and #define it
			dropper.pkt.weight *= (1.0F / CHANCE);
		}	
	}
	
	// Spinner work
	if(!spinner_stall) {
		spinner_rng = unit_rand(spinner_rng.seed);
		Spinner spinner_t = spin(spinner.pkt, spinner.mat, spinner_rng.rand);
		spinner.pkt.dir = spinner_t.dir;
		spinner.pkt.a = spinner_t.a;
		spinner.pkt.b = spinner_t.b;
	}
	
	// Matter work
	if(!matter_stall) {
		// TODO: Better to store adjMat ID's? Or to just look up the tetra and grab its mat ID?
		matter.mat2 = materials[matter.tetra.adjMats[matter.ifaceID]];
	}

	// Reflactor work
	if(!reflactor_stall) {
		// Refractive indices
		float ni = reflactor.mat.n;
		float nt = reflactor.mat2.n;

		float crit_cos = nt < ni ? 0.0F : GetCosCrit(ni,nt);
		float ca1 = -project_on_unit(reflactor.pkt.dir, reflactor.tetra.face[reflactor.ifaceID]); //ca1 is cos(theta incidence)

		// TODO: any special reason for doing in this order?
		// TODO: Waste of compute to do sqrt before check
		float sa1 = sqrt(1.0F-ca1*ca1); //sa1 is sin(theta incidence)
		if (ca1 > COSZERO) sa1 = 0.0F; 

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

		// Calculate reflected direction
		Point reflect_dir = reflect_packet(reflactor.pkt, ca1, reflactor.tetra.face[reflactor.ifaceID]);

		// Calculate refracted direction
		// refracted direction = n1/n2*original_direction - [n1/n2*cos(theta incidence) + sqrt(1-sin^2(theta transmitt))]*normal
		//pkt.dir = ni_nt*(pkt.dir) - fma(ni_nt,-ca1,ca2)*(float4)(normal[0],normal[1],normal[2],0);
		float4 temp = fma(ni_nt,-ca1,ca2)*(float4)(reflactor.tetra.face[reflactor.ifaceID][0],reflactor.tetra.face[reflactor.ifaceID][1],reflactor.tetra.face[reflactor.ifaceID][2],0);
		Point temp2;
		temp2.x = temp.x;
		temp2.y = temp.y;
		temp2.z = temp.z;
		Point refract_dir = sub_points(scalar_multiply(reflactor.pkt.dir, ni_nt), temp2);

		reflactor_rng = unit_rand(reflactor_rng.seed);
		bool reflect_refract = (ca1 <= crit_cos) || (rFresnel >= reflactor_rng.rand);

		if(ni != nt) {
			reflactor.pkt.dir = reflect_refract ? reflect_dir : refract_dir;

			// TODO: update_spin_helpers
		}

		reflactor.jump = (ni == nt || !reflect_refract);
		if(reflactor.jump) {
			reflactor.pkt.current_tid = reflactor.tetra.adjTetras[reflactor.ifaceID];
		}

		// Packet has left the model
		reflactor.bubble = (reflactor.jump && reflactor.tetra.adjTetras[reflactor.ifaceID] == 0);
		if(reflactor.bubble) {
			ATOMIC_ADD((transmittance[(reflactor.pkt.current_tid - 1) * 4 + reflactor.ifaceID]), (UINT32)(reflactor.pkt.weight * WEIGHT_SCALE));
		}
	}

	// Looker work
	if(!looker_stall) {
		if(looker.jump) {
			looker.tetra = mesh[looker.tetra.adjTetras[looker.ifaceID]];
		}
	}

	done = launcher.bubble && stepper.bubble && hopper.bubble && matter.bubble && reflactor.bubble && looker.bubble && dropper.bubble && spinner.bubble;
	}

	/*UINT32 seed_t = 80;
	UINT32* seed = &seed_t;
	Source launcher = *source;
	Packet launch_pkt = launch(launcher, seed);
	Packet step_pkt = draw_step(seed);
	Tetra tetra = mesh[pkt.current_tid];
	// Get the minimum distances from the pkt to each face and the index of which will hit first
	float dis[4];
	UINT32 minIndex = get_face_intercept_index(tetra, pkt, dis);
	TetraID nextTetraID = tetra.adjTetras[minIndex];

	float r = unit_rand(seed);
	Material mat = materials[tetra.matID];
	float travel_distance = pkt.step * mat.rmu_as;
	Spinner spinner = spin(pkt, mat, r, seed);
	
	bool hop = travel_distance >= dis[minIndex];
	pkt.step = hop ? pkt.step - dis[minIndex] * mat.mu_as : 0.0F;
	move_distance = hop ? dis[minIndex] : travel_distance;
	Tetra nextTetra = hop ? mesh[nextTetraID] : tetra;
	Material nextMat = hop ? materials[nextTetra.matID] : mat;

	pkt = move_packet(pkt, move_distance);



	float ni, nt; //refractive indices
	ni = mat.n;
	nt = nextMat.n;

	float crit_cos = nt < ni ? 0.0F : GetCosCrit(ni,nt);

	float ca1 = -project_on_unit(pkt.dir, tetra.face[minIndex]); //ca1 is cos(theta incidence)

	// TODO: any special reason for doing in this order?
	// TODO: Waste of compute to do sqrt before check
	float sa1 = sqrt(1.0F-ca1*ca1); //sa1 is sin(theta incidence)
	if (ca1 > COSZERO) sa1 = 0.0F; 
	

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
	// TODO: Should be 2.0F?
	float rFresnel = (Rs + Rp)/2;
	if (ca1 < COSNINETYDEG || sa2 == 1.0F) rFresnel = 1.0F;

	Point refract_dir;
	// refracted direction = n1/n2*original_direction - [n1/n2*cos(theta incidence) + sqrt(1-sin^2(theta transmitt))]*normal
	//pkt.dir = ni_nt*(pkt.dir) - fma(ni_nt,-ca1,ca2)*(float4)(normal[0],normal[1],normal[2],0);
	float4 temp = fma(ni_nt,-ca1,ca2)*(float4)(tetra.face[minIndex][0],tetra.face[minIndex][1],tetra.face[minIndex][2],0);
	Point temp2;
	temp2.x = temp.x;
	temp2.y = temp.y;
	temp2.z = temp.z;
	refract_dir = sub_points(scalar_multiply(pkt.dir, ni_nt), temp2);

	bool reflect_refract = (ca1 <= crit_cos) || (rFresnel >= r);
	Point hop_dir = (ni == nt) ? pkt.dir : (reflect_refract ? reflect_packet(pkt, ca1, tetra.face[minIndex]) : refract_dir);

	Point spin_dir = spinner.dir;
	Point spin_a = spinner.a;
	Point spin_b = spinner.b;

	Spinner spinner2 = update_spin_helpers(hop_dir);
	Point hop_a = spinner2.a;
	Point hop_b = spinner2.b;

	pkt.a = hop ? hop_a : stay_a;
	pkt.b = hop ? hop_b : stay_b;

	pkt.dir = hop ? hop_dir : spin_dir;

	if (((ni == nt) || (!reflect_refract)) && nextTetraID == 0) {
      		//store surface fluence TODO: find more efficient way of doing this!
		atomic_add(&(transmittance[(pkt.current_tid - 1) * 4 + minIndex]), (UINT32)(pkt.weight * WEIGHT_SCALE));
		// Kill the packet.
		pkt.weight = 0.0F;
	}

	float dw = pkt.weight*mat.absfrac;
	pkt.weight = hop ? pkt.weight : pkt.weight - dw;

	// Cons. of Energy guaranteed?
	// This pkt survives the roulette.
	pkt.weight = (pkt.weight < WEIGHT && pkt.weight != 0.0F && r < CHANCE) ? pkt.weight * (1.0F / CHANCE) : pkt.weight;

	if(!hop) {
		// This is a global memory access
		atomic_add(&(absorption[pkt.current_tid]), (UINT32)(dw*WEIGHT_SCALE));
	}

	pkt.current_tid = ((ni == nt) || (!reflect_refract)) ? nextTetraID : pkt.current_tid;*/
}

