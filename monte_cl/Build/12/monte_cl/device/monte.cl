
#include <misc.h>

#if defined(CACHE_TEST) && !defined(USE_CACHE)
#define USE_CACHE
#endif

#if defined(CACHE_TEST) || defined(SELF_TEST) || defined(DEBUG) || defined(LAPS)
#define ITERATION
#endif

#if defined(SELF_TEST) || defined(DEBUG)
#define PKT_ID
#endif

#define PI 3.1415926F
#define TWO_PI 6.2831853F
#define CHANCE 0.1F
#define CHANCE_INV 10.0F // 1.0F / CHANCE 
#define WEIGHT 1E-4F
#define WEIGHT_SCALE 65536//12000000

// TODO: Implementation specific to GPU for taylor series expansion approximation?
// TODO: Can FPGA give exact values?
#define COSNINETYDEG 1.0E-6F
#define COSZERO (1.0F - 1.0E-6F)

#define RAND_MAX 2147483647.0F

// Atomic is not required with only one pipeline
//#define ATOMIC_REQ
#ifdef ATOMIC_REQ
#define ATOMIC_ADD(a,b) atomic_add(&(a), (b))
#else
#define ATOMIC_ADD(a,b) (a += (b))
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

	Point a;
	Point b;

	float weight;
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
	Point p = pkt.dir;
	float temp = 2.0F*ca1;
	p.x += temp*normal[0];
	p.y += temp*normal[1];
	p.z += temp*normal[2];
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
	sinp = sincos(TWO_PI * rand, &cosp);

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
	Packet pkt;
	Tetra tetra;
	Material mat;
	UINT32 ifaceID;
#ifdef PKT_ID
	UINT32 payload_id;
#endif
	float step;
	bool bubble;
	bool hop;
	bool jump;
} Q_Payload;

typedef struct Q
{
	Q_Payload payload;
	UINT32 size;
	bool full;
} Q;

void q_push(Q* restrict q, Q_Payload payload)
{
	// TODO: Check not full
	q->payload = payload;
	q->full = true;
	q->size++;
}

Q_Payload q_pop(Q* restrict q)
{
	Q_Payload front = q->payload;
	q->size--;
	q->full = false;
#ifdef PKT_ID
	q->payload.payload_id = 0;
#endif
	return front;
}

typedef struct CacheBlock
{
	UINT32 element[CACHE_BLOCK_NUM_INTS];
} CacheBlock;

typedef UINT32 CacheTag;

void cache_load(__local CacheBlock* restrict blocks, __local CacheTag* restrict tags, UINT32 id, __global UINT32* restrict memory)
{
	UINT32 aligned_id = id & CACHE_BLOCK_MASK;
	UINT32 tag = ((id & CACHE_TAG_MASK) >> (CACHE_OFFSET_BITS + CACHE_INDEX_BITS));
	UINT32 index = ((id & CACHE_INDEX_MASK) >> CACHE_OFFSET_BITS);
	tags[index] = tag;

	// Memory accesses will be coalesced
	#pragma unroll
	for(int i = 0;i < CACHE_BLOCK_NUM_INTS;i++) {
		blocks[index].element[i] = memory[(aligned_id + i)];
	}
}

void cache_writeback(__local CacheBlock* restrict blocks, __local CacheTag* restrict tags, UINT32 id, __global UINT32* restrict memory)
{
	printf("Writeback\n");
	UINT32 aligned_id = id & CACHE_BLOCK_MASK;
	UINT32 tag = ((id & CACHE_TAG_MASK) >> (CACHE_OFFSET_BITS + CACHE_INDEX_BITS));
	UINT32 index = ((id & CACHE_INDEX_MASK) >> CACHE_OFFSET_BITS);

	// Memory accesses will be coalesced
	#pragma unroll
	for(int i = 0;i < CACHE_BLOCK_NUM_INTS;i++) {
		memory[(aligned_id + i)] = blocks[index].element[i];
	}
}

void cache_flush(__local CacheBlock* restrict blocks, __local CacheTag* restrict tags, __global UINT32* restrict memory)
{
	#pragma unroll
	for(int i = 0; i < CACHE_NUM_BLOCKS; i++) {
		if(tags[i] != 0xFFFFFFFF) {
			cache_writeback(blocks, tags, (tags[i] << (CACHE_OFFSET_BITS + CACHE_INDEX_BITS)) + (i << CACHE_OFFSET_BITS), memory);
		}
	}
}

UINT32 read_cache(__local CacheBlock* restrict blocks, __local CacheTag* restrict tags, UINT32 id, __global UINT32* restrict memory)
{
	UINT32 tag = ((id & CACHE_TAG_MASK) >> (CACHE_OFFSET_BITS + CACHE_INDEX_BITS));
	UINT32 index = ((id & CACHE_INDEX_MASK) >> CACHE_OFFSET_BITS);
	UINT32 offset = (id & CACHE_OFFSET_MASK);

	CacheTag cache_tag = tags[index];
	if(cache_tag == 0xFFFFFFFF) {
		// Miss
		cache_load(blocks, tags, id, memory);
		
	}
	else if(tag != cache_tag) {
		// Miss and Writeback
		cache_writeback(blocks, tags, id, memory);
		cache_load(blocks, tags, id, memory);
	}

	// Block is now loaded, Hit
	return blocks[index].element[offset];
}

void write_cache(__local CacheBlock* restrict blocks, __local CacheTag* restrict tags, UINT32 id, __global UINT32* restrict memory, UINT32 value)
{
	UINT32 tag = ((id & CACHE_TAG_MASK) >> (CACHE_OFFSET_BITS + CACHE_INDEX_BITS));
	UINT32 index = ((id & CACHE_INDEX_MASK) >> CACHE_OFFSET_BITS);
	UINT32 offset = (id & CACHE_OFFSET_MASK);

	CacheTag cache_tag = tags[index];
	if(cache_tag == 0xFFFFFFFF) {
		// Miss
		cache_load(blocks, tags, id, memory);
		
	}
	else if(tag != cache_tag) {
		// Miss and Writeback
		cache_writeback(blocks, tags, id, memory);
		cache_load(blocks, tags, id, memory);
	}

	// Block is now loaded, Hit
	blocks[index].element[offset] = value;

#ifdef CACHE_TEST
	// Writethrough
	memory[id] = value;
#endif
}

__attribute__((num_simd_work_items(1)))
__attribute__((num_compute_units(1)))
__attribute__((reqd_work_group_size(1,1,1)))
__kernel void entry (
	const UINT32 photon_pool,
	__global const Source* restrict source,
	__constant Material* restrict materials,
	__global const Tetra* restrict mesh,
//	const UINT32 num_tetras,
//	const UINT32 num_mats,
	__global UINT32* restrict absorption,
//	__global UINT32* restrict transmittance,
	const UINT32 launcher_rng_init,
	const UINT32 stepper_rng_init,
	const UINT32 dropper_rng_init,
	const UINT32 spinner_rng_init,
	const UINT32 reflactor_rng_init
#ifdef USE_CACHE
	,
	__attribute__((local_mem_size(CACHE_TOTAL_TAG_SIZE))) __local CacheTag* restrict cache_tags,
	__attribute__((local_mem_size(CACHE_TOTAL_DATA_SIZE))) __local CacheBlock* restrict cache_blocks
	#ifdef CACHE_TEST
	,
	__global UINT32* restrict absorption_cacheless
	#endif
#endif
#ifdef SELF_TEST
	,
	__attribute__((local_mem_size(MAX_LOCAL_MEM_SIZE))) __local UINT32* restrict deadpool
#endif
)
{
#ifdef SELF_TEST
	printf("Args: init_tid %u, pool_size %u\n", source->tid, photon_pool);
#endif

	// PIPELINE
	//
	// Launch
	// |
	// |<-------------------^
	// v			|
	// Step			|
	// |			|
	// |<-----------^	|
	// |	  Mat	|	|
	// v	  Look	|	|
	// Hop -> Ref--->	|
	// |			|
	// v			|
	// Drop			|
	// Spin ---------------->

	Q_Payload launcher;
	Q_Payload stepper;
	Q_Payload hopper, hopper2;
	Q_Payload matter;
	Q_Payload reflactor;
	Q_Payload looker;
	Q_Payload dropper;
	Q_Payload spinner;

	launcher.bubble = true;
	stepper.bubble = true;
	hopper.bubble = true;
	hopper2.bubble = true;
	matter.bubble = true;
	reflactor.bubble = true;
	looker.bubble = true;
	dropper.bubble = true;
	spinner.bubble = true;

	RNG launcher_rng;
	RNG stepper_rng;
	RNG dropper_rng;
	RNG spinner_rng;
	RNG reflactor_rng;

	// RNG initialization
	launcher_rng.seed = launcher_rng_init;
	stepper_rng.seed = stepper_rng_init;
	dropper_rng.seed = dropper_rng_init;
	spinner_rng.seed = spinner_rng_init;
	reflactor_rng.seed = reflactor_rng_init;

	// Launcher initialization
	Source l_source = *source;
	launcher.tetra = mesh[l_source.tid];
	launcher.mat = materials[launcher.tetra.matID];

	launcher.pkt.weight = 1.0F;
	launcher.pkt.current_tid = l_source.tid;
	launcher.pkt.pos = l_source.pos;

#ifdef USE_CACHE
	// Initialize cache blocks to invalid
	for(int i = 0;i < CACHE_NUM_BLOCKS;i++) {
		cache_tags[i] = 0xFFFFFFFF;
	}
#endif

#ifdef SELF_TEST
	float my_absorbed = 0.0;
	float my_transmitted = 0.0;
	for(UINT32 i = photon_pool/32; i > 0; i--) {
		deadpool[i] = 0;
	}
#endif

#ifdef ITERATION
	UINT32 iteration = 1;
#endif

	UINT32 launcher_photon_pool = photon_pool;
	bool done = false;
	bool pool_empty = false;
	while(!done) {
#ifdef DEBUG
		printf("-------------Iteration %d-------------\n", iteration);
#endif
		// Pipelining

		// Stall Note:
		// Stall only when work has been done on a stage
		// and the stage was not able to push or kill its packet
		// In this way work will occur on bubble stages to simplify logic

		// Flow direction based on previous state
		Q_Payload temp = hopper;
		Q_Payload temp2 = hopper2;

		bool hopper_stall = (!temp.bubble && temp.hop) && (!temp2.bubble && temp2.hop);
		bool hopper2_stall = (!temp.bubble && !temp.hop) && (!temp2.bubble && !temp2.hop);

		if(!hopper_stall) {
			hopper = stepper;
			stepper.bubble = true;
		}
		
		if(!hopper2_stall) {
			hopper2 = matter;
			matter.bubble = true;
		}

		// Drop path pipelining
		bool stepper_stall = !stepper.bubble;
		if(!stepper_stall) {
			stepper = spinner.bubble ? launcher : spinner;
			if(spinner.bubble) launcher.bubble = true;
			spinner.bubble = true; // Simplify logic, always true
		}
		bool spinner_stall = !spinner.bubble;
		if(!spinner_stall) {
			spinner = dropper;
			dropper.bubble = true;
		}
		bool dropper_stall = !dropper.bubble;
		if(!dropper_stall) dropper = (!temp.bubble && !temp.hop) ? temp : (!temp2.bubble && !temp2.hop) ? temp2 : dropper;
		bool launcher_stall = !launcher.bubble;

		// Iface path pipelining
		bool matter_stall = !matter.bubble;
		if(!matter_stall) {
			matter = looker;
			looker.bubble = true;
		}
		bool looker_stall = !looker.bubble;
		if(!looker_stall) {
			looker = reflactor;
			reflactor.bubble = true;
		}
		bool reflactor_stall = !reflactor.bubble;
		if(!reflactor_stall) reflactor = (!temp2.bubble && temp2.hop) ? temp2 : (!temp.bubble && temp.hop) ? temp : reflactor;

#ifdef PKT_ID
		if(launcher.bubble) launcher.payload_id = 0;
		if(stepper.bubble) stepper.payload_id = 0;
		if(hopper.bubble) hopper.payload_id = 0;
		if(hopper2.bubble) hopper2.payload_id = 0;
		if(dropper.bubble) dropper.payload_id = 0;
		if(spinner.bubble) spinner.payload_id = 0;
		if(reflactor.bubble) reflactor.payload_id = 0;
		if(looker.bubble) looker.payload_id = 0;
		if(matter.bubble) matter.payload_id = 0;
#endif

		// Work

		// Launcher work
		if(!launcher_stall) {
			// Pool empty avoids if statement around pool--
#ifdef PKT_ID
			launcher.payload_id = pool_empty ? 0 : launcher_photon_pool;
#endif
			pool_empty = launcher_photon_pool == 0 || pool_empty;
			launcher.bubble = pool_empty;

			// Simplify logic, do work anyway even if pool is empty
			launcher_photon_pool--;
			launcher_rng = unit_rand(launcher_rng.seed);
			float rand1 = launcher_rng.rand;
			launcher_rng = unit_rand(launcher_rng.seed);
			float rand2 = launcher_rng.rand;
			
			// Direction	
			float sint, sinp, cost, cosp;
			float theta = PI*rand1;
			float phi = TWO_PI*rand2;
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
		if(!stepper.bubble && !stepper_stall) {
			stepper_rng = unit_rand(stepper_rng.seed);
			stepper.step = -log(stepper_rng.rand);
		} 

		// Hopper work
		if(!hopper.bubble && !hopper_stall) {
			// Get the minimum distances from the pkt to each face and the index of which will hit first
			Interfacer ifacer  = get_face_intercept_index(hopper.tetra, hopper.pkt.pos, hopper.pkt.dir);
			hopper.ifaceID = ifacer.index;
			float iface_distance = ifacer.distance;

			float travel_distance = hopper.step * hopper.mat.rmu_as;
			hopper.hop = travel_distance >= iface_distance;
			float move_distance = hopper.hop ? iface_distance : travel_distance;
			hopper.pkt.pos = move_packet(hopper.pkt.pos, hopper.pkt.dir, move_distance);

			// Pipeline simplification
			// Always subtract from step. This is correct if it hops.
			// If it doesn't hop a new step will simply be drawn
			hopper.step -= iface_distance * hopper.mat.mu_as;
		}

		// Hopper2 work (should be identical to hopper1)
		if(!hopper2.bubble && !hopper2_stall) {
			// Get the minimum distances from the pkt to each face and the index of which will hit first
			Interfacer ifacer  = get_face_intercept_index(hopper2.tetra, hopper2.pkt.pos, hopper2.pkt.dir);
			hopper2.ifaceID = ifacer.index;
			float iface_distance = ifacer.distance;

			float travel_distance = hopper2.step * hopper2.mat.rmu_as;
			hopper2.hop = travel_distance >= iface_distance;
			float move_distance = hopper2.hop ? iface_distance : travel_distance;
			hopper2.pkt.pos = move_packet(hopper2.pkt.pos, hopper2.pkt.dir, move_distance);

			// Pipeline simplification
			// Always subtract from step. This is correct if it hops.
			// If it doesn't hop a new step will simply be drawn
			hopper2.step -= iface_distance * hopper2.mat.mu_as;
		}

		// Dropper work
		if(!dropper.bubble && !dropper_stall) {
			float dw = dropper.pkt.weight * dropper.mat.absfrac;
	
			// This is a global memory access
#ifdef USE_CACHE
			printf("Iteration %u\n", iteration);
			UINT32 absorb = read_cache(cache_blocks, cache_tags, dropper.pkt.current_tid, absorption);
			absorb += (UINT32)(dw*WEIGHT_SCALE);
			write_cache(cache_blocks, cache_tags, dropper.pkt.current_tid, absorption, absorb);

	#ifdef CACHE_TEST
			ATOMIC_ADD((absorption_cacheless[dropper.pkt.current_tid]), (UINT32)(dw*WEIGHT_SCALE));
			if(read_cache(cache_blocks, cache_tags, dropper.pkt.current_tid, absorption) != absorption_cacheless[dropper.pkt.current_tid]) {
				printf("Cache Divergence! Iteration: %u Cached: %u : Actual: %u\n", iteration, read_cache(cache_blocks, cache_tags, dropper.pkt.current_tid, absorption), absorption_cacheless[dropper.pkt.current_tid]);
			}
	#endif
#else
			ATOMIC_ADD((absorption[dropper.pkt.current_tid]), (UINT32)(dw*WEIGHT_SCALE));
#endif

			dropper.pkt.weight -= dw;
	
			// Cons. of Energy guaranteed? No, just monte carlo approximated.
			// This pkt survives the roulette.
			dropper_rng = unit_rand(dropper_rng.seed);
			bool roulette = (dropper.pkt.weight < WEIGHT && dropper.pkt.weight != 0.0F);
			bool survived = roulette && (dropper_rng.rand < CHANCE);
			dropper.bubble = (dropper.pkt.weight == 0.0F) || (roulette && !survived);
			if(roulette && survived) {
				// TODO: Precompute 1/chance and #define it
				dropper.pkt.weight *= (CHANCE_INV);
			}

#ifdef SELF_TEST
			my_absorbed += dw;
			if(dropper.bubble) {
				deadpool[dropper.payload_id/32] |= (1 << (dropper.payload_id % 32));
			}	
#endif	
		}

		// Spinner work
		if(!spinner.bubble && !spinner_stall) {
			spinner_rng = unit_rand(spinner_rng.seed);
			Spinner spinner_t = spin(spinner.pkt, spinner.mat, spinner_rng.rand);
			spinner.pkt.dir = spinner_t.dir;
			spinner.pkt.a = spinner_t.a;
			spinner.pkt.b = spinner_t.b;
		}

		// Reflactor work
		if(!reflactor.bubble && !reflactor_stall) {
			// Refractive indices
			float ni = reflactor.mat.n;
			float nt = reflactor.tetra.adjN[reflactor.ifaceID];
	
			float crit_cos = nt < ni ? 0.0F : GetCosCrit(ni,nt);
			float ca1 = -project_on_unit(reflactor.pkt.dir, reflactor.tetra.face[reflactor.ifaceID]); //ca1 is cos(theta incidence)

			float sa1 = sqrt(1.0F-ca1*ca1); //sa1 is sin(theta incidence)
			if (ca1 > COSZERO) sa1 = 0.0F; 

			// Calculate Fresnel
			// Rs = [(n1 cos(theta_i) - n2 cos(theta_t))/((n1 cos(theta_i) + n2 cos(theta_t)))] ^ 2
			// Rp = [(n1 cos(theta_t) - n2 cos(theta_i))/((n1 cos(theta_t) + n2 cos(theta_i)))] ^ 2

			float ni_nt = native_divide(ni, nt);
			float sa2 = min(ni_nt * sa1, 1.0F); //sa2 is sin(theta transmit)
			float ca2 = sqrt(1.0F - sa2*sa2); //ca2 is cos(theta transmit)
			float Rs = native_divide(ni*ca1 - nt*ca2, ni*ca1 + nt*ca2);
			Rs *= Rs;
			float Rp = native_divide(ni*ca2 - nt*ca1, ni*ca2 + nt*ca1);
			Rp *= Rp;
				
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

				Spinner spin_helper = update_spin_helpers(reflactor.pkt.dir);
				reflactor.pkt.a = spin_helper.a;
				reflactor.pkt.b = spin_helper.b;
			}

			reflactor.jump = (ni == nt || !reflect_refract);
			if(reflactor.jump) {
				reflactor.pkt.current_tid = reflactor.tetra.adjTetras[reflactor.ifaceID];
			}

			// Packet has left the model
			reflactor.bubble = (reflactor.pkt.current_tid == 0);
			if(reflactor.bubble) {
				// Broken line at run time
				//ATOMIC_ADD((transmittance[(reflactor.pkt.current_tid - 1) * 4 + reflactor.ifaceID]), (UINT32)(reflactor.pkt.weight * WEIGHT_SCALE));

#ifdef SELF_TEST
				my_transmitted += reflactor.pkt.weight;
				if(reflactor.bubble) {
					deadpool[reflactor.payload_id/32] |= (1 << (reflactor.payload_id % 32));
				}	
#endif
			}
		}

		// Looker work
		if(!looker.bubble && !looker_stall) {
			if(looker.jump) {
				looker.tetra = mesh[looker.tetra.adjTetras[looker.ifaceID]];
			}
		}

		// Matter work
		if(!matter.bubble && !matter_stall) {
			// TODO: store adj mat IDs in tetra so tetra and mat lookups
			// occur in parallel stages as opposed to serial?
			// This assumes the materials and tetras are stored in a way that
			// supports parallel memory lookups, perhaps they're in different blocks?
			// Experiment.
			matter.mat = materials[matter.tetra.matID];
		}

		done = launcher.bubble && stepper.bubble && hopper.bubble && hopper2.bubble && matter.bubble && reflactor.bubble && looker.bubble && dropper.bubble && spinner.bubble;

#ifdef DEBUG
		printf("%d-%d-%d %d-%d\n", launcher.payload_id, stepper.payload_id, hopper.payload_id, hopper2.payload_id, matter.payload_id);
		printf("  | |X| |\n");
		printf("  %d-%d %d-%d\n", spinner.payload_id, dropper.payload_id, reflactor.payload_id, looker.payload_id);
#endif

#ifdef LAPS
		done = (iteration == LAPS);
#endif

#ifdef ITERATION
		iteration++;
#endif
	}

#ifdef USE_CACHE
	cache_flush(cache_blocks, cache_tags, absorption);
#endif

#ifdef SELF_TEST
	printf("FPGA to Earth...Simulation Complete in %d iterations!\n", iteration - 1);
	printf("FPGA Absorbed: %f\n", my_absorbed);
	printf("FPGA Transmitted: %f\n", my_transmitted);
	printf("FPGA Total: %f\n", my_absorbed + my_transmitted);

	bool test_packets_dropped = false;
	for(UINT32 i = photon_pool; i > 0; i--) {
		if(((deadpool[i/32] >> (i % 32)) & 0x1) == 0) {
			printf("FPGA dropped packet %u %u\n", i, deadpool[i/32]);
			test_packets_dropped = true;
		}
	}

	if(test_packets_dropped) {
		printf("Test Failed: packets dropped\n");
	}
#endif

}

