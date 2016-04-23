
#ifndef MISC_H
#define MISC_H


#define MAX_LOCAL_MEM_SIZE 131072

#define CACHE_TAG_BITS 20
#define CACHE_TAG_MASK 0xFFFFF000

#define CACHE_INDEX_BITS 6
#define CACHE_NUM_BLOCKS 64
#define CACHE_INDEX_MASK 0x00000FC0

#define CACHE_OFFSET_BITS 6
#define CACHE_BLOCK_NUM_INTS 64
#define CACHE_OFFSET_MASK 0x0000003F

#define CACHE_BLOCK_MASK 0xFFFFFFC0

// DDR has wide load of 256 bytes
#define CACHE_BLOCK_SIZE 256
#define CACHE_TOTAL_TAG_SIZE 256 // sizeof(UINT32)*num_blocks
#define CACHE_TOTAL_DATA_SIZE 16384 // block_size*num_blocks

typedef unsigned long long UINT64;
typedef unsigned int UINT32;
typedef UINT32 TetraID;

typedef struct Tetra
{
	// Face i's plane is defined by equation ax + by + cz = d
	// The coefficients are found in face[i][0:3] respectively
	// Face i's normal vector face[i][0:2] has magnitude 1.
	// It always points into the tetrahedral.
	float face[4][4];

	// Adjacent tetrahedral ID's
	UINT32 adjTetras[4];

	// Adjacent material indices of refraction
	// Used to avoid two memory accesses if reflection occurs
	float adjN[4];

	// Material ID
	UINT32 matID;
} Tetra;

typedef struct Material
{
	float mu_a;	//absorption coefficient
	float mu_s;	//scatter coefficient
	float mu_as;	//attenuation coefficient: result of mu_a + mu_s
	float rmu_as;	//reciprocal of mu_as, store this to get rid of slow division arithmetic 
	float n;	//index of refraction
	float g;	//anisotropy constant
	float HGCoeff1;	// HGCoeff1 = (1+g^2)/(2g)
	float HGCoeff2;	// HGCoeff2 = (1-g^2)^2 /(2g). So cos(theta) = HGCoeff1 - HGCoeff2 / (1-g * rand(-1,1))
	float absfrac;	//absorb fraction = 1- albedo = 1 - mus / (mus+mua)
	int setMatched;
} Material;

#ifdef HOST
namespace FPGA
{
#endif

typedef struct Point
{
  float x;
  float y;
  float z;
} Point;

typedef struct Source
{
	Point pos;
	TetraID tid;
} Source;

#ifdef HOST
}
#endif

#endif // MISC_H

