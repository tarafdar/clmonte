
#ifndef MISC_H
#define MISC_H

#define MAX_LOCAL_MEM_SIZE 131072

typedef unsigned long long UINT64;
typedef unsigned int UINT32;
typedef UINT32 TetraID;

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

typedef struct Tetra
{
  // The face i's plane is defined by equation (face[i][0]) * x + (face[i][1]) * y + (face[i][2]) * z = face[i][3]
  // The face normal vectors (face[i][0],face[i][1],face[i][2]) always point into the tetrahedron. They are always unit vector.
  float face[4][4];
  
  UINT32 adjTetras[4];
  float adjN[4];
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
}
#endif

#endif // MISC_H

