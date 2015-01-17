#define NFLOATS 5
#define NINTS 5

#include <time.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <list>

#include <CL/cl.h>

#include "kernel.h"
#include "gpumcml.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////
//   Print Command Line Help - How to run program and pass in parameters 
//////////////////////////////////////////////////////////////////////////////
void usage(const char *prog_name)
{
  printf("\nUsage: %s [-A] [-S<seed>] [-G<num GPUs>] <input file>\n\n",
    prog_name);
  printf("  -A: ignore A detection\n");
  printf("  -S: seed for random number generation (MT only)\n");
  printf("  -G: set the number of GPUs this program uses\n");
  printf("\n");
  fflush(stdout);
}

//////////////////////////////////////////////////////////////////////////////
//   Parse command line arguments
//////////////////////////////////////////////////////////////////////////////
int interpret_arg(int argc, char* argv[], SimulationStruct *p_simulation)
{
  sscanf(argv[1],"%d",&(p_simulation->number_of_photons));
  strncpy(p_simulation->outp_filename, argv[2], STR_LEN);

  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//   Scale raw data and format data for file output 
//////////////////////////////////////////////////////////////////////////////
int Write_Simulation_Results(SimState* HostMem, SimulationStruct* sim, clock_t simulation_time)
{
  FILE* pFile_outp;
  pFile_outp = fopen (sim->outp_filename , "w");
  if (pFile_outp == NULL){perror ("Error opening output file");return 0;}
  fclose(pFile_outp);
  return 0;
}

void normalize(float &x, float &y, float &z)
{
  float magnitude = sqrt(x*x + y*y + z*z);
  x = x/magnitude;
  y = y/magnitude;
  z = z/magnitude;
}

void cross(float x1, float y1, float z1, float x2, float y2, float z2, float *p_xout, float *p_yout, float *p_zout)
{
  *p_xout = y1 * z2 - y2 * z1;
  *p_yout = z1 * x2 - z2 * x1;
  *p_zout = x1 * y2 - x2 * y1;
}

void Sort2Int(int *a, int *b)
{
  if (*a > *b)
  {
    int temp = *a;
    *a = *b;
    *b = temp;
  }
}

void Sort4Int(int *ids)
{
  //Get local mins and local maxes
  Sort2Int(ids, ids+1);
  Sort2Int(ids+2, ids+3);
  //Now min = ids[0] or ids[2]; max = ids[1] or ids[3];
  Sort2Int(ids, ids+2);
  Sort2Int(ids+1, ids+3);
  //Now min and max are properly located, only need to sort the middle two
  Sort2Int(ids+1, ids+2);
}

void PopulateFaceParameters(Tetra &tetra, Tetra &adjTetra, int tetraIndex, int adjTetraIndex, Point &p1, Point &p2, 
                            Point &p3, Point &p4)
{
  //1. Get vectors of any two sides of the face triangle
  float x1 = p1.x - p2.x;
  float y1 = p1.y - p2.y;
  float z1 = p1.z - p2.z;
  float x2 = p1.x - p3.x;
  float y2 = p1.y - p3.y;
  float z2 = p1.z - p3.z;
  //2. Take the cross product of those two sides to get the normal vector
  float nx,ny,nz;
  cross(x1,y1,z1,x2,y2,z2,&nx,&ny,&nz);
  //3. Normalize the normal vector to unit length
  normalize(nx,ny,nz);
  //4. Get the face constant
  float faceConstant = nx*p1.x + ny*p1.y + nz*p1.z;
  //5. Put nx,ny,nz and faceConstant to the correct tetrahedron by making sure (nx,ny,nz) always points into the tetra
  if(nx*p4.x+ny*p4.y+nz*p4.z > faceConstant)	//pointing into tetra
  {
    tetra.face[tetraIndex][0] = nx;
    tetra.face[tetraIndex][1] = ny;
    tetra.face[tetraIndex][2] = nz;
    tetra.face[tetraIndex][3] = faceConstant;
    adjTetra.face[adjTetraIndex][0] = -nx;
    adjTetra.face[adjTetraIndex][1] = -ny;
    adjTetra.face[adjTetraIndex][2] = -nz;
    adjTetra.face[adjTetraIndex][3] = -faceConstant;
  }
  else //pointing outside tetra, so into adjTetra
  {
    tetra.face[tetraIndex][0] = -nx;
    tetra.face[tetraIndex][1] = -ny;
    tetra.face[tetraIndex][2] = -nz;
    tetra.face[tetraIndex][3] = -faceConstant;
    adjTetra.face[adjTetraIndex][0] = nx;
    adjTetra.face[adjTetraIndex][1] = ny;
    adjTetra.face[adjTetraIndex][2] = nz;
    adjTetra.face[adjTetraIndex][3] = faceConstant;
  }
}

void HandleFace(int tetraID, Tetra *tetra_mesh, list<TwoPointIDsToTetraID> *faceToTetraMap, Point *points,
                          int lowID, int midID, int highID, int fourthID, int faceIndex)
{
  list<TwoPointIDsToTetraID> &nodeList = faceToTetraMap[lowID];
  TwoPointIDsToTetraID *node = NULL;
  int i;
  //Check if the nodeList already contains the adjacent tetra waiting to be paired.
  for(list<TwoPointIDsToTetraID>::iterator it = nodeList.begin();it != nodeList.end();it++)
  {
    if(it->lowerPointID == midID && it->higherPointID == highID)
    {
      node = &(*it);
      break;
    }
  }
  //TODO: may need mechanism to make sure no face is shared by more than two tetras.
  if(node == NULL)	//Did not find the node
  {
    //add the node to nodeList, indicating this Tetra is waiting for its adjacent one to be found later (if it exists).
    node = (TwoPointIDsToTetraID *)malloc(sizeof(TwoPointIDsToTetraID));
    node->lowerPointID = midID;
    node->higherPointID = highID;
    node->TetraID = tetraID;
    node->otherPointID = fourthID;
    nodeList.push_back(*node);
  }
  else	//the node is found
  {
    //set both tetras' adjacent tetra id as each other's
    Tetra &tetra = tetra_mesh[tetraID];
    tetra.adjTetras[faceIndex] = node->TetraID;
    Tetra &adjTetra = tetra_mesh[node->TetraID];
    //Look for the available place in face array
    for(i = 0; i < 4; i++)
    {
      if(adjTetra.face[i][0]==0 && adjTetra.face[i][1]==0 && adjTetra.face[i][2]==0 && adjTetra.adjTetras[i]==0) break;
    }
    //TODO: error checking for i==4?
	adjTetra.adjTetras[i] = tetraID;
	//Populate face parameters (normal vector and face constant)
	PopulateFaceParameters(tetra, adjTetra, faceIndex, i, points[lowID], points[midID], points[highID], points[fourthID]);
	//Delete the node
	nodeList.remove(*node);
  }
}


void PopulateTetraFromMeshFile(const char* filename, Tetra **p_tetra_mesh, int *p_Np, int *p_Nt)
{
  int i;
  int pointIDs[4];	//store the ids of 4 points of each tetrahedron
  FILE *pFile = fopen(filename , "r");
  if (pFile == NULL) 
  {
    printf ("Error Opening Mesh (.mesh) file.\n");
    fclose (pFile);
    return;
  }
  
  char line[64];
  fgets(line, 64, pFile);
  sscanf(line, "%d", p_Np);
  //make size+1 because the point ID refered by each tetrahedron in the input file starts from 1 not 0
  Point *points = (Point*)malloc(sizeof(Point)*(*p_Np+1));
  fgets(line, 64, pFile);
  sscanf(line, "%d", p_Nt);
  for(i=1; i<*p_Np+1; i++)
  {
    fgets(line, 64, pFile);
    sscanf(line, "%f %f %f", &(points[i].x), &(points[i].y), &(points[i].z));
  }
  
  *p_tetra_mesh = (Tetra*)malloc(sizeof(Tetra)*(*p_Nt+1));
  list<TwoPointIDsToTetraID> *faceToTetraMap = new list<TwoPointIDsToTetraID>[*p_Np+1];
  for(i=1; i<*p_Nt+1; i++)
  {
    fgets(line, 64, pFile);
    sscanf(line, "%d %d %d %d %d", &pointIDs[0], &pointIDs[1], &pointIDs[2], &pointIDs[3], &((*p_tetra_mesh)[i].matID));
    
    Sort4Int(pointIDs);
    //Find adjacent tetra IDs, and populate the face parameters into the tetra pair
    HandleFace(i, *p_tetra_mesh, faceToTetraMap, points, pointIDs[0], pointIDs[1], pointIDs[2], pointIDs[3], 0);
    HandleFace(i, *p_tetra_mesh, faceToTetraMap, points, pointIDs[0], pointIDs[1], pointIDs[3], pointIDs[2], 1);
    HandleFace(i, *p_tetra_mesh, faceToTetraMap, points, pointIDs[0], pointIDs[2], pointIDs[3], pointIDs[1], 2);
    HandleFace(i, *p_tetra_mesh, faceToTetraMap, points, pointIDs[1], pointIDs[2], pointIDs[3], pointIDs[0], 3);
  }
  fclose(pFile);
  //Now parameters of all the faces inside the mesh are populated, the rest are the mesh surfaces, which are the remaining nodes
  //in faceToTetraMap
  for(i=1; i<*p_Np+1; i++)
  {
    list<TwoPointIDsToTetraID> &nodeList = faceToTetraMap[i];
    while(!nodeList.empty())
    {
      TwoPointIDsToTetraID &firstNode = nodeList.front();

      Point p1 = points[i];
      Point p2 = points[firstNode.lowerPointID];
      Point p3 = points[firstNode.higherPointID];
      Point p4 = points[firstNode.otherPointID];
      Tetra &tetra = (*p_tetra_mesh)[firstNode.TetraID];
      int j;
      //Look for the available place in face array
      for(j=0;j<4;j++)
      {
        if(tetra.face[j][0]==0&&tetra.face[j][1]==0&&tetra.face[j][2]==0) break;
      }
      //TODO: error checking for j==4?
      //get normal vector
      float nx,ny,nz,faceConstant;
      cross(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z, p1.x-p3.x, p1.y-p3.y, p1.z-p3.z, &nx, &ny, &nz);
      normalize(nx,ny,nz);
      faceConstant = p1.x*nx+p1.y*ny+p1.z*nz;
      if (p4.x*nx+p4.y*ny+p4.z*nz>faceConstant)	//(nx,ny,nz) points into tetra
      {
        tetra.face[j][0] = nx;
        tetra.face[j][1] = ny;
        tetra.face[j][2] = nz;
        tetra.face[j][3] = faceConstant;
      }
      else
      {
        tetra.face[j][0] = -nx;
        tetra.face[j][1] = -ny;
        tetra.face[j][2] = -nz;
        tetra.face[j][3] = -faceConstant;
      }
      nodeList.remove(firstNode);
    }
  }
  free(points);
  delete[] faceToTetraMap;
  return;
}

void PopulateMaterialFromInput(const char* fileName, Material **p_material_spec, int *Nm)
{
  *Nm = 1;
  *p_material_spec = (Material *)malloc( sizeof(Material) * ((*Nm)+1) );
  Material &material = (*p_material_spec)[1];
  material.mu_as = (0.2+1);
  material.rmu_as = 1.0/(material.mu_as);
  material.n = 1.53;
  material.g = 0.9;
  material.HGCoeff1 = (1+material.g*material.g)/(2*material.g);
  material.HGCoeff2 = (1-material.g*material.g)*(1-material.g*material.g) / (2*material.g);
  material.absfrac = 1- 1000/(1000+0.2);
  
  Material &exterior = (*p_material_spec)[0];
  exterior.mu_as = (0.2+1000);
  exterior.rmu_as = 1.0/(exterior.rmu_as);
  exterior.n = 1;
  exterior.g = 0.9;
  exterior.HGCoeff1 = (1+exterior.g*material.g)/(2*exterior.g);
  exterior.HGCoeff2 = (1-exterior.g*material.g)*(1-exterior.g*exterior.g) / (2*exterior.g);
  exterior.absfrac = 1- 1000/(1000+0.2);
  return;
}

void ParseMaterial(const char* filename, Material** mat) 
{

  FILE *pFile;

  // Line 1 and 2 of the material (.opt) file
  // Line 1 -> Material Parameter Style -> Must be 1
  // Line 2 -> Number of materials present in the mesh
  char line1[100], line2[100], linemat[100];
  int paramstyle, numbermat, etype;
  float n_e;

  pFile = fopen(filename, "r");

  if (pFile == NULL) 
  {
    printf ("Error Opening Material (.opt) file.\n");
    fclose (pFile);
    return;
  }
  // Verify the material paramter style
  if (fgets(line1, 100, pFile) == NULL) {
    printf ("Error Reading Line 1 Of Material (.opt) File.\n");
    fclose (pFile);
    return;
  }
  sscanf (line1, "%d", &paramstyle);
  // Collect the number of materials present in the mesh
  if ( fgets (line1, 100, pFile) == NULL ) {
    printf ("Error Reading Line 2 Of Material (.opt) File.\n");
    fclose (pFile);
    return;
  }
  sscanf (line1, "%d", &numbermat);

  // printf("File opened successfully!\nMaterial parameter style: %d\nNumber of materials presnet:%d\n", paramstyle, numbermat);

  *mat = (Material*)malloc((numbermat+1)*sizeof(Material));

  /*
  if (materialproperties != NULL) {
  printf ("Structure created successfully\n");
  }*/

  for (int i = 1; i <= numbermat; i++) 
  {
    if ( fgets (linemat, 100, pFile) == NULL ) 
    {
      printf ("Error Reading Line %d Of Material (.opt) File.\n", i);
      fclose (pFile);
      return;
    }
    sscanf (linemat, "%f %f %f %f", &(mat[i]->mu_a), &(mat[i]->mu_s), &(mat[i]->g), &(mat[i]->n));
  }

  // Get E-Type
  if ( fgets (linemat, 100, pFile) == NULL ) 
  {
    printf ("Error Reading Line E-Type Of Material (.opt) File.\n");
    fclose (pFile);
    return;
  }

  sscanf (linemat, "%d", &etype);

  if (etype == 1) 
  {
    // Get n_e
    if ( fgets (linemat, 100, pFile) == NULL ) 
    {
      printf ("Error Reading Line n_e Of Material (.opt) File.\n");
      fclose (pFile);
      return;
    }
    sscanf (linemat, "%f", &n_e);
  }
  else if (etype == 2) 
  {
    // Set Matched
    mat[0]->setMatched = 1;

    // Get n_e
    if ( fgets (linemat, 100, pFile) == NULL ) 
    {
      printf ("Error Reading Line n_e Of Material (.opt) File.\n");
      fclose (pFile);
      return;
    }

    sscanf (linemat, "%f", &n_e);
  }
  else 
  {
    printf ("Invalid Environment Type\n");
    return;
  }

  mat[0]->mu_a = 0;
  mat[0]->mu_s = 0;
  mat[0]->g = 0;
  mat[0]->n = n_e;

  fclose (pFile);
}

//////////////////////////////////////////////////////////////////////////////
// Parse Source (.Source File)
//////////////////////////////////////////////////////////////////////////////

void ParseSource(char* filename, Source** sourcepoint) 
{
  FILE *sourcefile;

  // Line 1 and 2 of the material (.opt) file
  // Line 1 -> Material Parameter Style -> Must be 1
  // Line 2 -> Number of materials present in the mesh
  char line1[100], linesource[100];
  int numsources, sourcetype;

  sourcefile = fopen(filename, "r");

  if (sourcefile == NULL) 
  {
    printf ("Error Opening Source file.\n");
    fclose (sourcefile);
    return;
  }
  // Collect The Number Of Sources
  if ( fgets (line1, 100, sourcefile) == NULL ) 
  {
    printf ("Error Reading Line 1 Of Source File.\n");
    fclose (sourcefile);
    return;
  }
  
  sscanf (line1, "%d", &numsources);

  if (numsources == 0) 
  {
    printf("Error - Zero Sources Specified\n");
    return;
  }

  *sourcepoint = (Source*)malloc((numsources)*sizeof(Source));


  for (int i = 0; i < numsources; i++) 
  {
    if ( fgets (linesource, 100, sourcefile) == NULL ) 
    {
      printf ("Error Reading Line %d Of Source File.\n", i);
      fclose (sourcefile);
      return;
    }

    sscanf (linesource, "%d ", &(sourcepoint[i]->stype) );
    switch(sourcepoint[i]->stype)
    {
      // Point Position -> Cartesian Coordinates
      case 1:
        sscanf (linesource, "%d %f %f %f %d", &(sourcepoint[i]->stype), &(sourcepoint[i]->x), &(sourcepoint[i]->y), &(sourcepoint[i]->z), &(sourcepoint[i]->Np) );
        break;

      // IDt position
      case 2:
        sscanf (linesource, "%d %d %d", &(sourcepoint[i]->stype), &(sourcepoint[i]->IDt), &(sourcepoint[i]->Np) );
        break;

      // IDt position
      case 11:
        sscanf (linesource, "%d %d %f %f %f %f %f %f %d", &(sourcepoint[i]->stype), &(sourcepoint[i]->IDt), &(sourcepoint[i]->x), &(sourcepoint[i]->y), &(sourcepoint[i]->z), &(sourcepoint[i]->dx), &(sourcepoint[i]->dy), &(sourcepoint[i]->dz), &(sourcepoint[i]->Np) );
        break;

      default:
        printf ("Unrecognized Source Type.\n");
        break;
    }
  }
  // free(sourcepoint);

  fclose (sourcefile);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

