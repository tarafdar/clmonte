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
#include <fstream>
#include <iostream>
#include <CL/cl.h>
#include "host/h/gpumcml.h"
#include "AOCLUtils/aocl_utils.h"

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
  strncpy(p_simulation->inp_filename, argv[2], STR_LEN);
  strncpy(p_simulation->outp_filename, argv[3], STR_LEN);

  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//   Verify conservation of energy
//////////////////////////////////////////////////////////////////////////////
int Conservation_Of_Energy(SimState* HostMem, SimulationStruct* sim, TriNode* TriNodeList)
{
  float output_energy = 0;
  float absorption = 0;
  float transmittance = 0;
  // First sum the total absorption energy
  int i;
  for (i = 1; i <= sim->nTetras; i++)
    {
      output_energy += (double)(HostMem->absorption[i])/(double)WEIGHT_SCALE;
      absorption += (double)(HostMem->absorption[i])/(double)WEIGHT_SCALE;
    }

  // Next sum up the transmitted energy
  int TriNodeIndex = 0;
  for(int TetraID = 1; TetraID <= sim->nTetras; TetraID++)
  {
    for(int FaceID = 0; FaceID < 4; FaceID++)
    {
      TriNodeIndex = (TetraID-1)*4 + FaceID;

      if(TriNodeList[TriNodeIndex].N0!=0)
      {
        output_energy += (double) HostMem->transmittance[TriNodeIndex]/(double)(WEIGHT_SCALE);
        transmittance += (double) HostMem->transmittance[TriNodeIndex]/(double)(WEIGHT_SCALE);
      }else
      {
        //do nothing, not a surface element or just no light coming out of this part of the surface
      }
    }
  }

  printf ("Total Output Energy:	%f\n", output_energy);
  printf ("transmittance ratio: %f\n", transmittance/output_energy);
  printf ("absorption ratio: %f\n", absorption/output_energy);
  return 0;
}


//////////////////////////////////////////////////////////////////////////////
//   Scale raw data and format data for file output 
//////////////////////////////////////////////////////////////////////////////
int Write_Simulation_Results(SimState* HostMem, SimulationStruct* sim, double simulation_time, TriNode* TriNodeList, TetraNode* TetraNodeList, Material* material_spec, Tetra* tetra_mesh, char* output_filename)
{
  ofstream fout;
  
  fout.open(output_filename, ofstream::out);
  if(!fout.good()){
    cerr << "Could not write to output_filename." << endl;
    cerr << "Please check the write permission in this folder" << endl;
    fout.open("timos_tmp_result.dat", ofstream::out);
    if(fout.good()){
      cerr << "Current result will be write to: /tmp/timos_tmp_result.dat" << endl;
    }else{
      return 1;
    }
  }

  // Write simulation time
  fout <<"Wall clock simulation time: " << simulation_time << " seconds\n";

  // Calculate and write surface fluence!
  
  int TriNodeIndex = 0;
  list<TriNode> surfaceList;
  for(int TetraID = 1; TetraID <= sim->nTetras; ++TetraID)
  {
    for(int FaceID = 0; FaceID < 4; FaceID++)
    {
      TriNodeIndex = (TetraID-1)*4 + FaceID;
      //if (HostMem->transmittance[TriNodeIndex] > 0 ){
      if(TriNodeList[TriNodeIndex].N0!=0)
      {
      
        TriNodeList[TriNodeIndex].fluence = (double) HostMem->transmittance[TriNodeIndex]/WEIGHT_SCALE;
      
        surfaceList.push_back(TriNodeList[TriNodeIndex]);
      } 
      else
      {
      //do nothing, not a surface element
      }
      
    }
  }
  
  surfaceList.sort();
  fout << "Surface Fluence" << endl;
  for(list<TriNode>::iterator it=surfaceList.begin(); it!=surfaceList.end(); ++it)
  {
    fout << it->N0 << " \t"
         << it->N1 << " \t"
         << it->N2 << " \t"
         << it->area << " \t"
         << it->fluence << endl;
  }
  
  // Calculate and write absorption!
  fout << "Internal Fluence" << endl;
  for(int TetraID = 1; TetraID <= sim->nTetras; ++TetraID)
  {
      float fluence = (float) HostMem->absorption[TetraID]/WEIGHT_SCALE;
      //fluence = (double) HostMem->absorption[TetraID];
      fout << TetraNodeList[TetraID].N0 << " \t"
           << TetraNodeList[TetraID].N1 << " \t"
           << TetraNodeList[TetraID].N2 << " \t"
           << TetraNodeList[TetraID].N3 << " \t"
           << TetraNodeList[TetraID].volume << " \t"
           << fluence << endl;
  }

  fout.close();
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

float ComputeAreaFrom3Points(Point p1, Point p2, Point p3)
{
  //use Heron's formula
  float side1 = sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z));
  float side2 = sqrt((p1.x-p3.x)*(p1.x-p3.x)+(p1.y-p3.y)*(p1.y-p3.y)+(p1.z-p3.z)*(p1.z-p3.z));
  float side3 = sqrt((p3.x-p2.x)*(p3.x-p2.x)+(p3.y-p2.y)*(p3.y-p2.y)+(p3.z-p2.z)*(p3.z-p2.z));
  float s = (side1+side2+side3)/2; //semiperimeter
  return sqrt(s*(s-side1)*(s-side2)*(s-side3));
}

float ComputeVolume(Point p1, Point p2, Point p3, Point p4)
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
  //5. Get the height--distance from p4 to the face (p1,p2,p3)
  float height = fabs(nx*p4.x+ny*p4.y+nz*p4.z - faceConstant);
  //6. Get the bottom surface (p1,p2,p3) area
  float area = ComputeAreaFrom3Points(p1,p2,p3);
  //7. return the volume
  return (height*area)/3.0f;
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


void PopulateTetraFromMeshFile(const char* filename, Tetra **p_tetra_mesh, Point **p_points, TriNode **p_trinodes, TetraNode **p_tetranodes, int *p_Np, int *p_Nt)
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
  *p_points = points;
  fgets(line, 64, pFile);
  sscanf(line, "%d", p_Nt);
  for(i=1; i<*p_Np+1; i++)
  {
    fgets(line, 64, pFile);
    sscanf(line, "%f %f %f", &(points[i].x), &(points[i].y), &(points[i].z));
  }
  
  *p_tetra_mesh = (Tetra*)aocl_utils::alignedMalloc(sizeof(Tetra)*(*p_Nt+1));
  *p_trinodes = (TriNode*)malloc(sizeof(TriNode)*4*(*p_Nt));
  *p_tetranodes = (TetraNode*)malloc(sizeof(TetraNode)*(*p_Nt+1));
  for(i=0;i<4*(*p_Nt);i++)
  {
    (*p_trinodes)[i].N0 = 0;
    (*p_trinodes)[i].N1 = 0;
    (*p_trinodes)[i].N2 = 0;
    (*p_trinodes)[i].area = 0;
    (*p_trinodes)[i].fluence = 0;
  }
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
    (*p_tetranodes)[i].N0 = pointIDs[0];
    (*p_tetranodes)[i].N1 = pointIDs[1];
    (*p_tetranodes)[i].N2 = pointIDs[2];
    (*p_tetranodes)[i].N3 = pointIDs[3];
    (*p_tetranodes)[i].volume = ComputeVolume(points[pointIDs[0]], points[pointIDs[1]], points[pointIDs[2]], points[pointIDs[3]]);
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
      (*p_trinodes)[4*(firstNode.TetraID-1)+j].N0 = i;
      (*p_trinodes)[4*(firstNode.TetraID-1)+j].N1 = firstNode.lowerPointID;
      (*p_trinodes)[4*(firstNode.TetraID-1)+j].N2 = firstNode.higherPointID;
      (*p_trinodes)[4*(firstNode.TetraID-1)+j].area = ComputeAreaFrom3Points(p1, p2, p3);
    }
  }
  delete[] faceToTetraMap;

  return;
}

void ParseMaterial(const char* filename, Material** p_mats, int *p_Nm) 
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
  sscanf (line1, "%d", p_Nm);
  *p_mats = (Material*)aocl_utils::alignedMalloc(((*p_Nm)+1)*sizeof(Material));
  if (*p_mats == NULL)
  {
    printf("malloc failed\n");
  }
  for (int i = 1; i <= *p_Nm; i++) 
  {
    if ( fgets (linemat, 100, pFile) == NULL ) 
    {
      printf ("Error Reading Line %d Of Material (.opt) File.\n", i);
      fclose (pFile);
      return;
    }
    Material &mat = (*p_mats)[i];
    float mu_a, mu_s;
    sscanf (linemat, "%f %f %f %f", &mu_a, &mu_s, &(mat.g), &(mat.n));
    mat.mu_a = mu_a;
    mat.mu_s = mu_s;
    mat.mu_as = mu_a + mu_s;
    mat.rmu_as = 1.0f/mat.mu_as;
    mat.HGCoeff1 = (1+mat.g * mat.g)/(2*mat.g);
    mat.HGCoeff2 = (1-mat.g * mat.g)*(1-mat.g * mat.g) / (2*mat.g);
    mat.absfrac = 1 - mu_s / (mu_s+mu_a);
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
    (*p_mats)[0].setMatched = 1;

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

  (*p_mats)[0].g = 0;
  (*p_mats)[0].n = n_e;

  fclose (pFile);
}

int GetTetraIDFromCoordinates(const Tetra *tetra_mesh, const int Nt, const float x, const float y, const float z)
{
  for(int i = 1; i <= Nt; i++)
  {
    const Tetra &tetra = tetra_mesh[i];
    if(tetra.face[0][0]*x+tetra.face[0][1]*y+tetra.face[0][2]*z-tetra.face[0][3]<0)
    {
      continue;
    }
    if(tetra.face[1][0]*x+tetra.face[1][1]*y+tetra.face[1][2]*z-tetra.face[1][3]<0)
    {
      continue;
    }
    if(tetra.face[2][0]*x+tetra.face[2][1]*y+tetra.face[2][2]*z-tetra.face[2][3]<0)
    {
      continue;
    }
    if(tetra.face[3][0]*x+tetra.face[3][1]*y+tetra.face[3][2]*z-tetra.face[3][3]<0)
    {
      continue;
    }
    return i;
  }
  printf("Error: The source position is not inside the tissue.\n");
  return -1;
}

void PopulateCoordinatesFromTetraID(const Point *points, const TetraNode *tetranodes, Source *p_src)
{
  int IDt = p_src->IDt;
  Point p0=points[tetranodes[IDt].N0];
  Point p1=points[tetranodes[IDt].N1];
  Point p2=points[tetranodes[IDt].N2];
  Point p3=points[tetranodes[IDt].N3];
  p_src->x=(p0.x+p1.x+p2.x+p3.x)/4.0;
  p_src->y=(p0.y+p1.y+p2.y+p3.y)/4.0;
  p_src->z=(p0.z+p1.z+p2.z+p3.z)/4.0;
}

//////////////////////////////////////////////////////////////////////////////
// Parse Source (.Source File)
//////////////////////////////////////////////////////////////////////////////

void ParseSource(const char* filename, Source** sourcepoint, const Tetra *tetra_mesh, const int Nt, const Point *points, const TetraNode *tetranodes)
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
      exit(-1);
      fclose (sourcefile);
      return;
    }

    sscanf (linesource, "%d ", &(sourcepoint[i]->stype) );
    switch(sourcepoint[i]->stype)
    {
      // Point Position -> Cartesian Coordinates
      case 1:
        sscanf (linesource, "%d %f %f %f %d", &(sourcepoint[i]->stype), &(sourcepoint[i]->x), &(sourcepoint[i]->y), &(sourcepoint[i]->z), &(sourcepoint[i]->Np) );
        sourcepoint[i]->IDt = GetTetraIDFromCoordinates(tetra_mesh, Nt, sourcepoint[i]->x, sourcepoint[i]->y, sourcepoint[i]->z);
        printf("TetraID = %d\n",sourcepoint[i]->IDt );
        break;

      // IDt position
      case 2:
        sscanf (linesource, "%d %d %d", &(sourcepoint[i]->stype), &(sourcepoint[i]->IDt), &(sourcepoint[i]->Np) );
        PopulateCoordinatesFromTetraID(points, tetranodes, sourcepoint[i]);
        break;

      // IDt position
      case 11:
        sscanf (linesource, "%d %d %f %f %f %f %f %f %d", &(sourcepoint[i]->stype), &(sourcepoint[i]->IDt), &(sourcepoint[i]->x), &(sourcepoint[i]->y), &(sourcepoint[i]->z), &(sourcepoint[i]->dx), &(sourcepoint[i]->dy), &(sourcepoint[i]->dz), &(sourcepoint[i]->Np) );
        if (sourcepoint[i]->IDt != GetTetraIDFromCoordinates(tetra_mesh, Nt, sourcepoint[i]->x, sourcepoint[i]->y, sourcepoint[i]->z))
        {
          printf("Error: the source coordinate is outside the tetra ID in the input file.\n");
          fclose(sourcefile);
          exit(-1);
        }
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

