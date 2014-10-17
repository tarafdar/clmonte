#define NUM_STEPS 50000
#define NUM_BLOCKS 30
#define NUM_THREADS_PER_BLOCK 512
#define NUM_THREADS (NUM_BLOCKS * NUM_THREADS_PER_BLOCK)

#include <stdio.h>

int main()
{
  int num_photons;
  float time;
  float timeSum;
  float rate;
  int repeatTime = 5;
  FILE* fp_read = fopen("photon2time.csv", "r");
  FILE* fp_write = fopen("photon2AvgTime.csv", "w");
  int i;
  while(fscanf(fp_read,"%d,%f,%e",&num_photons,&timeSum,&rate) != EOF)
  {
    int j;
    for(j=0;j<repeatTime-1;j++)
    {
      fscanf(fp_read, "%d,%f,%e",&num_photons,&time,&rate);
      timeSum = timeSum+time;
    }
    float avgTime = timeSum/repeatTime;
    
    fprintf(fp_write, "%d,%f,%e\n", num_photons, avgTime, NUM_STEPS*NUM_THREADS/avgTime);
  }
  printf("%d,%f,%e\n",num_photons,time,rate);
}
