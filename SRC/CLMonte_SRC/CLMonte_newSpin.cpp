#include <CL/cl.h>
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "CLMonte.h"

unsigned int xtest[NUM_THREADS];
unsigned int ctest[NUM_THREADS];
unsigned int atest[NUM_THREADS];


// forward declaration of the device code

void check_return(cl_int ret, const char* message){
    
    if(ret!=CL_SUCCESS){
        printf("%d %s\n", ret,message);
        exit(-1);
    }

}
//#include "CLMonte_goldstandard.c"
// wrapper for device code
int MC(unsigned int* x,unsigned int* c,unsigned int* a)
{
	unsigned long long num_tot, hist_tot, num_l_tot;
	unsigned int i, j;
#ifdef SPATIAL_HISTOGRAM	
    unsigned int hist[RBUCKETSXTEMP];
    unsigned int hist_transpose[RBUCKETSXTEMP];
	unsigned int histh[RBUCKETSXTEMP];
#else
	unsigned int histh[TEMP];
	unsigned int hist[TEMP];
#endif
#ifdef EVENT_LOGGING
	unsigned int num[NUM_THREADS];
    unsigned int num_launched[NUM_THREADS]; 
#endif    
	unsigned int numh[NUM_THREADS];
	
    FILE *fp;
    char *source_str;
    size_t source_size;



#ifdef JORDAN
    fp = fopen("C:\\Users\\Jordan\\Documents\\GitHub\\ece496\\CLMonteTransport.cl", "r");
#endif
#ifdef LINUX
    //fp = fopen("SRC/CLMonte_SRC/CLMonteTransport.cl", "r");
    //fp = fopen("SRC/CLMonte_SRC/CLMonteTransport_n1.cl", "r");
    fp = fopen("SRC/CLMonte_SRC/CLMonteTransport_newSpin.cl", "r");

#endif
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source_str = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
    fclose( fp );


	FILE *file; 
	FILE *print_file;
	FILE *build_file;
	
	file = fopen("outp_newSpin.txt", "w");
	print_file = fopen("print_out.txt", "w");
	build_file = fopen("build.txt", "w");

	clock_t time1,time2,GPUtime,CPUtime;
    int size;
    size=NUM_THREADS*sizeof(unsigned int);
	
    
	// Get platform and device information
    cl_platform_id platform_id = NULL;
    cl_device_id device_id = NULL;   
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_GPU, 1, 
            &device_id, &ret_num_devices);


#ifdef INPUT_PARAMETERS
    float G = 0.9f;	
    float MUS_MAX = 90.0f;	//[1/cm]
    float V = 0.0214f;		//[cm/ps] (c=0.03 [cm/ps] v=c/n) here n=1.4
    float COS_CRIT = 0.6999f;	//the critical angle for total internal reflection at the border cos_crit=sqrt(1-(nt/ni)^2)
    float N = 1.4f;
#endif



	// Create an OpenCL context
    cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);

    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

	// Create memory buffers on the device for each vector 
    cl_mem xd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create xd buff fail\n");

    
    cl_mem cd_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create cd buff fail\n");

    cl_mem ad_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
    check_return(ret, "create ad buff fail\n");
#ifdef EVENT_LOGGING
    cl_mem numd_mem_obj = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size, NULL, &ret);
    check_return(ret, "create numd buff fail\n");
    
    cl_mem numl_mem_obj = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size, NULL, &ret);
    check_return(ret, "create numl buff fail\n");
#endif

#ifdef SPATIAL_HISTOGRAM
    cl_mem histd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, RBUCKETSXTEMP*sizeof(unsigned int), NULL, &ret);
#else
    cl_mem histd_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, TEMP*sizeof(unsigned int), NULL, &ret);
#endif
    
    check_return(ret, "create histd buff fail\n");
	
    ret = clEnqueueWriteBuffer(command_queue, xd_mem_obj, CL_TRUE, 0, size, xtest, 0, NULL, NULL);
    check_return(ret, "write buff xd fail\n");
	
    ret = clEnqueueWriteBuffer(command_queue, ad_mem_obj, CL_TRUE, 0, size, atest, 0, NULL, NULL);
    check_return(ret, "write buff ad fail\n");

	ret = clEnqueueWriteBuffer(command_queue, cd_mem_obj, CL_TRUE, 0, size, ctest, 0, NULL, NULL);
    check_return(ret, "write buff cd fail\n");
	
#ifdef SPATIAL_HISTOGRAM
	for(i=0;i<RBUCKETSXTEMP;i++)hist[i]=0;
	ret = clEnqueueWriteBuffer(command_queue, histd_mem_obj, CL_TRUE, 0, RBUCKETSXTEMP*sizeof(unsigned int), hist, 0, NULL, NULL);
#else
	for(i=0;i<TEMP;i++)hist[i]=0;
	ret = clEnqueueWriteBuffer(command_queue, histd_mem_obj, CL_TRUE, 0, TEMP*sizeof(unsigned int), hist, 0, NULL, NULL);
#endif
    check_return(ret, "write buff histd fail\n");

	// Create a program from the kernel source
    cl_program program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
    check_return(ret, "create program fail\n");

    // Build the program
//    ret = clBuildProgram(program, 1, &device_id, "-g -s \"C:\\Users\\Naif\\Documents\\Visual Studio 2010\\Projects\\CudaMC\\CudaMC\\CUDAMC\\CUDAMCtransport.cl\"", NULL, NULL);

	ret = clBuildProgram(program, 1, &device_id, "-I SRC/CLMonte_SRC/", NULL, NULL);


	if(ret!=CL_SUCCESS){
		char *build_log;
		size_t ret_val_size;
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size);
		build_log = new char[ret_val_size+1];
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL);
		build_log[ret_val_size] = '\0';
		fputs (build_log, build_file);
		delete build_log;
		fclose(build_file);
		printf("Build Program fail ad\n");
		return -1;
	}

    // Create the OpenCL kernel
    cl_kernel kernel = clCreateKernel(program, "MCd", &ret);
    check_return(ret, "create kernel MCd fail\n");

    int argval = 0;

    // Set the arguments of the kernel
    ret = clSetKernelArg(kernel, argval++, sizeof(cl_mem), (void *)&xd_mem_obj);
    check_return(ret, "set arg fail\n");

	ret = clSetKernelArg(kernel, argval++, sizeof(cl_mem), (void *)&cd_mem_obj);
    check_return(ret, "set arg fail\n");
    
	ret = clSetKernelArg(kernel, argval++, sizeof(cl_mem), (void *)&ad_mem_obj);
    check_return(ret, "set arg fail\n");
    
	ret = clSetKernelArg(kernel, argval++, sizeof(cl_mem), (void *)&histd_mem_obj);
    check_return(ret, "set arg fail\n");
#ifdef INPUT_PARAMETERS
	ret = clSetKernelArg(kernel, argval++, sizeof(float), (void *)&N);
    check_return(ret, "set arg fail\n");

	ret = clSetKernelArg(kernel, argval++, sizeof(float), (void *)&COS_CRIT);
    check_return(ret, "set arg fail\n");


	ret = clSetKernelArg(kernel, argval++, sizeof(float), (void *)&V);
    check_return(ret, "set arg fail\n");

	ret = clSetKernelArg(kernel, argval++, sizeof(float), (void *)&MUS_MAX);
    check_return(ret, "set arg fail\n");
	
    ret = clSetKernelArg(kernel, argval++, sizeof(float), (void *)&G);
    check_return(ret, "set arg fail\n");
#endif
#ifdef EVENT_LOGGING    
	ret = clSetKernelArg(kernel, argval++, sizeof(cl_mem), (void *)&numd_mem_obj);
    check_return(ret, "set arg fail\n");
	
	ret = clSetKernelArg(kernel, argval++, sizeof(cl_mem), (void *)&numl_mem_obj);
    check_return(ret, "set arg fail\n");
#endif
    
    // Execute the OpenCL kernel on the list
	cl_event kern_event;
    size_t global_item_size = NUM_THREADS; // Process the entire lists
    size_t local_item_size = NUM_THREADS_PER_BLOCK; // Process in groups of NUM_THREADS_PER_BLOCk
    printf("Running code on Accelerator in OpenCL\n");
	fprintf(print_file, "Running code on Accelerator in OpenCL\n");
	ret = clFlush(command_queue);
    ret = clFinish(command_queue);
	time1=clock();
	ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, &kern_event);
    check_return(ret, "enqueuendrange kernel fail\n");
	
	clWaitForEvents(1,&kern_event);
	
	time2=clock();

 #ifdef SPATIAL_HISTOGRAM
	ret = clEnqueueReadBuffer(command_queue, histd_mem_obj, CL_TRUE, 0,RBUCKETSXTEMP*sizeof(unsigned int), hist, 1, &kern_event, NULL);
#else	
    ret = clEnqueueReadBuffer(command_queue, histd_mem_obj, CL_TRUE, 0,TEMP*sizeof(unsigned int), hist, 1, &kern_event, NULL);
#endif
    check_return(ret, "read buff hist fail\n");

#ifdef EVENT_LOGGING
	ret = clEnqueueReadBuffer(command_queue, numd_mem_obj, CL_TRUE, 0,size, num, 1, &kern_event, NULL);
    check_return(ret, "read buff num fail\n");
	
	ret = clEnqueueReadBuffer(command_queue, numl_mem_obj, CL_TRUE, 0,size, num_launched, 1, &kern_event, NULL);
    check_return(ret, "read buff numl fail\n");

	num_tot=0;
	for(i=0;i<NUM_THREADS;i++){
		num_tot+=num[i];
	}
	
    num_l_tot=0;
    for(i=0;i<NUM_THREADS;i++){
        num_l_tot+=num_launched[i];
    }
#endif	
	

 #ifdef SPATIAL_HISTOGRAM
	hist_tot=0;
	for(i=0;i<RBUCKETSXTEMP;i++)hist_tot+=hist[i];
    for(i=0;i<RBUCKETS;i++){
	    printf("RADIUS BUCKET %d\n", i);
		for(j=0; j<TEMP; j++){
            printf("%d ",hist[i*TEMP + j]);
		    fprintf(print_file, "%d ", hist[i*TEMP + j]);
        }
        printf("\n\n");
        fprintf(print_file, "\n");
	}

	for(i=0;i<RBUCKETS;i++){
        for(j=0; j<TEMP; j++)
            //fprintf(file,"%d %d\n",j, hist[i*TEMP + j]);
            fprintf(file,"%d ", hist[i*TEMP + j]);
            fprintf(file, "\n");
    }

//TRANSPOSE
    for(i=0;i<RBUCKETS;i++)
        for(j=0; j<TEMP; j++)
            hist_transpose[j*RBUCKETS + i] = hist[i*TEMP + j];
#else	
	hist_tot=0;
	for(i=0;i<TEMP;i++)hist_tot+=hist[i];
    for(i=0;i<TEMP;i++){
		printf("%d ",hist[i]);
		fprintf(print_file, "%d ", hist[i]);
	}

   for(i=0; i<TEMP; i++)fprintf(file,"%d %d\n",i, hist[i]);
#endif
	fclose(file);
#ifdef EVENT_LOGGING
	printf("\nTotal number of photons terminated (i.e. full path simulated): %llu\n",num_tot);
    printf("Total number of photons launched: %llu\n",num_l_tot);
    printf("Average number of scattering events per photon: %e\n", ((double)NUM_THREADS*(double)NUMSTEPS_GPU)/num_l_tot);
#endif    
    printf("Number of photons contribution to the histogram: %llu\n",hist_tot);
	printf("Total number of photons steps simulated: %e\n",(double)NUM_THREADS*(double)NUMSTEPS_GPU);
    printf("time1=%u, time2=%u.\n",time1,time2);

	printf("Photon steps per sec: %e\n",((double)NUM_THREADS*(double)NUMSTEPS_GPU)/((double(time2-time1))/CLOCKS_PER_SEC));

#ifdef EVENT_LOGGING
	fprintf(print_file, "\nTotal number of photons terminated (i.e. full path simulated): %llu\n", num_tot);
    fprintf(print_file, "Total number of photons launched: %llu\n",num_l_tot);
    fprintf(print_file, "Average number of scattering events per photon: %e\n", ((double)NUM_THREADS*(double)NUMSTEPS_GPU)/num_l_tot);
#endif    
    fprintf(print_file,"Number of photons contribution to the histogram: %llu\n",hist_tot);
	fprintf(print_file, "Total number of photons steps simulated: %e\n",(double)NUM_THREADS*(double)NUMSTEPS_GPU);
    fprintf(print_file, "time1=%u, time2=%u.\n",time1,time2);

	fprintf(print_file, "Photon steps per sec: %e\n",((double)NUM_THREADS*(double)NUMSTEPS_GPU)/(((double(time2-time1)))/CLOCKS_PER_SEC));
	
	GPUtime=time2-time1;
	
	ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(xd_mem_obj);
	ret = clReleaseMemObject(cd_mem_obj);
	ret = clReleaseMemObject(ad_mem_obj);
#ifdef EVENT_LOGGING
	ret = clReleaseMemObject(numd_mem_obj);
	ret = clReleaseMemObject(numl_mem_obj);
#endif	
    ret = clReleaseMemObject(histd_mem_obj);
	ret = clReleaseCommandQueue(command_queue);

	ret = clReleaseContext(context);
#ifdef RUN_HOST
	printf("\n\nRunning CPU code (sequential C++)\n");
	fprintf(print_file, "\n\nRunning CPU code (sequential C++)\n");
	
#ifdef SPATIAL_HISTOGRAM
    for(i=0;i<RBUCKETSXTEMP;i++)histh[i]=0;
#else
    for(i=0;i<TEMP;i++)histh[i]=0;
#endif
    
	//run CPU code 
	time1=clock();
	MCh(xtest,ctest,atest,numh,histh, GPUtime);
    time2=clock();

	num_tot=0;
	for(i=0;i<NUM_THREADS;i++)num_tot+=numh[i];

	hist_tot=0;
#ifdef SPATIAL_HISTOGRAM
	for(i=0;i<RBUCKETSXTEMP;i++)hist_tot+=histh[i];
#else
    for(i=0;i<TEMP;i++)histh[i]=0;
#endif


file = fopen("outph.txt", "w");
#ifdef SPATIAL_HISTOGRAM
	for(i=0;i<RBUCKETS;i++){
        for(j=0; j<TEMP; j++){
		    printf("%d ", histh[i*TEMP + j]);
            fprintf(file,"%d ", histh[i*TEMP + j]);
            fprintf(print_file, "%d ", histh[i*TEMP + j]);
        }
        fprintf(file,"\n");
	}
#else

	for(i=0;i<TEMP;i++){
		printf("%d ", histh[i]);
		fprintf(file,"%d %d\n",i,histh[i]);
		fprintf(print_file, "%d ", histh[i]);
	}
#endif
	fclose(file);
	printf("\n\nTotal number of photons (i.e. full path simulated): %llu\nNumber of photons contribution to the histogram: %llu\n",num_tot,hist_tot);
	printf("Total number of photons steps simulated: %e\n",(double)NUM_THREADS*(double)NUMSTEPS_CPU);
    printf("time1=%u, time2=%u.\n",time1,time2);

	printf("Photon steps per sec: %e\n",((double)NUM_THREADS*(double)NUMSTEPS_CPU)/(double(time2-time1)/CLOCKS_PER_SEC));
	
	fprintf(print_file, "\n\nTotal number of photons (i.e. full path simulated): %llu\nNumber of photons contribution to the histogram: %llu\n",num_tot,hist_tot);
	fprintf(print_file, "Total number of photons steps simulated: %e\n",(double)NUM_THREADS*(double)NUMSTEPS_CPU);
    fprintf(print_file, "time1=%u, time2=%u.\n",time1,time2);

	fprintf(print_file, "Photon steps per sec: %e\n",((double)NUM_THREADS*(double)NUMSTEPS_CPU)/(double(time2-time1)/CLOCKS_PER_SEC));

	CPUtime=time2-time1;

	printf("\n\nSpeedup: %f\n",(NUMSTEPS_GPU*double(CPUtime))/(NUMSTEPS_CPU*double(GPUtime)));
	fprintf(print_file, "\n\nSpeedup: %f\n",(NUMSTEPS_GPU*double(CPUtime))/(NUMSTEPS_CPU*double(GPUtime)));
#endif	
    fclose(print_file);
}




void initialize(void)//Straight from Steven Gratton's code
{
    FILE *fp;
    unsigned int begin=0u;
    unsigned long long int xinit=1ull;
    unsigned int cinit=0u;
    unsigned int fora,tmp1,tmp2;
#ifdef JORDAN
    fp=fopen("C:\\Users\\Jordan\\Documents\\GitHub\\ece496\\safeprimes_base32.txt","r");//use an expanded list containing 50000 safeprimes instead of Steven's shorter list
#endif
#ifdef LINUX
	fp=fopen("SRC/HelperFiles/safeprimes_base32.txt","r");//use an expanded list containing 50000 safeprimes instead of Steven's shorter list
#endif

// use begin as a multiplier to generate the initial x's for 
// the other generators...
	if(fp!=NULL)
		fscanf(fp,"%u %u %u",&begin,&tmp1,&tmp2);
	else{
		printf("bro null\n");
		exit(1);
	}



    for (unsigned int i=0;i<NUM_THREADS;i++)
    {

	xinit=xinit*begin+cinit;
	cinit=xinit>>32;
	xinit=xinit&0xffffffffull;
	xtest[i]=(unsigned int) xinit;
	fscanf(fp,"%u %u %u",&fora,&tmp1,&tmp2);
	atest[i]=fora;

	xinit=xinit*begin+cinit;
	cinit=xinit>>32;
	xinit=xinit&0xffffffffull;
	ctest[i]=(unsigned int) ((((double)xinit)/UINT_MAX)*fora);

    }
    fclose(fp);
}



int main(int argc,char* argv[])
{
	//do all the initialization for the RNG's (one MWC per thread)
    initialize();
    MC(xtest,ctest,atest);
	return 0;
    
}
