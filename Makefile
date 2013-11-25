CLSOURCES=$(CLDIR)/CLMonte.cpp
CLSOURCES_newSpin=$(CLDIR)/CLMonte_newSpin.cpp
CUDASOURCES=$(CUDADIR)/CUDAMC.cu
CLOUT=CLMonte
CLOUT_newSpin=CLMonte_newSpin
CUDAOUT=CUDAMC
CLFLAGS =-O3 -o
CUDAFLAGS=-O3 -o
CLLIB=-l OpenCL
CUDAARCH=-arch=sm_11
SRCDIR=SRC
CUDADIR=$(SRCDIR)/CUDA_SRC
CLDIR=$(SRCDIR)/CLMonte_SRC
CLDIR_UNITTESTS=$(SRCDIR)/CLMonte_SRC/UnitTests
UNITTESTS = unitTests
UTSOURCES=$(CLDIR_UNITTESTS)/CLMonteUT.cpp

all: $(CLOUT) $(CUDAOUT)

$(CLOUT): $(CLDIR)/CLMonte.cpp $(CLDIR)/CLMonteTransport.cl $(CLDIR)/CLMonte_goldstandard.c
	g++ $(CLSOURCES) $(CLFLAGS) $(CLOUT) $(CLLIB)

$(CLOUT_newSpin): $(CLDIR)/CLMonte_newSpin.cpp $(CLDIR)/CLMonteTransport_newSpin.cl $(CLDIR)/CLMonte_goldstandard.c
	g++ $(CLSOURCES_newSpin) $(CLFLAGS) $(CLOUT_newSpin) $(CLLIB)

$(CUDAOUT): $(CUDADIR)/CUDAMC.cu $(CUDADIR)/CUDAMCtransport.cu $(CUDADIR)/CUDAMC_goldstandard.c
	nvcc $(CUDAARCH)  $(CUDASOURCES) $(CUDAFLAGS) $(CUDAOUT)

$(UNITTESTS): $(CLDIR_UNITTESTS)/CLMonteUT.cpp $(CLDIR_UNITTESTS)/CLMonteUT.cl
	g++ $(UTSOURCES) $(CLFLAGS) $(UNITTESTS) $(CLLIB)

clean:
	rm -rf *.txt $(CLOUT) $(CUDAOUT) $(UNITTESTS) $(CLOUT_newSpin)
