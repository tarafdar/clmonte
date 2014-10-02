CLSOURCES=main.cpp io.cpp mem.cpp 
CLOUT= gpumcml_opencl
CLFLAGS= -O3 -g -o
CLLIB=-l OpenCL
all: $(CLOUT) 

$(CLOUT): $(CLSOURCES) gpumcml.h 
	g++ $(CLSOURCES) $(CLFLAGS) $(CLOUT) $(CLLIB)

