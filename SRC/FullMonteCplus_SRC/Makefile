# base options only
GCC_OPTS=-Wall -mfpmath=sse -Wstrict-aliasing=2 -g -msse4 -DSVNVERSION=\"$(VERSION)\" -DSSE -mtune=core2 -march=core2 -fpermissive -std=gnu++0x -Wno-narrowing -fpermissive
LIBS=-lboost_program_options -lboost_timer -lcrypto -lboost_system
LIBDIRS=-L/usr/local/lib -L/usr/local/lib/boost
INCLDIRS=-I/usr/local/boost -I/usr/local/include -I. -I/usr/local/include/boost -I/usr/local/include/pgsql

# Switch debug on/off
#GCC_OPTS += -O0
GCC_OPTS += -DNDEBUG -O3

# Enable profiling
#GCC_OPTS += -pg

# Comment out if you want a count of how often each tetra/face is hit
GCC_OPTS += -DNOHIT

# Check for Mac OS (no POSIX timer)
OS:=$(shell uname)
VERSION:=$(shell svnversion -n 2>/dev/null )
ifeq ($(OS),Darwin)
GCC_OPTS += -DPLATFORM_DARWIN
endif


all: fullmonte

random.o: random.cpp random.hpp
	g++ -O1 -msse4 -g -Wall -msse3 -DNDEBUG -DPLATFORM_DARWIN $< -c -o $@

SFMT.o: SFMT.c SFMT*.h
	gcc -O1 -c -g -Wall -Wmissing-prototypes -Winline -finline-functions -fno-strict-aliasing --param max-inline-insns-single=1800 -msse2 -DHAVE_SSE2 -DSFMT_MEXP=19937 -DNDEBUG $< -o $@

%.o: %.cpp *.hpp
	g++ $(GCC_OPTS) $(INCLDIRS) -c $*.cpp -o $@

fullmonte: graph.o newgeom.o face.o helpers.o source.o fullmonte.o LoggerSurface.o io_timos.o progress.o linefile.o fluencemap.o mainloop.o sse.o random.o SFMT.o LoggerConservation.o LoggerEvent.o LoggerVolume.o
	g++ $(GCC_OPTS) $^ $(LIBS) $(LIBDIRS) -o $@

clean:
	rm -f *.o sse_int montecarlo compareflight blobmaint texwriter fullmonte
