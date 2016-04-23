
ifeq ($(MONTE_TOOLS),)
$(error MONTE_TOOLS Environment not set)
endif

all : arc mcl_test

arc : $(MONTE_TOOLS)/h/utils.h $(MONTE_TOOLS)/arc_build/arc.cpp
	g++ arc_build/arc.cpp -I h -o arc

mcl_test : $(MONTE_TOOLS)/h/utils.h $(MONTE_TOOLS)/test_build/mcl_test.cpp
	g++ test_build/mcl_test.cpp -I h -o mcl_test
