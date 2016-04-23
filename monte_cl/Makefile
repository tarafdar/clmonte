# Copyright (C) 2013-2015 Altera Corporation, San Jose, California, USA. All rights reserved.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to
# whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
# 
# This agreement shall be governed in all respects by the laws of the State of California and
# by the laws of the United States of America.
# This is a GNU Makefile.

# You must configure ALTERAOCLSDKROOT to point the root directory of the Altera SDK for OpenCL
# software installation.
# See http://www.altera.com/literature/hb/opencl-sdk/aocl_getting_started.pdf 
# for more information on installing and configuring the Altera SDK for OpenCL.

ifeq ($(VERBOSE),1)
ECHO := 
else
ECHO := @
endif

# Where is the Altera SDK for OpenCL software?
ifeq ($(wildcard $(ALTERAOCLSDKROOT)),)
$(error Set ALTERAOCLSDKROOT to the root directory of the Altera SDK for OpenCL software installation)
endif
ifeq ($(wildcard $(ALTERAOCLSDKROOT)/host/include/CL/opencl.h),)
$(error Set ALTERAOCLSDKROOT to the root directory of the Altera SDK for OpenCL software installation.)
endif

ifeq ($(MONTE_ROOT),)
$(error ARC Environment is not set)
endif

ifeq ($(MONTE_BIN),)
$(error ARC Environment is not set)
endif

ifeq ($(MONTE_SRC),)
$error ARC Environment is not set)
endif

# OpenCL compile and link flags.
AOCL_COMPILE_CONFIG := $(shell aocl compile-config )
AOCL_LINK_CONFIG := $(shell aocl link-config )

# Compilation flags
ifeq ($(DEBUG),1)
CXXFLAGS += -g
else
CXXFLAGS += -O2
endif

# Compiler
CXX := g++

# Target
TARGET := montecl
TARGET_KERNEL_L0 := monte.aoco
TARGET_KERNEL_L1 := monte_emulated.aocx
TARGET_KERNEL_L2 := monte_profiled.aocx
TARGET_KERNEL_L3 := monte_deployed.aocx
TARGET_HOST_DIR := $(MONTE_BIN)/host
TARGET_DEVICE_DIR := $(MONTE_BIN)/device

# Directories
KERNEL_INC_DIRS := $(MONTE_SRC)/h
INC_DIRS := $(MONTE_SRC) $(MONTE_SRC)/host/common/inc # $(MONTE_SRC)/host/h $(MONTE_SRC)/h
LIB_DIRS := 

# Files
INCS := $(wildcard $(MONTE_SRC)/h/*.h $(MONTE_SRC)/host/h/*.h)
SRCS := $(wildcard $(MONTE_SRC)/host/src/*.cpp $(MONTE_SRC)/host/common/src/AOCLUtils/*.cpp)
KERNEL_SRCS := $(MONTE_SRC)/device/monte.cl
AOCO := $(TARGET_DEVICE_DIR)/monte.aoco
LIBS := 

BOARD := p385_hpc_d5

DEBUG := 0
SELF_TEST := 0

# Set from command line
ifeq ($(DEBUG),1)
DEFINE_DEBUG := -D DEBUG=1
else
DEFINE_DEBUG :=
endif

ifeq ($(SELF_TEST),1)
DEFINE_SELF_TEST := -D SELF_TEST=1
else
DEFINE_SELF_TEST :=
endif

ifdef LAPS
DEFINE_LAPS := -D LAPS=$(LAPS)
else
DEFINE_LAPS :=
endif

ifeq ($(CACHE_TEST),1)
DEFINE_CACHE_TEST := -D CACHE_TEST=1
else
DEFINE_CACHE_TEST :=
endif

ifeq ($(USE_CACHE),1)
DEFINE_USE_CACHE := -D USE_CACHE=1
else
DEFINE_USE_CACHE :=
endif

DEFINES := $(DEFINE_SELF_TEST) $(DEFINE_DEBUG) $(DEFINE_LAPS) $(DEFINE_CACHE_TEST) $(DEFINE_USE_CACHE)

# Make it all!
all : host kernel

# Host executable target.
host : $(TARGET_HOST_DIR)/$(TARGET)

$(TARGET_HOST_DIR)/$(TARGET) : Makefile $(SRCS) $(INCS)
	echo "Building Host"
	@[ -d $(TARGET_HOST_DIR) ] || mkdir $(TARGET_HOST_DIR)
	$(ECHO)$(CXX) -DHOST $(CPPFLAGS) $(CXXFLAGS) -fPIC $(foreach D,$(INC_DIRS),-I$D) \
			$(AOCL_COMPILE_CONFIG) $(SRCS) $(AOCL_LINK_CONFIG) \
			$(foreach D,$(LIB_DIRS),-L$D) \
			$(foreach L,$(LIBS),-l$L) \
			-o $(TARGET_HOST_DIR)/$(TARGET)
	
kernel_l0 : $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L0)
kernel_l1 : $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L1)
kernel_l2 : $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L2)
kernel_l3 : $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L3)

$(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L0) : $(KERNEL_SRCS)
	echo "Building kernel L0"
	aoc -v -c -g --report $(DEFINES) --const-cache-bytes 1024 --board $(BOARD) -I $(KERNEL_INC_DIRS) $(KERNEL_SRCS) -o $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L0)

$(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L1) : $(KERNEL_SRCS)
	echo "Building kernel for Emulator"
	aoc -march=emulator -v -g $(DEFINES) --const-cache-bytes 1024 --board $(BOARD) -I $(KERNEL_INC_DIRS) $(KERNEL_SRCS) -o $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L1)

$(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L2) : $(KERNEL_SRCS)
	echo "Building kernel for Profiler"
	aoc -v --profile $(DEFINES) --const-cache-bytes 1024 --board $(BOARD) -I $(KERNEL_INC_DIRS) $(KERNEL_SRCS) -o $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L2)

$(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L3) : $(KERNEL_SRCS)
	echo "Building Full Deployment"
	#aoc -v --fp-relaxed --fpc $(DEFINES) --board $(BOARD) -I $(KERNEL_INC_DIRS) $(KERNEL_SRCS) -o $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L3)
	aoc -v --seed 5 $(DEFINES) --const-cache-bytes 1024 --board $(BOARD) -I $(KERNEL_INC_DIRS) $(KERNEL_SRCS) -o $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L3)

# Standard make targets
clean :
	echo "Clean All"
	rm -rf $(MONTE_BIN)
	mkdir $(MONTE_BIN)
	mkdir $(TARGET_HOST_DIR)
	mkdir $(TARGET_DEVICE_DIR)

clean_host :
	echo "Clean Host"
	rm -rf $(TARGET_HOST_DIR)
	mkdir $(TARGET_HOST_DIR)

clean_kernel_l0 :
	echo "Clean kernel_l0"
	rm -rf $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L0)

clean_kernel_l1 :
	echo "Clean kernel_l1"
	rm -rf $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L1)

clean_kernel_l3 :
	echo "Clean kernel_l3"
	rm -rf $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L3)


.PHONY : all clean
