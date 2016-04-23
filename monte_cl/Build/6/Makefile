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
#ifeq ($(wildcard $(ALTERAOCLSDKROOT)),)
#$(error Set ALTERAOCLSDKROOT to the root directory of the Altera SDK for OpenCL software installation)
#endif
#ifeq ($(wildcard $(ALTERAOCLSDKROOT)/host/include/CL/opencl.h),)
#$(error Set ALTERAOCLSDKROOT to the root directory of the Altera SDK for OpenCL software installation.)
#endif

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
INCS := $(wildcard $(MONTE_SRC)/h/*.h monte_cl/host/h/*.h)
SRCS := $(wildcard $(MONTE_SRC)/host/src/*.cpp $(MONTE_SRC)/host/common/src/AOCLUtils/*.cpp)
KERNEL_SRCS := $(MONTE_SRC)/device/monte.cl
AOCO := $(TARGET_DEVICE_DIR)/monte.aoco
LIBS := 

BOARD := p385_hpc_d5

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
	aoc -v -c -g --report --board $(BOARD) -I $(KERNEL_INC_DIRS) $(KERNEL_SRCS) -o $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L0)

$(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L1) : $(KERNEL_SRCS)
	echo "Building kernel for emulator"
	aoc -march=emulator -v -g --board $(BOARD) -I $(KERNEL_INC_DIRS) $(KERNEL_SRCS) -o $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L1)

$(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L2) : $(KERNEL_SRCS)
	echo "Building kernel for FPGA"
	aoc -v --profile --board $(BOARD) -I $(KERNEL_INC_DIRS) $(KERNEL_SRCS) -o $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L2)
	#aoc -v --profile --board $(BOARD) $(AOCO) -o $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L2)

$(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L3) : $(KERNEL_SRCS)
	echo "Building Full Deployment"
	aoc -v --board $(BOARD) -I $(KERNEL_INC_DIRS) $(KERNEL_SRCS) -o $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L3)
	#aoc -v --board $(BOARD) $(AOCO) -o $(TARGET_DEVICE_DIR)/$(TARGET_KERNEL_L3)

# Standard make targets
clean :
	echo "Clean All"
	rm -rf bin
	mkdir bin
	mkdir bin/host
	mkdir bin/device

.PHONY : all clean
