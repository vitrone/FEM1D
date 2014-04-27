#=============================================================================#
#                                                                             #
# Author: Vishal Vaibhav                                                      #
#                                                                             #
# Description: makefile for creating test targets                             #
#                                                                             #
# History:                                                                    #
# 22 Feb 2014 : created                                                       #
#                                                                             #
#=============================================================================#
SHELL = /bin/sh

CC = gcc
#
# MKL options
MKL_PATH   = "$(MKLROOT)/lib/intel64"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/intel64"
EXT = so

IFACE_LIB     = $(MKL_PATH)/libmkl_intel_lp64.$(EXT)
THREADING_LIB = $(MKL_PATH)/libmkl_intel_thread.$(EXT)
OMP_LIB       = -L$(CMPLR_PATH) -liomp5
CORE_LIB      = $(MKL_PATH)/libmkl_core.$(EXT)

MKL_LIBS = -Wl,--start-group $(IFACE_LIB) $(THREADING_LIB) $(CORE_LIB) -Wl,--end-group $(OMP_LIB)

INCLUDES = -I../INCLUDE  -I$(MKLROOT)/include

# CFLAG is used for the implicit rule to generate object files for C
# $(CC) -c $(CFLAGS)
# use -save-temps to see the preprocess effects in .i files
CFLAGS = -ansi -D_GNU_SOURCE  -fexceptions -fPIC  -O2 \
         -fno-omit-frame-pointer -pthread -std=c99 \
         $(INCLUDES) 

LDFLAGS = -pthread -lm -L../BUILD -lfem1d -L$(MKL_PATH) $(MKL_LIBS) -lcunit
SOURCES = CUnit_matrix.c   \
          CUnit_legendre.c \
          CUnit_jacobi.c   \
          CUnit_fem1d.c    \
          CUnit_solver.c   \
          speedup_pfem1d.c \
          CUnit_fem1d_Schroedinger.c   \
          CUnit_pfem1d_Schroedinger.c  \
          CUnit_pthpool_func.c         \
          CUnit_pfem1d.c       

TSTBIN = BIN
OBJECTS = $(SOURCES:%.c=$(TSTBIN)/%.o)
TARGETS = $(SOURCES:%.c=$(TSTBIN)/%)


.PHONY: all clean


all : $(TARGETS)

$(TSTBIN)/% : $(TSTBIN)/%.o
	echo "Building binary files..."
	$(CC) $(CFLAGS) $(LDFLAGS) $(CLIBS) $< -o $@ 

$(TSTBIN)/%.o : %.c
	echo "Building object files..."
	$(CC) $(CFLAGS) -c $< -o $@ 


clean :
	rm -f $(OBJECTS)

