#=============================================================================#
#                                                                             #
# Author: Vishal Vaibhav                                                      #
#                                                                             #
# Description: makefile for creating shared object libraries                  #
#                                                                             #
# History:                                                                    #
# 19 July 2013: created                                                       #
#                                                                             #
#=============================================================================#
SHELL = /bin/sh
CC = gcc
# Default locations where gcc searches for header files 
# /usr/local/include
# /usr/include
INC_DIR = ../INCLUDE
BUILDDIR = ../BUILD

# Implicit Rules: Some rules are known to "make". They need not be stated
# explicitly, for example, how to create .o files from .c files.
# For compiling C programs, following implicit rule is used 
# $(CC) -c $(CFLAGS)
# Therefore, CFLAG must be defined correctly.

# compiler options:
# -fPIC : generate position independent code
# -fexceptions: enable exception handling
# -fno-omit-frame-pointer: 

# MKL options
MKL_PATH   = "$(MKLROOT)/lib/intel64"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/intel64"
EXT = so

# OpenMP library : OMP_LIB for parallelism in MKL
IFACE_LIB      = $(MKL_PATH)/libmkl_intel_lp64.$(EXT)
THREADING_LIB  = $(MKL_PATH)/libmkl_intel_thread.$(EXT)
SEQUENTIAL_LIB = $(MKL_PATH)/libmkl_sequential.$(EXT)
OMP_LIB        = -L$(CMPLR_PATH) -liomp5
CORE_LIB       = $(MKL_PATH)/libmkl_core.$(EXT)

MKL_LIBS = -Wl,--start-group $(IFACE_LIB) $(THREADING_LIB) $(CORE_LIB) -Wl,--end-group $(OMP_LIB)
#MKL_LIBS = -Wl,--start-group $(IFACE_LIB) $(SEQUENTIAL_LIB) $(CORE_LIB) -Wl,--end-group

INCLUDES = -I$(INC_DIR) -I$(MKLROOT)/include

CFLAGS = -ansi -D_GNU_SOURCE -fexceptions -fPIC  -O2 \
         -fno-omit-frame-pointer -pthread -std=c99  \
         $(INCLUDES) 

LDFLAGS = -pthread -shared -lm      \
          -L$(MKL_PATH) $(MKL_LIBS)

# Default locations where gcc searches for libraries 
# /usr/local/lib
# /usr/lib

TARGET = $(BUILDDIR)/libfem1d.so

SOURCES = fem1d.c    \
          legendre.c \
          jacobi.c   \
          matlib.c   \
          matlib_io.c     \
          matlib_solver.c \
          pthpool.c  \
          pfem1d.c

OBJECTS = $(SOURCES:%.c=%.o)
# HEADERS = $(SOURCES:%.c=$(INC_DIR)/%.h)

HEADERS = $(INC_DIR)/fem1d.h    \
	  $(INC_DIR)/legendre.h \
	  $(INC_DIR)/jacobi.h   \
	  $(INC_DIR)/matlib.h   \
	  $(INC_DIR)/pthpool.h  \
	  $(INC_DIR)/pfem1d.h   \
	  $(INC_DIR)/ehandler.h \
	  $(INC_DIR)/debug.h


.PHONY: all clean

.SECONDARY: $(OBJECTS)

all : $(TARGET) 

$(TARGET): $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJECTS)


clean :
	rm -f $(OBJECTS)

