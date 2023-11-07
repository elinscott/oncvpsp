# System-dependent makefile options for ONCVPSP
# This must be carefully edited before executing "make" in src
#
# Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
 
##### Edit the following lines to correspond to your compilers ####

F77        = ifort
F90        = ifort
CC         = icx
FCCPP      = cpp

FLINKER     = $(F90)

FCCPPFLAGS = -ansi -DLIBXC_VERSION=221  #Use this for versions <301
#FCCPPFLAGS = -ansi -DLIBXC_VERSION=400  #Use this for versions >400

##### Edit the following optimization flags for your system ####

# FFLAGS     = -O3 -funroll-loops -check all -warn interfaces
FFLAGS     = -O3 -funroll-loops -traceback -g
CFLAGS     = -O3 -I/home/elinscott/code/libxc-2.2.1/src
FCCPPFLAGS     = -O3 -I/home/elinscott/code/libxc-2.2.1/src

##### Edit the following LAPACK and BLAS library paths for your system ####

LIBS = -lmkl_intel_lp64  -lmkl_sequential -lmkl_core

##### Edit the following for to use libxc if available #####

OBJS_LIBXC =	exc_libxc_stub.o

# oncvpsp is compatible with libxc
# To build oncvpsp with libxc, uncomment the e following lines and edit
# the paths to point to your libxc library and include directories
# make clean in src before rebuilding after changing this

LIBS += -L/usr/local/lib -L/home/elinscott/code/libxc-2.2.1/src/.libs/ -L/home/elinscott/code/libxc-2.2.1/src/ -lxcf90 -lxc
FFLAGS += -I/usr/local/include/ -I/home/elinscott/code/libxc-2.2.1/src/

#LIBS += -L/home/drh/abinit/fallbacks/exports/lib -lxcf90 -lxc
#FFLAGS += -I/home/drh/abinit/fallbacks/exports/include

OBJS_LIBXC =	functionals.o exc_libxc.o 

