# Makefile for compiling with personal F90 and
# link with MPICH and FFTW libraries
# To compile and link write simply "make"
#
SHELL = /bin/bash
FFLAG = -O3 -w -r8
FCOMP = mpif90 -c ${FFLAG}
LINK = mpif90 
LIBS = -lm

OBJ = var_inc.o main.o para.o collision.o streaming.o initial.o partlib.o saveload.o

.SUFFIXES: .o .f90

.f90.o:
	${FCOMP} $*.f90 ${LIBS}
main: ${OBJ}
	${LINK} ${FFLAG} ${OBJ} ${LIBS} -o main

clean:
	rm -rf *.o *.mod main core*

