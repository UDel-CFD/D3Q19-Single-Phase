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
#Defines file siffixes we'll use
.SUFFIXES: .o .f90

.f90.o:
	${FCOMP} $*.f90 ${LIBS}
main: ${OBJ}
	${LINK} ${FFLAG} ${OBJ} ${LIBS} -o main
#Delete all .o and other object files
	rm -rf ${OBJ} var_inc.mod
#Note this forces a recompile everytime

#Type make clean to clean all object and program files
clean:
	rm -rf *.o *.mod main core*

