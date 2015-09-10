# Makefile for compiling with personal F90 and
# link with MPI libraries
# To compile and link write simply "make"
#
SHELL = /bin/bash
FFLAG = -O3 -w -r8
FCOMP = mpif90 -c ${FFLAG}
LINK = mpif90 
LIBS = -lm
#List of files needed to compile
OBJ = var_inc.o main.o para.o collision.o initial.o partlib.o saveload.o benchmark.o
#Defines file siffixes we'll use
.SUFFIXES: .o .f90

.f90.o:
	${FCOMP} $*.f90 ${LIBS}
main: ${OBJ}
	${LINK} ${FFLAG} ${OBJ} ${LIBS} -o main
#Optional auto clean object files
	rm -rf *.o *.mod
# Debug build changes flags to assist in program debugging
debug: FFLAG = -g -r8 -C -CB -traceback -debug all -fp-model precise
debug: main
#Remove Old Object files, useful when recompiling
clean:
	rm -rf *.o *.mod core*
#Type make clean to clean all object and program files
fullclean:
	rm -rf *.o *.mod main core*

