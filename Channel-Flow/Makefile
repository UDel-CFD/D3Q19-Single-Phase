#=============================================================================
#University of Delaware LBM D3Q19 Single Phase Simulation 
#Copyright (C) 2017 Lian-Ping Wang

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#=============================================================================

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
OBJ = var_inc.o main.o para.o collision.o initial.o partlib.o saveload.o
#Defines file siffixes we'll use
.SUFFIXES: .o .f90

.f90.o:
	${FCOMP} $*.f90 ${LIBS}
main: ${OBJ}
	${LINK} ${FFLAG} ${OBJ} ${LIBS} -o main
#Optional auto clean object files
	rm -rf *.o *.mod
# Debug build changes flags to assist in program debugging
debug: FFLAG = -g -r8 -C -CB -check uninit -traceback -debug all -fp-model precise -fpe0 -fp-stack-check
debug: main
#Remove Old Object files, useful when recompiling
clean:
	rm -rf *.o *.mod core*
#Type make clean to clean all object and program files
fullclean:
	rm -rf *.o *.mod main core*

