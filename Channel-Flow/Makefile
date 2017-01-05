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
# Updated by Nicholas Geneva June 01, 2016
# make: normal compiling
# make debug: compile with various debugging flags (no optimization)
# make clean: remove object and module files along with error logs
# make fullclean: clean + remove main executable

SHELL = /bin/bash
FFLAG = -O3 -w -r8
FCOMP = mpif90 -c ${FFLAG}
LINK = mpif90 
LIBS = -lm

#Bin Folder Location (or where the object files get placed)
OBJDIR =./bin
#List of files needed to compile
OBJNAMES = var_inc.o main.o para.o collision.o initial.o saveload.o
MODS = var_inc.mod #Modules used in the program
OBJ = $(addprefix $(OBJDIR)/,${OBJNAMES}) #Add on bin folder to directory

#Defines file siffixes we'll use
.SUFFIXES: .o .f90

#Check if bin directory exists, create if needed
main ${OBJ}: | $(OBJDIR) #Make having the directory a prereq for the objest files
$(OBJDIR):
	mkdir -p ${OBJDIR}

#Compile module seperately to force entire program to recompile
${MODS} : ${MODS:.mod=.f90}
	${FCOMP} $< ${LIBS} -o $@
	
#Compile object files into bin
${OBJDIR}/%.o : %.f90 ${MODS}
	${FCOMP} $< ${LIBS} -o $@

#Create executable
main : ${OBJ}
	${LINK} ${FFLAG} ${OBJ} ${LIBS} -o main

# Debug build changes flags to assist in program debugging
.PHONY: debug
debug: FFLAG = -g -r8 -C -CB -check uninit -traceback -debug all -fp-model precise -fpe0 -fp-stack-check
debug: main
#Remove Old Object files, useful when recompiling
.PHONY: clean
clean:
	rm -rf bin/*.o *.mod *.o core*
#Type make clean to clean all object and program files
.PHONY: fullclean
fullclean:
	rm -rf *.o *.mod main core*