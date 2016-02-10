# Makefile for the Dense Core code
# Written by Aaron Lee, 2015
# $< :: matches the first dependency
# $^ :: matches all the targets, removes duplicates
# $@ :: matches the first target


# Variables you can change
# Uniform grid cell size (same in r as in z)
YESUNIFORM = 0
CYLINDERRADRAT = 1.0
DOONE = 1
DEBUG = 1
CPPFLAGS = -DDOONE=$(DOONE) -DDEBUG=$(DEBUG) -DUNIFORM=$(YESUNIFORM) -DRADIALRATIO=$(CYLINDERRADRAT)

# Code compiling related stuff
GCC = g++ 
CPPFLAGS += -O2 #-Wall
INCLUDES = -I ./eigen/

CPPFLAGS += $(INCLUDES)

OUTPUTNAME = FilamentCode

all: Filament_Main

Filament_Main : Filament_Main.o Filament_PrepInit.o Filament_SolveAll.o Filament_MagCyl.o Filament_Poisson.o Filament_Ampere.o
	$(GCC) $(CPPFLAGS) -o $(OUTPUTNAME) $^
	rm ./*.o

clean:
	rm $(OUTPUTNAME)
