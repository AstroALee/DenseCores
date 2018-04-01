


# Variables you can change
# Uniform grid cell size (same in r as in z)
YESUNIFORM = 0
# Ratio of filament radius to half height 
CYLINDERRADRAT = 1.0
# Max number of loops in Poisson/Ampere solve
LOOPMAX = 100
# Are we in debug mode? (1= IC cylinder ; 2= Messy IC, relaxes to cylinder )
DEBUG = 0
#
# Append to flags
CPPFLAGS = -DLOOPMAX=$(LOOPMAX) -DDEBUG=$(DEBUG) -DUNIFORM=$(YESUNIFORM) -DRADIALRATIO=$(CYLINDERRADRAT)

# Code compiling related stuff
GCC = g++ 
CPPFLAGS += -O2 #-Wall
INCLUDES = -I ./eigen/

CPPFLAGS += $(INCLUDES)

OUTPUTNAME = FilamentCode

all: Filament_Main

Filament_Main : Filament_Main.o Filament_PrepInit.o Filament_newSolveAll.o Filament_newMagCyl.o Filament_Poisson.o Filament_Ampere.o
	$(GCC) $(CPPFLAGS) -o $(OUTPUTNAME) $^
	rm ./*.o

clean:
	rm $(OUTPUTNAME) *.o
