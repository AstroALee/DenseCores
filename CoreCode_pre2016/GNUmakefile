# Makefile for the Dense Core code
# Written by Aaron Lee, 2015
# $< :: matches the first dependency
# $^ :: matches all the targets, removes duplicates
# $@ :: matches the first target

GCC = g++
CPPFLAGS = -O2 #-Wall
INCLUDES = -I ./eigen/

OUTPUTNAME = DenseCoreCode
OFILEFOLDER = OutFiles

all: DenseCore

DenseCore: DenseCoreMain.o SolveAll.o MagCylinder.o
	$(GCC) $(CPPFLAGS) $(INCLUDES) -o $(OUTPUTNAME) $^ 
	mv *.o $(OFILEFOLDER)/

DenseCoreMain.o : DenseCoreMain.cpp
	$(GCC) $(CPPFLAGS) $(INCLUDES) -c $<

SolveAll.o : SolveAll.cpp
	$(GCC) $(CPPFLAGS) $(INCLUDES) -c $<

MagCylinder.o : MagCylinder.cpp
	$(GCC) $(CPPFLAGS) $(INCLUDES) -c $<

clean:
	rm $(OUTPUTNAME)
	rm $(OFILEFOLDER)/*.o
