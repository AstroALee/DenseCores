GCC = g++
CPPFLAGS = -O2
INCLUDES = -I ./eigen/
all: DenseCore

DenseCore: DenseCoreMain.o SolveAll.o MagCylinder.o 
	$(GCC) $(CPPFLAGS) $(INCLUDES) -o DenseCoreCode DenseCoreMain.o SolveAll.o MagCylinder.o

DenseCoreMain.o : DenseCoreMain.cpp DenseCoreMain.H ErrorMessages.H
	$(GCC) $(CPPFLAGS) $(INCLUDES) -c DenseCoreMain.cpp

SolveAll.o : SolveAll.cpp SolveAll.H DenseCoreGlobals.H ErrorMessages.H
	$(GCC) $(CPPFLAGS) $(INCLUDES) -c SolveAll.cpp

MagCylinder.o : MagCylinder.cpp
	$(GCC) $(CPPFLAGS) $(INCLUDES) -c MagCylinder.cpp

clean:
	rm *.o DenseCoreCode
