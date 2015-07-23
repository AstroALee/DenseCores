GCC = g++
CPPFLAGS = -O2

all: DenseCore

DenseCore: DenseCoreMain.o SolveAll.o MagCylinder.o 
	$(GCC) $(CPPFLAGS) -o DenseCoreCode DenseCoreMain.o SolveAll.o MagCylinder.o

DenseCoreMain.o : DenseCoreMain.cpp DenseCoreMain.H ErrorMessages.H
	$(GCC) $(CPPFLAGS) -c DenseCoreMain.cpp

SolveAll.o : SolveAll.cpp SolveAll.H DenseCoreGlobals.H ErrorMessages.H
	$(GCC) $(CPPFLAGS) -c SolveAll.cpp

MagCylinder.o : MagCylinder.cpp
	$(GCC) $(CPPFLAGS) -c MagCylinder.cpp

clean:
	rm *.o DenseCoreCode
