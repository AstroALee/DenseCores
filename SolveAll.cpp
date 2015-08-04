#include "SolveAll.H"



void SolveAll(double ***curState)
{
    
    // First step is to solve Poisson's equation
    solvePoisson(curState);
    
    // Then you solve the induction equation
    //solveMagnetic(curState);
    
}




