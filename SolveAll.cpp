#include "SolveAll.H"



void SolveAll(double ***curState)
{
    
    // First step is to solve Poisson's equation
    solvePoisson(curState);
    
    // Then you solve the induction equation
    solveMagnetic(curState);
    
    // With the updated state, determine if we need to keep relaxing
    
    
    // Once relaxed, calculate the actual mass in the box
    
    cout << "asdfasdf " << NStates << endl;
    
    
}




