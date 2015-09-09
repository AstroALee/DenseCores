#include "SolveAll.H"


void updateQ(double ***curState)
{
    // We use the top row (RhoTop) and the updated Q
    // to make a means of mapping a value of phi to a value of Q.
    // Beyond the value of phi at the corner where the top boundary meets
    // the counter, the value of Q is fixed, since all flux tubes then
    // touch the outer boundary, where rho = 1 and V = constant.

    // How many cells are before the contour on the upper boundary?
    int i,j;
    int Ncells = 0;



    for(j=0;j<Z;j++) for(i=0;i<N;i++)
    {
        // current value of Phi
        double curPhi = cPos(i,DeltaR) * curState[Apot][i][j];

        // Find two indices where curPhi is between


    }


}


void SolveAll(double ***curState)
{

    // First step is to solve Poisson's equation
    //solvePoisson(curState);

    // Then you solve the induction equation
    //solveMagnetic(curState);

    // Then update Q based on the new state
    //updateQ(curState);
}
