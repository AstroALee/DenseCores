#include "SolveAll.H"


void updateQ(double ***curState)
{
    // We use the top row (RhoTop) and the updated Q
    // to make a means of mapping a value of phi to a value of Q.
    // Beyond the value of phi at the corner where the top boundary meets
    // the counter, the value of Q is fixed, since all flux tubes then
    // touch the outer boundary, where rho = 1 and V = constant.

    // How many cells are before the contour on the upper boundary?
    int i,j,k;
    int Ncells = 0;
    for(i=0;i<N;i++) if( cPos(i,DeltaR) > VContour[i] ) { Ncells = i; break; }

    // New values of Q on top row
    double NewQ[N];
    for(i=0;i<N;i++) NewQ[i] = RhoTop[i]*exp( curState[Vpot][i][Z-1] );

    // Top rows of Phi
    double TopPhi[N];
    for(i=0;i<N;i++) TopPhi[i] = cPos(i,DeltaR)*curState[Apot][i][Z-1];

    // Value of Phi at the contour
    double curB = sqrt(8.0*PI*1.0/betaedge);
    double PhiEnd = VContour[Z-1]*0.5*curB*pow(VContour[Z-1],1);

    // Loops over each cell, grabs its current value of Phi, and uses the top row to
    // interpolate a new value of Q. Assumes that Phi is monotonic.
    for(j=0;j<Z;j++) for(i=0;i<N;i++)
    {
        // current value of Phi
        double curPhi = cPos(i,DeltaR) * curState[Apot][i][j];

        // Find two indices where curPhi is between
        int idx = 0, done = 0;
        double rPhi = 0, lPhi = 0;
        for(k=0;k<Ncells;k++)
        {
            idx = k;

            if( curPhi <= TopPhi[k] )
            {
                done = 1;
                rPhi = TopPhi[k];
                break;
            }
            else
            {
                lPhi = rPhi;
            }
        }

        // If done = 1, we found index. If not, then this particular value of Phi is in the
        // force free zone.
        if(done)
        {
            // Interpolate
            double m = ( NewQ[idx] - NewQ[idx-1] ) / ( TopPhi[idx] - TopPhi[idx-1] ); // what if idx = 0 ???
            curState[Q][i][j] = NewQ[idx-1] + m*(curPhi - TopPhi[idx-1]);
        }
        else
        {
            curState[Q][i][j] = ContourQ;
        }
    }

// At this point, Q has been updated using the updated values of V and A.
}


void SolveAll(double ***curState)
{

    // First step is to solve Poisson's equation
    //solvePoisson(curState);

    // Then you solve the induction equation
    solveMagnetic(curState);

    // Then update Q based on the new state
    //updateQ(curState);

}
