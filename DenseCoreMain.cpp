/* Dense Core */
#include "DenseCoreMain.H"


void SolveAll(double ***curState);


// Fill the state with the initial guess
void FillState(double ***curState)
{
    int i,j,k;
    
    // Fills Vpot
    PointMassPotential(curState);
    

    // Fills Q
    double initbeta[N+1];
    InitialQ(curState,initbeta);
    
    
    // Fills Apot
    InitialA(curState,initbeta);
 
    
}


int main(int argc, char *argv[])
{
    // Usual counters
    int i,j,k;
    
    // Checks and balances on inputs
    WATERLOO_ArgCountCheck(argc);
    
    // Read in from command line
    N = atoi(argv[1]);
    Z = atoi(argv[2]);
    RhoCorner = atof(argv[3]);
    RhoCenter = atof(argv[4]);
    MassGuess = atof(argv[5]);
    betaedge  = atof(argv[6]);
    
    
    // Allocate Contour Array
    VContour = new double[Z+1];
    
    // Allocate State Data, Fill with initial guesses
    AllocateState(curState);
    FillState(curState);
    
    
    // Enter the Grand Solving Loop
    //SolveAll(curState);
        
    
    // Solution obtained, print
    //PrintScreen(curState);
    PrintFile(curState);
    
    // Publish!
    
    delete[] VContour;
    return 0;
}

