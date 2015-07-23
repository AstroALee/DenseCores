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

// Includes boundary (n=0; all z from 0 up to and including Z & z=Z; all n from 0 up to and including N)
void AllocateState(double ***&curState)
{
    int i,j,k;
    curState = new double **[NStates];
    
    
    for(i=0; i< NStates; i++)
    {
        curState[i] = new double *[N+1];
        for(j=0; j< N+1; j++) curState[i][j] = new double[Z+1];
    }
    
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
    pLength = atof(argv[3]);
    RhoCenter = atof(argv[4]);
    MassGuess = atof(argv[5]);
    betaedge  = atof(argv[6]);
    
    // Print out information
    TELLME_Inputs();
    
    
    // Determines grid cell size for a given ratio of Z to R
    double ratio = 1.0;
    DeltaR = pLength/((double)(N+1.0));
    DeltaZ = ratio*pLength/((double)(Z+1.0));
    cout << DeltaR << " " << DeltaZ << endl;
    
    
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

