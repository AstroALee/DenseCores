/* Dense Core */
#include "DenseCoreMain.H"


void SolveAll(double ***curState);


// Fill the state with the initial guess
void FillState(double ***curState)
{
    // Initial beta, used with initializing A
    // Then we don't need it again.
    double initbeta[N];
    

    // Fills Vpot
    PointMassPotential(curState);
    
    // Fills Q
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
        curState[i] = new double *[N];
        for(j=0; j< N; j++) curState[i][j] = new double[Z];
    }
    
}


int main(int argc, char *argv[])
{
    // Start the clock
    clock_t startTime = clock();
    
    
    
    cout << endl << endl;
    
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
    SHOWME_Inputs(argv);
    
    
    // Determines grid cell size for a given ratio of Z to R
    double ratio = 1.0;
    DeltaR = pLength/((double)(N));
    DeltaZ = ratio*pLength/((double)(Z));
    //cout << DeltaR << " " << DeltaZ << endl;
    
    // Allocate Contour Array
    VContour = new double[Z];
    RhoTop = new double[N];
    
#if DEBUG
    ContourValues = new double[Z];
#endif
    
    // Allocate State Data, Fill with initial guesses
    AllocateState(curState);
    FillState(curState);
    
    
    // Enter the Grand Solving Loop, given the curState,
    // solve for an updated V and A
    SolveAll(curState);
        
    
    // Solution obtained, print
    //PrintScreen(curState);
    PrintFile(curState);
    
    // Publish!
    delete[] VContour;
    delete[] RhoTop;
    
#if DEBUG
    delete[] ContourValues;
#endif
    
    // Stop the clock
    clock_t endTime = clock();
    cout << "Execution took " << double( endTime - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << endl;
    
    return 0;
}

