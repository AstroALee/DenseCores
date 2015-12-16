/* Dense Core */
#include "DenseCoreMain.H"



int main(int argc, char *argv[])
{
    // Start the clock
    clock_t startTime = clock();

    // Clean line
    cout << endl << endl;

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

    // Allocate Contour and Top Arrays
    VContour = new double[Z];
    RhoTop = new double[N];

#if DEBUG
    ContourValues = new double[Z];
    ContourATop = new double[N];
    ContourASide = new double[Z];
#endif

    // Allocate State Data
    AllocateState(curState);
    AllocateState(prevState);

    // Fill with initial guesses
    FillState(curState);
    CopyState(curState,prevState);


    // Enter the Grand Solving Loop, given the curState,
    // solve for an updated V and A, then updates Q
    SolveAll(curState);

    // Relax the solution
    //RelaxSolution(curState,prevState);

    // Determine if we have converged
    //ConvergeTest(curState,prevState);


    // Solution obtained, print
    //PrintScreen(curState);
    PrintFile(curState);


    // Publish!
    delete[] VContour;
    delete[] RhoTop;

#if DEBUG
    delete[] ContourValues;
    delete[] ContourATop;
    delete[] ContourASide;
#endif

    // Stop the clock
    clock_t endTime = clock();
    cout << "Execution took " << double( endTime - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << endl;

    return 0;
}

// Fill the state with the initial guess
void FillState(double ***curState)
{
    // Fills Vpot
    PointMassPotential(curState);

    // Initial beta, used with initializing A
    // Then we don't need it again.
    double initbeta[N];

    // Fills Q
    InitialQ(curState,initbeta);

    // Fills Apot
    InitialA(curState,initbeta);

}

// Includes boundary (n=0; all z from 0 up to and including Z & z=Z; all n from 0 up to and including N)
void AllocateState(double ***&State)
{
    int i,j,k;
    State = new double **[NStates];


    for(i=0; i< NStates; i++)
    {
        State[i] = new double *[N];
        for(j=0; j< N; j++) State[i][j] = new double[Z];
    }

}

// Copies first state array to the second state array
void CopyState(double ***State1, double ***State2)
{
    int i,j,k;
    for(i=0;i<NStates;i++) for(j=0;j<N;j++) for(k=0;k<Z;k++) State2[i][j][k] = State1[i][j][k];
}


void ConvergeTest(double ***curState, double ***prevState)
{
    // Sees how different the current state is from the previous state
    double Vtol = 1e-2;
    double Atol = 1e-2;
    double Qtol = 1e-2;




}

void RelaxSolution(double ***curState, double ***prevState)
{
    // Uses a weighting of the previous and current solution as the actual solution
    int i,j,k;

    double prevWeight = 0.5;

    for(i=0;i<NStates;i++) for(j=0;j<N;j++) for(k=0;k<Z;k++)
        curState[i][j][k] = (1.0-prevWeight)*curState[i][j][k] + prevWeight*prevState[i][j][k];



}
