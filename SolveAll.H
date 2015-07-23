#include "DenseCoreGlobals.H"
#include "ErrorMessages.H"






void solveMagnetic(double ***curState)
{
    
}


// Creates the matrix that is used
void createPoissonMatrix(double ***curState, double **PoMatrix, double *Source)
{
    int i,j,k;
    
    // First clear the matrix and source
    for(i=0;i<N+1;i++) for(j=0;j<Z+1;j++) PoMatrix[i][j] = 0.0;
    for(i=0;i<N*Z;i++) Source[i] = 0.0;
    
    // For each element, fill in the appropriate matrix elements based on the finite difference formulas
    
    
}

void solvePoisson(double ***curState)
{
    int i;
    
    /* We need to create the finite difference matrix and the source vector */
    
    // Allocate first
    double **PoMatrix;
    PoMatrix = new double *[N];
    for(i=0;i<N;i++) PoMatrix[i] = new double[Z];
    double Source[N*Z];
    double Soln[N*Z];
    
    // Then fill
    createPoissonMatrix(curState,PoMatrix,Source);
    
    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    
    
    
    // All done
    for(i=0;i<N;i++) delete[] PoMatrix[i];
    delete[] PoMatrix;
    
}


