// Eigen Routines (path in Makefile, needs to come first)
#include <Eigen/Dense>


#include "DenseCoreGlobals.H"
#include "ErrorMessages.H"

using namespace Eigen;



// Row major on interior cells : i = Ridx + Zidx * N-2
int Ridx(int i)
{
    return i%(N-2);
}

int Zidx(int i)
{
    return (i-Ridx(i))/(N-2);
}



// Creates the matrix that is used
void createAmpereMatrix(double ***curState, MatrixXd AmMatrix, VectorXd Source)
{
    int i,j,k;
    
    // For each element, fill in the appropriate matrix elements based on the finite difference formulas
    
    //   dQ      (dQ/dx)dx + (dQ/dy)dy
    //   --   =  ---------------------
    //  dPhi     (dP/dx)dx + (dP/dy)dx
    
}


// Creates the matrix that is used
void createPoissonMatrix(double ***curState, MatrixXd PoMatrix, VectorXd Source)
{
    int i,j,k;
    
    
    // Source Vector = 4 pi q exp(-V_previous)
    for(i=0;i<(Z-2)*(N-2);i++)
        Source(i) = pow(DeltaR,2)*4.0*PI*curState[Q][Ridx(i)+1][Zidx(i)+1]*exp(-curState[Vpot][Ridx(i)+1][Zidx(i)+1]);
    

    // Loops over every cell to be solved, each which has it's own finite difference
    // scheme and fills one particular row of the big matrix A in   A u = f
    for(i=0;i<(Z-2)*(N-2);i++)
        {
            // Each of these cells gets a finite difference scheme
            
            double wi = DeltaR/2.0/cPos(Ridx(i)+1,DeltaR);
            double rho = pow(DeltaR/DeltaZ,2);
            
            // If there's no boundary exceptions, here are the terms
            double rightV = 1.0 + wi;
            double leftV  = 1.0 - wi;
            double upV    = rho;
            double downV  = -rho;
            double midV   = -2.0*(1.0+rho);
            
            // Check to see if we're at the bottom
            if( Zidx(i)==0 )
            {
                // Then the V_(k-1) term gets added to the V_i term instead
                midV = midV + downV;
                downV = 0.0; // <-- Not really needed, there's no matrix entry
            }
            
            // Check to see if we're at the top
            if( Zidx(i) == (Z-2)-1 )
            {
                // Then the V_(k+1) term gets added to the V_i term instead
                midV = midV + upV;
                upV = 0.0; // <-- Not really needed, there's no matrix entry
            }
            
            // Check to see if we're at the left
            if( Ridx(i) == 0 )
            {
                // Then the V_(i-1) term gets added to the V_i term instead
                midV = midV + leftV;
                leftV = 0.0; // <-- Not really needed, there's no matrix entry
            }
            
            // Check to see if we're near the right boundary
            // If we're beyond it, who cares, overwrite source vector
            if( cPos(Ridx(i)+1,DeltaR) > VContour[Zidx(i)+1])
            {
                // One could check to see if this particular i value is never
                // access in any other row...
                rightV = 0;
                leftV = 0;
                upV = 0;
                downV = 0;
                midV = 1;
                Source(i) = VContour[Zidx(i)+1];
                
            }
            // Else if we're just shy of the boundary
            else if( cPos(Ridx(i)+1+1,DeltaR) >= VContour[Zidx(i)+1] )
            {
                // Then we need to linearly interpolate V based on the current
                // value of V_ik and the contour value.
                // (We only interpolate in the R direction at the moment,
                //  could also interpolate up and down as well)
                double leftpointV = curState[Vpot][Ridx(i)+1-1][Zidx(i)+1-1];
                double leftpos = cPos(Ridx(i)+1-1,DeltaR);
                double m = (RightContour - leftpointV)/(VContour[Zidx(i)+1] - leftpos);
                
                double midV = leftpointV + m*(cPos(Ridx(i)+1,DeltaR) - leftpos);
                
                // Now we want this midV to become the new value. This matrix
                // row becomes relatively simple. Overwrites the source vector
                rightV = 0;
                leftV = 0;
                upV = 0;
                downV = 0;
                midV = 1;
                Source(i) = midV;
            }
            
            // Now we fill the entries of the matrix with the values
            PoMatrix(i,i) = midV;
            if(Ridx(i)!=0) PoMatrix(i,i-1) = leftV;
            if(Ridx(i)!=N-2-1) PoMatrix(i,i+1) = rightV;
            if(Zidx(i)!=Z-2-1) PoMatrix(i,i+(N-2)) = upV;
            if(Zidx(i)!=0) PoMatrix(i,i-(N-2)) = downV;
            
        }
    
    
}

void solvePoisson(double ***curState)
{
    
    /* We need to create the finite difference matrix and the source vector */
    
    // Allocate first using the Eigen structures, initalize to zero
    // The boundaries are not solved for
    MatrixXd PoMatrix = MatrixXd::Zero((Z-2)*(N-2),(Z-2)*(N-2));
    VectorXd Source = VectorXd::Zero((Z-2)*(N-2));
    VectorXd Soln = VectorXd::Zero((Z-2)*(N-2));
    
    
    // Then fill
    createPoissonMatrix(curState,PoMatrix,Source);
    
    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    Soln = PoMatrix.colPivHouseholderQr().solve(Source);
    
    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (PoMatrix*Soln - Source).norm() / Source.norm();
    
    
}


void solveMagnetic(double ***curState)
{
    
    /* We need to create the finite difference matrix and the source vector */
    
    // Allocate first using the Eigen structures
    MatrixXd AmMatrix = MatrixXd::Zero(N+1,Z+1);
    VectorXd Source = VectorXd::Zero(N*Z);
    VectorXd Soln = VectorXd::Zero(N*Z);
    
    
    // Then fill
    createAmpereMatrix(curState,AmMatrix,Source);
    
    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    Soln = AmMatrix.colPivHouseholderQr().solve(Source);
    
    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (AmMatrix*Soln - Source).norm() / Source.norm();
}
