// Eigen Routines (path in Makefile, needs to come first)
#include <Eigen/Dense>


#include "DenseCoreGlobals.H"
#include "ErrorMessages.H"

using namespace Eigen;



// Row major on interior cells : i = Ridx + Zidx * N+1-2
// (N+1 total cells, ignore 2 edges)
int Ridx(int i)
{
    return i%(N+1-2);
}

int Zidx(int i)
{
    return (i-Ridx(i))/(N+1-2);
}



// Creates the matrix that is used
void createAmpereMatrix(double ***curState, MatrixXd AmMatrix, VectorXd Source)
{
    int i,j,k;
    
    // For each element, fill in the appropriate matrix elements based on the finite difference formulas
    
    //   dQ      (dQ/dx)dx + (dQ/dy)dy
    //   --   =  ---------------------
    //  dPhi     (dP/dx)dx + (dP/dy)dx
    
    // Source Vector = 4 pi q exp(-V_previous)
    double qdot[(Z+1-2)*(N+1-2)];
    
    for(i=0;i<(Z+1-2)*(N+1-2);i++)
        Source(i) = -4*PI*pow(DeltaR,2)*cPos(Ridx(i),DeltaR)*exp(-curState[Vpot][Ridx(i)][Zidx(i)])*qdot[i];
    
    
    
}


// Creates the matrix that is used
void createPoissonMatrix(double ***curState, MatrixXd& PoMatrix, VectorXd& Source)
{
    int i,j,k;
    
    
    // Source Vector = 4 pi q exp(-V_previous)
    for(i=0;i<(Z+1-2)*(N+1-2);i++)
        Source(i) = pow(DeltaR,2)*4.0*PI*curState[Q][Ridx(i)+1][Zidx(i)+1]*exp(-curState[Vpot][Ridx(i)+1][Zidx(i)+1]);
    

    // Loops over every cell to be solved, each which has it's own finite difference
    // scheme and fills one particular row of the big matrix A in   A u = f
    for(i=0;i<(Z+1-2)*(N+1-2);i++)
        {
            // Each of these cells gets a finite difference scheme
            
            double wi = DeltaR/2.0/cPos(Ridx(i)+1,DeltaR);
            double rho = pow(DeltaR/DeltaZ,2);
            
            // If there's no boundary exceptions, here are the terms
            double rightV = 1.0 + wi;
            double leftV  = 1.0 - wi;
            double upV    = rho;
            double downV  = rho;
            double midV   = -2.0*(1.0+rho);
            
            // Check to see if we're at the bottom
            if( Zidx(i)==0 )
            {
                // Then the V_(k-1) term gets added to the V_i term instead
                midV = midV + downV;
                downV = 0.0; // <-- Not really needed, there's no matrix entry
            }
            
            // Check to see if we're at the top
            if( Zidx(i) == (Z+1-2)-1 )
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
            
            //cout << "i=" << i << " " << cPos(Ridx(i)+1+1,DeltaR) << " " << VContour[Zidx(i)+1] << endl;
            
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
                double m = (ContourVpot - leftpointV)/(VContour[Zidx(i)+1] - leftpos);
                
                double midIV = leftpointV + m*(cPos(Ridx(i)+1,DeltaR) - leftpos);
                
                // Now we want this midV to become the new value. This matrix
                // row becomes relatively simple. Overwrites the source vector
                rightV = 0;
                leftV  = 0;
                upV    = 0;
                downV  = 0;
                midV   = 1;
                Source(i) = midIV;
            }
            
            // Now we fill the entries of the matrix with the values
            PoMatrix(i,i) = midV;
            if(Ridx(i)!=0) PoMatrix(i,i-1) = leftV;
            if(Ridx(i)!=N+1-2-1) PoMatrix(i,i+1) = rightV;
            if(Zidx(i)!=Z+1-2-1) PoMatrix(i,i+(N+1-2)) = upV;
            if(Zidx(i)!=0) PoMatrix(i,i-(N+1-2)) = downV;
            
        }
    
    
}

void solvePoisson(double ***curState)
{
    double L2tol = 1e-10;
    
    /* We need to create the finite difference matrix and the source vector */
    
    // Allocate first using the Eigen structures, initalize to zero
    // The boundaries are not solved for
    MatrixXd PoMatrix = MatrixXd::Zero((Z+1-2)*(N+1-2),(Z+1-2)*(N+1-2));
    VectorXd Source = VectorXd::Zero((Z+1-2)*(N+1-2));
    VectorXd Soln = VectorXd::Zero((Z+1-2)*(N+1-2));
    
    
    // Then fill
    createPoissonMatrix(curState,PoMatrix,Source);
    
    //cout << PoMatrix << endl;
    //cout << Source << endl;
    
    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    Soln = PoMatrix.colPivHouseholderQr().solve(Source);
    
    //cout << endl << Soln << endl;
    
    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (PoMatrix*Soln - Source).norm() / Source.norm();
    cout << "Gravity Matrix L2 error = " << L2error << endl;
    
    if(L2error < L2tol)
    {
        // Copy new results to curState
        int i,j;
        
        //Bottom-left corner
        curState[Vpot][0][0] = Soln[0];
        //Bottom row
        for(i=0;i<N+1-2;i++) curState[Vpot][i+1][0] = Soln[i];
        //Bottom-right corner (by construction should be ContourVpot)
        curState[Vpot][N][0] = ContourVpot;
        
        //Left side
        for(i=0;i<Z+1-2;i++) curState[Vpot][0][i+1] = Soln[ i*(N+1-2) ];
        
        //Top-left corner (should be able to look either right or below and get
        //                 the same answer... could check this?)
        curState[Vpot][0][Z] = curState[Vpot][0][Z-1];
        //Top row
        for(i=0;i<N+1-2;i++) curState[Vpot][i+1][Z] = Soln[ i + (Z+1-2-1)*(Z+1-2)];
        //Top-right corner (dV/dz = 0)
        curState[Vpot][N][Z] = curState[Vpot][N][Z-1];
        
        //Interior
        int radjust = 1;
        int zadjust = 1;
        for(j=0;j<(Z+1-2);j++)
            for(i=0;i<(N+1-2);i++)
                curState[Vpot][i+radjust][j+zadjust] = Soln[ i + j*(N+1-2) ];
        
        // Overall, who cares about what's right of the contour
        for(i=0;i<N+1;i++) for(j=0;j<Z+1;j++)
            if(cPos(i,DeltaR) >= VContour[j]) curState[Vpot][i][j]=ContourVpot;
            
            
    }
    
}


void solveMagnetic(double ***curState)
{
    
    /* We need to create the finite difference matrix and the source vector */
    
    // Allocate first using the Eigen structures
    MatrixXd AmMatrix = MatrixXd::Zero((Z+1-2)*(N+1-2),(Z+1-2)*(N+1-2));
    VectorXd Source = VectorXd::Zero((Z+1-2)*(N+1-2));
    VectorXd Soln = VectorXd::Zero((Z+1-2)*(N+1-2));
    
    
    // Then fill
    createAmpereMatrix(curState,AmMatrix,Source);
    
    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    Soln = AmMatrix.colPivHouseholderQr().solve(Source);
    
    cout << endl << Soln << endl;
    
    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (AmMatrix*Soln - Source).norm() / Source.norm();
    cout << "Magnetic Matrix L2 error = " << L2error << endl;
    
    if(L2error < 1)
    {
        // Copy new results to curState
        int i,j;
        
        //Bottom-left corner
        curState[Apot][0][0] = Soln[0];
        //Bottom row
        for(i=0;i<N+1-2;i++) curState[Apot][i+1][0] = Soln[i];
        //Bottom-right corner
        curState[Apot][N][0] = 0.0;
        
        //Left side
        for(i=0;i<Z+1-2;i++) curState[Apot][0][i+1] = Soln[ i*(N+1-2) ];
        
        //Top-left corner
        curState[Vpot][0][Z] = curState[Apot][0][Z-1];
        //Top row
        for(i=0;i<N+1-2;i++) curState[Apot][i+1][Z] = Soln[ i + (Z+1-2-1)*(Z+1-2)];
        //Top-right corner
        curState[Apot][N][Z] = curState[Apot][N-1][Z];
        
        //Interior
        int radjust = 1;
        int zadjust = 1;
        for(j=0;j<(Z+1-2);j++)
            for(i=0;i<(N+1-2);i++)
                curState[Apot][i+radjust][j+zadjust] = Soln[ i + j*(N+1-2) ];
        
        for(i=0;i<N+1;i++) for(j=0;j<Z+1;j++)
            if(cPos(i,DeltaR) > VContour[j]) curState[Apot][i][j]= 0.0;
        
        
    }
}
