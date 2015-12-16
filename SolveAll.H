// Eigen Routines (path in Makefile, needs to come first)
#include <Eigen/Dense>


#include "DenseCoreGlobals.H"
#include "ErrorMessages.H"

using namespace Eigen;



// Row major on interior cells : i = Ridx + Zidx * N
int Ridx(int i)
{
    return i%(N);
}

int Zidx(int i)
{
    return (i-Ridx(i))/(N);
}


// Creates the matrix that is used
void createAmpereMatrix(double ***curState, MatrixXd& AmMatrix, VectorXd& Source)
{

    // For each element, fill in the appropriate matrix elements based on the finite difference formulas

    //   dQ      (dQ/dx)dx + (dQ/dy)dy
    //   --   =  ---------------------
    //  dPhi     (dP/dx)dx + (dP/dy)dx


    int i,j,k;

    double dqdp[Z*N];

    for(i=0;i<N*Z;i++)
    {
        // At the left edge dq/dx = 0 , at top and bottom dq/dy = 0
        double dqdx = 0, dqdy = 0, dpdx = 0, dpdy = 0;

        // Phi is r*A

        // if left edge or right edge (which will ultimately be always right of contour)
        if(Ridx(i)==0 || Ridx(i)==N-1) { dqdx = 0; dpdx = 0; }
        else
        {
            dqdx = ( curState[Q][Ridx(i)+1][Zidx(i)] - curState[Q][Ridx(i)-1][Zidx(i)] ) / (2*DeltaR);
            dpdx = ( cPos(Ridx(i)+1,DeltaR)*curState[Apot][Ridx(i)+1][Zidx(i)] - cPos(Ridx(i)-1,DeltaR)*curState[Apot][Ridx(i)-1][Zidx(i)] ) / (2*DeltaR);
        }

        // if bottom or top
        if(Zidx(i)==0 || Zidx(i)==Z-1) { dqdy = 0; dpdy = 0; }
        else
        {
            dqdy = ( curState[Q][Ridx(i)][Zidx(i)+1] - curState[Q][Ridx(i)][Zidx(i)-1] ) / (2*DeltaZ);
            dpdy = cPos(Ridx(i),DeltaR)*( curState[Apot][Ridx(i)][Zidx(i)+1] - curState[Apot][Ridx(i)][Zidx(i)-1] ) / (2*DeltaZ); // same R position, so factored out
        }


        // Right-side can provide some exceptions
        // if beyond right contour, who cares, set the values so the rigth-hand side = 0
        if( cPos(Ridx(i),DeltaR) >= VContour[Zidx(i)])
        {
            dqdx = 0;
            dpdy = 1;
            dqdx = 0;
            dpdy = 1;
        }
        // else if right next to contour
        else if( cPos(Ridx(i)+1,DeltaR) > VContour[Zidx(i)] )
        {
            // Need to know position to the contour.
            double leftP = cPos(Ridx(i)-1,DeltaR);
            double rightP= VContour[Zidx(i)];
            double leftQ = curState[Q][Ridx(i)-1][Zidx(i)];

            dqdx = (ContourQ - leftQ)/( rightP - leftP );


            // The contour value of Phi assumes uniform field (rho=1.0 by assumption)
            double curB = sqrt(8.0*PI*1.0/betaedge);
            double rightPhi = VContour[Zidx(i)]*( 0.5*curB*pow(VContour[Zidx(i)],1) ); // r = Contour
            double leftPhi = cPos(Ridx(i)-1,DeltaR)*curState[Apot][Ridx(i)-1][Zidx(i)];

            dpdx = (rightPhi - leftPhi)/( rightP - leftP );

        }

        // NOTE: Need something to look if contour is passed above and below too.


        // If we're at a corner, we'll interpolate the value of A rather than deal with dqdPhi
        // so who cares about the value of dQdPhi. If we're beyond the contour, who cares.

        if( Ridx(i)==0 && ( firstASolve==1 || (Zidx(i)==0 || Zidx(i)==Z-1)) ) dqdp[i] = 0; // <-- yea?
        else dqdp[i] = ( dqdx*DeltaR + dqdy*DeltaZ ) / ( dpdx*DeltaR + dpdy*DeltaZ ) ;

    }

    // Source Vector = 4 pi r * exp(-V) * dq/dPhi
    for(i=0;i<(Z)*(N);i++)
        Source(i) = -4.0*PI*pow(DeltaR,2)*cPos(Ridx(i),DeltaR)*exp(-curState[Vpot][Ridx(i)][Zidx(i)])*dqdp[i];

#if DEBUG
    for(i=0;i<Z*N;i++) Source(i) = 0.0; // initial conditions should be a force-free field.
#endif


    // Loops over every cell to be solved, each which has it's own finite difference
    // scheme and fills one particular row of the big matrix A in   A u = f
    for(i=0;i<(Z)*(N);i++)
    {
        // Each of these cells gets a finite difference scheme

        double wi = DeltaR/2.0/cPos(Ridx(i),DeltaR);
        double rho = pow(DeltaR/DeltaZ,2);

        // If there's no boundary exceptions, here are the terms
        double rightV = 1.0 + wi;
        double leftV  = 1.0 - wi;
        double upV    = rho;
        double downV  = rho;
        double midV   = -2.0*(1.0+rho+2*pow(wi,2));

        // Check to see if we're at the bottom
        if( Zidx(i)==0 )
        {
            // Then the V_(k-1) term gets added to the V_k+1 term instead
            //midV = midV + downV;
            upV = upV + downV;
            downV = 0.0; // <-- Not really needed, there's no matrix entry

            // If debugging, still using periodic Bcond
        }


        // Check to see if we're at the top
        if( Zidx(i) == Z-1 )
        {
            // Then the V_(k+1) term gets added to the V_k-1 term instead
            //midV = midV + upV;
            downV = downV + upV;
            upV = 0.0; // <-- Not really needed, there's no matrix entry

#if DEBUG
            // We override the top with the contour value
            Source(i) = ContourATop[Ridx(i)];
            downV = 0;
            upV = 0;
            leftV = 0;
            rightV = 0;
            midV = 1;
#endif
        }

        // Check to see if we're at the left
        if( Ridx(i) == 0 )
        {
            // Then the V_(i-1) term gets added to the V_i+1 term instead
            //midV = midV + leftV;
            rightV = rightV + leftV;
            leftV = 0.0; // <-- Not really needed, there's no matrix entry

            // If debugging, still using periodic Bcond
        }

        //cout << "i=" << i << " " << cPos(Ridx(i)+1+1,DeltaR) << " " << VContour[Zidx(i)+1] << endl;

        // Check to see if we're near the right boundary
        // If we're beyond it, who cares, overwrite source vector
        if( cPos(Ridx(i),DeltaR) >= VContour[Zidx(i)])
        {


            // One could check to see if this particular i value is never
            // access in any other row...
            rightV = 0;
            leftV = 0;
            upV = 0;
            downV = 0;
            midV = 1;

            double curB = sqrt(8.0*PI*1.0/betaedge);
            Source(i) = 0.5*curB*pow(VContour[Zidx(i)],1); //VContour[Zidx(i)];

#if DEBUG
            Source(i) = ContourASide[Zidx(i)];
#endif


        }
        // Else if we're just shy of the boundary
        else if( cPos(Ridx(i)+1,DeltaR) > VContour[Zidx(i)] )
        {

            #if DEBUG
            cout << "SHOULD NOT BE HERE!" << endl;
            #endif


            // Then we need to linearly interpolate A based on the current
            // value of A_ik and the contour value.
            // (We only interpolate in the R direction at the moment,
            //  could also interpolate up and down as well)

            // B-field at the contour
            double curB = sqrt(8.0*PI*1.0/betaedge);

            double leftpointA = curState[Apot][Ridx(i)-1][Zidx(i)];
            double rightpointA = 0.5*curB*pow(VContour[Zidx(i)],1); // analytic
            double leftpos = cPos(Ridx(i)-1,DeltaR);

            double m = (rightpointA - leftpointA)/(VContour[Zidx(i)] - leftpos);

            double midIA = leftpointA + m*(cPos(Ridx(i),DeltaR) - leftpos);

            // Now we want this midV to become the new value. This matrix
            // row becomes relatively simple. Overwrites the source vector
            rightV = 0; leftV  = 0; upV  = 0; downV  = 0; midV   = 1;
            Source(i) = midIA;
        }

        // Exception if we're at the upper or bottom-left corner
        // For the first round, A will be set up to be vertical, so the entire left boundary
        // will have both x and y derivatives as zero. So we'll interpolate only in the first solve.
        if( Ridx(i)==0 && DEBUG==0 && ( firstASolve==1 || ( Zidx(i)==0 || Zidx(i)==Z-1 ) ) )
        {
            #if DEBUG
            cout << "SHOULD NOT BE HERE!" << endl;
            #endif

            // Here we just interpolate the value of A.
            leftV = 0; rightV = 0; upV = 0; downV = 0; midV = 1;

            // Interpolation (quadratic a*R^2)
            //(assuming uniform grid in the R-direction and A=0 on left boundary).
            //Source(i) = (1.0/9.0)*curState[Apot][Ridx(i)+1][Zidx(i)] ;

            // Interpolation (linear a*R)
            Source(i) = (1.0/3.0)*curState[Apot][1][Zidx(i)] ;

        }


        // Now we fill the entries of the matrix with the values
        AmMatrix(i,i) = midV;
        if(Ridx(i)!=0) AmMatrix(i,i-1) = leftV;
        if(Ridx(i)!=N-1) AmMatrix(i,i+1) = rightV;
        if(Zidx(i)!=Z-1) AmMatrix(i,i+N) = upV;
        if(Zidx(i)!=0) AmMatrix(i,i-N) = downV;

    }



 // After the first time, Interpolation only happens for the corners
 if(firstASolve) firstASolve=0;

}


// Creates the matrix that is used
void createPoissonMatrix(double ***curState, MatrixXd& PoMatrix, VectorXd& Source)
{
    int i,j,k;


    // Source Vector = 4 pi q exp(-V_previous)
    for(i=0;i<(Z)*(N);i++)
        Source(i) = pow(DeltaR,2)*4.0*PI*curState[Q][Ridx(i)][Zidx(i)]*exp(-curState[Vpot][Ridx(i)][Zidx(i)]);


    // Loops over every cell to be solved, each which has it's own finite difference
    // scheme and fills one particular row of the big matrix A in   A u = f
    for(i=0;i<(Z)*(N);i++)
        {
            // Each of these cells gets a finite difference scheme

            double wi = DeltaR/2.0/cPos(Ridx(i),DeltaR);
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
                // Then the V_(k-1) term gets added to the V_k+1 term instead
                //midV = midV + downV;
                upV = upV + downV;
                downV = 0.0; // <-- Not really needed, there's no matrix entry
            }


            // Check to see if we're at the top
            if( Zidx(i) == Z-1 )
            {
                // Then the V_(k+1) term gets added to the V_k-1 term instead
                //midV = midV + upV;
                downV = downV + upV;
                upV = 0.0; // <-- Not really needed, there's no matrix entry
            }

            // Check to see if we're at the left
            if( Ridx(i) == 0 )
            {
                // Then the V_(i-1) term gets added to the V_i+1 term instead
                //midV = midV + leftV;
                rightV = rightV + leftV;
                leftV = 0.0; // <-- Not really needed, there's no matrix entry
            }

            //cout << "i=" << i << " " << cPos(Ridx(i)+1+1,DeltaR) << " " << VContour[Zidx(i)+1] << endl;

            // Check to see if we're near the right boundary
            // If we're beyond it, who cares, overwrite source vector
            if( cPos(Ridx(i),DeltaR) > VContour[Zidx(i)])
            {


                // One could check to see if this particular i value is never
                // access in any other row...
                rightV = 0;
                leftV = 0;
                upV = 0;
                downV = 0;
                midV = 1;
                Source(i) = ContourVpot; //VContour[Zidx(i)];

#if DEBUG
                Source(i) = ContourValues[Zidx(i)];
#endif


            }
            // Else if we're just shy of the boundary
            else if( cPos(Ridx(i)+1,DeltaR) >= VContour[Zidx(i)] )
            {


                // Then we need to linearly interpolate V based on the current
                // value of V_ik and the contour value.
                // (We only interpolate in the R direction at the moment,
                //  could also interpolate up and down as well)
                double leftpointV = curState[Vpot][Ridx(i)-1][Zidx(i)];
                double leftpos = cPos(Ridx(i)-1,DeltaR);
                double m = (ContourVpot - leftpointV)/(VContour[Zidx(i)] - leftpos);

#if DEBUG
                m = (ContourValues[Zidx(i)] - leftpointV)/(VContour[Zidx(i)] - leftpos);
#endif

                double midIV = leftpointV + m*(cPos(Ridx(i),DeltaR) - leftpos);

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
            if(Ridx(i)!=N-1) PoMatrix(i,i+1) = rightV;
            if(Zidx(i)!=Z-1) PoMatrix(i,i+(N)) = upV;
            if(Zidx(i)!=0) PoMatrix(i,i-(N)) = downV;

        }


}

void solvePoisson(double ***curState)
{
    double L2tol = 1e-10;

    /* We need to create the finite difference matrix and the source vector */

    // Allocate first using the Eigen structures, initalize to zero
    // The boundaries are not solved for
    MatrixXd PoMatrix = MatrixXd::Zero((Z)*(N),(Z)*(N));
    VectorXd Source = VectorXd::Zero((Z)*(N));
    VectorXd Soln = VectorXd::Zero((Z)*(N));


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

        if(0)
        {
            //Bottom-left corner
            curState[Vpot][0][0] = Soln[0];
            //Bottom row
            for(i=0;i<N;i++) curState[Vpot][i][0] = Soln[i];
            //Bottom-right corner (by construction should be ContourVpot)
            curState[Vpot][N-1][0] = ContourVpot;

#if DEBUG
            curState[Vpot][N-1][0] = ContourValues[0];
#endif

            //Left side
            for(i=0;i<Z;i++) curState[Vpot][0][i] = Soln[ i*(N) ];

            //Top-left corner (should be able to look either right or below and get
            //                 the same answer... could check this?)
            curState[Vpot][0][Z-1] = curState[Vpot][0][Z-2];
            //Top row
            for(i=0;i<N;i++) curState[Vpot][i][Z-1] = Soln[ i + (Z-1)*(Z)];
            //Top-right corner (dV/dz = 0)
            curState[Vpot][N-1][Z-1] = curState[Vpot][N-1][Z-2];

        }


        for(j=0;j<(Z);j++)
            for(i=0;i<(N);i++)
                curState[Vpot][i][j] = Soln[ i + j*(N) ];

        // Overall, who cares about what's right of the contour
        for(i=0;i<N;i++) for(j=0;j<Z;j++)
            if(cPos(i,DeltaR) >= VContour[j]) curState[Vpot][i][j]=ContourVpot;


        // Test
#if DEBUG
        for(i=0;i<N;i++) for(j=0;j<Z;j++)
                if(cPos(i,DeltaR) >= VContour[j]) curState[Vpot][i][j]=ContourValues[j];
#endif

    }

}


void solveMagnetic(double ***curState)
{
    double L2tol = 1e-10;

    /* We need to create the finite difference matrix and the source vector */

    // Allocate first using the Eigen structures, initalize to zero
    // The boundaries are not solved for
    MatrixXd AmMatrix = MatrixXd::Zero((Z)*(N),(Z)*(N));
    VectorXd Source = VectorXd::Zero((Z)*(N));
    VectorXd Soln = VectorXd::Zero((Z)*(N));


    // Then fill
    createAmpereMatrix(curState,AmMatrix,Source);

    cout << AmMatrix << endl;
    cout << Source << endl;

    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    Soln = AmMatrix.colPivHouseholderQr().solve(Source);

    //cout << endl << Soln << endl;

    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (AmMatrix*Soln - Source).norm() / Source.norm();
    cout << "Ampere Matrix L2 error = " << L2error << endl;

    if(L2error < L2tol)
    {
        // Copy new results to curState
        int i,j;

        if(0)
        {
            //Bottom-left corner
            curState[Apot][0][0] = Soln[0];
            //Bottom row
            for(i=0;i<N;i++) curState[Apot][i][0] = Soln[i];
            //Bottom-right corner (by construction should be ContourVpot)
            curState[Apot][N-1][0] = 1.0;

            //Left side
            for(i=0;i<Z;i++) curState[Apot][0][i] = Soln[ i*(N) ];

            //Top-left corner (should be able to look either right or below and get
            //                 the same answer... could check this?)
            curState[Apot][0][Z-1] = curState[Apot][0][Z-2];
            //Top row
            for(i=0;i<N;i++) curState[Apot][i][Z-1] = Soln[ i + (Z-1)*(Z)];
            //Top-right corner (dV/dz = 0)
            curState[Apot][N-1][Z-1] = curState[Apot][N-1][Z-2];

        }


        for(j=0;j<(Z);j++)
            for(i=0;i<(N);i++)
                curState[Apot][i][j] = Soln[ i + j*(N) ];

        // Overall, who cares about what's right of the contour
        //for(i=0;i<N;i++) for(j=0;j<Z;j++)
        //    if(cPos(i,DeltaR) > VContour[j]) curState[Apot][i][j]=1.0;


    }


}
