// Matrix routines
#include <Eigen/Dense>
#include "Filament_Globals.H"

using namespace Eigen;

void createPoissonMatrix(MatrixXd& PoMatrix, VectorXd& Source);

void SolvePoisson()
{
    int i,j;

    // First allocate the matrix and vectors for the finite difference scheme, source, and solution
    MatrixXd PoMatrix = MatrixXd::Zero((M)*(N),(M)*(N));
    VectorXd Source = VectorXd::Zero((M)*(N));
    VectorXd Soln = VectorXd::Zero((M)*(N));

    // The finite difference scheme is applied
    createPoissonMatrix(PoMatrix,Source);

    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    //Soln = PoMatrix.householderQr().solve(Source);
    //Soln = PoMatrix.colPivHouseholderQr().solve(Source);
    Soln = PoMatrix.fullPivHouseholderQr().solve(Source);

    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (PoMatrix*Soln - Source).norm() / Source.norm();
    cout << "Gravity Matrix L2 error = " << L2error << endl;

    //cout << (PoMatrix*Soln - Source) << endl;

    // Copy the solution to the newState vector
    for(i=0;i<M;i++) for(j=0;j<N;j++) newState[Vpot][i][j] = Soln[ Sidx(i,j) ];

    // Vpot has been updated, re-normalize so V = 0 in the upper-left corner
    double Vcorner = newState[Vpot][0][N-1];
    for(i=0;i<M;i++) for(j=0;j<N;j++) newState[Vpot][i][j] = newState[Vpot][i][j] - Vcorner;
};


void createPoissonMatrix(MatrixXd& PoMatrix, VectorXd& Source)
{

    int s;

    // First, the source vector :  4*pi*dr^2*rho = 4*pi*dr^2*Q*exp(-V)
    for(s=0;s<M*N;s++)
    {
        if( cPos(Ridx(s),DeltaR) <= VContour[Zidx(s)] )  // In the filament
            Source(s) = 4.0*PI*pow(DeltaR,2)*curState[Q][Ridx(s)][Zidx(s)]*exp(-curState[Vpot][Ridx(s)][Zidx(s)]);
        else                                            // Outside the filament, rho = 0
            Source(s) = 0.0;
    }

    // Second, loop through and fill up the finite difference matrix
    for(s=0;s<M*N;s++)
    {
        double wi,f,Mleft,Mright,Mup,Mdown,Mcenter;

        if( cPos(Ridx(s),DeltaR) == 0 )
        {
            // Useful numbers
            wi = 0;
            f  = DeltaR/DeltaZ;

            // r=0 version, derived using L'Hopital on the original DiffEq and employing dV/dr=0
            Mleft = 2.0;
            Mright = 2.0;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(2.0+pow(f,2));
        }
        else
        {
            // Useful numbers
            wi = DeltaR/2.0/cPos(Ridx(s),DeltaR);
            f  = DeltaR/DeltaZ;

            // Original values of matrix entries (may change on boundaries)
            Mleft = 1.0 - wi;
            Mright = 1.0 + wi;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(1.0+pow(f,2));
        }

        // Now we will change some of these values to employ the boundary conditions

        // The top boundary employs a dV/dz = 0 condition so at z=N-1, V(z=N) = V(r=N-2)
        if( Zidx(s) == N-1 ) { Mdown = Mdown + Mup; Mup = 0;} // there is no up


        // The bottom boundary employs a dV/dz = 0 condition so at z=0, V(z=-1) = V(z=1)
        if( Zidx(s) == 0 ) { Mup = Mup + Mdown;  Mdown = 0;} // there is no down


        // The left boundary employs a dV/dr = 0 condition, so V(r=-1) = V(r=1)
        if( Ridx(s) == 0 ) { Mright = Mright + Mleft ; Mleft = 0;}  // there is no left


        // The right boundary employs a dV/dr = f(r,z) condition, so at r=M-1: V(M) = V(M-2) + 2*dr*f(r,z)
        if( Ridx(s) == M-1 )
        {
            // Uses Poisson equation and a ghost zone
            //Mleft = Mleft + Mright;
            //Mright = 0;
            //Source(s) = Source(s) - (1.0+wi)*2.0*DeltaR*VRight[Zidx(s)];

            // Employs a 2nd order approximation to the first derivative at M-1
            Mleft = -4.0;
            Mright = 0; // doesn't exist
            Mup = 0;
            Mdown = 0;
            Mcenter = 1;
            PoMatrix(s,s-2) = 1; // Need two neighbor points to the left
            Source(s) = 6.0*DeltaR*VRight[Zidx(s)];

        }

        // Now we fill the entries of the matrix with the values
        PoMatrix(s,s) = Mcenter;
        if(Ridx(s)!=0) PoMatrix(s,s-1) = Mleft;
        if(Ridx(s)!=M-1) PoMatrix(s,s+1) = Mright;
        if(Zidx(s)!=N-1) PoMatrix(s,s+M) = Mup;
        if(Zidx(s)!=0) PoMatrix(s,s-M) = Mdown;

    }

    //cout << PoMatrix << endl << endl;
    //cout << Source << endl;
};
