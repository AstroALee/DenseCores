// Matrix routines
#include <Eigen/Dense>
#include "Filament_Globals.H"

using namespace Eigen;

// What boundary condition are we using at the right edge
int BdyRight = 2; // 0 = value of A, 1 = value of dA/dr, 2 = value of Cross(A) = (1/r)*d(rA)/dr


void createAmpereMatrix(MatrixXd& AmMatrix, VectorXd& Source);

void SolveAmpere()
{
    int i,j;

    // First allocate the matrix and vectors for the finite difference scheme, source, and solution
    MatrixXd AmMatrix = MatrixXd::Zero((M)*(N),(M)*(N));
    VectorXd Source = VectorXd::Zero((M)*(N));
    VectorXd Soln = VectorXd::Zero((M)*(N));

    // The finite difference scheme is applied
    createAmpereMatrix(AmMatrix,Source);
    //cout << AmMatrix << endl;
    //cout << Source << endl;

    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    Soln = AmMatrix.colPivHouseholderQr().solve(Source);

    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (AmMatrix*Soln - Source).norm() / Source.norm();

    //cout << (AmMatrix) << endl;
    //cout << Source << endl;
    //cout << (AmMatrix*Soln - Source) << endl;

    cout << "Ampere Matrix L2 error = " << L2error << endl;

    // Copy the solution to the newState vector
    for(i=0;i<M;i++) for(j=0;j<N;j++) newState[Apot][i][j] = Soln(Sidx(i,j));

    // Apot has been updated, we don't need to renormalize
};


void createAmpereMatrix(MatrixXd& AmMatrix, VectorXd& Source)
{
    int s;

    double Binf = sqrt(8.0*PI/beta0)*pow(Rbdy,nCyl);
    double alpha = sqrt(8*PI/beta0); //Binf; // alpha^2 = 8*pi/beta_inf = 8*pi*(B^2/8*pi)/(rho_*c^2)_inft -> B^2
    double Ccst = PI*RhoTop[0]/2.0/(1.0+1.0/beta0); // RhoTop[0] should be 1

    // Source vector
    for(s=0;s<M*N;s++)
    {
        if( Ridx(s) == 0 )
        {
            // We will manually set A = 0 on the left hand side
            Source(s) = 0.0;
        }
        else if( Ridx(s) == M-1 )
        {
            if(BdyRight==0) // setting A
            {
                // We will manually set A = B_inf * r / 2
                Source(s) = Binf*cPos(M-1,DeltaR)/2.0;

                // Manually set to cylinder + empty space (the cylinder solution assumes nCyl = 0.5)
                double Redge = VContour[N-1];
                Source(s) = (0.5*alpha*sqrt(RhoTop[0])/Ccst)*log(Ccst*Redge*Redge+1.0)/rL + 0.5*Binf*(rL-Redge*Redge/rL);
            }
            else if(BdyRight==1) // setting dA/dr
            {

            }
            else // setting cros(A)
            {
                // We will impose a boundary condition (1/r)*d(rA)/dr = B_infty at this boundary
                Source(s) = 2.0*DeltaR*Binf;
            }
        }
        else
        {
            double dQdPhi = curState[dQdP][Ridx(s)][Zidx(s)];

            if( cPos(Ridx(s),DeltaR) <= VContour[Zidx(s)]) Source(s) = -4.0*PI*pow(DeltaR,2)*cPos(Ridx(s),DeltaR)*dQdPhi*exp(-curState[Vpot][Ridx(s)][Zidx(s)]);
            else Source(s) = 0;

        }
    }



    //cout << Source << endl;
    //cout << DeltaR << endl;

    // Finite Difference Matrix
    for(s=0;s<M*N;s++)
    {
        double wi,f,Mleft,Mright,Mup,Mdown,Mcenter;

        if( Ridx(s) == 0 )
        {
            // Left side will normalize A = 0, so all the entries here are 0
            Mleft = 0;
            Mright = 0;
            Mup = 0;
            Mdown = 0;
            Mcenter = 1;
        }
        else if(Ridx(s)==M-1)
        {
            if(BdyRight==0)
            {
                // Imposing value of A
                Mleft = 0;
                Mright = 0;
                Mup = 0;
                Mdown = 0;
                Mcenter = 1.0;
            }
            else if(BdyRight==1)
            {

            }
            else
            {
                // Right side implemented a Robin boundary condition
                Mleft = -4.0;
                Mright = 0; // doesn't exist
                Mup = 0;
                Mdown = 0;
                Mcenter = ( 2.0*DeltaR/cPos(M-1,DeltaR) + 3.0 );
                AmMatrix(s,s-2) = 1.0;
            }
        }
        else
        {
            // Base values
            wi = DeltaR/2.0/cPos(Ridx(s),DeltaR);
            f = DeltaR/DeltaZ;

            Mright = 1.0+wi;
            Mleft = 1.0-wi;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0-2.0*pow(f,2)-4.0*pow(wi,2);
        }

        // Now we will adjust to incorproate z-boundary conditions

        // at top and bottom, dA/dz = 0
        if( Zidx(s) == N-1)
        {
            Mdown = Mdown + Mup;
            Mup = 0;
        }

        if( Zidx(s) == 0 )
        {
            Mup = Mup + Mdown;
            Mdown = 0;
        }

        // Now we fill the entries of the matrix with the values
        AmMatrix(s,s) = Mcenter;
        if(Ridx(s)!=0) AmMatrix(s,s-1) = Mleft;
        if(Ridx(s)!=M-1) AmMatrix(s,s+1) = Mright;
        if(Zidx(s)!=N-1) AmMatrix(s,s+M) = Mup;
        if(Zidx(s)!=0) AmMatrix(s,s-M) = Mdown;
    }

};
