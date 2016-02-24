// Matrix routines
#include <Eigen/Dense>
#include "Filament_Globals.H"

using namespace Eigen;


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

    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    Soln = AmMatrix.colPivHouseholderQr().solve(Source);

    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (AmMatrix*Soln - Source).norm() / Source.norm();
    cout << "Ampere Matrix L2 error = " << L2error << endl;

    // Copy the solution to the newState vector
    for(i=0;i<M;i++) for(j=0;j<N;j++) newState[Apot][i][j] = Soln(Sidx(i,j));

    // Apot has been updated, we don't need to renormalize
};


void createAmpereMatrix(MatrixXd& AmMatrix, VectorXd& Source)
{
    int s;
    double Binf = sqrt(8.0*PI/betaInf);

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
            // We will manually set A = B_inf * r / 2
            Source(s) = Binf*cPos(M-1,DeltaR)/2.0;
        }
        else
        {
            // Need dQ/dPhi = Q(i+1)-Q(i-1)+Q(j+1)-Q(j-1) / Phi(i+1)-Phi(i-1)+Phi(j+1)-Phi(j-1)
            double Qright, Qleft, Qup, Qdown;
            double Phiright, Phileft, Phiup, Phidown;

            // Left will always be available.
            Qleft = curState[Q][Ridx(s)-1][Zidx(s)];
            if( Ridx(s) == M-1 ) Qright = Qleft;    // at or outside filament, Q is constant
            else Qright = curState[Q][Ridx(s)+1][Zidx(s)];

            // Since rho and V obey d/dz=0 boundary conditions, so does Q
            if( Zidx(s) == 0 ) { Qup = curState[Q][Ridx(s)][Zidx(s)+1]; Qdown = Qup;}
            else if(Zidx(s)==N-1) {Qdown = curState[Q][Ridx(s)][Zidx(s)-1]; Qup = Qdown;}
            else
            {
                Qup = curState[Q][Ridx(s)][Zidx(s)+1];
                Qdown = curState[Q][Ridx(s)][Zidx(s)-1];
            }


            // Left will always be available, right is not periodic, but uses analytic for straight uniform field
            // Up and down periodic
            Phileft = cPos(Ridx(s)-1,DeltaR) * curState[Apot][Ridx(s)-1][Zidx(s)];

            if( Zidx(s) == 0 )
            {
                Phiup = cPos(Ridx(s),DeltaR) * curState[Apot][Ridx(s)][Zidx(s)+1];
                Phidown = Phiup;
            }
            else if( Zidx(s) == N-1 )
            {
                Phidown = cPos(Ridx(s),DeltaR) * curState[Apot][Ridx(s)][Zidx(s)-1];
                Phiup = Phidown;
            }
            else
            {
                Phiup = cPos(Ridx(s),DeltaR) * curState[Apot][Ridx(s)][Zidx(s)+1];
                Phidown = cPos(Ridx(s),DeltaR) * curState[Apot][Ridx(s)][Zidx(s)-1];
            }

            if( Ridx(s) == M-1 )
            {
                // Analytic A = 0.5*Binf*r, beta = rho*cs2 / B^2/8pi -> beta = 8pi/B^2
                Phiright = 0.5*pow(cPos(Ridx(s),DeltaR),2)*Binf;
            }
            else
            {
                Phiright = cPos(Ridx(s)+1,DeltaR) * curState[Apot][Ridx(s)+1][Zidx(s)];
            }

            // dQ/dPhi
            double dQdPhi = (Qright-Qleft+Qup-Qdown)/(Phiright-Phileft+Phiup-Phidown);
            //cout << dQdPhi << endl;

            if( cPos(Ridx(s),DeltaR) <= VContour[Zidx(s)]) Source(s) = -0.5*pow(DeltaR,2)*betaInf*cPos(Ridx(s),DeltaR)*dQdPhi*exp(-curState[Vpot][Ridx(s)][Zidx(s)]);
            else Source(s) = 0;
        }
    }

    // Finite Difference Matrix
    for(s=0;s<M*N;s++)
    {
        double wi,f,Mleft,Mright,Mup,Mdown,Mcenter;

        if( Ridx(s) == 0 || Ridx(s) == M-1 )
        {
            // Left side will normalize A = 0, so all the entries here are 0
            // Right side normalize to A = Binf*r/2 (loaded into source vector)
            Mleft = 0;
            Mright = 0;
            Mup = 0;
            Mdown = 0;
            Mcenter = 1;
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
            Mcenter = -2-2*pow(f,2)-4*pow(wi,2);
        }

        // Now we will adjust to incorproate boundary conditions

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
