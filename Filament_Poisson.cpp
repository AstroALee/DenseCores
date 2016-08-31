// Matrix routines
#include <Eigen/Dense>
#include "Filament_Globals.H"

using namespace Eigen;

void createPoissonMatrix(MatrixXd& PoMatrix, VectorXd& Source);
void createPertPoissonMatrix(MatrixXd& PoMatrix, VectorXd& Source);

double UpdateDVDR(double r, double z);

void SolvePertPoisson(double** deltaV)
{
    int i,j;

    for(i=0;i<M;i++) for(j=0;j<N;j++) deltaV[i][j] = 0;

    // First allocate the matrix and vectors for the finite difference scheme, source, and solution
    MatrixXd PoMatrix = MatrixXd::Zero((M)*(N),(M)*(N));
    VectorXd Source = VectorXd::Zero((M)*(N));
    VectorXd Soln = VectorXd::Zero((M)*(N));

    // The finite difference scheme is applied
    createPertPoissonMatrix(PoMatrix,Source);

    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    //Soln = PoMatrix.householderQr().solve(Source);
    //Soln = PoMatrix.colPivHouseholderQr().solve(Source);
    Soln = PoMatrix.colPivHouseholderQr().solve(Source);

    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (PoMatrix*Soln - Source).norm() / Source.norm();
    cout << "Pert Gravity Matrix L2 error = " << L2error << endl;

    //cout << (PoMatrix*Soln - Source) << endl;

    // Copy the solution to the newState vector
    for(i=0;i<M;i++) for(j=0;j<N;j++) deltaV[i][j] = Soln[ Sidx(i,j) ];

    // Make the upper-left corner 0 again
    double Vnorm = deltaV[0][N-1];
    for(i=0;i<M;i++) for(j=0;j<N;j++) deltaV[i][j] = deltaV[i][j] - Vnorm;


}

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
    Soln = PoMatrix.colPivHouseholderQr().solve(Source);

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


void createPertPoissonMatrix(MatrixXd& PoMatrix, VectorXd& Source)
{
    int s;

    double w = pow(DeltaR/DeltaZ,2);
    double delta = 1;



    // First, the source vector : 4*Pi*exp(-V)*q
    for(s=0;s<M*N;s++)
    {
        if( cPos(Ridx(s),DeltaR) > VContour[Zidx(s)] ) delta = 0; else delta = 1;

        Source(s) = 0;// 4.0*PI*delta*pow(DeltaR,2)*curState[Q][Ridx(s)][Zidx(s)]*exp(-curState[Vpot][Ridx(s)][Zidx(s)]);

    }


    // Then the finite difference matrix
    for(s=0;s<M*N;s++)
    {
        if( cPos(Ridx(s),DeltaR) > VContour[Zidx(s)] ) delta = 0; else delta = 1;

        double alp,Mleft,Mright,Mup,Mdown,Mcenter;

        if( cPos(Ridx(s),DeltaR) == 0 )
        {
            Mleft = 0;
            Mright = 4.0;
            Mup = w;
            Mdown = w;
            Mcenter = 4.0*PI*delta*pow(DeltaR,2)*curState[Q][Ridx(s)][Zidx(s)]*exp(-curState[Vpot][Ridx(s)][Zidx(s)]) - 4.0 - 2.0*w;

        }
        else
        {
            alp = DeltaR/2.0/cPos(Ridx(s),DeltaR);

            Mleft = 1-alp;
            Mright = 1+alp;
            Mup = w;
            Mdown = w;
            Mcenter = 4.0*PI*delta*pow(DeltaR,2)*curState[Q][Ridx(s)][Zidx(s)]*exp(-curState[Vpot][Ridx(s)][Zidx(s)]) - 2.0 - 2.0*w ;

        }


        // Apply boundary conditions
        // --------------------------
        // The top boundary employs a dV/dz = 0 condition so at z=N-1, V(z=N) = V(r=N-2)
        if( Zidx(s) == N-1 ) { Mdown = Mdown + Mup; Mup = 0;} // there is no up
        // The bottom boundary employs a dV/dz = 0 condition so at z=0, V(z=-1) = V(z=1)
        if( Zidx(s) == 0 ) { Mup = Mup + Mdown;  Mdown = 0;} // there is no down
        // The left boundary employs a dV/dr = 0 condition, so V(r=-1) = V(r=1)
        if( Ridx(s) == 0 ) { Mright = Mright + Mleft ; Mleft = 0;}  // there is no left
        // The right boundary changes with every solve
        if( Ridx(s) == M-1 ) // use delta M = Mexcess_new - Mexcess_old
        {
            // Let's assume Mexcess isn't changing. Then the boundary condition
            // is just the cylinder.
            // Employs a 2nd order approximation to the first derivative at M-1

            Mleft = -4.0;
            Mright = 0; // doesn't exist
            Mup = 0;
            Mdown = 0;
            Mcenter = 1;
            PoMatrix(s,s-2) = 1; // Need two neighbor points to the left

            double rad = sqrt( pow(cPos(Ridx(s),DeltaR),2) + pow(cPos(Zidx(s),DeltaZ),2) );
            Source(s) = UpdateDVDR(cPos(Ridx(s),DeltaR),cPos(Zidx(s),DeltaZ)); //  ;testCst*mEx/pow(rad,2) ; //Ccst[0]/rL ;


        }



        // Now we fill the entries of the matrix with the values
        PoMatrix(s,s) = Mcenter;
        if(Ridx(s)!=0) PoMatrix(s,s-1) = Mleft;
        if(Ridx(s)!=M-1) PoMatrix(s,s+1) = Mright;
        if(Zidx(s)!=N-1) PoMatrix(s,s+M) = Mup;
        if(Zidx(s)!=0) PoMatrix(s,s-M) = Mdown;

    }


}

void createPoissonMatrix(MatrixXd& PoMatrix, VectorXd& Source)
{

    int s;

    int sMax = 0;
    double RhoMax = 0;

    // First, the source vector :  4*pi*dr^2*rho = 4*pi*dr^2*Q*exp(-V)
    for(s=0;s<M*N;s++)
    {
        if( cPos(Ridx(s),DeltaR) <= VContour[Zidx(s)] )  // In the filament
        {

                //Source(s) = 4.0*PI*pow(DeltaR,2)*curState[Q][Ridx(s)][Zidx(s)]*exp(-curState[Vpot][Ridx(s)][Zidx(s)]);
                Source(s) = 4.0*PI*pow(DeltaR,2)*curState[Rho][Ridx(s)][Zidx(s)];
                //cout << "Rho test:  Rho = " << curState[Rho][Ridx(s)][Zidx(s)] << " q e^-V = " << curState[Q][Ridx(s)][Zidx(s)]*exp(-curState[Vpot][Ridx(s)][Zidx(s)]) << endl;
                //if(fabs(curState[Rho][Ridx(s)][Zidx(s)]-curState[Q][Ridx(s)][Zidx(s)]*exp(-curState[Vpot][Ridx(s)][Zidx(s)]))>RhoMax) {sMax = s; RhoMax = fabs(curState[Rho][Ridx(s)][Zidx(s)]-curState[Q][Ridx(s)][Zidx(s)]*exp(-curState[Vpot][Ridx(s)][Zidx(s)]));}

        }
        else                                            // Outside the filament, rho = 0
            Source(s) = 0.0;
    }

    //cout << "sMax is = " << sMax << " or (r,z)=" << Ridx(sMax) << " , " << Zidx(sMax) << endl;
    //cout << "Where Delta Rho = " << RhoMax << endl;
    //cout << "Rho = " << curState[Rho][Ridx(s)][Zidx(s)] << " and qexpV = " << curState[Q][Ridx(s)][Zidx(s)]*exp(-curState[Vpot][Ridx(s)][Zidx(s)]) << endl;

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



double UpdateDVDR(double r, double z)
{
    double dm = mCurrent - mPrevious;
    if(dm<0) dm = 0;

    // Debug
    //dm = testCst;

    int LoopIdx;

    // finds dV/dr from point masses with mass equal to dm at location (r,z)
    // Uses same number of points as used in initial conditions.

    double dvdr = 0;
    double zP = 0;

    for(LoopIdx=0;LoopIdx<=PointLoopMax;LoopIdx++)
    {
            double X1 = sqrt( pow(r,2) + pow(z-zP,2) );
            double X2 = sqrt( pow(r,2) + pow(z+zP,2) );

            dvdr = dvdr + dm/pow(X1,2);
            if(zP>0) dvdr = dvdr + dm/pow(X2,2);

            zP = zP + 2.0*zL;
    }



    return dvdr;
};
