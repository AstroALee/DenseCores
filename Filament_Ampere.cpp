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
    cout << "Gravity Matrix L2 error = " << L2error << endl;

    // Copy the solution to the curState vector
    for(i=0;i<M;i++) for(j=0;j<N;j++) curState[Apot][i][j] = Soln[ Sidx(i,j) ];

    // Apot has been updated, we don't need to renormalize
};


void createAmpereMatrix(MatrixXd& AmMatrix, VectorXd& Source)
{

};
