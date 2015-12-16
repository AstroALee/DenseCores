

/* Global Variables */
int N=1;
int Z=1;
double pLength=1;
double RhoCenter=1;
double MassGuess=1;
double DeltaR = 1.0;
double DeltaZ = 1.0;
double ContourQ = 1.0;
double ContourVpot = 1.0;
double betaedge = 1.0;
int firstASolve = 1;

/* Global States */
double ***curState;
double ***prevState;

/* Needed In Solving */
double *VContour;
double *RhoTop;


// Header that extern defines the above variables. Since they are defined
// for the first time above, the include statement needs to be here, below
// the declarations.
#include "DenseCoreGlobals.H"
#include "ErrorMessages.H"

#if DEBUG
double *ContourValues;
double *ContourATop;
double *ContourASide;
#endif

// Functions
void MagCylinder(double*,double*,double*);
void SolveAll(double ***curState);
void FillState(double ***curState);
void AllocateState(double ***&State);

void ConvergeTest(double ***curState, double ***prevState);
void CopyState(double ***curState, double ***prevState);
void RelaxSolution(double ***curState, double ***prevState);
double EllipticComplete(double k, int which);


double EllipticComplete(double k, int which)
{
    // Checks and balances
    if (k < 0 or k > 1) return 0; // throw error
    if (which < 0 or which > 1) return 0; // throw error

    double A[100],B[100];
    double e = 1e-7;
    double e1,e2;
    int j, m; double pi;
    pi = 4*atan(1);
    A[0]=1.0+k ; B[0]=1-k;
    int n=0;

    e100: n++;
    // Generate improved values
    A[n]=(A[n-1]+B[n-1])/2.0;
    B[n]=sqrt(A[n-1]*B[n-1]);
    if (fabs(A[n]-B[n]) > e) goto e100;
    e1=pi/2.0/A[n];
    e2=2.0;
    m=1;
    for (j=1; j<n+1; j++) {
      e2=e2-m*(A[j]*A[j]-B[j]*B[j]);
      m=m*2;
    }
    e2 *= (e1/2.0);

    if(which==1) return e1;
    else return e2;
}


// Linearly Interpolate an array onto a point. Returns value.
inline double LInt(double* Xarray, double* Yarray, double point, int size)
{
    double value = 0;
    int found = 0;

    // find the two points surrounding 'point' in the X array
    // Assumes X array is arranged in ascending value
    int i;
    for(i=1;i<size;i++) if( Xarray[i] >= point)
    {
        double RightY = Yarray[i];
        double LeftY  = Yarray[i-1];
        double RightX = Xarray[i];
        double LeftX  = Xarray[i-1];

        double m = (RightY-LeftY)/(RightX-LeftX);

        value = LeftY + m*(point-LeftX);
        found = 1;
        break;
    }

    if(!found)
    {
        cout << "Interpolated value was not found!" << endl;
        WATERLOO_ShouldNotBeHere("LInt in DenseCoreMain");
    }

    return value;
}


// Prints out the data to the screen
void PrintScreen(double ***curState)
{
    int i,j,k;

    for(j=0; j< N; j++) for(k=0;k<Z;k++)
    {
        cout << j << "," << k;
        for(i = 0; i< NStates; i++) cout << "," << curState[i][j][k];
        cout << endl;
    }
}


// Prints out the data to a file
void PrintV(double V[],int j)
{
    int i;

    ofstream myfile;
    myfile.open("VContour.txt");
    for(i=0;i<j;i++)
    {
        myfile << V[i] << endl;
    }
    myfile.close();

}

// Prints out the data to a file
void PrintFile(double ***curState)
{
    int i,j,k;

    ofstream myfile;
    myfile.open("test.txt");
    for(j=0; j< N; j++) for(k=0;k<Z;k++)
    {
        myfile << j << "," << k;
        for(i = 0; i< NStates; i++) myfile << "," << curState[i][j][k];
        myfile << endl;
    }
    myfile.close();



}


// Creates the initial A values everywhere in the grid.
void InitialA(double ***curState,double* initbeta)
{
    // Given the initial density distribution, impose a vertical magnetic
    // field structure.
    // B_inf * r / 2 = A  :: for uniform field everywhere.
    // B = sqrt(8*pi*rho/beta)

    int i,j;
    for(j=0;j<Z;j++)
    {
        int idone = 0;

        for(i=0;i<N;i++)
        {
            // Only fills spots left of the right-most V contour.
            if(idone==0 && cPos(i,DeltaR) >= VContour[j] ) idone = 1;

            if(idone)
            {
                double curB = sqrt(8.0*PI*1.0/betaedge);
                curState[Apot][i][j] = curB*VContour[j]/2.0 ;
            }
            else
            {
                //double currho = curState[Q][i][j]*exp(-curState[Vpot][i][j]);
                double currho = RhoTop[i];
                //double curbeta = initbeta[i];
                double cenrho = RhoTop[0];
                double cenbeta = initbeta[0];
                //double curB = sqrt(8.0*PI*currho/curbeta);
                double cenB = sqrt(8.0*PI*cenrho/cenbeta);
                double curR = cPos(i,DeltaR); // always i , not idx

                curState[Apot][i][j] =  0.5*cenB*sqrt(currho/cenrho)*curR; //curB*pow(curR,2)/2.0;

            }
        }
    }

// Resets magnetic potential to a force free version, interior to a ring of charge
// Also sets the upper and right contours
#if DEBUG

    // Radius of ring is a, in our units (must be bigger than entire domain)
    // Current is icur, mu_0 = 1  (both overall constants, should not matter in finite difference scheme)
    for(i=0;i<N;i++) for(j=0;j<Z;j++)
    {
        double a = 0.60;
        double icur = 2.0;

        double k = 4.0*a*cPos(i,DeltaR) / ( pow(cPos(j,DeltaZ),2) + pow( a + cPos(i,DeltaR),2) );
        k = sqrt(k);

        curState[Apot][i][j] = icur*sqrt(a/cPos(i,DeltaR))/k/PI * ( (1-0.5*k*k)*EllipticComplete(k,1) - EllipticComplete(k,2));
        cout << "(" << cPos(i,DeltaR) << ", " << cPos(j,DeltaZ) << ") k = " << k << " " << EllipticComplete(k,1) << " " << EllipticComplete(k,2) << endl;
    }

    // Sets boundary values
    for(i=0;i<Z;i++) ContourASide[i] = curState[Apot][N-1][i];
    for(i=0;i<N;i++) ContourATop[i] = curState[Apot][i][Z-1];

#endif



    cout << "Apot values = " ;
    for(j=0;j<N;j++) cout << curState[Apot][j][Z-1] << " " ;
    cout << endl;
}


// Creates the initial Q values everywhere in the grid.
void InitialQ(double ***curState,double initbeta[])
{
    int i,j;

    // To get the initial Q, we need the initial density distribution.

    // Creates the array that will hold the density
    double MagDensity[N];
    fill_n(MagDensity,N,0);

    // If we're debugging, reset the location of the contour.
#if DEBUG
    int WhichDebug = 1;

    if(WhichDebug==1)
    {
        // Overwrites location of counter to be the outermost cell
        for(j=0;j<Z;j++) VContour[j] = pLength - 0.5*DeltaR;
    }
    else if(WhichDebug==2)
    {
        // Contour location is a slanted line
        for(j=0;j<Z;j++) VContour[j] = pLength - (0.5*DeltaR)*((double)j+2.0)/2.0;
    }
#endif


    // Integrates the fluid equations so that the density is that of
    // an isothermal magnetized cylinder. Iterates until beta equals a
    // particular value at the boundary where the top of the domain intersects
    // the V=cst. surface.
    double an_pot[N];
    MagCylinder(MagDensity,initbeta,an_pot);
    for(i=0;i<N;i++) RhoTop[i] = MagDensity[i];


    // Test, resets potential
#if DEBUG


        if(WhichDebug==1)
        {
            // Overwrites potential with analytic answer for cylinder
            // Note we've added a fixed constant.
            for(i=0;i<N;i++) for(j=0;j<Z;j++) curState[Vpot][i][j] = 1.0 + an_pot[i];

            // Counter value is fixed.
            ContourVpot = curState[Vpot][N-1][Z-1];
            // In this case, we have an array of contour values. But for the cylinder
            // it is fixed.
            for(j=0;j<Z;j++) ContourValues[j] = ContourVpot;

            cout << "ContourV = " << ContourVpot << endl;
            //Print VContour to file
            PrintV(VContour,Z);
        }
        else if(WhichDebug==2)
        {

            // Constant added (an_pot used several places below)
            for(i=0;i<N;i++) an_pot[i] = an_pot[i] + 1.0;

            // Rewrites the potential to be the cylinder solution
            for(i=0;i<N;i++) for(j=0;j<Z;j++) curState[Vpot][i][j] = an_pot[i];

            // Doesn't really mean much in this debugging case.
            ContourVpot = curState[Vpot][N-1][Z-1];

            // Interpolate Contour Value for the slanted line.
            double Xvalue[N];
            for(i=0;i<N;i++) Xvalue[i] = cPos(i,DeltaR);
            for(j=0;j<Z;j++) ContourValues[j] = LInt(Xvalue,an_pot,VContour[j],N);

            // Re-write values right of these locations
            for(j=0;j<Z;j++) for(i=0;i<N;i++)
                if(cPos(i,DeltaR) >= VContour[j]) curState[Vpot][i][j] = ContourValues[j];

            // Right of contour, who cares
            //for(j=0;j<Z;j++) for(i=0;i<N;i++) if(cPos(i,DeltaR) > VContour[j]) curState[Vpot][i][j] = 0;

            cout << "Contour Values = " << endl;
            cout << ContourValues[0];
            for(j=1;j<Z;j++) cout << " " << ContourValues[j];
            cout << endl;

            PrintV(VContour,Z);
        }
#endif


    // Now with the initial density distribution, along with the initial V
    // we can determine the initial Q everywhere
    // Q == rho*c^2 exp(V/c^2) --> rho*exp(V)

    // Exact value of Q based on rho = 1 (definition) and Vpotential of contour
    ContourQ = 1.0*exp(ContourVpot);
    cout << "Contour Q = " << ContourQ << endl;

    for(j=0;j<Z;j++)
    {
        int idone = 0;

        for(i=0;i<N;i++)
        {
            // If left of this contour, use cell data. Else use fixed value (we'll never need it)
            curState[Q][i][j] = MagDensity[i]*exp(curState[Vpot][i][j]);
            if( cPos(i,DeltaR) >= VContour[j]) curState[Q][i][j] = ContourQ;

#if DEBUG
            // Why do I have this?
            if( cPos(i,DeltaR) >= VContour[j]) curState[Q][i][j] = ContourQ;
#endif

        }
    }


}




// Creates the potential resulting from an array of point masses equally
// spaced along the R=0 axis. Fills the potential of the current state
// with these values.
void PointMassPotential(double ***curState)
{
    int i,j;

    // G*Mguess
    double GM = MassGuess * GCST ;

    // Tolerance (keeps adding point masses above and below until
    // potential changes by less than this percent
    double tol = 0.001;

    // Max number of times to iterate through loop
    // make double since double will be used in comparision
    double nMax = 10000.0;

    // Previous state array
    double **PotPrevious;
    PotPrevious = new double *[N];
    for(i=0;i<N;i++) PotPrevious[i] = new double[Z];

    // Clears arrays
    for(i=0;i<N;i++) for(j=0;j<Z;j++) PotPrevious[i][j]=0.0;
    for(i=0;i<N;i++) for(j=0;j<Z;j++) curState[Vpot][i][j]=0.0;

    double Ppos = 0.0; // holds particle position
    double di = 0.0; // holds distance to particle

    // Loop
    while(true)
    {

        // Updates particle position
        Ppos = di*2.0*Z*DeltaZ;

        // Updates potential
        for(i=0;i<N;i++) for(j=0;j<Z;j++)
        {
            // Get position of current cell
            double CposR = cPos(i,DeltaR); //(0.5+(double)i)*DeltaR;
            double CposZ = cPos(j,DeltaZ); //(0.5+(double)j)*DeltaZ;

            double RelPos = sqrt( pow(CposR,2) + pow(CposZ-Ppos,2) );
            double RelPosN = sqrt( pow(CposR,2) + pow(CposZ+Ppos,2) );

            // Adds potential to each cell
            curState[Vpot][i][j] += -GM / RelPos;
            if(di != 0.0) curState[Vpot][i][j] += -GM / RelPosN ;

        }

        // Checks to see how much potential has changed from previous state
        // Loops over all particles and calculates the max percent error
        double error = 1.0;
        if(di != 0.0 ) for(i=0;i<N;i++) for(j=0;j<Z;j++)
        {
            double curError = fabs( (curState[Vpot][i][j] - PotPrevious[i][j])/curState[Vpot][i][j]);
            if(curError<error) error = curError;
        }

        //cout << "Iteration " << di << " has an error of " << error << endl;

        // If error is low enough, break
        if(error < tol)
        {
            cout << "Potential: Tolerence achieved after " << di << " iterations." << endl;
            cout << "Potential: The error was " << error << "." << endl;
            break;
        }


        // Else add some more particles and start again

        // Copy current potential to previous for error comparison
        for(i=0;i<N;i++) for(j=0;j<Z;j++)
            PotPrevious[i][j] = curState[Vpot][i][j];

        di += 1.0;

        if(di > nMax) WATERLOO_InitPotNotConverging(nMax,tol);


    }


    // All done

    // Clear previous potential memory
    for(i=0;i<N;i++) delete[] PotPrevious[i];
    delete[] PotPrevious;

    // Adjusts so the center-most cell has V=0
    //double Vmid = curState[Vpot][0][0];
    //for(i=0;i<N+1;i++) for(j=0;j<Z+1;j++) curState[Vpot][i][j] -= Vmid;

    // Determines R location of V=ContourVpot surface for each value of Z
    // After which grid cell on top boundary is V > ContourVpot ?

    ContourVpot = curState[Vpot][N-1][0];
    cout << "ContourVpot = " << ContourVpot << endl;

    for(j=0;j<Z;j++) for(i=0;i<N;i++)
        if(curState[Vpot][i][j] >= ContourVpot)
         {
             if(i>0)
             {
                 double m = (curState[Vpot][i][j] - curState[Vpot][i-1][j])/DeltaR;
                 double x1 = cPos(i-1,DeltaR); //0.5 + ((double)(i-1))*DeltaR ;

                 VContour[j] = x1 + (ContourVpot - curState[Vpot][i-1][j])/m;
                 break;
             }
             else
             {
                 // Probabably should never be here...
                 WATERLOO_ShouldNotBeHere("PointMassPotential");
                 VContour[j] = 0;
             }

         }


    cout << "VContour values: " << VContour[0];
    for(j=1;j<Z;j++) cout << "," << VContour[j];
    cout << endl;




    //Print VContour to file
    PrintV(VContour,Z);

}
