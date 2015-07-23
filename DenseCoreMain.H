

/* Global Variables */
int N=1;
int Z=1;
double RhoCorner=1;
double RhoCenter=1;
double MassGuess=1;
double DeltaR = 1.0;
double DeltaZ = 1.0;
double RightContour = 1.0;
double betaedge = 1.0;

/* Global States */
double ***curState;

/* Needed In Solving */
double *VContour;


// Header that extern defines the above variables. Since they are defined
// for the first time above, the include statement needs to be here, below
// the declarations. 
#include "DenseCoreGlobals.H"
#include "ErrorMessages.H"


// Functions
void MagCylinder(double*,double*);
void PointMassPotential(double***);



// Includes boundary (n=0; all z from 0 up to and including Z & z=Z; all n from 0 up to and including N)
void AllocateState(double ***&curState)
{
    int i,j,k;
    curState = new double **[NStates];
    
    
    for(i=0; i< NStates; i++)
    {
        curState[i] = new double *[N+1];
        for(j=0; j< N+1; j++) curState[i][j] = new double[Z+1];
    }

}





// Prints out the data to the screen
void PrintScreen(double ***curState)
{
    int i,j,k;
    
    for(j=0; j< N+1; j++) for(k=0;k<Z+1;k++)
    {
        cout << j << "," << k;
        for(i = 0; i< NStates; i++) cout << "," << curState[i][j][k];
        cout << endl;
    }
}


// Prints out the data to a file
void PrintFile(double ***curState)
{
    int i,j,k;
    
    ofstream myfile;
    myfile.open("test.txt");
    for(j=0; j< N+1; j++) for(k=0;k<Z+1;k++)
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
    // B_0 * r / 2 = A
    // B_0 = sqrt(8*pi*rho/beta)
    
    int i,j;
    for(j=0;j<Z+1;j++)
    {
        int idone = 0;
        int iedge = 0;
        
        for(i=0;i<N+1;i++)
        {
            // Only fills spots left of the right-most V contour.
            if(idone==0) if( (0.5+(double)i)*DeltaR > VContour[j] )
            {
                idone = 1;
                iedge = i-1;
            }
            
            int idx = i;
            if(idone==1) idx = iedge;
            
            double currho = curState[Q][idx][j]*exp(-curState[Vpot][idx][j]);
            double curbeta = initbeta[idx];
            double curB = sqrt(8.0*PI*currho/curbeta);
            double curR = (0.5+(double)i)*DeltaR; // always i , not idx
            
            curState[Apot][i][j] = curB*pow(curR,2)/2.0;
            
        }
    }
    
    
}



// Creates the initial Q values everywhere in the grid.
void InitialQ(double ***curState,double initbeta[])
{
    // To get the initial Q, we need the initial density distribution.
    
    // Creates the array that will hold the density
    double MagDensity[N+1];
    fill_n(MagDensity,N+1,0);
    
    // Integrates the fluid equations so that the density is that of
    // an isothermal magnetized cylinder. Iterates until beta equals a
    // particular value at the boundary where the top of the domain intersects
    // the V=cst. surface.
    MagCylinder(MagDensity,initbeta);
    
    // Now with the initial density distribution, along with the initial V
    // we can determine the initial Q everywhere
    // Q == rho*c^2 exp(V/c^2)
    int i,j;
    for(j=0;j<Z+1;j++)
    {
        int idone = 0;
        int iedge = 0;
        
        for(i=0;i<N+1;i++)
        {
            // Only fills spots left of the right-most V contour.
            if(idone==0) if( (0.5+(double)i)*DeltaR > VContour[j] )
            {
               idone = 1;
               iedge = i-1;
            }
            
            // If left of this contour, use cell data. Else copy edge data.
            if(idone==0) curState[Q][i][j] = MagDensity[i]*exp(curState[Vpot][i][j]);
            else curState[Q][i][j] = MagDensity[iedge]*exp(curState[Vpot][iedge][j]);
            
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
    PotPrevious = new double *[N+1];
    for(i=0;i<N+1;i++) PotPrevious[i] = new double[Z+1];
    
    for(i=0;i<N+1;i++) for(j=0;j<Z+1;j++) PotPrevious[i][j]=0.0;
    for(i=0;i<N+1;i++) for(j=0;j<Z+1;j++) curState[Vpot][i][j]=0.0;
    
    // Loop
    
    double Ppos = 0.0;
    double di = 0.0;
    
    while(true)
    {
        
        // Updates particle position
        Ppos = di*2.0*(1+Z)*DeltaZ;
        
        // Updates potential
        for(i=0;i<N+1;i++) for(j=0;j<Z+1;j++)
        {
            // Get position of current cell
            double CposR = (0.5+(double)i)*DeltaR;
            double CposZ = (0.5+(double)j)*DeltaZ;
            
            double RelPos = sqrt( pow(CposR,2) + pow(CposZ-Ppos,2) );
            double RelPosN = sqrt( pow(CposR,2) + pow(CposZ+Ppos,2) );
            
            // Adds potential to each cell
            curState[Vpot][i][j] += -GM / RelPos;
            if(di != 0.0) curState[Vpot][i][j] += -GM / RelPosN ;
            
        }
        
        // Checks to see how much potential has changed from previous state
        // Loops over all particles and calculates the max percent error
        double error = 1.0;
        if(di != 0.0 ) for(i=0;i<N+1;i++) for(j=0;j<Z+1;j++)
        {
            double curError = fabs( (curState[Vpot][i][j] - PotPrevious[i][j])/curState[Vpot][i][j]);
            if(curError<error) error = curError;
        }
        
        
        // If error is low enough, break
        if(error < tol) {cout << "Tolerence achieved after " << di << " iterations." << endl; break;}
        
        // Else add some more particles and start again
        for(i=0;i<N+1;i++) for(j=0;j<Z+1;j++)
            PotPrevious[i][j] = curState[Vpot][i][j];
        
        
        di += 1.0;
        
        if(di > nMax) WATERLOO_InitPotNotConverging(nMax,tol);
    }
    
    
    // All done
    for(i=0;i<N;i++) delete[] PotPrevious[i];
    delete[] PotPrevious;
    
    // Adjusts so the center-most cell has V=0
    double Vmid = curState[Vpot][0][0];
    for(i=0;i<N+1;i++) for(j=0;j<Z+1;j++) curState[Vpot][i][j] -= Vmid;
    
    // Determines R location of V=RightContour surface for each value of Z
    // After which grid cell on top boundary is V > RightContour ?
    RightContour = curState[Vpot][N][0];
    for(j=0;j<Z+1;j++) for(i=0;i<N+1;i++)
        if(curState[Vpot][i][j] >= RightContour)
         {
             double m = (curState[Vpot][i][j] - curState[Vpot][i-1][j]);
             double x1 = (0.5 + (double)(i-1))*DeltaR;
             double x2 = (0.5 + (double) i)*DeltaR;
             m = m/(x2-x1);
             
             VContour[j] = x1 + (RightContour - curState[Vpot][i-1][j])/m;
             break;
         }
    
    cout << VContour[0];
    for(j=1;j<Z+1;j++) cout << "," << VContour[j];
    cout << endl;
    
}

