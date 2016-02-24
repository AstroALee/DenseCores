
/*
   The global variables need to be declared here, outside of main(). You can assign them a value
   in the code
 */
int M;
int N;
double rRatio;
double mExcess;
double betaInf;
double nCyl;

/* Declaration of global derived quantities */
double zL,rL;
double DeltaR;
double DeltaZ;

/* Filament Boundary and value */
double *VContour;
double Vbdy = 0;

/* Rho at the top boundary */
double *RhoTop;

/* Ambient background */
double Pbackground = 4.14194556e-13; //1.38064852e-12; // cgs Pressure (1000 cm-3 , 10 K)
double Cbackground = 18545.4363;     // cgs Sound speed (10 K, 2.4 mmw)

/* Boundary Condition Arrays for values or derivatives */
double *VLeft, *VRight;
double *ALeft, *ARight;

/* Global States */
double ***newState;
double ***curState;
double ***prevState;

/* Other values */

/* Unit coversion from parsecs to code length units (G=c_s=rho(edge)=1) */
double Pc2Code;
double Pc2Cm = 3.086e18;

/* Unit conversion from solar masses to code units */
double Sol2G = 1.989e33;
double Sol2Code;

/* Output File Name */
string OutFileName = "Output.out";

/* Variables related to the cylinder solve */
int LoopMAX = 100;
double CylTol = 1e-8;
double Ccst[2];

/* Variables related to the point potential solve */
int PointLoopMax = 100000;
double PointLoopTol = 1e-3;

/* Variables related to the convergence solve */
int ConvergeLoopMax = 1; // Set by LOOPMAX
double ConvergeTol = 1e-4;
double relaxFrac = 0.1;



// How long did all of this take?
void ExecutionTime(clock_t start, clock_t end )
{
    cout << "Execution took " << double( end - start ) / (double)CLOCKS_PER_SEC << " seconds." << endl;
};

void Tallies(clock_t startTime)
{


    // Stop the clock
    ExecutionTime(startTime,clock());

};


// Did you provide the right number of arguments at the command line?
void CheckArguments(int a)
{
    // The number should be correct
    if(a-1 < 7 || a-1 > 8)
    {
        WaterlooHeader("CheckArguments");
        cout << "You didn't provide the right number of arguments. You need seven or eight." << endl;
        cout << "(M,N,z-length (pc), Ratio r-length to z-length,Non-Dim Mass excess, beta_inf, n, (loopmax; opt))" << endl;
        cout << "You provided " << a-1 << " argument" << Plural(a-1,"s") << "." << endl;
        exit(1);
    }

};

// Read in the arguments and store them in the appropriate global variables.
void ReadArguments(int argc, char *a[])
{
    M = atoi(a[1]);
    N = atoi(a[2]);
    zL = atof(a[3]);
    rRatio = atof(a[4]);
    mExcess = atof(a[5]);
    betaInf = atof(a[6]);
    nCyl = atof(a[7]);
    if(argc-1==8) ConvergeLoopMax = atoi(a[8]);



    // Check that some numbers are reasonable.
    if(M < 3 || N < 3)
    {
        WaterlooHeader("ReadArguments");
        cout << "The number of grid cells is weird. You provided (M,N) = (" << M << "," << N << ")" << endl;
        cout << "You need at least 3 in each dimension" << endl;
        exit(1);
    }

    if(zL < 0 || rRatio < 1 || mExcess < 0)
    {
        WaterlooHeader("ReadArguments");
        cout << "Some values of zL, rRatio, or mExcess are weird. You provided zL, rRatio, mExcess = " << zL << ", " << rRatio << ", " << mExcess << endl;
        cout << "For example, the code assumes that edge of the filament is at most at the edge of the box (requires rRatio >= 1)" << endl;
        exit(1);
    }

    if(betaInf < 0 || nCyl < 0.5)
    {
        WaterlooHeader("ReadArguments");
        cout << "Your betaInf or nCyl values may be unphysical. You provided betaInf, nCyl = " << betaInf << ", " << nCyl << endl;
        exit(1);
    }

};

void WelcomeMessage()
{
    cout << endl;
    LongLine();
    cout << "Welcome to the Magnetized Dense Core Equilibrium code" << endl;
    cout << "Written by Aaron Lee, 2016" << endl;
    LongLine();
    cout << "Have an astronomical time!" << endl << endl;

};

void PrintArguments(int argc)
{
    cout << "Here are your input arguments:" << endl;
    cout << "M = " << M << endl;
    cout << "N = " << N << endl;
    cout << "zL = " << zL << " parsecs" << endl;
    cout << "rRatio = " << rRatio << endl;
    cout << "mExcess = " << mExcess << endl;
    cout << "betaInf = " << betaInf << endl;
    cout << "nCyl = " << nCyl << endl << endl;
    cout << "Here are your compiler options:" << endl;
    cout << "DEBUG = " << DEBUG << endl;
    cout << "LOOPMAX = " << LOOPMAX << endl;
    cout << "UNIFORM = " << UNIFORM << endl;
    cout << "RADIALRATIO = " << RADIALRATIO << endl;
    if(argc-1==8) cout << "Overwriting LOOPMAX to " << ConvergeLoopMax << endl;
};

void CalcDerived(int argc)
{
    // Unit conversion, multiply by this to convert parsecs to code units
    double cgsCodeLength = Cbackground*Cbackground/sqrt(Gcst*Pbackground);  // cm/code length
    double pcCodeLength = cgsCodeLength/Pc2Cm; // pc/code length
    Pc2Code = 1/pcCodeLength; // code/pc length

    // Unit conversion, multiply by this to convert solar masses to code units
    double cgsCodeMass = pow(Cbackground,4)/sqrt(Pbackground*pow(Gcst,3));
    double solCodeMass = cgsCodeMass/Sol2G;
    Sol2Code = 1/solCodeMass;

    cout << "Unit conversions:" << endl;
    cout << "Parsec to Code = " << Pc2Code << endl;
    cout << "Solar mass to Code = " << Sol2Code << endl;

    // Radial length
    rL = rRatio * zL;
    cout << "rL = " << rL << " parsecs" << endl;

    // Now converting to code units;
    zL = zL*Pc2Code;
    rL = rL*Pc2Code;

    // Cell sizes (in code units)
    if(UNIFORM) M = rRatio * N;   // Adjusts the value of M so DeltaR = DeltaZ
    if(UNIFORM) cout << "Uniform cell size adjustment: Adjusted the number of cells in the radial direction to " << M << "!" << endl;
    DeltaR = rL / (double(M)-1);
    DeltaZ = zL / (double(N)-1); // Need the -1 to include the edge as a grid point

    // Boundary condition arrays
    VLeft = new double[N];
    VRight = new double[N];
    ALeft = new double[N];
    ARight = new double[N];

    // Filament Boundary
    VContour = new double[N];

    // Rho at top
    RhoTop = new double[M];

    // Are we only doing one converge loop? (currently value potentially set by command line arg (CLA))
    if(DEBUG==1 && argc-1!=8) ConvergeLoopMax=1; // w/o CLA, override to only solve once
    if(LOOPMAX>=0 && argc-1!=8) ConvergeLoopMax = LOOPMAX; // w/o CLA, use GNUMake value for LoopMax
    if(ConvergeLoopMax<0 || DEBUG==-1) ConvergeLoopMax = 0; // Sanity or all-power DEBUG override

};


void AllocateState(double ***&State)
{
    int i,j,k;
    State = new double **[NStates];

    for(i=0; i< NStates; i++)
    {
        State[i] = new double *[M];
        for(j=0; j< M; j++) State[i][j] = new double[N];
    }

    // In case it's not initialized to zero
    for(i=0;i<NStates;i++) for(j=0;j<M;j++) for(k=0;k<N;k++) State[i][j][k] = 0.0;
};

void DeallocateState(double ***&State)
{
    delete[] State; // not sure if this is what I want it to to do...
};


// Prints out the data to a file
void PrintState(double ***curState)
{
    int i,j,k;

    ofstream myfile;
    myfile.open(OutFileName.c_str());

    // First print inputs as the first line
    myfile << M << "," << N << "," << zL << "," << rRatio << "," << mExcess << "," << betaInf << "," << nCyl << "," << Vbdy << endl;

    // Print VContour (converting to parsecs) as the second line
    myfile << VContour[0]/Pc2Code;
    for(j=1;j<N;j++) myfile << "," << VContour[j]/Pc2Code;
    myfile << endl;

    for(j=0; j< M; j++) for(k=0;k<N;k++)
    {
        double Rpc = cPos(j,DeltaR)/Pc2Code;
        double Zpc = cPos(k,DeltaZ)/Pc2Code;

        myfile << j << "," << k << "," << Rpc << "," << Zpc;
        for(i = 0; i< NStates; i++) myfile << "," << curState[i][j][k];
        myfile << endl;
    }
    myfile.close();
};

void CrudeIntegral(double& CylMass, double& TotMass)
{
    // Approximates the integral int(f dxdy) as DeltaX*DeltaY*Sum(f(xi,yi))

    int i,j;
    double sum=0;

    // Cylinder mass
    for(i=0;i<M;i++) for(j=0;j<N;j++) sum = sum + cPos(i,DeltaR)*RhoTop[i];
    CylMass = 4.0*PI*DeltaR*DeltaZ*sum;

    // Total mass
    sum = 0;
    for(i=0;i<M;i++) for(j=0;j<N;j++) sum = sum + cPos(i,DeltaR)*curState[Rho][i][j];
    TotMass = 4.0*PI*DeltaR*DeltaZ*sum;

};

void TrapIntegral(double& CylMass, double& TotMass)
{
    // Approximates the integral int(f dxdy) as 0.25*DeltaX*DeltaY*Sum(w*f(xi,yi))
    // where w = 1 for corners, 2 for edges, 4 for interior

    int i,j; double sum=0;

    //Cylinder
    sum = cPos(0,DeltaR)*RhoTop[0] + cPos(0,DeltaR)*RhoTop[0] + cPos(M-1,DeltaR)*RhoTop[M-1] + cPos(M-1,DeltaR)*RhoTop[M-1]; // corners
    for(i=1;i<M-1;i++) sum = sum + 2.0*(cPos(i,DeltaR)*RhoTop[i] + cPos(i,DeltaR)*RhoTop[i]); // top and bottom edges
    for(j=1;j<N-1;j++) sum = sum + 2.0*(cPos(0,DeltaR)*RhoTop[0] + cPos(M-1,DeltaR)*RhoTop[M-1]); // left and right edges
    for(i=1;i<M-1;i++) for(j=1;j<N-1;j++) sum = sum + 4.0*(cPos(i,DeltaR)*RhoTop[i]); // interiors
    CylMass = 4.0*PI*(DeltaR*DeltaZ/4.0)*sum;

    //Total Mass
    sum = cPos(0,DeltaR)*curState[Rho][0][0] + cPos(M-1,DeltaR)*curState[Rho][M-1][0] + cPos(0,DeltaR)*curState[Rho][0][N-1] + cPos(M-1,DeltaR)*curState[Rho][M-1][N-1]; // corners
    for(i=1;i<M-1;i++) sum = sum + 2.0*(cPos(i,DeltaR)*curState[Rho][i][0] + cPos(i,DeltaR)*curState[Rho][i][N-1]); // top and bottom edges
    for(j=1;j<N-1;j++) sum = sum + 2.0*(cPos(0,DeltaR)*curState[Rho][0][j] + cPos(M-1,DeltaR)*curState[Rho][M-1][j]); // left and right edges
    for(i=1;i<M-1;i++) for(j=1;j<N-1;j++) sum = sum + 4.0*(cPos(i,DeltaR)*curState[Rho][i][j]); // interiors
    TotMass = 4.0*PI*(DeltaR*DeltaZ/4.0)*sum;


};


void SimpIntegral(double& CylMass, double& TotMass)
{
    // Approximates the integral int(f dxdy) as (1/9)*DeltaX*DeltaY*Sum(w*f(xi,yi))
    // where w has a more involved formulation

    int i,j; double sum=0;

    // Since we require an odd number of points, we'll copy the rho data to an array that works for us
    int Mhere,Nhere;
    if(M%2==0) Mhere = M+1; else Mhere=M;
    if(N%2==0) Nhere = N+1; else Nhere=N;

    double RhoData[Mhere][Nhere];
    double RhoDataTop[Mhere];

    for(i=0;i<M;i++) for(j=0;j<N;j++) RhoData[i][j] = curState[Rho][i][j];
    if(M%2==0) for(j=0;j<N;j++) RhoData[M][j] = 0.0;
    if(N%2==0) for(i=0;i<M;i++) RhoData[i][N] = 0.0;
    if(M%2==0 && N%2==0) RhoData[M][N] = 0.0;

    for(i=0;i<M;i++) RhoDataTop[i] = RhoTop[i];
    if(M%2==0) RhoDataTop[M] = 0.0;

    //Cylinder
    sum = cPos(0,DeltaR)*RhoTop[0] + cPos(0,DeltaR)*RhoDataTop[0] + cPos(M-1,DeltaR)*RhoDataTop[M-1] + cPos(M-1,DeltaR)*RhoDataTop[M-1]; // corners
    for(i=1;i<Mhere-1;i++)
    {
        if(i%2==0) sum = sum + 2.0*(cPos(i,DeltaR)*RhoDataTop[i] + cPos(i,DeltaR)*RhoDataTop[i]);   // top and bottom, even
        else sum = sum + 4.0*(cPos(i,DeltaR)*RhoDataTop[i] + cPos(i,DeltaR)*RhoDataTop[i]);   // top and bottom, odd
    }
    for(j=1;j<Nhere-1;j++)
    {
        if(j%2==0) sum = sum + 2.0*(cPos(0,DeltaR)*RhoDataTop[0] + cPos(Mhere-1,DeltaR)*RhoDataTop[Mhere-1]);   // left and right, even  (assumes RhoTop=0 at right)
        else sum = sum + 4.0*(cPos(0,DeltaR)*RhoDataTop[0] + cPos(Mhere-1,DeltaR)*RhoDataTop[Mhere-1]);   // left and right, even  (assumes RhoTop=0 at right)
    }
    for(i=1;i<Mhere-1;i++) for(j=1;j<Nhere-1;j++)
    {
        if(j%2==0)  // even rows, weights 8 and 4
        {
            if(i%2==0) sum = sum + 4.0*cPos(i,DeltaR)*RhoDataTop[i];
            else sum = sum + 8.0*cPos(i,DeltaR)*RhoDataTop[i];
        }
        else        // odd rows, weights 16 and 8
        {
            if(i%2==0) sum = sum + 8.0*cPos(i,DeltaR)*RhoDataTop[i];
            else sum = sum + 16.0*cPos(i,DeltaR)*RhoDataTop[i];
        }
    }

    CylMass = 4.0*PI*(DeltaR*DeltaZ/9.0)*sum;

    // Total Mass
    sum = cPos(0,DeltaR)*RhoData[0][0] + cPos(M-1,DeltaR)*RhoData[M-1][0] + cPos(0,DeltaR)*RhoData[0][N-1] + cPos(M-1,DeltaR)*RhoData[M-1][N-1]; // corners
    for(i=1;i<Mhere-1;i++)
    {
        if(i%2==0) sum = sum + 2.0*(cPos(i,DeltaR)*RhoData[i][0] + cPos(i,DeltaR)*RhoData[i][Nhere-1]);   // top and bottom, even
        else sum = sum + 4.0*(cPos(i,DeltaR)*RhoData[i][0] + cPos(i,DeltaR)*RhoData[i][Nhere-1]);   // top and bottom, odd
    }
    for(j=1;j<Nhere-1;j++)
    {
        if(j%2==0) sum = sum + 2.0*(cPos(0,DeltaR)*RhoData[0][j] + cPos(Mhere-1,DeltaR)*RhoData[Mhere-1][j]);   // left and right, even  (assumes RhoTop=0 at right)
        else sum = sum + 4.0*(cPos(0,DeltaR)*RhoData[0][j] + cPos(Mhere-1,DeltaR)*RhoData[Mhere-1][j]);   // left and right, even  (assumes RhoTop=0 at right)
    }
    for(i=1;i<Mhere-1;i++) for(j=1;j<Nhere-1;j++)
    {
        if(j%2==0)  // even rows, weights 8 and 4
        {
            if(i%2==0) sum = sum + 4.0*cPos(i,DeltaR)*RhoData[i][j];
            else sum = sum + 8.0*cPos(i,DeltaR)*RhoData[i][j];
        }
        else        // odd rows, weights 16 and 8
        {
            if(i%2==0) sum = sum + 8.0*cPos(i,DeltaR)*RhoData[i][j];
            else sum = sum + 16.0*cPos(i,DeltaR)*RhoData[i][j];
        }
    }

    TotMass = 4.0*PI*(DeltaR*DeltaZ/9.0)*sum;

};


void CalcMass()
{

    // Calculates the total mass of the core by solving the axisymmetric integral
    // 2*int( rho dr dtheta dz ) for r = 0 to rL, theta = 0 to 2pi, z = 0 to zL (hence the 2 out front)
    //
    // Can do this in a few ways

    // Crude approximation
    double CylMassC, TotMassC;
    CrudeIntegral(CylMassC,TotMassC);

    // Trapezoid approximation
    double CylMassT, TotMassT;
    TrapIntegral(CylMassT,TotMassT);

    // Simpson's approximation
    double CylMassS, TotMassS;
    SimpIntegral(CylMassS,TotMassS);

    cout << "Mass calculations (non-dim)" << endl;
    cout << "Cylinder masses (C,T,S):" << endl;
    cout << CylMassC << ", " << CylMassT << ", " << CylMassS << endl;
    cout << "Total masses (C,T,S):" << endl;
    cout << TotMassC << ", " << TotMassT << ", " << TotMassS << endl;
    cout << "Excess masses (C,T,S):" << endl;
    cout << (TotMassC-CylMassC) << ", " << (TotMassT-CylMassT) << ", " << (TotMassS-CylMassS) << endl;
    cout << endl;

    cout << "Mass calculations (dim; solar masses)" << endl;
    cout << "Cylinder masses (C,T,S):" << endl;
    cout << CylMassC/Sol2Code << ", " << CylMassT/Sol2Code << ", " << CylMassS/Sol2Code << endl;
    cout << "Total masses (C,T,S):" << endl;
    cout << TotMassC/Sol2Code << ", " << TotMassT/Sol2Code << ", " << TotMassS/Sol2Code << endl;
    cout << "Excess masses (C,T,S):" << endl;
    cout << (TotMassC-CylMassC)/Sol2Code << ", " << (TotMassT-CylMassT)/Sol2Code << ", " << (TotMassS-CylMassS)/Sol2Code << endl;



};
