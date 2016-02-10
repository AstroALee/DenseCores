
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
int ConvergeLoopMax = 100000;
double ConvergeTol = 1e-4;



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
    if(a-1 != 7)
    {
        WaterlooHeader("CheckArguments");
        cout << "You didn't provide the right number of arguments. You need seven." << endl;
        cout << "(M,N,z-length (pc), Ratio r-length to z-length,Non-Dim Mass excess, beta_inf, n)" << endl;
        cout << "You provided " << a-1 << " argument" << Plural(a-1,"s") << "." << endl;
        exit(1);
    }

};

// Read in the arguments and store them in the appropriate global variables.
void ReadArguments(char *a[])
{
    M = atoi(a[1]);
    N = atoi(a[2]);
    zL = atof(a[3]);
    rRatio = atof(a[4]);
    mExcess = atof(a[5]);
    betaInf = atof(a[6]);
    nCyl = atof(a[7]);

    // Check that some numbers are reasonable.
    if(M < 3 || N < 3)
    {
        WaterlooHeader("ReadArguments");
        cout << "The number of grid cells is weird. You provided (M,N) = (" << M << "," << N << ")" << endl;
        cout << "You need at least 2 in each dimension" << endl;
        exit(1);
    }

    if(zL < 0 || rRatio < 0 || mExcess < 0)
    {
        WaterlooHeader("ReadArguments");
        cout << "Some values of zL, rRatio, or mExcess are weird. You provided zL, rRatio, mExcess = " << zL << ", " << rRatio << ", " << mExcess << endl;
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

void PrintArguments()
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
    cout << "DOONE = " << DOONE << endl;
    cout << "UNIFORM = " << UNIFORM << endl;
    cout << "RADIALRATIO = " << RADIALRATIO << endl;
};

void CalcDerived()
{
    // Unit conversion, multiply by this to convert parsecs to code units
    double cgsCodeLength = Cbackground*Cbackground/sqrt(Gcst*Pbackground);  // cm/code length
    double pcCodeLength = cgsCodeLength/Pc2Cm; // pc/code length
    Pc2Code = 1/pcCodeLength; // code/pc length

    // Unit conversion, multiply by this to convert solar masses to code units
    double cgsCodeMass = pow(Cbackground,4)/sqrt(Pbackground*pow(Gcst,3));
    double solCodeMass = cgsCodeMass/Sol2G;
    Sol2Code = 1/solCodeMass;

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

    // Are we only doing one converge loop?
    if(DOONE) ConvergeLoopMax = 1;
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
