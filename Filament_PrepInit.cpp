#include "Filament_PrepInit.H"

// Magnetized cylinder solution generator
void MagnetizedCylinder(double& Redge, double desiredLambda, int i);



void PrepareInitialState()
{
    // Calculate the baseline cylinder (what would be if mExcess = 0, and set units)
    CodeHeader("Baseline Magnetized Cylinder integrations");
    BaselineCylinder();

    // Given now that we know quantities of the background cylinder, we can evaluate the 'units'
    UnitConvert();

    // As of now, mExcess is actually the fraction of the baseline cylinder. Convert to non-dim mass
    mExcess = mExcess * (2.0*zL*lambda);

    // What is the total non-dim mass?
    cout << endl << "The total mass of the baseline cylinder is = " << 2.0*zL*lambda << endl;
    cout << "(" <<  2.0*zL*lambda/Sol2Code << " solar masses)" << endl;

    // Determine Vpot, Apot, RhoTop for reduced cylinder
    CodeHeader("Background Magnetized Cylinder integrations");
    InitCylinder();

    // What is the reduced non-dim mass?
    double desiredLambda = lambda - mExcess/(2.0*zL);
    cout << endl << "The total mass of the background cylinder is = " << 2.0*zL*desiredLambda << endl;
    cout << "(" <<  2.0*zL*desiredLambda/Sol2Code << " solar masses)" << endl;

    // Adds on Vpot from a chain of point masses
    CodeHeader("Determining Point Mass Potential");
    InitPoints();

    // We know the full initial potential. Get boundary conditions
    CodeHeader("Determining dVdR at the right edge");
    detDVDR();

    // Get boundary
    CodeHeader("Getting first filament boundary");
    getVbdy();

    // Adjust rho if mExcess > 0
    adjustRho();

    // With V and rho known everywhere, Q is easy
    CodeHeader("Initializing Q");
    InitQ();

    // Initial dQdPhi
    CodeHeader("Getting the first dQdPhi");
    CalcdQdP();

};


void BaselineCylinder()
{
    // Calculates the mExcess = 0 cylinder
    double Redge=0;

    MagnetizedCylinder(Redge,lambda,0);

    // This non-dim radius will be used to create our units.
    for(int i=0;i<N;i++) VContour[i] = Redge;

};


void InitCylinder()
{
    // Initialize Potentials from Cylinder
    // Given the value of lambda and mExcess, we can determine what the lambda value the reduced cylinder
    // must have.

    double desiredLambda = lambda - mExcess/(2.0*zL);
    if(desiredLambda < 0){ cout << endl << "ERROR! mExcess > totMass (" << mExcess << "," << lambda*2.0*zL << ")" << endl; exit(1);}

    double Redge=0; // don't actually need
    MagnetizedCylinder(Redge,desiredLambda,1); // 1 = populate curState

};


void InitPoints()
{

    // If there's an excess mass, impose the point potential
    if(mExcess<=0 || DEBUG==1) return;

    // Array for the point mass chain
    double PointV[M][N];
    double PointVold[M][N];

    int i,j;
    for(i=0;i<M;i++) for(j=0;j<N;j++) PointV[i][j] = 0.0;
    for(i=0;i<M;i++) for(j=0;j<N;j++) PointVold[i][j] = 0.0;

    // Loops, adding point masses to the chain. Breaks from additional points do not change potential
    int PPidx = 0; //1;// skips one particular because GM/r blows up for the idx=0 particle in the corner.
    double zP=0;//2.0*zL;

    // Mass in code units
    double Mcode = mExcess;
    cout << "MCode = " << Mcode << endl;



    while(true)
    {
        PPidx++;

        // Adds potential
        for(i=0;i<M;i++) for(j=0;j<N;j++)
        {
            // Distance locations of point particles
            double X1 = sqrt( pow(cPos(i,DeltaR),2) + pow(cPos(j,DeltaZ)-zP,2) );
            double X2 = sqrt( pow(cPos(i,DeltaR),2) + pow(cPos(j,DeltaZ)+zP,2) );

            if(zP==0)
            {
                // (GM/X blows up, use a softening parameter)
                double Xsoft = zL/2.0; //1.0*sqrt(pow(DeltaR,2)+pow(DeltaZ,2));
                X1 = X1 + Xsoft;
            }

            PointV[i][j] = PointV[i][j] + -Mcode/X1;
            if(zP > 0) PointV[i][j] = PointV[i][j] + -Mcode/X2;


        }


        // Calculates max error
        if(PPidx>2)
        {
            double err=0;

            for(i=0;i<M;i++) for(j=0;j<N;j++)
            {
                double locerr = fabs(PointV[i][j]-PointVold[i][j])/fabs(PointV[i][j]);
                if(locerr > err) err = locerr;
            }

            if(PPidx%1==0) cout << "Loop " << PPidx << " : Number of Points = " << 1+(PPidx-1)*2 << " (Max error = " << err << ")" << endl;

            // Is the max error small enough?
            if(err <= PointLoopTol) { PointLoopMax = PPidx; break; }
        }

        // Prepares for another iteration
        for(i=0;i<M;i++) for(j=0;j<N;j++) PointVold[i][j] = PointV[i][j];

        // New z location of point masses
        zP = zP + 2.0*zL;


        if(PPidx==PointLoopMax)
        {
            // Did not converge
            WaterlooHeader("PrepInit:InitPoints");
            cout << "Did not converge to a solution! (PointLoopMax = " << PointLoopMax << ")" << endl;
            exit(1);
        }

    }

    // If we're here, we have converged.
    // Normalize so that the upper-left corner has Vpot=0
    double Vnorm = PointV[0][N-1];
    for(i=0;i<M;i++) for(j=0;j<N;j++) PointV[i][j] = PointV[i][j] - Vnorm;

    // Add solution to the potential
    for(i=0;i<M;i++) for(j=0;j<N;j++) curState[Vpot][i][j] = curState[Vpot][i][j] + PointV[i][j];

};


void detDVDR()
{
    // We can calculate dV/dR on the right boundary, analytically from the cylinder + points

    // For points, we updated PointLoopMax to be the number of loops we used. We can use it again
    // to be consistent
    int i,j;
    int LoopIdx=0;
    double Mcode = mExcess;
    double zP = 0;

    // Just in case
    for(i=0;i<N;i++) VRight[i] = 0;

    // Potential grad from point masses
    for(LoopIdx=0;LoopIdx<=PointLoopMax;LoopIdx++)
    {
        if(DEBUG==2 || Mcode==0) break;

        for(j=0;j<N;j++)
        {
            double X1 = sqrt( pow(cPos(M-1,DeltaR),2) + pow(cPos(j,DeltaZ)-zP,2) );
            double X2 = sqrt( pow(cPos(M-1,DeltaR),2) + pow(cPos(j,DeltaZ)+zP,2) );

            VRight[j] = VRight[j] + Mcode/pow(X1,2);
            if(zP>0) VRight[j] = VRight[j] + Mcode/pow(X2,2);
        }
        zP = zP + 2.0*zL;
    }


    // The potential from the cylinder then adds the same value to each row
    // V has the form V = C1*ln(r/rEdge) + C2
    // so dV/dr = C1/r
    // C1 = Vdot(rEdge)*rEdge, C2 = V(rEdge)
    // This assumes the right boundary is at the filament boundary or outside the filament bdy

    // We're at r = rL (not necessarily rEdge = filament boundary)
    for(i=0;i<N;i++) VRight[i] = VRight[i] + Ccst[0]/rL;

    cout << "VRight = " << VRight[0];
    for(i=1;i<N;i++) cout << ", " << VRight[i];
    cout << endl;

};

void getVbdy()
{
    double Vrow[M];
    double rEdge = VContour[N-1];

    int i,j;

    // Vpot values we're interested in
    for(i=0;i<M;i++) Vrow[i] = curState[Vpot][i][N-1];

    //for(i=0;i<M;i++) cout <<  Vrow[i] << " ";
    //cout << endl;

    // Upate the global variable
    int lastIDX=0;
    if(rRatio>1) Vbdy = LIntY(Vrow,rEdge,DeltaR,M,lastIDX);
    else Vbdy = curState[Vpot][M-1][N-1];

    cout << "Vbdy is " << Vbdy << endl;

    // For each row, find where Vcontour is (N-1 -> skip the top)
    for(j=0;j<N-1;j++)
    {
        // Find first instance of two indices where Vval lie between
        // This should be fine if V(r,z) is monotonic in r for a given z
        int idx=0;
        for(i=1;i<M;i++)
        {
            idx = i;
            double Vleft = curState[Vpot][i-1][j];
            double Vright= curState[Vpot][i][j];
            //cout << "i=" << i << endl;

            if( (Vleft <= Vbdy && Vbdy <= Vright) || (Vleft >= Vbdy && Vbdy >= Vright) )
            {
                break;
            }
        }

        if(idx==0)
        {
            // never found one
            cout << "EEK! Never found a cell (z=" << j << ")" << endl;
        }

        double y1 = curState[Vpot][idx-1][j];
        double y2= curState[Vpot][idx][j];
        double x1 = cPos(i-1,DeltaR);
        double x2 = cPos(i,DeltaR);
        double m = (y2-y1)/(x2-x1);

        // contour location
        VContour[j] = x1 + (Vbdy-y1)/m;

    }

};


void adjustRho()
{
    int i,j;

    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        if( cPos(i,DeltaR) <= VContour[j] && curState[Rho][i][j]==0 ) curState[Rho][i][j]=1;
    }
};

void InitQ()
{
    // Q = rho * cs^2 exp(V/cs^2)

    // First, find Q value at boundary (where rho = 1 by definition)
    double Qbdy = 1.0 * exp(Vbdy);
    cout << "Init Qbdy = " << Qbdy << endl;

    // Loop over all cells, if left of boundary, evaluate Q, else assign value at bdy
    int i,j;
    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        if( cPos(i,DeltaR) <= VContour[j] )
        {
            curState[Q][i][j] = (curState[Rho][i][j])*exp(curState[Vpot][i][j]);
            //if(DEBUG==2) curState[Q][i][j] = curState[Q][i][j] + 1.0*cPos(i,DeltaR)*cPos(j,DeltaZ)/zL/rL;
        }
        else
        {
            curState[Q][i][j] = Qbdy;
        }
    }

    //for(i=0;i<M;i++) for(j=0;j<N;j++) cout << i << "," << j << " : " << curState[Q][i][j] << endl;


};



void UnitConvert()
{
    // VContour[i] is the non-dim radius of the filament
    double Rcyl = VContour[0];

    // Dimensional radius of filament (cgs)
    double Rdim = CYLINDERRADRAT*zL*Pc2Cm;

    // The pressure unit is Rcyl^2 * (c^4/R^2/G)
    Pbackground = pow(Rcyl,2)*( pow(Cbackground,4)/Gcst/pow(Rdim,2) );

    // Unit conversion, multiply by this to convert parsecs to code units
    double cgsCodeLength = Cbackground*Cbackground/sqrt(Gcst*Pbackground);  // cm/code length
    double pcCodeLength = cgsCodeLength/Pc2Cm; // pc/code length
    Pc2Code = 1/pcCodeLength; // code/pc length

    // Unit conversion, multiply by this to convert solar masses to code units
    double cgsCodeMass = pow(Cbackground,4)/sqrt(Pbackground*pow(Gcst,3)); // g/code
    double solCodeMass = cgsCodeMass/Sol2G; // sol/code
    Sol2Code = 1/solCodeMass;

    cout << "Unit conversions:" << endl;
    cout << "Parsec to Code = " << Pc2Code << endl;
    cout << "Solar mass to Code = " << Sol2Code << endl;

    // Non-dimensionalize the values
    zL = zL*Pc2Code;
    rL = rL*Pc2Code;

    // Define these values
    DeltaR = rL / (double(M)-1);
    DeltaZ = zL / (double(N)-1);

};
