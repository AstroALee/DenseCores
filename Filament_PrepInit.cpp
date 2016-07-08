#include "Filament_PrepInit.H"

// Magnetized cylinder solution generator
void MagnetizedCylinder(double& Redge, double desiredLambda, int i);
void FastMagnetizedCylinder(double& Redge, double desiredLambda, int i);

void PrepareInitialState()
{
    // We will always have a baseline cylinder with lambda equal to the user input value
    CodeHeader("Cylinder Integration to get units");
    BaselineCylinder(0);
    // Given now that we know quantities of the background cylinder, we can evaluate the 'units'
    UnitConvert(); // Sets DeltaR, DeltaZ, zL and rL.
    BaselineCylinder(1); // Sets RhoTop, VContour_init
    BaselineCylinder(2); // Sets V and A and rho from cylinder

    // Now we have the cylinder determined.
    // Let's adjust density so the corner density is correct,
    // Also track total mass added for boundary condition.
    NewBall();

    // Adds on Vpot from a chain of point masses
    CodeHeader("Determining Point Mass Potential");
    double mPoints = TrapInt() - 2*zL*lambda;
    InitPoints(mPoints);

    // We know the full initial potential. Get boundary conditions
    CodeHeader("Determining dVdR at the right edge");
    detDVDR(mPoints);

    // With V and rho known everywhere, Q is easy
    CodeHeader("Initializing Q");
    //InitQ();
    UpdateQ(0);
    UpdateQ(1);

    // Initial dQdPhi
    CodeHeader("Getting the first dQdPhi");
    CalcdQdP();

    cout << "End of initial conditions" << endl;
}

void NewBall()
{
    double DeltaRho = rhoC - curState[Rho][0][0]; // Need to add this to the region around 0,0

    int i,j;
    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        double r = cPos(i,DeltaR), z = cPos(j,DeltaZ);
        double rad = pow(r*r+z*z,0.5);

        double ballRad = min(zL/1.1 , VContour[0]/1.1);

        if(rad <= ballRad )
        {
            double wgt = (1.0 - rad/ballRad);
            curState[Rho][i][j] = curState[Rho][i][j] + wgt*DeltaRho;
        }

    }

}

void PrepareInitialState_old()
{
    // Calculate the baseline cylinder (what would be if mExcess = 0)
    CodeHeader("Baseline Magnetized Cylinder integrations");
    BaselineCylinder(0);

    // Given now that we know quantities of the background cylinder, we can evaluate the 'units'
    UnitConvert(); // Sets DeltaR, DeltaZ, zL and rL.
    BaselineCylinder(1); // Sets RhoTop

    // What is the total non-dim mass?
    cout << endl << "The total mass of the baseline cylinder is = " << 2.0*zL*lambda << endl;
    cout << "(" <<  2.0*zL*lambda/Sol2Code << " solar masses)" << endl;

    // Get boundary
    CodeHeader("Getting first filament boundary");
    firstRbdy();

    // Rho is then the top values stretched out toward the boundary for each row.
    // Sets Rho for everywhere in the box, also sets value of Rbdy = Rho at contour
    CodeHeader("Stretching Density");
    StretchRho();

    // Now the boundary is a rho contour, but we've introduced a lot of extra mass in the process
    // We want the total mass to be equal to the mass of the lambda = desiredLambda cylinder.
    // We actually want to remove more so to set up a centrally concentrated density
    // We will remove mass so this central concentration has mass (1-contR)*Mass(Lambda=desiredLambda)
    // Doing so will mean the 'cylinder' will have a mass contR*Mass(Lambda=desiredLambda)
    double mCylGuess = contR*(2.0*zL*lambda);
    double mPoints   = (1-contR)*(2.0*zL*lambda);

    // We use the weighting function, with the goal of finding hte value of A such that the mass
    // enclosed is equal to this mCylGuess.
    if(contR<1)
    {
        CodeHeader("Finding and Applying Weighting Function for Density");
        double Awf = reduceCyl(mCylGuess);

        // Replace rho now with the altered density
        for(int i=0;i<M;i++) for(int j=0;j<N;j++)
        {
            double r = cPos(i,DeltaR); double z = cPos(j,DeltaZ);
            curState[Rho][i][j] = (1.0-WFunc(Awf,r,z))*curState[Rho][i][j]; // Rho = 0 outside the boundary, so who cares what WFunc returns
        }

        // Now we want to introduce a centrally concentrated mass equal to mPoints, so the total mass equals 2*zL*lambda
        // Do so by introducing a blob that smoothly connects to the cylinder
        CodeHeader("Adding ball");
        BallDensity(mPoints);
    }
    else
    {
        // Add a ball to get corner density to be rhoC
        double DeltaRho = rhoC - curState[Rho][0][0];
        double ballRad = VContour[0]/1.5;

        int i,j;
        for(i=0;i<M;i++) for(j=0;j<N;j++)
        {
            double r = cPos(i,DeltaR), z = cPos(j,DeltaZ);
            double curRad = pow(r*r+z*z,0.5);
            if(curRad <= ballRad) curState[Rho][i][j] = curState[Rho][i][j] + DeltaR;
        }

        double ballMass = DeltaRho*(4./3.)*PI*pow(ballRad,3);
        if(ballMass <= 2.0*zL*lambda)
        {
            mPoints = ballMass;
            mCylGuess = 2.0*zL*lambda - ballMass;
        }
        else
        {
            mPoints = zL*lambda;
            mCylGuess = zL*lambda;
        }

    }

    // With density and contour determined, we can determine the remaining quantities for integration

    // Let's get V and A of the cylinder
    CodeHeader("Getting reduced cylinder potentials");
    InitCylinder();

    // Adds on Vpot from a chain of point masses
    CodeHeader("Determining Point Mass Potential");
    InitPoints(mPoints);

    // We know the full initial potential. Get boundary conditions
    CodeHeader("Determining dVdR at the right edge");
    detDVDR(mPoints);

    // With V and rho known everywhere, Q is easy
    CodeHeader("Initializing Q");
    //InitQ();
    UpdateQ(0);
    UpdateQ(1);

    // Initial dQdPhi
    CodeHeader("Getting the first dQdPhi");
    CalcdQdP();

    cout << "End of initial conditions" << endl;

};

void BallDensity(double mPoints)
{
    double drag = pow(CYLINDERRADRAT,2);
    double eta = 1.0-contR;
    double ballRad = ((1/pow(drag,2) + pow(eta/drag/drag,2))/pow(1+pow(eta/drag,2),2))*VContour[0];
    cout << "Maximum radius for ball is " << ballRad/VContour[0] << " times the radius of the filament" << endl;
    for(int i=0;i<M;i++) for(int j=0;j<N;j++)
    {
        double r = cPos(i,DeltaR), z = cPos(j,DeltaZ);
        double rad = pow(r*r+z*z,0.5);
        double ballDen = 3*mPoints/PI/pow(ballRad,3); //15.0*mPoints/8/PI/pow(ballRad,3); // quadratic   3*mPoints/PI/pow(ballRad,3); // linear
        if(rad<=ballRad) curState[Rho][i][j] = ballDen*(1.0-pow(rad/ballRad,1)) + curState[Rho][i][j];
    }

}

void BaselineCylinder(int idx)
{
    // Calculates the mExcess = 0 cylinder
    // Finds where the calculated lambda = the command line value 'lambda'
    double Redge=0;

    cout << "Searching for cylinder with lambda = " << lambda << endl;

    FastMagnetizedCylinder(Redge,lambda,idx);

    // This non-dim radius will be used to create our units.
    for(int i=0;i<N;i++) VContour[i] = Redge;

};


void InitCylinder()
{
    // Initialize Potentials from Cylinder

    double Redge=0; // don't actually need
    FastMagnetizedCylinder(Redge,contR*lambda,2); // 2 = populate curState

};


void InitPoints(double mPoints)
{

    // If there's an excess mass, impose the point potential
    if(mPoints <=0 || DEBUG==1) return;

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
    double Mcode = mPoints;
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

            if(PPidx%10==0 || err <= PointLoopTol) cout << "Loop " << PPidx << " : Number of Points = " << 1+(PPidx-1)*2 << " (Max error = " << err << ")" << endl;

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


void detDVDR(double mPoints)
{
    // We can calculate dV/dR on the right boundary, analytically from the cylinder + points

    // For points, we updated PointLoopMax to be the number of loops we used. We can use it again
    // to be consistent
    int i,j;
    int LoopIdx=0;
    double Mcode = mPoints;
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

void getRbdy()
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
    if(rRatio>1) Rbdy = LIntY(Vrow,rEdge,DeltaR,M,lastIDX);
    else Rbdy = curState[Vpot][M-1][N-1];

    cout << "Rbdy is " << Rbdy << endl;

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

            if( (Vleft <= Rbdy && Rbdy <= Vright) || (Vleft >= Rbdy && Rbdy >= Vright) )
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
        VContour[j] = x1 + (Rbdy-y1)/m;

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
    // Q = P exp(V/cs^2) = rho*c^2 exp(V/c^2)

    // First, find Q value at boundary
    double Qbdy = Rbdy; // Rbdy is now the value of rho at the contour.  // 1.0 * exp(Rbdy);
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

    // The Rho unit is Rcyl^2 * (c^2/R^2/G)
    RhoTopCenter = pow(Rcyl,2)*( pow(Cbackground,2)/Gcst/pow(Rdim,2) );

    // Unit conversion, multiply by this to convert parsecs to code units
    double cgsCodeLength = Cbackground/sqrt(Gcst*RhoTopCenter);  // cm/code length
    double pcCodeLength = cgsCodeLength/Pc2Cm; // pc/code length
    Pc2Code = 1/pcCodeLength; // code/pc length

    // Unit conversion, multiply by this to convert solar masses to code units
    double cgsCodeMass = pow(Cbackground,3)/sqrt(RhoTopCenter*pow(Gcst,3)); // g/code
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

void firstRbdy()
{
    // Draws a straight line with the right ratio contR.

    double r2 = VContour[0]; // still set to rEdge for lambda = desiredLambda cylinder

    for(int j = 0; j<N; j++ )
    {
        double curZ = cPos(j,DeltaZ);
        VContour[j] = InitVB(curZ); //r2 * ( 1.0 - (1.0-contR)*(curZ/zL) );
    }

};

void StretchRho()
{
    //if(contR==1) return;

    // The top rho inside the arbitrary contour choice will be stretched at every row out to InitVB(z)
    double rEdgeT = InitVB(zL);

    int i,j;
    for(j=0;j<N;j++)
    {
        double rEdgeCur = InitVB(cPos(j,DeltaZ));

        for(i=0;i<M;i++)
        {
            if(cPos(i,DeltaR) > rEdgeCur)
            {
                curState[Rho][i][j] = 0;
            }
            else
            {
                double adjRat = rEdgeT/rEdgeCur;
                double adjPos = adjRat*cPos(i,DeltaR);
                //cout << adjRat << endl;

                // This adjusted position is used to interpolate
                int lastIDX; // not used
                curState[Rho][i][j] = LIntY(RhoTop, adjPos, DeltaR, M, lastIDX);
            }
        }

    }

    // The value of rho on this contour needs to be saved too
    int idx; // not used
    Rbdy = LIntY(RhoTop, rEdgeT, DeltaR, M, idx);


};

double reduceCyl(double mCylGuess) // want to find A such that mass = mCylGuess
{
    // want to find root of function f() = curMass - mCylGuess as a function of A
    double Awf=0.5;

    // Newton Raphson to get the right mass
    double err = 1;
    double NRidx=0;

    double oldA = Awf;
    double oldM = MassInt(Awf);
    cout << "Loop " << 0 << ": Using " << Awf << " the reduced mass is " << oldM << " (want " << mCylGuess << "; " << err << ")" << endl;
    double curA = 1.1*Awf;
    double curM = MassInt(curA);
    cout << "Loop " << 0 << ": Using " << curA << " the reduced mass is " << curM << " (want " << mCylGuess << "; " << err << ")" << endl;

    while(true)
    {
        NRidx++;

        double derv = (curM-oldM)/(curA-oldA); // the mCylGuesses cancel out

        // current becomes old
        oldA = curA;
        oldM = curM;
        // new guess
        curA = oldA - (oldM-mCylGuess)/derv;
        // new mass
        curM = MassInt(curA);

        err = fabs(curM-mCylGuess)/mCylGuess;
        cout << "Loop " << NRidx << ": Using " << curA << " the reduced mass is " << curM << " (want " << mCylGuess << "; " << err << ")" << endl;

        if( err < CylTol) break;
        if(NRidx==LoopMAX) {cout << "Did not find cylinder fit!" << endl; break;}
    }


    return curA;
};

double WRho(int i, int j, double A)
{
    double r = cPos(i,DeltaR);
    double z = cPos(j,DeltaZ);
    double wr = (1.0-WFunc(A,r,z))*curState[Rho][i][j];
    return wr;
};


double MassInt(double A)
{
    // Approximates the integral int(f dxdy) as 0.25*DeltaX*DeltaY*Sum(w*f(xi,yi))
    // where w = 1 for corners, 2 for edges, 4 for interior

    int i,j; double sum=0;


    //Total Mass
    sum = cPos(0,DeltaR)*WRho(0,0,A) + cPos(M-1,DeltaR)*WRho(M-1,0,A) + cPos(0,DeltaR)*WRho(0,N-1,A) + cPos(M-1,DeltaR)*WRho(M-1,N-1,A); // corners
    for(i=1;i<M-1;i++) sum = sum + 2.0*(cPos(i,DeltaR)*WRho(i,0,A) + cPos(i,DeltaR)*WRho(i,N-1,A)); // top and bottom edges
    for(j=1;j<N-1;j++) sum = sum + 2.0*(cPos(0,DeltaR)*WRho(0,j,A) + cPos(M-1,DeltaR)*WRho(M-1,j,A)); // left and right edges
    for(i=1;i<M-1;i++) for(j=1;j<N-1;j++) sum = sum + 4.0*(cPos(i,DeltaR)*WRho(i,j,A)); // interiors
    double TotMass = 4.0*PI*(DeltaR*DeltaZ/4.0)*sum;

    return TotMass;
};
