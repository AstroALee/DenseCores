#include "Filament_SolveAll.H"

void SolvePoisson();
void SolveAmpere();


void Converge()
{
    int done = 0;

    // Loops
    int i;
    for(i=0;i<ConvergeLoopMax;i++)
    {
        cout << "Beginning iteration " << i << endl;

        // Solve Poisson
        SolvePoisson();

        // Solve Ampere
        SolveAmpere();

        // Relax V and A, use this to recalc boundary, Q, and Rho
        RelaxSoln();

        // Find the new filament boundary (Vbdy[])
        NewVBdy();

        // Update Q
        UpdateQ();

        // Update Rho
        UpdateRho();

        // Converged?
        ConvergeTest(i,done);
        if(done) { cout << "Converged! Breaking out." << endl; break; }
        
    }

};

// Relax
void RelaxSoln()
{
    int i,j;

    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        //cout << curState[Vpot][i][j] << " " << newState[Vpot][i][j] << endl;
        curState[Vpot][i][j] = relaxFrac*newState[Vpot][i][j] + (1.0-relaxFrac)*curState[Vpot][i][j];
        curState[Apot][i][j] = relaxFrac*newState[Apot][i][j] + (1.0-relaxFrac)*curState[Apot][i][j];
    }


};

// Find the new filament boundary
void NewVBdy()
{
    // The location at the top of the box will never change
    // Let's find the V value there

    double Vrow[M];
    double rEdge = RADIALRATIO*zL;

    int i,j;

    // Vpot value we're interested in
    for(i=0;i<M;i++) Vrow[i] = curState[Vpot][i][N-1];

    // Upate the global variable (boundary is always at rEdge on the top row)
    int lastIDX = 0;
    Vbdy = LIntY(Vrow,rEdge,DeltaR,M,lastIDX);
    cout << "Vbdy is now " << Vbdy << endl;

    // For each row, find where Vcontour (N-2 -> skip the top)
    // Do so by searching for indices near the last row's index (to avoid jumps for double-values)
    for(j=N-2;j<0;j--)
    {
        int idx=-1;

        // Find first instance of two instances where Vval falls between
        for(i=0;i<M;i++)
        {
            if(lastIDX-i-1 >= 0) // look left
            {
                double VleftVal = curState[Vpot][lastIDX-i-1][j];
                double VrightVal = curState[Vpot][lastIDX-i][j];

                if( (VleftVal <= Vbdy && Vbdy <= VrightVal) || (VleftVal >= Vbdy && Vbdy >= VrightVal) )
                {
                    idx = lastIDX-i; // right coordinate
                    break;
                }
            }
            else if(lastIDX+i+1 <= M-1) // look right
            {
                double VleftVal = curState[Vpot][lastIDX+i][j];
                double VrightVal = curState[Vpot][lastIDX+i+1][j];

                if( (VleftVal <= Vbdy && Vbdy <= VrightVal) || (VleftVal >= Vbdy && Vbdy >= VrightVal) )
                {
                    idx = lastIDX+i+1; // right coordinate
                    break;
                }
            }
        }


        if(idx==-1)
        {
            // Didn't find anything
            cout << "DID NOT FIND VALUE" << endl;
            VContour[j] = VContour[j+1];
        }
        else
        {
            // Interpolate
            double y1 = curState[Vpot][idx-1][j];
            double y2= curState[Vpot][idx][j];
            double x1 = cPos(idx-1,DeltaR);
            double x2 = cPos(idx,DeltaR);
            double m = (y2-y1)/(x2-x1);

            VContour[j] = x1 + (Vbdy-y1)/m;

            // Update lastIDX
            lastIDX = idx;
        }

    }

};

// Given a new V and A, update Q
// This is done differently than how we initialized Q in PrepInit
// Q is constant along field lines, so we use the tow row (whose density never changes) to create
// a map between the updated Phi = r*A and Q = rho*exp(V) values.
// Then everywhere else in the box, we update Q appropriately.
void UpdateQ()
{
    int i,j;

    double rEdge = VContour[N-1];

    // Make map
    double PhiMap[M];
    double QMap[M];
    for(i=0;i<M;i++) PhiMap[i] = cPos(i,DeltaR)*curState[Apot][i][N-1];


    // To get the new boundary value of Q, use the potential and the assumption that rho = 1
    double Qbdy = 1.0 * exp(Vbdy);

    // Since rho does not change on the upper row, we can update Q using its definition
    for(i=0;i<M;i++)
        if(cPos(i,DeltaR)<=VContour[N-1]) QMap[i] = RhoTop[i]*exp(curState[Vpot][i][N-1]);
        else QMap[i] = Qbdy;

    // Top row is done
    for(i=0;i<M;i++) curState[Q][i][N-1] = QMap[i];

    // This gives us a map. Now update the rest of the Q values.
    for(i=0;i<M;i++) for(j=0;j<N-1;j++) // skipping the top row
    {
        double localPhi = cPos(i,DeltaR)*curState[Apot][i][j];
        int idx;

        if(cPos(i,DeltaR)>VContour[j])
        {
            curState[Q][i][j] = Qbdy;
        }
        else
        {
            // Find the two indices localPhi lies between
            for(idx=1;idx<M;idx++)
            {
                double Pleft = cPos(idx-1,DeltaR)*curState[Apot][idx-1][j];
                double Pright= cPos(idx,DeltaR)*curState[Apot][idx][j];

                if( (Pleft <= localPhi && localPhi <= Pright) || (Pleft >= localPhi && localPhi >= Pright) )
                {
                    break;
                }
            }

            if(idx==M-1)
            {
                // Possibly didn't find, but we know Phi analytically?
            }

            // We now know which two Phi values the local Phi lies between. Intropolate a Q value

            double y1 = QMap[idx-1];
            double y2 = QMap[idx];
            double x1 = cPos(idx-1,DeltaR)*curState[Apot][idx-1][j];
            double x2 = cPos(idx,DeltaR)*curState[Apot][idx][j];
            double m = (y2-y1)/(x2-x1);

            curState[Q][i][j] = y1 + m*(localPhi-x1);
        }
    }

};

// Given a new Q and V, calculate the new rho
void UpdateRho()
{
    int i,j;
    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        if( cPos(i,DeltaR) > VContour[j] ) curState[Rho][i][j] = 0.0;
        else curState[Rho][i][j] = curState[Q][i][j]*exp(-curState[Vpot][i][j]);
    }

};


// Compares the previous state to the current updated state to assess if we have converged to a SS
void ConvergeTest(int loopnum, int& done)
{
    double errV=0, errA=0;

    int i,j;
    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        double locerrV=0,locerrA=0;

        // Since the actual value of V is not meaningful, we'll calculate error by esimating
        // dV = (dV/dr)dr + (dV/dz)*dz
        // using a 2nd order center difference
        // dV = ( V(i+1,j) - V(i-1,j) + V(i,j+1) - V(i,j-1) ) / 2

        double curDV=0,preDV=0;
        if(i>0 && i<M-1)
        {
            curDV = 0.5*(curState[Vpot][i+1][j] + curState[Vpot][i-1][j] + curState[Vpot][i][j+1] - curState[Vpot][i][j-1]);
            preDV = 0.5*(prevState[Vpot][i+1][j] + prevState[Vpot][i-1][j] + prevState[Vpot][i][j+1] - prevState[Vpot][i][j-1]);

            if(curDV!=0) locerrV = fabs( curDV - preDV ) / fabs( curDV );

        }

        if(curState[Apot][i][j] != 0) {
                locerrA = fabs( curState[Apot][i][j] - prevState[Apot][i][j] )
                        / fabs(curState[Apot][i][j]); }

        //cout << i << "," << j << " : " << locerrV << " " << curDV-preDV << " " << curState[Vpot][i][j] - prevState[Vpot][i][j] << " " << curState[Apot][i][j] - prevState[Apot][i][j] << endl;

        // Tabulates the max error
        if(locerrV > errV) errV = locerrV;
        if(locerrA > errA) errA = locerrA;
    }

    // Can we go home yet?
    if(errV < ConvergeTol && errA < ConvergeTol) done = 1;

    string a;
    if(done) a = "SUCCESS!"; else a = " ";
    cout << "Converge iteration " << loopnum << ": " << a << " Max errors V,A = " << errV << "," << errA << endl;

};
