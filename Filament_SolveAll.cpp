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

        // Solve Poisson
        SolvePoisson();

        // Solve Ampere
        //SolveAmpere();

        // Find the new filament boundary (Vbdy[])
        //NewVBdy();

        // Update Q
        //UpdateQ();

        // Update Rho
        //UpdateRho();

        // Converged?
        ConvergeTest(i,done); if(done || DOONE) break;

        // We never converged, sad face
        if(i+1==ConvergeLoopMax && DOONE==0)
        {
            cout << "We never converged message" << endl;
        }
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

    // Upate the global variable
    Vbdy = LIntY(Vrow,rEdge,DeltaR,M);
    cout << "Vbdy is now " << Vbdy << endl;

    // For each row, find where Vcontour (N-1 -> skip the top)
    for(j=0;j<N-1;j++)
    {
        // Find first instance of two indices where Vval lie between
        // This should be fine if V(r,z) is monotonic in r for a given z
        int idx=0;
        for(i=1;i<M;i++)
        {
            double Vleft = curState[Vpot][i-1][j];
            double Vright= curState[Vpot][i][j];
            //cout << "i=" << i << endl;

            if( (Vleft <= Vbdy && Vbdy <= Vright) || (Vleft >= Vbdy && Vbdy >= Vright) )
            {
                idx = i;
                break;
            }
        }

        if(idx==0)
        {
            // never found one
            cout << "Vbdy update problem! Never found a cell (z=" << j << ")" << endl;
        }

        double y1 = curState[Vpot][idx-1][j];
        double y2= curState[Vpot][idx][j];
        double x1 = cPos(i-1,DeltaR);
        double x2 = cPos(i,DeltaR);
        double m = (y2-y1)/(x2-x1);

        VContour[j] = x1 + (Vbdy-y1)/m;

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
        if(cPos(i,DeltaR)<=VContour[i]) QMap[i] = curState[Rho][i][N-1]*exp(curState[Vpot][i][N-1]);
        else QMap[i] = Qbdy;

    // Top row is done
    for(i=0;i<M;i++) curState[Q][i][N-1] = QMap[i];

    // This gives us a map. Now update the rest of the Q values.
    for(i=0;i<M;i++) for(j=0;j<N-1;j++) // skipping the top row
    {
        double localPhi = cPos(i,DeltaR)*curState[Apot][i][j];
        int idx;

        // Find the two indices localPhi lies between
        for(idx=1;idx<M;idx++)
        {
            double Pleft = cPos(i-1,DeltaR)*curState[Apot][i-1][j];
            double Pright= cPos(i,DeltaR)*curState[Apot][i][j];

            if( (Pleft <= localPhi && localPhi <= Pright) || (Pleft >= localPhi && localPhi >= Pright) )
            {
                break;
            }
        }

        if(idx==M-1)
        {
            // Possibly didn't find
        }

        // We now know which two Phi values the local Phi lies between. Intropolate a Q value

        double y1 = QMap[idx-1];
        double y2 = QMap[idx];
        double x1 = cPos(idx-1,DeltaR)*curState[Apot][idx-1][j];
        double x2 = cPos(idx,DeltaR)*curState[Apot][idx][j];
        double m = (y2-y1)/(x2-x1);

        curState[Q][i][j] = y1 + m*(localPhi-x1);

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

        cout << i << "," << j << " : " << locerrV << " " << curDV-preDV << " " << curState[Vpot][i][j] - prevState[Vpot][i][j] << endl;

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
