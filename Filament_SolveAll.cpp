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
        string NumStr = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
        PrintState(curState,"_",NumStr);

        cout << endl << "Beginning iteration " << i << endl;

        // Solve Poisson -- solution stored in newState
        SolvePoisson();

        // Solve Ampere -- solution stored in newState
        SolveAmpere();

        // About to use new solution to relax, make data in curState the prevState
        CopyState(curState,prevState);

        // Relax V and A using newState and prevState, use this to recalc boundary, Q, and Rho
        RelaxSoln();

        // Find the new filament boundary (Vbdy[])
        NewVBdy();

        // Update Q
        UpdateQ();

        // Update Rho
        UpdateRho();

        // Update dQdPhi
        CalcdQdP();

        // Converged?
        ConvergeTest(i,done);
        if(done) { cout << "Converged! Breaking out." << endl; break; }

    }

};

// Relax your newState
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
    double rEdge = VContour[N-1]; //RADIALRATIO*zL;

    int i,j;

    // Vpot values we're interested in
    for(i=0;i<M;i++) Vrow[i] = curState[Vpot][i][N-1];

    // Upate the global variable (boundary is always at rEdge on the top row)
    int lastIDX = 0;
    if(rRatio>1) Vbdy = LIntY(Vrow,rEdge,DeltaR,M,lastIDX);
    else Vbdy = curState[Vpot][M-1][N-1];

    cout << "Vbdy is now " << Vbdy << " at index " << lastIDX << endl;

    // For each row, find where Vcontour (N-2 <-> skipping the top)
    // Do so by searching for indices near the last row's index (to avoid jumps for double-values)
    for(j=N-2;j>=0;j--)
    {
        // Index where this cell and the cell to the left has the desired Vpot value
        int idx=-1;

        // Find first instance where Vval occurs
        // Does so by looking left and right from lastIDX, the location of the contour for
        // the row above this one.
        for(i=0;i<M;i++)
        {
            // can we look left?
            if(lastIDX-i-1 >= 0)
            {
                double Vleft = curState[Vpot][lastIDX-i-1][j];
                double Vright = curState[Vpot][lastIDX-i][j];
                if( (Vleft<=Vbdy && Vbdy<=Vright) || (Vleft>=Vbdy && Vbdy>=Vright) ) { idx = lastIDX-i; lastIDX = idx; break; }
            }

            // else look right?
            if(lastIDX+i+1 <= M-1)
            {
                double Vright = curState[Vpot][lastIDX+i+1][j];
                double Vleft = curState[Vpot][lastIDX+i][j];
                if( (Vleft<=Vbdy && Vbdy<=Vright) || (Vleft>=Vbdy && Vbdy>=Vright) ) { idx = lastIDX+i+1; lastIDX = idx; break; }
            }
        }

        if(idx==-1)
        {
            cout << "NEVER FOUND Vbdy in row " << j << endl;
            VContour[j] = VContour[j+1]; // manually set equal to value above (but something likely went wrong)
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
            //cout << "Contour value at " << j << " is " << VContour[j] << endl;
        }
    }

    cout << "New array of Contour radii is " << VContour[0];
    for(i=0;i<N;i++) cout  <<", " << VContour[i];
    cout << endl;

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

    // Test if this map is monotonic or not
    cout << "Updated Q at top: " ;
    for(i=0;i<M-1;i++) cout << QMap[i] << ", ";
    cout << QMap[M-1] << endl;

    // This gives us a map. Now update the rest of the Q values.
    for(i=0;i<M;i++) for(j=0;j<N-1;j++) // skipping the top row
    {
        // Local phi value we want to use with PhiMap and QMap
        double localPhi = cPos(i,DeltaR)*curState[Apot][i][j];
        int idx;

        // Find the two indices localPhi lies between (if not monotonic, could be weird)
        for(idx=1;idx<M;idx++)
        {
            double Pleft = PhiMap[idx-1]; // cPos(idx-1,DeltaR)*curState[Apot][idx-1][j];
            double Pright= PhiMap[idx];   // cPos(idx,DeltaR)*curState[Apot][idx][j];

            if( (Pleft <= localPhi && localPhi <= Pright) || (Pleft >= localPhi && localPhi >= Pright) )
            {
                break;
            }
        }

        if(idx==M-1)
        {
            // Possibly didn't find, but we know Phi analytically
            // Assumes we are beyond the boundary
        }

        // We now know which two Phi values the local Phi lies between. Intropolate a Q value

        double y1 = QMap[idx-1];
        double y2 = QMap[idx];
        double x1 = PhiMap[idx-1]; //cPos(idx-1,DeltaR)*curState[Apot][idx-1][j];
        double x2 = PhiMap[idx];   //cPos(idx,DeltaR)*curState[Apot][idx][j];
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
    int Vi,Vj,Ai,Aj;
    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        double locerrV=0,locerrA=0;

        // Since the actual value of V is not meaningful, we'll calculate error by esimating
        // dV = (dV/dr)dr + (dV/dz)*dz
        // using a 2nd order center difference
        // dV = ( V(i+1,j) - V(i-1,j) + V(i,j+1) - V(i,j-1) ) / 2

        double curDV=0,preDV=0;
        double curDA=0,preDA=0;

        if(i>1 && i<M-1)
        {
            curDV = 0.5*(curState[Vpot][i+1][j] + curState[Vpot][i-1][j] + curState[Vpot][i][j+1] - curState[Vpot][i][j-1]);
            preDV = 0.5*(prevState[Vpot][i+1][j] + prevState[Vpot][i-1][j] + prevState[Vpot][i][j+1] - prevState[Vpot][i][j-1]);
            if(preDV!=0) locerrV = fabs(curDV-preDV)/fabs(preDV);

            curDA = 0.5*(curState[Apot][i+1][j] + curState[Apot][i-1][j] + curState[Apot][i][j+1] - curState[Apot][i][j-1]);
            preDA = 0.5*(prevState[Apot][i+1][j] + prevState[Apot][i-1][j] + prevState[Apot][i][j+1] - prevState[Apot][i][j-1]);
            if(preDA!=0) locerrA = fabs(curDA-preDA)/fabs(preDA);
        }


        //cout << i << "," << j << " : " << locerrV << " " << curDV-preDV << " " << curState[Vpot][i][j] - prevState[Vpot][i][j] << " " << curState[Apot][i][j] - prevState[Apot][i][j] << endl;

        // Tabulates the max error
        if(locerrV > errV)
        {
            errV = locerrV;
            Vi = i;
            Vj = j;
        }
        if(locerrA > errA)
        {
            errA = locerrA;
            Ai = i;
            Aj = j;
        }
    }

    // Location of maximum error.
    cout << "The largest error for V occurs at (i,j)=(" << Vi << "," << Vj << ")" << endl;
    cout << "The largest error for A occurs at (i,j)=(" << Ai << "," << Aj << ")" << endl;

    // Print out
    PrintError(loopnum,errV,errA);

    // Can we go home yet?
    if(errV < ConvergeTol && errA < ConvergeTol) done = 1;

    string a;
    if(done) a = "SUCCESS!"; else a = " ";
    cout << "Converge iteration " << loopnum << " : " << a << " Max errors V,A = " << errV << " , " << errA << endl;

};

void PrintError(int i, double v, double a)
{
    ofstream myfile;
    if(i==0) myfile.open("Errors.out",ofstream::out);
    else myfile.open("Errors.out",ofstream::app);

    myfile << "Iteration " << i << " " << v << " " << a << endl;

    myfile.close();

};
