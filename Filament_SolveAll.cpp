#include "Filament_SolveAll.H"

void SolvePoisson();
void SolveAmpere();
//void UpdateQ(int type);

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

        // Update Q
        UpdateQ(0);

        // Update Rho
        int success=0;
        //success = UpdateRho(0);
        success = 1;

        if(success)
        {
            cout << "Found a rho contour!" << endl;

            // Update Q (also updates global value of Rbdy)
            UpdateQ(1);

            // Update Rho right of contour
            success=UpdateRho(1);

            // Update dQdPhi
            CalcdQdP();

            // Update boundary conditions
            UpdateBC(0);

            // Converged?
            ConvergeTest(i,done);
            if(done) { cout << "Converged! Breaking out." << endl; break; }
        }
        else
        {
            cout << "Did not find a rho contour, reverting back..." << endl;

            // Copy previous state back to current state
            CopyState(prevState,curState);

            // Adjust the boundary conditions
            UpdateBC(1);

        }


    }
    UpdateQ(0);

};

// Update Boundary conditions
void UpdateBC(int type)
{
    // type = 0, recalc dVdr
    // type = 1, reverting back, upping Malter
    if(type==0)
    {
        cout << "Updating boundary conditions!" << endl;

        // Need cylinder with radius equal to VContour[N-1]
        double Crad = VContour[N-1];
        // Find last index
        int i=0; while(true) {i++; if(cPos(i,DeltaR)>Crad) break; }
        if(i>M-1) i=M-1;

        // 1D trap integral to get reduced cylinder mass
        double Mredcy = OneDMassInt(RhoTop, i);
        cout << "Cylinder mass = " << Mredcy << " (2zlambda = " << 2*zL*lambda << ")" << endl;

        // Total mass = totMass
        if(Mredcy>2*zL*lambda)
        {
            cout << "ISSUE!: Reduced cylinder has more mass than total (2*zL*lambda) mass! WHAT?!" << endl;
            Mredcy = totMass;
        }

        // compute 2d mass integral to get mass in the box
        double Mbox = TrapInt();
        cout << "Mbox = " << Mbox << endl;

        if(Mredcy>Mbox)
        {
            cout << "ISSUE!: Reduced cylinder has more mass than mass in the box... weird as f*@$ !" << endl;
            exit(1);
        }


        if(Mredcy>2*zL*lambda) Mredcy = 2*zL*lambda;

        double Mexcess = 2*zL*lambda - Mredcy; //Mbox - Mredcy;
        if(Mexcess<0)
        {
            cout << "Mexcess < 0!" << endl;
            Mexcess = 0;
        }
        double Medgecy = Mredcy; //2*zL*lambda - Mexcess;
        if(Medgecy<0)
        {
            cout << "Medgecy < 0!" << endl;
            Medgecy = 0;
            Mexcess  = 2*zL*lambda;
        }

        // Over-rides
        Medgecy = 0.50*2*zL*lambda;
        Mexcess = 0.50*2*zL*lambda;

        // Boundary conditions will use Mexcess for the points and totMass - Mexcess for the cylinder
        double CyldVdr = Medgecy/zL/rL;
        for(i=0;i<N;i++) VRight[i]=CyldVdr; // overwrites

        // Now adds in point masses dVdr = G*Mexcess/Rad^2 (same number of points as in initial conditions)
        // Potential grad from point masses
        double zP=0;
        cout << "Max points = " << PointLoopMax << endl;
        for(i=0;i<=PointLoopMax;i++)
        {
            //if(DEBUG==2 || Mcode==0) break;
            int j;
            for(j=0;j<N;j++)
            {
                double X1 = sqrt( pow(cPos(M-1,DeltaR),2) + pow(cPos(j,DeltaZ)-zP,2) );
                double X2 = sqrt( pow(cPos(M-1,DeltaR),2) + pow(cPos(j,DeltaZ)+zP,2) );

                VRight[j] = VRight[j] + Mexcess/pow(X1,2);
                if(zP>0) VRight[j] = VRight[j] + Mexcess/pow(X2,2);
            }
            zP = zP + 2.0*zL;
        }

        // dVdR is now updated
        cout << "VRight values = ";
        for(int j=0;j<N;j++) cout << VRight[j] << ", ";
        cout << endl;

    }
    else
    {
        cout << "Reverting boundary conditions!" << endl;

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

// Find the Rho contour
void UpdateBdy()
{
    // We want to find a contour where the

}

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
    if(rRatio>1) Rbdy = LIntY(Vrow,rEdge,DeltaR,M,lastIDX);
    else Rbdy = curState[Vpot][M-1][N-1];

    cout << "Rbdy is now " << Rbdy << " at index " << lastIDX << endl;

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
                if( (Vleft<=Rbdy && Rbdy<=Vright) || (Vleft>=Rbdy && Rbdy>=Vright) ) { idx = lastIDX-i; lastIDX = idx; break; }
            }

            // else look right?
            if(lastIDX+i+1 <= M-1)
            {
                double Vright = curState[Vpot][lastIDX+i+1][j];
                double Vleft = curState[Vpot][lastIDX+i][j];
                if( (Vleft<=Rbdy && Rbdy<=Vright) || (Vleft>=Rbdy && Rbdy>=Vright) ) { idx = lastIDX+i+1; lastIDX = idx; break; }
            }
        }

        if(idx==-1)
        {
            cout << "NEVER FOUND Rbdy in row " << j << endl;
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

            VContour[j] = x1 + (Rbdy-y1)/m;
            //cout << "Contour value at " << j << " is " << VContour[j] << endl;
        }
    }

    cout << "New array of Contour radii is " << VContour[0];
    for(i=0;i<N;i++) cout  <<", " << VContour[i];
    cout << endl;

};



// Finds the ratio of a given contour. Input the upper-row's starting radius.
double FindContRatio(double topRad, int type)
{
    double curCont;

    int i=0,j=0;

    // We know the top radius. Get the density at that location
    while(true){i++; if(cPos(i,DeltaR)>=topRad) break;}

    // Interpolate to get density
    double Dl = RhoTop[i-1]; double Dr = RhoTop[i];
    double x = (topRad - cPos(i-1,DeltaR))/(DeltaR);
    double desiredRho = (1-x)*Dl + x*Dr;
    //cout << "Desired rho = " << desiredRho << endl;

    if(type)
    {
        VContour[N-1] = topRad;
    }

    // Now loop over every row and find the location of this density value
    int iLast = i;

    for(j=N-2;j>=0;j--)
    {
        //cout << "Trying for j = " << j << endl;

        int iLeft = iLast; int iRight = iLast;
        curCont = 0;

        while(true)
        {
            // Search around the last location of the contour and span out from there.

            // look left
            double Dl = curState[Rho][iLeft-1][j];
            double Dr = curState[Rho][iLeft][j];
            int Check1 = ( Dl<=desiredRho && desiredRho<=Dr );
            int Check2 = ( Dl>=desiredRho && desiredRho>=Dr );
            if(Check1 || Check2)
            {
                iLast = iLeft;
                double m = DeltaR/(Dr-Dl);
                curCont = cPos(iLeft-1,DeltaR) + m*(desiredRho-Dl);
                //cout << "Contour location for j=" << j << " is " << curCont << endl;
                if(type) VContour[j] = curCont;
                break;
            }

            // else look right
            Dl = curState[Rho][iRight-1][j];
            Dr = curState[Rho][iRight][j];
            Check1 = ( Dl<=desiredRho && desiredRho<=Dr );
            Check2 = ( Dl>=desiredRho && desiredRho>=Dr );

            if( Check1 || Check2 )
            {
                iLast = iRight;
                double m = DeltaR/(Dr-Dl);
                curCont = cPos(iRight-1,DeltaR) + m*(desiredRho-Dl);
                //cout << "Contour location for j=" << j << " is " << curCont << endl;
                if(type) VContour[j] = curCont;
                break;
            }

            // If here, didn't find it
            if(iLeft==1 && iRight==M-1)
            {
                cout << "Did not find the density value!" << endl;
                exit(1);
            }

            // Prepares to start again
            if(iLeft > 1) iLeft--;
            if(iRight < M-2) iRight++;
        }
    }

    double rat = 0;
    if(curCont!=0) rat = topRad/curCont;
    //cout << "Contour ratio is = " << rat << endl;

    return rat;
}


// Given a new Q and V, calculate the new rho
int UpdateRho(int type)
{
    int success=0;
    int i,j;

    if(type)
    {
        success=1;

        // We know where the contour is
        for(i=0;i<M;i++)for(j=0;j<N;j++)
        {
            if(cPos(i,DeltaR)>VContour[j]) curState[Rho][i][j] = 0;
        }

        return success;
    }


    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        if(j==N-1) curState[Rho][i][j] = RhoTop[i]; // top row is always the cylinder
        else curState[Rho][i][j] = curState[Q][i][j]*exp(-curState[Vpot][i][j]); // everywhere else uses the Q update
    }

    if(contR==1) return 1;

    // Now find the contour that has ratio rCont
    // Need to find the two values of Rho at the top the contour resides between

    // which two indices does the current boundary lie between? Gives first guess for where to start
    i=0; while(true) {i++; if( cPos(i,DeltaR) >= VContour[N-1] ) break;}
    int iLeft = i; int iRight = i;

    while(true)
    {
        // look left
        double contLe = FindContRatio(cPos(iLeft-1,DeltaR),0);
        double contRi = FindContRatio(cPos(iLeft,DeltaR),0);

        //cout << "Left(" << iLeft << ") Left and Right Contour Ratios are " << contLe << " and " << contRi << endl;

        // Was it there?
        if( (contLe <= contR && contR <= contRi) || (contLe >= contR && contR >= contRi) )
        {
            // Linearly interpolate to find the exact radius
            double m = DeltaR/(contRi-contLe);
            double contRadius = cPos(iLeft-1,DeltaR) + m*(contR-contLe); // FindContRatio( contLe + m*(contR-contRi) , 1 ); // also override VContour with new values
            double contFinal = FindContRatio(contRadius,1); // 1 = Overwrites VContour
            cout << "Found a contour with radius " << contRadius << " with ratio " << contFinal << " (wanted contR= " << contR << ")" << endl;
            success=1;
            break;
        }

        // else look right
        contLe = FindContRatio(cPos(iRight-1,DeltaR),0);
        contRi = FindContRatio(cPos(iRight,DeltaR),0);
        //cout << "Right(" << iRight << "): Left and Right Contour Ratios are " << contLe << " and " << contRi << endl;

        // Was it there?
        if( (contLe <= contR && contR <= contRi) || (contLe >= contR && contR >= contRi) )
        {
            // Linearly interpolate to find the exact radius
            double m = DeltaR/(contRi-contLe);
            double contRadius = cPos(iRight-1,DeltaR) + m*(contR-contLe); // FindContRatio( contLe + m*(contR-contRi) , 1 ); // also override VContour with new values
            double contFinal = FindContRatio(contRadius,1); // 1 = Overwrites VContour
            cout << "Found a contour with radius " << contRadius << " with ratio " << contFinal << " (wanted contR= " << contR << ")" << endl;
            success=1;
            break;
        }

        // Are we out of options?
        if(iLeft==1 and iRight==M-1) { cout << "Did not find contour!!!" << endl; break;}

        // else Adjust indices
        if(iRight < M-1) iRight++;
        if(iLeft > 1) iLeft--;
        cout << "New iLeft and iRight are = " << iLeft << " and " << iRight << endl;
    }

    //exit(0);
    return success;

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
