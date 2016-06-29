#include "Filament_newSolveAll.H"

// Functions from other files
void SolvePoisson();
void SolveAmpere();



void Converge()
{

    int done = 0;
    int dummy = 0;

    // Outermost Loop -- Seeks convergence in V, A, Q, and Rho
    int i;
    for(i=0;i<ConvergeLoopMax;i++)
    {
        cout << endl << "Outer loop: Beginning iteration " << i << endl;

        // Inner loop -- Seeks convergence in V, A, Q, fixes Rho
        int j;
        for(j=0;j<ConvergeLoopMax;j++)
        {
            cout << endl << "Inner loop: Beginning iteration " << j << endl;

            // Solve Poisson -- solution stored in newState
            SolvePoisson();

            // Solve Ampere -- solution stored in newState
            SolveAmpere();

            // About to use new solution to relax, make data in curState the prevState
            CopyState(curState,prevState);

            // Relax V and A using newState and prevState, use this to recalc boundary, Q, and Rho
            RelaxSoln();

            // Update Q via mapping
            UpdateQ(0);
            //UpdateQ(1);

            // After this inner loop, we keep Rho fixed
            // and return to solving Poisson and Ampere
            if(ConvergeTest(i,1)) { cout << "Converged! Breaking out." << endl; break; }

        } // end innermost loop

        //NewUpdateRho(0);


        if(UpdateRho(0)) dummy=UpdateRho(1);
        else { cout << "Did not find contour!" << endl; }
        UpdateQ(1);

        // debug
        //break;

        if(ConvergeTest(i,2)) { cout << "Converged in Rho! Breaking out." << endl; break; }

    } // end outermost loop


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


// Compares the previous state to the current updated state to assess if we have converged to a SS
int ConvergeTest(int loopnum,int type)
{
    int done = 0;

    double errV=0, errA=0, errR=0;

    int i,j;
    int Vi,Vj,Ai,Aj;

    if(type==1) {

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
    //PrintError(loopnum,errV,errA);

    // Can we go home yet?
    if(errV < ConvergeTol && errA < ConvergeTol) done = 1;

    string a;
    if(done) a = "SUCCESS!"; else a = " ";
    cout << "Converge iteration " << loopnum << " : " << a << " Max errors V,A = " << errV << " , " << errA << endl;

    }
    else
    {
        cout << "Test for convergence in Rho!" << endl;

        double locerrR=0;

        for(i=0;i<M;i++) for(j=0;j<N;j++)
        {
            locerrR = 0;

            if( cPos(i,DeltaR) < VContour[j]-DeltaR )
            {
                if( fabs( prevState[Rho][i][j] ) > 0 )
                    locerrR = fabs( curState[Rho][i][j] - prevState[Rho][i][j] )/( prevState[Rho][i][j] );
                else
                    locerrR = 0.0;
            }

            if(locerrR > errR) errR = locerrR;

        }
        cout << "Rho err was " << errR << endl;
        if(errR < ConvergeTol) done = 1;

    }



    return done;
};


void NewUpdateRho(int type)
{
    int i,j;

    if(type)
    {
        // We know where the contour is
        for(i=0;i<M;i++)for(j=0;j<N;j++)
        {
            if(cPos(i,DeltaR)>VContour[j]) curState[Rho][i][j] = 0;
        }
        return;
    }


    // Use the new V and Q to get a new rho everywhere
    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        if(j==N-1) curState[Rho][i][j] = RhoTop[i]; // top row is always the cylinder
        else curState[Rho][i][j] = curState[Q][i][j]*exp(-curState[Vpot][i][j]); // everywhere else uses the Q update
    }

    // Find the density contour the encloses the right amount of mass.
    double tempCont[N];

    double prevMass = 0, currMass = 0;
    // Find the contour that starts at a partiular top location.
    for(i=1;i<M;i++)
    {
        // Find the contour the starts at cPos(i)
        findDenCont(tempCont,cPos(i,DeltaR),0);
        // Given this contour, find the mass enclosed
        currMass = findMass(tempCont);
        // As i increases, mass will increase. Have we enclosed the right amount of mass?
        if(currMass >= 2.0*zL*lambda)
        {
            // The right amount of mass is somewhere between i and i-1
            // Linearly interpolate to get the guess for the top radius
            double m = DeltaR/(currMass-prevMass);
            double rightRad = cPos(i-1,DeltaR) + m*(2.0*zL*lambda - prevMass);
            // This is the contour we want.
            findDenCont(tempCont,rightRad,0);
            // Save in our permanent array.
            for(j=0;j<N;j++) VContour[j] = tempCont[j];
            break;
        }
        else
        {
            // Not enough mass yet. Need to keep going
            prevMass = currMass;
        }
    }

};

void findDenCont(double *rhoContour, double topRad, int type)
{
    // What is the density at topRad?
    int i,j;
    i=0; while(true) { i++; if(cPos(i,DeltaR)>=topRad) break; }

    double rRight = curState[Rho][i][N-1];
    double rLeft = curState[Rho][i-1][N-1];
    double m = (rRight-rLeft)/DeltaR;

    double curDen = rLeft + m*(topRad-cPos(i-1,DeltaR));

    // Top is done
    rhoContour[N-1] = topRad;

    // For each other row, find the location where rho = curDen
    // Assumes rho monotonically decreases
    for(j=0;j<N-1;j++)
    {
        i=0; while(true){ i++; if(curState[Rho][i][j] <= curDen) break; }
        rRight = curState[Rho][i][j];
        rLeft = curState[Rho][i-1][j];
        m = (rRight-rLeft)/DeltaR;

        double locDen = rLeft + m*(topRad-cPos(i-1,DeltaR));

        // This row is done
        rhoContour[j] = locDen;
    }

};

// Do an integral of rho, ignoring the density beyond rhoContour
double findMass(double *rhoContour)
{
        double sum=0;

        // finds the largest index
        int idxMax = 0;
        int i,j;


        for(j=0;j<N;j++) for(i=0;i<M;i++)
        {
            if( cPos(i,DeltaR) >= rhoContour[j] )
            {
                if(i>idxMax) idxMax=i;
                break;
            }
        }

        // Now do a trapezoidal integral out to idxMax
        sum = cPos(0,DeltaR)*curState[Rho][0][0] + cPos(idxMax,DeltaR)*curState[Rho][idxMax][0] + cPos(0,DeltaR)*curState[Rho][0][N-1] + cPos(idxMax,DeltaR)*curState[Rho][idxMax][N-1]; // corners
        for(i=1;i<idxMax;i++) sum = sum + 2.0*(cPos(i,DeltaR)*curState[Rho][i][0] + cPos(i,DeltaR)*curState[Rho][i][N-1]); // top and bottom edges
        for(j=1;j<N-1;j++) sum = sum + 2.0*(cPos(0,DeltaR)*curState[Rho][0][j] + cPos(idxMax,DeltaR)*curState[Rho][idxMax][j]); // left and right edges
        for(i=1;i<idxMax;i++) for(j=1;j<N-1;j++) sum = sum + 4.0*(cPos(i,DeltaR)*curState[Rho][i][j]); // interiors
        sum = 4.0*PI*(DeltaR*DeltaZ/4.0)*sum;

        return sum;

};


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
};
