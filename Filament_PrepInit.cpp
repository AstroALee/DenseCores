#include "Filament_PrepInit.H"

// Magnetized cylinder solution generator
void MagnetizedCylinder(double Redge);



void PrepareInitialState()
{
    // Determine Vpot, Apot, RhoTop from cylinder solution
    CodeHeader("Magnetized Cylinder integrations");
    InitCylinder();

    // Adds on Vpot from a chain of point masses
    CodeHeader("Determining Point Mass Potential");
    InitPoints();

    // We know the full initial potential. Get boundary conditions
    CodeHeader("Determining dVdR at the right edge");
    detDVDR();

    // Get boundary
    CodeHeader("Getting first filament boundary");
    getVbdy();

    // With V and rho known everywhere, Q is easy
    CodeHeader("Initializing Q");
    InitQ();

};






void InitCylinder()
{
    // Initialize Potentials from Cylinder
    // Given the value of rL and zL, we can find the cylinder solution that goes to unity at r=RADIALRATIO*zL,
    // while calcuating the full potential solution out to rL (assuming P=1 outside r=zL)

    // Radial edge of the cylinder
    double Redge = RADIALRATIO*zL;

    // We want the solution where rho = 1 at Redge
    MagnetizedCylinder(Redge);


};


void InitPoints()
{

    // If there's an excess mass, impose the point potential
    if(mExcess<=0) return;

    // Array for the point mass chain
    double PointV[M][N];
    double PointVold[M][N];

    int i,j;
    for(i=0;i<M;i++) for(j=0;j<N;j++) PointV[i][j] = 0.0;
    for(i=0;i<M;i++) for(j=0;j<N;j++) PointVold[i][j] = 0.0;

    // Loops, adding point masses to the chain. Breaks from additional points do not change potential
    int PPidx = 0;
    double zP=0;

    // Mass in code units
    double Mcode = mExcess * Sol2Code;

    while(true)
    {
        PPidx++;

        // Adds potential
        for(i=0;i<M;i++) for(j=0;j<N;j++)
        {
            // Distance locations of point particles
            double X1 = sqrt( pow(cPos(i,DeltaR),2) + pow(cPos(j,DeltaZ)-zP,2) );
            double X2 = sqrt( pow(cPos(i,DeltaR),2) + pow(cPos(j,DeltaZ)+zP,2) );

            if(i==0 && j==0 && zP==0)
            {
                // Do nothing (GM/X blows up)
            }
            else if(i==0 && j==1 && zP==0)
            {
                PointV[0][0] = PointV[0][0] + -1.1*Mcode/X1; // Since we skipped this earlier
                                                             // 1.1 factor so potential min at 0,0
                PointV[i][j] = PointV[i][j] + -Mcode/X1;
            }
            else
            {
                PointV[i][j] = PointV[i][j] + -Mcode/X1;
                if(zP > 0) PointV[i][j] = PointV[i][j] + -Mcode/X2;
            }

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
    double Mcode = mExcess * Sol2Code;
    double zP = 0;

    // Just in case
    for(i=0;i<N;i++) VRight[i] = 0;

    // Potential grad from point masses
    for(LoopIdx=0;LoopIdx<=PointLoopMax;LoopIdx++)
    {
        if(Mcode==0) break;

        for(j=0;j<N;j++)
        {
            double X1 = sqrt( pow(cPos(i,DeltaR),2) + pow(cPos(j,DeltaZ)-zP,2) );
            double X2 = sqrt( pow(cPos(i,DeltaR),2) + pow(cPos(j,DeltaZ)+zP,2) );

            VRight[j] = VRight[j] + Mcode/pow(X1,2);
            if(zP>0) VRight[j] = VRight[j] + Mcode/pow(X2,2);
        }
        zP = zP + 2.0*zL;
    }


    // The potential from the cylinder then adds the same value to each row
    // V has the form V = C1*ln(r) + C2
    // so dV/dr = C1/r
    // C1 = Vdot*rEdge, C2 = V(rEdge) - C1*ln(rEdge)
    // This assumes the right boundary is outside the filament

    double rEdge = RADIALRATIO*zL;

    for(i=0;i<N;i++) VRight[i] = VRight[i] + Ccst[0]/rL;

    for(i=0;i<N;i++) cout << VRight[i] << endl;

};

void getVbdy()
{
    double Vrow[M];
    double rEdge = RADIALRATIO*zL;

    int i,j;

    // Vpot value we're interested in
    for(i=0;i<M;i++) Vrow[i] = curState[Vpot][i][N-1];

    // Upate the global variable
    Vbdy = LIntY(Vrow,rEdge,DeltaR,M);
    cout << "Vbdy is " << Vbdy << endl;

    VContour[N-1] = rEdge; // Should have already done this

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
            cout << "EEK! Never found a cell (z=" << j << ")" << endl;
        }

        double y1 = curState[Vpot][idx-1][j];
        double y2= curState[Vpot][idx][j];
        double x1 = cPos(i-1,DeltaR);
        double x2 = cPos(i,DeltaR);
        double m = (y2-y1)/(x2-x1);

        VContour[j] = x1 + (Vbdy-y1)/m;

    }

    // Technically, if Mexcess > 0, rho does not drop to zero along this boundary. We technically
    // do not use rho in the finite difference solve, so while this is a little inconsistent, whatev.



};


void InitQ()
{
    // Q = rho * cs^2 exp(V/cs^2)

    // First, find Q value at boundary (where rho = 1 by definition)
    double Qbdy = 1.0 * exp(Vbdy);

    // Loop over all cells, if left of boundary, evaluate Q, else assign value at bdy
    int i,j;
    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        if( cPos(i,DeltaR) <= VContour[i] )
        {
            curState[Q][i][j] = curState[Rho][i][j]*exp(curState[Vpot][i][j]);
        }
        else
        {
            curState[Q][i][j] = Qbdy;
        }
    }


};
