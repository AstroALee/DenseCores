#include "Filament_newMagCyl.H"


// Integration works in units of cs = G = 1.
// We are going to find the set of (rho_c,r) such that rho(r)=1 and lambda(code) = lambda(input)

// Derivatives for OD solver
struct derivs {

    // Derivatives
    void operator() (const double t, VecDoub_I &y, VecDoub_O &dydx) {

        // Fixed parameters?
        double ncst= nCyl;  //0.5;
        double alp= sqrt(8.0*PI/beta0); // alp^2 = 8*pi*rho^(1-2n)_intfy / beta_infty

        // Derived quantities
        double a2 = alp*alp;
        double gm = a2*ncst/4.0/PI;

        // Use r instead of t (for sanity)
        double r = t;
        double rho = y[0]; // we'll use this often
        double psi = y[1];
        double dpsi = y[2];

        // Second derivative place holder
        double ddpsi;

        // First derivative of r*A
        double drA;

        if(r==0) // Likely derived using L'Hopital
        {
            ddpsi = 2.0*PI*rho;
        }
        else
        {
            ddpsi = 4.0*PI*rho - dpsi/r;
        }


        drA = r*(alp*pow(rho,ncst));

        double drho = -rho*dpsi/(1.0+gm*pow(rho,2.*ncst-1));


        dydx[0] = drho; //y[1]; // d/dr(rho) = rho'
        //dydx[1] = ddrho; // d/dr(rho') = rho''

        dydx[1] = dpsi; // d/dr(psi) = psi'
        dydx[2] = ddpsi; // d/dr(psi') = psi''

        dydx[3] = drA;   // d/dr(r*A) = (r*A)'

        // Mass/length integral of cylinder
        dydx[4] = 2.0*PI*r*rho; // d/dr(lambda) = 2*pi*r*rho

    }

};

void FastMagnetizedCylinder(double& rEdge, double desiredLambda, int dooutput)
{
    //Set for particular problem (independent variable = t)
    const int Nvar = 5;
    const double t_start = 0.0; // t here is the radius
    const double t_end   = max(3.0,1.5*rL); // this should be large enough so that lambda = desiredLambda somewhere
                                            // also should be larger than rL

    //ODE tolerances and such
    const double atol = 1.0e-8;
    const double rtol = atol;
    const double h1 = 0.001;
    const double hmin = 1.0e-6;
    derivs derv;  // Only need one declaration of this

    int NRout = 10*M;  // need a good number for accurate interpolations

    // Initialize intitial conditions
    VecDoub ystart(Nvar);
    // Rho start
    ystart[0] = 1.0;
    // Psi and Psi' start
    ystart[1] = 0;
    ystart[2] = 0;
    // A start
    ystart[3] = 0;
    // Mass/length start
    ystart[4] = 0;

    // Output
    Output out(NRout);

    // Integrate
    Odeint< StepperDopr5<derivs> > ode(ystart,t_start,t_end,atol,rtol,h1,hmin,out,derv);
    ode.integrate();

    //for(int i=0;i<out.count;i++) cout << "output " << i << " : " << out.ysave[4][i] << endl;
    // What is the r where lambda=desiredLambda?
    rEdge = LIntY(out.xsave, out.ysave, desiredLambda, out.count, 4); // 4 = lambda, returns radius
    cout << "Value of rEdge is " << rEdge << endl;

    if(dooutput==2)
    {
        NewUseOutput(out,rEdge);
    }
    else if(dooutput==1)
    {
        // The cylinder solution is good to go.
        for(int i=0 ; i<M; i++) RhoTop[i] = LInt(out.xsave,out.ysave,cPos(i,DeltaR),out.count,0); // 0 = Rho
        for(int i=0 ; i<M; i++) cout << "RhoTop(" << i << ")=" << RhoTop[i] << endl;
        //for(int i=0 ; i<M; i++) cout << "Lambda(" << i << ")=" << LInt(out.xsave,out.ysave,cPos(i,DeltaR),out.count,4) << endl;
    }


};

void MagnetizedCylinder(double& rEdge, double desiredLambda, int dooutput)
{
    //Set for particular problem (independent variable = t)
    const int Nvar = 5;
    const double t_start = 0.0; // t here is the radius
    const double t_end   = 2.0; // this should be large enough so that lambda = desiredLambda somewhere

    //ODE tolerances and such
    const double atol = 1.0e-8;
    const double rtol = atol;
    const double h1 = 0.001;
    const double hmin = 1.0e-6;
    derivs derv;  // Only need one declaration of this

    // going to Newton-Raphson until we find the value of rho_center where lambda = desiredLambda
    int NRidx = 0, NRmax = LoopMAX; // set in Filament_Main.H
    int NRout = 1000; //10*M;  // need a good number for accurate interpolations


    double rhoNG;
    // old:: new (first) guess (drawn from averages of typical cases)
    //if(beta0 < 0.5) rhoNG = 4.0;  // 'intuition'
    //else rhoNG = 1.21;
    rhoNG = 1.0;

    double rhoOG = 0;   // old guess
    double NRerrNG = 1;   // error in rho
    double NRerrOG = 1;
    double NRtol = CylTol; // when do we claim success? ; set in Filament_Main.H

    // Newton Rhapson routine
    while(true)
    {
        NRidx++;

        // Initialize intitial conditions
        VecDoub ystart(Nvar);
        // Rho start
        ystart[0] = rhoNG;
        // Psi and Psi' start
        ystart[1] = 0;
        ystart[2] = 0;
        // A start
        ystart[3] = 0;
        // Mass/length start
        ystart[4] = 0;

        // Output
        Output out(NRout);

        // Integrate
        Odeint< StepperDopr5<derivs> > ode(ystart,t_start,t_end,atol,rtol,h1,hmin,out,derv);
        ode.integrate();

        // Let's see how we did

        // What is the r where lambda=desiredLambda?
        rEdge = LIntY(out.xsave, out.ysave, desiredLambda, out.count, 4); // 4 = lambda

        // What is lambda at this radius?
        double curLam = LInt(out.xsave,out.ysave,rEdge,out.count,4); // 4 = lambda

        // What's the error?
        NRerrNG = fabs(curLam-desiredLambda)/desiredLambda;
        cout << "Num " << NRidx << " : lambda(rEdge) = " << curLam << " (error = " << NRerrNG << ")" << endl;

        if(NRerrNG < NRtol)
        {
            cout << "Converged: rEdge = " << rEdge << " and lambda = " << curLam << endl;

            // If we're using this to set up initial conditions, do it here
            if(dooutput) UseOutput(out,rEdge);

            // get out of the loop
            break;
        }
        else
        {
            // Prepare for next integration

            // Is this our first time here, set next guess, else Newton Rhapson
            if(rhoOG==0)
            {
                rhoOG = rhoNG;
                rhoNG = 1.1*rhoOG;
                NRerrOG = NRerrNG;
            }
            else
            {
                //slope
                double m = (NRerrNG - NRerrOG) / (rhoNG - rhoOG);
                rhoOG = rhoNG;
                rhoNG = rhoOG - NRerrNG/m; // rhoOG is the old rhoNG... a little confusing
                NRerrOG = NRerrNG;

                // If rhoNG < 1, yikes.
                if(rhoNG < 1)
                {
                    int rn = rand() % 30 + 1;
                    rhoNG = 1.0 + 0.01 * ((double) rn );
                    cout << "      Had to use random numbers to make a rhoNewGuess(0)" << endl;
                }
            }

            cout << "      : rhoNewGuess(0) = " << rhoNG << endl;

        }

        // Has something gone wrong?
        if(NRidx == NRmax)
        {
            WaterlooHeader("MagCyl");
            cout << "Did not converge to a cylinder solution! (NRmax = " << NRmax << ")" << endl;
            exit(1);
        }

    } // End of NR loop

};


void NewUseOutput(Output out, double rEdge)
{
    // We have already done Rho, we only need the values for V and A here

    // Let's fill up the state vector
    int i,j;

    // Rho at top edge
    double RTopEdge = LInt(out.xsave,out.ysave,rEdge,out.count,0); // Rho

    // Outside the filament, we use the analytic form for the potential
    double VEdge = LInt(out.xsave,out.ysave,rEdge,out.count,1); // Vpot
    double VEdgeD = LInt(out.xsave,out.ysave,rEdge,out.count,2); // d(Vpot)/dR

    // V outside has the form  V = C1*ln(r/rEdge) + C2
    // C1 and C2 chosen so V is smooth at rEdge (V(rE) and DV(rE) match)
    // C1 = rE*DV(rE)     C2 = V(rE)
    double C1 = VEdgeD*rEdge;
    double C2 = VEdge;

    // We will need these constants later
    Ccst[0] = C1;
    Ccst[1] = C2;

    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        // Assign Rho
        curState[Rho][i][j] = RhoTop[i];

        double rPos = cPos(i,DeltaR);
        double zPos = cPos(j,DeltaZ);

        if(rPos <= rEdge)
        {
            curState[Vpot][i][j] = LInt(out.xsave,out.ysave,rPos,out.count,1); // 1 = Vpot

            // We solved for r*A, have to divide by r
            if(i>0) curState[Apot][i][j] = LInt(out.xsave,out.ysave,rPos,out.count,3)/rPos; // 3 = Apot*r

        }
        else
        {
            // Rho is zero
            curState[Rho][i][j] = 0.0;
            
            // Analytic form for V
            curState[Vpot][i][j] = Ccst[0]*log(rPos/rEdge)+Ccst[1] ;

            // Analytic form for A*r outside cylinder = Cst*ln() + 0.5*Binf*(r^2-Redge^2)
            double Binf = sqrt(8.0*PI/beta0)*pow(RTopEdge,nCyl);
            double Acyl = LInt(out.xsave,out.ysave,rEdge,out.count,3) + 0.5*Binf*( pow(rPos,2) - pow(rEdge,2) );
            Acyl = Acyl/rPos;
            curState[Apot][i][j] = Acyl;

        }
    }


    // Want Vpot = 0 at top left corner
    // Should have that by default at this point
    double Vnorm = curState[Vpot][0][N-1];
    for(i=0;i<M;i++) for(j=0;j<N;j++) curState[Vpot][i][j] = curState[Vpot][i][j] - Vnorm;

};

void UseOutput(Output out,double rEdge)
{
    // The passed in output is the solution we like. Let's use it!

    // Clear stuff
    for(int i=0;i<M;i++) RhoTop[i]=0.0;

    for(int i=0;i<M;i++)
        for(int j=0;j<N;j++)
            for(int k=0;k<NStates;k++) curState[k][i][j] = 0.0;


    // We know Rho, so let's fill in RhoTop (this will never change)
    for(int i=0 ; i<M; i++)
    {
        double curPos = cPos(i,DeltaR);
        if(curPos <= rEdge) RhoTop[i] = LInt(out.xsave,out.ysave,curPos,out.count,0); // 0 = Rho
        else RhoTop[i] = 0.0;
    }

    if(0)
    {
    cout << "RhoTop:" << endl;
    for(int i=0 ; i<M; i++) cout << RhoTop[i] << ",";
    cout << endl << "Radius:" << endl;
    for(int i=0 ; i<M; i++) cout << cPos(i,DeltaR) << ",";
    cout << endl;
    }


    // Let's fill up the state vector

    int i,j;

    // Outside the filament, we use the analytic form for the potential
    double VEdge = LInt(out.xsave,out.ysave,rEdge,out.count,1); // Vpot
    double VEdgeD = LInt(out.xsave,out.ysave,rEdge,out.count,2); // d(Vpot)/dR

    // V outside has the form  V = C1*ln(r/rEdge) + C2
    // C1 and C2 chosen so V is smooth at rEdge (V(rE) and DV(rE) match)
    // C1 = rE*DV(rE)     C2 = V(rE)
    double C1 = VEdgeD*rEdge;
    double C2 = VEdge;

    // We will need these constants later
    Ccst[0] = C1;
    Ccst[1] = C2;

    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        double rPos = cPos(i,DeltaR);
        double zPos = cPos(j,DeltaZ);

        if(rPos <= rEdge)
        {
            curState[Vpot][i][j] = LInt(out.xsave,out.ysave,rPos,out.count,1); // 1 = Vpot
            //if(rPos==rEdge) cout << "Cyl Vpot " << curState[Vpot][i][j] << endl;

            curState[Rho][i][j]  = LInt(out.xsave,out.ysave,rPos,out.count,0); // 0 = Rho

            // We solved for r*A, have to divide by r
            if(i>0) curState[Apot][i][j] = LInt(out.xsave,out.ysave,rPos,out.count,3)/rPos; // 3 = Apot*r

            //if(i>0 && DEBUG==2) curState[Apot][i][j] = curState[Apot][i][j] + 0.3*fabs(rPos*zPos)*(1.0/zL/rL);

        }
        else
        {
            // Analytic form for V
            curState[Vpot][i][j] = Ccst[0]*log(rPos/rEdge)+Ccst[1] ;

            // Analytic form for A*r outside cylinder = Cst*ln() + 0.5*Binf*(r^2-Redge^2)
            double Binf = sqrt(8.0*PI/beta0); // assumes rho(edge)=1
            double Acyl = LInt(out.xsave,out.ysave,rEdge,out.count,3) + 0.5*Binf*( pow(rPos,2) - pow(rEdge,2) );
            Acyl = Acyl/rPos;
            curState[Apot][i][j] = Acyl;

            // Rho is easy in this region
            curState[Rho][i][j] = 0.0;
        }
    }


    // Want Vpot = 0 at top left corner
    // Should have that by default at this point
    double Vnorm = curState[Vpot][0][N-1];
    for(i=0;i<M;i++) for(j=0;j<N;j++) curState[Vpot][i][j] = curState[Vpot][i][j] - Vnorm;

    // Also, the exact location of the filament boundary is fixed at the top
    VContour[N-1] = rEdge;


    // Debug
    if(DEBUG==2)
    {
        // Makes the potentials wrong. The boundary conditions should make the solution settle
        // back to the analytic solution above.

        for(i=0;i<M;i++) for(j=0;j<N;j++)
        {
            double rloc = cPos(i,DeltaR);
            double zloc = cPos(j,DeltaZ);

            curState[Vpot][i][j] = curState[Vpot][i][j] + 1.0*sin(1.0*PI/rL*rloc)*sin(1.0*PI/zL*zloc);

            curState[Apot][i][j] = curState[Apot][i][j] + 1.0*sin(1.0*PI/rL*rloc)*sin(1.0*PI/zL*zloc);
        }


    }


};
