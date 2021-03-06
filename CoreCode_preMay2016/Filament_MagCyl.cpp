#include "Filament_MagCyl.H"


// Derivatives for OD solver
struct derivs {

    // Derivatives
    void operator() (const double t, VecDoub_I &y, VecDoub_O &dydx) {

        // Fixed parameters?
        double ncst= nCyl;  //0.5;
        double alp= sqrt(8.0*PI/betaInf); // alp^2 = 8*pi*rho^(1-2n)_intfy / beta_infty

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
            //ddrho = -2.0*PI*pow(rho,3)/(rho+gm*pow(rho,2.*ncst));
            ddpsi = 2.0*PI*rho;

        }
        else
        {
            //ddrho = 4*PI*pow(rho,3) - (1+gm*pow(rho,2.*ncst-1))*pow(rhop,2) + gm*(2*ncst-1)*pow(rho,2*ncst-1)*pow(rhop,2) + (rho+gm*pow(rho,2.*ncst))*rhop/r;
            //ddrho = -ddrho/(rho+gm*pow(rho,2.*ncst));

            //dpsi = -rhop/rho - gm*rhop*pow(rho,2*ncst-2);
            ddpsi = 4.0*PI*rho - dpsi/r;
        }


        drA = r*(alp*pow(rho,ncst));

        double rhop = -rho*dpsi/(1.0+gm*pow(rho,2.*ncst-1));


        dydx[0] = rhop; //y[1]; // d/dr(rho) = rho'
        //dydx[1] = ddrho; // d/dr(rho') = rho''

        dydx[1] = dpsi; // d/dr(psi) = psi'
        dydx[2] = ddpsi; // d/dr(psi') = psi''

        dydx[3] = drA;   // d/dr(r*A) = (r*A)'

        // Mass/length integral of cylinder
        dydx[4] = 2.0*PI*r*rho; // d/dr(lambda) = 2*pi*r*rho

    }

};

void MagnetizedCylinder(double rEdge,double& curMass,int dooutput)
{
    //Set for particular problem (independent variable = t)
    const int Nvar = 5;
    const double t_start = 0.0;
    const double t_end   = 1.2*rEdge;

    //ODE tolerances and such
    const double atol = 1.0e-8;
    const double rtol = atol;
    const double h1 = 0.001;
    const double hmin = 1.0e-6;
    derivs derv;  // Only need one declaration of this

    // going to Newton-Raphson until we find the value of rho_center where rho=1 at rEdge
    int NRidx=0, NRmax = LoopMAX; // set in Filament_Main.H
    int NRout=10*M;  // need a good number for accurate interpolations

    // new (first) guess (drawn from averages of typical cases)
    double rhoNG = 10.0*rEdge;
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

        // Linearly interpolate the value of rho at r = rEdge from our solution
        double rhoSoln = LInt(out.xsave,out.ysave,rEdge,out.count,0);
        // Update error
        NRerrNG = fabs(rhoSoln-1.0)/1.0; // want rhoSoln to be = 1.0
        cout << "Num " << NRidx << " : rho(rEdge) = " << rhoSoln << " (error = " << NRerrNG << ")" << endl;


        if(NRerrNG < NRtol)
        {
            cout << "rEdge = " << rEdge << endl;

            // We succeeded. Use the output to initialize things
            UseOutput(out,rEdge,curMass);

            // break out of this loop
            break;
        }
        else
        {
            // Prepare for next integration

            // If this is our first time here?
            if(rhoOG == 0)
            {
                // We have to get a second guess so we can calculate derivatives
                rhoOG = rhoNG;
                rhoNG = rhoOG*1.1;

                NRerrOG = NRerrNG;
            }
            else
            {
                // Use the POWER of Newton Rhapson

                // slope
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

            // Rho start
            ystart[0] = rhoNG;
            // Psi and Psi' start
            ystart[1] = 0;
            ystart[2] = 0;
            // A start
            ystart[3] = 0;
            // Mass/length start
            ystart[4] = 0;

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

}



void UseOutput(Output out,double rEdge,double& curMass)
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

    if(1)
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

    // V outside has the form  V = C1*ln(r) + C2
    // C1 and C2 chosen so V is smooth at rEdge (V(rE) and DV(rE) match)
    // C1 = DV(rE)*rEdge, C2 = V(rE) - C1*ln(rE)
    double C1 = VEdgeD*rEdge;
    double C2 = VEdge - C1*log(rEdge);

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

            if(i>0 && DEBUG==2) curState[Apot][i][j] = curState[Apot][i][j] + 0.3*fabs(rPos*zPos)*(1.0/zL/rL);

        }
        else
        {
            // Analytic form for V
            curState[Vpot][i][j] = Ccst[0]*log(rPos)+Ccst[1] ;

            // Analytic form for A
            double Binf = sqrt(8.0*PI/betaInf);
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

    // The total mass of the filament
    curMass = 2.0*zL*LInt(out.xsave,out.ysave,rEdge,out.count,4) / Sol2Code;
    if(totMass == -1) totMass = curMass;
    cout << "Cylinder calculation total mass = " << totMass << endl;

};
