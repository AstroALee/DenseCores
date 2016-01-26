// Dense Core files
#include "DenseCoreGlobals.H"
#include "ErrorMessages.H"

// Numerical Recipes
#include "nr3.h"
#include "stepper.h"
#include "odeint.h"
#include "stepperdopr5.h"

double errbeta(double beta_top) { return fabs(beta_top-betaedge)/betaedge; }
double errrho(double density_edge) { return fabs(density_edge-1.0)/1.0; }


// Linearly Interpolate an array onto a point. Returns value.
inline double LInt(VecDoub Xarray, MatDoub Yarray, double point, int size,int elem)
{
    double value = 0;
    int found = 0;

    // find the two points surrounding 'point' in the X array
    // Assumes X array is arranged in ascending value
    int i;
    for(i=0;i<size-1;i++) if( Xarray[i+1] >= point)
    {
        double RightY = Yarray[elem][i+1];
        double LeftY  = Yarray[elem][i];
        double RightX = Xarray[i+1];
        double LeftX  = Xarray[i];

        double m = (RightY-LeftY)/(RightX-LeftX);

        value = LeftY + m*(point-LeftX);
        found = 1;
        break;
    }

    if(!found) cout << "Interpolated value was not found!" << endl;

    return value;
}





struct derivs {

    void setbeta(double b){ beta = b;}
    void setalph(double a){ alp = a;}
    void setncst(double n){ ncst = n;}

    double ret_alpha(){return alp;}
    double ret_beta(){return beta;}
    double ret_n(){return ncst;}

    double beta;
    double alp;
    double ncst;

    // Derivatives
    void operator() (const double t, VecDoub_I &y, VecDoub_O &dydx) {

        // Fixed parameters?
        double cs=1.0;
        double Gcst=1.0;

        // Derived quantities
        double alp2 = alp*alp;
        double nalpcst = -ncst*alp2/4.0/PI;

        // d(rho) = Rhodot
        dydx[0] = y[1];

        //Magnetostatic

        double r = t + 1e-20;
        double dPsi = nalpcst*pow(y[0],2.*ncst-2.)*y[1]-pow(cs,2)*y[1]/y[0];
        double ddPsi = nalpcst*(2.*ncst-2.)*pow(y[0],2.*ncst-3.)*pow(y[1],2) + pow(cs*y[1]/y[0],2);

        //dydx[1] = dPsi/r + ddPsi - 4.0*PI*Gcst*y[0] ;
        //dydx[1] = dydx[1] / ( pow(cs,2)/y[0] - nalpcst*pow(y[0],2.*ncst-2.) );

        // d(Rhodot) = RhoDotDot
        dydx[1] = ( 4*PI*Gcst*pow(y[0],2) - dPsi*y[0]/r - nalpcst*(2*ncst-1)*pow(y[0],2*ncst-2)*y[1] + (nalpcst*pow(y[0],2*ncst-1)-cs*cs)*pow(y[1],2)/y[0] ) / ( nalpcst*pow(y[0],2*ncst-1)-cs*cs) ;


        // Tracks potential as well (V = entry 2, Vdot = entry 3)

        // d(V) = Vdot
        dydx[2] = y[3];
        // Vdotdot
        dydx[3] = 4*PI*y[0] - y[3]/r;



    }

};



void MagCylinder(double* density, double* initbeta, double* an_pot)
{
    // Number of integrations
    int Nint = 0;
    int Nmax = 50;

    // Set for particular problem (independent variable = t)
    // Integrate to the contour
    const int Nvar = 4;
    const double t_start = 0.0;
    const double t_end   = pLength - 0.5*DeltaR + 0.0001*DeltaR; //VContour[Z];


    //ODE tolerances and such
    const double atol = 1.0e-5;
    const double rtol = atol;
    const double h1 = 0.001;
    const double hmin = 1.0e-15;
    derivs derv;


    derv.setncst(0.5);                      // <- PARAMETER: either 1/2 or 2/3
    derv.setbeta(betaedge);
    derv.setalph( sqrt(8.0*PI/betaedge) );  // <- determined at edge where rho=1


    // Output (NR version)
    int Nouts = 2*N;
    double DeltaOut = t_end / ( (double) Nouts);
    double VpotRho[N];

    //Initial guess, use analytic relationship (when n = 1/2, this will be exact)
    VecDoub ystart(Nvar);
    double yfirstguess = 1.86;
    double ysecondguess;
    double errorfirst;
    ystart[0] = yfirstguess;
    ystart[1] = 0.0;
    ystart[2] = 0.0;
    ystart[3] = 0.0;

    //Get ready to integrate

    //Odeint< StepperRoss<derivs> > ode(ystart,t_start,t_end,atol,rtol,h1,hmin,out,derv);

    int i=0,j,k;
    double betatol = 0.0001;
    double rhotol = 0.0001;

    while(true)
    {

        // INTEGRATE!
        Output out(Nouts);
        Odeint< StepperDopr5<derivs> > ode(ystart,t_start,t_end,atol,rtol,h1,hmin,out,derv);
        cout << "Integrating." << endl;
        ode.integrate();


        // Interpolate density onto our grid
        for(i=0;i<N;i++) density[i] = LInt(out.xsave, out.ysave, cPos(i,DeltaR), out.count,0);

        // Interpolate Potential onto our grid
        for(i=0;i<N;i++) an_pot[i] = LInt(out.xsave,out.ysave,cPos(i,DeltaR),out.count,2);


        // Now Interpolated, check the value of Beta and rho at the upper-right edge
        // First need to find the density value at R = VContour[Z];
        double density_edge;
        density_edge = LInt(out.xsave,out.ysave,VContour[Z-1],out.count,0);

        // Test
#if DEBUG
            //density_edge = LInt(out.xsave,out.ysave,VContour[Z-1],out.count,0);
            //density_edge = LInt(out.xsave,out.ysave,cPos(N-1,DeltaR),out.count,0);
            //ContourVpot = an_pot[N-1];
#endif


        // Get beta by analytic relationships
        double beta_top = derv.ret_beta() * pow(density_edge,1-2*derv.ret_n());
        cout << "density edge = " << density_edge << " and beta edge = " << beta_top << endl;


        // If these match betaedge and P/c^2, we're done. Else iterate until we converge.
        // Units of P_0 = c_s = 1 implies P/c^2 = 1.
        if( errbeta(beta_top) < betatol  &&  errrho(density_edge) < rhotol)
        {
            // Success?
            cout << "Success with error = " << errrho(density_edge) << endl;
            break;
        }


        // Create a new guess and start again
        if(Nint==0)
        {
            errorfirst = errrho(density_edge);
            ysecondguess = yfirstguess*1.1;
            ystart[0] = ysecondguess;
            ystart[1] = 0.0;
            ystart[2] = 0.0;
            ystart[3] = 0.0;
        }
        else
        {
            //Newton-Rhapson
            double dennew;

            double m = (errrho(density_edge) - errorfirst)/(ysecondguess-yfirstguess);
            dennew = ysecondguess - (errrho(density_edge) / m);

            if(dennew < 0) dennew = ysecondguess*1.1;

            ystart[0] = dennew;
            ystart[1] = 0.0;
            ystart[2] = 0.0;
            ystart[3] = 0.0;

            // Updates for next Newton-Rhapson, if needed
            yfirstguess = ysecondguess;
            ysecondguess = dennew;
            errorfirst = errrho(density_edge);
        }


        cout << "Trying new integration with beta = " << derv.ret_beta() << " and density = " << ystart[0] << endl;
        cout << "error = " << errrho(density_edge) << endl;


        Nint++;
        if(Nint == Nmax) WATERLOO_MagCylNoConverge(Nint);


    }


    // Interpolate beta onto grid
    for(j=0;j<N;j++) initbeta[j] = derv.ret_beta() * pow(density[j],1-2*derv.ret_n());


    cout << "Integration Done!" << endl;
    //cout << "Number good steps: " << ode.nok << endl;
    //cout << "Number bad steps:  " << ode.nbad << endl << endl;



    cout << "Density values: " << density[0];
    for(j=1;j<N;j++) cout << ", " << density[j];
    cout << endl;


    cout << "Potential values: " << an_pot[0];
    for(j=1;j<N;j++) cout << ", " << an_pot[j];
    cout << endl;





}
