#include "DenseCoreGlobals.H"
#include "ErrorMessages.H"


#include "nr3.h"
#include "stepper.h"
#include "odeint.h"
#include "stepperdopr5.h"


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
        setncst(0.5);   // <--- either 1/2 or 2/3
        
        setalph( sqrt(8.0*PI/beta) );
        double cs=1.0;
        double Gcst=1.0;
        
        // Derived quantities
        double alp2 = alp*alp;
        double nalpcst = -ncst*alp2/4.0/PI;
        
        dydx[0] = y[1];
        
        //Magnetostatic
        
        double r = t + 1e-20;
        double dPsi = nalpcst*pow(y[0],2.*ncst-2.)*y[1]-pow(cs,2)*y[1]/y[0];
        double ddPsi = nalpcst*(2.*ncst-2.)*pow(y[0],2.*ncst-3.)*pow(y[1],2) + pow(cs*y[1]/y[0],2);
        
        //dydx[1] = dPsi/r + ddPsi - 4.0*PI*Gcst*y[0] ;
        //dydx[1] = dydx[1] / ( pow(cs,2)/y[0] - nalpcst*pow(y[0],2.*ncst-2.) );
        
        dydx[1] = ( 4*PI*Gcst*pow(y[0],2) - dPsi*y[0]/r - nalpcst*(2*ncst-1)*pow(y[0],2*ncst-2)*y[1] + (nalpcst*pow(y[0],2*ncst-1)-cs*cs)*pow(y[1],2)/y[0] ) / ( nalpcst*pow(y[0],2*ncst-1)-cs*cs) ;
        
        
        //Hydrostatic
        //dydx[1] = pow(y[1],2)/y[0] - 4*PI*Gcst*pow(y[0],2)/(cs*cs) - y[1]/(t+0.00000001);
        
        //cout << dydx[1] << endl;
    }
    
};



void MagCylinder(double* density, double* initbeta)
{
    // Number of integrations
    int Nint = 0;
    int Nmax = 100;
    
    //Set for particular problem (independent variable = t)
    const int Nvar = 2;
    const double t_start = 0.0;
    const double t_end   = ((double)N+1)*DeltaR;
    
    
    //ODE tolerances and such
    const double atol = 1.0e-4;
    const double rtol = atol;
    const double h1 = 0.001;
    const double hmin = 1.0e-15;
    derivs derv;
    
    
    derv.setbeta(betaedge);
    
    // Output (NR version)
    int Nouts = 2*N;
    Output out(Nouts);
    double DeltaOut = t_end / ( (double) Nouts);
    
    
    
    
    //Initialize
    VecDoub ystart(Nvar);
    ystart[0] = 1.0;
    ystart[1] = 0.0;
    
    //Get ready to integrate
    Odeint< StepperDopr5<derivs> > ode(ystart,t_start,t_end,atol,rtol,h1,hmin,out,derv);
    //Odeint< StepperRoss<derivs> > ode(ystart,t_start,t_end,atol,rtol,h1,hmin,out,derv);
    
    int i,j,k;
    double betatol = 0.1;
    while(true)
    {
        
        // INTEGRATE!
        ode.integrate();
        
        // Interpolate onto our grid
        for(i=0;i<N+1;i++)
        {
            // for grid cell i
            // find the two output cells it lies between
            for(j=1;j<out.count;j++)
            {
                if( out.xsave[j] > (0.5+(double)i)*DeltaR )
                {
                    double Xsave = (0.5+(double)i)*DeltaR;
                    
                    double Yright = out.ysave[0][j];
                    double Yleft = out.ysave[0][j-1];
                    double Xright = out.xsave[j];
                    double Xleft = out.xsave[j-1];
                    double m = (Yright-Yleft)/(Xright-Xleft);
                    
                    density[i] = Yleft + m*(Xsave-Xleft);
                    break;
                    
                }
            }
            
        }
    
        // Now Interpolated, check the value of Beta at the upper-right edge
        // First need to find the density value at R = VContour[Z];
        double density_edge;
        
        for(j=1;j<out.count;j++)
        {
            if( out.xsave[j] > VContour[Z] )
            {
                double Xsave = VContour[Z];
                
                double Yright = out.ysave[0][j];
                double Yleft = out.ysave[0][j-1];
                double Xright = out.xsave[j];
                double Xleft = out.xsave[j-1];
                double m = (Yright-Yleft)/(Xright-Xleft);
                
                density_edge = Yleft + m*(Xsave-Xleft);
                cout << "density edge " << density_edge << endl;
                break;
                
            }
        }
        
        // Get beta by analytic relationships
        double beta_top = derv.ret_beta() * pow(density_edge,1-2*derv.ret_n());
        
        
        // If this matches betaedge, we're done. Else iterate until we converge.
        double errbeta = fabs(beta_top-betaedge)/betaedge;
        
        if(errbeta < betatol) break;
        
        // Create a new guess and start again
        derv.setbeta( betaedge * pow(density_edge,2*derv.ret_n()-1) );
        
        Nint++;
        if(Nint == Nmax) WATERLOO_MagCylNoConverge(Nint);
        
            
    }
    
    // Interpolate beta onto grid
    for(j=0;j<N+1;j++) initbeta[j] = derv.ret_beta() * pow(density[j],1-2*derv.ret_n());
   
    
    
    cout << "Integration Done!" << endl;
    cout << "Number good steps: " << ode.nok << endl;
    cout << "Number bad steps:  " << ode.nbad << endl << endl;
 
    
    
    
    
    // Now interpolate the answer onto the density array and beta onto the density array
    
    

    
}