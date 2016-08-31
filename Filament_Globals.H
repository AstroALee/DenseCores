/*
   Global Variable Declarations for the Filament Code by Aaron T. Lee
   Project is with Steven Stahler, 2015-2016
*/

#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <math.h>
#include <time.h>


extern double testCst;
extern double mPrevious;
extern double mCurrent;
// The code will assume that c_s = G = P_ambient = 1

// Constants (wrap in parentheses to avoid subtle arithmetic issues)
#define PI (4.0*atan(1))      // Delicious with whipped cream
#define Gcst (6.67408e-8)     // Gravitational constant, cgs

// Global Inputs
extern int M;               // grid poitns in radial direciton; command line input
extern int N;               // grid points in veritcal direction; command line input
extern double rRatio;       // ratio of vertical extent to entire domain's radial extent (not just cylinder radius; compiler input)
extern double mEx;          // Mass excess value from command line

// Cylinder knobs
extern double beta0;        // value of plasma beta at r=0
extern double nCyl;         // power index relating rho to B (B prop \rho^n); command line input
extern double lambda;       // mass/length value (e.g., hydro crit = 2*c^2/G); command line input
extern double rhoC;         // non-dim central density (must be > 1); command line input

// Global Derived Values
extern double zL;           // Vertical length of the box, input is in parsecs; command line input
extern double rL;           // Radial length of the box; derived
extern double DeltaR;       // cell size in radial direction; derived
extern double DeltaZ;       // cell size in vertical direction; derived

// (exact) Filament Boundary and Rho value
extern double *VContour;
extern double Rbdy;

// Rho at the top boundary
extern double *RhoTop;

// Total mass of cylinder without point masses
extern double totMass;

// Boundary Condition Arrays for values or derivatives
//extern double *VLeft;
extern double *VRight;
extern double *ALeft, *ARight;

// Global States
extern double ***newState;
extern double ***curState;
extern double ***prevState;

// Unit Conversions
extern double Pc2Code;
extern double Pc2Cm;
extern double Sol2Code;
extern double Sol2G;
extern double Cbackground;
extern double RhoTopCenter;

// Altered the boundary condition?
extern double Malter;


// Other preprocessors
#ifndef UNIFORM
   #define UNIFORM 0
#endif

#ifndef DEBUG
   #define DEBUG 0
#endif

#ifndef RADIALRATIO
    #define RADIALRATIO 1.0
#endif

#ifndef CYLINDERRADRAT
    #define CYLINDERRADRAT 1.0
#endif

// Enuemrate Declaration
#define NStates 5
enum states { Vpot, Apot, Q, Rho, dQdP };

// No'bdy gonna make me type std:: over and over again
using namespace std;

// Name of the output file
extern string OutFileName;

// Variables related to the cylinder Loop
extern int LoopMAX;
extern double CylTol;
extern double Ccst[2];

// Variables related to the point potential Loop
extern int PointLoopMax;
extern double PointLoopTol;

// Variables related to the convergence solve
extern int ConvergeLoopMax;
extern double ConvergeTol;
extern double relaxFrac;




/* --------------- Functions ------------------- */

inline void LongLine()
{
    cout << "-----------------------------------------------------" << endl;
}

inline void CodeHeader(string f)
{
    cout << endl;
    LongLine();
    cout << f << endl;
    LongLine();
};

inline void WaterlooHeader(string f)
{
    LongLine();
    cout << "                    WATERLOO                         " << endl;
    LongLine();
    cout << "In function: " << f << endl;
};

// Copying state 1 to state 2
inline void CopyState(double ***&State1, double ***&State2)
{
    int i,j,k;
    for(i=0;i<NStates;i++) for(j=0;j<M;j++) for(k=0;k<N;k++) State2[i][j][k] = State1[i][j][k];
};


// Plural function, because grammar matters
inline string Plural(int i,string s)
{
    if(i==1)
    {
        string blank = "";
        return blank;
    }
    else
    {
        return s;
    }
}

// Debugging
inline void Here()
{
    cout << "HERE!!" << endl;
}

// Index Generators
// Row major on interior cells : i = Ridx + Zidx * M
inline int Ridx(int i)
{
    return i%(M);
}

inline int Zidx(int i)
{
    return (i-Ridx(i))/(M);
}

inline int Sidx(int r, int z)
{
    return ( r + ( M * z ) );
};

// Calculates the position of the cell in one of the directions.
// Feed in the cell size and the cell index. Cell points are at edges (i=0 corresponds to r=0 or z=0)
inline double cPos(int i, double dx) {return ((double) i)*dx;}


// Input a Yarray and a particular X value, interpolates to find Y(X)
inline double LIntY(double* Yarray, double x, double Delta, int size, int& lastIDX)
{
    double value = 0;

    // Search for the indices 'x' lies between
    int idx=0,i;
    for(i=1;i<size;i++)
    {
        double curPos = cPos(i,Delta);
        if(x < curPos)
        {
            idx = i;
            lastIDX = idx;
            break;
        }
    }

    if(idx==0)
    {
        // Never found a point, throw error
        WaterlooHeader("LIntY: Failed to find point");
    }

    double y2 = Yarray[idx];
    double y1 = Yarray[idx-1];
    double x2 = cPos(idx,Delta);
    double x1 = cPos(idx-1,Delta);
    double m = (y2-y1)/(x2-x1);

    // Linearly interpolation
    value = y1 + m*(x-x1);

    return value;
};

inline void CalcRho()
{
    int i,j;
    for(i=0;i<M;i++) for(j=0;j<N;j++)
    {
        if( cPos(i,DeltaR) < VContour[j] )
        {
            curState[Rho][i][j] = curState[Q][i][j]*exp(-curState[Vpot][i][j]);
        }
        else
        {
            curState[Rho][i][j] = 0.0;
        }
    }
};

// Given a Yarray Y and a particular y value, finds the x such Y(x)=y
inline double LIntX(double* Yarray, double y, double Delta, int size)
{
    double value = 0;

    // Search for the indices 'y' lies between
    int idx=0,i;
    for(i=1;i<size;i++)
    {


        double curPos = cPos(i,Delta);
        if(y < curPos)
        {
            idx = i;
            break;
        }
    }

    if(idx==0)
    {
        // Never found a point, throw error
        WaterlooHeader("LIntX: Failed to find point");
    }

    double y2 = Yarray[i];
    double y1 = Yarray[i-1];
    double x2 = cPos(i,Delta);
    double x1 = cPos(i-1,Delta);
    double m = (y2-y1)/(x2-x1);

    // Linearly interpolation
    value = x1 + (y-y1)/m;

    return value;
};

// Prints out the data to a file
inline void PrintState(double ***curState,string divide, string num)
{
    int i,j,k;

    string OutFile;
    OutFile = OutFileName + divide + num;

    //OutFile = OutFileName.c_str() + num.c_str();

    ofstream myfile;
    myfile.open(OutFile.c_str());

    // First print inputs as the first line
    myfile << M << "," << N << "," << zL << "," << rRatio << ","
        << mEx << "," << beta0 << "," << nCyl << ","
        << Rbdy << "," << lambda << "," << Pc2Code << "," << Sol2Code << endl;

    // Print VContour as the second line
    myfile << VContour[0];
    for(j=1;j<N;j++) myfile << "," << VContour[j];
    myfile << endl;

    for(j=0; j< M; j++) for(k=0;k<N;k++)
    {
        double Rval = cPos(j,DeltaR);
        double Zval = cPos(k,DeltaZ);

        myfile << j << "," << k << "," << Rval << "," << Zval;
        for(i = 0; i< NStates; i++) myfile << "," << curState[i][j][k];
        myfile << endl;
    }
    myfile.close();
};


// Trap Integral
inline double TrapInt()
{
    double sum=0;

    // finds the largest index
    int idxMax = 0;
    int i,j;

    for(j=0;j<N;j++) for(i=0;i<M;i++)
    {
        if( cPos(i,DeltaR) >= VContour[j] )
        {
            if(i>idxMax) idxMax=i;
            break;
        }
    }

    idxMax = M-1;

    // Now do a trapezoidal integral out to idxMax
    sum = cPos(0,DeltaR)*curState[Rho][0][0] + cPos(idxMax,DeltaR)*curState[Rho][idxMax][0] + cPos(0,DeltaR)*curState[Rho][0][N-1] + cPos(idxMax,DeltaR)*curState[Rho][idxMax][N-1]; // corners
    for(i=1;i<idxMax;i++) sum = sum + 2.0*(cPos(i,DeltaR)*curState[Rho][i][0] + cPos(i,DeltaR)*curState[Rho][i][N-1]); // top and bottom edges
    for(j=1;j<N-1;j++) sum = sum + 2.0*(cPos(0,DeltaR)*curState[Rho][0][j] + cPos(idxMax,DeltaR)*curState[Rho][idxMax][j]); // left and right edges
    for(i=1;i<idxMax;i++) for(j=1;j<N-1;j++) sum = sum + 4.0*(cPos(i,DeltaR)*curState[Rho][i][j]); // interiors
    sum = 4.0*PI*(DeltaR*DeltaZ/4.0)*sum;

    return sum;
}

// 1D mass tradezoid integral (used for cylinder stuff)
inline double OneDMassInt(double* Rho, int Midx)
{
    double sum = 0;

    // Input rho, calculates  int( rho dV ) = 4*pi*zL int( r*rho dr ) using trapezoidal rule

    for(int i = 0; i < Midx-1; i++)
    {
        double fpo = cPos(i+1,DeltaR)*Rho[i+1];
        double f = cPos(i,DeltaR)*Rho[i];
        double LOCdeltaR = (cPos(i+1,DeltaR) - cPos(i,DeltaR));
        sum = sum + 0.5*LOCdeltaR*( fpo + f );
    }

    sum = 4.0*PI*zL*sum;

    return sum;
};

inline void CalcdQdP()
{
    //                      (dQ/dr)*DeltaR + (dQ/dz)*DeltaZ       Q(r+1)-Q(r-1)+Q(z+1)-Q(z-1)
    // Estimates dQ/dPhi =  ------------------------------- = ----------------------------------
    //                      (dP/dr)*DeltaR + (dP/dz)*DeltaZ   Phi(r+1)-Phi(r-1)+Phi(z+1)-Phi(z-1)
    //
    // Note: The above is true only if you perform a centered-difference for each derivative
    // (else factors of 1/2 and such would not cancel out entirely)

    // Boundary conditions:
    // top and bottom have reflective
    // left has reflective
    // right will utilize a three stencil estimate for the derivative


    // find Q at filament boundary
    int i=0; while(true){ i++; if(cPos(i,DeltaR)>=VContour[N-1]) break;}
    double Qright=curState[Q][i][N-1];
    double Qleft=curState[Q][i-1][N-1];
    double m = (Qright-Qleft)/DeltaR;
    double Qedge = Qleft + m*( VContour[N-1] - cPos(i,DeltaR) );

    int s;
    for(s=0;s<M*N;s++)
    {
        Qright = 0; Qleft = 0;
        double Qup=0, Qdown=0;
        double Phiright=0, Phileft=0, Phiup=0, Phidown=0;
        double dQdr=0,dQdz=0,dPdr=0,dPdz=0;




        if(Ridx(s)==0) // Left side
        {
            // We use one-sided derivatives to estimate these quantities.
            // We are dividing by DeltaR so that it cancels out the DeltaR at the end
            double Q0 = curState[Q][Ridx(s)][Zidx(s)];
            double Q1 = curState[Q][Ridx(s)+1][Zidx(s)];
            double Q2 = curState[Q][Ridx(s)+2][Zidx(s)];
            double Q3 = curState[Q][Ridx(s)+3][Zidx(s)];
            double P0 = cPos(Ridx(s),DeltaR)*curState[Apot][Ridx(s)][Zidx(s)]; // = 0
            double P1 = cPos(Ridx(s)+1,DeltaR)*curState[Apot][Ridx(s)+1][Zidx(s)];
            double P2 = cPos(Ridx(s)+2,DeltaR)*curState[Apot][Ridx(s)+2][Zidx(s)];
            double P3 = cPos(Ridx(s)+3,DeltaR)*curState[Apot][Ridx(s)+3][Zidx(s)];

            dQdr=(2*Q0-5*Q1+4*Q2-Q3)/DeltaR;
            dPdr=(2*P0-5*P1+4*P2-P3)/DeltaR;
        }
        else if(Ridx(s)==M-1) // Right side
        {
            // Instead of saying what Q(r+1) and Q(r-1) is, we'll
            // specify a value for Q(r+1) and set Q(r-1), where Q(r+1)
            // will estimate (dQ/dr)*DeltaR

            if(cPos(Ridx(s),DeltaR) > VContour[Zidx(s)]) // on right side of contour
            {
                dQdr = 0;
            }
            else
            {
                double Qr = curState[Q][Ridx(s)][Zidx(s)] ;
                double Qrm1 = curState[Q][Ridx(s)-1][Zidx(s)] ;
                double Qrm2 = curState[Q][Ridx(s)-2][Zidx(s)] ;
                dQdr = ( 3*Qr - 4*Qrm1 + 1*Qrm2 )/(2*DeltaR);
            }

            // Phi is done similarly, but dPhi/dr = d(r*A)/dr = A + r*(dA/dr)
            // Phiright will estimate = (dPhi/dr)*DeltaR = A*DeltaR + r*DeltaR*(dA/dr)
            double Ar = curState[Apot][Ridx(s)][Zidx(s)] ;
            double Arm1 = curState[Apot][Ridx(s)-1][Zidx(s)] ;
            double Arm2 = curState[Apot][Ridx(s)-2][Zidx(s)] ;
            double dAdr = (3*Ar - 4*Arm1 + 1*Arm2 )/(2*DeltaR); // continuous across filament, no special treatment needed

            dPdr = Ar + cPos(Ridx(s),DeltaR)*dAdr;

        }
        else // r derivatives -- interior
        {
            Qleft = curState[Q][Ridx(s)-1][Zidx(s)];
            Qright = curState[Q][Ridx(s)+1][Zidx(s)];
            dQdr = (Qright-Qleft)/(2*DeltaR); // this is true if both points are left of the contour

            Phileft  = cPos(Ridx(s)-1,DeltaR)*curState[Apot][Ridx(s)-1][Zidx(s)];
            Phiright = cPos(Ridx(s)+1,DeltaR)*curState[Apot][Ridx(s)+1][Zidx(s)];
            dPdr = (Phiright-Phileft)/(2*DeltaR);  // this is true if both points are left of the contour

            if(cPos(Ridx(s),DeltaR) >= VContour[Zidx(s)]) dQdr = 0; // cell is right of contour
            else if( cPos(Ridx(s)+1,DeltaR) >= VContour[Zidx(s)] )  // cell is near the contour
            {
                // Center difference overlaps the contour. Instead interpolate

                // Need the value of Q at the contour: We need the value of V at the contour
                //double Vr = curState[Vpot][1+Ridx(s)][Zidx(s)];
                //double Vl = curState[Vpot][Ridx(s)][Zidx(s)];
                double x  = (VContour[Zidx(s)] - cPos(Ridx(s),DeltaR))/(DeltaR);
                //double Vbdy = (1-x)*Vl + x*Vr;
                //double Qbdy = Rbdy * exp(Vbdy);

                double dr = (VContour[Zidx(s)] - cPos(Ridx(s)-1,DeltaR)) ;
                double dQ = Qedge - Qleft;
                dQdr = dQ/dr;



                // Do the same for Phi
                Phileft = cPos(Ridx(s)-1,DeltaR)*curState[Apot][Ridx(s)-1][Zidx(s)];
                double PhiC = cPos(Ridx(s),DeltaR)*curState[Apot][Ridx(s)][Zidx(s)];
                double Pbdy = (1-x)*PhiC + x*Phiright;
                dPdr = (Pbdy - Phileft)/(VContour[Zidx(s)] - cPos(Ridx(s)-1,DeltaR)) ;
            }
        }

        if(Zidx(s)==0 || Zidx(s)==N-1) // top or bottom
        {
            // Qup - Qdown = 0, similar for Phi.
            dQdz = 0;
            dPdz = 0;
        }
        else // interior
        {
            Qdown = curState[Q][Ridx(s)][Zidx(s)-1];
            Qup   = curState[Q][Ridx(s)][Zidx(s)+1];
            dQdz = (Qup-Qdown)/(2*DeltaZ); // Assumes both points are interior to contour

            Phidown  = cPos(Ridx(s),DeltaR)*curState[Apot][Ridx(s)][Zidx(s)-1];
            Phiup    = cPos(Ridx(s),DeltaR)*curState[Apot][Ridx(s)][Zidx(s)+1];
            dPdz = (Phiup-Phidown)/(2*DeltaZ); // Assumes both points are interior to contour

            if(cPos(Ridx(s),DeltaR) >= VContour[Zidx(s)]) dQdz = 0; // point is right of contour
            else if( cPos(Ridx(s),DeltaR) >= VContour[Zidx(s)+1] ) // point is next to contour
            {
                // Need z location of contour at cPos(Ridx(s),DeltaR) radius
                double m = DeltaZ/(VContour[Zidx(s)+1]-VContour[Zidx(s)]);
                double zCt = cPos(Zidx(s),DeltaZ) + m*(cPos(Ridx(s),DeltaR)-VContour[Zidx(s)]);

                // Now need to get Q at that location; to do so needs V at that location
                double Vu = curState[Vpot][Ridx(s)][Zidx(s)+1];
                double Vd = curState[Vpot][Ridx(s)][Zidx(s)];
                double x = (zCt-cPos(Zidx(s),DeltaZ))/DeltaZ;
                double Vbdy = (1-x)*Vd + x*Vu;
                double Qbdy = Rbdy * exp(Vbdy);

                double dQ = Qbdy - Qdown;
                double dz = (zCt - cPos(Zidx(s)-1,DeltaZ));// + DeltaZ;
                dQdz = dQ/dz;

                // Do the same for Phi
                Phidown  = cPos(Ridx(s),DeltaR)*curState[Apot][Ridx(s)][Zidx(s)];
                double Pbdy = (1-x)*Phidown + x*Phiup;
                dPdz = (Pbdy - Phidown)/(zCt - cPos(Zidx(s),DeltaZ));

            }


        }

        double dQ = dQdr*DeltaR + dQdz*DeltaZ; //Qup - Qdown + Qright - Qleft;
        double dP = dPdr*DeltaR + dPdz*DeltaZ; //Phiup - Phidown + Phiright - Phileft;

        double dQdPhi = 0;
        if(dP!=0) dQdPhi = dQ/dP;

        if(isnan(dQdPhi))
        {
            cout << "dQdPhi has a nan! (i,j)=(" << Ridx(s) << "," << Zidx(s) << ")" << endl;
            cout << "    dQ = " << dQ << " , dP = " << dP << endl;
            cout << "    dQdr = " << dQdr << " , dPdr = " << dPdr << endl;
            cout << "    dQdz = " << dQdz << " , dPdz = " << dPdz << endl;
            cout << "    cPos = " << cPos(Ridx(s),DeltaR) << " , VCont+1 = " << VContour[Zidx(s)+1] << endl;
            cout << "    cPos = " << cPos(Ridx(s),DeltaR) << " , VCont = " << VContour[Zidx(s)] << endl;
            if(cPos(Ridx(s),DeltaR) >= VContour[Zidx(s)+1]) cout << "Near a contour!" << endl;
            exit(1);
        }

        //if(isnan(dQdPhi)) dQdPhi = 0;

        //if(Zidx(s)==N-1) cout << dQdPhi << endl;

        if(cPos(Ridx(s),DeltaZ) <= VContour[Zidx(s)]) curState[dQdP][Ridx(s)][Zidx(s)] = dQdPhi;
        //if(dQdPhi > 0.3) cout << "Large derivative! at " << Ridx(s) << " and " << Zidx(s) << "! (" << dQdPhi << ")" << endl;
        else curState[dQdP][Ridx(s)][Zidx(s)] = 0; // dQdPhi should have been set to zero anyway, but just to be safe.

    }

};


// Given a new V and A, update Q
// This is done differently than how we initialized Q in PrepInit
// Q is constant along field lines, so we use the tow row (whose density never changes) to create
// a map between the updated Phi = r*A and Q = rho*exp(V) values.
// Then everywhere else in the box, we update Q appropriately.
// if type = 0, we update Q everywhere without caring where the boundary is
// if type = 1, we update Q by finding the value of Q right of the boundary and then setting all the right-of-boundary
// cells to this value (=rho)
inline void UpdateQ(int type)
{
    cout << "Updating Q (type=" << type << ")" << endl;

    int i,j;

    for(i=1;i<M-1;i++) { if(cPos(i,DeltaR)>VContour[N-1]) break;}
    double Rr = RhoTop[i];
    double Rl = RhoTop[i-1];
    double m2 = (Rr-Rl)/DeltaR;
    Rbdy = Rl + m2*(VContour[N-1]-cPos(i-1,DeltaR)); // updates global Rbdy
    cout << "New Rbdy value = " << Rbdy << endl; // shouldn't ever change...


    if(type==1) // we know where the new boundary is and Rho has been updated
    {
        double rEdge = VContour[N-1];

        // Interpolate Rho to rEdge
        for(i=1;i<M-1;i++) { if(cPos(i,DeltaR)>rEdge) break;}
        double Rr = RhoTop[i];
        double Rl = RhoTop[i-1];
        double m = (Rr-Rl)/DeltaR;
        Rbdy = Rl + m*(rEdge-cPos(i-1,DeltaR)); // updates global Rbdy
        cout << "New Rbdy value = " << Rbdy << endl;

        for(i=0;i<M;i++) for(j=0;j<N;j++)
        {
            if(cPos(i,DeltaR)>VContour[j]) curState[Q][i][j] = Rbdy;
        }

        return; // we don't have to make maps and stuff
    }

    // Make map
    double PhiMap[M];
    double QMap[M];
    for(i=0;i<M;i++) PhiMap[i] = cPos(i,DeltaR)*curState[Apot][i][N-1]; // Phi at top
    for(i=0;i<M;i++)   QMap[i] = RhoTop[i] * exp( curState[Vpot][i][N-1] ); // Q at top using RhoTop

    // What is Q at the filament boundary?
    i=0; while(true){ i++; if(cPos(i,DeltaR)>=VContour[N-1]) break; }
    double m = ( QMap[i] - QMap[i-1] )/DeltaR;
    double Qedge = QMap[i-1] + m*(VContour[N-1] - cPos(i-1,DeltaR));

    // Beyond the boundary, we use this Q value.
    for(i=0;i<M;i++) if(cPos(i,DeltaR)>=VContour[N-1]) QMap[i] = Qedge;

    // Top row is done
    for(i=0;i<M;i++) curState[Q][i][N-1] = QMap[i];

    // Now for every other cell, use map to interpolate a value of Q
    for(i=0;i<M;i++) for(j=0;j<N-1;j++) // skipping the top row
    {
        // If right of the contour, who cares
        if(cPos(i,DeltaR) >= VContour[j]){ curState[Q][i][j]=Qedge; continue; }

        // did we succeed?
        int yay=0;

        // Local phi value we want to use with PhiMap and QMap
        double localPhi = cPos(i,DeltaR)*curState[Apot][i][j];
        int idx;

        // Find the two indices localPhi lies between (if not monotonic, could be weird)
        for(idx=1;idx<M;idx++)
        {
            double Pleft = PhiMap[idx-1]; // cPos(idx-1,DeltaR)*curState[Apot][idx-1][j];
            double Pright= PhiMap[idx];   // cPos(idx,DeltaR)*curState[Apot][idx][j];

            if( (Pleft <= localPhi && localPhi <= Pright) || (Pleft >= localPhi && localPhi >= Pright) ) {yay=1; break;}
        }

        if(idx==M-1 && yay==0)
        {
            cout << "Possibly didn't find a Phi value to interpolate in updateQ (i,j)=(" << i << "," << j <<")" << endl;
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

        if(isnan(curState[Q][i][j]))
        {
            cout << "Q has a nan! (i,j) = (" << i << "," << j << ")" << endl;
            exit(1);
        }
    }

};
