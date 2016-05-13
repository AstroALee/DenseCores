/* Global Variable Declarations for the Filament Code by Aaron T. Lee
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

// The code will assume that c_s = G = P_ambient = 1

// Constants (wrap in parentheses to avoid subtle arithmetic issues)
#define PI (4.0*atan(1))      // Delicious with whipped cream
#define Gcst (6.67408e-8)     // Gravitational constant, cgs

// Global Inputs
extern int M;               // grid poitns in radial direciton; command line input
extern int N;               // grid points in veritcal direction; command line input
extern double rRatio;       // ratio of vertical extent to entire domain's radial extent (not just cylinder radius; compiler input)
extern double mExcess;      // mass excess beyond cylinder's mass; command line input

// Cylinder knobs
extern double betaInf;      // value of plasma beta in the ambient region; command line input
extern double nCyl;         // power index relating rho to B (B prop \rho^n); command line input
extern double lambda;       // mass/length value (e.g., hydro crit = 2*c^2/G); command line input

// Global Derived Values
extern double zL;           // Vertical length of the box, input is in parsecs; command line input
extern double rL;           // Radial length of the box; derived
extern double DeltaR;       // cell size in radial direction; derived
extern double DeltaZ;       // cell size in vertical direction; derived

// (exact) Filament Boundary and V value
extern double *VContour;
extern double Vbdy;

// Rho at the top boundary
extern double *RhoTop;

// Total mass of cylinder without point masses
extern double totMass;

// Boundary Condition Arrays for values or derivatives
extern double *VLeft, *VRight;
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
extern double Pbackground;


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
        << mExcess << "," << betaInf << "," << nCyl << ","
        << Vbdy << "," << lambda << "," << Pc2Code << "," << Sol2Code << endl;

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

    // Boundary conditions:
    // top and bottom have reflective
    // left has reflective
    // right will utilize a three stencil estimate for the derivative

    int s;
    for(s=0;s<M*N;s++)
    {
        double Qright=0, Qleft=0, Qup=0, Qdown=0;
        double Phiright=0, Phileft=0, Phiup=0, Phidown=0;

        if(Ridx(s)==0) // Left side
        {
            // Qright-Qleft = 0, similar for Phi
        }
        else if(Ridx(s)==M-1) // Right side
        {
            // Instead of saying what Q(r+1) and Q(r-1) is, we'll
            // specify a value for Q(r+1) and set Q(r-1), where Q(r+1)
            // will estimate (dQ/dr)*DeltaR

            double Qr = curState[Q][Ridx(s)][Zidx(s)] ;
            double Qrm1 = curState[Q][Ridx(s)-1][Zidx(s)] ;
            double Qrm2 = curState[Q][Ridx(s)-2][Zidx(s)] ;
            double dQdr = 0.5*( 3*Qr - 4*Qrm1 + 1*Qrm2 );

            Qright = dQdr;

            // Phi is done similarly, but dPhi/dr = d(r*A)/dr = A + r*(dA/dr)
            // Phiright will estimate = (dPhi/dr)*DeltaR = A*DeltaR + r*DeltaR*(dA/dr)
            double Ar = curState[Apot][Ridx(s)][Zidx(s)] ;
            double Arm1 = curState[Apot][Ridx(s)-1][Zidx(s)] ;
            double Arm2 = curState[Apot][Ridx(s)-2][Zidx(s)] ;
            double dAdr = 0.5*(3*Ar - 4*Arm1 + 1*Arm2 );

            Phiright = DeltaR*Ar + cPos(Ridx(s),DeltaR)*dAdr;

        }
        else // interior
        {
            Qleft = curState[Q][Ridx(s)-1][Zidx(s)];
            Qright = curState[Q][Ridx(s)+1][Zidx(s)];

            Phileft  = cPos(Ridx(s)-1,DeltaR)*curState[Apot][Ridx(s)-1][Zidx(s)];
            Phiright = cPos(Ridx(s)+1,DeltaR)*curState[Apot][Ridx(s)+1][Zidx(s)];
        }

        if(Zidx(s)==0 || Zidx(s)==N-1) // top or bottom
        {
            // Qup - Qdown = 0, similar for Phi.
        }
        else // interior
        {
            Qdown = curState[Q][Ridx(s)][Zidx(s)-1];
            Qup   = curState[Q][Ridx(s)][Zidx(s)+1];

            Phidown  = cPos(Ridx(s),DeltaR)*curState[Apot][Ridx(s)][Zidx(s)-1];
            Phiup    = cPos(Ridx(s),DeltaR)*curState[Apot][Ridx(s)][Zidx(s)+1];
        }

        double dQ = Qup - Qdown + Qright - Qleft;
        double dP = Phiup - Phidown + Phiright - Phileft;

        double dQdPhi = 0;
        if(dP!=0) dQdPhi = dQ/dP;

        curState[dQdP][Ridx(s)][Zidx(s)] = dQdPhi;

    }

};
