/* Global Variable Declarations for the Filament Code by Aaron T. Lee
   Project is with Steven Stahler, 2015-2016
*/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// The code will assume that c_s = G = rho_ambient = 1

// Constants
#define PI 4.0*atan(1)      // Delicious with whipped cream
#define Gcst 6.67408e-8     // Gravitational constant, cgs

// Global Inputs
extern int M;               // grid poitns in radial direciton
extern int N;               // grid points in veritcal direction
extern double rRatio;       // ratio of vertical extent to entire domain's radial extent (not just cylinder radius)
extern double mExcess;      // mass excess beyond cylinder's mass

// Cylinder knobs
extern double betaInf;      // value of plasma beta in the ambient region
extern double nCyl;         // power index relating rho to B (B prop \rho^n)

// Global Derived Values
extern double zL;           // Vertical length of the box, input is in parsecs
extern double rL;           // Radial length of the box
extern double DeltaR;       // cell size in radial direction
extern double DeltaZ;       // cell size in vertical direction

// (exact) Filament Boundary and V value
extern double *VContour;
extern double Vbdy;

// Ambient number density
extern double nBackground;

// Rho at the top boundary
extern double *RhoTop;

// Boundary Condition Arrays for values or derivatives
extern double *VLeft, *VRight;
extern double *ALeft, *ARight;

// Global States
extern double ***curState;
extern double ***prevState;

// Unit Conversions
extern double Pc2Code;
extern double Sol2Code;

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



// Enuemrate Declaration
#define NStates 4
enum states { Vpot, Apot, Q, Rho };

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

// Copying state
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
inline double LIntY(double* Yarray, double x, double Delta, int size)
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
            break;
        }
    }

    if(idx==0)
    {
        // Never found a point, throw error
        WaterlooHeader("LIntY: Failed to find point");
    }

    double y2 = Yarray[i];
    double y1 = Yarray[i-1];
    double x2 = cPos(i,Delta);
    double x1 = cPos(i-1,Delta);
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
