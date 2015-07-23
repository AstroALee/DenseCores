#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* Enuemrate Declaration */
#define NStates 3

// Constants (cgs)
#define PI 3.14159265359
#define GCST 6.67259e-8;

// Global variables, so they can be global across different files
extern int N;
extern int Z;
extern double RhoCorner;
extern double RhoCenter;
extern double MassGuess;
extern double betaedge;

extern double DeltaR;
extern double DeltaZ;

extern double* VContour;
extern double RightContour;



// Enumeration Definitions (should be the same number as NStates above)
enum states { Vpot, Apot, Q };

// No'bdy gonna make me type std:: over and over again
using namespace std;



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
