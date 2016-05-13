#include "Filament_Globals.H"

// NR routines
#include "NRRoutines/nr3.h"
#include "NRRoutines/stepper.h"
#include "NRRoutines/odeint.h"
#include "NRRoutines/stepperdopr5.h"
//#include "stepperross.h"



void UseOutput(Output out,double rEdge,double& curMass);


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
        // hi

        double m = (RightY-LeftY)/(RightX-LeftX);

        value = LeftY + m*(point-LeftX);
        found = 1;
        break;
    }

    if(!found) cout << "Interpolated value was not found!" << endl;

    return value;
};
