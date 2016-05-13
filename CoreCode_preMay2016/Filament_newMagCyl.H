#include "Filament_Globals.H"

// NR routines
#include "NRRoutines/nr3.h"
#include "NRRoutines/stepper.h"
#include "NRRoutines/odeint.h"
#include "NRRoutines/stepperdopr5.h"
//#include "stepperross.h"



void UseOutput(Output out, double rEdge);


// Linearly interpolate on the output to find the particular r value where it occurs
inline double LIntY(VecDoub Xarray, MatDoub Yarray, double point, int size, int elem)
{
    double value = 0;
    int found = 0;
    // Find the two points surrounding 'point' in Yarray[elem][:]
    for(int i=0;i<size-1;i++)
    {
        if( (Yarray[elem][i] <= point && Yarray[elem][i+1] >= point) || (Yarray[elem][i] >= point && Yarray[elem][i+1] <= point)  )
        {
            found = i;

            double y1 = Yarray[elem][i];
            double y2 = Yarray[elem][i+1];
            double x1 = Xarray[i];
            double x2 = Xarray[i+1];
            double m = (y2-y1)/(x2-x1);

            value = x1 + (point-y1)/m;
        }
    }

    if(!found) cout << "Interpolated value was not found!" << endl;
    return value;
};


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
