#include "Filament_Globals.H"

double InitVB(double z) { return ( VContour[0]*(1-(1-contR)*(z/zL)) ); };
double WFunc(double A,double r, double z)
{
    double WF = A*(1-pow(z/zL,4))*(1-pow(r/InitVB(z),4));

    return WF;
};

void PrepareInitialState();
void PrepareInitialState_old();
void InitCylinder();
void InitPoints(double m);
void detDVDR(double m);
void getRbdy();
void adjustRho();
void InitQ();
void UnitConvert();
void BaselineCylinder(int i);
void firstRbdy();
void StretchRho();
void BallDensity(double m);
double reduceCyl(double x);
double MassInt(double A);
void NewBall();
