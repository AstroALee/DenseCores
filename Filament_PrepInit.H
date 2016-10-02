#include "Filament_Globals.H"

double InitVB(double z) { return ( VContour[0]*(1-(1-0.8)*(z/zL)) ); }; // 0.8 = ratio of top to bottom radius
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
void getVbdy();
void adjustRho();
void truncateRho();
void InitQ();
void UnitConvert();
void BaselineCylinder(int i);
void firstRbdy();
void stretchRho();
void StretchRho_old();
void BallDensity(double m);
void BallExcess(double m);
double reduceCyl(double x);
double MassInt(double A);
void NewBall(double m);
void SoftPotential();
