#include "Filament_Globals.H"

void ConvergeNew();
void Converge();
void RelaxSoln();
int ConvergeTest(int i,int j);
int UpdateRho(int type);
void NewUpdateRho(int t);
double FindContRatio(double topRad, int type);

void findDenCont(double *a, double r, int t);
double findMass(double *a);
void NewRhoBdy();
void NewBdyCond();

double getMaxDV(double **dv,int debug);

void SingleUpdate();
