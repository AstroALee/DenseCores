#include "Filament_Globals.H"

void Converge();
void RelaxSoln();
void UpdateBdy();
void NewVBdy();
int UpdateRho(int i);
void UpdateBC(int i);

void ConvergeTest(int l, int& done);
void PrintError(int i, double errV, double errA);
