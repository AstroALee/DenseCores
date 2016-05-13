#include "Filament_Globals.H"

void Converge();
void RelaxSoln();
void NewVBdy();
void UpdateQ();
void UpdateRho();
void CalcdQdP();

void ConvergeTest(int l, int& done);
void PrintError(int i, double errV, double errA);
