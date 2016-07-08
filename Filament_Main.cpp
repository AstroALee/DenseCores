/* Dense Core */
#include "Filament_Globals.H"
#include "Filament_Main.H"



void PrepareInitialState();
void Converge();


int main(int argc, char *argv[])
{
    // Set precision (give me all the numbers)
    cout << std::setprecision(16);

    // Start the clock
    clock_t startTime = clock();

    // In case we need random numbers
    srand(time(NULL));

    // Check to make sure your arguments are good
    // If so, read them in
    CheckArguments(argc);
    ReadArguments(argc,argv);

    // Welcome message
    WelcomeMessage();

    // Print arguments to screen
    PrintArguments(argc);

    // Compute various derived quantities
    // Allocates remaining quantities
    CalcDerived(argc);

    // Allocate State Arrays
    AllocateState(newState);
    AllocateState(curState);
    AllocateState(prevState);

    // Prepare initial conditions
    PrepareInitialState();

    // Copy initial state (I think this is redundant with what's in SolveAll.cpp)
    //CopyState(curState,prevState);

    // Enter the Solving Routine
    CodeHeader("Convergence Loop");
    Converge();

    // Calculate the total mass
    //CodeHeader("Calculating Final Masses");
    //CalcMass();

    // Print the final state
    PrintState(curState,".","out");

    // Deallocate State Arrays
    DeallocateState(newState);
    DeallocateState(curState);
    DeallocateState(prevState);

    // Stop the Clock and Final Tallies
    CodeHeader("Timing and Tallies");
    Tallies(startTime);

    // See ya!
    CodeHeader("Goodbye!");
    return 0;
}
