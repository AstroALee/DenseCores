/* Dense Core */
#include "Filament_Globals.H"
#include "Filament_Main.H"



void PrepareInitialState();
void Converge();


int main(int argc, char *argv[])
{
    // Start the clock
    clock_t startTime = clock();

    // In case we need random numbers
    srand(time(NULL));

    // Check to make sure your arguments are good
    // If so, read them in
    CheckArguments(argc);
    ReadArguments(argv);

    // Welcome message
    WelcomeMessage();

    // Print arguments to screen
    PrintArguments();

    // Compute various derived quantities
    // Allocates remaining quantities
    CalcDerived();

    // Allocate State Arrays
    AllocateState(curState);
    AllocateState(prevState);

    // Prepare initial conditions
    PrepareInitialState();

    // Copy initial state
    CopyState(curState,prevState);

    // Enter the Solving Routine
    CodeHeader("Entering Convergence Loop");
    Converge();

    // Print the final state
    PrintState(curState);

    // Deallocate State Arrays
    DeallocateState(curState);
    DeallocateState(prevState);

    // Stop the Clock and Final Tallies
    CodeHeader("Timing and Tallies");
    Tallies(startTime);

    // See ya!
    CodeHeader("Goodbye!");
    return 0;
}
