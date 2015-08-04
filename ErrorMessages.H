/* Error Messages and Info Messages */


inline void SHOWME_Inputs(char** argv)
{
    cout << "You had N = " << argv[1] << endl;
    cout << "You had Z = " << argv[2] << endl;
    cout << "You had pLength = " << argv[3] << endl;
    cout << "You had rhoCenter = " << argv[4] << endl;
    cout << "You had MassCore = " << argv[5] << endl;
    cout << "You had beta = " << argv[6] << endl;
    cout << endl << endl;
  
}


inline void WATERLOO_ArgCountCheck(int i)
{
    if(i < 7)
    {
        cout << "WATERLOO: Not enough arguments (N,Z,Non-Dim Width,RhoCenter,MassGuess,betaedge)" << endl;
        cout << "You provided " << i-1 << " argument" << Plural(i-1,"s") << "." << endl;
        exit(1);
    }
    if(i > 7)
    {
        cout << "WATERLOO: Too many arguements? (N,Z,Non-Dim Width,RhoCenter,MassGuess,betaedge) " << endl;
        cout << "You provided " << i-1 << " argument" << Plural(i-1,"s") << "." << endl;
        exit(1);
    }
    
    
}

inline void WATERLOO_InitPotNotConverging(double i, double tol)
{
    cout << "WATERLOO: Calculating the initial gravitational potential from a "
         << "set of point masses is not converging. Failed to meet a tolerance "
         << "of " << tol << " after " << i << " iterations. Sad face." << endl;
    exit(1);
}


inline void WATERLOO_MagCylNoConverge(int i)
{
    cout << "Magnetized Cylinder integrations failed to converge after " << i << " attempts.";
    cout << endl;
    exit(1);
}