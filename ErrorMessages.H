/* Error Messages */



inline void WATERLOO_ArgCountCheck(int i)
{
    if(i < 7)
    {
        cout << "WATERLOO: Not enough arguments (N,Z,RhoCorner,RhoCenter,MassGuess,betaedge)" << endl;
        cout << "You provided " << i-1 << " argument" << Plural(i-1,"s") << "." << endl;
        exit(1);
    }
    if(i > 7)
    {
        cout << "WATERLOO: Too many arguements? (N,Z,RhoCorner,RhoCenter,MassGuess,betaedge) " << endl;
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