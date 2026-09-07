// Momentum predictor

tmp<fvVectorMatrix> tUEqn // expressed in terms of superficial velocity, following "https://github.com/csoulain/porousMedia4Foam/blob/master/solvers/dbsFoam/UEqn.H"
(
    
    1./sqr(porosity+1.0e-8)*fvm::div(phi, U)
    - fvm::laplacian(nu/(porosity+1.0e-8), U) // was "//+ turbulence->divDevReff(U)"
    +  fvm::Sp(nu*invPermeability, U) // Darcy drag term
    + fvm::Sp(forchCoeff*mag(U), U) // Forchheimer drag term (observed not to do much)

    //removed "== fvOptions(U)", as it is not used
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);
