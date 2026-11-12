// of14-adjoint: the objective and the wall shape sensitivity.
//
// J = -sum over the non-wall, non-constraint patches of (p + |U|^2/2) (U . Sf): the total-pressure
// power lost between the inlets and outlets (kinematic, m^5/s^3), the objective the adjoint
// boundary conditions encode.
//
// shapeSensitivity, on each wall face: s = nuEff (dU/dn) . (dUa/dn), with n the face normal out of
// the fluid and the normal derivatives OpenFOAM's own boundary snGrad. It is the change of J per
// unit area and per unit displacement of the wall along n (Othmer 2008, frozen turbulence). The
// sign was fixed against finite differences of J (Stokes ADR-0038): moving a wall out of the fluid
// lowers the loss, so s < 0 there.

{
    scalar J = 0;
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if (isA<wallFvPatch>(patch) || patch.constraint())
        {
            continue;
        }
        const vectorField& Up = U.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        J -= sum((pp + 0.5*magSqr(Up))*(Up & patch.Sf()));
    }
    reduce(J, sumOp<scalar>());
    if (objectiveFile.valid())
    {
        objectiveFile() << runTime.userTimeName() << tab << J << endl;
    }
    Info<< "Objective J = " << J << nl;

    // Updated every iteration (walls only), so every write carries it, the final write of a
    // run that stops on residualControl included.
    {
        const volScalarField nuEff(turbulence->nuEff());
        volScalarField::Boundary& sb = shapeSensitivity.boundaryFieldRef();
        forAll(mesh.boundary(), patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                sb[patchi] =
                    nuEff.boundaryField()[patchi]
                   *(U.boundaryField()[patchi].snGrad() & Ua.boundaryField()[patchi].snGrad());
            }
        }
    }
}
