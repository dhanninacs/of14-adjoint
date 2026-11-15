// of14-adjoint: the objective and the wall shape sensitivity.
//
// J = -sum over the non-wall, non-constraint patches of (p + |U|^2/2) (U . Sf): the total-pressure
// power lost between the inlets and outlets (kinematic, m^5/s^3), the objective the adjoint
// boundary conditions encode.
//
// shapeSensitivity, on each wall face: s = nuEff (dU/dn) . (dUa/dn), with n the face normal out of
// the fluid and the normal derivatives OpenFOAM's own boundary snGrad by default
// (sensitivityWallGradient: corrected and leastSquares below). It is the change of J per
// unit area and per unit displacement of the wall along n (Othmer 2008, frozen turbulence). The
// sign was fixed against finite differences of J (Stokes ADR-0038): moving a wall out of the fluid
// lowers the loss, so s < 0 there.

//
// With the flowRatePartition objective (objective.H), J is the flow-split J of the last
// objectiveUpdate.H (the same phi: the primal is not solved again in between), and objective.dat
// also carries the outlet fractions and dJ/dQ_i. The wall formula is the same for both.

{
    scalar J = 0;
    if (flowRatePartition)
    {
        J = objectiveJ;
        if (objectiveFile.valid())
        {
            objectiveFile() << runTime.userTimeName() << tab << J;
            forAll(objectiveFractions, i)
            {
                objectiveFile() << tab << objectiveFractions[i];
            }
            forAll(objectiveDJdQ, i)
            {
                objectiveFile() << tab << objectiveDJdQ[i];
            }
            objectiveFile() << endl;
        }
    }
    else
    {
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
    }
    Info<< "Objective J = " << J << nl;

    // Updated every iteration (walls only), so every write carries it, the final write of a
    // run that stops on residualControl included.
    {
        const volScalarField nuEff(turbulence->nuEff());
        volScalarField::Boundary& sb = shapeSensitivity.boundaryFieldRef();
        if (sensitivityLeastSquares)
        {
            // As "corrected" below, with the wall cell's gradient from least squares over its
            // face neighbours and its boundary faces at their true centres (wallCellGrad.H):
            // exact for a linear field on any cell shape, where the Gauss gradient of a snapped
            // wall cell is not (tests/flowsplit.md, Phase C).
            forAll(mesh.boundary(), patchi)
            {
                const fvPatch& patch = mesh.boundary()[patchi];
                if (isA<wallFvPatch>(patch))
                {
                    const vectorField nf(patch.nf());
                    const vectorField d(patch.Cf() - patch.Cn());
                    const scalarField nd(max(nf & d, small*mag(d)));
                    const vectorField k(d - nf*(nf & d));
                    const fvPatchVectorField& Up = U.boundaryField()[patchi];
                    const fvPatchVectorField& Uap = Ua.boundaryField()[patchi];
                    const vectorField snU
                    (
                        (Up - Up.patchInternalField() - (k & wallCellGrad(U, patchi)))/nd
                    );
                    const vectorField snUa
                    (
                        (Uap - Uap.patchInternalField() - (k & wallCellGrad(Ua, patchi)))/nd
                    );
                    sb[patchi] = nuEff.boundaryField()[patchi]*(snU & snUa);
                }
            }
        }
        else if (!sensitivityCorrected)
        {
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
        else
        {
            // The wall gradient along the face normal. The boundary snGrad is
            // (X_w - X_c)/(n & d) with d from the cell centre to the face centre: on a
            // non-orthogonal wall cell (snapped, or a body-fitted slope) the centre lies off
            // the normal by k = d - n (n & d), so the difference is not taken along n. Here
            // the cell value is moved to the foot of the normal with the cell gradient first
            // (adjointControls.H).
            const volTensorField gradU(fvc::grad(U));
            const volTensorField gradUa(fvc::grad(Ua));
            forAll(mesh.boundary(), patchi)
            {
                const fvPatch& patch = mesh.boundary()[patchi];
                if (isA<wallFvPatch>(patch))
                {
                    const vectorField nf(patch.nf());
                    // Not patch.delta(): Foundation's is already the normal part,
                    // nf (nf & (Cf - Cn)), so boundary snGrad divides by n & d and
                    // ignores the tangential offset k.
                    const vectorField d(patch.Cf() - patch.Cn());
                    const scalarField nd(max(nf & d, small*mag(d)));
                    const vectorField k(d - nf*(nf & d));
                    const fvPatchVectorField& Up = U.boundaryField()[patchi];
                    const fvPatchVectorField& Uap = Ua.boundaryField()[patchi];
                    const vectorField snU
                    (
                        (
                            Up - Up.patchInternalField()
                          - (k & gradU.boundaryField()[patchi].patchInternalField())
                        )/nd
                    );
                    const vectorField snUa
                    (
                        (
                            Uap - Uap.patchInternalField()
                          - (k & gradUa.boundaryField()[patchi].patchInternalField())
                        )/nd
                    );
                    sb[patchi] = nuEff.boundaryField()[patchi]*(snU & snUa);
                }
            }
        }
    }
}
