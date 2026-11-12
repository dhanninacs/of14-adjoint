/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Adapted for OpenFOAM 14 (Foundation) by of14-adjoint, 2026.
-------------------------------------------------------------------------------
License
    This file is part of of14-adjoint, derived from OpenFOAM (ESI-OpenCFD,
    src/optimisation/adjointOptimisation/adjoint/ATCModel: ATCModel,
    ATCstandard, ATCUaGradU, cancelATC, zeroATCcells/faceCells).

    It is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version. It is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details: <http://www.gnu.org/licenses/>.

Description
    The adjoint transpose convection (ATC) term of the adjoint momentum
    equation, with ESI's choice of form and limiter, as one header for the
    solver (ESI's run-time selection classes collapsed into a switch):

        ATC
        {
            type                UaGradU;   // UaGradU | standard | cancel
            zeroATCPatchTypes   ();        // e.g. (wall patch): zero ATC next to these
            nSmooth             0;         // smoothing passes of the limiter
            extraConvection     0;         // implicit adjoint convection added
        }

    - UaGradU:  ATC = -(grad(Ua) & U), the Foundation solver's form (default)
    - standard: ATC =  (grad(U) & Ua), ESI's default form
    - cancel:   no ATC term (robust, less accurate)

    The term is added to the left-hand side of the adjoint momentum equation,
    multiplied by the limiter: 1 everywhere except the cells next to the given
    patch types and the inlets, where it is 0, then smoothed nSmooth times by
    a local-minimum interpolation and averaging (ESI's computeLimiter).
    extraConvection > 0 adds that multiple of the adjoint convection
    implicitly, for diagonal dominance, and corrects the right-hand side by
    the same amount (ESI's extraConvection).

\*---------------------------------------------------------------------------*/

// The settings, read once.
const dictionary ATCDict
(
    adjointDict.found("ATC") ? adjointDict.subDict("ATC") : dictionary()
);
const word ATCType(ATCDict.lookupOrDefault<word>("type", "UaGradU"));
if (ATCType != "UaGradU" && ATCType != "standard" && ATCType != "cancel")
{
    FatalIOErrorInFunction(ATCDict)
        << "Unknown ATC type " << ATCType
        << "; use UaGradU, standard or cancel" << exit(FatalIOError);
}
const wordList zeroATCPatchTypes
(
    ATCDict.lookupOrDefault<wordList>("zeroATCPatchTypes", wordList())
);
const label ATCnSmooth(ATCDict.lookupOrDefault<label>("nSmooth", 0));
const scalar ATCextraConvection
(
    ATCDict.lookupOrDefault<scalar>("extraConvection", 0)
);
Info<< "ATC: " << ATCType << ", zeroed next to " << zeroATCPatchTypes
    << " and the inlets, " << ATCnSmooth << " smoothing passes, extra convection "
    << ATCextraConvection << nl << endl;

// The limiter (ESI ATCModel::computeLimiter with zeroATCcells::faceCells).
volScalarField ATClimiter
(
    IOobject
    (
        "ATClimiter",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 1),
    zeroGradientFvPatchScalarField::typeName
);
{
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if (findIndex(zeroATCPatchTypes, word(patch.type())) != -1)
        {
            for (const label celli : patch.faceCells())
            {
                ATClimiter[celli] = 0;
            }
        }
    }
    // The Foundation solver zeroed the term on the inlet cells: kept.
    for (const label celli : inletCells)
    {
        ATClimiter[celli] = 0;
    }
    ATClimiter.correctBoundaryConditions();
    const localMin<scalar> limiterScheme(mesh);
    for (label pass = 0; pass < ATCnSmooth; ++pass)
    {
        surfaceScalarField faceLimiter(limiterScheme.interpolate(ATClimiter));
        ATClimiter = fvc::average(faceLimiter);
        ATClimiter.correctBoundaryConditions();
    }
}
