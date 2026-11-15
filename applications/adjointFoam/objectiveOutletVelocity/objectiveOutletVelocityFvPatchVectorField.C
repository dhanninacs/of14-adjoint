/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2023 OpenCFD Ltd.
    Adapted for OpenFOAM 14 (Foundation) by of14-adjoint, 2026.
-------------------------------------------------------------------------------
License
    This file is part of of14-adjoint (GPL-3.0-or-later); see the header file
    and NOTICE.md for its provenance.

\*---------------------------------------------------------------------------*/

#include "objectiveOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectiveOutletVelocityFvPatchVectorField::
objectiveOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    tangentialRobin_(true),
    nuEff_(p.size(), 0)
{}


Foam::objectiveOutletVelocityFvPatchVectorField::
objectiveOutletVelocityFvPatchVectorField
(
    const objectiveOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    tangentialRobin_(ptf.tangentialRobin_),
    nuEff_(p.size(), 0)
{}


Foam::objectiveOutletVelocityFvPatchVectorField::
objectiveOutletVelocityFvPatchVectorField
(
    const objectiveOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, fvMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    tangentialRobin_(pivpvf.tangentialRobin_),
    nuEff_(pivpvf.nuEff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::objectiveOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& magSf = patch().magSf();
    const vectorField nf(patch().nf());

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phia");

    vectorField value(nf*phiap/magSf);

    if (tangentialRobin_ && nuEff_.size() == patch().size())
    {
        // ESI: Ua_t = (nd Ua_t,c - explicit deviatoric flux - source)/(U_n + nd);
        // the flow-rate objective has no tangential source, the explicit
        // deviatoric flux is left out, and U_n >= 0 keeps the denominator
        // positive (backflow is refused upstream).
        const scalarField Un(max(phip/magSf, scalar(0)));
        const scalarField nd(nuEff_*patch().deltaCoeffs());
        const vectorField Uac(patchInternalField());
        value += (Uac - nf*(nf & Uac))*nd/max(Un + nd, vSmall);
    }

    vectorField::operator=(value);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::objectiveOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        objectiveOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
