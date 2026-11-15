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

#include "objectiveOutletPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectiveOutletPressureFvPatchScalarField::
objectiveOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    source_(dict.lookupOrDefault<scalar>("source", 0)),
    addUaU_(true),
    viscous_(false),
    nuEff_(p.size(), 0)
{}


Foam::objectiveOutletPressureFvPatchScalarField::
objectiveOutletPressureFvPatchScalarField
(
    const objectiveOutletPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    source_(ptf.source_),
    addUaU_(ptf.addUaU_),
    viscous_(false),
    nuEff_(p.size(), 0)
{}


Foam::objectiveOutletPressureFvPatchScalarField::
objectiveOutletPressureFvPatchScalarField
(
    const objectiveOutletPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    source_(tppsf.source_),
    addUaU_(tppsf.addUaU_),
    viscous_(tppsf.viscous_),
    nuEff_(tppsf.nuEff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::objectiveOutletPressureFvPatchScalarField::updateCoeffs()
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

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("Ua");

    // ESI: pa = Ua_n U_n + nuEff snGrad(Ua).n + source (+ Ua.U for UaGradU);
    // here the objective source enters with the opposite sign (pa = -q).
    scalarField value((phiap/magSf)*(phip/magSf) - source_);

    if (addUaU_)
    {
        value += Up & Uap;
    }

    if (viscous_ && nuEff_.size() == patch().size())
    {
        value += nuEff_*(Uap.snGrad() & nf);
    }

    operator==(value);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::objectiveOutletPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "source", source_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        objectiveOutletPressureFvPatchScalarField
    );
}

// ************************************************************************* //
