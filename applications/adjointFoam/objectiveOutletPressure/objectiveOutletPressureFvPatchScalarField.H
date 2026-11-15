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
    This file is part of of14-adjoint. The class layout is the Foundation's
    adjointOutletPressure; the boundary value follows ESI-OpenCFD's
    adjointOutletPressure (src/optimisation/adjointOptimisation/adjoint/
    adjointBoundaryConditions), with the sign of the objective source flipped
    for this solver's convention (Ua = -u, pa = -q of Othmer). See NOTICE.md.

    It is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version. It is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details: <http://www.gnu.org/licenses/>.

Class
    Foam::objectiveOutletPressureFvPatchScalarField

Description
    Adjoint pressure at an outlet for an objective that depends on the
    outlet's flow rate (the flow-rate partition objective of objective.H):

        pa = Ua_n U_n + nuEff d(Ua_n)/dn - source  [+ Ua.U with the UaGradU ATC]

    with U_n = phi/|Sf|, Ua_n = phia/|Sf| and source = dJ/dQ for this outlet,
    pushed by the solver every iteration (objectiveUpdate.H), as are nuEff,
    whether the viscous term is kept and whether the Ua.U term belongs (the
    UaGradU form of the adjoint transpose convection only).

    Usage:
        outlet1 { type objectiveOutletPressure; value uniform 0; }

SourceFiles
    objectiveOutletPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveOutletPressureFvPatchScalarField_H
#define objectiveOutletPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objectiveOutletPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- dJ/dQ of this outlet (uniform), set by the solver
        scalar source_;

        //- Add Ua.U (the UaGradU ATC form)
        bool addUaU_;

        //- Add nuEff d(Ua_n)/dn
        bool viscous_;

        //- Effective viscosity on the patch, set by the solver
        scalarField nuEff_;


public:

    //- Runtime type information
    TypeName("objectiveOutletPressure");


    // Constructors

        objectiveOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const dictionary&
        );

        objectiveOutletPressureFvPatchScalarField
        (
            const objectiveOutletPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const fieldMapper&
        );

        objectiveOutletPressureFvPatchScalarField
        (
            const objectiveOutletPressureFvPatchScalarField&,
            const DimensionedField<scalar, fvMesh>&
        );

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new objectiveOutletPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Setters (the solver pushes these every iteration)

            void source(const scalar s)
            {
                source_ = s;
            }

            void addUaU(const bool b)
            {
                addUaU_ = b;
            }

            void viscous(const bool b)
            {
                viscous_ = b;
            }

            void nuEff(const scalarField& nu)
            {
                nuEff_ = nu;
            }

            scalar source() const
            {
                return source_;
            }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        using fixedValueFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
