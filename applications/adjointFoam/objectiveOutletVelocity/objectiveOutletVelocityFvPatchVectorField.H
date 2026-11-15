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
    adjointOutletVelocity; the boundary value follows ESI-OpenCFD's
    adjointOutletVelocity (src/optimisation/adjointOptimisation/adjoint/
    adjointBoundaryConditions), without its explicit deviatoric flux and with
    no tangential objective source. See NOTICE.md.

    It is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version. It is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details: <http://www.gnu.org/licenses/>.

Class
    Foam::objectiveOutletVelocityFvPatchVectorField

Description
    Adjoint velocity at an outlet for a flow-rate objective (objective.H):

        normal:     Ua_n = phia/|Sf|
        tangential: Ua_t = Ua_t,cell nd/(max(U_n, 0) + nd), nd = nuEff/delta
                    (U_n Ua_t + nuEff d(Ua_t)/dn = 0, discretised), or zero

    nuEff and the choice of the tangential part are pushed by the solver
    every iteration (objectiveUpdate.H).

    Usage:
        outlet1 { type objectiveOutletVelocity; value uniform (0 0 0); }

SourceFiles
    objectiveOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveOutletVelocityFvPatchVectorField_H
#define objectiveOutletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objectiveOutletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- The tangential Robin condition (else Ua_t = 0)
        bool tangentialRobin_;

        //- Effective viscosity on the patch, set by the solver
        scalarField nuEff_;


public:

    //- Runtime type information
    TypeName("objectiveOutletVelocity");


    // Constructors

        objectiveOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const dictionary&
        );

        objectiveOutletVelocityFvPatchVectorField
        (
            const objectiveOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const fieldMapper&
        );

        objectiveOutletVelocityFvPatchVectorField
        (
            const objectiveOutletVelocityFvPatchVectorField&,
            const DimensionedField<vector, fvMesh>&
        );

        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new objectiveOutletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Setters (the solver pushes these every iteration)

            void tangentialRobin(const bool b)
            {
                tangentialRobin_ = b;
            }

            void nuEff(const scalarField& nu)
            {
                nuEff_ = nu;
            }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        using fixedValueFvPatchVectorField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
