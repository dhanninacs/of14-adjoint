// of14-adjoint: the solver's own controls, from the optional system/adjointDict.
//
//     porosity    false;   // the Foundation solver's topology update; off for shape sensitivity
//     ATC         { type standard; }   // standard | UaGradU | cancel (see ATC.H)
//
// With porosity on, lambda and alphaMax are read from physicalProperties as before.

IOdictionary adjointDict
(
    IOobject
    (
        "adjointDict",
        runTime.system(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

const bool porosity = adjointDict.lookupOrDefault<Switch>("porosity", false);

// The inlets: every patch that is not a wall or a constraint and on which U is fixed by a plain
// fixed-value condition. The Foundation solver looked for a patch called "inlet" by name.
// **Mixed conditions are not inlets**: an outlet's inletOutlet reports fixesValue() too, and
// counting it zeroed the ATC on the outlet cells as well (found by bisection against the
// Foundation solver; zeroing there is ESI's zeroATCPatchTypes (patch), an explicit choice).
labelList inletCells;
{
    DynamicList<label> cells;
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if
        (
            !isA<wallFvPatch>(patch)
         && !patch.constraint()
         && U.boundaryField()[patchi].fixesValue()
         && !isA<mixedFvPatchVectorField>(U.boundaryField()[patchi])
        )
        {
            cells.append(patch.faceCells());
        }
    }
    inletCells.transfer(cells);
}
Info<< "Inlet cells (ATC zeroed, porosity kept off): "
    << returnReduce(inletCells.size(), sumOp<label>()) << nl << endl;

// The objective's output: postProcessing/adjoint/0/objective.dat, one row per iteration (its
// header is written by objective.H, which knows the objective).
autoPtr<OFstream> objectiveFile;
if (Pstream::master())
{
    const fileName dir(runTime.path()/"postProcessing"/"adjoint"/runTime.name());
    mkDir(dir);
    objectiveFile.reset(new OFstream(dir/"objective.dat"));
}

// The wall shape sensitivity (shapeSensitivity.H): zero inside, written with the other fields.
volScalarField shapeSensitivity
(
    IOobject
    (
        "shapeSensitivity",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(sqr(dimLength)/pow3(dimTime), 0),
    calculatedFvPatchScalarField::typeName
);
