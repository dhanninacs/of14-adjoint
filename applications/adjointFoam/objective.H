// of14-adjoint: the objective the adjoint boundary conditions encode, from the optional
// `objective` sub-dictionary of system/adjointDict.
//
//     objective
//     {
//         type             flowRatePartition;   // ptLosses (default) | flowRatePartition
//         outletPatches    (out1 out2);         // flowRatePartition: >= 2 patch names
//         targetFractions  ();                  // optional, default 1/N each, patch-index order
//         viscousTerm      true;                // nuEff d(Ua_n)/dn in the outlet pa
//         tangential       robin;               // robin | zero  (outlet Ua_t)
//         startIteration   0;                   // the outlet source is held at 0 before this
//         convection       fvSchemes;           // or a scheme, e.g. bounded Gauss linearUpwind
//                                               // grad(Ua); fvSchemes: the case's div(-phi,Ua)
//         ATCInlets        legacy;              // legacy | strict (ATC kept in outlet cells)
//     }
//
// ptLosses (or no block): the Foundation solver's total-pressure loss, unchanged; the case uses
// adjointOutletPressure/adjointOutletVelocity with Ua = -U_in at the inlets.
//
// flowRatePartition: J = sum_i (Q_i/Q - T_i)^2 over the chosen outlets, Q_i = outlet flux, Q their
// sum (Stokes' flow_split; ESI's objectiveFlowRatePartition with J doubled and Q in place of the
// inlet flux). The case uses objectiveOutletPressure/objectiveOutletVelocity on those outlets and
// Ua = 0 at the inlets; objectiveUpdate.H sets the outlet sources every iteration.

const dictionary objectiveDict
(
    adjointDict.found("objective") ? adjointDict.subDict("objective") : dictionary()
);
const word objectiveType(objectiveDict.lookupOrDefault<word>("type", "ptLosses"));
if (objectiveType != "ptLosses" && objectiveType != "flowRatePartition")
{
    FatalIOErrorInFunction(objectiveDict)
        << "Unknown objective type " << objectiveType
        << "; use ptLosses or flowRatePartition" << exit(FatalIOError);
}
const bool flowRatePartition = (objectiveType == "flowRatePartition");

labelList objectiveOutlets;
scalarField objectiveTargets;
bool objectiveViscous = true;
bool objectiveRobin = true;
label objectiveStart = 0;

// The adjoint convection scheme of the objective path (adjointFoam.C), e.g.
// "bounded Gauss linearUpwind grad(Ua)"; the default, "fvSchemes", is the case's own
// div(-phi,Ua). linearUpwind lowered the snapped T-junction's error by 3 to 10 % but did not
// make R2 converge, so it is not the default (tests/flowsplit.md, Phase B). The pressure-loss
// path always uses fvSchemes.
string objectiveConvection("fvSchemes");

// The wall gradient of shapeSensitivity (adjointControls.H): uncorrected unless asked for, on
// both paths. Corrected fixed the body-fitted T-junction but made the snapped one worse
// (tests/flowsplit.md, Phase B). leastSquares is the most accurate on snapped and body-fitted
// walls against exact Poiseuille gradients, but with linearUpwind it still leaves the manifold
// at 0.80 / 0.84 and the 1 mm T-junction's R1 at 1.12 against finite differences, so it is
// not the default either (Phase C).
if (sensitivityWallGradient == "default")
{
    sensitivityWallGradient = "uncorrected";
}
const bool sensitivityCorrected = (sensitivityWallGradient == "corrected");
const bool sensitivityLeastSquares = (sensitivityWallGradient == "leastSquares");
if (flowRatePartition || sensitivityCorrected || sensitivityLeastSquares)
{
    Info<< "shapeSensitivity wall gradient: " << sensitivityWallGradient << nl << endl;
}

// The values of the last update (objectiveUpdate.H), logged by shapeSensitivity.H.
scalar objectiveJ = 0;
scalarField objectiveFractions;
scalarField objectiveDJdQ;

if (flowRatePartition)
{
    const wordList names(objectiveDict.lookup<wordList>("outletPatches"));
    {
        DynamicList<label> ids;
        forAll(names, i)
        {
            const label patchi = mesh.boundary().findIndex(names[i]);
            if (patchi < 0)
            {
                FatalIOErrorInFunction(objectiveDict)
                    << "outletPatches: no patch named " << names[i]
                    << exit(FatalIOError);
            }
            if (findIndex(ids, patchi) == -1)
            {
                ids.append(patchi);
            }
        }
        objectiveOutlets.transfer(ids);
    }
    // Patch-index order, as ESI (targetFractions follow it).
    sort(objectiveOutlets);

    const label nOut = objectiveOutlets.size();
    if (nOut < 2)
    {
        FatalIOErrorInFunction(objectiveDict)
            << "flowRatePartition needs at least two outletPatches, got " << names
            << exit(FatalIOError);
    }

    objectiveTargets = scalarField(nOut, 1.0/scalar(nOut));
    if (objectiveDict.found("targetFractions"))
    {
        const scalarList targets(objectiveDict.lookup<scalarList>("targetFractions"));
        if (targets.size())
        {
            if (targets.size() != nOut)
            {
                FatalIOErrorInFunction(objectiveDict)
                    << "targetFractions has " << targets.size() << " entries for "
                    << nOut << " outletPatches" << exit(FatalIOError);
            }
            if (mag(sum(targets) - 1) > 1e-6)
            {
                FatalIOErrorInFunction(objectiveDict)
                    << "targetFractions " << targets << " do not sum to 1"
                    << exit(FatalIOError);
            }
            objectiveTargets = targets;
        }
    }

    objectiveViscous = objectiveDict.lookupOrDefault<Switch>("viscousTerm", true);
    const word tangential(objectiveDict.lookupOrDefault<word>("tangential", "robin"));
    if (tangential != "robin" && tangential != "zero")
    {
        FatalIOErrorInFunction(objectiveDict)
            << "tangential " << tangential << "; use robin or zero" << exit(FatalIOError);
    }
    objectiveRobin = (tangential == "robin");
    objectiveStart = objectiveDict.lookupOrDefault<label>("startIteration", 0);
    // Written as in fvSchemes (convection bounded Gauss linearUpwind grad(Ua);), quoted or not.
    if (objectiveDict.found("convection"))
    {
        ITstream& is = objectiveDict.lookup("convection");
        if (is.size() == 1 && is[0].isWord())
        {
            objectiveConvection = is[0].wordToken();
        }
        else if (is.size() == 1 && is[0].isString())
        {
            objectiveConvection = is[0].stringToken();
        }
        else
        {
            OStringStream os;
            forAll(is, i)
            {
                os << (i ? " " : "") << is[i];
            }
            objectiveConvection = os.str();
        }
    }
    if (objectiveConvection != "fvSchemes")
    {
        // Fail at start-up on a scheme that does not parse.
        const surfaceScalarField phiNeg("-phi", -phi);
        IStringStream schemeData(objectiveConvection);
        fv::convectionScheme<vector>::New(mesh, phiNeg, schemeData);
    }

    // The outlets must carry the two conditions this objective drives.
    forAll(objectiveOutlets, i)
    {
        const label patchi = objectiveOutlets[i];
        if
        (
            !isA<objectiveOutletPressureFvPatchScalarField>(pa.boundaryField()[patchi])
         || !isA<objectiveOutletVelocityFvPatchVectorField>(Ua.boundaryField()[patchi])
        )
        {
            FatalIOErrorInFunction(objectiveDict)
                << "outlet " << mesh.boundary()[patchi].name()
                << " needs pa of type objectiveOutletPressure and Ua of type "
                << "objectiveOutletVelocity; it has "
                << pa.boundaryField()[patchi].type() << " and "
                << Ua.boundaryField()[patchi].type() << exit(FatalIOError);
        }
    }

    // The inlets must have Ua = 0: J does not depend on p. The rule of adjointControls.H, without
    // the objective's outlets and without directionMixed conditions (an outlet's
    // pressureInletOutletVelocity reports fixesValue() and is not mixed; adjointControls.H keeps
    // counting it, for the ATC limiter, so that the pressure-loss path is unchanged).
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if
        (
            !isA<wallFvPatch>(patch)
         && !patch.constraint()
         && findIndex(objectiveOutlets, patchi) == -1
         && U.boundaryField()[patchi].fixesValue()
         && !isA<mixedFvPatchVectorField>(U.boundaryField()[patchi])
         && !isA<directionMixedFvPatchVectorField>(U.boundaryField()[patchi])
        )
        {
            const scalar maxUa =
                returnReduce
                (
                    Ua.boundaryField()[patchi].size()
                  ? max(mag(Ua.boundaryField()[patchi]))
                  : scalar(0),
                    maxOp<scalar>()
                );
            if (maxUa > small)
            {
                FatalIOErrorInFunction(objectiveDict)
                    << "inlet " << patch.name() << " has max|Ua| = " << maxUa
                    << "; the flow-rate objective wants Ua = 0 at the inlets"
                    << exit(FatalIOError);
            }
        }
    }

    // The ATC limiter's inlet cells (ATC.H): "legacy" keeps adjointControls.H's rule, which also
    // counts an outlet's pressureInletOutletVelocity (directionMixed) as an inlet and so zeroes
    // the ATC in the outlet cells; "strict" zeroes it only next to the inlets found above.
    const word ATCInlets(objectiveDict.lookupOrDefault<word>("ATCInlets", "legacy"));
    if (ATCInlets != "legacy" && ATCInlets != "strict")
    {
        FatalIOErrorInFunction(objectiveDict)
            << "ATCInlets " << ATCInlets << "; use legacy or strict" << exit(FatalIOError);
    }
    if (ATCInlets == "strict")
    {
        ATClimiter = dimensionedScalar(dimless, 1);
        label nZeroed = 0;
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& patch = mesh.boundary()[patchi];
            const bool strictInlet =
                !isA<wallFvPatch>(patch)
             && !patch.constraint()
             && findIndex(objectiveOutlets, patchi) == -1
             && U.boundaryField()[patchi].fixesValue()
             && !isA<mixedFvPatchVectorField>(U.boundaryField()[patchi])
             && !isA<directionMixedFvPatchVectorField>(U.boundaryField()[patchi]);
            if
            (
                strictInlet
             || findIndex(zeroATCPatchTypes, word(patch.type())) != -1
            )
            {
                for (const label celli : patch.faceCells())
                {
                    ATClimiter[celli] = 0;
                    ++nZeroed;
                }
            }
        }
        ATClimiter.correctBoundaryConditions();
        const localMin<scalar> limiterScheme(mesh);
        for (label pass = 0; pass < ATCnSmooth; ++pass)
        {
            surfaceScalarField faceLimiter(limiterScheme.interpolate(ATClimiter));
            ATClimiter = fvc::average(faceLimiter);
            ATClimiter.correctBoundaryConditions();
        }
        Info<< "ATC limiter, strict inlets: zeroed in "
            << returnReduce(nZeroed, sumOp<label>()) << " face cells" << nl << endl;
    }

    objectiveFractions = scalarField(nOut, 0);
    objectiveDJdQ = scalarField(nOut, 0);

    Info<< "Objective: flowRatePartition over";
    forAll(objectiveOutlets, i)
    {
        Info<< " " << mesh.boundary()[objectiveOutlets[i]].name();
    }
    Info<< ", targets " << objectiveTargets
        << ", viscousTerm " << objectiveViscous
        << ", tangential " << tangential
        << ", startIteration " << objectiveStart
        << ", adjoint convection "
        << (
               objectiveConvection == "fvSchemes"
             ? string("fvSchemes div(-phi,Ua)")
             : objectiveConvection
           )
        << nl << endl;
}

// The objective's output: postProcessing/adjoint/0/objective.dat, one row per iteration.
if (objectiveFile.valid())
{
    if (flowRatePartition)
    {
        // Finite-difference checks difference J between runs: more than the default 6 digits.
        objectiveFile().precision(12);
        objectiveFile() << "# Time  J (flow-rate partition, dimensionless)";
        forAll(objectiveOutlets, i)
        {
            objectiveFile() << "  f_" << mesh.boundary()[objectiveOutlets[i]].name();
        }
        forAll(objectiveOutlets, i)
        {
            objectiveFile() << "  dJdQ_" << mesh.boundary()[objectiveOutlets[i]].name();
        }
        objectiveFile() << endl;
    }
    else
    {
        objectiveFile() << "# Time  J (total-pressure power loss, kinematic: m^5/s^3)" << endl;
    }
}
