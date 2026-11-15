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
        << ", startIteration " << objectiveStart << nl << endl;
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
