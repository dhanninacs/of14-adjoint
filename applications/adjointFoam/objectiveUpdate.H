// of14-adjoint: the flow-rate partition objective's outlet sources, from the current primal flux.
// Included after the primal SIMPLE corrector and before the adjoint one (objective.H).
//
// J = sum_i e_i^2, e_i = f_i - T_i, f_i = Q_i/Q, Q = sum_k Q_k; with the solver's convention
// (Ua = -u, pa = -q of Othmer) the outlet adjoint pressure carries -dJ/dQ_i, and
//
//     dJ/dQ_i = (2/Q) (e_i - sum_k e_k f_k)
//
// J is homogeneous of degree 0 in the Q_i, so sum_i Q_i dJ/dQ_i = 0 (logged as "euler").

if (flowRatePartition)
{
    const label nOut = objectiveOutlets.size();
    scalarField Qout(nOut, 0);
    scalar Qsum = 0;
    scalar phiaSum = 0;
    scalar phiaAbs = 0;
    forAll(objectiveOutlets, i)
    {
        const label patchi = objectiveOutlets[i];
        Qout[i] = gSum(phi.boundaryField()[patchi]);
        Qsum += Qout[i];
        phiaSum += gSum(phia.boundaryField()[patchi]);
        phiaAbs += mag(gSum(phia.boundaryField()[patchi]));
    }

    if (mag(Qsum) > vSmall)
    {
        objectiveFractions = Qout/Qsum;
    }
    else
    {
        objectiveFractions = 0;
    }
    const scalarField e(objectiveFractions - objectiveTargets);
    objectiveJ = sum(sqr(e));
    const scalar ef = sum(e*objectiveFractions);
    if (mag(Qsum) > vSmall)
    {
        objectiveDJdQ = (2.0/Qsum)*(e - ef);
    }
    else
    {
        objectiveDJdQ = 0;
    }
    const scalar euler = sum(objectiveDJdQ*Qout);

    const bool active = (runTime.timeIndex() >= objectiveStart);
    const bool addUaU = (ATCType == "UaGradU");

    volScalarField::Boundary& pab = pa.boundaryFieldRef();
    volVectorField::Boundary& Uab = Ua.boundaryFieldRef();
    forAll(objectiveOutlets, i)
    {
        const label patchi = objectiveOutlets[i];
        const scalarField nuEffp(turbulence->nuEff(patchi));

        objectiveOutletPressureFvPatchScalarField& pap =
            refCast<objectiveOutletPressureFvPatchScalarField>(pab[patchi]);
        pap.source(active ? objectiveDJdQ[i] : 0);
        pap.addUaU(addUaU);
        pap.viscous(objectiveViscous);
        pap.nuEff(nuEffp);

        objectiveOutletVelocityFvPatchVectorField& Uap =
            refCast<objectiveOutletVelocityFvPatchVectorField>(Uab[patchi]);
        Uap.tangentialRobin(objectiveRobin);
        Uap.nuEff(nuEffp);
    }

    Info<< "flowRatePartition: J " << objectiveJ
        << " fractions " << objectiveFractions
        << " dJdQ " << objectiveDJdQ
        << " euler " << euler
        << " adjointFlux " << phiaSum << " (sum |port| " << phiaAbs << ")"
        << (active ? "" : " [source held at 0]") << endl;
}
