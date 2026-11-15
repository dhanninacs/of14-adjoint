// of14-adjoint: the least-squares gradient of a field in the cells next to one patch, for the
// wall gradient of shapeSensitivity (sensitivityWallGradient leastSquares).
//
// For each face of the patch, the gradient in its cell from a weighted least-squares fit over
// every face of the cell: the neighbour cell across an internal face, the neighbour across a
// coupled face, and the face value at the true face centre (Cf - C, not Foundation's
// normal-projected patch delta) on any other boundary face, weights 1/|d|^2. It is exact for a
// linear field whatever the cell's shape, which the Gauss gradient of a snapped (cut,
// non-orthogonal, skewed) wall cell is not. Measured against plane and pipe Poiseuille flow on
// snapped meshes and against finite differences in tests/flowsplit.md (Phase C).

template<class Type>
tmp<Field<typename outerProduct<vector, Type>::type>> wallCellGrad
(
    const VolField<Type>& vf,
    const label patchi
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vf.mesh();
    const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const volVectorField& C = mesh.C();
    const cellList& cells = mesh.cells();
    const polyBoundaryMesh& pbm = static_cast<const polyMesh&>(mesh).boundary();

    // Across coupled patches: the cell-to-cell vector and the neighbour value, once per patch
    // that the cells next to patchi touch.
    PtrList<vectorField> coupledDelta(mesh.boundary().size());
    PtrList<Field<Type>> coupledValue(mesh.boundary().size());

    // The directions a 2-D (empty) mesh does not solve: closed with a unit entry so dd is
    // invertible and the gradient has no component there.
    symmTensor emptyDirections(Zero);
    {
        const Vector<label>& sd = mesh.solutionD();
        if (sd.x() == -1) emptyDirections.xx() = 1;
        if (sd.y() == -1) emptyDirections.yy() = 1;
        if (sd.z() == -1) emptyDirections.zz() = 1;
    }

    tmp<Field<GradType>> tgrad(new Field<GradType>(faceCells.size(), Zero));
    Field<GradType>& grad = tgrad.ref();

    forAll(faceCells, i)
    {
        const label celli = faceCells[i];
        const cell& cFaces = cells[celli];

        symmTensor dd(emptyDirections);
        GradType dX(Zero);

        forAll(cFaces, j)
        {
            const label facei = cFaces[j];
            vector d;
            Type delta;
            if (mesh.isInternalFace(facei))
            {
                const label other = (own[facei] == celli) ? nei[facei] : own[facei];
                d = C[other] - C[celli];
                delta = vf[other] - vf[celli];
            }
            else
            {
                const label pi = pbm.whichPatch(facei);
                const fvPatch& p = mesh.boundary()[pi];
                if (isA<emptyFvPatch>(p))
                {
                    continue;
                }
                const label lf = facei - p.start();
                const fvPatchField<Type>& pf = vf.boundaryField()[pi];
                if (pf.coupled())
                {
                    if (!coupledDelta.set(pi))
                    {
                        coupledDelta.set(pi, p.delta().ptr());
                        coupledValue.set(pi, pf.patchNeighbourField().ptr());
                    }
                    d = coupledDelta[pi][lf];
                    delta = coupledValue[pi][lf] - vf[celli];
                }
                else
                {
                    d = p.Cf()[lf] - C[celli];
                    delta = pf[lf] - vf[celli];
                }
            }
            const scalar w = 1.0/max(magSqr(d), vSmall);
            dd += w*sqr(d);
            dX += w*(d*delta);
        }

        grad[i] = inv(dd) & dX;
    }

    return tgrad;
}
