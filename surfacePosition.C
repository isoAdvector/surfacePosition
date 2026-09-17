/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "surfacePosition.H"
#include "faceOnlySet.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "processorFvPatchFields.H"
#include "processorPolyPatch.H"
#include "PstreamBuffers.H"
#include "mapPolyMesh.H"
#include "Pair.H"
#include "FixedList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(surfacePosition, 0);
    addToRunTimeSelectionTable(functionObject, surfacePosition, dictionary);
}
}


const Foam::Enum
<
    Foam::functionObjects::surfacePosition::methodType
>
Foam::functionObjects::surfacePosition::methodTypeNames_
({
    { methodType::AUTO, "auto" },
    { methodType::RECONSTRUCTED, "reconstructed" },
    { methodType::PLIC, "plic" },
    { methodType::ALPHA, "alpha" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

using namespace Foam;

//- Volume of a tetrahedron
inline scalar tetVolume
(
    const point& a,
    const point& b,
    const point& c,
    const point& d
)
{
    return mag((b - a) & ((c - a) ^ (d - a)))/6;
}


//- Point where the plane crosses the edge p-q, from the signed distances
inline point edgeCut
(
    const point& p,
    const point& q,
    const scalar fp,
    const scalar fq
)
{
    return p + fp/(fp - fq)*(q - p);
}


//- Volume of the part of a tetrahedron with negative signed distance
scalar tetVolumeBelow
(
    const FixedList<point, 4>& p,
    const FixedList<scalar, 4>& f
)
{
    FixedList<label, 4> below;
    FixedList<label, 4> above;
    label nBelow = 0;
    label nAbove = 0;

    for (label i = 0; i < 4; ++i)
    {
        if (f[i] < 0)
        {
            below[nBelow++] = i;
        }
        else
        {
            above[nAbove++] = i;
        }
    }

    if (nBelow == 0)
    {
        return 0;
    }

    const scalar vTet = tetVolume(p[0], p[1], p[2], p[3]);

    if (nBelow == 4)
    {
        return vTet;
    }

    if (nBelow == 1 || nBelow == 3)
    {
        // Small tetrahedron at the single vertex on its own side
        const label i = (nBelow == 1 ? below[0] : above[0]);
        const label j0 = (nBelow == 1 ? above[0] : below[0]);
        const label j1 = (nBelow == 1 ? above[1] : below[1]);
        const label j2 = (nBelow == 1 ? above[2] : below[2]);

        const scalar vSmall = tetVolume
        (
            p[i],
            edgeCut(p[i], p[j0], f[i], f[j0]),
            edgeCut(p[i], p[j1], f[i], f[j1]),
            edgeCut(p[i], p[j2], f[i], f[j2])
        );

        return (nBelow == 1 ? vSmall : vTet - vSmall);
    }

    // Two vertices on each side: triangular prism between the edge a-b and
    // the cut, split into three tetrahedra
    const label a = below[0];
    const label b = below[1];
    const label c = above[0];
    const label e = above[1];

    const point pac(edgeCut(p[a], p[c], f[a], f[c]));
    const point pae(edgeCut(p[a], p[e], f[a], f[e]));
    const point pbc(edgeCut(p[b], p[c], f[b], f[c]));
    const point pbe(edgeCut(p[b], p[e], f[b], f[e]));

    return
        tetVolume(p[a], pac, pae, p[b])
      + tetVolume(pac, pae, p[b], pbc)
      + tetVolume(pae, p[b], pbc, pbe);
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::surfacePosition::crossesLine
(
    const label celli,
    const label li,
    scalar& sIn,
    scalar& sOut,
    label& entryFace
) const
{
    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& Cf = mesh_.faceCentres();
    const labelUList& own = mesh_.faceOwner();
    const point& start = starts_[li];
    const vector& dir = dirs_[li];

    // Clip the line parameter against the half-spaces of the cell faces
    scalar sMin = 0;
    scalar sMax = lengths_[li];
    entryFace = -1;

    for (const label facei : mesh_.cells()[celli])
    {
        const vector n = (own[facei] == celli ? Sf[facei] : -Sf[facei]);
        const scalar den = n & dir;
        const scalar num = n & (Cf[facei] - start);

        const scalar magN = mag(n);
        const scalar tol = 1e-9*magN*Foam::sqrt(magN);

        if (mag(den)*lengths_[li] <= tol)
        {
            // Line parallel to the face (its distance from the face plane
            // changes by less than 1e-9 of the face size along the whole
            // line): inside or outside for all s. A line on the face plane is
            // inside, so that a line on a face is not lost to round-off, which
            // would otherwise make num/den an arbitrary distance.
            if (num < -tol)
            {
                return false;
            }
            continue;
        }

        const scalar s = num/den;

        if (den > 0)
        {
            sMax = min(sMax, s);
        }
        else if (s > sMin)
        {
            sMin = s;
            entryFace = facei;
        }

        if (sMax <= sMin)
        {
            return false;
        }
    }

    sIn = sMin;
    sOut = sMax;

    return (sMax - sMin > SMALL*lengths_[li]);
}


void Foam::functionObjects::surfacePosition::searchSeeds
(
    const label li,
    labelList& seeds
)
{
    if (!searchPtr_)
    {
        searchPtr_.reset(new meshSearch(mesh_));
        searchStale_ = false;
    }
    else if (searchStale_)
    {
        searchPtr_->correct();
        searchStale_ = false;
    }

    const faceOnlySet set
    (
        "surfacePosition",
        mesh_,
        *searchPtr_,
        "xyz",
        starts_[li],
        starts_[li] + lengths_[li]*dirs_[li]
    );

    // The cells of the set seed a flood fill, which also finds cells crossed
    // by the line that the set misses. The set can miss all cells of a line
    // lying on faces of the mesh, so the cells nearest to points of the line
    // inside the bounds of the local mesh are added as seeds.
    DynamicList<label> found(set.cells());

    const point& start = starts_[li];
    const vector& dir = dirs_[li];
    scalar s0 = 0;
    scalar s1 = lengths_[li];

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        const scalar lo = localBounds_.min()[cmpt];
        const scalar hi = localBounds_.max()[cmpt];

        if (mag(dir[cmpt]) < VSMALL)
        {
            if (start[cmpt] < lo || start[cmpt] > hi)
            {
                s1 = -1;    // Outside the local bounds
            }
            continue;
        }

        scalar a = (lo - start[cmpt])/dir[cmpt];
        scalar b = (hi - start[cmpt])/dir[cmpt];
        if (a > b)
        {
            std::swap(a, b);
        }
        s0 = max(s0, a);
        s1 = min(s1, b);
    }

    if (s1 >= s0)
    {
        const label nPoints = 5;
        for (label i = 0; i < nPoints; ++i)
        {
            const point p(start + (s0 + (i + 0.5)*(s1 - s0)/nPoints)*dir);
            found.push_back(searchPtr_->findNearestCell(p));
        }
    }

    seeds.transfer(found);
}


void Foam::functionObjects::surfacePosition::floodLineCells
(
    const UList<labelList>& seeds,
    List<labelList>& cells
) const
{
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const cellList& meshCells = mesh_.cells();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const label nLines = seeds.size();

    // Candidate cells carry the number of rings of neighbours still to be
    // searched around them if they do not cross the line. The seeds start
    // with maxRings; the neighbours of cells crossing the line with none.
    const label maxRings = 10;

    typedef FixedList<label, 2> cellRings;

    // Per line: cells tested, cells found, the rings of neighbours still to
    // search, and the found cells whose neighbours are still to be tested
    List<labelHashSet> visited(nLines);
    List<DynamicList<label>> found(nLines);
    List<DynamicList<cellRings>> rings(nLines);
    List<DynamicList<label>> front(nLines);

    // Candidates for the neighbouring processors, per processor patch:
    // line, face index in the patch, rings
    List<DynamicList<FixedList<label, 3>>> send(pbm.size());

    // Send the cell across a processor face (other boundary faces are
    // ignored)
    auto sendAcross = [&](const label li, const label facei, const label n)
    {
        const label patchi = pbm.whichPatch(facei);
        if (isA<processorPolyPatch>(pbm[patchi]))
        {
            send[patchi].push_back
            (
                FixedList<label, 3>({li, facei - pbm[patchi].start(), n})
            );
        }
    };

    forAll(seeds, li)
    {
        for (const label celli : seeds[li])
        {
            if (celli >= 0 && celli < mesh_.nCells())
            {
                rings[li].push_back(cellRings({celli, maxRings}));
            }
        }
    }

    while (true)
    {
        for (label li = 0; li < nLines; ++li)
        {
            labelHashSet& vis = visited[li];
            DynamicList<label>& fnd = found[li];
            DynamicList<label>& frt = front[li];

            // Rings of neighbours, until a ring holds cells crossing the line
            DynamicList<cellRings>& ring = rings[li];

            while (!ring.empty())
            {
                const bool wasEmpty = fnd.empty();
                DynamicList<cellRings> next;
                DynamicList<cellRings> across;  // Boundary face, rings

                for (const cellRings& cr : ring)
                {
                    const label celli = cr[0];

                    if (!vis.insert(celli))
                    {
                        continue;
                    }

                    if (crossesLine(celli, li))
                    {
                        fnd.push_back(celli);
                        frt.push_back(celli);
                    }
                    else if (wasEmpty && cr[1] > 0)
                    {
                        for (const label facei : meshCells[celli])
                        {
                            if (mesh_.isInternalFace(facei))
                            {
                                next.push_back
                                (
                                    cellRings
                                    ({
                                        own[facei] == celli
                                      ? nei[facei]
                                      : own[facei],
                                        cr[1] - 1
                                    })
                                );
                            }
                            else
                            {
                                across.push_back(cellRings({facei, cr[1] - 1}));
                            }
                        }
                    }
                }

                if (!fnd.empty())
                {
                    ring.clear();
                }
                else
                {
                    ring.transfer(next);
                    for (const cellRings& fr : across)
                    {
                        sendAcross(li, fr[0], fr[1]);
                    }
                }
            }

            // Flood fill from the cells crossing the line
            while (!frt.empty())
            {
                const label celli = frt.back();
                frt.pop_back();

                for (const label facei : meshCells[celli])
                {
                    if (!mesh_.isInternalFace(facei))
                    {
                        sendAcross(li, facei, 0);
                        continue;
                    }

                    const label nbr =
                        (own[facei] == celli ? nei[facei] : own[facei]);

                    if (vis.insert(nbr) && crossesLine(nbr, li))
                    {
                        fnd.push_back(nbr);
                        frt.push_back(nbr);
                    }
                }
            }
        }

        if (!UPstream::parRun())
        {
            break;
        }

        label nSend = 0;
        for (const auto& list : send)
        {
            nSend += list.size();
        }

        if (returnReduce(nSend, sumOp<label>()) == 0)
        {
            break;
        }

        // Exchange the candidates, in the same order of the processor patches
        // on both sides
        PstreamBuffers pBufs;

        forAll(pbm, patchi)
        {
            const auto* ppp = isA<processorPolyPatch>(pbm[patchi]);
            if (ppp)
            {
                UOPstream toNbr(ppp->neighbProcNo(), pBufs);
                toNbr << send[patchi];
                send[patchi].clear();
            }
        }

        pBufs.finishedSends();

        forAll(pbm, patchi)
        {
            const auto* ppp = isA<processorPolyPatch>(pbm[patchi]);
            if (!ppp)
            {
                continue;
            }

            List<FixedList<label, 3>> received;
            {
                UIPstream fromNbr(ppp->neighbProcNo(), pBufs);
                fromNbr >> received;
            }

            const labelUList& faceCells = ppp->faceCells();

            for (const FixedList<label, 3>& r : received)
            {
                const label li = r[0];
                const label celli = faceCells[r[1]];

                if (r[2] > 0)
                {
                    if (!visited[li].found(celli))
                    {
                        rings[li].push_back(cellRings({celli, r[2]}));
                    }
                }
                else if
                (
                    visited[li].insert(celli)
                 && crossesLine(celli, li)
                )
                {
                    // Neighbour of a cell crossing the line on the other side
                    found[li].push_back(celli);
                    front[li].push_back(celli);
                }
            }
        }
    }

    cells.resize(nLines);
    forAll(cells, li)
    {
        cells[li].transfer(found[li]);
    }
}


void Foam::functionObjects::surfacePosition::setLineCells
(
    const label li,
    const labelUList& cells
)
{
    const vectorField& C = mesh_.cellCentres();

    // Entry and exit distances and entry faces
    scalarList sIn(cells.size());
    scalarList sOut(cells.size());
    labelList entry(cells.size());

    forAll(cells, i)
    {
        crossesLine(cells[i], li, sIn[i], sOut[i], entry[i]);
    }

    // Order the cells along the line
    labelList order;
    sortedOrder(sIn, order);

    lineCells_[li] = labelList(UIndirectList<label>(cells, order));
    lineIn_[li] = scalarList(UIndirectList<scalar>(sIn, order));
    lineOut_[li] = scalarList(UIndirectList<scalar>(sOut, order));
    lineEntry_[li] = labelList(UIndirectList<label>(entry, order));

    scalarList& dist = lineDist_[li];
    dist.resize(cells.size());
    forAll(dist, i)
    {
        dist[i] = (C[lineCells_[li][i]] - starts_[li]) & dirs_[li];
    }
}


void Foam::functionObjects::surfacePosition::calcLineCells()
{
    const label nLines = starts_.size();

    bool seeded = (seeds_.size() == nLines);

    if (seeded)
    {
        ++nRebuilds_;

        if (fullSearchInterval_ > 0 && nRebuilds_ % fullSearchInterval_ == 0)
        {
            seeded = false;
        }
    }

    if (!seeded)
    {
        // Bounds of the mesh on this processor, for the full search
        localBounds_ = boundBox(mesh_.points(), false);
        localBounds_.inflate(1e-6);
    }

    lineCells_.resize(nLines);
    lineDist_.resize(nLines);
    lineIn_.resize(nLines);
    lineOut_.resize(nLines);
    lineEntry_.resize(nLines);

    label nCells = 0;
    label nMissing = 0;
    label nExtra = 0;

    // Seeds: the old line cells, or cells found with the search engine
    List<labelList> seeds;

    if (seeded)
    {
        seeds.transfer(seeds_);
    }
    else
    {
        seeds.resize(nLines);
        for (label li = 0; li < nLines; ++li)
        {
            searchSeeds(li, seeds[li]);
        }
    }

    List<labelList> lineCells;
    floodLineCells(seeds, lineCells);

    for (label li = 0; li < nLines; ++li)
    {
        const labelList& cells = lineCells[li];

        if (seeded && checkLineCells_)
        {
            // Reference: all cells crossed by the line
            DynamicList<label> reference;
            for (label celli = 0; celli < mesh_.nCells(); ++celli)
            {
                if (crossesLine(celli, li))
                {
                    reference.push_back(celli);
                }
            }

            const labelHashSet found(cells);
            const labelHashSet expected(reference);

            for (const label celli : reference)
            {
                if (!found.found(celli)) ++nMissing;
            }
            for (const label celli : cells)
            {
                if (!expected.found(celli)) ++nExtra;
            }
        }

        nCells += cells.size();
        setLineCells(li, cells);
    }

    seeds_.clear();
    update_ = false;

    if (seeded && checkLineCells_)
    {
        reduce(nCells, sumOp<label>());
        reduce(nMissing, sumOp<label>());
        reduce(nExtra, sumOp<label>());

        Info<< type() << ' ' << name() << ": rebuild " << nRebuilds_
            << " at time " << mesh_.time().timeName()
            << ": line cells " << nCells
            << ", missing " << nMissing
            << ", extra " << nExtra << endl;
    }
    else if (!seeded)
    {
        reduce(nCells, sumOp<label>());

        Info<< type() << ' ' << name() << ": searched " << nLines
            << " lines at time " << mesh_.time().timeName()
            << ", line cells " << nCells << endl;
    }
}


bool Foam::functionObjects::surfacePosition::predecessor
(
    const label li,
    const label k,
    const volScalarField& alpha,
    scalar& value,
    scalar* distPtr
) const
{
    const label facei = lineEntry_[li][k];

    if (facei < 0)
    {
        return false;   // The line starts in the cell
    }

    const label celli = lineCells_[li][k];

    if (mesh_.isInternalFace(facei))
    {
        const label own = mesh_.faceOwner()[facei];
        const label nbr = (own == celli ? mesh_.faceNeighbour()[facei] : own);
        value = alpha[nbr];

        if (distPtr)
        {
            *distPtr = (mesh_.C()[nbr] - starts_[li]) & dirs_[li];
        }
        return true;
    }

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const label patchi = pbm.whichPatch(facei);
    const fvPatchScalarField& pf = alpha.boundaryField()[patchi];

    if (!pf.coupled() || pf.empty())
    {
        return false;   // The line enters the domain here
    }

    // Processor patch values are the neighbour cell values
    const label i = facei - pbm[patchi].start();
    value =
    (
        isA<processorFvPatchScalarField>(pf)
      ? pf[i]
      : pf.patchNeighbourField()()[i]
    );

    if (distPtr)
    {
        // Neighbour cell centre from the cell centre and the patch delta
        const point nbrCentre
        (
            mesh_.C()[celli] + mesh_.boundary()[patchi].delta()()[i]
        );
        *distPtr = (nbrCentre - starts_[li]) & dirs_[li];
    }

    return true;
}


void Foam::functionObjects::surfacePosition::firstInterface
(
    const label li,
    const volScalarField& alpha,
    const volVectorField* normalPtr,
    const volVectorField* centrePtr,
    const pointScalarField* pointAlphaPtr,
    scalar& event,
    scalar& fallback
) const
{
    const point& start = starts_[li];
    const vector& dir = dirs_[li];
    const labelList& cells = lineCells_[li];
    const scalarList& sIn = lineIn_[li];

    event = GREAT;
    fallback = GREAT;

    forAll(cells, k)
    {
        if (sIn[k] > event)
        {
            break;  // The cells further on cannot hold an earlier event
        }

        const label celli = cells[k];
        const scalar a = alpha[celli];

        // Sharp jump of alpha at the entry face, between two full cells
        if (a <= jumpTol_ || a >= 1 - jumpTol_)
        {
            scalar aPrev = 0;

            if
            (
                predecessor(li, k, alpha, aPrev)
             && (aPrev <= jumpTol_ || aPrev >= 1 - jumpTol_)
             && (aPrev - 0.5)*(a - 0.5) < 0
            )
            {
                event = min(event, sIn[k]);
                continue;
            }
        }

        // Interface plane in the cell
        vector n(Zero);
        scalar d = 0;

        if (normalPtr)
        {
            n = (*normalPtr)[celli];
            const scalar magN = mag(n);

            if (magN < VSMALL)
            {
                continue;   // Not an interface cell
            }

            n /= magN;
            d = n & (*centrePtr)[celli];
        }
        else
        {
            if (a <= alphaTol_ || a >= 1 - alphaTol_)
            {
                continue;   // Not an interface cell
            }

            // Normal pointing out of the phase alpha
            n = -alphaGradient(celli, *pointAlphaPtr);
            const scalar magN = mag(n);

            if (magN < VSMALL)
            {
                continue;
            }

            n /= magN;

            // Plane (n & x) = d with the phase alpha on the side (n & x) < d
            d = plicPlane(celli, a, n);
        }

        const scalar dirN = dir & n;

        if (mag(dirN) <= SMALL)
        {
            continue;   // Line parallel to the interface plane
        }

        const scalar s = (d - (n & start))/dirN;

        if (insideCell(li, k, s))
        {
            event = min(event, s);
        }
        else if (fallback == GREAT)
        {
            fallback = s;
        }
    }
}


Foam::vector Foam::functionObjects::surfacePosition::alphaGradient
(
    const label celli,
    const pointScalarField& pointAlpha
) const
{
    const vectorField& Sf = mesh_.faceAreas();
    const labelUList& own = mesh_.faceOwner();
    const faceList& faces = mesh_.faces();

    // Gauss gradient with face values from the point values (Youngs)
    vector grad(Zero);

    for (const label facei : mesh_.cells()[celli])
    {
        const face& f = faces[facei];

        scalar af = 0;
        for (const label pointi : f)
        {
            af += pointAlpha[pointi];
        }
        af /= f.size();

        grad += (own[facei] == celli ? af : -af)*Sf[facei];
    }

    return grad;
}


Foam::scalar Foam::functionObjects::surfacePosition::volumeBelow
(
    const label celli,
    const vector& n,
    const scalar d
) const
{
    const pointField& points = mesh_.points();
    const faceList& faces = mesh_.faces();
    const vectorField& Cf = mesh_.faceCentres();
    const point& cc = mesh_.cellCentres()[celli];
    const scalar fc = (n & cc) - d;

    scalar volume = 0;

    // Tetrahedra from the cell centre, a face centre and a face edge
    for (const label facei : mesh_.cells()[celli])
    {
        const face& f = faces[facei];
        const scalar ff = (n & Cf[facei]) - d;

        forAll(f, fp)
        {
            const point& p1 = points[f[fp]];
            const point& p2 = points[f.nextLabel(fp)];

            volume += tetVolumeBelow
            (
                FixedList<point, 4>({cc, Cf[facei], p1, p2}),
                FixedList<scalar, 4>
                ({
                    fc,
                    ff,
                    (n & p1) - d,
                    (n & p2) - d
                })
            );
        }
    }

    return volume;
}


Foam::scalar Foam::functionObjects::surfacePosition::plicPlane
(
    const label celli,
    const scalar alpha1,
    const vector& n
) const
{
    const pointField& points = mesh_.points();
    const faceList& faces = mesh_.faces();

    scalar dMin = GREAT;
    scalar dMax = -GREAT;

    for (const label facei : mesh_.cells()[celli])
    {
        for (const label pointi : faces[facei])
        {
            const scalar d = n & points[pointi];
            dMin = min(dMin, d);
            dMax = max(dMax, d);
        }
    }

    const scalar volume = volumeBelow(celli, n, dMax + (dMax - dMin));
    const scalar target = alpha1*volume;
    const scalar tol = 1e-12*(dMax - dMin);

    // Bisection on the monotonic volume below the plane
    for (label iter = 0; iter < 100 && dMax - dMin > tol; ++iter)
    {
        const scalar d = 0.5*(dMin + dMax);

        if (volumeBelow(celli, n, d) < target)
        {
            dMin = d;
        }
        else
        {
            dMax = d;
        }
    }

    return 0.5*(dMin + dMax);
}


Foam::scalar Foam::functionObjects::surfacePosition::alphaDistance
(
    const label li,
    const volScalarField& alpha
) const
{
    const labelList& cells = lineCells_[li];
    const scalarList& dist = lineDist_[li];
    const scalarList& sIn = lineIn_[li];

    scalar best = GREAT;

    // Crossings of alpha = 0.5 between each cell and the cell before it
    // along the line, each found only by the later cell
    forAll(cells, k)
    {
        if (sIn[k] > best)
        {
            break;
        }

        scalar aPrev = 0;
        scalar sPrev = 0;

        if (!predecessor(li, k, alpha, aPrev, &sPrev))
        {
            continue;
        }

        const scalar a0 = aPrev - 0.5;
        const scalar a1 = alpha[cells[k]] - 0.5;

        if (a0*a1 < 0 || (a1 == 0 && a0 != 0))
        {
            best = min(best, sPrev + a0/(a0 - a1)*(dist[k] - sPrev));
        }
    }

    return best;
}


void Foam::functionObjects::surfacePosition::writeFileHeader(Ostream& os)
{
    writeHeader(os, "Interface position along lines");

    writeHeaderValue(os, "Method", methodTypeNames_[usedMethod_]);

    if (useLocations_)
    {
        writeHeaderValue(os, "Upward direction", up_);
    }

    forAll(starts_, li)
    {
        const string line("Line " + Foam::name(li));
        writeHeaderValue(os, line + " start", starts_[li]);
        writeHeaderValue(os, line + " direction", dirs_[li]);
        writeHeaderValue(os, line + " length", lengths_[li]);

        if (useLocations_)
        {
            writeHeaderValue(os, line + " location", locations_[li]);
        }
    }

    writeCommented(os, "Time");

    forAll(starts_, li)
    {
        const word n(Foam::name(li));

        writeTabbed(os, "s" + n);

        if (useLocations_)
        {
            writeTabbed(os, "h" + n);
        }

        writeTabbed(os, "x" + n);
        writeTabbed(os, "y" + n);
        writeTabbed(os, "z" + n);
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::surfacePosition::surfacePosition
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    alphaName_("alpha.water"),
    method_(methodType::AUTO),
    usedMethod_(methodType::AUTO),
    alphaTol_(1e-8),
    checkLineCells_(false),
    fullSearchInterval_(0),
    nRebuilds_(0),
    useLocations_(true),
    locations_(),
    up_(Zero),
    starts_(),
    dirs_(),
    lengths_(),
    lineCells_(),
    lineDist_(),
    seeds_(),
    searchPtr_(nullptr),
    searchStale_(false),
    update_(true),
    headerDone_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::surfacePosition::read(const dictionary& dict)
{
    if (!(fvMeshFunctionObject::read(dict) && writeFile::read(dict)))
    {
        return false;
    }

    alphaName_ = dict.getOrDefault<word>("alpha", "alpha.water");
    method_ = methodTypeNames_.getOrDefault("method", dict, methodType::AUTO);
    alphaTol_ = dict.getOrDefault<scalar>("alphaTol", 1e-8);
    checkLineCells_ = dict.getOrDefault("checkLineCells", false);
    fullSearchInterval_ = dict.getOrDefault<label>("fullSearchInterval", 0);

    if (dict.found("lines"))
    {
        useLocations_ = false;

        const List<Pair<point>> lines(dict.get<List<Pair<point>>>("lines"));

        starts_.resize(lines.size());
        dirs_.resize(lines.size());
        lengths_.resize(lines.size());

        forAll(lines, li)
        {
            const vector d = lines[li].second() - lines[li].first();
            lengths_[li] = mag(d);

            if (lengths_[li] < VSMALL)
            {
                FatalIOErrorInFunction(dict)
                    << "Line " << li << " has zero length"
                    << exit(FatalIOError);
            }

            starts_[li] = lines[li].first();
            dirs_[li] = d/lengths_[li];
        }
    }
    else
    {
        useLocations_ = true;
        locations_ = dict.get<List<point>>("locations");

        vector direction(Zero);

        if (!dict.readIfPresent("direction", direction))
        {
            const auto* gPtr =
                mesh_.time().cfindObject<uniformDimensionedVectorField>("g");

            if (!gPtr)
            {
                FatalIOErrorInFunction(dict)
                    << "No direction given and no gravity (g) registered"
                    << exit(FatalIOError);
            }

            direction = -gPtr->value();
        }

        if (mag(direction) < VSMALL)
        {
            FatalIOErrorInFunction(dict)
                << "The direction is zero"
                << exit(FatalIOError);
        }

        up_ = direction/mag(direction);

        // Lines from above to below through the locations, across the whole
        // mesh wherever the location is, with a margin for mesh motion
        const boundBox& bb = mesh_.bounds();
        const tmp<pointField> tcorners(bb.points());
        scalar hMin = GREAT;
        scalar hMax = -GREAT;

        for (const point& p : tcorners())
        {
            hMin = min(hMin, up_ & p);
            hMax = max(hMax, up_ & p);
        }

        const scalar margin = bb.mag();

        starts_.resize(locations_.size());
        dirs_ = List<vector>(locations_.size(), -up_);
        lengths_ = List<scalar>(locations_.size(), hMax - hMin + 2*margin);

        forAll(locations_, li)
        {
            const scalar hLoc = up_ & locations_[li];
            starts_[li] = locations_[li] + (hMax + margin - hLoc)*up_;
        }
    }

    seeds_.clear();
    update_ = true;

    return true;
}


bool Foam::functionObjects::surfacePosition::execute()
{
    return true;
}


bool Foam::functionObjects::surfacePosition::write()
{
    if (update_)
    {
        calcLineCells();
    }

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    const word normalName
    (
        IOobject::groupName("interfaceNormal", alpha.group())
    );
    const word centreName
    (
        IOobject::groupName("interfaceCentre", alpha.group())
    );

    const volVectorField* normalPtr = nullptr;
    const volVectorField* centrePtr = nullptr;

    usedMethod_ = method_;

    if (method_ == methodType::AUTO || method_ == methodType::RECONSTRUCTED)
    {
        normalPtr = mesh_.cfindObject<volVectorField>(normalName);
        centrePtr = mesh_.cfindObject<volVectorField>(centreName);

        if (normalPtr && centrePtr)
        {
            usedMethod_ = methodType::RECONSTRUCTED;
        }
        else if (method_ == methodType::RECONSTRUCTED)
        {
            FatalErrorInFunction
                << "Method reconstructed needs the fields " << normalName
                << " and " << centreName << " of an interface"
                << " reconstruction, which are not registered"
                << exit(FatalError);
        }
        else
        {
            usedMethod_ = methodType::PLIC;
        }
    }

    // Point values of alpha for the PLIC normals. The interpolation covers the
    // whole mesh and is collective, so it is done on all processors.
    tmp<pointScalarField> tpointAlpha;

    if (usedMethod_ == methodType::PLIC)
    {
        tpointAlpha = volPointInterpolation::New(mesh_).interpolate(alpha);
    }

    // First interface along each line: genuine events (interface inside its
    // cell, or a sharp jump), and as a fallback an interface plane cutting
    // the line outside its cell. One reduction over all lines and both.
    const label nLines = starts_.size();
    scalarList found(2*nLines, GREAT);

    for (label li = 0; li < nLines; ++li)
    {
        if (usedMethod_ == methodType::ALPHA)
        {
            found[li] = alphaDistance(li, alpha);
        }
        else
        {
            firstInterface
            (
                li,
                alpha,
                normalPtr,
                centrePtr,
                tpointAlpha.get(),
                found[li],
                found[nLines + li]
            );
        }
    }

    Pstream::listReduce(found, minOp<scalar>());

    scalarList dist(nLines);
    forAll(dist, li)
    {
        dist[li] = (found[li] < 0.5*GREAT ? found[li] : found[nLines + li]);
    }

    if (Pstream::master())
    {
        if (!headerDone_)
        {
            writeFileHeader(file());
            headerDone_ = true;
        }

        writeCurrentTime(file());

        forAll(dist, li)
        {
            if (dist[li] < 0.5*GREAT)
            {
                const point p = starts_[li] + dist[li]*dirs_[li];

                file() << tab << dist[li];

                if (useLocations_)
                {
                    file() << tab << ((p - locations_[li]) & up_);
                }

                file() << tab << p.x() << tab << p.y() << tab << p.z();
            }
            else
            {
                const label nColumns = useLocations_ ? 5 : 4;

                for (label i = 0; i < nColumns; ++i)
                {
                    file() << tab << "nan";
                }
            }
        }

        file() << endl;
    }

    return true;
}


void Foam::functionObjects::surfacePosition::updateMesh
(
    const mapPolyMesh& mpm
)
{
    if (&mpm.mesh() != &mesh_)
    {
        return;
    }


    // Cells before this change: pending seeds, or the current line cells
    const List<labelList> oldCells
    (
        (update_ && seeds_.size() == lineCells_.size()) ? seeds_ : lineCells_
    );

    const labelList& reverseCellMap = mpm.reverseCellMap();

    seeds_.resize(oldCells.size());

    forAll(oldCells, li)
    {
        DynamicList<label> newCells(oldCells[li].size());

        for (const label oldCelli : oldCells[li])
        {
            if (oldCelli >= reverseCellMap.size())
            {
                continue;
            }

            label newCelli = reverseCellMap[oldCelli];

            if (newCelli < -1)
            {
                newCelli = -newCelli - 2;   // Merged into this cell
            }

            if (newCelli >= 0)
            {
                newCells.push_back(newCelli);
            }
        }

        seeds_[li].transfer(newCells);
    }

    update_ = true;
    searchStale_ = true;
}


void Foam::functionObjects::surfacePosition::movePoints(const polyMesh& mesh)
{
    if (&mesh != &mesh_)
    {
        return;
    }


    // Cell labels are unchanged: the current line cells are the seeds
    if (!(update_ && seeds_.size() == lineCells_.size()))
    {
        seeds_ = lineCells_;
    }

    update_ = true;
    searchStale_ = true;
}


// ************************************************************************* //
