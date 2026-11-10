/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Johan Roenby
-------------------------------------------------------------------------------
License
    This file is an extension of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "surfacePosition.H"
#include "faceOnlySet.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "surfaceIteratorIso.H"
#include "processorPolyPatch.H"
#include "PstreamBuffers.H"
#include "mapPolyMesh.H"
#include "FixedList.H"
#include "addToRunTimeSelectionTable.H"

#include <type_traits>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// Compatibility with OpenFOAM.com v2006 and later: DynamicList::push_back came
// with v2212 (append before), and Pstream::listGatherValues gathers
// non-contiguous types only from v2406

namespace
{

using namespace Foam;

//- Append an element to a list
template<class T, int SizeMin>
inline void pushBack
(
    DynamicList<T, SizeMin>& list,
    const typename std::common_type<T>::type& val
)
{
#if (OPENFOAM >= 2212)
    list.push_back(val);
#else
    list.append(val);
#endif
}


//- The last element of a list
template<class T, int SizeMin>
inline T& lastItem(DynamicList<T, SizeMin>& list)
{
    return list[list.size() - 1];
}


//- Remove the last element of a list and return it
template<class T, int SizeMin>
inline T popBack(DynamicList<T, SizeMin>& list)
{
    T val(list[list.size() - 1]);
    list.resize(list.size() - 1);
    return val;
}


//- The values of all processors, on the master (empty on the others)
template<class T>
List<T> gatherToMaster(const T& local)
{
    List<T> all;

    if (!UPstream::parRun())
    {
        all.resize(1);
        all[0] = local;
        return all;
    }

    PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);

    if (!UPstream::master())
    {
        UOPstream toMaster(UPstream::masterNo(), pBufs);
        toMaster << local;
    }

    pBufs.finishedSends();

    if (UPstream::master())
    {
        all.resize(UPstream::nProcs());

        forAll(all, proci)
        {
            if (proci == UPstream::masterNo())
            {
                all[proci] = local;
            }
            else
            {
                UIPstream fromProc(proci, pBufs);
                fromProc >> all[proci];
            }
        }
    }

    return all;
}

} // End anonymous namespace


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
    Foam::functionObjects::surfacePosition::modeType
>
Foam::functionObjects::surfacePosition::modeTypeNames_
({
    { modeType::ACCUMULATED, "accumulated" },
    { modeType::HIGHEST, "highest" },
    { modeType::LOWEST, "lowest" },
});



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::surfacePosition::setLines()
{
    const auto* gPtr =
        mesh_.time().cfindObject<uniformDimensionedVectorField>("g");

    if (mag(direction_) > VSMALL)
    {
        down_ = direction_/mag(direction_);

        // The direction points downwards, as in interfaceHeight, whose
        // documentation calls it the direction of the interface: a direction
        // against gravity is most likely an upward one given by mistake
        if (gPtr && (down_ & gPtr->value()) < 0)
        {
            WarningInFunction
                << type() << ' ' << name() << ": direction " << direction_
                << " points against gravity " << gPtr->value() << nl
                << "    The direction points downwards, along gravity, as in"
                << " interfaceHeight; the heights are measured against it."
                << endl;
        }
    }
    else
    {
        if (!gPtr || mag(gPtr->value()) < VSMALL)
        {
            FatalErrorInFunction
                << "No direction given and no gravity (g) registered"
                << exit(FatalError);
        }

        down_ = gPtr->value()/mag(gPtr->value());
    }

    // Lines from above to below through the locations, across the whole mesh
    // wherever the location is, with a margin for mesh motion
    const boundBox& bb = mesh_.bounds();
    const tmp<pointField> tcorners(bb.points());
    scalar hMin = GREAT;
    scalar hMax = -GREAT;

    for (const point& p : tcorners())
    {
        hMin = min(hMin, -down_ & p);
        hMax = max(hMax, -down_ & p);
    }

    const scalar margin = bb.mag();

    length_ = hMax - hMin + 2*margin;
    starts_.resize(locations_.size());

    forAll(locations_, li)
    {
        const scalar hLoc = -down_ & locations_[li];
        starts_[li] = locations_[li] - (hMax + margin - hLoc)*down_;
    }
}


bool Foam::functionObjects::surfacePosition::crossesLine
(
    const label celli,
    const label li,
    scalar& sIn,
    scalar& sOut
) const
{
    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& Cf = mesh_.faceCentres();
    const labelUList& own = mesh_.faceOwner();
    const point& start = starts_[li];
    const vector& dir = down_;

    // Clip the line parameter against the half-spaces of the cell faces
    scalar sMin = 0;
    scalar sMax = length_;

    for (const label facei : mesh_.cells()[celli])
    {
        const vector n = (own[facei] == celli ? Sf[facei] : -Sf[facei]);
        const scalar den = n & dir;
        const scalar num = n & (Cf[facei] - start);

        const scalar magN = mag(n);
        const scalar tol = 1e-9*magN*Foam::sqrt(magN);

        if (mag(den)*length_ <= tol)
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
        else
        {
            sMin = max(sMin, s);
        }

        if (sMax <= sMin)
        {
            return false;
        }
    }

    sIn = sMin;
    sOut = sMax;

    return (sMax - sMin > SMALL*length_);
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
        starts_[li] + length_*down_
    );

    // The cells of the set seed a flood fill, which also finds cells crossed
    // by the line that the set misses. The set can miss all cells of a line
    // lying on faces of the mesh, so the cells nearest to points of the line
    // inside the bounds of the local mesh are added as seeds.
    DynamicList<label> found(set.cells());

    const point& start = starts_[li];
    const vector& dir = down_;
    scalar s0 = 0;
    scalar s1 = length_;

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
            pushBack(found, searchPtr_->findNearestCell(p));
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
            pushBack
            (
                send[patchi],
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
                pushBack(rings[li], cellRings({celli, maxRings}));
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
                        pushBack(fnd, celli);
                        pushBack(frt, celli);
                    }
                    else if (wasEmpty && cr[1] > 0)
                    {
                        for (const label facei : meshCells[celli])
                        {
                            if (mesh_.isInternalFace(facei))
                            {
                                pushBack
                                (
                                    next,
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
                                pushBack(across, cellRings({facei, cr[1] - 1}));
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
                const label celli = popBack(frt);

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
                        pushBack(fnd, nbr);
                        pushBack(frt, nbr);
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
        PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);

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
                        pushBack(rings[li], cellRings({celli, r[2]}));
                    }
                }
                else if
                (
                    visited[li].insert(celli)
                 && crossesLine(celli, li)
                )
                {
                    // Neighbour of a cell crossing the line on the other side
                    pushBack(found[li], celli);
                    pushBack(front[li], celli);
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
    scalarList sIn(cells.size());
    scalarList sOut(cells.size());

    forAll(cells, i)
    {
        crossesLine(cells[i], li, sIn[i], sOut[i]);
    }

    // Order the cells along the line
    labelList order;
    sortedOrder(sIn, order);

    lineCells_[li] = labelList(UIndirectList<label>(cells, order));
    lineIn_[li] = scalarList(UIndirectList<scalar>(sIn, order));
    lineOut_[li] = scalarList(UIndirectList<scalar>(sOut, order));
}


void Foam::functionObjects::surfacePosition::calcCover()
{
    const label nLines = starts_.size();

    // The segments of the line cells of all processors
    List<List<Pair<scalar>>> local(nLines);

    forAll(local, li)
    {
        List<Pair<scalar>>& segments = local[li];
        segments.resize(lineIn_[li].size());

        forAll(segments, k)
        {
            segments[k] = Pair<scalar>(lineIn_[li][k], lineOut_[li][k]);
        }
    }

    const List<List<List<Pair<scalar>>>> all
    (
        gatherToMaster(local)
    );

    if (!Pstream::master())
    {
        return;
    }

    breaks_.resize(nLines);
    cover_.resize(nLines);

    for (label li = 0; li < nLines; ++li)
    {
        // The ends of all segments, sorted and unique, divide the line into
        // parts, each covered by a fixed number of cells
        DynamicList<scalar> ends;

        for (const auto& proc : all)
        {
            for (const Pair<scalar>& seg : proc[li])
            {
                pushBack(ends, seg.first());
                pushBack(ends, seg.second());
            }
        }

        Foam::sort(ends);

        scalarList& breaks = breaks_[li];
        breaks.resize(ends.size());
        label n = 0;

        for (const scalar s : ends)
        {
            if (n == 0 || s > breaks[n - 1])
            {
                breaks[n++] = s;
            }
        }
        breaks.resize(n);

        labelList& cover = cover_[li];
        cover.resize(max(n - 1, 0));
        cover = 0;

        for (const auto& proc : all)
        {
            for (const Pair<scalar>& seg : proc[li])
            {
                const label first = findLower(breaks, seg.first()) + 1;

                for
                (
                    label j = first;
                    j < cover.size() && breaks[j] < seg.second();
                    ++j
                )
                {
                    ++cover[j];
                }
            }
        }
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
    lineIn_.resize(nLines);
    lineOut_.resize(nLines);

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
                    pushBack(reference, celli);
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

    calcCover();

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


void Foam::functionObjects::surfacePosition::liquidParts
(
    const label li,
    const volScalarField& alpha,
    List<Pair<scalar>>& parts
)
{
    const point& start = starts_[li];
    const labelList& cells = lineCells_[li];
    const scalarList& sIn = lineIn_[li];
    const scalarList& sOut = lineOut_[li];

    DynamicList<Pair<scalar>> found;

    // Add the part [s0, s1] of the line, merged with the last part if they
    // touch (which leaves the covered length unchanged)
    auto add = [&found](const scalar s0, const scalar s1)
    {
        if (s1 <= s0)
        {
            return;
        }
        if (found.size() && lastItem(found).second() == s0)
        {
            lastItem(found).second() = s1;
        }
        else
        {
            pushBack(found, Pair<scalar>(s0, s1));
        }
    };

    forAll(cells, k)
    {
        const label celli = cells[k];
        const scalar a = alpha[celli];

        if (a >= 1 - surfCellTol_)
        {
            add(sIn[k], sOut[k]);
            continue;
        }
        if (a <= surfCellTol_)
        {
            continue;
        }

        // Interface cell: the part on the liquid side of the isosurface. The
        // status returned for cells within isoFaceTol of full or empty is not
        // that of the cut, so those are decided by their volume fraction.
        const label status =
            isoPtr_->vofCutCell(celli, a, isoFaceTol_, maxIter_);

        const vector& n = isoPtr_->surfaceArea();
        const scalar magN = mag(n);

        if (status != 0 || magN < VSMALL)
        {
            if (a >= 0.5)
            {
                add(sIn[k], sOut[k]);
            }
            continue;
        }

        // Normal of the isosurface into the liquid, i.e. towards the centre of
        // the submerged part of the cell (the orientation of the isosurface
        // area is not relied on). The liquid is where (x - centre) & m > 0,
        // with x = start + s*down.
        const point& centre = isoPtr_->surfaceCentre();
        const vector m
        (
            ((isoPtr_->subCellCentre() - centre) & n) > 0 ? n : -n
        );
        const scalar dm = down_ & m;
        const scalar f0 = (start - centre) & m;

        if (mag(dm) <= SMALL*magN)
        {
            // Line parallel to the isosurface
            if (f0 + 0.5*(sIn[k] + sOut[k])*dm > 0)
            {
                add(sIn[k], sOut[k]);
            }
        }
        else
        {
            const scalar sCut = -f0/dm;

            if (dm > 0)
            {
                add(max(sIn[k], sCut), sOut[k]);    // Liquid below the cut
            }
            else
            {
                add(sIn[k], min(sOut[k], sCut));    // Liquid above the cut
            }
        }
    }

    parts.transfer(found);
}


void Foam::functionObjects::surfacePosition::lineHeight
(
    const label li,
    const UList<List<List<Pair<scalar>>>>& allParts,
    scalar& height,
    scalar& meshLength
) const
{
    const scalarList& breaks = breaks_[li];
    const labelList& cover = cover_[li];

    height = 0;
    meshLength = 0;
    forAll(cover, j)
    {
        if (cover[j])
        {
            meshLength += breaks[j + 1] - breaks[j];
        }
    }

    if (meshLength <= 0)
    {
        return;
    }

    if (mode_ == modeType::ACCUMULATED)
    {
        // Each part in the alpha phase adds its overlap with each piece of
        // the line divided by the number of cells covering that piece
        scalar alphaLength = 0;

        for (const auto& proc : allParts)
        {
            for (const Pair<scalar>& part : proc[li])
            {
                label j = max(findLower(breaks, part.first()), 0);

                for (; j < cover.size() && breaks[j] < part.second(); ++j)
                {
                    if (cover[j])
                    {
                        const scalar overlap =
                            min(part.second(), breaks[j + 1])
                          - max(part.first(), breaks[j]);

                        if (overlap > 0)
                        {
                            alphaLength += overlap/cover[j];
                        }
                    }
                }
            }
        }

        alphaLength = min(alphaLength, meshLength);
        height = liquid_ ? alphaLength : meshLength - alphaLength;
        return;
    }

    // Highest or lowest: the line is divided at the ends of the pieces and of
    // the parts in the alpha phase. Each interval is in the liquid if at least
    // half of the cells covering it are. Sweeping down the line (increasing
    // s), the number of alpha parts covering an interval is updated from the
    // part ends.
    DynamicList<scalar> ends(breaks);
    DynamicList<Pair<scalar>> events;   // (s, +1 or -1)

    for (const auto& proc : allParts)
    {
        for (const Pair<scalar>& part : proc[li])
        {
            pushBack(ends, part.first());
            pushBack(ends, part.second());
            pushBack(events, Pair<scalar>(part.first(), 1));
            pushBack(events, Pair<scalar>(part.second(), -1));
        }
    }

    Foam::sort(ends);
    Foam::sort
    (
        events,
        [](const Pair<scalar>& x, const Pair<scalar>& y)
        {
            return x.first() < y.first();
        }
    );

    const scalar sLow = breaks.last();
    scalar sTopFirst = -1;     // Top of the highest liquid run
    scalar sTopLast = -1;      // Top of the lowest liquid run
    bool inLiquid = false;
    label nAlpha = 0;
    label ei = 0;

    for (label i = 0; i + 1 < ends.size(); ++i)
    {
        const scalar s0 = ends[i];
        const scalar s1 = ends[i + 1];

        while (ei < events.size() && events[ei].first() <= s0)
        {
            nAlpha += label(events[ei].second());
            ++ei;
        }

        if (s1 <= s0)
        {
            continue;
        }

        const label j = findLower(breaks, 0.5*(s0 + s1));
        const label n = (j >= 0 && j < cover.size()) ? cover[j] : 0;

        bool liquid = false;
        if (n > 0)
        {
            const label nLiquid = liquid_ ? nAlpha : n - nAlpha;
            liquid = (2*nLiquid >= n);
        }

        if (liquid && !inLiquid)
        {
            if (sTopFirst < 0)
            {
                sTopFirst = s0;
            }
            sTopLast = s0;
        }
        inLiquid = liquid;
    }

    const scalar sTop =
        (mode_ == modeType::HIGHEST) ? sTopFirst : sTopLast;

    if (sTop >= 0)
    {
        height = sLow - sTop;
    }
}


void Foam::functionObjects::surfacePosition::writeFileHeader
(
    const fileID fid
)
{
    // As interfaceHeight, with the mode
    forAll(locations_, li)
    {
        writeHeaderValue
        (
            files(label(fid)),
            "Location " + Foam::name(li),
            locations_[li]
        );
    }

    writeHeaderValue(files(label(fid)), "Mode", modeTypeNames_[mode_]);

    switch (fid)
    {
        case fileID::heightFile:
        {
            writeHeaderValue
            (
                files(label(fid)),
                "hB",
                "Interface height above the boundary"
            );
            writeHeaderValue
            (
                files(label(fid)),
                "hL",
                "Interface height above the location"
            );
            break;
        }
        case fileID::positionFile:
        {
            writeHeaderValue(files(label(fid)), "p", "Interface position");
            break;
        }
    }

    const Foam::Omanip<int> w = valueWidth(1);

    writeCommented(files(label(fid)), "Location");
    forAll(locations_, li)
    {
        switch (fid)
        {
            case fileID::heightFile:
                files(label(fid)) << w << li << w << ' ';
                break;
            case fileID::positionFile:
                files(label(fid)) << w << li << w << ' ' << w << ' ' << "  ";
                break;
        }
    }
    files(label(fid)).endl();

    writeCommented(files(label(fid)), "Time");
    forAll(locations_, li)
    {
        switch (fid)
        {
            case fileID::heightFile:
                files(label(fid)) << w << "hB" << w << "hL";
                break;
            case fileID::positionFile:
                files(label(fid)) << w << "p" << w << ' ' << w << ' ' << "  ";
                break;
        }
    }
    files(label(fid)).endl();
}


void Foam::functionObjects::surfacePosition::writePositions()
{
    if (update_)
    {
        calcLineCells();
    }

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    // Point values of alpha for the isosurface cuts. The interpolation covers
    // the whole mesh and is collective, so it is done on all processors.
    const tmp<pointScalarField> tpointAlpha
    (
        volPointInterpolation::New(mesh_).interpolate(alpha)
    );
    ap_ = tpointAlpha().primitiveField();

    if (!isoPtr_)
    {
        isoPtr_.reset(new surfaceIteratorIso(mesh_, ap_, surfCellTol_));
    }

    const label nLines = starts_.size();
    List<List<Pair<scalar>>> local(nLines);

    for (label li = 0; li < nLines; ++li)
    {
        liquidParts(li, alpha, local[li]);
    }

    const List<List<List<Pair<scalar>>>> allParts
    (
        gatherToMaster(local)
    );

    if (!Pstream::master())
    {
        return;
    }

    files(label(fileID::heightFile)) << mesh_.time().timeName() << tab;
    files(label(fileID::positionFile)) << mesh_.time().timeName() << tab;

    const vector up(-down_);

    for (label li = 0; li < nLines; ++li)
    {
        scalar height = 0;
        scalar meshLength = 0;
        lineHeight(li, allParts, height, meshLength);

        // Height of the location above the lower boundary, where the line
        // leaves the mesh (as interfaceHeight: -GREAT if it misses the mesh)
        scalar hLB = -GREAT;

        if (meshLength > 0)
        {
            const point pLow(starts_[li] + breaks_[li].last()*down_);
            hLB = up & (locations_[li] - pLow);
        }

        // Interface heights above the boundary and location
        const scalar hIB = height;
        const scalar hIL = hIB - hLB;

        // Position of the interface
        const point p = locations_[li] + up*hIL;

        const Foam::Omanip<int> w = valueWidth(1);

        files(label(fileID::heightFile)) << w << hIB << w << hIL;
        files(label(fileID::positionFile))
            << '(' << w << p.x() << w << p.y()
            << valueWidth() << p.z() << ") ";
    }

    files(label(fileID::heightFile)).endl();
    files(label(fileID::positionFile)).endl();
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
    logFiles(obr_, name),
    alphaName_("alpha"),
    liquid_(true),
    mode_(modeType::ACCUMULATED),
    direction_(Zero),
    locations_(),
    surfCellTol_(1e-8),
    isoFaceTol_(1e-8),
    checkLineCells_(false),
    fullSearchInterval_(0),
    nRebuilds_(0),
    down_(Zero),
    starts_(),
    length_(0),
    lineCells_(),
    lineIn_(),
    lineOut_(),
    breaks_(),
    cover_(),
    seeds_(),
    searchPtr_(nullptr),
    searchStale_(false),
    update_(true),
    ap_(),
    isoPtr_(nullptr)
{
    read(dict);
    resetNames({"height", "position"});

    if (Pstream::master())
    {
        writeFileHeader(fileID::heightFile);
        writeFileHeader(fileID::positionFile);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::surfacePosition::~surfacePosition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::surfacePosition::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    alphaName_ = dict.getOrDefault<word>("alpha", "alpha");
    liquid_ = dict.getOrDefault("liquid", true);
    mode_ = modeTypeNames_.getOrDefault("mode", dict, modeType::ACCUMULATED);
    dict.readEntry("locations", locations_);
    direction_ = dict.getOrDefault<vector>("direction", Zero);
    surfCellTol_ = dict.getOrDefault<scalar>("surfCellTol", 1e-8);
    isoFaceTol_ = dict.getOrDefault<scalar>("isoFaceTol", 1e-8);
    checkLineCells_ = dict.getOrDefault("checkLineCells", false);
    fullSearchInterval_ = dict.getOrDefault<label>("fullSearchInterval", 0);

    setLines();

    seeds_.clear();
    update_ = true;
    isoPtr_.reset(nullptr);

    return true;
}


bool Foam::functionObjects::surfacePosition::execute()
{
    return true;
}


bool Foam::functionObjects::surfacePosition::write()
{
    logFiles::write();

    writePositions();

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
                pushBack(newCells, newCelli);
            }
        }

        seeds_[li].transfer(newCells);
    }

    update_ = true;
    searchStale_ = true;
    isoPtr_.reset(nullptr);
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
