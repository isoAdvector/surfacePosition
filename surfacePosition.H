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

Class
    Foam::functionObjects::surfacePosition

Group
    grpFieldFunctionObjects

Description
    Interface height and position at a set of locations, for VoF solvers,
    with the same entries and output as interfaceHeight, at a cost per write
    that does not depend on the mesh size.

    Through each location a line runs along the direction (default: gravity)
    across the whole mesh, wherever the location is. The height of the liquid
    along the line is the sum of the lengths of the line inside the liquid:
    the whole segment in full cells, nothing in empty cells, and in interface
    cells (surfCellTol < alpha < 1 - surfCellTol) the part on the liquid side
    of the isosurface of alpha, interpolated to the mesh points, that cuts the
    cell's volume fraction (surfaceIteratorIso::vofCutCell, as in the isoAlpha
    reconstruction). Where several cells cover the same part of the line (a
    line in a face plane or on an edge), their liquid lengths are averaged.

    The mode sets which height is reported:
    - accumulated (default): the length of the line in the liquid, as
      interfaceHeight;
    - highest: the height of the highest point of the liquid along the line,
      e.g. including droplets above the surface;
    - lowest: the height of the top of the lowest liquid part along the line,
      i.e. the first transition from liquid to gas going up the line.
    For highest and lowest, a part of the line covered by several cells is in
    the liquid if at least half of them are.

    Written, as by interfaceHeight, for each location:
    - hB: the height (of the mode) above the lower boundary, i.e. where the
      line leaves the mesh along the direction
    - hL: the height of the interface above the location, hB minus the
      height of the location above the lower boundary
    - p: the position of the interface, the location moved by hL against the
      direction
    A dry line gives hB = 0 (the interface at the lower boundary), a line
    missing the mesh hB = 0 and hL = GREAT.

    The cells crossed by each line are found once with a cached meshSearch.
    After a mesh change they are found again by a flood fill from the old
    line cells (mapped to the new mesh after a topology change), testing each
    candidate cell for intersection with the line, so that no search tree is
    rebuilt. The flood fill continues across processor faces, so a line that
    moves into a neighbouring processor subdomain is followed. The segments of
    the line in its cells, and how many cells cover each part of the line,
    are computed only then, so that full cells cost nothing at a write.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    surfacePosition1
    {
        // Mandatory entries
        type            surfacePosition;
        libs            (surfacePosition);
        locations       ((0 0 0) (10 0 0) (20 0 0));

        // Optional entries, as for interfaceHeight
        alpha           alpha.water;    // default: alpha
        liquid          true;           // alpha is the liquid's
        direction       (0 0 -1);       // default: gravity

        // Optional entries
        mode            accumulated;    // accumulated, highest or lowest
        surfCellTol     1e-8;           // interface cells
        isoFaceTol      1e-8;           // tolerance of the isosurface cut

        // Optional checks of the line cells after mesh changes
        checkLineCells      false;  // compare with all cells (tests only)
        fullSearchInterval  0;      // full search every N rebuilds (0: never)
    }
    \endverbatim

    Output files, as for interfaceHeight:
    \table
      Operand        | Type | Location
      output file 1  | dat  | postProcessing/\<FO\>/\<time\>/height
      output file 2  | dat  | postProcessing/\<FO\>/\<time\>/position
    \endtable
    The headers also state the mode (# Mode : ...).

Note
    - The point interpolation of alpha covers the whole mesh once per write.
    - The flood fill follows face neighbours and assumes convex cells. If no
      old line cell crosses the line any more, rings of neighbours (up to 10)
      are searched, also across processor faces. A line that moves further
      than that, or into a part of the mesh connected to its old cells only
      through cyclic faces, is only found by a full search
      (fullSearchInterval).
    - A line crossing the mesh in several parts (e.g. below an overhang)
      counts the liquid of all parts from the lowest exit, as
      interfaceHeight does.
    - The library links OpenFOAM's geometricVoF. With a solver linking a
      modified copy of it, both are loaded.

SourceFiles
    surfacePosition.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfacePosition_H
#define functionObjects_surfacePosition_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "meshSearch.H"
#include "volFieldsFwd.H"
#include "boundBox.H"
#include "scalarField.H"
#include "Pair.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;
class surfaceIteratorIso;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class surfacePosition Declaration
\*---------------------------------------------------------------------------*/

class surfacePosition
:
    public fvMeshFunctionObject,
    public logFiles
{
public:

    // Public Data Types

        //- Which height along a line is reported
        enum class modeType
        {
            ACCUMULATED,    //!< Length of the line in the liquid
            HIGHEST,        //!< Highest point of the liquid
            LOWEST          //!< Top of the lowest liquid part
        };

        //- Names of the modes
        static const Enum<modeType> modeTypeNames_;


private:

    // Private Enumerations

        //- Output files, as for interfaceHeight
        enum class fileID
        {
            heightFile = 0,
            positionFile = 1
        };


    // Private Data

        //- Name of the volume fraction field
        word alphaName_;

        //- Is the alpha field that of the liquid?
        bool liquid_;

        //- Which height is reported
        modeType mode_;

        //- Direction of the lines, downwards (zero: gravity)
        vector direction_;

        //- Locations to report the height at
        List<point> locations_;

        //- Tolerance for interface cells:
        //- surfCellTol < alpha < 1 - surfCellTol
        scalar surfCellTol_;

        //- Tolerance on the volume fraction of the isosurface cut
        scalar isoFaceTol_;

        //- Maximum number of iterations of the isosurface cut
        static constexpr label maxIter_ = 100;

        //- Compare the line cells after each rebuild with all cells
        bool checkLineCells_;

        //- Full search every this many rebuilds (0: never)
        label fullSearchInterval_;

        //- Number of rebuilds after mesh changes
        label nRebuilds_;

        //- Unit downward direction of the lines
        vector down_;

        //- Line start points, above the mesh
        List<point> starts_;

        //- Line length, across the mesh with a margin
        scalar length_;

        //- Cells crossed by each line on this processor, ordered along it
        List<labelList> lineCells_;

        //- Distance along the line where it enters and leaves those cells
        List<scalarList> lineIn_;
        List<scalarList> lineOut_;

        //- Parts of each line between the ends of the segments of all its
        //- cells, on all processors (master only)
        List<scalarList> breaks_;

        //- Number of cells covering each part (master only)
        List<labelList> cover_;

        //- Seed cells for the next rebuild (empty: full search)
        List<labelList> seeds_;

        //- Cached search engine
        autoPtr<meshSearch> searchPtr_;

        //- Bounds of the mesh on this processor (for full searches)
        boundBox localBounds_;

        //- The search engine needs correcting after a mesh change
        bool searchStale_;

        //- The line cells must be (re)built before the next write
        bool update_;

        //- Values of alpha at the mesh points, for the isosurface cuts
        scalarField ap_;

        //- Isosurface cuts of cells, on ap_
        autoPtr<surfaceIteratorIso> isoPtr_;


    // Private Member Functions

        //- Set the downward direction and the lines
        void setLines();

        //- Does the line cross the (convex) cell? Also returns the distances
        //- along the line where it enters and leaves the cell.
        bool crossesLine
        (
            const label celli,
            const label li,
            scalar& sIn,
            scalar& sOut
        ) const;

        //- Does the line cross the (convex) cell?
        bool crossesLine(const label celli, const label li) const
        {
            scalar sIn, sOut;
            return crossesLine(celli, li, sIn, sOut);
        }

        //- Seed cells for a full search of the line, with the search engine
        void searchSeeds(const label li, labelList& seeds);

        //- Cells crossed by each line, by a flood fill from seed cells that
        //- continues across processor faces
        void floodLineCells
        (
            const UList<labelList>& seeds,
            List<labelList>& cells
        ) const;

        //- Order the cells along the line and store their segments
        void setLineCells(const label li, const labelUList& cells);

        //- Parts of the lines and their cover from the segments of the line
        //- cells on all processors
        void calcCover();

        //- Find the cells crossed by the lines
        void calcLineCells();

        //- The parts of line li in the alpha phase in its cells on this
        //- processor, touching parts merged
        void liquidParts
        (
            const label li,
            const volScalarField& alpha,
            List<Pair<scalar>>& parts
        );

        //- Height of the mode above the lower boundary along line li, from
        //- the parts of the line in the alpha phase on all processors, and
        //- the length of the line in the mesh (master only)
        void lineHeight
        (
            const label li,
            const UList<List<List<Pair<scalar>>>>& allParts,
            scalar& height,
            scalar& meshLength
        ) const;

        //- Output file header, as for interfaceHeight
        void writeFileHeader(const fileID fid);

        //- Write the heights and positions
        void writePositions();

        //- No copy construct
        surfacePosition(const surfacePosition&) = delete;

        //- No copy assignment
        void operator=(const surfacePosition&) = delete;


public:

    //- Runtime type information
    TypeName("surfacePosition");


    // Constructors

        //- Construct from name, Time and dictionary
        surfacePosition
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfacePosition();


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary& dict);

        //- Do nothing
        virtual bool execute();

        //- Write the heights and positions
        virtual bool write();

        //- Rebuild the line cells after a topology change
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Rebuild the line cells after mesh motion
        virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
