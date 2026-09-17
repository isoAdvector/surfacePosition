/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::functionObjects::surfacePosition

Group
    grpFieldFunctionObjects

Description
    Position of the interface along user-specified lines, for VoF solvers, at
    a cost per write that is independent of the mesh size.

    The cells crossed by each line are found once with a cached meshSearch.
    After a mesh change they are found again by a flood fill from the old
    line cells (mapped to the new mesh after a topology change), testing each
    candidate cell for intersection with the line, so that no search tree is
    rebuilt. The flood fill continues across processor faces, so a line that
    moves into a neighbouring processor subdomain is followed.

    At each write, the cells of each line are visited in order along the line,
    and the first interface cell whose interface meets the line inside the
    cell is used (or, if there is none, the first interface cell). The
    position is found by one of the methods:
    - reconstructed: intersection of the line with the interface plane
      (fields interfaceNormal.\<phase\> and interfaceCentre.\<phase\>) of the
      reconstruction scheme of a geometric VoF solver.
    - plic: intersection with a plane cutting the cell at the volume fraction
      of the cell (piecewise linear interface calculation), with the Youngs
      normal from alpha interpolated to the cell points. Needs no
      reconstruction and works with any solver; the point interpolation of
      alpha covers the whole mesh once per write.
    - alpha: where alpha = 0.5, interpolated linearly between the centres of
      each cell and the cell before it along the line (also across processor
      and cyclic faces).
    - auto (default): reconstructed if the fields are registered, else plic.

    An interface cell is a cell with a non-zero reconstructed normal
    (reconstructed), or with alphaTol < alpha < 1 - alphaTol (plic). An
    interface lying on a face between two full cells (alpha within 1e-6 of 0
    and 1) is found as a jump of alpha where the line enters the cell, also
    across processor and cyclic faces.

    Several cells, or several processors, may find the same point, for
    example when the interface or the line lies on a face, edge or point of
    the mesh: only the distance along the line is reduced, with a minimum,
    so the point is reported once. Interfaces meeting the line inside their
    cell take precedence on all processors over interface planes that cut
    the line outside the cell.

    In parallel each processor finds the first interface along its part of
    each line, and a single reduction over all lines selects the first one
    along the whole line.

    Lines are given either as locations and a direction, or as explicit start
    and end points. With locations, each line runs through the location along
    the direction, from above to below across the whole mesh (wherever the
    location is, also outside the mesh), so the first interface found is the
    uppermost one, and the height of the interface above the location is also
    written. The direction defaults to the opposite of gravity.

    Output file postProcessing/\<name\>/\<time\>/surfacePosition.dat with, per
    line, the distance s from the line start, the height h above the location
    (locations only), and the position x y z of the interface. A line without
    an interface is written as nan.

Usage
    \verbatim
    surfacePosition1
    {
        type            surfacePosition;
        libs            (surfacePosition);

        // Lines through locations along a direction
        locations       ((1.0001 0.0001 0.0001));
        direction       (0 0 1);            // optional, default -g

        // Or explicit lines; the first interface from the start is used
        // lines        (((1 0 0) (1 0 1)));

        alpha           alpha.water;        // optional
        method          auto;               // optional
        alphaTol        1e-8;               // optional

        // Optional checks of the line cells after mesh changes
        checkLineCells      false;  // compare with all cells (tests only)
        fullSearchInterval  0;      // full search every N rebuilds (0: never)

        writeControl    timeStep;
        writeInterval   1;
    }
    \endverbatim

Note
    - Only the first interface along each line is reported.
    - The flood fill follows face neighbours and assumes convex cells. If no
      old line cell crosses the line any more, rings of neighbours (up to 10)
      are searched, also across processor faces. A line that moves further
      than that, or into a part of the mesh connected to its old cells only
      through cyclic faces, is only found by a full search
      (fullSearchInterval).

SourceFiles
    surfacePosition.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfacePosition_H
#define functionObjects_surfacePosition_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "meshSearch.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "Enum.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class surfacePosition Declaration
\*---------------------------------------------------------------------------*/

class surfacePosition
:
    public fvMeshFunctionObject,
    public writeFile
{
public:

    // Public Data Types

        //- Method for the interface position in a cell
        enum class methodType
        {
            AUTO,
            RECONSTRUCTED,
            PLIC,
            ALPHA
        };

        //- Names of the methods
        static const Enum<methodType> methodTypeNames_;


private:

    // Private Data

        //- Name of the volume fraction field
        word alphaName_;

        //- Requested method
        methodType method_;

        //- Method used (method_ with AUTO resolved)
        methodType usedMethod_;

        //- Volume fraction tolerance for interface cells (plic)
        scalar alphaTol_;

        //- Volume fraction tolerance for full cells at a sharp jump
        static constexpr scalar jumpTol_ = 1e-6;

        //- Compare the line cells after each rebuild with all cells
        bool checkLineCells_;

        //- Full search every this many rebuilds (0: never)
        label fullSearchInterval_;

        //- Number of rebuilds after mesh changes
        label nRebuilds_;

        //- Lines are given as locations and a direction
        bool useLocations_;

        //- Locations (locations mode only)
        List<point> locations_;

        //- Unit upward direction (locations mode only)
        vector up_;

        //- Line start points
        List<point> starts_;

        //- Unit line directions, from start to end
        List<vector> dirs_;

        //- Line lengths
        List<scalar> lengths_;

        //- Cells crossed by each line on this processor, ordered along it
        List<labelList> lineCells_;

        //- Distance of those cell centres along the line
        List<scalarList> lineDist_;

        //- Distance along the line where it enters and leaves those cells
        List<scalarList> lineIn_;
        List<scalarList> lineOut_;

        //- Face where the line enters those cells (-1: the line starts
        //- in the cell)
        List<labelList> lineEntry_;

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

        //- The file header has been written
        bool headerDone_;


    // Private Member Functions

        //- Does the line cross the (convex) cell? Also returns the distances
        //- along the line where it enters and leaves the cell.
        bool crossesLine
        (
            const label celli,
            const label li,
            scalar& sIn,
            scalar& sOut,
            label& entryFace
        ) const;

        //- Does the line cross the (convex) cell?
        bool crossesLine(const label celli, const label li) const
        {
            scalar sIn, sOut;
            label entryFace;
            return crossesLine(celli, li, sIn, sOut, entryFace);
        }

        //- Is distance s along line li within the k-th cell of the line?
        bool insideCell(const label li, const label k, const scalar s) const
        {
            const scalar sIn = lineIn_[li][k];
            const scalar sOut = lineOut_[li][k];
            const scalar tol = 1e-6*(sOut - sIn);
            return (s >= sIn - tol && s <= sOut + tol);
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

        //- Order the cells along the line and store them
        void setLineCells(const label li, const labelUList& cells);

        //- Find the cells crossed by the lines
        void calcLineCells();

        //- Alpha (and optionally the distance of the centre along the line)
        //- of the cell before the k-th cell of line li, across the face
        //- where the line enters; false at the line start or where the line
        //- enters the domain
        bool predecessor
        (
            const label li,
            const label k,
            const volScalarField& alpha,
            scalar& value,
            scalar* distPtr = nullptr
        ) const;

        //- First interface along line li on this processor: the distance to
        //- the first interface meeting the line inside its cell or sharp
        //- jump of alpha (event), and to the first interface plane cutting
        //- the line outside its cell (fallback); GREAT if none. With
        //- normalPtr the reconstructed planes are used, else PLIC.
        void firstInterface
        (
            const label li,
            const volScalarField& alpha,
            const volVectorField* normalPtr,
            const volVectorField* centrePtr,
            const pointScalarField* pointAlphaPtr,
            scalar& event,
            scalar& fallback
        ) const;

        //- Distance along line li to the first alpha = 0.5 crossing
        scalar alphaDistance
        (
            const label li,
            const volScalarField& alpha
        ) const;

        //- Gradient of alpha in a cell (Youngs), times the cell volume
        vector alphaGradient
        (
            const label celli,
            const pointScalarField& pointAlpha
        ) const;

        //- Volume of the part of a cell with (n & x) < d
        scalar volumeBelow
        (
            const label celli,
            const vector& n,
            const scalar d
        ) const;

        //- Plane position d with (n & x) < d filling the fraction alpha1
        //- of the cell
        scalar plicPlane
        (
            const label celli,
            const scalar alpha1,
            const vector& n
        ) const;

        //- Write the output file header
        void writeFileHeader(Ostream& os);

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
    virtual ~surfacePosition() = default;


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary& dict);

        //- Do nothing
        virtual bool execute();

        //- Find the interface positions and write them
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
