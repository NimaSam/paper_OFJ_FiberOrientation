/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::functionObjects::fiberOrientationModeling

SourceFiles
    fiberOrientationModeling.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fiberOrientationModeling_H
#define functionObjects_fiberOrientationModeling_H

#include "functionObject.H"
#include "volFields.H"
#include "fiberOrientationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class fiberOrientationModeling Declaration
\*---------------------------------------------------------------------------*/

class fiberOrientationModeling
:
    public functionObject
{
    // Private data
        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        const fvMesh& mesh_;

        //- Name of phase field (optional)
        word phaseName_;

        //- Name of flux field (optional)
        word phiName_;
        
        // Reference to the velocity field
        const volVectorField& U_;

         // Reference to the flux
        const surfaceScalarField& phi_;

        // Fiber Orientation model
        autoPtr<fiberOrientation::fiberOrientationModel> fiberOrientationModel_;

    // Private Member Functions

        //- No copy construct
        fiberOrientationModeling(const fiberOrientationModeling&) = delete;

        //- No copy assignment
        void operator=(const fiberOrientationModeling&) = delete;


public:

    //- Runtime type information
    TypeName("fiberOrientationModeling");


    // Constructors

        //- Construct from Time and dictionary
        fiberOrientationModeling
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~fiberOrientationModeling();

    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- Calculate the fiberOrientationModeling
        virtual bool execute(const bool forceWrite);

        //- Read the fiberOrientationModeling data
        virtual bool read(const dictionary&);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for mesh point-motion
        virtual void movePoints(const pointField&)
        {}

         //- Do nothing.
        //  The volScalarField is registered and written automatically
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
