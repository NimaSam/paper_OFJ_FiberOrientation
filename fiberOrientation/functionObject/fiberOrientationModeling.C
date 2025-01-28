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

\*---------------------------------------------------------------------------*/

#include "fiberOrientationModeling.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fiberOrientationModeling, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fiberOrientationModeling,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fiberOrientationModeling::fiberOrientationModeling
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    regionName_(polyMesh::defaultRegion),
    mesh_(time_.lookupObject<fvMesh>(regionName_)),
    phaseName_(dict.lookupOrDefault<word>("phase", "none")),
    phiName_(   phaseName_ == "none" 
                ? 
                dict.lookupOrDefault<word>("phi", "phi"):
                dict.lookupOrDefault<word>("phasePhiCompressed", "alphaPhiUn") ),
    U_( mesh_.lookupObject<volVectorField>("U") ),
    phi_( mesh_.lookupObject<surfaceScalarField>(phiName_) ),
    fiberOrientationModel_(fiberOrientation::fiberOrientationModel::New(dict, mesh_, U_, phi_))
{

    
    read(dict);

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        if (isA<emptyPolyPatch>(pp))
        {
            FatalErrorInFunction
                << "Fiber Orientation modeling is not ready for 2D cases!" << exit(FatalError);
        }
    }

    Info << nl
         << "FO is created" << nl
         << "Phase name is:\t" << phaseName_ << nl
         << "Phi name is:\t" << phiName_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fiberOrientationModeling::~fiberOrientationModeling()
{
    Info << "FiberOrientation FO was destroyed" << endl;
} 


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::fiberOrientationModeling::start()
{
    return true;
}

bool Foam::fiberOrientationModeling::read(const dictionary& dict)
{
    

    dict.readIfPresent("phase", phaseName_);

    if (phaseName_ != "none")
    {
        phiName_ = dict.lookupOrDefault<word>("phasePhiCompressed", "alphaPhiUn");
    }
    else
    {
        phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    }
    
    //functionObject::read(dict);

    return false;
}


bool Foam::fiberOrientationModeling::execute(const bool forceWrite)
{
    Info << type() << ":\tSolving fiber orientation:" << endl;

    fiberOrientationModel_->solve();

    Info << endl;
    
    return true;
}

bool Foam::fiberOrientationModeling::write()
{
    return true;
}


// ************************************************************************* //
