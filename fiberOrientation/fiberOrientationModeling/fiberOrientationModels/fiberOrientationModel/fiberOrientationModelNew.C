/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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

#include "fiberOrientationModel.H"

// ************************************************************************* //

Foam::autoPtr<Foam::fiberOrientation::fiberOrientationModel> 
Foam::fiberOrientation::fiberOrientationModel::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const volVectorField& U,
    const surfaceScalarField& phi
)
{
    const word modelType(dict.subDict("fiberOrientationProperties").lookup("model"));

    Info << nl << "Selected fiber orientation model: " << modelType << endl;

    //auto* ctorPtr = dictionaryConstructorTable(modelType);
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "fiberOrientation::New"
            "(const dictionary&, const fvMesh&, const volVectorField&, const surfaceScalarField&)"
        )
            << "Unknown fiberOrientation type "
            << modelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid closureModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<fiberOrientation::fiberOrientationModel >(cstrIter()(dict, mesh, U, phi));

}

// ************************************************************************* //
