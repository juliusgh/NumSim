/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "PatchFunction1.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 0f9d725e99189ec8a8cae2c9f53d1a1053c4d270
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void parabolicVelocity_0f9d725e99189ec8a8cae2c9f53d1a1053c4d270(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    parabolicVelocityFixedValueFvPatchVectorField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct parabolicVelocity : patch/DimensionedField");
    }
}


Foam::
parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const parabolicVelocityFixedValueFvPatchVectorField& rhs,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct parabolicVelocity : patch/DimensionedField/mapper");
    }
}


Foam::
parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (false)
    {
        printMessage("Construct parabolicVelocity : patch/dictionary");
    }
}


Foam::
parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const parabolicVelocityFixedValueFvPatchVectorField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct parabolicVelocity");
    }
}


Foam::
parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const parabolicVelocityFixedValueFvPatchVectorField& rhs,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct parabolicVelocity : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
parabolicVelocityFixedValueFvPatchVectorField::
~parabolicVelocityFixedValueFvPatchVectorField()
{
    if (false)
    {
        printMessage("Destroy parabolicVelocity");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
parabolicVelocityFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs parabolicVelocity");
    }

//{{{ begin code
    #line 32 "/home/julius/Documents/NumSim/NumSim/task3_openfoam/karmanvortexstreet/0/U.boundaryField.inlet"
const vectorField& Cf = patch().Cf();
            vectorField& field = *this;

            const scalar Umax = 1.0;
            forAll(Cf, faceI)
            {
                const scalar y = Cf[faceI][1];
                field[faceI] = vector(Umax*(1-(((y-0.5)/0.5)*((y-0.5)/0.5))),0,0);
            }
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

